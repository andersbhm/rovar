import os, time
from math import radians, cos, sin, asin, sqrt, degrees
from calendar import timegm

import numpy as np
import datetime
import pickle
from scipy.constants import codata
from netCDF4 import Dataset

from pathlib import Path
from spacetrack import SpaceTrackClient
import spacetrack.operators as op
from astropy.utils.iers import IERS_A
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from astropy.time import Time
from astropy import units as u

import pyproj
geo_WGS84 = pyproj.Geod(ellps='WGS84')
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')




aepoch = time.strptime('2004 1 1 0 0 0', '%Y %m %d %H %M %S') # epoch of AGILE satellite clock
tepoch = timegm(aepoch)
################################################################################
R_earth = 6378137. # meters, radius of the Earth
speed_of_light = codata.value('speed of light in vacuum')
################################################################################
def check_if_folder_exist_or_make_it(directory):
    """If directory does not exist - > make it"""
    if not os.path.exists(directory):
        os.makedirs(directory)
    return 0

def check_order_of_magnitude(number):
    '''Return order of magnitude of number.'''
    return int(np.log10(number))

def get_AGILE_obt_from_date(date):
    '''Convert UTC time to onboardtime AGILE satellite'''
    return timegm(date.timetuple()) - tepoch

def haversine(lat1, lon1, lat2, lon2):
    '''Calculate great circle distance [m] and angle [radians] between coordinates in latitude, and longitude, ((lat1, lon1) and (lat2, lon2)), using haversine'''
    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    theta = 2*asin(sqrt(sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2))
    dist = R_earth*theta
    return dist, theta

def get_propagation_time_haversine(lat1, lon1, altitude1, lat2, lon2, altitude2):
    '''
    Get propagation time, assuming speed of light, from two different (lat, lon, altitude) positions
    Latitude and longitude in degrees. Altitude in meters.
    Return propagation time in seconds'''
    _, theta = haversine(lat1, lon1, lat2, lon2)
    distance_lightning_sat = sqrt((R_earth+altitude1)**2 + (R_earth+altitude2)**2 - 2*(R_earth+altitude1)*(R_earth+altitude2)*cos(theta));
    return distance_lightning_sat/speed_of_light

def get_propagation_time_and_distance(lat1, lon1, altitude1, lat2, lon2, altitude2, calculation_method="WGS84"):
    ''' More accurate than haversine and get_propagation_time functions. Do not assume Earth as a perfect circle.

    ###Inputs Can be arrays.###

    Input units in degrees and meters.
    Calculation method can be "WGS84" (faster and more accurate) or "HAVERSINE"

    Return propagation time between position1 and position2 in seconds
    and great circle distance along surface of Earth, between position1 and position2 in meters.
    '''
    if calculation_method == "WGS84": # faster an more accurate
        x_sat, y_sat, z_sat = pyproj.transform(lla, ecef, lon1, lat1, altitude1, radians=False)
        x_lightning, y_lightning, z_lightning = pyproj.transform(lla, ecef, lon2, lat2, altitude2, radians=False)
        x_sat, y_sat, z_sat = np.array(x_sat), np.array(y_sat), np.array(z_sat)
        x_lightning, y_lightning, z_lightning = np.array(x_lightning), np.array(y_lightning), np.array(z_lightning)
        distance_between_lightning_sat = np.sqrt( (x_lightning-x_sat)**2 + (y_lightning-y_sat)**2 + (z_lightning-z_sat)**2)
        propagation_time_lightning_sat = distance_between_lightning_sat/speed_of_light
        _, _, distance = geo_WGS84.inv(lon1,lat1,lon2,lat2)

    if calculation_method == "HAVERSINE":
        propagation_time_lightning_sat = []
        distance = []
        for lats, lons, alt, latl, lonl, altl in zip(lat1,lon1,altitude1,lat2,lon2, altitude2):
            propagation_time_lightning_sat_0 = get_propagation_time_haversine(lats, lons, alt, latl, lonl, altl)
            propagation_time_lightning_sat.append(propagation_time_lightning_sat_0)
            dist, theta = haversine(lats, lons, latl, lonl)
            distance.append(dist)
        propagation_time_lightning_sat = np.array(propagation_time_lightning_sat)
        distance = np.array(distance)

    return propagation_time_lightning_sat, distance

def distance_to_shoreline(lon_array, lat_array, path_to_GSHHG = "./useful_files/dist_to_GSHHG_v2.3.7_1m.nc"):
    '''Get distance to coastline based on the GSHHG precalculated distance to shore grid: dist_to_GSHHG_v2.3.7_1m.nc
    The datset is available from https://www.soest.hawaii.edu/pwessel/gshhg/
    Input: longitude and latitude arrays.
    Return: array of distances to coastline in meters, where negative distances are over ocean.
    '''
    lon_array, lat_array = np.array(lon_array), np.array(lat_array)
    lon_lat_distance_to_coast_dataset = Dataset(path_to_GSHHG)
    lon_grid = np.array(lon_lat_distance_to_coast_dataset.variables['lon']) # 0 - 360
    lon_grid = ((lon_grid+180) % 360) - 180 # convert from 360 to +-180
    lat_grid = np.array(lon_lat_distance_to_coast_dataset.variables['lat'])
    dist_grid = np.array(lon_lat_distance_to_coast_dataset.variables['dist'])
    dist_coast = []
    for lon, lat in zip(lon_array, lat_array):
        index_lon = (np.abs(lon_grid - lon)).argmin()
        index_lat = (np.abs(lat_grid - lat)).argmin()
        distance = dist_grid[index_lat, index_lon] # negative over ocean
        dist_coast.append(distance)
    return np.array(dist_coast)*1000




def get_lat_lon_h_satellite(year, month, day, hour, minute, second, microsecond, satellite_name,
    printwarning=True, updateTLE=True,
    update_TLE_information_array = [25544, "./useful_files/spacetrackclient_usernamepw.txt"],
    file_finals2000A_path = "./useful_files/finals2000A.all.txt"):

    '''
    year month, day hour minute second microsecond are integers.
    (2018, 2, 16, 20, 57, 49, 123456, "ISS")
    Satellite_name can be "ISS", "RHESSI", "AGILE" or "Fermi", depends on
    The function return latitude, longitude in degrees, altitude in meters,
    and error in radial distance in km between the TLE before and after the closest TLE.
    The last parameter is a quality check if ISS has done some maneuvers.

    Download Earth orientation data from: https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
    Filename = finals2000A.all

    Run update_TLEs() to get TLE data for this calculation.
    '''
    sat_ID, spacetrack_usernamepassword_file = update_TLE_information_array[0], update_TLE_information_array[1]

    def get_spacetrackclient_username_and_password(spacetrack_usernamepassword_file):
        f = open(spacetrack_usernamepassword_file).read().split(",")
        username, password = f[0], f[1]
        return username, password

    def update_TLEs(username, password, f, satid, printwarning=True):
        '''
        Automatic download TLE information from space-track.org
        Username: your username to space-track.org
        Password: your password to space-track.org
        printwarning: Print information
        file_list: your preferred filesave name, typically ["ISS_orbital_info.txt", "Fermi_orbital_info.txt", "AGILE_orbital_info.txt"]
        sat_ID_list: the corresponding satellited IDs: [25544, 33053, 31135, 27370]
        '''
        # Update only if the day has changed
        now = datetime.datetime.now()

        if printwarning:
            print(' ')
            print('Trying to update TLE data')

            print("Current year: %d" % now.year)
            print("Current month: %d" % now.month)
            print("Current day: %d" % now.day)
        today={}
        today['year'] = now.year
        today['month'] = now.month
        today['day'] = now.day

        directory = "dataFiles_satellite_TLE/"
        check_if_folder_exist_or_make_it(directory)

        TLE_nb_lim = 200000

        st = SpaceTrackClient(username, password)

        f = "%s%s" % (directory,f)
        if os.path.exists(f):
            os.remove(Path(f))

        lines = st.tle(norad_cat_id=[satid], orderby='epoch asc', format='tle',limit=TLE_nb_lim)

        with open(Path(f), 'w') as fp:
            for line in lines:
                fp.write(line)

        if printwarning:
            print('TLE data successfully updated.')
            print(' ')
        return 0

    def cbrt(x):
        if x >= 0:
            return pow(x, 1.0/3.0)
        else:
            return -pow(abs(x), 1.0/3.0)


    def ecef2geodetic(x, y, z):
        """Convert ECEF coordinates to geodetic.
        J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates \
        to geodetic coordinates," IEEE Transactions on Aerospace and \
        Electronic Systems, vol. 30, pp. 957-961, 1994."""

        # fprm: https://code.google.com/p/pysatel/source/browse/trunk/coord.py?r=22

        # Constants defined by the World Geodetic System 1984 (WGS84)
        a = 6378.137
        b = 6356.7523142
        esq = 6.69437999014 * 0.001
        e1sq = 6.73949674228 * 0.001
        f = 1 / 298.257223563

        r = sqrt(x * x + y * y)
        Esq = a * a - b * b
        F = 54 * b * b * z * z
        G = r * r + (1 - esq) * z * z - esq * Esq
        C = (esq * esq * F * r * r) / (pow(G, 3))
        S = cbrt(1 + C + sqrt(C * C + 2 * C))
        P = F / (3 * pow((S + 1 / S + 1), 2) * G * G)
        Q = sqrt(1 + 2 * esq * esq * P)
        r_0 =  -(P * esq * r) / (1 + Q) + sqrt(0.5 * a * a*(1 + 1.0 / Q) - \
            P * (1 - esq) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r)
        U = sqrt(pow((r - esq * r_0), 2) + z * z)
        V = sqrt(pow((r - esq * r_0), 2) + (1 - esq) * z * z)
        Z_0 = b * b * z / (a * V)
        h = U * (1 - b * b / (a * V))
        lat = np.arctan((z + e1sq * Z_0) / r)
        lon = np.arctan2(y, x)
        return degrees(lat), degrees(lon), h

    def rotation_z(a):
        # Rotation matrix around z axis
        # of an angle a (in rad)
        r = np.eye(3,3)
        cos_a = np.cos(a)
        sin_a = np.sin(a)
        r[0,0] = cos_a
        r[0,1] = sin_a
        r[1,0] = -sin_a
        r[1,1] = cos_a
        return r

    def get_lat_lon_h(l5,l6,wgs72, d):
        satellite0 = twoline2rv(l5, l6, wgs72)
        xsatellite, vsatellite = satellite0.propagate(d.year, d.month, d.day, d.hour, d.minute, d.second + d.microsecond * 1.e-6)
        x0 = np.array(xsatellite)
        X = x0[0]
        Y = x0[1]
        Z = x0[2]

        t = Time(d)
        t.delta_ut1_utc = iers_a.ut1_utc(t)
        s = t.sidereal_time('mean', 'greenwich');
        s.wrap_angle=180 * u.deg
        r = rotation_z(s)
        x0_ecef = r.dot(x0)
        lat, lon, h = ecef2geodetic(x0_ecef[0], x0_ecef[1], x0_ecef[2])
        return lat, lon, h

    ##########################
    if updateTLE:
        username, password = get_spacetrackclient_username_and_password(spacetrack_usernamepassword_file)
        f = "%s_orbital_info.txt" % satellite_name
        update_TLEs(username, password, f, sat_ID, printwarning=True)
    ##########################

    iers_a = IERS_A.open(file_finals2000A_path)

    dini = datetime.datetime(year, month, day, hour, minute, second, microsecond) # start date for the analysis

    # load TLE
    file_path = "./dataFiles_satellite_TLE/%s_orbital_info.txt" % satellite_name
    file_path = Path(file_path)

    satellite_tle_file = open(file_path, 'r')
    satellite_tle_list = satellite_tle_file.readlines()

    TLE_index_epoch = []
    l5_list, l6_list = [], []
    for i in range(0, len(satellite_tle_list), 2):
        l5 = satellite_tle_list[i]
        l6 = satellite_tle_list[i+1]
        satellite0 = twoline2rv(l5, l6, wgs72) # if this line failes try update sgp4
        TLE_index_epoch.append(get_AGILE_obt_from_date(satellite0.epoch))
        l5_list.append(l5)
        l6_list.append(l6)

    TLE_index_epoch = np.array(TLE_index_epoch)
    l5_list = np.array(l5_list)
    l6_list = np.array(l6_list)
    d = dini

    indx = np.argmin(abs(get_AGILE_obt_from_date(d) - TLE_index_epoch))
    lat_closest, lon_closest, h_closest = get_lat_lon_h(l5_list[indx],l6_list[indx],wgs72, d)
    lat_before, lon_before, h_before = get_lat_lon_h(l5_list[indx-1],l6_list[indx-1],wgs72, d)
    lat_after, lon_after, h_after = get_lat_lon_h(l5_list[indx+1],l6_list[indx+1],wgs72, d)


    dist1, _ = haversine(lat_closest, lon_closest, lat_before, lon_before)
    dist2, _ = haversine(lat_closest, lon_closest, lat_after, lon_after)
    dist1, dist2 = dist1/1000, dist2/1000 # meter to km
    diff = abs(dist1 - dist2)
    if (dist1 >= 5 or dist2 > 5) and printwarning:
        print("\n############## LARGE UNCERTAINTY IN POSITION. CHECK WHICH TLE TO USE ############")
        satellite0 = twoline2rv(l5_list[indx-1], l6_list[indx-1], wgs72)
        print("BEFORE: ", satellite0.epoch)
        satellite0 = twoline2rv(l5_list[indx+1], l6_list[indx+1], wgs72)
        print("AFTER:  ", satellite0.epoch)
        satellite0 = twoline2rv(l5_list[indx], l6_list[indx], wgs72)
        print("CLOSEST:", satellite0.epoch, " (THIS IS THE ONE USED)")
        print("Distance between TLE closest and TLE before closest: %.0f km" % (dist1))
        print("Distance between TLE closest and TLE after closest : %.0f km" % (dist2))

    lat, lon, h = lat_closest, lon_closest, h_closest
    return lat, lon, h*1000, dist1, dist2



if __name__ == '__main__':
    print("")
