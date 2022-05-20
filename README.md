# Useful functions
Useful python functions used for my PhD in Space Physics

# Description of functions


## get_propagation_time_and_distance(lat1, lon1, altitude1, lat2, lon2, altitude2, calculation_method="WGS84")
Get propagation time, assuming speed of light, and great circle distance from two different positions

***More accurate than haversine() and get_propagation_time() functions. Does not assume Earth as a perfect circle.***




* Parameters:
  * lat1 : latitude of first coordinate in degrees. Can be array.
  * lon1 : longitude of first coordinate in degrees. Can be array.
  * altitude1 : altitude of first coordinate in meters. Can be array.
  * lat2 : latitude of second coordinate in degrees. Can be array.
  * lon2 : longitude of second coordinate in degrees. Can be array.
  * altitude2 : altitude of second coordinate in meters. Can be array.
  * calculation_method : can be "WGS84" (faster and more accurate) or "HAVERSINE"

* Returns:    
  * prop_time: propagation time between position1 and position2 in seconds
  * distance: great circle distance along surface of Earth, between position1 and position2 in meters.

```python
>>> lat_lightning, lon_lightning, alt_lightning = 0.8033, 60.8386, 13000
>>> lat_satellite, lon_satellite, alt_satellite = 2.2478, 58.6766, 476210
>>> get_propagation_time_and_distance(lat_lightning, lon_lightning, alt_lightning, lat_satellite, lon_satellite, alt_satellite, calculation_method="WGS84")
(0.001840, 288776.470869)
>>>
>>>
>>> lat_lightning, lon_lightning, alt_lightning = [0.8033, 9.8014], [60.8386, -73.3565], [13000, 13000]
>>> lat_satellite, lon_satellite, alt_satellite = [2.2478, 1.9268], [58.6766, -74.3703], [476210, 477958]
>>> get_propagation_time_and_distance(lat_lightning, lon_lightning, alt_lightning, lat_satellite, lon_satellite, alt_satellite)
(array([0.001840, 0.003411]), [288776.470869, 878028.471848])
```






## haversine(lat1, lon1, lat2, lon2)
Calculate great circle distance and angle between two coordinates using the haversine formula
* Parameters:
  * lat1 : latitude of first coordinate in degrees
  * lon1 : longitude of first coordinate in degrees
  * lat2 : latitude of second coordinate in degrees
  * lon2 : longitude of second coordinate in degrees

* Returns:    
  * distance: great circle distance in meters
  * theta: angle in radians
```python
>>> haversine(0.51, -73.31, 0.15, -76.26)
(330822.772, 0.052)
```


## get_propagation_time_haversine(lat1, lon1, altitude1, lat2, lon2, altitude2)
Get propagation time, assuming speed of light, from two different positions


Latitude and longitude in degrees. Altitude in meters.
Return propagation time in seconds'''

* Parameters:
  * lat1 : latitude of first coordinate in degrees
  * lon1 : longitude of first coordinate in degrees
  * altitude1 : altitude of first coordinate in meters
  * lat2 : latitude of second coordinate in degrees
  * lon2 : longitude of second coordinate in degrees
  * altitude2 : altitude of second coordinate in meters

* Returns:    
  * prop_time: propagation time in seconds
```python
>>> lat_lightning, lon_lightning, alt_lightning = 0.8033, 60.8386, 13000
>>> lat_satellite, lon_satellite, alt_satellite = 2.2478, 58.6766, 476210
>>> get_propagation_time_haversine(lat_lightning, lon_lightning, alt_lightning, lat_satellite, lon_satellite, alt_satellite)
0.001841
```


## distance_to_shoreline(lon_array, lat_array, path_to_GSHHG)

Get distance to coastline based on the GSHHG precalculated distance to shore grid

Filename: dist_to_GSHHG_v2.3.7_1m.nc

The datset is available from https://www.soest.hawaii.edu/pwessel/gshhg/

* Parameters:
  * lon_array : longitude array in degrees
  * lat_array : latitude array in degrees
  * path_to_GSHHG : path to "dist_to_GSHHG_v2.3.7_1m.nc"

* Returns:    
  * distance : array of distances to coastline in meters. Negative distances are over ocean.
```python
>>> distance_to_shoreline([102.66, 116.81], [3.06, 0.72], path_to_GSHHG = "./useful_files/dist_to_GSHHG_v2.3.7_1m.nc")
[85280.41 88413.38]
```



## get_lat_lon_h_satellite(year, month, day, hour, minute, second, microsecond, satellite_name, printwarning, updateTLE, update_TLE_information_array = [], file_finals2000A_path = "")
Get latitude, longitude and altitude of a satellite. The function automatically download TLE files from spacetrack provided given the correct satellite ID. The function also checks if there are some discrepancies between the closest in time TLE measurements. This prints a warning if for example ISS has completed a maneuver between two TLEs leading to large errors in the calculated position. The functions does not decide which is the correct TLE in these cases.

| Satellite name      | Satellite ID |
| ----------- | ----------- |
| ISS      | 25544       |
| Fermi   | 33053        |
| AGILE   | 31135        |
| RHESSI   | 27370        |
| ...   | ...        |





* Parameters:
  * year : year as integer
  * month : month as integer
  * day : day as integer
  * hour : hour as integer
  * minute : minute as integer
  * second : second as integer
  * microsecond : microsecond as integer
  * satellite_name : satellite name, for example "ISS"
  * printwarning : True or False
  * updateTLE : True or False
    * TLE = TWO LINE ELEMENT file from SpaceTrack.
    * if True the function download and save a TLE file for the satellite of interest.
    ***This only needs to be done once for calculating satellite positions before the current date.***
  * update_TLE_information_array : a list where first element is satellite ID and the second element is a path to a txt file where you have written your SpaceTrack username and password
  * file_finals2000A_path : finals2000A.all file with Earth orientation data, available from https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html


* Returns:    
  * latitude : latitude of satellite in degrees
  * longitude : longitude of satellite in degrees
  * altitude : altitude of satellite in meters
  * radial distance : safety check. Radial distance [km] between the satellite positions calculated with the closest TLE and the TLE before the closest TLE.
  * radial distance : safety check. Radial distance [km] between the satellite positions calculated with the closest TLE and the TLE after the closest TLE

```python
>>> get_lat_lon_h_satellite(2018, 2, 16, 20, 57, 49, 123456, "ISS",
      printwarning=True, updateTLE=True,
      update_TLE_information_array = [25544, "./useful_files/spacetrackclient_usernamepw.txt"],
      file_finals2000A_path = "./useful_files/finals2000A.all.txt")

Trying to update TLE data
Current year: 2022
Current month: 5
Current day: 19
TLE data successfully updated.
(-48.085368276723145, 38.4916975322021, 423033.19557148905, 0.08541662669725056, 0.1604117499451488)
>>>
>>>
>>> get_lat_lon_h_satellite(2018, 8, 3, 16, 21, 43, 0, "AGILE",
      printwarning=True, updateTLE=True,
      update_TLE_information_array = [31135, "./useful_files/spacetrackclient_usernamepw.txt"],
      file_finals2000A_path = "./useful_files/finals2000A.all.txt")

Trying to update TLE data
Current year: 2022
Current month: 5
Current day: 19
TLE data successfully updated.
(-0.8903378881572883, 137.2636161735792, 471921.8884056839, 0.3568265953470963, 0.40153284482341883)
>>>
>>>
>>> get_lat_lon_h_satellite(2018, 8, 30, 6, 23, 16, 0, "AGILE",
      printwarning=True, updateTLE=False,
      update_TLE_information_array = [31135, "./useful_files/spacetrackclient_usernamepw.txt"],
      file_finals2000A_path = "./useful_files/finals2000A.all.txt")
(-2.340697039821966, 102.90778847628414, 472517.01810010563, 0.5807744896745226, 0.8704856077922032)
```





## check_if_folder_exist_or_make_it(directory)
Creates a new directory if the directory does not exist.

* Parameters:
  * directory : name of directory you want to create

* Returns:    
  * 0
```python
check_if_folder_exist_or_make_it("newdirectoryname")
```



## check_order_of_magnitude(number)
Returns the order of magnitude of the number

* Parameters:
  * number : number you want to check the order of magnitude of

* Returns:    
  * order of magnitude of number
```python
>>> check_order_of_magnitude(10)
1
>>> check_order_of_magnitude(100)
2
```
