import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])



list = ["numpy", "matplotlib", "pandas", "progress", "astropy", "sgp4", "scipy", "pyproj", "netCDF4", "spacetrack"]


for package in list:
    install(package)
