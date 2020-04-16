import sys
import os
import numpy as np

from numpy import cos,pi,sin,pi,arccos, tan, arctan
from datetime import datetime, timedelta
from pathlib import Path

vgnss_path=Path('../') # Get the path to the folder containing the mycode and mydata folders

print(vgnss_path.resolve())
sys.path.append(str(vgnss_path.resolve())) # add the folder to the list of paths 

from mycode.gnss import GNSS
from mycode.sp3 import SP3

print('TAADAA')
# Step 2
lat=30.25*pi/180
lon=(-89-37/60-12/3600)*pi/180
height=10
epoch = datetime(2005,8,29,15,0,0)
# Worked example (bias equivalent to ~1000m)
c = 299792458
rx_dt = 1000/c
mydata_path = os.path.join(vgnss_path,"mydata")
sys.path.append(mydata_path)
rnd_path = os.path.join(mydata_path,"random_worked_example.txt")
print(rnd_path)
rnd = list()
rnd_file = open(rnd_path)
data = rnd_file.read()
rnd_file.close()
data = data.splitlines()
rnd = [float(i) for i in data]
k = 2.31
C = 2.5*pi/180
blunder = 2000
bv = [1, 2]    # GPS (1) vehicle 2, GLONASS (2) vehicle 23 would be [2 23]
bias = [1000,1000,1000,60]
# Step 3
gnss=GNSS( mydata_path)
gnss.add_next_epoch_ephemeris(epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, epoch.second, epoch.microsecond)
# Use the ephemeris from the SP3 file interpolation
eph_gps=gnss.eph_gps_sp3[-1]
eph_gln=gnss.eph_gln_sp3[-1]
