import sys
import os
import numpy as np

from numpy import cos,pi,sin,pi,arccos, tan, arctan, sqrt, arctan2
from datetime import datetime, timedelta
from pathlib import Path
import matplotlib.pyplot as plt

vgnss_path=Path('../') # Get the path to the folder containing the mycode and mydata folders

print(vgnss_path.resolve())
sys.path.append(str(vgnss_path.resolve())) # add the folder to the list of paths 

from mycode.gnss import GNSS
from mycode.sp3 import SP3

verbose = False
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
eph_gps=gnss.eph_gps_sp3[-1]
eph_gln=gnss.eph_gln_sp3[-1]

# 4.1.0 The parameters for WGS-84
a = 6378137.0000
b = 6356752.3142

## Calculate f, e2 and N
f = (a-b)/a
e2=1-b**2/a**2

## 4.1.1 Prime radius of curvature 
N = a/sqrt(1-e2*sin(lat)**2)

## 4.1.2 To create a 3x1 array use xyz = np.array([[1],[2],[3]])
xyz = np.array([[(N+height)*cos(lat)*cos(lon)], \
              [(N+height)*cos(lat)*sin(lon)], \
              [(N*b**2/a**2+height)*sin(lat)]]  )

if verbose:
    print('xyz = ')
    for i in range(len(xyz)):
        print('\t%15.3f'%(xyz[i]))

## 4.2.0 Calculate the Receiver to GPS SV baselines
Dxyz_gps = eph_gps[:,1:4].T-xyz
n_gps = Dxyz_gps.shape[ 1]

if verbose:
    print('n_gps: %d\n'%n_gps)
    print('Dxyz_gps = ')
    for i in range(n_gps):
        print('\t%15.3f %15.3f %15.3f'%(Dxyz_gps[0,i],Dxyz_gps[1,i],Dxyz_gps[2,i]))

## 4.2.1 Calculate the Receiver to GLONASS SV baselines

Dxyz_gln = eph_gln[:,1:4].T-xyz
n_gln = Dxyz_gln.shape[1]

if verbose:
    print('n_gln: %d\n'%n_gln)
    print('Dxyz_gln = ')
    for i in range(n_gln):
        print('\t%15.3f %15.3f %15.3f'%(Dxyz_gln[0,i],Dxyz_gln[1,i],Dxyz_gln[2,i]))
        
## 4.3.0 Calculate the range from the receiver to each GPS satellite
# First allocate memory by creating a n_gpsx1 2D array
rho_gps=np.zeros((n_gps,1))

for i in range(n_gps):
    rho_gps[i] = sqrt(Dxyz_gps[:,i].T@Dxyz_gps[:,i])
    
if verbose:
    print('rho_gps = ')
    for i in range(n_gps):
        print('\t%15.3f'%(rho_gps[i]))
        
## 4.3.1 Calculate the range from the receiver to each GPS satellite
# First allocate memory by creating a n_gpsx1 2D array
rho_gln=np.zeros((n_gln,1))

for i in range(n_gln):
    rho_gln[i] = sqrt(Dxyz_gln[:,i].T@Dxyz_gln[:,i])
    
if verbose:
    print('rho_gln = ')
    for i in range(n_gln):
        print('\t%15.3f'%(rho_gln[i]))
        
## 4.4.0.0 Calculate Reflection Matrix P2
P2 = np.eye(3)
P2[1,1]=-P2[1,1]

if verbose:
    print('P2:')
    print(P2)
    
## 4.4.0.1 Calculate Reflection Matrix R2
R2 = np.asarray([[cos(lat-pi/2), 0, -sin(lat-pi/2)],\
                 [0,1,0],\
                 [sin(lat-pi/2),0,cos(lat-pi/2)]])

if verbose:
    print('R2:')
    print(R2)
    
#### 4.4.0.2 Calculate Reflection Matrix R3
R3 = np.asarray([[cos(lon-pi),sin(lon-pi),0],\
                 [-sin(lon-pi),cos(lon-pi),0],\
                 [0,0,1]])

if verbose:
    print('R3:')
    print(R3)
    
#### 4.4.0.3 Calculate Reflection Matrix Q

Q = P2@R2@R3

if verbose:
    print('Q:')
    print(Q)
    
## 4.4.0.4 Calculate Topocentric Cartesian (North East Up) coordinates for GPS satellites 
# Allocate memory by using the shape of Dxyz_gps (Dxyz_gps.shape)
NEU_gps = np.zeros((Dxyz_gps.shape))

# Cycle through all GPS satellites
for i in range(n_gps):
    NEU_gps[:,i:i+1] = Q@Dxyz_gps[:,i:i+1]
    
if verbose:
    print('NEU_gps = ')
    for i in range(n_gps):
        print('\t%15.3f %15.3f %15.3f'%(NEU_gps[0,i],NEU_gps[1,i],NEU_gps[2,i]))
        
## 4.4.0.5 Calculate Topocentric Cartesian (North East Up) coordinates for 
# GLONASS satellites
# Allocate memory by using the shape of Dxyz_gps (Dxyz_gps.shape)
NEU_gln = np.zeros((Dxyz_gln.shape))

# Cycle through all GLONASS satellites
for i in range(n_gln):
    NEU_gln[:,i:i+1] = Q@Dxyz_gln[:,i:i+1]
    
if verbose:
    print('NEU_gln = ')
    for i in range(n_gln):
        print('\t%15.3f %15.3f %15.3f'%(NEU_gln[0,i],NEU_gln[1,i],NEU_gln[2,i]))
        
## 4.5.0 Determination of the Azimuth from the Receiver to the GPS satellites
# Allocate memory 
p_az_elev_gps = np.zeros((n_gps,4))

# Cycle through all GPS satellites
for i in range(n_gps):
    p_az_elev_gps[i,0] = eph_gps[i,0]
    p_az_elev_gps[i,1] = sqrt(NEU_gps[1,i]**2+NEU_gps[0,i]**2)
    p_az_elev_gps[i,2] = arctan2(NEU_gps[1,i],NEU_gps[0,i])
    p_az_elev_gps[i,3] = arctan(NEU_gps[2,i]/sqrt(NEU_gps[1,i]**2+NEU_gps[0,i]**2))
    
if verbose:
    print('p_az_elev_gps = ')
    for i in range(n_gps):
        print('\t%3d %15.3f %15.3f %15.3f'%(p_az_elev_gps[i,0],p_az_elev_gps[i,1],p_az_elev_gps[i,2],p_az_elev_gps[i,3]))
        
## 4.5.1 Determination of the Azimuth from the Receiver to the GLONASS satellites
# Allocate memory 
p_az_elev_gln = np.zeros((n_gln,4))

# Cycle through all GLONASS satellites
for i in range(n_gln):
    p_az_elev_gln[i,0] = eph_gln[i,0]
    p_az_elev_gln[i,1] = sqrt(NEU_gln[1,i]**2+NEU_gln[0,i]**2)
    p_az_elev_gln[i,2] = arctan2(NEU_gln[1,i],NEU_gln[0,i])
    p_az_elev_gln[i,3] = arctan(NEU_gln[2,i]/sqrt(NEU_gln[1,i]**2+NEU_gln[0,i]**2))
    
if verbose:
    print('p_az_elev_gln = ')
    for i in range(n_gln):
        print('\t%3d %15.3f %15.3f %15.3f'% (p_az_elev_gln[i,0],p_az_elev_gln[i,1],p_az_elev_gln[i,2],p_az_elev_gln[i,3]))
        
## 4.5.1 Determination of the Azimuth from the Receiver to the GLONASS satellites
# Allocate memory 
p_az_elev_gln = np.zeros((n_gln,4))

# Cycle through all GLONASS satellites
for i in range(n_gln):
    p_az_elev_gln[i,0] = eph_gln[i,0]
    p_az_elev_gln[i,1] = sqrt(NEU_gln[1,i]**2+NEU_gln[0,i]**2)
    p_az_elev_gln[i,2] = arctan2(NEU_gln[1,i],NEU_gln[0,i])
    p_az_elev_gln[i,3] = arctan(NEU_gln[2,i]/sqrt(NEU_gln[1,i]**2+NEU_gln[0,i]**2))
    
if verbose:
    print('p_az_elev_gln = ')
    for i in range(n_gln):
        print('\t%3d %15.3f %15.3f %15.3f'%(p_az_elev_gln[i,0],p_az_elev_gln[i,1],p_az_elev_gln[i,2],p_az_elev_gln[i,3]))

## 4.6.0 Determine which GPS Satellites are Visible

visible_gps = p_az_elev_gps[:,3]>0

print('GPS       PRN        Range             Azimuth        Elevation ')
print('______________________________________________________________ ')
for prn in p_az_elev_gps[visible_gps,:]:
    print('\t%5d %15.3f %15.3f %15.3f'%(prn[0],prn[1],prn[2],prn[3]))
print('')
print('')
print('')     
## 4.6.1 Determine which GLONASS Satellites are Visible

visible_gln = p_az_elev_gln[:,3]>0

print('GLONASS   PRN        Range             Azimuth        Elevation ')
print('______________________________________________________________ ')
for prn in p_az_elev_gln[visible_gln,:]:
    print('\t%5d %15.3f %15.3f %15.3f'%(prn[0],prn[1],prn[2],prn[3]))
    
## 4.7.0 Create the arrays t, x, and z
t = np.arange(0,2.01*pi,0.01*pi)
x = 1/2*pi*cos(t)
y = 1/2*pi*sin(t)

# 4.7.0 Plot the Circle 
plt.figure(figsize=(10, 10))
plt.plot(x,y,'k',linewidth=2)
# 4.7.1 Plot the zenith Dot
plt.plot(0,0,'.k',linewidth=3)
# 4.7.2 Make sure that the scaling for the axes is the same and the axes are invisible
plt.axis('equal')
plt.axis('off')
# 4.7.3 Add the Suptitle 'Sky Plot' in font size 24
plt.suptitle('Sky Plot', fontsize = 24)
# 4.7.4 Create the position string
str_pos = r'$\phi$: '
pos = lat*180/pi
str_pos += str(int(pos))
pos = (pos-int(pos))*60
str_pos += '$^{\circ}$'+str(int(pos))
pos = (pos-int(pos))*60
str_pos += '\'%.2f'%pos+'\"'
str_pos += '  $\lambda$: '
pos = lon*180/pi
str_pos += str(int(pos))
pos = (pos-int(pos))*60
str_pos += '$^{\circ}$'+str(int(pos))
pos = (pos-int(pos))*60
str_pos += '\'%.2f'%pos+'\"'
# 4.7.5 Add epoch string as the title
plt.title('For pos: ' + str_pos + ' WGS-84 at epoch: ' + epoch.strftime('%c') + ' UTC')

# 4.7.6 Transform the polar coordinates to plot coordiniates for the GPS satellites
x_plot = (pi/2-p_az_elev_gps[visible_gps,3:4])*sin(p_az_elev_gps[visible_gps,2:3])
y_plot = (pi/2-p_az_elev_gps[visible_gps,3:4])*cos(p_az_elev_gps[visible_gps,2:3])
# 4.7.7 Plot the visible GPS satellites, make sure that we may add a legend
plt.plot(x_plot,y_plot,'r.',label=u'GPS')
# 4.7.8 Label the visible GPS satellites with their PRN Codes
for i, prn in enumerate(p_az_elev_gps[visible_gps,:]):
    plt.text(x_plot[i]+.01+.01,y_plot[i]+.01,str(prn[0]), fontsize = 12);
# 4.7.9 Transform the polar coordinates to plot coordiniates for the GLONASS satellites
x_plot = (pi/2-p_az_elev_gln[visible_gln,3:4])*sin(p_az_elev_gln[visible_gln,2:3])
y_plot = (pi/2-p_az_elev_gln[visible_gln,3:4])*cos(p_az_elev_gln[visible_gln,2:3])
# 4.7.10 Plot the visible GLONASS satellites, make sure that we may add a legend
plt.plot(x_plot,y_plot,'b.',label=u'GLONASS')
# 4.7.11 Label the visible GLONASS satellites with their PRN Codes
for i, prn in enumerate(p_az_elev_gln[visible_gln,:]):
    plt.text(x_plot[i]+.01,y_plot[i]+.01,str(prn[0]), fontsize = 12);

    
plt.legend(ncol=2, loc='upper right')    
plt.show();

       