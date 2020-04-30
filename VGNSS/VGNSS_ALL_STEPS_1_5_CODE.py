import sys
import os
import numpy as np

from numpy import cos,pi,sin,pi,arccos, tan, arctan, sqrt, arctan2
from numpy.linalg import inv
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

if verbose:
    print('GPS       PRN        Range             Azimuth        Elevation ')
    print('______________________________________________________________ ')
    for prn in p_az_elev_gps[visible_gps,:]:
        print('\t%5d %15.3f %15.3f %15.3f'%(prn[0],prn[1],prn[2],prn[3]))
    print('')
    print('')
    print('')     
## 4.6.1 Determine which GLONASS Satellites are Visible

visible_gln = p_az_elev_gln[:,3]>0

if verbose:
    print('GLONASS   PRN        Range             Azimuth        Elevation ')
    print('______________________________________________________________ ')
    for prn in p_az_elev_gln[visible_gln,:]:
        print('\t%5d %15.3f %15.3f %15.3f'%(prn[0],prn[1],prn[2],prn[3]))
    
## 4.7.0 Create the arrays t, x, and z
t = np.arange(0,2.01*pi,0.01*pi)
x = 1/2*pi*cos(t)
y = 1/2*pi*sin(t)

if not 'drawn' in locals():
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
    drawn = True

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

## 5.0.0 re-assign the satellite identified by the list bv if necesarry

if bv[0] == 1:
    if any(p_az_elev_gps[visible_gps,0]==bv[1]):
        if verbose:
            print( 'GPS vehicle '+str(bv[1])+ ' is visible')
    else:
        if verbose:
            print( 'GPS PRN '+str(bv[1])+ ' is not visible; re-assigning bv to GPS PRN '\
                   +str(int(p_az_elev_gps[visible_gps,0][0])))
        bv = [1,int(p_az_elev_gps[visible_gps,0][0])]
elif bv[0] == 2:
    if any(p_az_elev_gln[visible_gln,0]==bv[1]):
        if verbose:
            print( 'GLONASS vehicle '+str(bv[1])+ ' is visible')
    else:
        if verbose:
            print( 'GLONASS PRN '+str(bv[1])+ ' is not visible; re-assigning bv to GPS PRN '\
              +str(int( [visible_gln,0][0])))
        bv = [1,int(p_az_elev_gln[visible_gps,0][0])]
        
## 5.0.1 Calculate Total Uncertainty to be Applied to each GPS Satellite Range
# Assing the first n_gps random numbers to the column vector rnd_gps
rnd_gps=np.array(rnd[0:n_gps]).reshape(n_gps,1)

# Determine the number of visible GPS satellites
n_gps_vis = sum(visible_gps)

# Allocate memory
error_gps = np.zeros((n_gps_vis,7))

# Set the PRNs of the GPS satellites in the first column
error_gps[:,0:1]=p_az_elev_gps[visible_gps,0:1]

# Add the elevations in the 2nd column
error_gps[:,1:2]=p_az_elev_gps[visible_gps,3:4]

# Add the tropospheric error from the model defined in step 1.3.2 in the 3d column
if verbose:
    print('Tropo delay model: %4.2f/sin(sqrt(elev^2+%4.2f^2))'%(k,C*180/pi))
error_gps[:,2:3]=k/sin(sqrt(error_gps[:,1:2]**2+C**2))

# Add the random noise in the 4th column
error_gps[:,3:4]=rnd_gps[visible_gps]

# Add the clockbias (mapped to meters) in the 5th column
error_gps[:,4:5]=rx_dt*c

# If the blunder vehicle is part of the GPS constellation add the blunder to the vehicle identifies by bv
if bv[0]==1:
    error_gps[error_gps[:,0]==bv[1],5]=blunder
    
# Finally determine the total uncertainty
error_gps[:,-1]=np.sum(error_gps[:,2:],1)


if verbose:
    print('rnd_gps = ')
    for r in rnd_gps:
        print('\t%15.3f'%r)
    print('\nn_gps_vis = %d\n'%n_gps_vis)
    print('error_gps = ')
    print('\t  PRN    elev  Tr. Delay Rnd Noise   Clock Bias  Blunder   Tot. Uncertainty')
    for i in range(n_gps_vis):
        print('\t%5d %7.3f %10.3f %9.3f %12.3f %8.3f   %16.3f'% \
            (error_gps[i,0],error_gps[i,1],error_gps[i,2],error_gps[i,3], \
             error_gps[i,4],error_gps[i,5],error_gps[i,6]))
        
## 5.0.2 Calculate Total Uncertainty to be Applied to each GLONASS Satellite Range
# Assing the first n_gps random numbers to the column vector rnd_gps
rnd_gln=np.array(rnd[n_gps:]).reshape(n_gln,1)

# Determine the number of visible GPS satellites
n_gln_vis = sum(visible_gln)

# Allocate memory for matrix error_gln
error_gln = np.zeros((n_gln_vis,7))

# Set the PRNs of the GLONASS satellites in the first column
error_gln[:,0:1]=p_az_elev_gln[visible_gln,0:1]

# Add the elevations in the 2nd column
error_gln[:,1:2]=p_az_elev_gln[visible_gln,3:4]

# Add the tropospheric error from the model defined in step 1.3.2 in the 3d column
if verbose:
    print('Tropo delay model: %4.2f/sin(sqrt(elev^2+%4.2f^2))'%(k,C*180/pi))
error_gln[:,2:3]=k/sin(sqrt(error_gln[:,1:2]**2+C**2))

# Add the random noise in the 4th column
error_gln[:,3:4]=rnd_gln[visible_gln]

# Add the clockbias (mapped to meters) in the 5th column
error_gln[:,4:5]=rx_dt*c

# If the blunder vehicle is part of the GLONASS constellation add the blunder to the vehicle identifies by bv
if bv[0]==2:
    error_gln[error_gln[:,0]==bv[1],5]=blunder
    
# Finally determine the total uncertainty
error_gln[:,-1]=np.sum(error_gln[:,2:],1)


if verbose:
    print('rnd_gln = ')
    for r in rnd_gln:
        print('\t%15.3f'%r)
    print('\nn_gln_vis = %d\n'%n_gln_vis)
    print('error_gln = ')
    print('\t  PRN    elev  Tr. Delay Rnd Noise   Clock Bias  Blunder   Tot. Uncertainty')
    for i in range(n_gln_vis):
        print('\t%5d %7.3f %10.3f %9.3f %12.3f %8.3f   %16.3f'% \
            (error_gln[i,0],error_gln[i,1],error_gln[i,2],error_gln[i,3], \
             error_gln[i,4],error_gln[i,5],error_gln[i,6]))
##5.0.3 Calculate the GPS Pseudoranges

# Allocate memory
obs_gps = np.zeros((n_gps_vis,2))

obs_gps[:,0]=error_gps[:,0]
obs_gps[:,1:2]=rho_gps[visible_gps]+error_gps[:,6:7]

if verbose:
    for i in range(n_gps_vis):
        print( 'Observed pseudorange: %15.3f to GPS PRN %2d '%(obs_gps[i,1],obs_gps[i,0]))

##5.0.4 Calculate the GLONASS Pseudoranges

# Allocate memory
obs_gln = np.zeros((n_gln_vis,2))

obs_gln[:,0]=error_gln[:,0]
obs_gln[:,1:2]=rho_gln[visible_gln]+error_gln[:,6:7]

if verbose:
    for i in range(n_gln_vis):
        print( 'Observed pseudorange: %15.3f to GLONASS PRN %2d '%(obs_gln[i,1],obs_gln[i,0]))
    
## 5.0.5 GPS A-Priori Measurement Uncertainty

# Allocate memory
var_gps = np.zeros((n_gps_vis,1))
# Determine the variance of each of the observations
var_gps = error_gps[:,2:3]**2+error_gps[:,3:4]**2
# Determine the average variance and use this as the estimated a-priori msmt uncertainty
var_a_priori_gps = sum(var_gps)/n_gps_vis

if verbose:
    print( 'var_gps = ')
    for i in range(n_gps_vis):
        print('%12.3f'%var_gps[i])
        
    print('\nvar_a_priori_gps  = %12.3f'%var_a_priori_gps)
    
## 5.0.5 GPS A-Priori Measurement Uncertainty

# Allocate memory
var_gps = np.zeros((n_gps_vis,1))
# Determine the variance of each of the observations
var_gps = error_gps[:,2:3]**2+error_gps[:,3:4]**2
# Determine the average variance and use this as the estimated a-priori msmt uncertainty
var_a_priori_gps = sum(var_gps)/n_gps_vis

if verbose:
    print( 'var_gps = ')
    for i in range(n_gps_vis):
        print('%12.3f'%var_gps[i])
        
    print('\nvar_a_priori_gps  = %12.3f'%var_a_priori_gps)
   
## 5.0.6 GLONASS A-Priori Measurement Uncertainty

# Allocate memory
var_gln = np.zeros((n_gln_vis,1))
# Determine the variance of each of the observations
var_gln = error_gln[:,2:3]**2+error_gln[:,3:4]**2
# Determine the average variance and use this as the estimated a-priori msmt uncertainty
var_a_priori_gln = sum(var_gln)/n_gln_vis

if verbose:
    print( 'var_gln = ')
    for i in range(n_gln_vis):
        print('%12.3f'%var_gln[i])
        
    print('\nvar_a_priori_gln  = %12.3f'%var_a_priori_gln)

## 5.0.7 A-priori measurement uncertainty for the combined constellations

var_a_priori = (var_a_priori_gps*n_gps_vis+var_a_priori_gln*n_gln_vis)/(n_gps_vis+n_gln_vis)

if verbose:
    print('var_a_priori = '+str(var_a_priori))

## 5.0.8 Estimate of Initial Unknown Parameter Values

#Allocate memory
xyzt_0 = np.zeros((4,1))
xyzt_0[0:3] = xyz
xyzt_0 += np.array(bias).reshape(4,1)

if verbose:
    print('xyzt_0 = ')
    for i in range(4):
        print('%12.3f'%xyzt_0[i])
        
### 5.1.0 Creating the lists x_est, w and A
x_est = [xyzt_0] # The first (and currently last) estimate of receiver location x


### 5.1.1 Combining the Observations and Ephemeris

# The total number of satellite vehicles
n_svs = n_gps_vis+n_gln_vis

# Augment a vector of ones (GPS) with the GPS PRNs and Pseudoranges (note the use of the hstack() function)
l_obs = np.hstack((np.ones((n_gps_vis,1)),obs_gps))

# Augment l_obs with GLONASS data (note the use of both the hstack and vstack() functions)
l_obs = np.vstack((l_obs,np.hstack((2*np.ones((n_gln_vis,1)),obs_gln))))

# Similary create the array eph holding the constellation identifier, PRN and ephemeri
eph = np.hstack((np.ones((n_gps_vis,1)), eph_gps[visible_gps,:]))
eph = np.vstack((eph,np.hstack((2*np.ones((n_gln_vis,1)),eph_gln[visible_gln,:]))))

if verbose:
    print('n_svs = %d'%n_svs)
    
    print('\nl_obs = ')
    for pr in l_obs:
        print( '%7d %3d %15.3f'%(pr[0],pr[1],pr[2]))    
    
    print('\neph = ')
    for sv in eph:
        print( '%7d %3d %15.3f %15.3f %15.3f'%(sv[0],sv[1],sv[2],sv[3],sv[4]))
        
### 5.1.2 Determining the Misclosures and Design Matrix

def get_design_misclosure_range(x_est, eph, l_obs, c, verbose):
    
    # 5.1.2 Allocate memory for the misclosure vector w and design matrix A
    n_svs = eph.shape[0]
    
    w = np.zeros((n_svs,1))
    A = np.zeros((n_svs,4))
    
    ### 5.1.3 Create the model ranges rho_mod and observations l_mod
    
    rho_mod = sqrt((eph[:,2:3]-x_est[0])**2+(eph[:,3:4]-x_est[1])**2+(eph[:,4:5]-x_est[2])**2)   
    l_mod = rho_mod + c*x_est[3]  
    
    ### 5.1.4 Using the model observation and the pseudoranges determine the misclosures
    
    w = l_mod - l_obs[:,2:3]
    
    ### 5.1.5 Using the model Observations determine design natrix A

    A[:,0:1] = -(eph[:,2:3]-x_est[0])/rho_mod
    A[:,1:2] = -(eph[:,3:4]-x_est[1])/rho_mod
    A[:,2:3] = -(eph[:,4:5]-x_est[2])/rho_mod
    A[:,3:4] = c
    
    if verbose:
        print('l_mod = ')
        for l in l_mod:
            print('%15.3f'%l)
        print('rho_mod = ')
        for rho in rho_mod:
            print('%15.3f'%rho)


    return w, A
result = get_design_misclosure_range(x_est[-1],eph,l_obs, c, False)
w = [result[0]]
A = [result[1]]
