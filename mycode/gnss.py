import os
import urllib.request
from datetime import datetime, timezone, timedelta

import numpy as np

from mycode.sp3 import SP3

# from mycode.ephemeris import Ephemeris

## Created by Semme J. Dijkstra 12/24/2019
    #  

class GNSS:
    """A Class for handling GNSS Data"""
    
    def __init__(self, datapath):
        
        # Check to see whether the datapath exists
        
        if os.path.exists(datapath): 
            self.datapath = datapath
        else:  # raise a meaningful error if the data folder does not exist
            raise RuntimeError("GNSS: Unable to locate the data folder: " + data_folder)
            
        c=299792458       # Speed of light in vaccuum
        
        # Time management
        self.epochs = list()
        self.gnss_weeks=list()
        self.gnss_times=list()
        self.gnss_week_rollovers=list()
        
        # Ephemeris data file management
        self.eph_gps_sp3_filenames = list()
        self.eph_gln_sp3_filenames = list()
        self.eph_gps_nav_filenames = list()
        self.eph_gln_nav_filenames = list()
        self.eph = list()
        
        # Ephemeris data structures 
        
        self.gps_sp3 = list()
        self.gln_sp3 = list()
        self.gps_nav = list()
        self.gln_nav = list()
        
        
        
    
    def add_epoch(self, year, month, day, hour, minute, second):
        # Add any epoch - not necesarily in order
        
        # Make sure to insert the data in the proper location!!
    
        error( 'GNSS.add_epoch; not yet implemented')
        
    def add_next_epoch_ephemeris(self, year, month, day, hour, minute, second, micro_second):
        # This is a simple version that will only allow you to add epochs that are 
        # successive. Using this prevent you from having to search for already read 
        # ephemeris data
        
        t = datetime(year,month,day,hour,minute,second,micro_second)
        epoch = t - datetime(1980,1,6,0,0,0)
        
        if (not len( self.epochs)) or (epoch > self.epochs[-1]):
            # Add the time
            self.epochs.append( epoch)
        else:
             print( 'Current epoch occured before or on last epoch - use add_epoch instead') 
             return # do nothing
        
        # Determine the gnns week and time for this epoch

        self.gnss_weeks.append(np.floor(epoch.days / 7))
        self.gnss_times.append((epoch.days - self.gnss_weeks[-1] * 7)*24*3600 + epoch.seconds)
        self.gnss_week_rollovers.append(np.fix(self.gnss_weeks[-1]/1024))
        
        # Deal with GNSS first
        
        # Form the ephemeris file name for the current epoch in GPS time, best quality by default
        
        eph_gps_filename='igs'
        eph_gps_filename += '%04d'%self.gnss_weeks[-1] 
        eph_gps_filename += '%1d'%np.fix(self.gnss_times[-1]/24/60/60)
        eph_gps_filename += '.sp3'
        
        # Check to see whether the file has already been processed
        
        if not eph_gps_filename in self.eph_gps_sp3_filenames:
            
            datapath = os.path.join(self.datapath, eph_gps_filename)
            
            # Only download the data if we do not already have it
            
            if not os.path.exists(datapath):
                print( "Downloading ephemeris file: " + eph_gps_filename)
                # Form the URL to download the file from

                eph_url = 'ftp://cddis.gsfc.nasa.gov/gnss/products/' 
                eph_url += '%04d'%self.gnss_weeks[-1] + '/' + eph_gps_filename + '.Z';

                # Download the data file
        
                weburl=urllib.request.urlretrieve(eph_url, \
                                                  os.path.join(self.datapath, eph_gps_filename + '.Z'))
        
                # Unfortunately Python does not have a Unix uncompress available, which is 
                # needed to decompress a .Z file on the ePom server. We will use a sidestep to a 
                # system call assuming that gunzip is available (true for the ePOM server - 
                # but really should test)
        
                os.system('gunzip -f ' + os.path.join(self.datapath, eph_gps_filename + '.Z'))
            
            # Load the data into an sp3 object and add it to the list
            
            gps_sp3 = SP3( os.path.join(self.datapath, eph_gps_filename),'gps')
            self.gps_sp3.append( gps_sp3)
            self.eph_gps_sp3_filenames.append( eph_gps_filename) 
       
        # The Glonass filename will be almost the same
        eph_gln_filename = eph_gps_filename[:2]+'l'+eph_gps_filename[3:]   
        
        if not eph_gln_filename in self.eph_gln_sp3_filenames:
            datapath = os.path.join(self.datapath, eph_gln_filename)
            
            # Only download the data if we do not already have it
            
            if not os.path.exists(datapath):
                print( "Downloading ephemeris file: " + eph_gln_filename)
                # Form the URL to download the file from

                eph_url = 'ftp://cddis.gsfc.nasa.gov/glonass/products/' 
                eph_url += '%04d'%self.gnss_weeks[-1] + '/' + eph_gln_filename + '.Z';

                # Download the data file
        
                weburl=urllib.request.urlretrieve(eph_url, \
                                                  os.path.join(self.datapath, eph_gln_filename + '.Z'))
        
                # Unfortunately Python does not have a Unix uncompress available, which is 
                # needed to decompress a .Z file on the ePom server. We will use a sidestep to a 
                # system call assuming that gunzip is available (true for the ePOM server - 
                # but really should test)
        
                os.system('gunzip -f ' + os.path.join(self.datapath, eph_gln_filename + '.Z'))
            
            # Load the data into an sp3 object and add it to the list
            
            glonass_sp3 = SP3( os.path.join(self.datapath, eph_gln_filename),'glonass')
            self.gln_sp3.append( glonass_sp3)
            self.eph_gln_sp3_filenames.append( eph_gln_filename)
                
        
        

    def get_ephemeris(self, epoch):
        if len( self.epochs) != len( self.ephemeris):
            get_ephemeris = get_ephemeris( epoch)
            
            
        return get_ephemeris
            
            
                   
              
