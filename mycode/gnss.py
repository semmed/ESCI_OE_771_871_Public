import os
from ftplib import FTP
import urllib
from ftplib import FTP_TLS
from datetime import datetime, timezone, timedelta
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial
from numpy import sqrt, sin, cos, arctan2
import numpy as np
import scipy as sp
from mycode.sp3 import SP3
from mycode.rinex_nav import RINEX_nav

# from mycode.ephemeris import Ephemeris

## Created by Semme J. Dijkstra 12/24/2019
# Updated 2/10/2022 FTP Server interfacing has changed (once again...)

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
        
        # Ephemeris data structures 
        self.gps_sp3 = list()
        self.gln_sp3 = list()
        self.gps_nav = list()
        self.gln_nav = list()
        
        # Associated ephemeris data
        self.eph_gps_nav = list()
        self.eph_gln_nav = list()
        self.eph_gps_sp3 = list()
        self.eph_gln_sp3 = list()
        
    
    def add_epoch(self, year, month, day, hour, minute, second):
        # Add any epoch - not necesarily in order
        # Make sure to insert the data in the proper location!!
    
        error( 'GNSS.add_epoch; not yet implemented')
        
    def add_next_epoch_ephemeris(self, year, month, day, hour, minute, second, micro_second, email='student@ccom.unh.edu'):

        # This is a simple version that will only allow you to add single epochs that are 
        # successive. Using this prevents you from having to search for already read 
        # ephemeris data, but it is not very efficient for adding multiple epochs
        
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
        
        # Check to see whether the final ephemeris file has already been processed
        


        if not eph_gps_filename in self.eph_gps_sp3_filenames:
            # Create a data path to the GPS file
            datapath = os.path.join(self.datapath, eph_gps_filename)
            
            # Only download the data if we do not already have it
            if not os.path.exists(datapath):
                print( 'Attempting to download precise ephemerides file: ' + eph_gps_filename)
                ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
                ftps.login(user='anonymous', passwd=email)
                ftps.prot_p()
                ftps.cwd('pub/gps/products/%04d'%(self.gnss_weeks[-1]))
                ftps.retrbinary("RETR " + eph_gps_filename+'.Z', open(os.path.join(datapath + '.Z'), 'wb').write)
                os.system('gunzip -f ' + os.path.join(datapath + '.Z'))

            # Create an SP3 object and add it to the list
            gps_sp3 = SP3( os.path.join(self.datapath, eph_gps_filename),'gps')
            self.gps_sp3.append( gps_sp3)
            self.eph_gps_sp3_filenames.append( eph_gps_filename)

            # Now that we are guaranteed to have the data calculate the ephemeris for the epoch
            self.eph_gps_sp3.append(self.get_single_epoch_ephemeris_from_sp3(self.epochs[-1], self.gps_sp3[-1]))



        # The Glonass filename will be almost the same
        eph_gln_filename = eph_gps_filename[:2]+'l'+eph_gps_filename[3:]  
        
        if not eph_gln_filename in self.eph_gln_sp3_filenames:
            datapath = os.path.join(self.datapath, eph_gln_filename)
            
            # Only download the data if we do not already have it
            if not os.path.exists(datapath):
                print( 'Attempting to download precise ephemerides file: ' + eph_gln_filename)
                ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
                ftps.login(user='anonymous', passwd=email)
                ftps.prot_p()
                ftps.cwd('pub/glonass/products/%04d'%(self.gnss_weeks[-1]))
                ftps.retrbinary("RETR " + eph_gln_filename+'.Z', open(os.path.join(datapath + '.Z'), 'wb').write)
                os.system('gunzip -f ' + os.path.join(datapath + '.Z'))

            # Create an SP3 object and add it to the list
            glonass_sp3 = SP3( os.path.join(self.datapath, eph_gln_filename),'glonass')
            self.gln_sp3.append( glonass_sp3)
            self.eph_gln_sp3_filenames.append( eph_gln_filename)

            # Now that we are guaranteed to have the data calculate the ephemeris for the epoch
            self.eph_gln_sp3.append(self.get_single_epoch_ephemeris_from_sp3(self.epochs[-1], self.gln_sp3[-1]))
                
       # We also will need the broadcast ephemeris files - we will only read the GPS broadcast data for now
       # (the other systems have different data formats - too much work for this assignment)
      
    
       # First the GPS nav data files
            
        dt = t - datetime(year,1,1,0,0,0)
        brdc_filename='brdc'
        brdc_filename += '%03d0.'%(dt.days + 1) 
        brdc_filename += t.strftime('%y')
        brdc_filename += 'n'
            
        datapath = os.path.join(self.datapath, brdc_filename)
            
        # Only attempt to load if it is not already loaded
           
        if not brdc_filename in self.eph_gps_nav_filenames:
            
            # If the data file does not already exist try to download it
            
            if not os.path.exists(datapath):
                print( "Downloading ephemeris file: " + brdc_filename)

                ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
                ftps.login(user='anonymous', passwd=email)
                ftps.prot_p()
                print('pub/gps/data/daily/' + t.strftime('%Y') + '/brdc/')
                ftps.cwd('pub/gps/data/daily/' + t.strftime('%Y') + '/brdc/')
                try:
                    ftps.retrbinary("RETR " + brdc_filename+'.gz', open(os.path.join(datapath + '.gz'), 'wb').write)
                    os.system('gunzip -f ' + os.path.join(self.datapath, brdc_filename + '.gz'))
                except:
                    ftps.retrbinary("RETR " + brdc_filename+'.Z', open(os.path.join(datapath + '.Z'), 'wb').write)
                    os.system('gunzip -f ' + os.path.join(self.datapath, brdc_filename + '.Z'))
                    


                # Also unzip this data file
           
                
                    
            
            gps_nav = RINEX_nav( os.path.join(self.datapath, brdc_filename),'gps')
            self.gps_nav.append( gps_nav)
            self.eph_gps_nav_filenames.append( brdc_filename) 
            
            # Now that we are guaranteed to have the data calculate the ephemeris for the epoch
            
            self.eph_gps_nav.append(self.get_single_epoch_ephemeris_from_gps_nav(self.gnss_times[-1], self.gps_nav[-1]))

        # We are all done

        if 'ftp' in locals():
            ftp.quit()

    # Return the ephemeris of the last epoch in the object  
    # This method is for the VGNSS assignment only as it deals poorly with epochs close to the borders of the span
    # covered by the sp3 files.
    # This requires significant extra coding that is not helpful to the assignment
    
    def get_single_epoch_ephemeris_from_sp3(self, epoch, gps_sp3):
       
        # Chose an order of 7 for the Lagrange Polynomials, this should provide cm level accuracy
        poly_fit_order = 5

        # Determine the time differences between the sp3 starting epoch and the other sp3 epochs
        
        dt = np.zeros(len(gps_sp3.epoch_times))
        
        i = 0
        for t in gps_sp3.epoch_times:
            dt[i] = (t - gps_sp3.epoch_times[0]).seconds
            i += 1
            
        # Determine the time difference between the current epoch and the sp3 starting epoch
        
        t = (epoch - (gps_sp3.epoch_times[0] - datetime(1980,1,6,0,0,0))).seconds
        
        # Find the start index of the epochs to use for the interpolation 
        
        index = sum(t>dt)-poly_fit_order//2
        
        # Make sure that the polynomial fit order does not exceed the number of available epochs
        # This is where extra effort is needed to deal with epochs close to the start end end epochs contained within
        # the sp3 file

        while poly_fit_order+1 > len(gps_sp3.epoch_times):
            print( 'Order of Lagrange polynomial fitting exceeds number of samples!')
            print( 'Reducing the order from: '+str(poly_fit_order)+' to: ' + str(len(gps_sp3.epoch_times)))
            poly_fit_order = len(gps_sp3.epoch_times)-1
        
        # Prevent negative indexes and index over-runs
        if index < 0:
            index = 0
        elif index + poly_fit_order+1 >= len(gps_sp3.epoch_times):
            index = len(gps_sp3.epoch_times) - poly_fit_order-1 
        
        # Create a 3D array holding the epoch data in the first dimension, the sv data in the 2nd dimension
        # and the individual data fields in the 3d dimension
        # So that all the data associated to sv n for all the epochs of interest is contained within sp_3[:,n,:] 
        eph_epochs = np.asarray(gps_sp3.epoch_pos_t_data[index:index+poly_fit_order+1])
        
        # Example, for all the epochs in the interpolation window the satellite PRN for the 2nd vehicle is given by
        # sp_3[:,1,1])
        
        eph=np.zeros([gps_sp3.nr_sats,4])
 
        # Produce ephemeris for all the satellite vehicles
        for i in range(gps_sp3.nr_sats):
            # Determine the Lagrange polynomials for the ECEF X, Y and Z coordinates
            eph[i,0] = eph_epochs[0,i,1]
            poly = lagrange(dt[index:index+poly_fit_order+1],eph_epochs[:,i,2])
            p = Polynomial(poly).coef
            eph[i,1] = np.polyval(p,t)
            poly = lagrange(dt[index:index+poly_fit_order+1],eph_epochs[:,i,3])
            p = Polynomial(poly).coef
            eph[i,2] = np.polyval(p,t)
            poly = lagrange(dt[index:index+poly_fit_order+1],eph_epochs[:,i,4])
            p = Polynomial(poly).coef
            eph[i,3] = np.polyval(p,t)

        return eph

            
    def get_single_epoch_ephemeris_from_gps_nav(self, gps_time, gps_nav):
        
        # This method goes along with the paper 'Computing satellite velocity using the broadcast ephemeris'
        # By Remondi, Benjamin W - Ben Remondi is one of the authorative figures in the GNSS literature of whom
        # You should be aware.

        # Calculate ephemeris using Keplerian model analysis
        
        ### 3.1.0 Establishment of the WGS84 System Parameters
        
        GM = 3.986005*10**14      # m3/s2 WGS-84 value for the product of gravitational constant G and the mass of the Earth M
        w_e = 7.2921151467*10**-5 # rad/s WGS-84 value of the Earth’s rotation rate
        

        # 3.1.1 Allocation of Memory for Ephemeris Data
        eph=np.zeros([len( gps_nav.sat_prns),4])
        
        # There is no ephemeris data for every satellite in every record - loop through the satellites and then the 
        # epochs - at least the gps_time of interest may predate the t_oe, which avoids all kinds of trouble
        row = 0
        for prn in gps_nav.sat_prns:
            
            # The subset of records for this PRN
            records = gps_nav.records[ gps_nav.records['prn'] == prn]
            
            # Calculate the index of the time of epoch t_oe
            i_oc = 0    # Ensures that the index is -1 if there is no data available for the epoch of interest 
            i_oc = sum( records['t_oe'] <= gps_time) - 1          
                
            # Time since reference epoch (t_oe). Where gps_time is the time of measurement by the receiver 
            t_k = (gps_time - records[i_oc]['t_oe'])
    
            ## Keplerian orbit parameter modeling
            # Satellite orbital period
            n_0 = sqrt(GM / records[i_oc]['a_rt']**6)       # Computed mean motion,also known as greek mu
            n = n_0 + records[i_oc]['D_n']
            M_k = records[i_oc]['M_0'] + n * t_k
            M_k_dot = n
            E_k = M_k
            for i in range(7):                              # iterate 7 times just to be sure (2 should be enough)
                E_k = M_k + records[i_oc]['e']*sin(E_k)
                
            E_k_dot = M_k_dot/(1.0 - records[i_oc]['e']*cos(E_k))
            
            # In the line, below, ta_k is the true anomaly (which is nu in the ICD-200).
            ta_k = arctan2( sqrt(1.0-records[i_oc]['e']**2)*sin(E_k), cos(E_k)-records[i_oc]['e'])

            phi_k = ta_k + records[i_oc]['w']
            
            corr_u = records[i_oc]['c_us']*sin(2.0*phi_k) + records[i_oc]['c_uc']*cos(2.0*phi_k)
            corr_r = records[i_oc]['c_rs']*sin(2.0*phi_k) + records[i_oc]['c_rc']*cos(2.0*phi_k)
            corr_i = records[i_oc]['c_is']*sin(2.0*phi_k) + records[i_oc]['c_ic']*cos(2.0*phi_k)
            
            u_k = phi_k + corr_u
            r_k = (records[i_oc]['a_rt']**2)*(1.0-records[i_oc]['e']*cos(E_k)) + corr_r
            i_k = records[i_oc]['i_0'] + records[i_oc]['i_dot']*t_k + corr_i
            
            # Polar coordinates of satellite vehicle prn

            xp_k = r_k*cos(u_k)
            yp_k = r_k*sin(u_k)
            
            W_k = records[i_oc]['W_0'] + (records[i_oc]['W_dot']-w_e)*t_k - w_e*records[i_oc]['t_oe']
            
            # Coordinate transformation to ECEF coordinates.
            
            eph[row, 0] = prn
            eph[row, 1] = xp_k*cos(W_k) - yp_k*sin(W_k)*cos(i_k);
            eph[row, 2] = xp_k*sin(W_k) + yp_k*cos(W_k)*cos(i_k);
            eph[row, 3] =                 yp_k*sin(i_k);
            row += 1         


        return eph

                   
              
