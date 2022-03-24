import os
import urllib.request
from datetime import datetime, timezone, timedelta
import numpy as np

    ## semme J. Dijkstra  2/10/2020 Created rinex_nav class with rinex_nav_a reader

# Create a data type for the epoch data in the RINEX nav files in a structured array
# Note that the fields are ordered by the sequence in which they are defined in the RINEX 2.11 format
# Explanation as per Seeber(2003) - here the data type for GPS boroadcast ephemeris


# d_type = 'prn', np.uint8), \     # SV PRN code
#           't_oc'      # Epoch - Time of clock (s)
#           'a_0'       # Clock correction bias coefficient
#           'a_1'       # Clock correction drift coefficient
#           'a_2'       # Clock correction drift rate coefficient
#           'iode'      # Issue of data ephemeris (s)
#           'c_rs'      # Amplitude of the sine harmonic correction term to the orbit radius (m)
#           'D_n'       # Delta n: Mean motion difference from computed value (rad/s)
#           'M_0'       # Mean anomaly at reference time (rad)
#           'c_uc'      # Amplitude of the cosine harmonic correction term to the argument of latitude(rad)
#           'e'         # Eccentricity (dimensionless)
#           'c_us'      # Amplitude of the sine harmonic correction term to the argument of latitude (rad)
#           'a_rt'      # Square root of the semi-major axis (m**1/2)
#           't_oe'      # Reference time, ephemeris parameters (s)
#           'c_ic'      # Amplitude of the cosine harmonic correction term to the angle of inclination(rad)
#           'W_0'       # Longitude of ascending node at reference time (rad)
#           'c_is'      # Amplitude of the sine harmonic correction term to the angle of inclination (rad)
#           'i_0'       # Inclination angle at reference time (rad)
#           'c_rc'      # Amplitude of the cosine harmonic correction term to the orbit radius (m)
#           'w'         # Argument of perigee (rad)
#           'W_dot'     # Rate of change of right ascension (rad/s)
#           'i_dot'     # Rate of change of inclination (rad/s)
#           'L2_code'   # Codes on L2 Channel
#           'toe_week'  # GPS Week of toe - not modulus(1024)!
#           'L2_flags'  # L2 P data flags
#           'sv_acc'    # SV accuracy (RMS?) in (m)
#           'sv_health' # Rate of change of right ascension (rad/s)
#           't_gd'      # L2 group phase delay (s)
#           'iode'      # Issue of data Clock (s)

gps_d_type = [('prn', np.uint8), \
          ('t_oc', datetime), \
          ('a_0', np.float64), \
          ('a_1', np.float64), \
          ('a_2', np.float64), \
          ('iode', np.float64), \
          ('c_rs', np.float64), \
          ('D_n', np.float64), \
          ('M_0', np.float64), \
          ('c_uc', np.float64), \
          ('e', np.float64), \
          ('c_us', np.float64), \
          ('a_rt', np.float64), \
          ('t_oe', np.float64), \
          ('c_ic', np.float64), \
          ('W_0', np.float64), \
          ('c_is', np.float64), \
          ('i_0', np.float64), \
          ('c_rc', np.float64), \
          ('w', np.float64), \
          ('W_dot', np.float64), \
          ('i_dot', np.float64), \
          ('L2_code', np.float64), \
          ('toe_week', np.float64), \
          ('L2_flags', np.float64), \
          ('sv_acc', np.float64), \
          ('sv_health', np.float64), \
          ('t_gd', np.float64), \
          ('iodc', np.float64)]

class RINEX_nav:
    """A Class for handling RINEX Nav Ephemeris Data"""
    
    # An object of this class will hold the data of one RINEX nav file
    
    def __init__(self, data_path, file_type = ' '):  
        
        if not (file_type.lower() == ' ' or  \
                file_type.lower() == 'gps' or  \
                file_type.lower() == 'glonass' or \
                file_type.lower() == 'beidou' or \
                file_type.lower() == 'galileo' ):
            raise RuntimeError('rinex_nav.__init__: Not currently implemented for : ' + file_type)
        
        self.file_type = file_type.lower()
        self.data_path = data_path
        self.sat_prns = []
        self.nr_sats = len( self.sat_prns) # Derived quantity, but so commonly used it is worthwhile 
        
        self.records = []      # List of records contained within the file        
        self.read()            # Read the data file
        
        
    def read( self):
        print( "Reading rinex nav file: " + self.data_path)
        
        rinex_nav_file = open(self.data_path)
        rinex_nav_content = rinex_nav_file.read()
        rinex_nav_file.close()
        
        # Separate the records
        rinex_nav_lines = rinex_nav_content.splitlines() 
        
        # Parse the header data
        
        self.version = float(rinex_nav_lines[0][0:9])
        self.type = rinex_nav_lines[0][20:60]
        
        if self.version == 2.0 and self.file_type.lower() == 'gps':
            self.read_rinex_nav_gps_2( rinex_nav_lines)
        elif self.version == 2.01 and self.file_type.lower() == 'gps':
            self.read_rinex_nav_gps_2( rinex_nav_lines)
        elif self.version == 2.01 and self.file_type.lower() == 'glonass':
            raise RuntimeError('RINEX_nav.read_sp3_a: NOt yet implemented for GLONASS data')
        else:
            print( 'RINEX_nav.read: rinex_nav_' + str(self.version) + ' parser not yet implemented!')
    
    def read_rinex_nav_gps_2( self, rinex_nav_lines):

        # rinex_nav_a format only - this function can also read rinex_nav_b files that are mislabeled as rinex_nav_a
        # as long as they are labelled correctly

        # Read the file header (just skip by everything but the leap seconds for now
       
        n_lines = 0 
        for rinex_nav_line in rinex_nav_lines:    
            if rinex_nav_line[60:73].lower() == "leap seconds":
                self.leap_sec = int(rinex_nav_line[0:6])
            elif rinex_nav_line[60:73].lower() == "end of header":
                n_lines+=1
                break
             
            n_lines+=1
            
        n_records = (len(rinex_nav_lines)-n_lines)/8
        
        if not np.floor( n_records) == n_records:
            raise  RuntimeError('rinex_nav.read: Incorrect file size, records corrupted')
        
        n_records = int(n_records)
                
        # Allocate the memory needed to hold the structured array
        self.records = np.zeros(n_records, dtype= gps_d_type)
        
        # Read the records into the structured array
        for i in range(n_records):
            self.records[i]['prn'] = int(rinex_nav_lines[n_lines][0:2])
            
            # Read the satellite data
            prn = rinex_nav_line[0:2]
            
            # Determine the year of the t_oc
            if datetime.now().year >= 2078:
                # If this ever gets raised then I'll eat my socks, bit we'll test anyway
                # (It would mean this code is still in use 58 years after creation)
                raise  RuntimeError('rinex_nav.read: Ambiguity in epoch year')
            
            # Create a full representation of the year
            year = int(rinex_nav_lines[n_lines][3:5])
            if  year < 78:
                year += 2000    
            else: 
                year += 1900
                
            # Get the date and time as a timedelta in GPS second
                
            self.records[i]['t_oc'] = datetime(year, \
                             int(rinex_nav_lines[n_lines][6:8]), \
                             int(rinex_nav_lines[n_lines][9:11]), \
                             int(rinex_nav_lines[n_lines][12:14]), \
                             int(rinex_nav_lines[n_lines][15:17]), \
                             int(rinex_nav_lines[n_lines][17:20]), \
                             int(rinex_nav_lines[n_lines][21:22])*100000) - datetime(1980,1,6,0,0,0)
            
            if not self.records[i]['prn'] in self.sat_prns:
                self.sat_prns.append( self.records[i]['prn'])
        
            # The clock correction polynomial coefficients

            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['a_0']       = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['a_1']       = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['a_2']       = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 1
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['iode']      = np.float64(rinex_nav_lines[n_lines][ 3:22])
            self.records[i]['c_rs']      = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['D_n']       = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['M_0']       = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 2
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['c_uc']      = np.float64(rinex_nav_lines[n_lines][3:22])
            self.records[i]['e']         = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['c_us']      = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['a_rt']      = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 3
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['t_oe']      = np.float64(rinex_nav_lines[n_lines][3:22])
            self.records[i]['c_ic']      = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['W_0']       = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['c_is']      = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 4
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['i_0']       = np.float64(rinex_nav_lines[n_lines][ 3:22])
            self.records[i]['c_rc']      = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['w']         = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['W_dot']     = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 5
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['i_dot']     = np.float64(rinex_nav_lines[n_lines][ 3:22])
            self.records[i]['L2_code']   = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['toe_week']  = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['L2_flags']  = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 6
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['sv_acc']    = np.float64(rinex_nav_lines[n_lines][ 3:22])
            self.records[i]['sv_health'] = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['t_gd']      = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['iodc']      = np.float64(rinex_nav_lines[n_lines][60:79])
            # BROADCAST ORBIT - 7
            n_lines+=1
            rinex_nav_lines[n_lines] = rinex_nav_lines[n_lines].replace("D", "E")
            self.records[i]['sv_acc']    = np.float64(rinex_nav_lines[n_lines][ 3:22])
            self.records[i]['sv_health'] = np.float64(rinex_nav_lines[n_lines][22:41])
            self.records[i]['t_gd']      = np.float64(rinex_nav_lines[n_lines][41:60])
            self.records[i]['iodc']      = np.float64(rinex_nav_lines[n_lines][60:79])

            n_lines+=1


