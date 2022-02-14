import os
import urllib.request
from datetime import datetime, timezone, timedelta
import numpy as np

    ## semme J. Dijkstra  1/31/2020 Created SP3 class with SP3_a reader


class SP3:
    """A Class for handling SP3 Ephemeris Data"""
    
    # An object of this class will hold the data of one sp3 file
    
    def __init__(self, data_path, file_type = ' '):  
        
        if not (file_type.lower() == ' ' or  \
                file_type.lower() == 'gps' or  \
                file_type.lower() == 'glonass' or \
                file_type.lower() == 'beidou' or \
                file_type.lower() == 'galileo' ):
            raise RuntimeError('SP3.__init__: Not currently implemented for : ' + file_type)
        
        self.file_type = file_type
        self.data_path = data_path
        
        # Header variables
        # initial velaues are from the .sp3_c definition example
        
        self.version = 'c'
        self.pos_or_vel = 'V'
        self.start_date = datetime(2001,8,8,0,0,0) # Time to the nearest microsecond
        self.start_s = 0  # Needed to get to the 10 second resolution used for GNSS time intervals
        self.nr_epochs = 0
        self.data_used = 'ORBIT'
        self.coord_sys = 'IGS97'
        self.orbit_type = 'FIT'
        self.agency = 'NGS'
        self.gps_week = 1126
        self.gps_seconds = 259200.00000000
        self.epoch_interval = 900
        self.mod_julian_date = 52129  # the modified Julian Day Start (where 44244 
                                      #represents GPS zero time -- January 6, 1980)
            
        self.fractional_day = 0.0000000000000
        self.nr_sats = 12
        self.sat_ids = list()
        self.sat_accuracies = list()
        self.time_system = 'GPS'
        self.pos_vel_base = 1.2500000
        self.clock_rate_base = 1.025000000
        self.header_comments = list()
        
        # Epoch data variables
        
        self.epoch_times = list()
        self.epoch_seconds = list()
        self.epoch_pos_t_data = list()
        self.epoch_vel_data = list()
        
        self.read()
        
        
    def read( self):
        print( "Reading SP3 file: " + self.data_path)
        
        sp3_file = open(self.data_path)
        sp3_content = sp3_file.read()
        sp3_file.close()
        
        # Separate the records
        sp3_lines = sp3_content.splitlines() 
        
        # Parse the header data
        
        i = 0
        self.version = sp3_lines[i][1].lower()
        
        if self.version == 'a':
            self.read_sp3_a( sp3_lines)
        elif self.version == 'b':
            self.read_sp3_b( sp3_lines)
        elif self.version == 'c':
            self.read_sp3_c( sp3_lines)
        else:
            print( 'SP3: SP3_' + self.version + ' parser not yet implemented!')
    
    def read_sp3_a( self, sp3_lines):
        
        # sp3_a format only - this function can also read sp3_b files that are mislabeled as sp3_a
        
        if (self.file_type.lower() == ' '):
            s_error = 'SP3.__init__: For sp3_a files the system type must be provided \n'
            s_error += 'Allowable system types\n'
            s_error += 'G or GPS\n'
            s_error += 'R or GLONASS\n'
            s_error += 'L or LEO or "Low Earth"\n'
            s_error += 'E or GALILEO\n'
            s_error += 'C or COMPASS\n'
            s_error += 'J or QZSS\n'
            s_error += 'B or BEIDOU'
            s_error += 'M or MIXED'
            raise RuntimeError(s_error)
        
        #### REMEMBER TO UPDATE ALL read_sp3_XX METHODS IF NECESARRY ####

        self.pos_or_vel      = sp3_lines[0][2].lower()     # Position or Velocity File
              
        self.start_date      = datetime(int(sp3_lines[0][3:7]), \
                                   int(sp3_lines[0][8:10]), \
                                   int(sp3_lines[0][11:13]), \
                                   int(sp3_lines[0][14:16]), \
                                   int(sp3_lines[0][17:19]), \
                                   int(sp3_lines[0][20:22]), \
                                   int(sp3_lines[0][23:31]))
        
        self.start_second    = sp3_lines[0][20:31]
        self.nr_epochs       = int(sp3_lines[0][32:39])
        self.data_used       = sp3_lines[0][40:45]
        self.coord_sys       = sp3_lines[0][46:51]
        self.orbit_type      = sp3_lines[0][52:55]
        self.agency          = sp3_lines[0][56:60]
        self.gps_week        = int(sp3_lines[1][3:7])
        self.gps_seconds     = float(sp3_lines[1][8:23])
        self.epoch_interval  = float(sp3_lines[1][24:38])
        self.mod_julian_date = int(sp3_lines[1][39:44])
        self.fractional_day  = float(sp3_lines[1][45:60])
        self.nr_sats         = int(sp3_lines[2][4:6])
        
        # This is kind of useless, as it is repeated in the epochs, but oh well
        # It would address the issue of a satellite disappearing during the day, but that 
        # has not ever happened as far as I know
 
        self.sat_ids.append( sp3_lines[2][ 9:12])
        self.sat_ids.append( sp3_lines[2][12:15])
        self.sat_ids.append( sp3_lines[2][15:18])
        self.sat_ids.append( sp3_lines[2][18:21])
        self.sat_ids.append( sp3_lines[2][21:24])
        self.sat_ids.append( sp3_lines[2][24:27])
        self.sat_ids.append( sp3_lines[2][27:30])
        self.sat_ids.append( sp3_lines[2][30:33])
        self.sat_ids.append( sp3_lines[2][33:36])
        self.sat_ids.append( sp3_lines[2][36:39])
        self.sat_ids.append( sp3_lines[2][39:42])
        self.sat_ids.append( sp3_lines[2][42:45])
        self.sat_ids.append( sp3_lines[2][45:48])
        self.sat_ids.append( sp3_lines[2][48:51])
        self.sat_ids.append( sp3_lines[2][51:54])
        self.sat_ids.append( sp3_lines[2][54:57])
        self.sat_ids.append( sp3_lines[2][57:60])
        self.sat_ids.append( sp3_lines[3][ 9:12])
        self.sat_ids.append( sp3_lines[3][12:15])
        self.sat_ids.append( sp3_lines[3][15:18])
        self.sat_ids.append( sp3_lines[3][18:21])
        self.sat_ids.append( sp3_lines[3][21:24])
        self.sat_ids.append( sp3_lines[3][24:27])
        self.sat_ids.append( sp3_lines[3][27:30])
        self.sat_ids.append( sp3_lines[3][30:33])
        self.sat_ids.append( sp3_lines[3][33:36])
        self.sat_ids.append( sp3_lines[3][36:39])
        self.sat_ids.append( sp3_lines[3][39:42])
        self.sat_ids.append( sp3_lines[3][42:45])
        self.sat_ids.append( sp3_lines[3][45:48])
        self.sat_ids.append( sp3_lines[3][48:51])
        self.sat_ids.append( sp3_lines[3][51:54])
        self.sat_ids.append( sp3_lines[3][54:57])
        self.sat_ids.append( sp3_lines[3][57:60])
        self.sat_ids.append( sp3_lines[4][ 9:12])
        self.sat_ids.append( sp3_lines[4][12:15])
        self.sat_ids.append( sp3_lines[4][15:18])
        self.sat_ids.append( sp3_lines[4][18:21])
        self.sat_ids.append( sp3_lines[4][21:24])
        self.sat_ids.append( sp3_lines[4][24:27])
        self.sat_ids.append( sp3_lines[4][27:30])
        self.sat_ids.append( sp3_lines[4][30:33])
        self.sat_ids.append( sp3_lines[4][33:36])
        self.sat_ids.append( sp3_lines[4][36:39])
        self.sat_ids.append( sp3_lines[4][39:42])
        self.sat_ids.append( sp3_lines[4][42:45])
        self.sat_ids.append( sp3_lines[4][45:48])
        self.sat_ids.append( sp3_lines[4][48:51])
        self.sat_ids.append( sp3_lines[4][51:54])
        self.sat_ids.append( sp3_lines[4][54:57])
        self.sat_ids.append( sp3_lines[4][57:60])
        self.sat_ids.append( sp3_lines[5][ 9:12])
        self.sat_ids.append( sp3_lines[5][12:15])
        self.sat_ids.append( sp3_lines[5][15:18])
        self.sat_ids.append( sp3_lines[5][18:21])
        self.sat_ids.append( sp3_lines[5][21:24])
        self.sat_ids.append( sp3_lines[5][24:27])
        self.sat_ids.append( sp3_lines[5][27:30])
        self.sat_ids.append( sp3_lines[5][30:33])
        self.sat_ids.append( sp3_lines[5][33:36])
        self.sat_ids.append( sp3_lines[5][36:39])
        self.sat_ids.append( sp3_lines[5][39:42])
        self.sat_ids.append( sp3_lines[5][42:45])
        self.sat_ids.append( sp3_lines[5][45:48])
        self.sat_ids.append( sp3_lines[5][48:51])
        self.sat_ids.append( sp3_lines[5][51:54])
        self.sat_ids.append( sp3_lines[5][54:57])
        self.sat_ids.append( sp3_lines[5][57:60])
        self.sat_ids.append( sp3_lines[6][ 9:12])
        self.sat_ids.append( sp3_lines[6][12:15])
        self.sat_ids.append( sp3_lines[6][15:18])
        self.sat_ids.append( sp3_lines[6][18:21])
        self.sat_ids.append( sp3_lines[6][21:24])
        self.sat_ids.append( sp3_lines[6][24:27])
        self.sat_ids.append( sp3_lines[6][27:30])
        self.sat_ids.append( sp3_lines[6][30:33])
        self.sat_ids.append( sp3_lines[6][33:36])
        self.sat_ids.append( sp3_lines[6][36:39])
        self.sat_ids.append( sp3_lines[6][39:42])
        self.sat_ids.append( sp3_lines[6][42:45])
        self.sat_ids.append( sp3_lines[6][45:48])
        self.sat_ids.append( sp3_lines[6][48:51])
        self.sat_ids.append( sp3_lines[6][51:54])
        self.sat_ids.append( sp3_lines[6][54:57])
        self.sat_ids.append( sp3_lines[6][57:60])

        self.sat_accuracies.append( 2**float(sp3_lines[7][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][57:60])/1000) 
        
        # Replace all unknown accuracies with nan
        
        for i in range( len( self.sat_accuracies)):
            if self.sat_accuracies[i] == 0.001:
                self.sat_accuracies[i] = np.nan
                
                
        if not sp3_lines[12][9:12].lower() == 'ccc':
            raise RuntimeError('SP3.read_sp3_a: File format indicated as SP3_a, but at least SP3_b')
        
        self.time_system = 'GPS'  # Always for sp3_a files 
        self.pos_vel_base = 2   # Always for sp3_a files
        self.clock_rate_base = 2  # Always for sp3_a files
        self.header_comments.append(sp3_lines[18][3:60])
        self.header_comments.append(sp3_lines[19][3:60])
        self.header_comments.append(sp3_lines[20][3:60])
        self.header_comments.append(sp3_lines[21][3:60])
        
        # Predict file size from number of records assuming that there are no correclation data blocks
        # For now we will only read the clock and position records
        
        if 22+self.nr_epochs*(1*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 1
        elif 22+self.nr_epochs*(2*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 2
        elif 22+self.nr_epochs*(3*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 3
        else:
            print(len(sp3_lines))
            raise RuntimeError('SP3.read_sp3_a: Predicted file size does not match actual file size')
                  
        # Cycle through the epoch data
        
        epoch_indexes = np.arange(22, \
                                  22+self.nr_epochs*(nr_recordtypes*self.nr_sats+1), \
                                  nr_recordtypes*self.nr_sats+1)
        

        for i in epoch_indexes:     # Cycle through all epochs
            self.epoch_pos_t_data.append(np.zeros([self.nr_sats,34]))  # 34 SP3_c has up to 34 fields 
            
            # Parse epoch header record
            self.epoch_times.append( datetime(int(sp3_lines[i][3:7]), \
                                   int(sp3_lines[i][8:10]), \
                                   int(sp3_lines[i][11:13]), \
                                   int(sp3_lines[i][14:16]), \
                                   int(sp3_lines[i][17:19]), \
                                   int(sp3_lines[i][20:22]), \
                                   int(sp3_lines[i][23:31])))
            
            self.epoch_seconds.append( float(sp3_lines[i][20:31]))
        
            k = 0
            for j in range(i+1,i+nr_recordtypes*self.nr_sats + 1):
                if sp3_lines[j][0:1] == 'P':
                    # Parse position and clock record
                    # GPS = 0 is default
                    if sp3_lines[j][1:2] == ' ' or sp3_lines[j][1:2] == 'G':
                        self.epoch_pos_t_data[-1][k,0]=0 #GPS
                    elif sp3_lines[j][1:2] == 'R':
                        self.epoch_pos_t_data[-1][k,0]=1 #Glonass
                    else:
                         raise RuntimeError( 'SP3.read_sp3_a: Not GPS or GLONASS')
                    self.epoch_pos_t_data[-1][k,1] = sp3_lines[j][2:4]
                    self.epoch_pos_t_data[-1][k,2] = float(sp3_lines[j][4:18])*1000
                    self.epoch_pos_t_data[-1][k,3] = float(sp3_lines[j][18:32])*1000
                    self.epoch_pos_t_data[-1][k,4] = float(sp3_lines[j][32:46])*1000
                    self.epoch_pos_t_data[-1][k,5] = float(sp3_lines[j][46:60])
                elif  sp3_lines[j][0:1] == 'V':
                    # GPS = 0 is default
                    if sp3_lines[j][1:2] == ' ' and not self.epoch_pos_t_data[-1][k,0]==0:
                        raise RuntimeError( 'SP3.read_sp3_a: Mixing satellite type in records')
                    elif sp3_lines[j][1:2] == 'G'  and not self.epoch_pos_t_data[-1][k,0]==1:
                        raise RuntimeError( 'SP3.Mixing satellite types in records')
                    elif not (sp3_lines[j][1:2] == ' ' or sp3_lines[j][1:2] == 'G'):
                        raise RuntimeError( 'SP3.read_sp3_a: Not GPS or GLONASS')
                    if not self.epoch_pos_t_data[-1][k,1] == float(sp3_lines[j][2:4]):
                        raise RuntimeError( 'SP3.read_sp3_a: satellite id mismatch')
                    self.epoch_pos_t_data[-1][k,10] = float(sp3_lines[j][4:18])*1000
                    self.epoch_pos_t_data[-1][k,11] = float(sp3_lines[j][18:32])*1000
                    self.epoch_pos_t_data[-1][k,12] = float(sp3_lines[j][32:46])*1000
                    self.epoch_pos_t_data[-1][k,13] = float(sp3_lines[j][46:60])
                else: 
                    raise RuntimeError( 'SP3.read_sp3_a: Not a known sp3_a record')
                
    
                k += 1


    def read_sp3_b( self, sp3_lines):
        
        #### REMEMBER TO UPDATE ALL read_sp3_XX METHODS IF NECESARRY ####
        
        ## The only difference with sp3_a is how the satellites are defined

        self.pos_or_vel      = sp3_lines[0][2].lower()     # Position or Velocity File
              
        self.start_date      = datetime(int(sp3_lines[0][3:7]), \
                                   int(sp3_lines[0][8:10]), \
                                   int(sp3_lines[0][11:13]), \
                                   int(sp3_lines[0][14:16]), \
                                   int(sp3_lines[0][17:19]), \
                                   int(sp3_lines[0][20:22]), \
                                   int(sp3_lines[0][23:31]))
        
        self.start_second    = sp3_lines[0][20:31]
        self.nr_epochs       = int(sp3_lines[0][32:39])
        self.data_used       = sp3_lines[0][40:45]
        self.coord_sys       = sp3_lines[0][46:51]
        self.orbit_type      = sp3_lines[0][52:55]
        self.agency          = sp3_lines[0][56:60]
        self.gps_week        = int(sp3_lines[1][3:7])
        self.gps_seconds     = float(sp3_lines[1][8:23])
        self.epoch_interval  = float(sp3_lines[1][24:38])
        self.mod_julian_date = int(sp3_lines[1][39:44])
        self.fractional_day  = float(sp3_lines[1][45:60])
        self.nr_sats         = int(sp3_lines[2][4:6])
        
        # This is kind of useless, as it is repeated in the epochs, but oh well
        # It would address the issue of a satellite disappearing during the day, but that 
        # has not ever happened as far as I know
 
        self.sat_ids.append( sp3_lines[2][ 9:12])
        self.sat_ids.append( sp3_lines[2][12:15])
        self.sat_ids.append( sp3_lines[2][15:18])
        self.sat_ids.append( sp3_lines[2][18:21])
        self.sat_ids.append( sp3_lines[2][21:24])
        self.sat_ids.append( sp3_lines[2][24:27])
        self.sat_ids.append( sp3_lines[2][27:30])
        self.sat_ids.append( sp3_lines[2][30:33])
        self.sat_ids.append( sp3_lines[2][33:36])
        self.sat_ids.append( sp3_lines[2][36:39])
        self.sat_ids.append( sp3_lines[2][39:42])
        self.sat_ids.append( sp3_lines[2][42:45])
        self.sat_ids.append( sp3_lines[2][45:48])
        self.sat_ids.append( sp3_lines[2][48:51])
        self.sat_ids.append( sp3_lines[2][51:54])
        self.sat_ids.append( sp3_lines[2][54:57])
        self.sat_ids.append( sp3_lines[2][57:60])
        self.sat_ids.append( sp3_lines[3][ 9:12])
        self.sat_ids.append( sp3_lines[3][12:15])
        self.sat_ids.append( sp3_lines[3][15:18])
        self.sat_ids.append( sp3_lines[3][18:21])
        self.sat_ids.append( sp3_lines[3][21:24])
        self.sat_ids.append( sp3_lines[3][24:27])
        self.sat_ids.append( sp3_lines[3][27:30])
        self.sat_ids.append( sp3_lines[3][30:33])
        self.sat_ids.append( sp3_lines[3][33:36])
        self.sat_ids.append( sp3_lines[3][36:39])
        self.sat_ids.append( sp3_lines[3][39:42])
        self.sat_ids.append( sp3_lines[3][42:45])
        self.sat_ids.append( sp3_lines[3][45:48])
        self.sat_ids.append( sp3_lines[3][48:51])
        self.sat_ids.append( sp3_lines[3][51:54])
        self.sat_ids.append( sp3_lines[3][54:57])
        self.sat_ids.append( sp3_lines[3][57:60])
        self.sat_ids.append( sp3_lines[4][ 9:12])
        self.sat_ids.append( sp3_lines[4][12:15])
        self.sat_ids.append( sp3_lines[4][15:18])
        self.sat_ids.append( sp3_lines[4][18:21])
        self.sat_ids.append( sp3_lines[4][21:24])
        self.sat_ids.append( sp3_lines[4][24:27])
        self.sat_ids.append( sp3_lines[4][27:30])
        self.sat_ids.append( sp3_lines[4][30:33])
        self.sat_ids.append( sp3_lines[4][33:36])
        self.sat_ids.append( sp3_lines[4][36:39])
        self.sat_ids.append( sp3_lines[4][39:42])
        self.sat_ids.append( sp3_lines[4][42:45])
        self.sat_ids.append( sp3_lines[4][45:48])
        self.sat_ids.append( sp3_lines[4][48:51])
        self.sat_ids.append( sp3_lines[4][51:54])
        self.sat_ids.append( sp3_lines[4][54:57])
        self.sat_ids.append( sp3_lines[4][57:60])
        self.sat_ids.append( sp3_lines[5][ 9:12])
        self.sat_ids.append( sp3_lines[5][12:15])
        self.sat_ids.append( sp3_lines[5][15:18])
        self.sat_ids.append( sp3_lines[5][18:21])
        self.sat_ids.append( sp3_lines[5][21:24])
        self.sat_ids.append( sp3_lines[5][24:27])
        self.sat_ids.append( sp3_lines[5][27:30])
        self.sat_ids.append( sp3_lines[5][30:33])
        self.sat_ids.append( sp3_lines[5][33:36])
        self.sat_ids.append( sp3_lines[5][36:39])
        self.sat_ids.append( sp3_lines[5][39:42])
        self.sat_ids.append( sp3_lines[5][42:45])
        self.sat_ids.append( sp3_lines[5][45:48])
        self.sat_ids.append( sp3_lines[5][48:51])
        self.sat_ids.append( sp3_lines[5][51:54])
        self.sat_ids.append( sp3_lines[5][54:57])
        self.sat_ids.append( sp3_lines[5][57:60])
        self.sat_ids.append( sp3_lines[6][ 9:12])
        self.sat_ids.append( sp3_lines[6][12:15])
        self.sat_ids.append( sp3_lines[6][15:18])
        self.sat_ids.append( sp3_lines[6][18:21])
        self.sat_ids.append( sp3_lines[6][21:24])
        self.sat_ids.append( sp3_lines[6][24:27])
        self.sat_ids.append( sp3_lines[6][27:30])
        self.sat_ids.append( sp3_lines[6][30:33])
        self.sat_ids.append( sp3_lines[6][33:36])
        self.sat_ids.append( sp3_lines[6][36:39])
        self.sat_ids.append( sp3_lines[6][39:42])
        self.sat_ids.append( sp3_lines[6][42:45])
        self.sat_ids.append( sp3_lines[6][45:48])
        self.sat_ids.append( sp3_lines[6][48:51])
        self.sat_ids.append( sp3_lines[6][51:54])
        self.sat_ids.append( sp3_lines[6][54:57])
        self.sat_ids.append( sp3_lines[6][57:60])

        self.sat_accuracies.append( 2**float(sp3_lines[7][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][57:60])/1000)
        
        # Replace all unknown accuracies with nan
        
        for i in range( len( self.sat_accuracies)):
            if self.sat_accuracies[i] == 0.001:
                self.sat_accuracies[i] = np.nan

        self.file_type = sp3_lines[12][3:5]
        self.time_system = sp3_lines[12][9:12]
        
        # Warn if the time used is not GPS time
        
        if not self.time_system == 'GPS':
            print('SP3.read_sp3_b: ' + self.time_system+' time used in file: '+self.data_path)
            
        self.pos_vel_base = sp3_lines[14][3:13]
        self.clock_rate_base = sp3_lines[14][14:26]
        self.header_comments.append(sp3_lines[18][3:60])
        self.header_comments.append(sp3_lines[19][3:60])
        self.header_comments.append(sp3_lines[20][3:60])
        self.header_comments.append(sp3_lines[21][3:60])

        # Predict file size from number of records
        
        if 22+self.nr_epochs*(1*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 1
        elif 22+self.nr_epochs*(2*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 2
        elif 22+self.nr_epochs*(3*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 3
        elif 22+self.nr_epochs*(4*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 4
        else:
            print(len(sp3_lines))
            raise RuntimeError('SP3.read_sp3_b: Predicted file size does not match actual file size')
        
        # Cycle through the epoch data
        
        epoch_indexes = np.arange(22, \
                                  22+self.nr_epochs*(nr_recordtypes*self.nr_sats+1), \
                                  nr_recordtypes*self.nr_sats+1)

        for i in epoch_indexes:
            self.epoch_pos_t_data.append(np.zeros([self.nr_sats,34]))  # 34 SP3_c has up to 34 fields 
            
            # Parse epoch header record
            self.epoch_times.append( datetime(int(sp3_lines[i][3:7]), \
                                   int(sp3_lines[i][8:10]), \
                                   int(sp3_lines[i][11:13]), \
                                   int(sp3_lines[i][14:16]), \
                                   int(sp3_lines[i][17:19]), \
                                   int(sp3_lines[i][20:22]), \
                                   int(sp3_lines[i][23:31])))
            
            self.epoch_seconds.append( float(sp3_lines[i][20:31]))
            
            k = 0
            for j in range(i+1,i+nr_recordtypes*self.nr_sats + 1):
                
                if sp3_lines[j][0:1] == 'P':
                    # Parse position and clock record
                    # GPS = 0 is default
                    if sp3_lines[j][1:2] == ' ' or sp3_lines[j][1:2] == 'G':
                        self.epoch_pos_t_data[-1][k,0]=0 #GPS
                    elif sp3_lines[j][1:2] == 'R':
                        self.epoch_pos_t_data[-1][k,0]=1 #Glonass
                    else:
                        print( 'file line : '+str(j+1))
                        print(sp3_lines[j])
                        raise RuntimeError( 'SP3.read_sp3_b: Not GPS or GLONASS')
                    self.epoch_pos_t_data[-1][k,1] = sp3_lines[j][2:4]
                    self.epoch_pos_t_data[-1][k,2] = float(sp3_lines[j][4:18])*1000
                    self.epoch_pos_t_data[-1][k,3] = float(sp3_lines[j][18:32])*1000
                    self.epoch_pos_t_data[-1][k,4] = float(sp3_lines[j][32:46])*1000
                    self.epoch_pos_t_data[-1][k,5] = float(sp3_lines[j][46:60])
                elif  sp3_lines[j][0:1] == 'V':
                    # GPS = 0 is default
                    if sp3_lines[j][1:2] == ' ' and not self.epoch_pos_t_data[-1][k,0]==0:
                        raise RuntimeError( 'SP3.read_sp3_b: Mixing satellite type in records')
                    elif sp3_lines[j][1:2] == 'G'  and not self.epoch_pos_t_data[-1][k,0]==1:
                        raise RuntimeError( 'SP3.read_sp3_b: Mixing satellite types in records')
                    elif not (sp3_lines[j][1:2] == ' ' or sp3_lines[j][1:2] == 'G'):
                        raise RuntimeError( 'SP3.read_sp3_b: Not GPS or GLONASS')
                    if not self.epoch_pos_t_data[-1][k,1] == float(sp3_lines[j][2:4]):
                        raise RuntimeError( 'SP3.read_sp3_b: satellite id mismatch')
                    self.epoch_pos_t_data[-1][k,10] = float(sp3_lines[j][4:18])*1000
                    self.epoch_pos_t_data[-1][k,11] = float(sp3_lines[j][18:32])*1000
                    self.epoch_pos_t_data[-1][k,12] = float(sp3_lines[j][32:46])*1000
                    self.epoch_pos_t_data[-1][k,13] = float(sp3_lines[j][46:60])
                elif  sp3_lines[j][0:2] == 'EP':
                    pass
                elif  sp3_lines[j][0:2] == 'EV':
                    pass
                else: 
                    raise RuntimeError( 'SP3.read_sp3_b: ' + \
                                       sp3_lines[j][0:1] + \
                                       ' Not a known sp3_b record')
                k += 1          
        
    def read_sp3_c( self, sp3_lines):
        
        #### REMEMBER TO UPDATE ALL read_sp3_XX METHODS IF NECESARRY ####

        self.pos_or_vel      = sp3_lines[0][2].lower()     # Position or Velocity File
              
        self.start_date      = datetime(int(sp3_lines[0][3:7]), \
                                   int(sp3_lines[0][8:10]), \
                                   int(sp3_lines[0][11:13]), \
                                   int(sp3_lines[0][14:16]), \
                                   int(sp3_lines[0][17:19]), \
                                   int(sp3_lines[0][20:22]), \
                                   int(sp3_lines[0][23:31]))
        
        self.start_second    = sp3_lines[0][20:31]
        self.nr_epochs       = int(sp3_lines[0][32:39])
        self.data_used       = sp3_lines[0][40:45]
        self.coord_sys       = sp3_lines[0][46:51]
        self.orbit_type      = sp3_lines[0][52:55]
        self.agency          = sp3_lines[0][56:60]
        self.gps_week        = int(sp3_lines[1][3:7])
        self.gps_seconds     = float(sp3_lines[1][8:23])
        self.epoch_interval  = float(sp3_lines[1][24:38])
        self.mod_julian_date = int(sp3_lines[1][39:44])
        self.fractional_day  = float(sp3_lines[1][45:60])
        self.nr_sats         = int(sp3_lines[2][4:6])
        
        # This is kind of useless, as it is repeated in the epochs, but oh well
        # It would address the issue of a satellite disappearing during the day, but that 
        # has not ever happened as far as I know
 
        self.sat_ids.append( sp3_lines[2][ 9:12])
        self.sat_ids.append( sp3_lines[2][12:15])
        self.sat_ids.append( sp3_lines[2][15:18])
        self.sat_ids.append( sp3_lines[2][18:21])
        self.sat_ids.append( sp3_lines[2][21:24])
        self.sat_ids.append( sp3_lines[2][24:27])
        self.sat_ids.append( sp3_lines[2][27:30])
        self.sat_ids.append( sp3_lines[2][30:33])
        self.sat_ids.append( sp3_lines[2][33:36])
        self.sat_ids.append( sp3_lines[2][36:39])
        self.sat_ids.append( sp3_lines[2][39:42])
        self.sat_ids.append( sp3_lines[2][42:45])
        self.sat_ids.append( sp3_lines[2][45:48])
        self.sat_ids.append( sp3_lines[2][48:51])
        self.sat_ids.append( sp3_lines[2][51:54])
        self.sat_ids.append( sp3_lines[2][54:57])
        self.sat_ids.append( sp3_lines[2][57:60])
        self.sat_ids.append( sp3_lines[3][ 9:12])
        self.sat_ids.append( sp3_lines[3][12:15])
        self.sat_ids.append( sp3_lines[3][15:18])
        self.sat_ids.append( sp3_lines[3][18:21])
        self.sat_ids.append( sp3_lines[3][21:24])
        self.sat_ids.append( sp3_lines[3][24:27])
        self.sat_ids.append( sp3_lines[3][27:30])
        self.sat_ids.append( sp3_lines[3][30:33])
        self.sat_ids.append( sp3_lines[3][33:36])
        self.sat_ids.append( sp3_lines[3][36:39])
        self.sat_ids.append( sp3_lines[3][39:42])
        self.sat_ids.append( sp3_lines[3][42:45])
        self.sat_ids.append( sp3_lines[3][45:48])
        self.sat_ids.append( sp3_lines[3][48:51])
        self.sat_ids.append( sp3_lines[3][51:54])
        self.sat_ids.append( sp3_lines[3][54:57])
        self.sat_ids.append( sp3_lines[3][57:60])
        self.sat_ids.append( sp3_lines[4][ 9:12])
        self.sat_ids.append( sp3_lines[4][12:15])
        self.sat_ids.append( sp3_lines[4][15:18])
        self.sat_ids.append( sp3_lines[4][18:21])
        self.sat_ids.append( sp3_lines[4][21:24])
        self.sat_ids.append( sp3_lines[4][24:27])
        self.sat_ids.append( sp3_lines[4][27:30])
        self.sat_ids.append( sp3_lines[4][30:33])
        self.sat_ids.append( sp3_lines[4][33:36])
        self.sat_ids.append( sp3_lines[4][36:39])
        self.sat_ids.append( sp3_lines[4][39:42])
        self.sat_ids.append( sp3_lines[4][42:45])
        self.sat_ids.append( sp3_lines[4][45:48])
        self.sat_ids.append( sp3_lines[4][48:51])
        self.sat_ids.append( sp3_lines[4][51:54])
        self.sat_ids.append( sp3_lines[4][54:57])
        self.sat_ids.append( sp3_lines[4][57:60])
        self.sat_ids.append( sp3_lines[5][ 9:12])
        self.sat_ids.append( sp3_lines[5][12:15])
        self.sat_ids.append( sp3_lines[5][15:18])
        self.sat_ids.append( sp3_lines[5][18:21])
        self.sat_ids.append( sp3_lines[5][21:24])
        self.sat_ids.append( sp3_lines[5][24:27])
        self.sat_ids.append( sp3_lines[5][27:30])
        self.sat_ids.append( sp3_lines[5][30:33])
        self.sat_ids.append( sp3_lines[5][33:36])
        self.sat_ids.append( sp3_lines[5][36:39])
        self.sat_ids.append( sp3_lines[5][39:42])
        self.sat_ids.append( sp3_lines[5][42:45])
        self.sat_ids.append( sp3_lines[5][45:48])
        self.sat_ids.append( sp3_lines[5][48:51])
        self.sat_ids.append( sp3_lines[5][51:54])
        self.sat_ids.append( sp3_lines[5][54:57])
        self.sat_ids.append( sp3_lines[5][57:60])
        self.sat_ids.append( sp3_lines[6][ 9:12])
        self.sat_ids.append( sp3_lines[6][12:15])
        self.sat_ids.append( sp3_lines[6][15:18])
        self.sat_ids.append( sp3_lines[6][18:21])
        self.sat_ids.append( sp3_lines[6][21:24])
        self.sat_ids.append( sp3_lines[6][24:27])
        self.sat_ids.append( sp3_lines[6][27:30])
        self.sat_ids.append( sp3_lines[6][30:33])
        self.sat_ids.append( sp3_lines[6][33:36])
        self.sat_ids.append( sp3_lines[6][36:39])
        self.sat_ids.append( sp3_lines[6][39:42])
        self.sat_ids.append( sp3_lines[6][42:45])
        self.sat_ids.append( sp3_lines[6][45:48])
        self.sat_ids.append( sp3_lines[6][48:51])
        self.sat_ids.append( sp3_lines[6][51:54])
        self.sat_ids.append( sp3_lines[6][54:57])
        self.sat_ids.append( sp3_lines[6][57:60])

        self.sat_accuracies.append( 2**float(sp3_lines[7][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[7][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[8][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[9][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[10][57:60])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][ 9:12])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][12:15])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][15:18])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][18:21])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][21:24])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][24:27])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][27:30])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][30:33])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][33:36])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][36:39])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][39:42])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][42:45])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][45:48])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][48:51])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][51:54])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][54:57])/1000)
        self.sat_accuracies.append( 2**float(sp3_lines[11][57:60])/1000)     

        # Replace all unknown accuracies with nan
        
        for i in range( len( self.sat_accuracies)):
            if self.sat_accuracies[i] == 0.001:
                self.sat_accuracies[i] = np.nan
                 
        self.file_type = sp3_lines[12][3:5]
        self.time_system = sp3_lines[12][9:12]
        
        # Warn if the time used is not GPS time
        
        if not self.time_system == 'GPS':
            print('SP3.read_sp3_c: ' + self.time_system+' time used in file: '+self.data_path)
            
        self.pos_vel_base = float(sp3_lines[14][3:13])
        self.clock_rate_base = float(sp3_lines[14][14:26])
        self.header_comments.append(sp3_lines[18][3:60])
        self.header_comments.append(sp3_lines[19][3:60])
        self.header_comments.append(sp3_lines[20][3:60])
        self.header_comments.append(sp3_lines[21][3:60])
        

        # Predict file size from number of records
        
        if   22 + self.nr_epochs*(1*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 1
        elif 22 + self.nr_epochs*(2*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 2
        elif 22 + self.nr_epochs*(3*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 3
        elif 22 + self.nr_epochs*(4*self.nr_sats+1)+1 == len(sp3_lines):
            nr_recordtypes = 4
        else:
            print(len(sp3_lines))
            raise RuntimeError('SP3.read_sp3_c: Predicted file size does not match actual file size')
        
        # Cycle through the epoch data
        
        epoch_indexes = np.arange(22, \
                                  22+self.nr_epochs*(nr_recordtypes*self.nr_sats+1), \
                                  nr_recordtypes*self.nr_sats+1)

        for i in epoch_indexes:
            self.epoch_pos_t_data.append(np.zeros([self.nr_sats,34]))  # 34 SP3_c has up to 34 fields 
            
            # Parse epoch header record
            self.epoch_times.append( datetime(int(sp3_lines[i][3:7]), \
                                   int(sp3_lines[i][8:10]), \
                                   int(sp3_lines[i][11:13]), \
                                   int(sp3_lines[i][14:16]), \
                                   int(sp3_lines[i][17:19]), \
                                   int(sp3_lines[i][20:22]), \
                                   int(sp3_lines[i][23:31])))

            
            self.epoch_seconds.append( float(sp3_lines[i][20:31]))

            
            k = 0
            for j in range(i+1,i+nr_recordtypes*self.nr_sats + 1):
                if sp3_lines[j][0:1] == 'P':
                    # Parse position and clock record
                    # GPS = 0 is default
                    if sp3_lines[j][1:2] == ' ' or sp3_lines[j][1:2] == 'G':
                        self.epoch_pos_t_data[-1][k,0]=0 #GPS
                    elif sp3_lines[j][1:2] == 'R':
                        self.epoch_pos_t_data[-1][k,0]=1 #Glonass
                    else:
                         raise RuntimeError( 'SP3.read_sp3_c: Not GPS or GLONASS')
                    self.epoch_pos_t_data[-1][k,1] = sp3_lines[j][2:4]
                    self.epoch_pos_t_data[-1][k,2] = float(sp3_lines[j][4:18])*1000
                    self.epoch_pos_t_data[-1][k,3] = float(sp3_lines[j][18:32])*1000
                    self.epoch_pos_t_data[-1][k,4] = float(sp3_lines[j][32:46])*1000
                    self.epoch_pos_t_data[-1][k,5] = float(sp3_lines[j][46:60])
                    
                    if not self.epoch_pos_t_data[-1][k,5] == 999999.999999:
                        if not sp3_lines[j][61:63] == '  ':
                            self.epoch_pos_t_data[-1][k,6] = self.pos_vel_base**float(sp3_lines[j][61:63])/1000
                        else:
                            self.epoch_pos_t_data[-1][k,6] = np.nan
                        if not sp3_lines[j][64:66] == '  ':
                            self.epoch_pos_t_data[-1][k,7] = self.pos_vel_base**float(sp3_lines[j][64:66])/1000
                        else:
                            self.epoch_pos_t_data[-1][k,7] = np.nan
                        if not sp3_lines[j][67:69] == '  ':
                            self.epoch_pos_t_data[-1][k,8] = self.pos_vel_base**float(sp3_lines[j][67:69])/1000
                        else:
                            self.epoch_pos_t_data[-1][k,8] = np.nan
                        if not sp3_lines[j][70:73] == '   ':
                            self.epoch_pos_t_data[-1][k,9] = self.clock_rate_base**float(sp3_lines[j][70:73])/1000
                        else:
                            self.epoch_pos_t_data[-1][k,9] = np.nan
                        
                        
                    else:
                        self.epoch_pos_t_data[-1][k,6] = np.nan
                        self.epoch_pos_t_data[-1][k,7] = np.nan
                        self.epoch_pos_t_data[-1][k,8] = np.nan
                        self.epoch_pos_t_data[-1][k,9] = np.nan
                elif  sp3_lines[j][0:1] == 'V':
                    # GPS = 0 is default
                    if sp3_lines[j][1:2] == ' ' and not self.epoch_pos_t_data[-1][k,0]==0:
                        raise RuntimeError( 'SP3.read_sp3_c: Mixing satellite type in records')
                    elif sp3_lines[j][1:2] == 'G'  and not self.epoch_pos_t_data[-1][k,0]==0:
                        raise RuntimeError( 'SP3.read_sp3_c: Mixing satellite types  on line '+str(j+1))
                    elif not (sp3_lines[j][1:2] == ' ' or sp3_lines[j][1:2] == 'G'):
                        raise RuntimeError( 'SP3.read_sp3_c: Not GPS or GLONASS on line '+str(j+1))
                    if not self.epoch_pos_t_data[-1][k,1] == float(sp3_lines[j][2:4]):
                        raise RuntimeError( 'SP3.read_sp3_c: satellite id mismatch')
                    self.epoch_pos_t_data[-1][k,10] = float(sp3_lines[j][4:18])*1000
                    self.epoch_pos_t_data[-1][k,11] = float(sp3_lines[j][18:32])*1000
                    self.epoch_pos_t_data[-1][k,12] = float(sp3_lines[j][32:46])*1000
                    self.epoch_pos_t_data[-1][k,13] = float(sp3_lines[j][46:60])
                elif  sp3_lines[j][0:2] == 'EP':
                    self.epoch_pos_t_data[-1][k,14] = float(sp3_lines[j][4:8])/1000
                    self.epoch_pos_t_data[-1][k,15] = float(sp3_lines[j][9:13])/1000
                    self.epoch_pos_t_data[-1][k,16] = float(sp3_lines[j][14:18])/1000
                    self.epoch_pos_t_data[-1][k,17] = float(sp3_lines[j][19:26])
                    self.epoch_pos_t_data[-1][k,18] = float(sp3_lines[j][27:35])
                    self.epoch_pos_t_data[-1][k,19] = float(sp3_lines[j][36:44])
                    self.epoch_pos_t_data[-1][k,20] = float(sp3_lines[j][45:53])
                    self.epoch_pos_t_data[-1][k,21] = float(sp3_lines[j][54:62])
                    self.epoch_pos_t_data[-1][k,22] = float(sp3_lines[j][63:71])
                    self.epoch_pos_t_data[-1][k,23] = float(sp3_lines[j][72:80])
                elif  sp3_lines[j][0:2] == 'EV':
                    self.epoch_pos_t_data[-1][k,24] = float(sp3_lines[j][4:8])/10000
                    self.epoch_pos_t_data[-1][k,25] = float(sp3_lines[j][9:13])/10000
                    self.epoch_pos_t_data[-1][k,26] = float(sp3_lines[j][14:18])/10000
                    self.epoch_pos_t_data[-1][k,27] = float(sp3_lines[j][19:26])
                    self.epoch_pos_t_data[-1][k,28] = float(sp3_lines[j][27:35])
                    self.epoch_pos_t_data[-1][k,29] = float(sp3_lines[j][36:44])
                    self.epoch_pos_t_data[-1][k,30] = float(sp3_lines[j][45:53])
                    self.epoch_pos_t_data[-1][k,31] = float(sp3_lines[j][54:62])
                    self.epoch_pos_t_data[-1][k,32] = float(sp3_lines[j][63:71])
                    self.epoch_pos_t_data[-1][k,33] = float(sp3_lines[j][72:80])

                else: 
                    raise RuntimeError( 'SP3.read_sp3_c: ' + \
                                       sp3_lines[j][0:1] + \
                                       ' Not a known sp3_c record')

                k += 1