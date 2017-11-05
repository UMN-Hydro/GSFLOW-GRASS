#!/usr/bin/python2.7

# To run within Spyder:
# runfile('/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/settings_test.py', args='settings_new.ini',  wdir='/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts')
# (make sure current directory contains the input file)
#
# To run at command line:


from ConfigParser import SafeConfigParser
import platform
import os
import sys # to read in command line arguments

class Settings(object):

    def __init__(self, input_file):
        """
        Parses settings.ini for GSFLOW
        """
        
        if platform.system() == 'Linux':
            slashstr = '/'
        else:
            slashstr = '\\'

        parser = SafeConfigParser()
        
        """
        # Set input file
        if len(sys.argv) < 2:
            input_file = 'settings.ini'
            print 'Using default input file: ' + input_file
        else:
            input_file = sys.argv[1]
            print 'Using specified input file: ' + input_file
        """

        # Read in directory information
        #parser.read('settings_test.ini')
        parser.read(input_file)
        self.PROJ_NAME = parser.get('settings', 'proj_name')

        # command-line executable for GSFLOW (just used to print message)
        self.GSFLOW_exe = parser.get('settings', 'gsflow_exe')
        #DEM = parser.get('settings', 'DEM') # name of file w/ topography data (in GIS data directory)
        self.GISinput_dir = parser.get('settings', 'GISinput_dir')
        self.climate_data_file = parser.get('settings', 'climate_data_file')

        self.PROJ_CODE=self.PROJ_NAME.replace(" ", "") # remove blank spaces

        self.gsflow_simdir = parser.get('settings', 'gsflow_simdir')

        # 1: for spinup (starts with steady-state run), 2: for restart (run AFTER spinup)
        self.sw_1spinup_2restart = int(parser.get('settings', 'sw_1spinup_2restart'))
        if self.sw_1spinup_2restart == 2:
            # point to files created from spinup run
            self.restart_PRMSfil = parser.get('settings', 'restart_PRMSfil')
            self.restart_MODfil = parser.get('settings', 'restart_MODfil')


        # for relative pathname
        self.PRMSinput_dir_rel = 'inputs' + slashstr + 'PRMS_GSFLOW' 
        self.MODFLOWinput_dir_rel = 'inputs' + slashstr + 'MODFLOW_NWT'
        self.PRMSoutput_dir_rel = 'outputs' + slashstr + 'PRMS_GSFLOW' 
        self.MODFLOWoutput_dir_rel = 'outputs' + slashstr + 'MODFLOW_NWT'

        # full pathnames
        self.control_dir = self.gsflow_simdir + slashstr + 'control' 
        self.PRMSinput_dir = self.gsflow_simdir + slashstr + self.PRMSinput_dir_rel 
        self.MODFLOWinput_dir = self.gsflow_simdir + slashstr + self.MODFLOWinput_dir_rel
        self.PRMSoutput_dir = self.gsflow_simdir + slashstr + self.PRMSoutput_dir_rel
        self.MODFLOWoutput_dir = self.gsflow_simdir + slashstr + self.MODFLOWoutput_dir_rel


        # create directories if they do not exist:
        if not os.path.isdir(self.control_dir):
            os.makedirs(self.control_dir)   
        if not os.path.isdir(self.PRMSinput_dir):
            os.makedirs(self.PRMSinput_dir)
        if not os.path.isdir(self.MODFLOWinput_dir):
            os.makedirs(self.MODFLOWinput_dir)
        if not os.path.isdir(self.PRMSoutput_dir):
            os.makedirs(self.PRMSoutput_dir)
        if not os.path.isdir(self.MODFLOWoutput_dir):
            os.makedirs(self.MODFLOWoutput_dir)

            

        # -- problem-specifc variables

        # either single value for constant K or name of file with array [m/d]
        # 
        hydcond0 = parser.get('custom_params', 'hydcond') 
        # either single value for spatially constant finf or name of file with array [m/d]
        self.finf0 = parser.get('custom_params', 'finf') 

        self.START_DATE = parser.get('domain', 'start_date')
        self.END_DATE = parser.get('domain', 'end_date')

        self.NLAY = int(parser.get('domain', 'NLAY'))
        DZ_str = parser.get('domain', 'DZ')  # for NLAY>1, this is comma-separated list (e.g., dz=50, 100), spaces don't matter
        value = DZ_str.split(',')
        self.DZ = []
        for ii in range(self.NLAY):
            try:
              self.DZ.append(int(value[ii]))
            except:
              self.DZ.append(float(value[ii]))


        hydcond0_arr_str = hydcond0
        value = hydcond0_arr_str.split(',')
        self.hydcond0 = []
        for ii in range(self.NLAY):
            try:
              self.hydcond0.append(float(value[ii]))
            except:
              self.hydcond0 = hydcond0_arr_str

        ## move to MODFLOW_NWT_lib_Shullcas_test.py
        #import numpy as np
        #NROW = 25
        #NCOL = 14
        #NLAY = 2
        #hydcond = np.ones((NROW,NCOL,NLAY))
        #try:
        #   float(hydcond0)
        #   hydcond = float(hydcond0) * hydcond
        #except ValueError:
        #   for ii in range(NLAY):
        #        hydcond[:,:,ii] = np.genfromtxt(hydcond0, skip_header=1+ii*(NROW+1), \
        #        max_rows=NROW, dtype=float)
        #        
        #finf = np.ones((NROW,NCOL))
        #try:
        #   float(finf0)
        #   finf = float(finf0) * finf
        #except ValueError:
        #   for ii in range(NLAY):
        #        finf[:,:,ii] = np.genfromtxt(finf0, skip_header=1+ii*(NROW+1), \
        #        max_rows=NROW, dtype=float)



        #hydcond_fil = parser.get('custom_params', 'hydcond_array_fil') # file with array, NROW x NCOL
        #finf_const = parser.get('custom_params', 'finf') # constant finf (infiltration for MODFLOW's steady-state stress period
        #finf_fil = parser.get('custom_params', 'finf_array_fil') # file with array, NROW x NCOL
        #NLAY = parser.get('custom_params', 'NLAY_MODFLOW') # number of MODFLOW layers
        #dz = parser.get('custom_params', 'NLAY_MODFLOW') # 

