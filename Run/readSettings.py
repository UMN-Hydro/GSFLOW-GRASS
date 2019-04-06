#!/usr/bin/python2.7

# To run within Spyder:
# runfile('/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/settings_test.py', args='settings_new.ini',  wdir='/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts')
# (make sure current directory contains the input file)
#
# To run at command line:


from ConfigParser import SafeConfigParser
import platform
import os
import shutil
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
        

        # Read in directory information
        #parser.read('settings_test.ini')
        parser.read(input_file)
        self.PROJ_NAME = parser.get('paths', 'proj_name')
        self.PROJ_CODE=self.PROJ_NAME.replace(" ", "") # remove blank spaces

        # command-line executable for GSFLOW (just used to print message)
        self.GSFLOW_exe = parser.get('paths', 'gsflow_exe')
        self.GSFLOW_ver = parser.get('paths', 'gsflow_ver')
        
        self.gsflow_path_simdir = parser.get('paths', 'gsflow_path_simdir')
        self.gsflow_simdir = self.gsflow_path_simdir + slashstr + self.PROJ_CODE
        self.GISinput_dir = self.gsflow_simdir + slashstr + 'GIS'
        
        # 1: for spinup (starts with steady-state run), 2: for restart (run AFTER spinup)
        self.sw_1spinup_2restart = int(parser.get('run_mode', 'sw_1spinup_2restart'))
        if self.sw_1spinup_2restart == 2:
            # point to files created from spinup run
            self.init_PRMSfil = parser.get('run_mode', 'init_PRMSfil')
            self.init_MODfil = parser.get('run_mode', 'init_MODfil')


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
        
        # GRASS GIS core variables
        self.GIS_output_rootdir = self.gsflow_simdir + slashstr + 'GIS'
        
        # GRASS GIS DEM input raster map
        self.DEM_input = parser.get('elevation_inputs', 'DEM_file_path_to_import')

        # GRASS GIS LAND-SURFACE raster maps
        self.LAND_COVER_file_path_to_import = parser.get('land-surface_inputs', 'LAND_COVER_file_path_to_import')
        self.SOIL_file_path_to_import = parser.get('land-surface_inputs', 'SOIL_file_path_to_import')
        self.cov_type_uniform = parser.get('land-surface_inputs', 'cov_type')
        self.soil_type_uniform = parser.get('land-surface_inputs', 'soil_type')

        # GRASS GIS drainage variables
        self.drainage_threshold = parser.get('GRASS_drainage', 'threshold_drainage_area_meters2')
        self.flow_weights = parser.get('GRASS_drainage', 'flow_weights')
        self.MODFLOW_grid_resolution = parser.get('GRASS_drainage', 'MODFLOW_grid_resolution_meters')
        self.outlet_point_x = parser.get('GRASS_drainage', 'outlet_point_x')
        self.outlet_point_y = parser.get('GRASS_drainage', 'outlet_point_y')
        
        # GRASS GIS hydraulics variables
        self.icalc = parser.get('GRASS_hydraulics', 'icalc')
        self.channel_Mannings_n = parser.get('GRASS_hydraulics', 'channel_Mannings_n')
        self.channel_Mannings_n_grid = parser.get('GRASS_hydraulics', 'channel_Mannings_n_grid')
        self.channel_Mannings_n_vector = parser.get('GRASS_hydraulics', 'channel_Mannings_n_vector')
        self.channel_Mannings_n_vector_col = parser.get('GRASS_hydraulics', 'channel_Mannings_n_vector_col')
        self.overbank_Mannings_n = parser.get('GRASS_hydraulics', 'overbank_Mannings_n')
        self.channel_width = parser.get('GRASS_hydraulics', 'channel_width')
        self.channel_width_vector = parser.get('GRASS_hydraulics', 'channel_width_vector')
        self.channel_width_vector_col = parser.get('GRASS_hydraulics', 'channel_width_vector_col')
        self.floodplain_width = parser.get('GRASS_hydraulics', 'floodplain_width')
        self.floodplain_width_vector = parser.get('GRASS_hydraulics', 'floodplain_width_vector')
        self.floodplain_width_vector_col = parser.get('GRASS_hydraulics', 'floodplain_width_vector_col')

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

        self.fl_print_climate_hru = int(parser.get('climate_inputs', 'fl_print_climate_hru'))
        if self.fl_print_climate_hru == 1:        
            self.climate_data_file = parser.get('climate_inputs', 'climate_data_file')
        elif self.fl_print_climate_hru == 0:
            self.climate_hru_dir = parser.get('climate_inputs', 'climate_hru_dir')            
            climatehru_fils = ['tmin.day', 'tmax.day', 'precip.day', 'empty.day']
            for item in climatehru_fils:
                source_fil = self.climate_hru_dir + slashstr + item
                shutil.copy2(source_fil, self.PRMSinput_dir)
        else:
            sys.exit("Must choose 1 or 0 for 'fl_print_climate_hru'")

        # either single value for constant K or name of file with array [m/d]
        self.fl_create_hydcond = int(parser.get('hydrogeologic_inputs', 'fl_create_hydcond'))
        hydcond0 = parser.get('hydrogeologic_inputs', 'hydcond') 
        # either single value for spatially constant finf or name of file with array [m/d]
        self.finf0 = parser.get('hydrogeologic_inputs', 'finf') 

        self.START_DATE = parser.get('time', 'start_date')
        self.END_DATE = parser.get('time', 'end_date')
        if self.sw_1spinup_2restart == 2:
            self.INIT_START_DATE = parser.get('time', 'init_start_date')

        self.NLAY = int(parser.get('hydrogeologic_inputs', 'NLAY'))
        DZ_str = parser.get('hydrogeologic_inputs', 'DZ')  # for NLAY>1, this is comma-separated list (e.g., dz=50, 100), spaces don't matter
        value = DZ_str.split(',')
        self.DZ = []
        for ii in range(self.NLAY):
            try:
              self.DZ.append(int(value[ii]))
            except:
              self.DZ.append(float(value[ii]))


        self.hydcond0 = hydcond0.split(',')

