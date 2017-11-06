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

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

parser = SafeConfigParser()

# Set input file
if len(sys.argv) < 2:
    input_file = 'settings.ini'
    print 'Using default input file: ' + input_file
else:
    input_file = sys.argv[1]
    print 'Using specified input file: ' + input_file

# Read in directory information
#parser.read('settings_test.ini')
parser.read(input_file)
PROJ_NAME = parser.get('settings', 'proj_name')

# command-line executable for GSFLOW (just used to print message)
GSFLOW_exe = parser.get('settings', 'gsflow_exe')
#DEM = parser.get('settings', 'DEM') # name of file w/ topography data (in GIS data directory)
GISinput_dir = parser.get('settings', 'GISinput_dir')
climate_data_file = parser.get('settings', 'climate_data_file')

PROJ_CODE=PROJ_NAME.replace(" ", "") # remove blank spaces

gsflow_path_simdir = parser.get('settings', 'gsflow_path_simdir')
gsflow_simdir = gsflow_path_simdir + slashstr + PROJ_CODE

# 1: for spinup (starts with steady-state run), 2: for restart (run AFTER spinup)
sw_1spinup_2restart = int(parser.get('settings', 'sw_1spinup_2restart'))
if sw_1spinup_2restart == 2:
    # point to files created from spinup run
    restart_PRMSfil = parser.get('settings', 'restart_PRMSfil')
    restart_MODfil = parser.get('settings', 'restart_MODfil')


# for relative pathname
PRMSinput_dir_rel = 'inputs' + slashstr + 'PRMS_GSFLOW' 
MODFLOWinput_dir_rel = 'inputs' + slashstr + 'MODFLOW_NWT'
PRMSoutput_dir_rel = 'outputs' + slashstr + 'PRMS_GSFLOW' 
MODFLOWoutput_dir_rel = 'outputs' + slashstr + 'MODFLOW_NWT'

# full pathnames
control_dir = gsflow_simdir + slashstr + 'control' 
PRMSinput_dir = gsflow_simdir + slashstr + PRMSinput_dir_rel 
MODFLOWinput_dir = gsflow_simdir + slashstr + MODFLOWinput_dir_rel
PRMSoutput_dir = gsflow_simdir + slashstr + PRMSoutput_dir_rel
MODFLOWoutput_dir = gsflow_simdir + slashstr + MODFLOWoutput_dir_rel


# create directories if they do not exist:
if not os.path.isdir(control_dir):
    os.makedirs(control_dir)   
if not os.path.isdir(PRMSinput_dir):
    os.makedirs(PRMSinput_dir)
if not os.path.isdir(MODFLOWinput_dir):
    os.makedirs(MODFLOWinput_dir)
if not os.path.isdir(PRMSoutput_dir):
    os.makedirs(PRMSoutput_dir)
if not os.path.isdir(MODFLOWoutput_dir):
    os.makedirs(MODFLOWoutput_dir)

    

# -- problem-specifc variables

# either single value for constant K or name of file with array [m/d]
# 
hydcond0 = parser.get('custom_params', 'hydcond') 
# either single value for spatially constant finf or name of file with array [m/d]
finf0 = parser.get('custom_params', 'finf') 

START_DATE = parser.get('domain', 'start_date')
END_DATE = parser.get('domain', 'end_date')

NLAY = int(parser.get('domain', 'NLAY'))
DZ_str = parser.get('domain', 'DZ')  # for NLAY>1, this is comma-separated list (e.g., dz=50, 100), spaces don't matter
value = DZ_str.split(',')
DZ = []
for ii in range(NLAY):
    try:
      DZ.append(int(value[ii]))
    except:
      DZ.append(float(value[ii]))


hydcond0_arr_str = hydcond0
value = hydcond0_arr_str.split(',')
hydcond0 = []
for ii in range(NLAY):
    try:
      hydcond0.append(float(value[ii]))
    except:
      hydcond0 = hydcond0_arr_str

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

