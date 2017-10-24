#!/usr/bin/python2.7

from ConfigParser import SafeConfigParser
import platform
import os

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

parser = SafeConfigParser()

# Read in directory information
parser.read('settings_test.ini')
PROJ_NAME = parser.get('settings', 'proj_name')

# command-line executable for GSFLOW (just used to print message)
GSFLOW_exe = parser.get('settings', 'gsflow_exe') + slashstr + 'gsflow'
DEM = parser.get('settings', 'DEM') # name of file w/ topography data
GISinput_dir = parser.get('settings', 'GISinput_dir')
climate_data_file = parser.get('settings', 'climate_data_file')

PROJ_CODE=PROJ_NAME.replace(" ", "") # remove blank spaces

gsflow_simdir = parser.get('settings', 'gsflow_simdir')

# 1: for spinup, 2: for restart
sw_1spinup_2restart = int(parser.get('settings', 'sw_1spinup_2restart'))
if sw_1spinup_2restart == 2:
    restart_PRMSfil = parser.get('settings', 'restart_PRMSfil')
    restart_MODfil = parser.get('settings', 'restart_MODfil')


# for relative pathname
PRMSinput_dir_rel = 'inputs' + slashstr + 'PRMS_GSFLOW' 
MODFLOWinput_dir_rel = 'inputs' + slashstr + 'MODFLOW_NWT'
PRMSoutput_dir_rel = 'outputs' + slashstr + 'PRMS_GSFLOW' 
MODFLOWoutput_dir_rel = 'outputs' + slashstr + 'MODFLOW_NWT'

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
# only ones of these will be read in
parser.read('custom_params.ini')
hydcond0 = parser.get('custom_params', 'hydcond') # either single value for constant K or name of file with array [m/d]
finf0 = parser.get('custom_params', 'finf') # either single value for spatially constant finf or name of file with array [m/d]

START_DATE = parser.get('domain', 'start_date')
END_DATE = parser.get('domain', 'end_date')

NLAY = int(parser.get('domain', 'NLAY'))
DZ_str = parser.get('domain', 'DZ')  # for NLAY>1, this is comma-separated array
value = DZ_str.split(',')
DZ = []
for ii in range(NLAY):
    try:
      DZ.append(int(value[ii]))
    except:
      DZ.append(float(value[ii]))


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

