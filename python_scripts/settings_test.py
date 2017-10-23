#!/usr/bin/python2.7

from ConfigParser import SafeConfigParser
import platform

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
DEM = parser.get('settings', 'DEM')
GISinput_dir = parser.get('settings', 'GISinput_dir')

START_DATE = parser.get('settings', 'start_date')
END_DATE = parser.get('settings', 'end_date')

PROJ_CODE=PROJ_NAME.replace(" ", "") # remove blank spaces

gsflow_simdir = parser.get('settings', 'gsflow_simdir')
control_dir = gsflow_simdir + slashstr + 'control' 
PRMSinput_dir = gsflow_simdir + slashstr + 'inputs' + slashstr + 'PRMS'
MODFLOWinput_dir = gsflow_simdir + slashstr + 'inputs' + slashstr + 'MODFLOW_NWT'
PRMSoutput_dir = gsflow_simdir + slashstr + 'outputs' + slashstr + 'PRMS' # eventually rename?  It's really GSFLOW outputs
MODFLOWoutput_dir = gsflow_simdir + slashstr + 'outputs' + slashstr + 'MODFLOW_NWT'

#LOCAL_DIR = parser.get('settings', 'local_dir')
#control_dir = parser.get('settings', 'control_dir')
#PRMSinput_dir = parser.get('settings', 'PRMSinput_dir')
#MODFLOWinput_dir = parser.get('settings', 'MODFLOWinput_dir')
#PRMSoutput_dir = parser.get('settings', 'PRMSoutput_dir')
#MODFLOWoutput_dir = parser.get('settings', 'MODFLOWoutput_dir')



# -- problem-specifc variables
# only ones of these will be read in
parser.read('custom_params.ini')
hydcond0 = parser.get('custom_params', 'hydcond') # either single value for constant K or name of file with array [m/d]
finf0 = parser.get('custom_params', 'finf') # either single value for spatially constant finf or name of file with array [m/d]


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

