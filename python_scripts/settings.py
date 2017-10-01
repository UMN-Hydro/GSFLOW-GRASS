#!/usr/bin/python2.7

from ConfigParser import SafeConfigParser

parser = SafeConfigParser()
parser.read('settings.ini')
LOCAL_DIR = parser.get('settings', 'local_dir')
PROJ_NAME = parser.get('settings', 'proj_name')
control_dir = parser.get('settings', 'control_dir')
PRMSinput_dir = parser.get('settings', 'PRMSinput_dir')
MODFLOWinput_dir = parser.get('settings', 'MODFLOWinput_dir')
PRMSoutput_dir = parser.get('settings', 'MODFLOWinput_dir')
# command-line executable for GSFLOW (just used to print message)
GSFLOW_exe = parser.get('settings', 'gsflow_exe') + '/gsflow'
DEM = parser.get('settings', 'DEM')

PROJ_CODE=PROJ_NAME.replace(" ", "")
