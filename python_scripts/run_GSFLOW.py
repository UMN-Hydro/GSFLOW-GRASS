# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 07:43:05 2017

@author: gcng
"""

import os
import platform
from readSettings import Settings
import sys

# Set input file
if len(sys.argv) < 2:
    settings_input_file = 'settings.ini'
    print 'Using default input file: ' + settings_input_file
else:
    settings_input_file = sys.argv[1]
    print 'Using specified input file: ' + settings_input_file
Settings = Settings(settings_input_file)

if platform.system() == 'Linux':
    slashstr = '/'
elif platform.system() == 'Windows':
    slashstr = '\\'

model_mode = 'GSFLOW'
cmd_str = Settings.control_dir + slashstr + Settings.PROJ_CODE + '_' + model_mode

if platform.system() == 'Linux':
    cmd_str = cmd_str + '.sh > out.txt'
elif platform.system() == 'Windows':
    cmd_str = cmd_str + '.bat'

os.system(cmd_str)


