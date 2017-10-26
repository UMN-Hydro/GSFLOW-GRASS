# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 07:43:05 2017

@author: gcng
"""

import os
import platform
import settings_test

if platform.system() == 'Linux':
    slashstr = '/'
elif platform.system() == 'Windows':
    slashstr = '\\'

model_mode = 'GSFLOW'
cmd_str = settings_test.control_dir + slashstr + settings_test.PROJ_CODE + '_' + model_mode

if platform.system() == 'Linux':
    cmd_str = cmd_str + '.sh > out.txt'
elif platform.system() == 'Windows':
    cmd_str = cmd_str + '.bat'

os.system(cmd_str)


