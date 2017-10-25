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
else:
    slashstr = '\\'

cmd_str = settings_test.control_dir + slashstr + settings_test.PROJ_CODE + '_GSFLOW.sh'

if platform.system() == 'Linux':
    cmd_str = cmd_str + ' > out.txt'

os.system(cmd_str)

