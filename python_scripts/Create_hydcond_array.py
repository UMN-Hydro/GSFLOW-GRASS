# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 18:33:11 2017

@author: gcng
"""

# Create_hydcond_array

import numpy as np
import pandas as pd
import platform
import settings_test

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

surfz_fil = settings_test.GISinput_dir + slashstr + settings_test.DEM + '.asc'

NLAY = 1

hydcond_fil = 'hydcond_test.txt'

sw_scheme = 2

# %%

f = open(surfz_fil, 'r')
sdata = {}
for i in range(6):
    line = f.readline()
    line = line.rstrip() # remove newline characters
    key, value = line.split(': ')
    try:
      value = int(value)
    except:
      value = float(value)
    sdata[key] = value
f.close()
    
NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
NROW = sdata['rows'] 
NCOL = sdata['cols']

TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)


hydcond = np.ones([NROW, NCOL, NLAY]) * 0.1 # default

#%%
# ----- Based on elevation -----

if sw_scheme == 1:
    # - domain dimensions, maybe already in surfz_fil and botm_fil{}?
    # NLAY = 2;
    # surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
    
    hydcond0 = np.copy(hydcond[:,:,0])    
    hydcond0[TOP<4500] = 0.25
    hydcond0[TOP<4200] = 0.4
    hydcond0[TOP<4000] = 0.5
    hydcond[:,:,0] = np.copy(hydcond0)
    
    
    

#%%
# ----- Based on stream channel -----

if sw_scheme == 2:
    reach_fil = settings_test.GISinput_dir + slashstr + 'reaches_tmp.txt'
    reach_data_all = pd.read_csv(reach_fil)       # 

    hydcond[reach_data_all.loc[:,'IRCH']-1, reach_data_all.loc[:,'JRCH']-1, 0] = 0.6


#%% Write to File
    
fobj = open(hydcond_fil, 'w+')
for ii in range(NLAY):
    fobj.write('Layer %d \n' % (ii+1))
    np.savetxt(fobj, hydcond[:,:,ii], delimiter=' ', fmt='%10g')        

fobj.close()
