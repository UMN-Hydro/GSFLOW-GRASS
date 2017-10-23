# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 18:33:11 2017

@author: gcng
"""

# Create_hydcond_array

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # matlab-like plots
import platform
import settings_test

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

surfz_fil = settings_test.GISinput_dir + slashstr + settings_test.DEM + '.asc'

NLAY = settings_test.NLAY

hydcond_fil = 'hydcond_test.txt'

sw_scheme = 2 # 1: based on elev, 2: based on streams, 3: based on both

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

# - space discretization
DELR = (NSEW[2]-NSEW[3])/NCOL # width of column [m]
DELC = (NSEW[0]-NSEW[1])/NROW # height of row [m]

TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)


hydcond = np.ones([NROW, NCOL, NLAY]) * 0.1 # default

#%%
# ----- Based on elevation -----

if sw_scheme == 1 or sw_scheme == 3:
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

if sw_scheme == 2 or sw_scheme == 3:
    reach_fil = settings_test.GISinput_dir + slashstr + 'reaches_tmp.txt'
    reach_data_all = pd.read_csv(reach_fil)       # 

    # Set stream cell only
#    hydcond[reach_data_all.loc[:,'IRCH']-1, reach_data_all.loc[:,'JRCH']-1, 0] = 0.6

    # hydcond depends on distance from stream
    row_ind = reach_data_all.loc[:,'IRCH']-1
    col_ind = reach_data_all.loc[:,'JRCH']-1
    xcoord = DELR * np.arange(NCOL)
    ycoord = DELC * np.arange(NROW)
    xstrm = xcoord[col_ind]
    ystrm = ycoord[row_ind]
    xcoord_ar = np.kron(np.ones((NROW,1)), xcoord)
    ycoord_ar = np.kron(np.ones((NCOL,1)), ycoord)
    ycoord_ar = ycoord_ar.transpose()    
    
    dx = np.ceil(np.maximum(DELR, DELC))
    npix_stream_buffer = 1
    buffer_dist = np.arange(1,npix_stream_buffer+1) * dx  # up to npix pixels from stream
    buffer_hydcond = np.array([0.4, 0.2])
    
    ind = np.argsort(buffer_dist)[::-1]
    buffer_dist = np.copy(buffer_dist[ind])
    buffer_hydcond = np.copy(buffer_hydcond[ind])
    
    strm_hydcond = 0.6
    
    hydcond0 = np.copy(hydcond[:,:,0])    
    
    # buffer distances from stream:
    for d_i in range(np.size(buffer_dist)):
        for strm_i in range(np.size(xstrm)): 
            dist = ((xcoord_ar-xstrm[strm_i])**2 + (ycoord_ar-ystrm[strm_i])**2)**0.5
            hydcond0[dist <= buffer_dist[d_i]] = buffer_hydcond[d_i]
                
    hydcond0[row_ind, col_ind] = strm_hydcond # stream
    
    hydcond[:,:,0] = hydcond0

# %% Plot
#ax = fig.add_subplot(2,2,1)
#im = ax.imshow(TOP_to_plot)

for ilay in range(NLAY):
    fig = plt.figure(figsize=(12,12))
#    plt.subplot(2,2,ilay+1)
    im = plt.imshow(hydcond[:,:,ilay])
#    im.set_clim(3800, 6200)
    fig.colorbar(im, orientation='horizontal')
    plt.title('BOTM lay' + str(ilay+1));


#%% Write to File
    
fobj = open(hydcond_fil, 'w+')
for ii in range(NLAY):
    fobj.write('Layer %d \n' % (ii+1))
    np.savetxt(fobj, hydcond[:,:,ii], delimiter=' ', fmt='%10g')        

fobj.close()
