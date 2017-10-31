# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 19:04:13 2017

Will plot:
  IBOUND, TOP, BOTM, and Hyd Cond 
  (can easily be expanded to also plot vertical hyd cond, Ss, and Sy)

@author: gcng
"""


import platform
import numpy as np
from matplotlib import pyplot as plt
import settings_test

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'


# *** Change file names as needed
dis_fil = settings_test.MODFLOWinput_dir + slashstr + settings_test.PROJ_CODE + '.dis'
ba6_fil = settings_test.MODFLOWinput_dir + slashstr + settings_test.PROJ_CODE + '.ba6'
flo_fil = settings_test.MODFLOWinput_dir + slashstr + settings_test.PROJ_CODE + '.upw'


#%%
f = open(dis_fil, 'r')
for i in range(3): # first 2 lines are comments
    line = f.readline().rstrip()
line = line.split()
NLAY = int(line[0]) 
NROW = int(line[1]) 
NCOL = int(line[2]) 
NPER = int(line[3]) 
ITMUNI = int(line[4]) 
LENUNI = int(line[5])    

print(NROW, NCOL)

line = f.readline().rstrip()
LAYCBD = np.array(line.split(), float)

line = f.readline().rstrip()
line = line.split()
DELR = float(line[1])
line = f.readline().rstrip()
line = line.split()
DELC = float(line[1])
f.close() 

TOP = np.genfromtxt(dis_fil, skip_header=7, max_rows=NROW, dtype=float)
BOTM = np.zeros((NROW, NCOL, NLAY),float);
for ii in range(NLAY):
    BOTM[:,:,ii] = np.genfromtxt(dis_fil, skip_header=7+(ii+1)*(NROW+1), \
    max_rows=NROW, dtype=float)

IBOUND = np.genfromtxt(ba6_fil, skip_header=3, max_rows=NROW, dtype=float)

# -- find boundary cells
IBOUND0 = np.copy(IBOUND)
IBOUND0.astype(int)
IBOUND0[IBOUND0>0] = 1 # active cells
IBOUND0[IBOUND0<0] = 0 # constant head cells
# (IBOUND0 = 0 for no flow)

IBOUNDin = IBOUND0[1:-1,1:-1]
IBOUNDu = IBOUND0[0:-2,1:-1] # up
IBOUNDd = IBOUND0[2:,1:-1] # down
IBOUNDl = IBOUND0[1:-1,0:-2] # left
IBOUNDr = IBOUND0[1:-1,2:] # right

# - inner boundary (of inner grid domain,i.e. domain[1:-1,1:-1])
ind_bound_in = np.logical_and(IBOUNDin==1, np.logical_or(np.logical_or(np.logical_or(IBOUNDin-IBOUNDu==1, IBOUNDin-IBOUNDd==1), \
IBOUNDin-IBOUNDl==1), IBOUNDin-IBOUNDr==1))
# - outer boundary (of inner grid domain, i.e. domain[1:-1,1:-1])
ind_bound_out = np.logical_and(IBOUNDin==0, np.logical_or(np.logical_or(np.logical_or(IBOUNDin-IBOUNDu==-1, IBOUNDin-IBOUNDd==-1), \
IBOUNDin-IBOUNDl==-1), IBOUNDin-IBOUNDr==-1))

ind_bound_out = np.logical_and(IBOUNDin==0, np.logical_or(np.logical_or(np.logical_or(IBOUNDin-IBOUNDu==-1, IBOUNDin-IBOUNDd==-1), \
IBOUNDin-IBOUNDl==-1), IBOUNDin-IBOUNDr==-1))

HY = np.zeros((NROW, NCOL, NLAY),float); # hyd conductivity
VKA = np.zeros((NROW, NCOL, NLAY),float); # vertical hyd conductivity
Ss = np.zeros((NROW, NCOL, NLAY),float); # specific storage
Sy = np.zeros((NROW, NCOL, 1),float); # specific yield
ctr = 7+1
for ii in range(NLAY):
    HY[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
    max_rows=NROW, dtype=float)
    ctr = ctr + NROW + 1
    VKA[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
    max_rows=NROW, dtype=float)
    ctr = ctr + NROW + 1
    Ss[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
    max_rows=NROW, dtype=float)
    ctr = ctr + NROW + 1
    if ii == 0:
        Sy[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
        max_rows=NROW, dtype=float)
        ctr = ctr + NROW + 1


# -- plot IBOUND [BA6] for active cells
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(2,2,1)
ncmap = np.max(IBOUND) - np.min(IBOUND) + 1
cmap = plt.get_cmap('Set1',ncmap)
im = ax.imshow(IBOUND, interpolation='none',cmap=cmap)
im.set_clim(np.min(IBOUND)-0.5, np.max(IBOUND)+0.5)
fig.colorbar(im, orientation='horizontal')
plt.title('IBOUND (active cells)')


# -- plot domain discretization [DIS]
TOP_to_plot = TOP.copy()
TOP_to_plot[TOP <= 0] == np.nan
BOTM_to_plot = BOTM.copy()
BOTM_to_plot[BOTM <= 0] == np.nan

min_to_plot = np.min(BOTM_to_plot)
max_to_plot = np.max(TOP_to_plot)

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(2,2,1)
x = np.copy(TOP_to_plot)
x2 = x[1:-1,1:-1]
x2[ind_bound_out] = 0
x[1:-1,1:-1] = x2
im = ax.imshow(x, interpolation='none')
im.set_clim(min_to_plot, max_to_plot)
im.set_cmap(plt.cm.terrain)
# use im.get_clim() to get the limits and generalize
fig.colorbar(im, orientation='horizontal')
plt.title('TOP')
for ilay in range(NLAY):
    x = np.copy(BOTM_to_plot[:,:,ilay])
    x2 = x[1:-1,1:-1]
    x2[ind_bound_out] = 0
    x[1:-1,1:-1] = x2
    plt.subplot(2,2,2+ilay)
    im = plt.imshow(x, interpolation='none')
    im.set_clim(min_to_plot, max_to_plot)
    fig.colorbar(im, orientation='horizontal')
    plt.set_cmap(plt.cm.terrain)
    plt.title('BOTM lay' + str(ilay+1));
    

# -- Hydraulic conductivity [UPW]   
for ilay in range(NLAY):
    x = np.copy(HY[:,:,ilay])
    x2 = x[1:-1,1:-1]
    x2[ind_bound_out] = 0
    x[1:-1,1:-1] = x2   
    fig = plt.figure(figsize=(12,12))
    plt.subplot(2,2,ilay+1)
    im = plt.imshow(x, interpolation='none')
    plt.set_cmap(plt.cm.cool)
#    im.set_clim(3800, 6200)
    fig.colorbar(im, orientation='horizontal')
    plt.title('HydCond in lay' + str(ilay+1));
    
    