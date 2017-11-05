# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 23:36:53 2017

@author: gcng
"""
import sys
import platform
import struct
import numpy as np
from matplotlib import pyplot as plt
import settings_test

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'


# ***Select which variable to plot:
# 0: 'UZF RECHARGE'
# 1: 'SURFACE LEAKAGE' (negative)
# 2: 'HORT+DUN'
# 3: 'SURFACE LEAKAGE' (positive)
# 4: 'STORAGE CHANGE'
sw_PlotVar = 0

# ***Optional: set plot variable limits, set to empty list for default of min / max values
clim = [0, 250]
#clim = [] # empty list for default min, max plot limits


#%% from Settings file 

# *** Change file names as needed
uzf_file = settings_test.MODFLOWoutput_dir + slashstr + settings_test.PROJ_CODE + '_uzf.dat'  # head data
surfz_fil = settings_test.GISinput_dir + slashstr + 'DEM.asc'
ba6_fil = settings_test.MODFLOWinput_dir + slashstr + settings_test.PROJ_CODE + '.ba6'

print '\n******************************************'
print 'Plotting results from: ', uzf_file
print ' (domain data in: ' + surfz_fil + ')'
print ' (active cells info in: ' + ba6_fil + ')'
print '******************************************\n'

#%% In general: don't change below here

# Only ONE can be 1, others 0
if sys.platform[:3] == 'win':
    nread = 0
elif (platform.linux_distribution()[0] == 'Ubuntu') or (platform.linux_distribution()[0] == 'debian'):
    nread = 1
elif platform.linux_distribution()[0][:3] == 'Red':
    # Hope this works; haven't tried Red Hat here
    nread = 2
else:
    sys.exit("You should add your OS binary formatting to this script!")

#uzf_file = settings_test.MODFLOWoutput_dir + slashstr + 'Shullcas_uzf_win.dat'  # head data
#nread = 0

# -- get surface elevations [m] (to plot WTD)
# function for parsing ASCII grid header in GIS data files
def read_grid_file_header(fname):
    f = open(fname, 'r')
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

    return sdata
    
sdata = read_grid_file_header(surfz_fil)
    
NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
NROW = sdata['rows'] 
NCOL = sdata['cols']

# - space discretization
DELR = (NSEW[2]-NSEW[3])/NCOL # width of column [m]
DELC = (NSEW[0]-NSEW[1])/NROW # height of row [m]


# -- get active cells
IBOUND = np.genfromtxt(ba6_fil, skip_header=3, max_rows=NROW, dtype=float)
# =========================================================================

# -- find boundary cells
IBOUND0 = np.copy(IBOUND)
IBOUND0.astype(int)
IBOUND0[IBOUND0>0] = 1 # active cells
IBOUND0[IBOUND0<0] = 1 # constant head cells
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



fl_binary = 1;  # 1 for binary
fl_dble = 0;  # 1 for dble prec, 0 for single prec


# Save precision to variables;
if fl_dble:
    prec = 8
else:
    prec = 4


# -- Get head data and plot it as contour image plots
fid = open(uzf_file, 'rb')

#fid.read(4)
#struct.unpack('i', fid.read(4)) 
#fid.close()

def binbuild(nitems, nbytes, typecode, infile):
    outdata = []
    for i in range(nitems):
        x = infile.read(nbytes)
        if len(x) == 0:
            return ''
        outdata.append( struct.unpack(typecode, x) )
    return np.squeeze(np.array(outdata))

time_info = np.zeros([2,0])
all_label = []
nvar = 5;
ii = 0
t_i = 0
while True:
    # NOTE: for some reason, int w/ bit-length info is trailed by 0!!!
    a_info = binbuild( nitems=nread, nbytes=prec, typecode='i', infile=fid )
    kstp = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    if not kstp:
        break
    kper = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
#    pertim = binbuild(nitems=1, nbytes=prec, typecode='f', infile=fid ) # FIX THESE TO SCALAR
#    totim = binbuild(nitems=1, nbytes=prec, typecode='f', infile=fid ) # FIX THESE TO SCALAR
    label = binbuild(nitems=16, nbytes=1, typecode='c', infile=fid ) # FIX TO CHAR ARRAY
    ncol = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    nrow = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    ilay = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    a_info = binbuild(nitems=nread, nbytes=prec, typecode='i', infile=fid )
    
    a_data = binbuild(nitems=nread, nbytes=prec, typecode='i', infile=fid )
    nn = ncol*nrow*ilay
    
#    if nread == 0:
#        nn = ncol*nrow*ilay
#    else:
#        nn = a_data/prec # is floor divide OK? Also, shouldn't it just be nlay?
    
    data = binbuild(nitems=nn, nbytes=prec, typecode='f', infile=fid)
    a_data = binbuild(nitems=nread, nbytes=prec, typecode='i', infile=fid )

#    print data.shape
        
#    if nread == 0:
#        all_data = np.reshape(data, (nrow,ncol), order='C') 
#    else:
#        all_data = np.reshape(data, (ncol,nrow,ilay), order='F') 
#        all_data = np.transpose(all_data, axes=(1,0,2)) 
    all_data = np.reshape(data, (ncol,nrow,ilay), order='F') 
    all_data = np.transpose(all_data, axes=(1,0,2)) 
#        print 'here1'
#        all_data = np.reshape(data[nrow,ncol,0], (nrow,ncol), order='C') 
#        print 'here'

#        all_data = reshape(data,ncol,nrow,ilay);
#        all_data = permute(all_data, [2 1 3]);

    if ii == 0:
#        if nread == 0:
#            all_data_all = np.zeros([NROW,NCOL,nvar,0])
#        elif nread == 1:
        all_data_all = np.zeros([NROW,NCOL,ilay,nvar,0])
        
    var_i = ii % nvar;
    if ii % 100 == 0:  # mod 100
#        if nread == 0:
#            all_data_all2 = np.zeros([NROW,NCOL,nvar,ii+100])
#            all_data_all2[:,:,:,:ii] = all_data_all
#            all_data_all = all_data_all2
#        elif nread == 1:
        all_data_all2 = np.zeros([NROW,NCOL,ilay,nvar,ii+100])
        all_data_all2[:,:,:,:,:ii] = all_data_all
        all_data_all = all_data_all2
            
        time_info2 = np.zeros([2,ii+100])
        time_info2[:,:ii] = time_info
        time_info = time_info2
    
#    if nread == 0:
#        all_data_all[:,:,var_i,t_i] = all_data
#    elif nread == 1:
    all_data_all[:,:,:,var_i,t_i] = all_data
    time_info[:,t_i] = [kstp, kper]
    if ii < 5:
        all_label.append(str.strip(''.join(label)))
    
    if (ii % nvar) == 0: 
        t_i = t_i + 1
    ii = ii + 1

fid.close()    

#if nread == 0:
#    all_data_all = all_data_all[:,:,:,:t_i]
#elif nread == 1:
all_data_all = all_data_all[:,:,:,:,:t_i]
time_info = time_info[:,:t_i]
ntimes = t_i

#if nread == 0:
#    data_all = all_data_all[:,:,sw_PlotVar,:]
#elif nread == 1:
data_all = all_data_all[:,:,:,sw_PlotVar,:]
#    dd = all_data_all[:,:,1,sw_PlotVar,:] - all_data_all[:,:,0,sw_PlotVar,:]
    

# -- Plot
x = np.arange(DELR/2., DELR*NCOL+DELR/2., DELR)
y = np.arange(DELC/2., DELC*NROW+DELC/2., DELC)
X, Y = np.meshgrid(x,y)


#data_head_all_NaN = data_head_all
#data_head_all_NaN[data_head_all_NaN > 1e29] = np.nan
#data_head_all_NaN[data_head_all_NaN <= 999] = np.nan

# mask with just outline (outer boundary) of watershed
outline = np.ones((NROW,NCOL))*np.nan
outline2 = outline[1:-1,1:-1]
outline2[ind_bound_out] = 1
outline[1:-1,1:-1] = outline2

lay_i0 = 0

# head plot movie
ti = all_label[sw_PlotVar]
plt.figure()
for ii in range(ntimes):
#for ii in range(2):
    for lay_i in [lay_i0]:
#    for lay_i in range(ilay):
        
#        if nread == 0:              
#            data = data_all[:,:,ii]    
#        elif nread == 1:
        data = data_all[:,:,lay_i,ii]   
#            data = dd[:,:,ii]   
        
        # only show active domain, with outline around it
        data[IBOUND == 0] = np.nan
        data2 = data[1:-1,1:-1]
        data2[ind_bound_out] = 0
        data[1:-1,1:-1] = data2
        
        if ii == 0:
            plt.subplot(2,2,1)
#            p = plt.imshow(data, extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto', interpolation='none')
            p = plt.imshow(data, interpolation='none')
            p.set_cmap(plt.cm.rainbow)
            plt.colorbar(p)
#            plt.clim()
            x = data_all[~np.isnan(data_all)]
#            if len(clim) == 0:
#                p.set_clim(vmin=np.min(x), vmax=np.max(x))
#            else:
#                p.set_clim(vmin=clim[0], vmax=clim[1])
            try:
                clim
                p.set_clim(vmin=clim[0], vmax=clim[1])
                del clim                
            except NameError:
                p.set_clim(vmin=np.min(x), vmax=np.max(x))
            plt.xlabel('[m]', fontsize=16)
            plt.ylabel('[m]', fontsize=16)
        else:
            p.set_data(data)        
            str0 = ti + ' ' + str(int(time_info[0,ii])) 
            plt.title(str0)
        plt.tight_layout()
        im2 = plt.imshow(outline, interpolation='none')
#        im2 = plt.imshow(outline, interpolation='none')
        im2.set_clim(0, 1)
        cmap = plt.get_cmap('binary',2)
        im2.set_cmap(cmap)   
           
    #    plt.show()
        plt.pause(0.5)
            
    #plt.savefig("myplot.png", dpi = 300)

#try:
#    clim
#    clim = None # set for next run
#except:
#    

