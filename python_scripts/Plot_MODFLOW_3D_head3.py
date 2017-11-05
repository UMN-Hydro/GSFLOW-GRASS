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
from readSettings import Settings

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
else:
    slashstr = '\\'


# ***Select which option:
# (1: head, 2: water table depth, 3: change in head per print-time increment)
sw_head_WTD_dhead = 2


#%% from Settings file
head_file = Settings.MODFLOWoutput_dir + slashstr + Settings.PROJ_CODE + '_head.bhd'  # head data
surfz_fil = Settings.GISinput_dir + slashstr + 'DEM.asc'
ba6_fil = Settings.MODFLOWinput_dir + slashstr + Settings.PROJ_CODE + '.ba6'

print '\n******************************************'
print 'Plotting results from: ', head_file
print ' (WTD=TOP-HEAD calculated using topo data in: ' + surfz_fil + ')'
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

TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)


# =========================================================================

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'


fl_binary = 1;  # 1 for binary
fl_dble = 0;  # 1 for dble prec, 0 for single prec

dry_cell = 1e30;
inactive_cell = -999.99;


# Save precision to variables;
if fl_dble:
    prec = 8
else:
    prec = 4


# -- Get head data and plot it as contour image plots
fid = open(head_file, 'rb')

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

data_head_all = np.zeros([NROW,NCOL,0])
time_info = np.zeros([4,0])
lay_info = np.zeros([1,0])
ii = 0
while True:
    # NOTE: for some reason, int w/ bit-length info is trailed by 0!!!
    a_info = binbuild( nitems=nread, nbytes=prec, typecode='i', infile=fid )
    kstp = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    if not kstp:
        break
    kper = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    pertim = binbuild(nitems=1, nbytes=prec, typecode='f', infile=fid ) # FIX THESE TO SCALAR
    totim = binbuild(nitems=1, nbytes=prec, typecode='f', infile=fid ) # FIX THESE TO SCALAR
    label = binbuild(nitems=16, nbytes=1, typecode='c', infile=fid ) # FIX TO CHAR ARRAY
    ncol = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    nrow = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    ilay = binbuild(nitems=1, nbytes=prec, typecode='i', infile=fid ) # FIX THESE TO SCALAR
    a_info = binbuild(nitems=nread, nbytes=prec, typecode='i', infile=fid )
    
    a_data = binbuild(nitems=nread, nbytes=prec, typecode='i', infile=fid )
    if nread == 0:
        nn = ncol*nrow
    else:
        nn = a_data/prec # is floor divide OK? Also, shouldn't it just be nlay?
    
    data = binbuild(nitems=nn, nbytes=prec, typecode='f', infile=fid)
    a_data = binbuild(nitems=nread, nbytes=prec, typecode='i', infile=fid )
    
    if ii % 100 == 0:  # mod 100
        data_head_all2 = np.zeros([NROW,NCOL,ii+100])
        data_head_all2[:,:,:ii] = data_head_all
        data_head_all = data_head_all2
        time_info2 = np.zeros([4,ii+100])
        time_info2[:,:ii] = time_info
        time_info = time_info2
        lay_info2 = np.zeros([1,ii+100],int)
        lay_info2[:,:ii] = lay_info
        lay_info = lay_info2        
    
    data_head_all[:,:,ii] = np.reshape(data, (nrow,ncol), order='C') 
    time_info[:,ii] = [kstp, kper, pertim, totim]
    lay_info[0,ii] = int(ilay)
    
    ii = ii + 1

fid.close()    

data_head_all = data_head_all[:,:,:ii]
time_info = time_info[:,:ii]
lay_info = lay_info[:,:ii]
ntimeslay = ii # ntimes x nlay

NLAY = np.max(lay_info)
ntimes = ntimeslay / NLAY

x = np.arange(DELR/2., DELR*NCOL+DELR/2., DELR)
y = np.arange(DELC/2., DELC*NROW+DELC/2., DELC)
X, Y = np.meshgrid(x,y)


data_head_all_NaN = data_head_all
data_head_all_NaN[data_head_all_NaN > 1e29] = np.nan # dry cell
data_head_all_NaN[data_head_all_NaN <= -999] = np.nan


# use this to plot WTD:
TOP2 = np.tile(TOP[:,:,np.newaxis], (1,1,ntimeslay))
WTD_all = TOP2 - data_head_all_NaN

# use this to plot change in head:        
dhead_all = np.zeros((NROW,NCOL,ntimeslay))
dhead_all[:,:,1:] = data_head_all_NaN[:,:,1:] - data_head_all_NaN[:,:,:-1]

# -- get active cells
IBOUND = np.genfromtxt(ba6_fil, skip_header=3, max_rows=NROW, dtype=float)

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


# mask with just outline (outer boundary) of watershed
outline = np.ones((NROW,NCOL))*np.nan
outline2 = outline[1:-1,1:-1]
outline2[ind_bound_out] = 1
outline[1:-1,1:-1] = outline2

# head plot movie
fig = plt.figure()
ctr = 0
for ii in range(ntimes):
    for lay_i in range(NLAY):
        
        if sw_head_WTD_dhead == 1:
            # head:
            ti = 'head [m], '
            data_all = data_head_all_NaN
        elif sw_head_WTD_dhead == 2:        
            # WTD:
            ti = 'WTD=TOP-HEAD [m], '
            data_all = WTD_all
        elif sw_head_WTD_dhead == 3:
            # change in head:
            ti = 'Change in head [m], '
            data_all = dhead_all
               
#        data = data_all[:,:,ii]    
        data = data_all[:,:,ctr]    
        
        if ii == 0:
            print ii
            if lay_i == 0:
                av = []
                pv = []
                
            av.append(plt.subplot(2,2,lay_info[0,ctr]))
            pv.append(av[lay_i].imshow(data, interpolation='nearest'))
            pv[lay_i].set_cmap(plt.cm.cool)
            plt.colorbar(pv[lay_i])
#            plt.clim()
            x = data_all[:,:,lay_i::2]
            x = x[~np.isnan(x)]
            pv[lay_i].set_clim(vmin=np.min(x), vmax=np.max(x))
        
#            plt.subplot(2,2,lay_info[0,ii])
#            p = plt.imshow(data, extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto', interpolation='none')
#            p.set_cmap(plt.cm.hsv)
#            plt.colorbar(p)
##            plt.clim()
#            x = data_all[~np.isnan(data_all)]
#            p.set_clim(vmin=np.min(x), vmax=np.max(x))
#            plt.xlabel('[m]', fontsize=16)
#            plt.ylabel('[m]', fontsize=16)
            av[lay_i].set_xlabel('[m]', fontsize=16)
            av[lay_i].set_ylabel('[m]', fontsize=16)
        else:
#            p.set_data(data)        
            pv[lay_i].set_data(data)        
        str0 = ti + str(time_info[0,ctr]) + 'd, lay ' + str(int(lay_info[0,ctr]))
#            plt.title(str0)
        av[lay_i].set_title(str0)
        plt.tight_layout()

        im2 = av[lay_i].imshow(outline, interpolation='nearest')
#        im2 = plt.imshow(outline, interpolation='none')
        im2.set_clim(0, 1)
        cmap = plt.get_cmap('binary',2)
        im2.set_cmap(cmap)   
                      
        ctr = ctr + 1
    #    plt.show()
    plt.pause(0.5)
            
    #plt.savefig("myplot.png", dpi = 300)


