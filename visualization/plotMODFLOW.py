# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 23:36:53 2017

@author: gcng
@author: awickert
"""
import sys
import platform
import struct
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as manimation
import matplotlib as mpl
import argparse

###############
## FUNCTIONS ##
###############


# DUMMY CLASS TO DECIDE ON MOVIE MAKING 
#########################################

# https://stackoverflow.com/questions/22226708/can-a-with-statement-be-used-conditionally

class Dummysink(object):
    def write(self, data):
        pass # ignore the data
    def __enter__(self): return self
    def __exit__(*x): pass

def datasink(writer, fig, moviefile_name=None):
    if moviefile_name:
        return writer.saving(fig, moviefile_name, 100)
    else:
        return Dummysink()

# function for parsing ASCII grid header in GIS data files
############################################################

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
    
# Truncate colormap
#####################

# https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


#########################
## COMMAND-LINE PARSER ##
#########################

parser = argparse.ArgumentParser(description= \
        'Plot animated time-series of groundwater conditions from GSFLOW.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

requiredArgs = parser.add_argument_group('required arguments')

# REQUIRED
requiredArgs.add_argument('-i', '--infile', type=str, default=argparse.SUPPRESS,
                    help='input <settings>.ini fiile for GSFLOW',
                    required = True)
requiredArgs.add_argument('-p', '--plot', type=str, default='wtd',
                          choices=['topo', 'head', 'wtd', 'dhead', 'hydcond',
                                   'hydcond_vertical', 'ss', 'sy'],
                          help='Plot variable selector: topography (topo), \
                                head, water table depth (wtd), change in head \
                                (dhead), hydraulic conductivity (hydcond), \
                                vertical hydraulic conductivity (hydcond_vert) \
                                specific storage (ss), or specific yield (sy)')

# OPTIONAL
parser.add_argument('-o', '--outmovie', type=str, default=None,
                    help='Output file (mp4) for movie, if desired')

args = parser.parse_args()
args = vars(args)

settings_input_file = args['infile']
plotvar = args['plot']
moviefile_name = args['outmovie']

##################
## MAIN PROGRAM ##
##################

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

# add path containing readSettings.py
sys.path.append('..' + slashstr + 'Run')

# Read in user-specified settings
from readSettings import Settings
Settings = Settings(settings_input_file)

# Entiries from Settings File
###############################

head_file = Settings.MODFLOWoutput_dir + slashstr + Settings.PROJ_CODE + '_head.bhd'
surfz_fil = Settings.GISinput_dir + slashstr + 'DEM.asc'
ba6_fil = Settings.MODFLOWinput_dir + slashstr + Settings.PROJ_CODE + '.ba6'
dis_fil = Settings.MODFLOWinput_dir + slashstr + Settings.PROJ_CODE + '.dis'
flo_fil = Settings.MODFLOWinput_dir + slashstr + Settings.PROJ_CODE + '.upw'


print '\n******************************************'
print 'Plotting results from: ', head_file
print ' (WTD=TOP-HEAD calculated using topo data in: ' + surfz_fil + ')'
print ' (active cells info in: ' + ba6_fil + ')'
print '******************************************\n'

# Information on domain
#########################
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

sdata = read_grid_file_header(surfz_fil)
NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]

#NROW = sdata['rows'] 
#NCOL = sdata['cols']

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


# - space discretization
_N, _S, _E, _W = NSEW[0], NSEW[1], NSEW[2], NSEW[3]

dx = DELR = (NSEW[2]-NSEW[3])/NCOL # width of column [m]
dy = DELC = (NSEW[0]-NSEW[1])/NROW # height of row [m]

x = np.arange(_W + dx/2., _E, dx)
y = np.arange(_S + dy/2., _N, dy)
X, Y = np.meshgrid(x,y)


# VARIABLES
#############

HY = np.zeros((NROW, NCOL, NLAY),float); # hyd conductivity
VKA = np.zeros((NROW, NCOL, NLAY),float); # vertical hyd conductivity
Ss = np.zeros((NROW, NCOL, NLAY),float); # specific storage
Sy = np.zeros((NROW, NCOL, 1),float); # specific yield
ctr = 7+1
for ii in range(NLAY):
    HY[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
    max_rows=NROW, dtype=float)

# Head
data_head_all_NaN = data_head_all
data_head_all_NaN[data_head_all_NaN > 1e29] = np.nan # dry cell
data_head_all_NaN[data_head_all_NaN <= -999] = np.nan

# WTD:
WTD_all = np.tile(TOP[:,:,np.newaxis], (1,1,ntimeslay)) - data_head_all_NaN

# change in head (i - (i-1)): 0's at first time
dhead_all = np.zeros((NROW,NCOL,ntimeslay))
dhead_all[:,:,1:] = data_head_all_NaN[:,:,1:] - data_head_all_NaN[:,:,:-1]

# Active cells
IBOUND = np.genfromtxt(ba6_fil, skip_header=3, max_rows=NROW, dtype=float)

# Topography in basin only
TOP_in_basin = TOP * (IBOUND == 1)
TOP_in_basin[TOP_in_basin == 0] = np.nan

# Hydraulic conductivity, specific storage, and specific yield
hydraulic_conductivity = np.zeros((NROW, NCOL, NLAY), float)
hydraulic_conductivity__vertical = np.zeros((NROW, NCOL, NLAY), float)
specific_storage = np.zeros((NROW, NCOL, NLAY), float)
specific_yield = np.zeros((NROW, NCOL, 1), float)
_ctr = 7+1
for ii in range(NLAY):
    hydraulic_conductivity[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
        max_rows=NROW, dtype=float)
    _ctr = _ctr + NROW + 1
    hydraulic_conductivity__vertical[:,:,ii] = np.genfromtxt(flo_fil, \
        skip_header=ctr, max_rows=NROW, dtype=float)
    _ctr = _ctr + NROW + 1
    specific_yield[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
        max_rows=NROW, dtype=float)
    _ctr = _ctr + NROW + 1
    if ii == 0:
        specific_yield[:,:,ii] = np.genfromtxt(flo_fil, skip_header=ctr, \
            max_rows=NROW, dtype=float)
        _ctr = _ctr + NROW + 1

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

# Axis labels without offsets
y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
x_formatter = mpl.ticker.ScalarFormatter(useOffset=False)

# Plot extent
_extent = (_W/1000., _E/1000., _S/1000., _N/1000.)
_extent_countour = (_W/1000., _E/1000., _N/1000., _S/1000.)

# Number of panels
if NLAY < 2:
    nrows = 1
else:
    nrows = 2
ncols = (NLAY+1)/2

# Multiple times or no
if plotvar in ['topo', 'hydcond', 'hydcond', 'ss', 'sy']:
    static_plot = True
else:
    static_plot = False

# Plot
fig = plt.figure()

if moviefile_name:
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='GSFLOW Movie', artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps=10, metadata=metadata)
else:
    FFMpegWriter = None
    metadata = None
    writer = None

if not static_plot:
    with datasink(writer=writer, fig=fig, moviefile_name=moviefile_name):
        ctr = 0
        for ii in range(ntimes):
            for lay_i in range(NLAY):
                if plotvar == 'head':
                    # head:
                    cbl = 'Hydraulic head [m]'
                    #ti = 'head [m], '
                    data_all = data_head_all_NaN
                elif plotvar == 'wtd':        
                    # WTD:
                    cbl = 'Water table depth [m]'
                    data_all = WTD_all
                elif plotvar == 'dhead':
                    # change in head:
                    cbl = 'Change in hydraulic head [m]'
                    data_all = dhead_all

                data = data_all[:,:,ctr]    
                
                if ii == 0:
                    print ii
                    if lay_i == 0:
                        av = []
                        pv = []
                        cv = []
                    av.append(plt.subplot(nrows, ncols, lay_info[0,ctr]))
                    pv.append(av[lay_i].imshow(data, interpolation='nearest', 
                                               extent=_extent))
                    pv[lay_i].set_cmap(plt.cm.cool)
                    cv.append(plt.colorbar(pv[lay_i]))
                    _col = data_all[:,:,lay_i::2]
                    _col = _col[~np.isnan(_col)]
                    cv[lay_i].set_label(cbl, fontsize=24)
                    pv[lay_i].set_clim(vmin=np.min(_col), vmax=np.max(_col))
                    av[lay_i].set_xlabel('E [km]', fontsize=24)
                    av[lay_i].set_ylabel('N [km]', fontsize=24)
                    av[lay_i].yaxis.set_major_formatter(y_formatter)
                    av[lay_i].xaxis.set_major_formatter(x_formatter)
                    cs = av[lay_i].contour(TOP_in_basin, colors='k', 
                                           extent=_extent_countour)
                    plt.clabel(cs, inline=1, fontsize=16, fmt='%d')
                    av[lay_i].set_aspect('equal', 'datalim')
                else:
                    pv[lay_i].set_data(data)        
                titlestr = '%d' %time_info[0,ctr] + ' days; layer ' + \
                               str(int(lay_info[0,ctr])) + \
                               '\nwith topographic contours [m]'
                av[lay_i].set_title(titlestr, fontsize=24)
                im2 = av[lay_i].imshow(outline, interpolation='nearest',
                                       extent=_extent)
                im2.set_clim(0, 1)
                cmap = plt.get_cmap('binary',2)
                im2.set_cmap(cmap)   

                ctr = ctr + 1
            #    plt.show()
            #plt.tight_layout()
            plt.pause(0.5)
            
            # Optional write movie frame
            if moviefile_name:
                writer.grab_frame()
                
            # Optional write figure to file
            
            #for _axis in av:
            #    _axis.cla()
                    
            #plt.savefig("myplot.png", dpi = 300)
else:
    av = []
    pv = []
    cv = []
    if plotvar == 'topo':
        cbl = 'Topographic elevation [m]'
        data_all = data_head_all_NaN
    elif plotvar == 'hydcond':
        cbl = 'Hydraulic conductivity [m/d]'
        data_all = hydraulic_conductivity
    elif plotvar == 'hydcond_vert':
        cbl = 'Hydraulic conductivity (vertical) [m/d]'
        data_all = hydraulic_conductivity__vertical
    elif plotvar == 'ss':
        cbl = 'Specific storage [1/m]'
        data_all = specific_storage
    elif plotvar == 'sy':
        cbl = 'Specific yield [-]'
        data_all = specific_yield
    lay_i = 0
    data = data_all = TOP_in_basin
    av.append(plt.subplot(nrows, ncols, 1))
    pv.append(av[lay_i].imshow(data, interpolation='nearest', 
                               extent=_extent))
    if plotvar == 'topo':
        pv[lay_i].set_cmap(truncate_colormap(plt.cm.terrain, 0.25, 1.0))
    else:
        pv[lay_i].set_cmap(plt.cm.cool)
    cv.append(plt.colorbar(pv[lay_i]))
    _col = data_all[:]
    _col = _col[~np.isnan(_col)]
    cv[lay_i].set_label(cbl, fontsize=24)
    pv[lay_i].set_clim(vmin=np.min(_col), vmax=np.max(_col))
    av[lay_i].set_xlabel('E [km]', fontsize=24)
    av[lay_i].set_ylabel('N [km]', fontsize=24)
    av[lay_i].yaxis.set_major_formatter(y_formatter)
    av[lay_i].xaxis.set_major_formatter(x_formatter)
    av[lay_i].set_aspect('equal', 'datalim')
    im2 = av[lay_i].imshow(outline, interpolation='nearest',
                           extent=_extent)
    im2.set_clim(0, 1)
    cmap = plt.get_cmap('binary',2)
    im2.set_cmap(cmap)   

