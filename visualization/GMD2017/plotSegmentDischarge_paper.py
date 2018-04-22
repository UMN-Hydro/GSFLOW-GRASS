#! /usr/bin/env python

from osgeo import ogr
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.animation as manimation
import platform
import sys
import re

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

# add path containing readSettings.py
sys.path.append('..' + slashstr + 'Run')

# Read in user-specified settings
from readSettings import Settings
# Set input file
if len(sys.argv) < 2:
    settings_input_file = 'settings.ini'
    print 'Using default input file: ' + settings_input_file
else:
    settings_input_file = sys.argv[1]
    print 'Using specified input file: ' + settings_input_file
Settings = Settings(settings_input_file)

#%% *** SET THE FOLLOWING *****************************************************

# *** Save movie to following file
moviefile_name = 'SR_strmseg.mp4'


###########################################
## Specify time to plot (instead of movie) 
###########################################
# default is all, for movie
ptime_ind = [] # movie
p_dates = []


# - only set ONE of the below (ptime_ind OR p_dates)
#ptime_ind = [-1] # plot only this time index, starts at 0 (-1 for last, but last entry isn't read in properly, so use -2 instead)
#p_dates = ['1993-01-09']

## - Shullcas 
## run plotSegmentDischarge_paper.py /home/gcng/workspace/ProjectFiles/GSFLOW-GRASS_ms/examples4ms/Shullcas_gcng.ini
### run plotSegmentDischarge_paper.py /media/gcng/STORAGE3A/ANDY/GSFLOW/Shullcas_gcng.ini
#figName = 'Shullcas_Q'
#p_dates = ['2015-02-15']
##figsize0 = (7.5,5.5) # default (8W,6H) [inches]
#plot_pos = (0,1) # row 0, col 1 
#xlim = [480, 498]
#ylim = [8665, 8689]
##plot_ti_ltr = 'C) '


# - Santa Rosa 
# run plotSegmentDischarge_paper.py /home/gcng/workspace/ProjectFiles/GSFLOW-GRASS_ms/examples4ms/SantaRosa_WaterCanyon_gcng.ini
## run plotSegmentDischarge_paper.py /media/gcng/STORAGE3A/ANDY/GSFLOW/SantaRosa_WaterCanyon_gcng.ini
figName = 'SR_Q'
p_dates = ['2017-02-17']
plot_pos = (0,1) # row 0, col 1 
xlim = [213, 220]
ylim = [3760, 3766]
#plot_ti_ltr = 'C) '

## - Cannon River - 2 layer
## run plotSegmentDischarge_paper.py /media/gcng/STORAGE3A/ANDY/GSFLOW/CannonRiver_2layer_gcng.ini
#figName = 'Cannon_Q'
#p_dates = ['1943-07-06']

# font sizes
FS_lab = 10
FS_cvtick = 8
FS_xylab = 10
FS_clab = 10
FS_ti = 10



#%% *** CHANGE FILE NAMES AS NEEDED *******************************************
# (default is to use entries from Settings File) 
segout_fil = Settings.PRMSoutput_dir + slashstr + Settings.PROJ_CODE + '.ani.nsegment'
projdir_GIS = Settings.GISinput_dir
segshp_fil = projdir_GIS + slashstr + "shapefiles/segments/segments.shp"

print '\n**********************************************************************'
print 'Plotting results from:'
print segout_fil
print '  (segments shapefile:'
print '  ', segshp_fil + ')'
print '************************************************************************'
print

#%% In general: don't change below here

plotting_variable = 'streamflow_sfr'

#for filename in [HRUout_fil, segment_filename]:
#for filename in [segout_fil]:
#    infile = file(filename, 'r')
#    outfile = file(filename + '.corrected', 'w')
#    for line in infile:
#        if line[:2] == '  ':
#            p = re.compile("\-[0-9]{2}\-")
#            for m in p.finditer(line):
#                if m.start():
#                    break
#            _start = m.start() - 4 # space for the year
#            outfile.write(line[_start:])
#        else:
#            outfile.write(line)
            
# Correct the file with the streamflow data
for filename in [segout_fil]:
    infile = file(filename, 'r')
    outfile = file(filename + '.corrected', 'w')
    for line in infile:
        if line[0] == ' ':
            p = re.compile("\-[0-9]{2}\-")
            for m in p.finditer(line):
                if m.start():
                    break
            _start = m.start() - 4 # space for the year
            outfile.write(line[_start:])
        else:
            outfile.write(line)              
outfile.close()

segout_fil2 = segout_fil + '.corrected'

_shapefile = ogr.Open(segshp_fil)
_shape = _shapefile.GetLayer(0)

segment_outputs = pd.read_csv(segout_fil2, comment='#', delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False, skiprows=[8])

dates = sorted(list(set(list(segment_outputs.timestamp))))

fl_plot_single = 1
if (len(ptime_ind) == 0) and (len(p_dates) == 0):
    p_dates0 = dates
    fl_plot_single = 0
elif len(ptime_ind) > 0:
    p_dates0 = [dates[ptime_ind[0]]]
elif len(p_dates) > 0:
    p_dates0 = p_dates


plotting_variable = 'streamflow_sfr'

# for color axis
if fl_plot_single == 1: # max only for that date
    _segment_outputs_on_date = segment_outputs.loc[segment_outputs['timestamp'] == p_dates0[0]] 
    data = _segment_outputs_on_date[plotting_variable].values* 0.0283168466    
else:
    data = segment_outputs[plotting_variable] * 0.0283168466
#data = segment_outputs[plotting_variable] * 0.0283168466
_min = np.min(data[data>0])
_max = np.max(data)
if fl_plot_single == 1:
    dd = np.log10(_max) - np.log10(_min)
    dd2 = dd*1.1
    lmax = np.log10(_min) + dd2
    _max = 10.**lmax


nrows = 2
ncols = 3
fig = plt.figure(1)
#gridspec.GridSpec(nrows,ncols)
#plt.ion()

#ax = plt.subplot(111)
ax = plt.subplot2grid((nrows, ncols), plot_pos)


cmap = plt.get_cmap('RdYlBu')
_min = np.max((_max * 1e-2, _min))
print "colorbar min, max: ", _min, _max
cax, _ = mpl.colorbar.make_axes(ax, location='right')
cbar = mpl.colorbar.ColorbarBase(cax, cmap=cm.jet,
               norm=mpl.colors.LogNorm(vmin=_min, vmax=_max))


y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
x_formatter = mpl.ticker.ScalarFormatter(useOffset=False)

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=10, metadata=metadata)

print "total number of dates", len(dates)    

#print segment_outputs
#with writer.saving(fig, moviefile_name, 100):
for date in p_dates0:
#        print date
    _segment_outputs_on_date = segment_outputs.loc[segment_outputs['timestamp'] == date]
    _values = []
    for i in range(_shape.GetFeatureCount()):
        _feature = _shape.GetFeature(i)
        _n = _feature['id']
        #print _nhru
        _row = _segment_outputs_on_date.loc[_segment_outputs_on_date['nsegment'] == _n]
        try:
            # cfs to m3/s
            _values.append(float(_row[plotting_variable].values))
        except:
            _values.append(np.nan)
            print _n
            continue
    _values = np.array(_values) * 0.0283168466
    # Floating colorbar
    colors = cm.jet(plt.Normalize( np.log10(_min), np.log10(_max)) 
                                   (np.log10(_values)) )
    ax.cla()
    _lines = []
    for i in range(_shape.GetFeatureCount()):
        _feature = _shape.GetFeature(i)
        #feature = shape.GetFeature(0) # how to get it otherwise
        _geometry = _feature.geometry()
        _line_points = np.array(_geometry.GetLinearGeometry().GetPoints())
        _x = _line_points[:,0]/1000.
        _y = _line_points[:,1]/1000.
        scale = .03
        _lines.append( ax.plot(_x, _y, '-', color=colors[i], linewidth=(scale*_values[i]/0.0283168466)**.5+.25) )
    #ax.set_title(plotting_variable+': '+date)
#    ax.set_title(date, fontsize=FS_ti, fontweight='bold')
#    ax.set_title(date, fontsize=FS_ti)
    ax.set_title('Streamflow [m$^3$/s]', fontsize=FS_ti)
#    cbar.set_label(r'Streamflow [m$^3$/s]', fontsize=FS_lab)
    ax.set_xlabel('E [km]', fontsize=FS_xylab)
    ax.set_ylabel('N [km]', fontsize=FS_xylab)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.tick_params(axis='both', which='major', labelsize=FS_cvtick)
    cax.tick_params(labelsize=FS_cvtick) 
                        
#    ax.set_aspect('equal', 'datalim')
    
#    plt.tight_layout()
#    writer.grab_frame()
    plt.pause(0.01)
    #plt.waitforbuttonpress()

#ax.set_aspect('equal')
#xlim = ax.get_xlim()
#ax.set_xticks(np.arange(np.floor(xlim[0]), np.floor(xlim[-1]), np.floor((xlim[-1]-xlim[0])/3)))  

#ax.set_aspect('equal', 'datalim')
ax.set_aspect('equal')
ax.set_xlim(xlim)  
ax.set_ylim(ylim)  
ax.set_xticks(2+np.arange(np.floor(xlim[0]), np.floor(xlim[-1]), np.ceil((xlim[-1]-xlim[0])/3)))  
ax.set_yticks(2+np.arange(np.floor(ylim[0]), np.floor(ylim[-1]), np.ceil((ylim[-1]-ylim[0])/3)))  

#plt.savefig(figName+'.png', dpi=100)
#plt.savefig(figName+'.svg')