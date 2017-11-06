from osgeo import ogr
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.animation as manimation
from readSettings import Settings
import platform
import sys
import re

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

#%% User-specified settings here:

# *** Enter PRMS variable name to plot (must be specified under 'aniOutVar_names' in control file, 
# see GSFLOW manual Table A1-2 for list of output variable names)
plotting_variable = 'hru_ppt'

# *** Save movie to following file
moviefile_name = 'testmovie_HRUs.mp4'


#%% from Settings file, change to plot something else
HRUout_fil = Settings.PRMSoutput_dir + slashstr + Settings.PROJ_CODE + '.ani.nhru'
projdir_GIS = Settings.GISinput_dir
HRUshp_fil = projdir_GIS + slashstr + "shapefiles/HRUs/HRUs.shp"

print '\n******************************************'
print 'Plotting results from: ', HRUout_fil
print ' (HRU shapefile: ' + HRUshp_fil + ')'
print '******************************************\n'

#%% In general: don't change below here


for filename in [HRUout_fil]:
    infile = file(filename, 'r')
    outfile = file(filename + '.corrected', 'w')
    for line in infile:
        if line[:2] == '  ':
            p = re.compile("\-[0-9]{2}\-")
            for m in p.finditer(line):
                if m.start():
                    break
            _start = m.start() - 4 # space for the year
            outfile.write(line[_start:])
        else:
            outfile.write(line)

HRUout_fil2 = HRUout_fil + '.corrected'
HRU_outputs = pd.read_csv(HRUout_fil2, comment='#', delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False, skiprows=[11])

_shapefile = ogr.Open(HRUshp_fil)
_shape = _shapefile.GetLayer(0)

dates = sorted(list(set(list(HRU_outputs.timestamp))))


_min = np.min(HRU_outputs[plotting_variable])
_max = np.max(HRU_outputs[plotting_variable])

fig = plt.figure(figsize=(8,6))
#plt.ion()

ax = plt.subplot(111)

cax, _ = mpl.colorbar.make_axes(ax, location='right')
cbar = mpl.colorbar.ColorbarBase(cax, cmap=cm.jet,
               norm=mpl.colors.Normalize(vmin=_min, vmax=_max))

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=10, metadata=metadata)

with writer.saving(fig, moviefile_name, 100):
    for date in dates:
        print date
        _HRU_outputs_on_date = HRU_outputs.loc[HRU_outputs['timestamp'] == date]
        _values = []
        for i in range(_shape.GetFeatureCount()):
            _feature = _shape.GetFeature(i)
            _nhru = _feature['id']
            #print _nhru
            _row = _HRU_outputs_on_date.loc[_HRU_outputs_on_date['nhru'] == _nhru]
            try:
                _values.append(float(_row[plotting_variable].values))
            except:
                #_values.append(np.nan)
                continue
        _values = np.array(_values)
        # Floating colorbar
        #colors = cm.jet(plt.Normalize( min(_values), max(_values)) (_values) )
        colors = cm.jet(plt.Normalize( _min, _max) (_values) )
        ax.cla()
        _polys = []
        for i in range(_shape.GetFeatureCount()):
            _feature = _shape.GetFeature(i)
            _geometry = _feature.geometry()
            _boundary = _geometry.GetBoundary()
            _boundary_points = np.array(_boundary.GetPoints())
            _x = _boundary_points[:,0]/1000.
            _y = _boundary_points[:,1]/1000.
            #plt.plot(_x, _y, 'k-')
            _polys.append(ax.fill(_x, _y, facecolor=colors[i], edgecolor='k')[0])
        ax.set_title(plotting_variable+': '+date)
        ax.set_xlabel('E [km]', fontsize=16)
        ax.set_ylabel('N [km]', fontsize=16)
        #plt.tight_layout()
        writer.grab_frame()
        plt.pause(0.1)
    
