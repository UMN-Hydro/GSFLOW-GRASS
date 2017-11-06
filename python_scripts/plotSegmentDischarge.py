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


# *** Save movie to following file
moviefile_name = 'testmovie_strmseg.mp4'


#%% from Settings file, change to plot something else
segout_fil = Settings.PRMSoutput_dir + slashstr + Settings.PROJ_CODE + '.ani.nsegment'
projdir_GIS = Settings.GISinput_dir
segshp_fil = projdir_GIS + slashstr + "shapefiles/segments/segments.shp"

print '\n******************************************'
print 'Plotting results from: ', segout_fil
print ' (segments shapefile: ' + segshp_fil + ')'
print '******************************************\n'

#%% In general: don't change below here

plotting_variable = 'streamflow_sfr'

#for filename in [HRUout_fil, segment_filename]:
for filename in [segout_fil]:
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

segout_fil2 = segout_fil + '.corrected'

_shapefile = ogr.Open(segshp_fil)
_shape = _shapefile.GetLayer(0)

segment_outputs = pd.read_csv(segout_fil2, comment='#', delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False, skiprows=[8])

dates = sorted(list(set(list(segment_outputs.timestamp))))
cmap = plt.get_cmap('RdYlBu')

plotting_variable = 'streamflow_sfr'
_min = np.min(segment_outputs[plotting_variable])
_max = np.max(segment_outputs[plotting_variable])

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
        _segment_outputs_on_date = segment_outputs.loc[segment_outputs['timestamp'] == date]
        _values = []
        for i in range(_shape.GetFeatureCount()):
            _feature = _shape.GetFeature(i)
            _n = _feature['id']
            #print _nhru
            _row = _segment_outputs_on_date.loc[_segment_outputs_on_date['nsegment'] == _n]
            try:
                _values.append(float(_row[plotting_variable].values))
            except:
                _values.append(np.nan)
                print _n
                continue
        _values = np.array(_values)
        # Floating colorbar
        colors = cm.jet(plt.Normalize( _min, _max) (_values) )
        ax.cla()
        _lines = []
        for i in range(_shape.GetFeatureCount()):
            _feature = _shape.GetFeature(i)
            #feature = shape.GetFeature(0) # how to get it otherwise
            _geometry = _feature.geometry()
            _line_points = np.array(_geometry.GetLinearGeometry().GetPoints())
            _x = _line_points[:,0]/1000.
            _y = _line_points[:,1]/1000.
            _lines.append( ax.plot(_x, _y, '-', color=colors[i], linewidth=_values[i]**.5+.25) )
        ax.set_title(plotting_variable+': '+date)
        ax.set_xlabel('E [km]', fontsize=16)
        ax.set_ylabel('N [km]', fontsize=16)
        #plt.tight_layout()
        writer.grab_frame()
        plt.pause(0.01)
        #plt.waitforbuttonpress()
    
