from osgeo import ogr
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from readSettings import Settings
import platform
import sys
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

# *** Enter list of HRU numbers to identify byhighlighting in plot
hru_to_ID = np.array([1, 23])

#%% from Settings file, change to plot something else
projdir_GIS = Settings.GISinput_dir
HRUshp_fil = projdir_GIS + slashstr + "shapefiles/HRUs/HRUs.shp"
segshp_fil = projdir_GIS + slashstr + "shapefiles/segments/segments.shp"

print '\n******************************************'
print 'Plotting from HRU shapefile: ' + HRUshp_fil
print '******************************************\n'


#%%

plt.ion()

fig = plt.figure()
ax = plt.subplot(111)
ax.set_xlabel('E [km]', fontsize=16)
ax.set_ylabel('N [km]', fontsize=16)

_shapefile = ogr.Open(HRUshp_fil)
_shape = _shapefile.GetLayer(0)
#first feature of the shapefile
for i in range(_shape.GetFeatureCount()):
    _feature = _shape.GetFeature(i)
    #feature = shape.GetFeature(0) # how to get it otherwise
    _geometry = _feature.geometry()
    _boundary = _geometry.GetBoundary()
    _boundary_points = np.array(_boundary.GetPoints())
    _x = _boundary_points[:,0]/1000.
    _y = _boundary_points[:,1]/1000.
    if  (hru_to_ID == _feature['id']).any():
        plt.fill(_x, _y, facecolor='.8')
    else:
        plt.plot(_x, _y, 'k-')
ti = 'HRU '
for ii in hru_to_ID:
    ti = ti + str(ii)
    if ii != hru_to_ID[-1]:
        ti = ti + ', '
plt.title(ti)

_shapefile = ogr.Open(segshp_fil)
_shape = _shapefile.GetLayer(0)
for i in range(_shape.GetFeatureCount()):
    _feature = _shape.GetFeature(i)
    #feature = shape.GetFeature(0) # how to get it otherwise
    _geometry = _feature.geometry()
    _line_points = np.array(_geometry.GetLinearGeometry().GetPoints())
    _x = _line_points[:,0]/1000.
    _y = _line_points[:,1]/1000.
    plt.plot(_x, _y, 'b-', linewidth=2)
