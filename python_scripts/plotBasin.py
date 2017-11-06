from osgeo import ogr
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec

projdir = '/home/awickert/dataanalysis/GRASS-fluvial-profiler/Shullcas_2lay/'
hru_to_ID = np.array([1, 23])

plt.ion()

fig = plt.figure()
ax = plt.subplot(111)
ax.set_xlabel('E [km]', fontsize=16)
ax.set_ylabel('N [km]', fontsize=16)

_shapefile = ogr.Open(projdir + "shapefiles/HRUs/HRUs.shp")
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

_shapefile = ogr.Open(projdir + "shapefiles/segments/segments.shp")
_shape = _shapefile.GetLayer(0)
for i in range(_shape.GetFeatureCount()):
    _feature = _shape.GetFeature(i)
    #feature = shape.GetFeature(0) # how to get it otherwise
    _geometry = _feature.geometry()
    _line_points = np.array(_geometry.GetLinearGeometry().GetPoints())
    _x = _line_points[:,0]/1000.
    _y = _line_points[:,1]/1000.
    plt.plot(_x, _y, 'b-', linewidth=2)
