#!/usr/bin/env python

# See: https://grasswiki.osgeo.org/wiki/Working_with_GRASS_without_starting_it_explicitly

##################
# IMPORT MODULES #
##################
# PYTHON
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import subprocess
from ConfigParser import ConfigParser

#########################
# READ CONFIG VARIABLES #
#########################
config = ConfigParser()
config.read('/home/awickert/dataanalysis/GRASS-fluvial-profiler/settings.ini')

# Global input variables
# project_name becomes GRASS location name
project_name = config.get('settings', 'proj_name')
DEM_input = config.get('GRASS', 'DEM_file_path_to_import')
A_threshold = config.get('GRASS', 'threshold_drainage_area_meters2')
MODFLOW_grid_resolution = config.get('GRASS', 'MODFLOW_grid_resolution_meters')
outlet_point_x = config.get('GRASS', 'outlet_point_x')
outlet_point_y = config.get('GRASS', 'outlet_point_y')
icalc = config.get('GRASS', 'icalc')
gisdb = config.get('GRASS', 'gisdb')
version = config.get('GRASS', 'version')

########################
# RUN GRASS IMPLICITLY #
########################

# path to the GRASS GIS launch script
# MS Windows
grass7bin_win = r'C:\OSGeo4W\bin\grass'+version+'svn.bat'
# uncomment when using standalone WinGRASS installer
# grass7bin_win = r'C:\Program Files (x86)\GRASS GIS 7.2.0\grass72.bat'
# Linux
grass7bin_lin = 'grass'+version
# Mac OS X
# this is TODO
grass7bin_mac = '/Applications/GRASS/GRASS-'+version[0]+'.'+version[1]+'.app/'

# specify (existing) location and mapset
location = project_name
mapset   = "PERMANENT"
location_path = os.path.join(gisdb, location)

########### SOFTWARE
if sys.platform.startswith('linux'):
    # we assume that the GRASS GIS start script is available and in the PATH
    # query GRASS 7 itself for its GISBASE
    grass7bin = grass7bin_lin
elif sys.platform.startswith('win'):
    grass7bin = grass7bin_win
else:
    raise OSError('Platform not configured.')

# Using existing location
if os.path.isdir(location_path):
    # query GRASS 7 itself for its GISBASE
    startcmd = [grass7bin, '--config', 'path']
    p = subprocess.Popen(startcmd, shell=False,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
        sys.exit(-1)
    gisbase = out.strip('\n\r')
# Making new location
else:
    startcmd = grass7bin + ' --config path'
    p = subprocess.Popen(startcmd, shell=True, 
					     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
  	    print >>sys.stderr, 'ERROR: %s' % err
  	    print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
  	    sys.exit(-1)
    if sys.platform.startswith('linux'):
  	    gisbase = out.strip('\n')
    elif sys.platform.startswith('win'):
        if out.find("OSGEO4W home is") != -1:
    		    gisbase = out.strip().split('\n')[1]
        else:
    		    gisbase = out.strip('\n')
        os.environ['GRASS_SH'] = os.path.join(gisbase, 'msys', 'bin', 'sh.exe')

# Set GISBASE environment variable
os.environ['GISBASE'] = gisbase
#os.environ['PATH'] += os.pathsep + os.path.join(gisbase, 'extrabin') # not needed with trunk
os.environ['PATH'] += os.pathsep + os.path.join(os.getenv('HOME'), '.grass7', 'addons', 'scripts')

# define GRASS-Python environment
gpydir = os.path.join(gisbase, "etc", "python")
sys.path.append(gpydir)

########
# define GRASS DATABASE
if sys.platform.startswith('win'):
    gisdb = os.path.join(os.getenv('APPDATA'), 'grassdata')
else:
    gisdb = os.path.join(os.getenv('HOME'), 'grassdata')

if os.path.isdir(location_path):
    pass
else:
    #  from georeferenced data file
    startcmd = grass7bin + ' -c ' + DEM_input + ' -e ' + location_path
    p = subprocess.Popen(startcmd, shell=True, 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >>sys.stderr, 'ERROR: %s' % err
        print >>sys.stderr, 'ERROR: Cannot generate location (%s)' % startcmd
        sys.exit(-1)
    else:
        print 'Created location %s' % location_path

# Set GISDBASE environment variable
os.environ['GISDBASE'] = gisdb

# Linux: Set path to GRASS libs (TODO: NEEDED?)
# THE ISSUE LIES HERE!!!!!!
path = os.getenv('LD_LIBRARY_PATH')
dir  = os.path.join(gisbase, 'lib')
if path:
    path = dir + os.pathsep + path
else:
    path = dir
os.environ['LD_LIBRARY_PATH'] = path

# language
os.environ['LANG'] = 'en_US'
os.environ['LOCALE'] = 'C'

# import GRASS Python bindings (see also pygrass)
import grass.script.setup as gsetup
 
###########
# launch session
gsetup.init(gisbase,
            gisdb, location, mapset)

# GRASS
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.gis import region
from grass.pygrass import vector # Change to "v"?
from grass.script import vector_db_select
from grass.pygrass.vector import Vector, VectorTopo
from grass.pygrass.raster import RasterRow
from grass.pygrass import utils
from grass import script as gscript
from grass.pygrass.vector.geometry import Point

# Internal variables: set names
DEM_original_import = 'DEM_original_import'   # Raw DEM
DEM                 = 'DEM'                   # DEM after offmap flow removed
DEM_MODFLOW         = 'DEM_MODFLOW'           # DEM for MODFLOW
cellArea_meters2    = 'cellArea_meters2'      # Grid cell size
accumulation        = 'accumulation'          # Flow accumulation
accumulation_onmap  = 'accumulation_onmap'    # Flow accum: no off-map flow
draindir            = 'draindir'              # Drainage direction
streams_all         = 'streams_all'           # Streams on full map
streams_inbasin     = 'streams_inbasin'       # Streams in the study basin
streams_MODFLOW     = 'streams_MODFLOW'       # Streams on MODFLOW grid
basins_all          = 'basins_all'            # All watershed subbasins
basins_inbasin      = 'basins_inbasin'        # Subbasins in the study basin
segments            = 'segments'              # Stream segments
reaches             = 'reaches'               # Stream reaches (for MODFLOW)
MODFLOW_grid        = 'grid'                  # MODFLOW grid vector
slope               = 'slope'                 # Topographic slope
aspect              = 'aspect'                # Topographic aspect
HRUs                = 'HRUs'                  # Hydrologic response units
gravity_reservoirs  = 'gravity_reservoirs'    # Connect HRUs to MODFLOW grid
basin_mask          = 'basin_mask'            # Mask out non-study-basin cells
pour_point          = 'pour_point'            # Outlet pour point
bc_cell             = 'bc_cell'               # Grid cell for MODFLOW b.c.

DEM_orig=DEM_orig
DEM=DEM
DEM_coarse=DEM_coarse # MODFLOW resolution
accumulation=accumulation_tmp
streams=streams_tmp
streams_MODFLOW=streams_MODFLOW
streams_onebasin=${streams}_onebasin
basins=basins_tmp
basins_onebasin=${basins}_onebasin
segments=segments_tmp
reaches=reaches_tmp
threshold=1000000 #m2 - drainage area
grid_res=500 #150 #1000 #m2 - for MODFLOW
grid=grid_tmp
slope=slope_tmp
aspect=aspect_tmp
HRUs=HRUs_tmp
gravity_reservoirs=gravity_reservoirs_tmp
basin_mask=basin_mask_tmp
basin_mask_out=basin_mask
pour_point=pp_tmp
bc_cell=bc_cell
icalc=1 # how to compute hydraulic geometry
x_outlet=482452.076604
y_outlet=8672978.0598

# Import DEM if required
# And perform the standard starting tasks.
# These take time, so skip if not needed
if DEM_input != '':
    # Import DEM and set region
    r.in_gdal(input=DEM_input, output=DEM_original_import)
    g.region(raster=DEM_original_import)
    # Build flow accumulation with only fully on-map flow
    # Cell areas
    r.cell_area(output=cellArea_meters2, units='m2', overwrite=True)
    # Hydrologic correction
    r.hydrodem(input=DEM_original_import, output=DEM, flags='a', overwrite=True)
    # No offmap flow
    r.watershed(elevation=DEM, flow=cellArea_meters2, accumulation=accumulation, flags='m', overwrite=True)
    r.mapcalc(accumulation_onmap+' = '+accumulation+' * ('+accumulation+' > 0)', overwrite=True)
    r.null(map=accumulation_onmap, setnull=0)
    r.mapcalc(DEM+' = if(isnull('+accumulation_onmap+'),null(),'+DEM+')', overwrite=True)
    # Ensure that null cells are shared
    r.mapcalc(accumulation_onmap+' = if(isnull('+DEM+'),null(),'+accumulation_onmap+')', overwrite=True)
    # Repeat is sometimes needed
    r.mapcalc(DEM+' = if(isnull('+accumulation_onmap+'),null(),'+DEM+')', overwrite=True)
    r.mapcalc(accumulation_onmap+' = if(isnull('+DEM+'),null(),'+accumulation_onmap+')', overwrite=True)


# Set region
g.region(raster=DEM_original_import)

# Build streams and sub-basins
r.stream_extract(elevation=DEM, accumulation=accumulation_onmap, stream_raster=streams_all, stream_vector=streams_all, threshold=A_threshold, direction=draindir, d8cut=0, overwrite=True)
r.stream_basins(direction=draindir, stream_rast=streams_all, basins=basins_all, overwrite=True)
r.to_vect(input=basins_all, output=basins_all, type='area', flags='v', overwrite=True)

# Build stream network
v.stream_network(map=streams_all)

# Restrict to a single basin
v.stream_inbasin(input_streams=streams_all, input_basins=basins_all, output_streams=streams_inbasin, output_basin=basins_inbasin, x_outlet=outlet_point_x, y_outlet=outlet_point_y, output_pour_point=pour_point, overwrite=True)

# GSFLOW segments: sections of stream that define subbasins
v.gsflow_segments(input=streams_inbasin, output=segments, icalc=icalc, overwrite=True)

# MODFLOW grid & basin mask (1s where basin exists and 0 where it doesn't)
v.gsflow_grid(basin=basins_inbasin, pour_point=pour_point, raster_input=DEM, dx=MODFLOW_grid_resolution, dy=MODFLOW_grid_resolution, output=MODFLOW_grid, mask_output=basin_mask, bc_cell=bc_cell, overwrite=True)

# Hydrologically-correct DEM for MODFLOW
r.gsflow_hydrodem(dem=DEM, grid=MODFLOW_grid, streams=streams, resolution=MODFLOW_grid_resolution, streams_modflow=streams_MODFLOW, dem_modflow=DEM_MODFLOW, overwrite=True)

# GSFLOW reaches: intersection of segments and grid
v.gsflow_reaches(segment_input=segments, grid_input=MODFLOW_grid, elevation=DEM, output=reaches, overwrite=True)

# GSFLOW HRU parameters
r.slope_aspect(elevation=DEM, slope=slope, aspect=aspect, format='percent', zscale=0.01, overwrite=True)
v.gsflow_hruparams(input=basins_inbasin, elevation=DEM, output=HRUs, slope=slope, aspect=aspect, overwrite=True)

# GSFLOW gravity reservoirs
v.gsflow_gravres(hru_input=HRUs, grid_input=MODFLOW_grid, output=gravity_reservoirs, overwrite=True)

# Export DEM with MODFLOW resolution
# Also export basin mask -- 1s where basin exists and 0 where it doesn't
# And make sure it is in an appropriate folder
if os.path.split(os.getcwd())[-1] != project_name:
    try:
        os.mkdir(project_name)
    except:
        pass
os.chdir(project_name)
g.region(raster=DEM_MODFLOW)
r.out_ascii(input=DEM_MODFLOW, output='DEM.asc', null_value=0, overwrite=True)
r.out_ascii(input=basin_mask, output=basin_mask+'.asc', null_value=0, overwrite=True)
g.region(raster=DEM)
#g.region vect=$basins_onebasin

# Export tables and discharge point
v.gsflow_export(reaches_input=reaches,
                segments_input=segments,
                gravres_input=gravity_reservoirs,
                hru_input=HRUs,
                pour_point_input=pour_point,
                bc_cell_input=bc_cell,
                reaches_output=reaches,
                segments_output=segments,
                gravres_output=gravity_reservoirs,
                hru_output=HRUs,
                pour_point_output=pour_point,
                bc_cell_output=bc_cell,
                overwrite=True)
                
