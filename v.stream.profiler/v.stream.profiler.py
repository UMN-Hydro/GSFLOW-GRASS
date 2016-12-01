#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.profiler
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Build long profiles and slope--area diagrams of a river network
#
# COPYRIGHT:    (c) 2016 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  uses inputs from r.stream.extract
 
# More information
# Started 14 October 2016
#%module
#% description: Build a linked stream network: each link knows its downstream link
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#%end
#%option
#%  key: cat
#%  label: Starting line segment category
#%  required: yes
#%  guidependency: layer,column
#%end
#%option G_OPT_V_INPUT
#%  key: streams
#%  label: Vector input of stream network created by r.stream.extract
#%  required: yes
#%end
#%option G_OPT_V_OUTPUT
#%  key: outstream
#%  label: Vector output stream
#%  required: yes
#%end
#%option
#%  key: direction
#%  type: string
#%  label: Which directon to march: up or down
#%  options: upstream,downstream
#%  answer: downstream
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: elevation
#%  label: Topography (DEM)
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: accumulation
#%  label: Flow accumulation
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: slope
#%  label: Map of slope created by r.slope.area
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: window
#%  label: Distance over which to average in river profiles and slope-area plots
#%  required: no
#%end
#%option
#%  key: units
#%  type: string
#%  label: Flow accumulation units
#%  options: m2, km2, cumecs, cfs
#%  required: no
#%end
#%option
#%  key: plots
#%  type: string
#%  label: Plots to generate
#%  options: LongProfile,SlopeArea,SlopeDistance,AreaDistance
#%  required: no
#%end

##################
# IMPORT MODULES #
##################
# PYTHON
import numpy as np
from matplotlib import pyplot as plt
import sys
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

###############
# MAIN MODULE #
###############

def main():
    """
    Links each river segment to the next downstream segment in a tributary 
    network by referencing its category (cat) number in a new column. "0"
    means that the river exits the map.
    """

    options, flags = gscript.parser()
    
    print options
    print flags

    # Attributes of streams
    colNames = np.array(vector_db_select(options['streams'])['columns'])
    colValues = np.array(vector_db_select(options['streams'])['values'].values())
    tostream = colValues[:,colNames == 'tostream'].astype(int).squeeze()
    cats = colValues[:,colNames == 'cat'].astype(int).squeeze() # = "fromstream"

    # We can loop over this list to get the shape of the full river network.
    selected_cats = []
    segment = int(options['cat'])
    selected_cats.append(segment)
    if options['direction'] == 'downstream':
        while selected_cats[-1] != 0:
            selected_cats.append(int(tostream[cats == selected_cats[-1]]))
        selected_cats = selected_cats[:-1] # remove 0 at end
    elif options['direction'] == 'upstream':
        print "Not yet active!"
        """
        # Add new lists for each successive upstream river
        river_is_upstream =
        while
        full_river_cats
        """
    
    selected_cats_str = list(np.array(selected_cats).astype(str))
    selected_cats_csv = ','.join(selected_cats_str)

    v.extract(input=options['streams'], output=options['outstream'], cats=selected_cats_csv, overwrite=True)

if __name__ == "__main__":
    main()

