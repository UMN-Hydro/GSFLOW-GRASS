#!/usr/bin/env python
############################################################################
#
# MODULE:       r.stream.network
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Attach IDs of upstream and downstream nodes as well as the
#               category value of the next downstream stream segment
#               (0 if the stream exits the map)
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
#%option G_OPT_V_INPUT
#%  key: streams
#%  label: Vector input of stream network created by r.stream.extract
#%  required: yes
#%  guidependency: layer,column
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
    streams = options['streams']

    streamsTopo = VectorTopo(streams)
    streamsTopo.build()

    # 1. Get vectorTopo
    streamsTopo.open(mode='rw')
    points_in_streams = []
    cat_of_line_segment = []

    # 2. Get coordinates
    for row in streamsTopo:
        cat_of_line_segment.append(row.cat)
        if type(row) == vector.geometry.Line:
            points_in_streams.append(row)

    # 3. Coordinates of points: 1 = start, 2 = end
    try:
        streamsTopo.table.columns.add('x1','double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add('y1','double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add('x2','double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add('y2','double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add('tostream','int')
    except:
        pass

    # Is this faster than v.to.db?
    # v.to_db(map='streams', option='start', columns='x1,y1')
    # v.to_db(map='streams', option='end', columns='x2,y2')
    cur = streamsTopo.table.conn.cursor()
    for i in range(len(points_in_streams)):
        cur.execute("update streams set x1="+str(points_in_streams[i][0].x)+" where cat="+str(cat_of_line_segment[i]))
        cur.execute("update streams set y1="+str(points_in_streams[i][0].y)+" where cat="+str(cat_of_line_segment[i]))
        cur.execute("update streams set x2="+str(points_in_streams[i][-1].x)+" where cat="+str(cat_of_line_segment[i]))
        cur.execute("update streams set y2="+str(points_in_streams[i][-1].y)+" where cat="+str(cat_of_line_segment[i]))
    streamsTopo.table.conn.commit()
    streamsTopo.build()

    colNames = np.array(vector_db_select('streams')['columns'])
    colValues = np.array(vector_db_select('streams')['values'].values())
    cats = colValues[:,colNames == 'cat'].astype(int).squeeze() # river number
    xy1 = colValues[:,(colNames == 'x1') + (colNames == 'y1')].astype(float) # upstream
    xy2 = colValues[:,(colNames == 'x2') + (colNames == 'y2')].astype(float) # downstream

    # Build river network
    tocat = []
    for i in range(len(cats)):
        tosegment_mask = np.prod(xy1 == xy2[i], axis=1)
        if np.sum(tosegment_mask) == 0:
            tocat.append(0)
        else:
            tocat.append(tosegment_mask.nonzero()[0][0])
    tocat = np.asarray(tocat).astype(int)

    # This gives us a set of downstream-facing adjacencies.
    # We will update the database with it.
    streamsTopo.build()
    streamsTopo.open('rw')
    cur = streamsTopo.table.conn.cursor()
    for i in range(len(tocat)):
        cur.execute("update streams set tostream="+str(tocat[i])+" where cat="+str(cats[i]))
    streamsTopo.table.conn.commit()
    streamsTopo.build()

    print ""
    print "Done."
    print ""

if __name__ == "__main__":
    main()

