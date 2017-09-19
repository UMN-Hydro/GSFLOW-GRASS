#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.network
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Build stream segments for the USGS models MODFLOW or GSFLOW.
#
# COPYRIGHT:    (c) 2016-2017 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  uses inputs from v.stream.network
 
# More information
# Started December 2016

#%module
#% description: Prepares stream segments for PRMS and GSFLOW
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#%end

#%option G_OPT_V_INPUT
#%  key: input
#%  label: Vector stream network from r.stream.extract
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_OUTPUT
#%  key: output
#%  label: Segments: stream segments for GSFLOW / PRMS
#%  required: yes
#%  guidependency: layer,column
#%end

#%option
#%  key: ICALC
#%  type: integer
#%  description: Stream depth option: 0-const; 1,2-Manning; 3-aQ^b
#%  answer: 3
#%  required: no
#%end

#%option
#%  key: CDPTH
#%  type: double
#%  description: Flow depth coefficient; used if ICALC=3
#%  answer: 0.25
#%  required: no
#%end

#%option
#%  key: FDPTH
#%  type: double
#%  description: Flow depth exponent; used if ICALC=3
#%  answer: 0.4
#%  required: no
#%end

#%option
#%  key: AWDTH
#%  type: double
#%  description: Flow width coefficient; used if ICALC=3
#%  answer: 8
#%  required: no
#%end

#%option
#%  key: BWDTH
#%  type: double
#%  description: Flow width exponent; used if ICALC=3
#%  answer: 0.5
#%  required: no
#%end

#%option
#%  key: IUPSEG
#%  type: string
#%  description: Category of upstream diversion segment (from_cat,to_cat,...)
#%  answer: 0,0
#%  required: no
#%end

#%option
#%  key: FLOW
#%  type: string
#%  description: Streamflow entering the upstream-most segments (cat,Q,cat,Q,...)
#%  answer: 0,0
#%  required: no
#%end

#%option
#%  key: RUNOFF
#%  type: string
#%  description: Diffuse runoff entering each segment (cat,Q,cat,Q,...)
#%  answer: 0,0
#%  required: no
#%end

#%option
#%  key: ETSW
#%  type: string
#%  description: Direct removal of in-channel water by ET (cat,Q,cat,Q)
#%  answer: 0,0
#%  required: no
#%end

#%option
#%  key: PPTSW
#%  type: string
#%  description: Direct precipitation on the stream
#%  answer: 0,0
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
from grass.pygrass.modules.shortcuts import miscellaneous as m
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
    Builds river segments for input to the USGS hydrologic models
    PRMS and GSFLOW.
    """

    ##################
    # OPTION PARSING #
    ##################

    options, flags = gscript.parser()
    
    # I/O
    streams = options['input']
    segments = options['output']
    
    # Hydraulic geometry
    ICALC = options['ICALC']
    
    # ICALC=0: Constant depth - NOT IMPLEMENTED
    
    # ICALC=0: Manning - NOT IMPLEMENTED

    # ICALC=0: Manning - NOT IMPLEMENTED

    # ICALC=3: Power-law relationships (following Leopold and others)
    CDPTH = options['CDPTH']
    FDPTH = options['FDPTH']
    AWDTH = options['AWDTH']
    BWDTH = options['BWDTH']
    
    ##################################################
    # CHECKING DEPENDENCIES WITH OPTIONAL PARAMETERS #
    ##################################################
    
    if ICALC == 3:
        if CDPTH and FDPTH and AWDTH and BWDTH:
            pass
        else:
            grass.fatal('Missing CDPTH, FDPTH, AWDTH, and/or BWDTH. \
                         These are required when ICALC = 3.')

    ###########
    # RUNNING #
    ###########

    # New Columns for Segments
    segment_columns = []
    # Self ID
    segment_columns.append('id integer') # segment number
    segment_columns.append('ISEG integer') # segment number
    segment_columns.append('NSEG integer') # segment number
    # Streamflow and lake routing
    #segment_columns.append('tosegment integer') # Index of downstream segment to which a segment
                                                #   flows (thus differentiating it from hru_segment,
                                                #   which is for HRU's, though segment and HRU ID's
                                                #   are the same when HRU's are sub-basins
    # for GSFLOW
    segment_columns.append('ICALC integer') # 3 for power function
    segment_columns.append('OUTSEG integer') # downstream segment -- tostream, renumbered
    segment_columns.append('CDPTH double precision') # depth coeff
    segment_columns.append('FDPTH double precision') # depth exp
    segment_columns.append('AWDTH double precision') # width coeff
    segment_columns.append('BWDTH double precision') # width exp
    # The below will be all 0
    segment_columns.append('IUPSEG integer') # upstream segment ID number, for diversions
    segment_columns.append('FLOW integer')
    segment_columns.append('RUNOFF integer')
    segment_columns.append('ETSW integer')
    segment_columns.append('PPTSW integer')

    segment_columns = ",".join(segment_columns)

    # CONSIDER THE EFFECT OF OVERWRITING COLUMNS -- WARN FOR THIS
    # IF MAP EXISTS ALREADY?

    # Create a map to work with
    g.copy(vector=(streams, segments), overwrite=gscript.overwrite())
    # and add its columns
    v.db_addcolumn(map=segments, columns=segment_columns)

    # Produce the data table entries
    ##################################
    colNames = np.array(gscript.vector_db_select(segments, layer=1)['columns'])
    colValues = np.array(gscript.vector_db_select(segments, layer=1)['values'].values())
    number_of_segments = colValues.shape[0]
    cats = colValues[:,colNames == 'cat'].astype(int).squeeze()

    nseg = np.arange(1, len(cats)+1)
    nseg_cats = []
    for i in range(len(cats)):
        nseg_cats.append( (nseg[i], cats[i]) )

    # Default OUTSEG to 0: will stay that way if there is no valid update value
    # (works if there is a stream that goes somewhere else that is numbered but 
    # that isn't in the current basin)
    v.db_update(map=segments, column='OUTSEG', value=0)

    segment = VectorTopo(segments)
    segment.open('rw')
    cur = segment.table.conn.cursor()

    # id = cat (as does ISEG and NSEG)
    cur.executemany("update segments set id=? where cat=?", nseg_cats)
    cur.executemany("update segments set ISEG=? where cat=?", nseg_cats)
    cur.executemany("update segments set NSEG=? where cat=?", nseg_cats)

    # outseg = tostream
    cur.executemany("update segments set OUTSEG=? where tostream=?", nseg_cats)

    # Discharge and hydraulic geometry
    cur.execute("update segments set ICALC="+str(ICALC))
    cur.execute("update segments set CDPTH="+str(CDPTH))
    cur.execute("update segments set FDPTH="+str(FDPTH))
    cur.execute("update segments set AWDTH="+str(AWDTH))
    cur.execute("update segments set BWDTH="+str(BWDTH))

    gscript.message('')
    gscript.message('NOTICE: not currently used:')
    gscript.message('IUPSEG, FLOW, RUNOFF, ETSW, and PPTSW.')
    gscript.message('All set to 0.')
    gscript.message('')

    # values that are 0
    cur.execute("update segments set IUPSEG="+str(0))
    cur.execute("update segments set FLOW="+str(0))
    cur.execute("update segments set RUNOFF="+str(0))
    cur.execute("update segments set ETSW="+str(0))
    cur.execute("update segments set PPTSW="+str(0))

    segment.table.conn.commit()
    segment.close()

if __name__ == "__main__":
    main()

