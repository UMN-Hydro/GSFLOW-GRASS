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

# Segments
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

# Create a map to work with
g.copy(vector=('streams','segments'), overwrite=True)

v.db_addcolumn(map='segments', columns=segment_columns)



# Produce the data table entries
##################################

colNames = np.array(gscript.vector_db_select('segments', layer=1)['columns'])
colValues = np.array(gscript.vector_db_select('segments', layer=1)['values'].values())
number_of_segments = colValues.shape[0]
cats = colValues[:,colNames == 'cat'].astype(int).squeeze()

nseg = np.arange(1, len(cats)+1)
nseg_cats = []
for i in range(len(cats)):
  nseg_cats.append( (nseg[i], cats[i]) )

# Default OUTSEG to 0: will stay that way if there is no valid update value
# (works if there is a stream that goes somewhere else that is numbered but 
# that isn't in the current basin)
v.db_update(map='segments', column='OUTSEG', value=0)
v.db_update(map='streams', column='OUTSEG', value=0)

# Somehow only works after I v.clean, not right after v.overlay
segment = VectorTopo('segments')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segments set id=? where cat=?", nseg_cats)
cur.executemany("update segments set OUTSEG=? where tostream=?", nseg_cats)
segment.table.conn.commit()
segment.close()

# Do the same for streams, at least for the ID
v.db_addcolumn(map='streams', columns='id integer, OUTSEG integer')
v.db_update(map='streams', column='OUTSEG', value=0)
segment = VectorTopo('streams')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update streams set id=? where cat=?", nseg_cats)
cur.executemany("update streams set OUTSEG=? where tostream=?", nseg_cats)
segment.table.conn.commit()
segment.close()

# Power law: COEFFS IN METERS -- ASK CRYSTAL IF OK
v.db_update(map='segments', column='ICALC', value=3)
v.db_update(map='segments', column='CDPTH', value=0.25)
v.db_update(map='segments', column='FDPTH', value=0.4)
v.db_update(map='segments', column='AWDTH', value=8.)
v.db_update(map='segments', column='BWDTH', value=0.5)

# Segment IDs
v.db_update(map='segments', column='ISEG', value='id')
v.db_update(map='segments', column='NSEG', value='id')

# values that are 0
v.db_update(map='segments', column='IUPSEG', value=0)
v.db_update(map='segments', column='FLOW', value=0)
v.db_update(map='segments', column='RUNOFF', value=0)
v.db_update(map='segments', column='ETSW', value=0)
v.db_update(map='segments', column='PPTSW', value=0)
"""
# Streamflow and lake routing
# tosegment

tosegment_cats = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
tosegment = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
# From outlet segment
for i in range(len(xy2)):
  # to outlet segment
  outlets = np.prod(xy2 == xy1[i], axis=1)
  # Update outlet segments with ID of inlets
  tosegment[outlets.nonzero()] = nhru[i]
  tosegment_cats[outlets.nonzero()] = cats[i]

# Now, just update tosegment (segments) and hru_segment (hru's)
# In this case, they are the same.
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general
# Tuple for upload to SQL
# 0 is the default value if it doesn't go into any other segment (i.e flows
# off-map)
tosegmentt = []
tosegment_cats_t = []
for i in range(len(nsegment)):
  tosegmentt.append( (tosegment[i], nsegment[i]) )
  tosegment_cats_t.append( (tosegment_cats[i], cats[i]) )
# Once again, special case
hru_segmentt = tosegmentt

# Loop check!
# Weak loop checker - will only detect direct ping-pong.
loops = []
tosegmenta = np.array(tosegmentt)
for i in range(len(tosegmenta)):
  for j in range(len(tosegmenta)):
    if (tosegmenta[i] == tosegmenta[j][::-1]).all():
      loops.append(tosegmenta[i])

segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set tosegment=? where id=?", tosegmentt)
segment.table.conn.commit()
segment.close()


# In study basin?
grass.run_command('v.db.addcolumn', map='segment', columns='in_study_basin int')
grass.run_command('v.db.addcolumn', map='HRU', columns='in_study_basin int')
grass.run_command('v.what.vect', map='segment', column='in_study_basin', query_map='studyBasin', query_column='value')
grass.run_command('v.what.vect', map='HRU', column='in_study_basin', query_map='segment', query_column='in_study_basin')


# Save global segment
grass.run_command('g.rename', vect='segment,segment_all')


# Output
colNames = np.array(grass.vector_db_select('segment')['columns'])
colValues = np.array(grass.vector_db_select('segment')['values'].values())


# Same for segments
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general

# Somehow only works after I v.clean, not right after v.overlay
segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set id=? where cat=?", nsegmentt)
segment.table.conn.commit()
segment.close()


tosegment_cats = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
tosegment = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
# From outlet segment
for i in range(len(xy2)):
  # to outlet segment
  outlets = np.prod(xy2 == xy1[i], axis=1)
  # Update outlet segments with ID of inlets
  tosegment[outlets.nonzero()] = nhru[i]
  tosegment_cats[outlets.nonzero()] = cats[i]

# Now, just update tosegment (segments) and hru_segment (hru's)
# In this case, they are the same.
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general
# Tuple for upload to SQL
# 0 is the default value if it doesn't go into any other segment (i.e flows
# off-map)
tosegmentt = []
tosegment_cats_t = []
for i in range(len(nsegment)):
  tosegmentt.append( (tosegment[i], nsegment[i]) )
  tosegment_cats_t.append( (tosegment_cats[i], cats[i]) )
# Once again, special case
hru_segmentt = tosegmentt

# Loop check!
# Weak loop checker - will only detect direct ping-pong.
loops = []
tosegmenta = np.array(tosegmentt)
for i in range(len(tosegmenta)):
  for j in range(len(tosegmenta)):
    if (tosegmenta[i] == tosegmenta[j][::-1]).all():
      loops.append(tosegmenta[i])


segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set tosegment=? where id=?", tosegmentt)
segment.table.conn.commit()
segment.close()
"""


#v.db.select segments sep=comma > segment.csv

