import numpy as np
from matplotlib import pyplot as plt
import sys
import warnings
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
reach_columns = []
# Self ID
reach_columns.append('KRCH integer')
reach_columns.append('IRCH integer')
reach_columns.append('JRCH integer')
reach_columns.append('ISEG integer')
reach_columns.append('IREACH integer')
reach_columns.append('RCHLEN integer')
reach_columns.append('STRTOP double precision')
reach_columns.append('SLOPE double precision')
reach_columns.append('STRTHICK double precision')

reach_columns = ",".join(reach_columns)

# Create a map to work with
v.extract(input='streams', output='tmp2', type='line', overwrite=True)
v.overlay(ainput='tmp2', atype='line', binput='grid', output='reaches', operator='and', overwrite=True)

v.db_addcolumn(map='reaches', columns=reach_columns)

# Rename a,b columns
v.db_renamecolumn(map='reaches', column=('a_x1', 'x1'))
v.db_renamecolumn(map='reaches', column=('a_x2', 'x2'))
v.db_renamecolumn(map='reaches', column=('a_y1', 'y1'))
v.db_renamecolumn(map='reaches', column=('a_y2', 'y2'))
v.db_renamecolumn(map='reaches', column=('a_stream_type', 'stream_type'))
v.db_renamecolumn(map='reaches', column=('a_type_code', 'type_code'))
v.db_renamecolumn(map='reaches', column=('a_cat', 'rnum_cat'))
v.db_renamecolumn(map='reaches', column=('a_tostream', 'tostream'))
v.db_renamecolumn(map='reaches', column=('a_id', 'segment_id'))
v.db_renamecolumn(map='reaches', column=('a_OUTSEG', 'OUTSEG'))
v.db_renamecolumn(map='reaches', column=('b_row', 'row'))
v.db_renamecolumn(map='reaches', column=('b_col', 'col'))
v.db_renamecolumn(map='reaches', column=('b_id', 'cell_id'))

# Drop some unnecessary columns
v.db_dropcolumn(map='reaches', columns='b_area_m2')

# Update some columns that can be done now
v.db_update(map='reaches', column='KRCH', value=1)
v.db_update(map='reaches', column='IRCH', value='row')
v.db_update(map='reaches', column='JRCH', value='col')
v.db_update(map='reaches', column='ISEG', value='segment_id')
v.to_db(map='reaches', columns='RCHLEN', option='length')
v.db_update(map='reaches', column='STRTHICK', value=0.1) # 10 cm, prescribed

# Still to go after these:
# STRTOP (added with slope)
# IREACH (whole next section dedicated to this)
# SLOPE (need z_start and z_end)

# Now, the light stuff is over: time to build the reach order
v.db_addcolumn(map='reaches', columns='xr1 double precision, yr1 double precision, xr2 double precision, yr2 double precision')
v.to_db(map='reaches', option='start', columns='xr1,yr1')
v.to_db(map='reaches', option='end', columns='xr2,yr2')

# Now just sort by category, find which stream has the same xr1 and yr1 as
# x1 and y1 (or a_x1, a_y1) and then find where its endpoint matches another 
# starting point and move down the line.
# v.db.select reaches col=cat,a_id,xr1,xr2 where="a_x1 = xr1"

# First, get the starting coordinates of each stream segment
# and a set of river ID's (ordered from 1...N)
colNames = np.array(gscript.vector_db_select('segments', layer=1)['columns'])
colValues = np.array(gscript.vector_db_select('segments', layer=1)['values'].values())
number_of_segments = colValues.shape[0]
segment_x1s = colValues[:,colNames == 'x1'].astype(float).squeeze()
segment_y1s = colValues[:,colNames == 'y1'].astype(float).squeeze()
segment_ids = colValues[:,colNames == 'id'].astype(float).squeeze()

# Then move back to the reaches map to produce the ordering
colNames = np.array(gscript.vector_db_select('reaches', layer=1)['columns'])
colValues = np.array(gscript.vector_db_select('reaches', layer=1)['values'].values())
reach_cats = colValues[:,colNames == 'cat'].astype(int).squeeze()
reach_x1s = colValues[:,colNames == 'xr1'].astype(float).squeeze()
reach_y1s = colValues[:,colNames == 'yr1'].astype(float).squeeze()
reach_x2s = colValues[:,colNames == 'xr2'].astype(float).squeeze()
reach_y2s = colValues[:,colNames == 'yr2'].astype(float).squeeze()
segment_ids__reach = colValues[:,colNames == 'segment_id'].astype(float).squeeze()

for segment_id in segment_ids:
  reach_order_cats = []
  downstream_directed = []
  ssel = segment_ids == segment_id
  rsel = segment_ids__reach == segment_id # selector
  # Find first segment: x1y1 first here, but not necessarily later
  downstream_directed.append(1)
  _x_match = reach_x1s[rsel] == segment_x1s[ssel]
  _y_match = reach_y1s[rsel] == segment_y1s[ssel]
  _i_match = _x_match * _y_match
  x1y1 = True # false if x2y2
  # Find cat
  _cat = int(reach_cats[rsel][_x_match * _y_match])
  reach_order_cats.append(_cat)
  # Get end of reach = start of next one
  reach_x_end = float(reach_x2s[reach_cats == _cat])
  reach_y_end = float(reach_y2s[reach_cats == _cat])
  # Make flexible for the answer being in 1 or 2: not directional
  # BUT IT LOOKS LIKE SEGMENTS ARE DIRECTIONAL! THIS CODE WILL STAY FLEXIBLE
  # BUT WILL ISSUE A WARNING
  while _i_match.any():
    _x_match = reach_x1s[rsel] == reach_x_end
    _y_match = reach_y1s[rsel] == reach_y_end
    _i_match = _x_match * _y_match
    if _i_match.any():
      _cat = int(reach_cats[rsel][_x_match * _y_match])
      reach_x_end = float(reach_x2s[reach_cats == _cat])
      reach_y_end = float(reach_y2s[reach_cats == _cat])
      #downstream_directed.append(1)
      reach_order_cats.append(_cat)
    """
    # Does not seem necessary :)
    else:
      _x_match = reach_x2s[rsel] == int(reach_x2s[reach_cats == _cat])
      _y_match = reach_y2s[rsel] == int(reach_y2s[reach_cats == _cat])
      _cat_valid = np.ones(_x_match.shape) - (reach_cats[rsel] == _cat)
      _i_match = _x_match * _y_match * _cat_valid
      if _i_match.any():
        _cat = int(reach_cats[rsel][_x_match * _y_match])
        reach_x_end = float(reach_x1s[reach_cats == _cat])
        reach_y_end = float(reach_y1s[reach_cats == _cat])
        downstream_directed.append(0)
        warnings.warn('Stream reach oriented upstream!')
    if _i_match.any():
      reach_order_cats.append(_cat)
    """
  print len(reach_order_cats), len(reach_cats[rsel])
    
  # Reach order to database table
  reach_order_cats_cats = []
  for i in range(len(reach_order_cats)):
    reach_order_cats_cats.append( (reach_order_cats[i], reach_cats[rsel][i]) )
  vname = 'reaches'
  vect = VectorTopo(vname)
  vect.open('rw')
  cur = vect.table.conn.cursor()
  cur.executemany("update "+vname+" set IREACH=? where cat=?", reach_order_cats_cats)
  vect.table.conn.commit()
  vect.close()
  

# TOP AND BOTTOM ARE OUT OF ORDER: SOME SEGS ARE BACKWARDS. UGH!!!!
# NEED TO GET THEM IN ORDER TO GET THE Z VALUES AT START AND END

# Compute slope and starting elevations from the elevations at the start and 
# end of the reaches and the length of each reach
v.db_addcolumn(map='reaches', columns='zr1 double precision, zr2 double precision')
zr1 = []
zr2 = []
for i in range(len(reach_cats)):
  _x = reach_x1s[i]
  _y = reach_y1s[i]
  _z = float(gscript.parse_command('r.what', map='srtm_local_filled', coordinates=str(_x)+','+str(_y)).keys()[0].split('|')[-1])
  zr1.append(_z)
  _x = reach_x2s[i]
  _y = reach_y2s[i]
  _z = float(gscript.parse_command('r.what', map='srtm_local_filled', coordinates=str(_x)+','+str(_y)).keys()[0].split('|')[-1])
  zr2.append(_z)

zr1_cats = []
zr2_cats = []
for i in range(len(reach_cats)):
  zr1_cats.append( (zr1[i], reach_cats[i]) )
  zr2_cats.append( (zr2[i], reach_cats[i]) )

# There is no reason that I have to upload both of these to the attribute 
# table, but why not?
vname = 'reaches'
vect = VectorTopo(vname)
vect.open('rw')
cur = vect.table.conn.cursor()
cur.executemany("update "+vname+" set zr1=? where cat=?", zr1_cats)
cur.executemany("update "+vname+" set zr2=? where cat=?", zr2_cats)
cur.executemany("update segments set OUTSEG=? where tostream=?", nseg_cats)
vect.table.conn.commit()
vect.close()

v.db_update(map='reaches', column='STRTOP', value='zr1')
v.db_update(map='reaches', column='SLOPE', value='(zr1 - zr2)/RCHLEN')

