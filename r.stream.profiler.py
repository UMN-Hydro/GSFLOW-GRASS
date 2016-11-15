# PyGRASS river profile analysis
# By Andrew Wickert
# Started 14 October, 2016

from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.gis import region

from grass.pygrass import vector
import numpy as np
from matplotlib import pyplot as plt

from grass.script import vector_db_select

"""
Required inputs:
- Flow accumulation map (i.e. drainage area -- can also scale to discharge)
- Slope from r.slope.aspect
- Vector map of river channels (look at r.gsflow)

Outputs:
- Channel networks as sets of points along a line that can be queried
- Slope, area, normalized steepness
- Locations of likely channel heads
- If run in "full landscape" mode, gives hillslope vs. channel via drainage direction (will I ever get this far?)
- For each channel:
    * k_s
    * chi?
    * remember, these are meant to work ASSUMING STREAM POWER
        --> and therefore probably assuming bedrock, though all quite
            fuzzy in stream-power land anyway.

First run r.stream.extract
Then run r.slope.aspect
"""

##########
# INPUTS #
##########

# Full set of commands:
elevation = 'topoclip'
#slope = 'srtm.slope'
#streams = 'srtm.stream
# or THRESH if not STREAMS -- in km**2
thresh = 0.2 # 10 km2 = minimum catchment size
thresh = 10000 # m2

# Made to work on a projected coordinate system PROJECTED!
reg = region.Region()

##############
# BEFOREHAND #
##############
# These also will generate inputs: code to run beforehand

# Slope first, easy
print "Computing slope with r.slope.aspect. Slope in unitless decimal values."
r.slope_aspect(elevation=elevation, slope='tmp', format='percent', overwrite=True)
r.mapcalc('slope = tmp/100.', overwrite=True)

# Assuming you work in meters. ASSUMING!
# Have to include output name or write to a temporary file, FOR ALL!
r.mapcalc('cellArea_meters2 = '+str(reg.nsres)+' * '+str(reg.ewres), overwrite=True)
r.mapcalc("cellArea_km2 = cellArea_meters2 / 10^6", overwrite=True)

print "Running r.watershed"
r.watershed(elevation=elevation, flow='cellArea_meters2', accumulation='drainageArea_m2', drainage='drainageDirection', stream='streams', threshold=thresh, flags='s', overwrite=True)
# for irregular outlines on topo
r.mapcalc("drainageArea_m2 = drainageArea_m2 + 0*"+elevation, overwrite=True)
#r.watershed(elevation=elevation, flow='cellArea_km2', accumulation='drainageArea_km2', drainage='drainageDirection', stream='streams', threshold=thresh, flags='s', overwrite=True)
# Remove areas of negative (offmap) accumulation
#r.mapcalc('drainageArea_km2 = drainageArea_km2 * (drainageArea_km2 > 0)', overwrite=True)
#r.null(map='drainageArea_km2', setnull=0)
#r.mapcalc(elevation+" = "+elevation+"*drainageArea_km2*0", overwrite=True)

# Get watershed
print "Building drainage network"
r.stream_extract(elevation=elevation, accumulation='drainageArea_m2', threshold=thresh, d8cut=0, mexp=0, stream_raster='streams', stream_vector='streams', direction='draindir', overwrite=True)

"""
# Get slope and area
v.db_addcolumn(map='streams_points', columns=('slope double precision, area_km2 double precision'))
v.what_rast(map='streams_points', type='point', raster='slope', column='slope')
v.what_rast(map='streams_points', type='point', raster='drainageArea_km2', column='area_km2')
"""

"""
~~~~~~~~~~~~~~~~~~~~~
awickert@dakib:~$ Now that I have a good ordering scheme, find adjacency between units based on starting and ending points and then convert each of these to points and get the slopes and areas between them. Concatenate these all into a single line and see what kind of averaging is needed -- of position, slope, and area. Poisition should include a base downstream distance as well, so we can keep averaged lines on the same course.
06 NOV
~~~~~~~~~~~~~~~~~~~~~
"""

####################
# COMPUTE NETWORKS #
####################

# The following is not following upstream/downstream conventions.
# Time to do this the hard way.
# 1. Get vectorTopo
# 2. Get coordinates
# 3. Get areas at coordinates
# 4. Sort points to go from small A to large A -- but not yet uploading
# 5. Upload small area as x1, y1; large area as x2, y2

from grass.pygrass.vector import Vector, VectorTopo
from grass.pygrass.raster import RasterRow
from grass.pygrass import utils
streamsTopo = VectorTopo('streams', overwrite=True)
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
# 3. Get areas at coordinates
drainageArea_km2 =  RasterRow('drainageArea_km2')
drainageArea_km2.open('r')
"""
##########
streamsTopo.open(mode='rw')
for row in streamsTopo:
  row = a
streamsTopo.build()
streamsTopo.close()

vector.libvect.Vect_cat_set(geo_obj.c_cats, 1, 1)
result = libvect.Vect_rewrite_line(self.c_mapinfo,
                                   cat, geo_obj.gtype,
                                   geo_obj.c_points,
                                   geo_obj.c_cats)


streamsTopo.open(mode='rw')
for i in range(1,len(streamsTopo)+1):
  cat = i
  geometry = streamsTopo.cat(cat,'lines')[0]
  _A_point1 = drainageArea_km2.get_value(points_in_streams[i][0])
  _A_point2 = drainageArea_km2.get_value(points_in_streams[i][-1])
  if _A_point1 > _A_point2:
    streamsTopo.rewrite(geometry.reverse(), cat=i, attrs=None)
    #row = vector.geometry.Line(row[::-1])
streamsTopo.build()
streamsTopo.write()c
streamsTopo.close()
##########
"""
streamsTopo.table.columns.add('drainageArea_km2_1','double precision')
streamsTopo.table.columns.add('drainageArea_km2_2','double precision')
streamsTopo.table.columns.add('x1','double precision')
streamsTopo.table.columns.add('y1','double precision')
streamsTopo.table.columns.add('x2','double precision')
streamsTopo.table.columns.add('y2','double precision')

cur = streamsTopo.table.conn.cursor()
for i in range(len(points_in_streams)):
  _A_point1 = drainageArea_km2.get_value(points_in_streams[i][0])
  _A_point2 = drainageArea_km2.get_value(points_in_streams[i][-1])
  # 4. Sort points to go from small A to large A
  #streamsTopo[i+1] = points_in_streams[i].reverse()
  # 5. Upload small area as x1, y1; large area as x2, y2
  # CATS had better be in unbroken ascending order!
  # Areas
  cur.execute("update streams set drainageArea_km2_1="+str(_A_point1)+" where cat="+str(cat_of_line_segment[i]))
  cur.execute("update streams set drainageArea_km2_2="+str(_A_point2)+" where cat="+str(cat_of_line_segment[i]))
  if _A_point1 > _A_point2:
    # Points
    cur.execute("update streams set x1="+str(points_in_streams[i][-1].x)+" where cat="+str(cat_of_line_segment[i]))
    cur.execute("update streams set y1="+str(points_in_streams[i][-1].y)+" where cat="+str(cat_of_line_segment[i]))
    cur.execute("update streams set x2="+str(points_in_streams[i][0].x)+" where cat="+str(cat_of_line_segment[i]))
    cur.execute("update streams set y2="+str(points_in_streams[i][0].y)+" where cat="+str(cat_of_line_segment[i]))
  else:
    print "WARNING!!!!"
    # Points
    cur.execute("update streams set x1="+str(points_in_streams[i][0].x)+" where cat="+str(cat_of_line_segment[i]))
    cur.execute("update streams set y1="+str(points_in_streams[i][0].y)+" where cat="+str(cat_of_line_segment[i]))
    cur.execute("update streams set x2="+str(points_in_streams[i][-1].x)+" where cat="+str(cat_of_line_segment[i]))
    cur.execute("update streams set y2="+str(points_in_streams[i][-1].y)+" where cat="+str(cat_of_line_segment[i]))
  print i
#streamsTopo.write()
streamsTopo.table.conn.commit()
streamsTopo.build()

# THEN:
# 5. Every line should have 

"""
# CLOSE TO BEING DONE WITH OLD RIVER NUMBERS,
# NOW THAT I JUST THRESHOLD DRAINAGE AREA
colNames = np.array(vector_db_select('streams')['columns'])
colValues = np.array(vector_db_select('streams')['values'].values())
cats = colValues[:,colNames == 'cat'].astype(int).squeeze()
river_numbers = colValues[:,colNames == 'river_number'].astype(int).squeeze()
drainageArea_km2_1 = colValues[:,colNames == 'drainageArea_km2_1'].astype(float).squeeze() # area at upstream end
xy1 = colValues[:,(colNames == 'x1') + (colNames == 'y1')].astype(float) # upstream
xy2 = colValues[:,(colNames == 'x2') + (colNames == 'y2')].astype(float) # downstream
xy  = np.vstack((xy1, xy2))
"""

"""
# new river numbers: ascending order from headwaters downstream.
# Sort first by area, and then by river number.
# Raise a warning if both are the same at some point, and I will have to rethink
# this approach.
# ... actually, some work to break ties
# #river_numbers_list = list(
# see http://stackoverflow.com/questions/19643099/sorting-a-list-of-tuples-with-multiple-conditions
# so not much work, but some!
# But just sorting by area for now
river_numbers = list(cats.copy()) # Just to have an ascending integer list [1, n]
cats_list = list(cats.copy()) # Just to have an ascending integer list [1, n]
area1_list = list(drainageArea_km2_1)
keydict = dict(zip(cats, area1_list))
cats_list.sort(key=keydict.get)

# THIS IS HOW YOU UPDATE A TABLE!!!!
v.db_addcolumn(map='streams', columns='stream_number_ascending int')
#streams.table.columns.add('stream_number_ascending','int')
streams = Vector('streams')
streams.open('rw')
cur = streams.table.conn.cursor()
for i in range(len(river_numbers)):
  print i
  cur.execute("update streams set stream_number_ascending="+str(river_numbers[i])+" where cat="+str(cats_list[i]))
streams.table.conn.commit()
streams.build()
"""

# CAT IS RIVER NUMBER FOR R.STREAM....
colNames = np.array(vector_db_select('streams')['columns'])
colValues = np.array(vector_db_select('streams')['values'].values())
river_numbers = colValues[:,colNames == 'cat'].astype(int).squeeze()
# Could also just use our river_numbers vector from before

# For points:
#slope = colValues[:,colNames == 'slope'].astype(int).squeeze()
#area = colValues[:,colNames == 'area_m2'].astype(int).squeeze()


# NEW #######################
#####################
# OR -- USE V.TO.DB # ALL THAT IS NEEDED!
#####################
v.to_db(map='streams', option='start', columns='x1,y1')
v.to_db(map='streams', option='end', columns='x2,y2')
# And v.what.rast to find a way to bring in ID's
# NOT NEEDED -- CATS ARE ALL IN ORDER WITH V.STREAM....
#v.db_addcolumn(map='streams',columns="Id int")
#v.what_rast(map='streams', raster='streams', column='Id')

# CLOSE TO BEING DONE WITH OLD RIVER NUMBERS,
# NOW THAT I JUST THRESHOLD DRAINAGE AREA
colNames = np.array(vector_db_select('streams')['columns'])
colValues = np.array(vector_db_select('streams')['values'].values())
cats = colValues[:,colNames == 'cat'].astype(int).squeeze()
# CAT IS RIVER NUMBER FOR R.STREAM....
river_numbers = colValues[:,colNames == 'cat'].astype(int).squeeze()
"""
drainageArea_km2_1 = colValues[:,colNames == 'drainageArea_km2_1'].astype(float).squeeze() """
# area at upstream end
xy1 = colValues[:,(colNames == 'x1') + (colNames == 'y1')].astype(float) # upstream
xy2 = colValues[:,(colNames == 'x2') + (colNames == 'y2')].astype(float) # downstream
xy  = np.vstack((xy1, xy2))

############################

# Need ascending river numbers

# So now can use this information to find headwaters and mouths
tosegment_river_numbers = np.nan * np.zeros(colValues.shape[0]).astype(int) # default to 0 if they do not flow to another segment
tosegment_cats = tosegment_river_numbers.copy()
# From outlet segment
for i in range(len(xy2)):
  from_river_number = river_numbers[i]
  # to outlet segment
  outlets = np.prod(xy2[i] == xy1, axis=1)
  # Update outlet segments with ID of inlets
  # in GRASS, segment numbers get larger downstream, so pick the biggest one
  if np.sum(outlets) == 0:
    tosegment_river_numbers[i] = 0
  else:
    # IMPORTANT: this use of np.max is why we had these rivers ordered
    # in terms of their basin area
    tosegment_river_numbers[i] = np.max(cats[outlets.nonzero()])
    tosegment_cats[i] = np.max(cats[outlets.nonzero()])

tosegment_river_numbers = tosegment_river_numbers.astype(int)

# This gives us a set of downstream-facing adjacencies.
# First let us update the database table with it
v.db_addcolumn(map='streams', columns='tostream int')
streams = Vector('streams')
streams.open('rw')
cur = streams.table.conn.cursor()
for i in range(len(tosegment_river_numbers)):
  print i
  cur.execute("update streams set tostream="+str(tosegment_river_numbers[i])+" where cat="+str(cats[i])) # CATS[I] WAS THE FIX
streams.table.conn.commit()
streams.build()



colNames = np.array(vector_db_select('streams')['columns'])
colValues = np.array(vector_db_select('streams')['values'].values())
stream_number_ascending = colValues[:,colNames == 'stream_number_ascending'].astype(int).squeeze()
tostream = colValues[:,colNames == 'tostream'].astype(int).squeeze()

# THIS MUST BE A SEPARATE MODULE
#################################

#####################################################
# GET EVERYTHING DOWNSTREAM FROM THIS RIVER SEGMENT #
#####################################################

# We can loop over this list to get the shape of the full river network.
full_river_cats = []
#segment = 406
#segment = 712 # Machai?
# 92, 597X, 114, 64
segment = 64 # CHANGEABLE  -- MAKE INPUT TO A DIFFERENT FUNCTION!!!
#segment=1
full_river_cats.append(segment)
while full_river_cats[-1] != 0:
  full_river_cats.append(int(tostream[cats == full_river_cats[-1]]))
full_river_cats = full_river_cats[:-1] # remove 0 at end

full_river_cats_str = list(np.array(full_river_cats).astype(str))
full_river_cats_csv = ','.join(full_river_cats_str)

v.extract(input='streams', output='specific_stream', cats=full_river_cats_csv, overwrite=True)



#v.dissolve(input='specific_stream_segmented', output='specific_stream', column='river_number', overwrite=True)
#v.db_addtable('specific_stream')
#v.category('specific_stream')
v.to_points(input='specific_stream', output='specific_stream_points', use='vertex', dmax=90, flags='i', overwrite=True)
v.db_addcolumn(map='specific_stream_points', layer=2, columns=('slope double precision, area_km2 double precision'))
v.what_rast(map='specific_stream_points', type='point', raster='slope', column='slope', layer=2)
v.what_rast(map='specific_stream_points', type='point', raster='drainageArea_km2', column='area_km2', layer=2)

# Elevation
v.db_addcolumn(map='specific_stream_points', layer=2, columns=('elevation double precision'))
v.what_rast(map='specific_stream_points', type='point', raster=elevation, column='elevation', layer=2)

########################################
# THIS CREATES THE RIVER PROFILE PLOTS #
########################################

SA = vector_db_select(map='specific_stream_points', columns='lcat,along,slope,area_km2,elevation', layer=2)
SA = np.asarray(SA.items()[0][1].values())
SA = SA[(SA[:,2] != '') * (SA[:,3] != '')] # no empty values
SA = SA.astype(float)
SA = SA[SA[:,3] > 0] # Only positive areas allowed
SA = np.vstack((SA[0], SA[1:][np.sum(np.diff(SA, axis=0), axis=1) != 0])) # Remove duplicates at nodes -- duplicates from same segment, in fact!
x = SA[:,1].copy()
for i in range(1, len(x)):
  if SA[:,1][i] == 0:
    x[i:] += SA[i-1,1]
# Now delete duplicates between columns
SA = np.vstack((SA[0], SA[1:][np.diff(x) > 1E-9]))
x = np.hstack((x[0], x[1:][np.diff(x) > 1E-9]))
S = SA[:,2].copy()
A = SA[:,3].copy()
z = SA[:,4].copy()

# Averaged
_x = []
_S = []
_A = []
_z = []
# window
L = 50 # m window
dL = L # moving steps
l = 0
r = L
rmax = x[-1]
while r < rmax:
  _criterion = [(x>=l)*(x<r)]
  _x.append(np.mean(x[_criterion]))
  _S.append(np.mean(S[_criterion]))
  _A.append(np.mean(A[_criterion]))
  _z.append(np.mean(z[_criterion]))
  l += dL
  r += dL

plt.ion()
plt.figure('Long profile')
plt.title('Long profile', fontsize=16)
plt.xlabel('$x$ [m]', fontsize=16)
plt.ylabel('$z$ [m]', fontsize=16)
plt.plot(_x,_z)
plt.savefig('Long profile.png')
plt.figure('Slope--distance')
plt.plot(_x,_S, 'ko')
plt.title('Slope--distance', fontsize=16)
plt.xlabel('$x$ [m]', fontsize=16)
plt.ylabel('$S$ [--]', fontsize=16)
plt.savefig('Sx.png')
plt.figure('Slope--area')
plt.title('Slope--area', fontsize=16)
plt.xlabel('$A$ [m^2]', fontsize=16)
plt.ylabel('$S$ [--]', fontsize=16)
plt.loglog(np.array(_A)*1E6, _S, 'k.')
plt.ylim((1E-2, 10))
plt.savefig('SA.png')
plt.figure('Slope--area tall', figsize=(8,12))
plt.title('Slope--area', fontsize=16)
plt.xlabel('$A$ [m^2]', fontsize=16)
plt.ylabel('$S$ [--]', fontsize=16)
plt.loglog(np.array(_A)*1E6, _S, 'k.')
plt.ylim((1E-2, 100))
plt.savefig('SAtall.png')
plt.figure('Area--distance')
plt.title('Area--distance', fontsize=16)
plt.xlabel('$x$ [m]', fontsize=16)
plt.ylabel('$A$ [m^2]', fontsize=16)
plt.plot(_x, np.array(_A)*1E6, 'k.')
plt.savefig('Ax.png')
plt.show()

# Get segment connectivity

###########################
# Make a slope--area plot #
###########################

#SA = vector_db_select(map='streams_points', columns='river_number,slope,area_m2')
#SA = np.asarray(SA.items()[0][1].values()).astype(float)

#SA2 = SA[SA[:,0] == 11814]
plt.ion()
plt.loglog(SA[:,1]/1E6, SA[:,0], 'k.')
plt.show()

# NOW GET SLOPE AT ALL POINTS ALONG DRAINAGES, ... 


# THIS CAUSES PROBLEMS HERE W/ NEARBY AREAS W/ SAME ACCUM, NOT TAKING DRAINGE DIRECTION INTO ACCOUNT. SO USED STREAM(ABOVE)
# BUT NOW HOW TO GET TOPOLOGY? LET'S SEE....
# AH, PROBLEM (BELOW) MORE ABOUT BASINS. ABOVE PROBLEM STILL A PROBLEM THERE
# BUT MY NETWORK SHOULD BE FINE.
"""
# Manually create streams from accumulation -- threshold should be provided by user.
# The one funny step is the cleaning w/ snap, because r.thin allows cells that are
# diagonal to each other to be next to each other -- creating boxes along the channel
# that are not 
r.mapcalc('streams_unthinned = flowAccum > '+str(thresh), overwrite=True)
r.null(map='streams_unthinned', setnull=0)
r.thin(input='streams_unthinned', output='streams', overwrite=True)
r.to_vect(input='streams', output='streams_raw', type='line', overwrite=True)
v.clean(input='streams_raw', output='streams', tool='snap', threshold=1.42*(reg.nsres + reg.ewres)/2., flags='c', overwrite=True) # threshold is one cell
v.to_rast(input='streams', output='streams_unthinned', use='val', val=1, overwrite=True)
r.thin(input='streams_unthinned', output='streams', overwrite=True)
r.to_vect(input='streams', output='streams', type='line', overwrite=True)
v.to_rast(input='streams', output='streams', use='cat', overwrite=True)
# Create drainage basins -- skipping for this tool (all borrowed from r.gsflow)
#grass.run_command('r.stream.basins', direction='drainageDirection', stream_rast='streams', basins='basins', overwrite=True)
# If there is any more need to work with nodes, I should check the code I wrote for Kelly Monteleone's paper -- this has river identification and extraction, including intersection points.
"""

