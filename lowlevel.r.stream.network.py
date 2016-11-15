#!/usr/bin/env python
############################################################################
#
# MODULE:       v.flexure
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Calculate flexure of the lithosphere under a specified
#               set of loads and with a given elastic thickness (scalar)
#
# COPYRIGHT:    (c) 2014, 2015 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  gFlex: http://csdms.colorado.edu/wiki/gFlex
#         (should be downloaded automatically along with the module)
#         github repository: https://github.com/awickert/gFlex
 
# More information
# Started 20 Jan 2015 to add GRASS GIS support for distributed point loads
# and their effects on lithospheric flexure
#%module
#% description: Lithospheric flexure: gridded deflections from scattered point loads
#% keyword: vector
#% keyword: stream
#% keyword: geomorphometry
#% keyword: hydrology
#% keyword: geomorphology
#%end
#%option G_OPT_R_INPUT
#%  key: elevation
#%  description: Vector map of loads (thickness * area * density * g) [N]
#%  guidependency: layer,column
#%end
#%option G_OPT_R_INPUT
#%  key: cellsize
#%  description: Vector map of loads (thickness * area * density * g) [N]
#%  guidependency: layer,column
#%end
#%option G_OPT_V_FIELD
#%  key: layer
#%  description: Layer containing load values
#%  guidependency: column
#%end
#%option G_OPT_DB_COLUMNS
#%  key: column
#%  description: Column containing load values [N]
#%  required : yes
#%end
#%option
#%  key: te
#%  type: double
#%  description: Elastic thicnkess: scalar; unis chosen in "te_units"
#%  required : yes
#%end
#%option
#%  key: te_units
#%  type: string
#%  description: Units for elastic thickness
#%  options: m, km
#%  required : yes
#%end
#%option G_OPT_V_OUTPUT
#%  key: output
#%  description: Output vector points map of vertical deflections [m]
#%  required : yes
#%end
#%option G_OPT_R_OUTPUT
#%  key: raster_output
#%  description: Output raster map of vertical deflections [m]
#%  required : no
#%end
#%option
#%  key: g
#%  type: double
#%  description: gravitational acceleration at surface [m/s^2]
#%  answer: 9.8
#%  required : no
#%end
#%option
#%  key: ym
#%  type: double
#%  description: Young's Modulus [Pa]
#%  answer: 65E9
#%  required : no
#%end
#%option
#%  key: nu
#%  type: double
#%  description: Poisson's ratio
#%  answer: 0.25
#%  required : no
#%end
#%option
#%  key: rho_fill
#%  type: double
#%  description: Density of material that fills flexural depressions [kg/m^3]
#%  answer: 0
#%  required : no
#%end
#%option
#%  key: rho_m
#%  type: double
#%  description: Mantle density [kg/m^3]
#%  answer: 3300
#%  required : no
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

#############
# WORKSPACE #
#############

# Unsorted arrays: what is at start and what is at end?
streamsTopo = VectorTopo('streams', overwrite=True)
streamsTopo.build()
streamsTopo.open(mode='rw')

# Order everything
cats_full = [] # points and lines
for row in streamsTopo:
  cats_full.append(row.cat)
streamsTopo.rewind()
cats_full_order = np.argsort(cats)
streamsTopo_a = np.array(streamsTopo)
streamsTopo_a = streamsTopo_a[cats_full_order] # Place in order of ascending category

# NEW STRATEGY -- NO NEED TO ORGANIZE. SIMPLY FIND START AND END POINTS!
# Those points with id 0 will occur only at the start of the whole system.
# These are the headwaters; march downstream from here.
points = []
point_type_codes = [] # 0 = headwaters; 1 = downstream
point_cats = []
lines = []
line_type_codes = [] # 0 = headwaters; 1 = downstream
line_cats = []
streamsTopo = VectorTopo('streams', overwrite=True)
streamsTopo.build()
streamsTopo.open(mode='rw')
for i in range(1, len(streamsTopo)+1):
  if type(streamsTopo[i]) == vector.geometry.Line:
    lines.append(streamsTopo[i].tolist())
    line_type_codes.append(streamsTopo[i].attrs['type_code'])
    line_cats.append(streamsTopo[i].cat)
  elif type(streamsTopo[i]) == vector.geometry.Point:
    points.append(streamsTopo[i].coords())
    point_type_codes.append(streamsTopo[i].attrs['type_code'])
    point_cats.append(streamsTopo[i].cat)
  else:
    sys.exit("Non-line, non-point, in stream network. Why?")
points = np.asarray(points)
point_type_codes = np.asarray(point_type_codes)
point_cats = np.asarray(point_cats)
lines = np.asarray(lines)
line_type_codes = np.asarray(line_type_codes)
line_cats = np.asarray(line_cats)

# Order everything -- at least will make thinking easier
point_cats_order = np.argsort(point_cats)
line_cats_order = np.argsort(line_cats)
points = points[point_cats_order]
point_type_codes = point_type_codes[point_cats_order]
point_cats = point_cats[point_cats_order]
lines = lines[line_cats_order]
line_type_codes = line_type_codes[line_cats_order]
line_cats = line_cats[line_cats_order]

# Gives repeat points, but not line segments. Still, be cautious with both
_cat, index = np.unique(point_cats, return_index=True)
points = points[index]
point_type_codes = point_type_codes[index]
point_cats = point_cats[index]
# ADD MATCHING LINE SECTION HERE IF IT BECOMES NECESSARY!

# Use the type codes as a mask: 
# Define those points and lines at the headwaters (start) and downstream
headwaters_points = points[point_type_codes == 0]
headwaters_point_cats = point_cats[point_type_codes == 0]
lower_points = points[point_type_codes]
lower_point_cats =  point_cats[point_type_codes]
headwaters_lines = lines[line_type_codes == 0]
headwaters_line_cats = line_cats[line_type_codes == 0]
lower_lines = lines[line_type_codes]
lower_line_cats =  line_cats[line_type_codes]

# THIS SHOWS THAT THEY ARE ALL IN THE PROPER ORDER!
for i in range(len(headwaters_lines)):
  _mask = np.prod(headwaters_points == headwaters_lines[i][0], axis=1)
  if np.sum(_mask):
    pass
  else:
    print i
    
    
# SO ASSUMING EVERYTHING IS IN ORDER, THEN WE CAN SIMPLIFY THINGS:
# LET'S SEE IF THIS WORKS
# Lower in-network points: up- and down-stream
lower_upstream_points = []
lower_downstream_points = []
for line in lower_lines:
  lower_upstream_points.append(line[0])
  lower_downstream_points.append(line[-1])

upstream_points = []
downstream_points = []
for line in lines:
  upstream_points.append(line[0])
  downstream_points.append(line[-1])
upstream_points = np.asarray(upstream_points)
downstream_points = np.asarray(downstream_points)

tocat = []
for i in range(len(lines)):
  tosegment_mask = np.prod(upstream_points == downstream_points[i], axis=1)
  if np.sum(tosegment_mask) == 0:
    tocat.append(0)
  else:
    tocat.append(tosegment_mask.nonzero()[0][0])
  
  
# Number of networks = headwaters
toline_lol = []
for i in range(len(headwaters_lines)):
  topoint = headwaters_lines[i][-1]
  nextline = downstream_lines[
  toline_list = []
  toline_list.append()



  print np.prod(headwaters_points == headwaters_lines[i][0], axis=1).nonzero()[0][0]
  print (headwaters_points == headwaters_lines[i][0]).nonzero()[0][0]












# Lists for lines in proper order
#lines_in_order = lines.copy()


# Points and start/end of lines -- and orient lines properly!
P1 = []
P2 = []
P1a = []
P2a = []
cats = []
for row in streamsTopo:
  if type(row) == vector.geometry.Line:
    #print row.cat
    line = row
    P1.append(line[0])
    P2.append(line[-1])
    P1a.append(line[0].coords())
    P2a.append(line[-1].coords())
    cats.append(line.cat)









# Points and start/end of lines
P1 = []
P2 = []
P1a = []
P2a = []
cats = []
for row in streamsTopo:
  if type(row) == vector.geometry.Line:
    #print row.cat
    line = row
    P1.append(line[0])
    P2.append(line[-1])
    P1a.append(line[0].coords())
    P2a.append(line[-1].coords())
    cats.append(line.cat)
    
# Indices for categories in order
cats_order = np.argsort(cats)

# Create arrays to make comparisons easier, and put in order.
P1a = np.array(P1a)[cats_order]
P2a = np.array(P2a)[cats_order]

# Check if this river is either the start or end of the network
isterminal = []
isterminal1 = []
isterminal2 = []
for i in range(P1a.shape[0]):
  _is_not_terminal1 = (P1a[i,:] == P1a[:i]).any() or \
                      (P1a[i,:] == P1a[i+1:]).any() or \
                      (P1a[i,:] == P2a).any()
  _is_not_terminal2 = (P2a[i,:] == P2a[:i]).any() or \
                      (P2a[i,:] == P2a[i+1:]).any() or \
                      (P2a[i,:] == P1a).any()
  isterminal1.append(_is_not_terminal1 == False)
  isterminal2.append(_is_not_terminal2 == False)
  if isterminal1[-1] and isterminal2[-1]:
    isterminal.append(-1) # start and end
  elif isterminal1[-1] or isterminal2[-1]:
    isterminal.append(1) # start or end
  else:
    isterminal.append(0) # intermediate
isterminal1 = np.array(isterminal1)
isterminal2 = np.array(isterminal2)
isterminal = np.array(isterminal)
  
  if :
     # Potential to be terminal -- 
    if (P1a[i,:] == P1a[:i]).any() or (P1a[i,:] == P1a[i+1:]).any() or \
       (P1a[i,:] == P2a).any():
    isterminal.append(False)
  else:
    isterminal.append(True)

#####################
# UTILITY FUNCTIONS #
#####################

def add_upstream_downstream_points(streams,)
    """
    Add points at the upstream and downstream ends of river segments into the
    attribute table: postiion, drainage area
    x1, y1: upstream
    x2, y2: downstream
    """
    streamsTopo = VectorTopo('streams', overwrite=True)
    streamsTopo.build()
    # 1. Get vectorTopo
    streamsTopo.open(mode='rw')
    points_in_streams = []
    cat_of_line = []
    # 2. Get coordinates
    for row in streamsTopo:
      cat_of_line.append(row.cat)
      if type(row) == vector.geometry.Line:
        points_in_streams.append(row)
    # 3. Get areas at coordinates
    drainageArea_km2 =  RasterRow('drainageArea_km2')
    drainageArea_km2.open('r')
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
      if _A_point1 > _A_point2:
        # Areas
        cur.execute("update streams set drainageArea_km2_1="+\
                    str(_A_point2)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set drainageArea_km2_2="+str(_A_point1)+\
                    " where cat="+str(cat_of_line[i]))
        # Points
        cur.execute("update streams set x1="+str(points_in_streams[i][-1].x)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set y1="+str(points_in_streams[i][-1].y)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set x2="+str(points_in_streams[i][0].x)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set y2="+str(points_in_streams[i][0].y)+\
                    " where cat="+str(cat_of_line[i]))
      else:
        cur.execute("update streams set drainageArea_km2_1="+str(_A_point1)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set drainageArea_km2_2="+str(_A_point2)+\
                    " where cat="+str(cat_of_line[i]))
        # Points
        cur.execute("update streams set x1="+str(points_in_streams[i][0].x)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set y1="+str(points_in_streams[i][0].y)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set x2="+str(points_in_streams[i][-1].x)+\
                    " where cat="+str(cat_of_line[i]))
        cur.execute("update streams set y2="+str(points_in_streams[i][-1].y)+\
                    " where cat="+str(cat_of_line[i]))
      print i
    #streamsTopo.write()
    streamsTopo.table.conn.commit()
    streamsTopo.build()


############################
# PASS VARIABLES AND SOLVE #
############################
def main():
    """
    Links each river segment to the next downstream segment in a tributary 
    network by referencing its category (cat) number in a new column. "0"
    means that the river exits the map.
    """
    
    options, flags = grass.parser()

    ##########
    # SET-UP #
    ##########
   
    # This code is for 2D flexural isostasy
    flex = gflex.F2D()
    # And show that it is coming from GRASS GIS
    flex.grass = True
   
    # Method
    flex.Method = 'SAS_NG'
   
    # Parameters that are often changed for the solution
    ######################################################
   
    # x, y, q
    flex.x, flex.y = get_points_xy(options['input'])
    # xw, yw: gridded output
    if len(grass.parse_command('g.list', type='vect', pattern=options['output'])):
        if not grass.overwrite():
            grass.fatal("Vector map '" + options['output'] + "' already exists. Use '--o' to overwrite.")
    # Just check raster at the same time if it exists
    if len(grass.parse_command('g.list', type='rast', pattern=options['raster_output'])):
        if not grass.overwrite():
            grass.fatal("Raster map '" + options['raster_output'] + "' already exists. Use '--o' to overwrite.")
    grass.run_command('v.mkgrid', map=options['output'], type='point', overwrite=grass.overwrite(), quiet=True)
    grass.run_command('v.db.addcolumn', map=options['output'], columns='w double precision', quiet=True)
    flex.xw, flex.yw = get_points_xy(options['output']) # gridded output coordinates
    vect_db = grass.vector_db_select(options['input'])
    col_names = np.array(vect_db['columns'])
    q_col = (col_names == options['column'])
    if np.sum(q_col):
        col_values = np.array(vect_db['values'].values()).astype(float)
        flex.q = col_values[:, q_col].squeeze() # Make it 1D for consistency w/ x, y
    else:
        grass.fatal("provided column name, "+options['column']+" does not match\nany column in "+options['q0']+".")
    # Elastic thickness
    flex.Te = float(options['te'])
    if options['te_units'] == 'km':
        flex.Te *= 1000
    elif options['te_units'] == 'm':
        pass
    else:
        grass.fatal("Inappropriate te_units. How? Options should be limited by GRASS.")
    flex.rho_fill = float(options['rho_fill'])
   
    # Parameters that often stay at their default values
    ######################################################
    flex.g = float(options['g'])
    flex.E = float(options['ym']) # Can't just use "E" because reserved for "east", I think
    flex.nu = float(options['nu'])
    flex.rho_m = float(options['rho_m'])
    # Set verbosity
    if grass.verbosity() >= 2:
        flex.Verbose = True
    if grass.verbosity() >= 3:
        flex.Debug = True
    elif grass.verbosity() == 0:
        flex.Quiet = True
   
    # Check if lat/lon and let user know if verbosity is True
    if grass.region_env()[6] == '3':
        flex.latlon = True
        flex.PlanetaryRadius = float(grass.parse_command('g.proj', flags='j')['+a'])
        if flex.Verbose:
            print "Latitude/longitude grid."
            print "Based on r_Earth = 6371 km"
            print "Computing distances between load points using great circle paths"
    ##########
    # SOLVE! #
    ##########
    flex.initialize()
    flex.run()
    flex.finalize()
   
    # Now to use lower-level GRASS vector commands to work with the database
    # table and update its entries
    # See for help:
    # http://nbviewer.ipython.org/github/zarch/workshop-pygrass/blob/master/02_Vector.ipynb
    w = vector.VectorTopo(options['output'])
    w.open('rw') # Get ready to read and write
    wdb = w.dblinks[0]
    wtable = wdb.table()
    col = int((np.array(wtable.columns.names()) == 'w').nonzero()[0]) # update this column
    for i in range(1, len(w)+1):
        # ignoring 1st column: assuming it will be category (always true here)
        wnewvalues = w[i].attrs.values()[1:col] + tuple([flex.w[i-1]]) + w[i].attrs.values()[col+1:]
        wtable.update(key=i, values=wnewvalues)
    wtable.conn.commit() # Save this
    w.close(build=False) # don't build here b/c it is always verbose
    grass.run_command('v.build', map=options['output'], quiet=True)
   
    # And raster export
    # "w" vector defined by raster resolution, so can do direct v.to.rast
    # though if this option isn't selected, the user can do a finer-grained
    # interpolation, which shouldn't introduce much error so long as these
    # outputs are spaced at << 1 flexural wavelength.
    if options['raster_output']:
        grass.run_command('v.to.rast', input=options['output'], output=options['raster_output'], use='attr', attribute_column='w', type='point', overwrite=grass.overwrite(), quiet=True)
        # And create a nice colormap!
        grass.run_command('r.colors', map=options['raster_output'], color='differences', quiet=True)
def install_dependencies():
    print "PLACEHOLDER"
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--install-dependencies':
        install_dependencies()
    else:
        main()


