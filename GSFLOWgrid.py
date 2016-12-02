# "value" is the cat for the 
r.stream.basins direction=draindir stream_rast=streams basins=basins --o
r.to.vect in=basins out=basins type=area --o

g.region
v.mkgrid grid=50,50 map=grid --o

# basins
v.overlay ainput=basins binput=grid output=tmp op=and --o
# clean up with some v.db.renamecolumn commands

# streams
v.extract in=streams out=tmp2 type=line
v.overlay ainput=tmp2 atype=line binput=boxtmp output=tmp3 op=and --o
# And then sort from top to bottom by:
# (1) finding segments with proper a_cat
# (2) finding which one ends with x_1, y_1
# (3) following up downriver until reaching x_2, y_2
# (4) while doing (3), make sure that all reaches are organized
#     upstream to downstream
# (5) Once all reaches are ordered, give them numbers. This may require
#     adding a column as step (1.5).



# UPDATE DECEMBER 1ST

# Basins 3x resolution
g.region rast=topo -p
g.region -p nsres=92.44186047 ewres=92.87925695999999
# Actually, if it isn't going to be perfect -- just 90x90
# Or no -- 200
g.region -pas res=200
r.mapcalc "topogrid = topo" --o # topogrid is output topo on output grid

# Now make grid
v.mkgrid map=grid --o

# And output
r.out.ascii in=topogrid out=topo.asc
r.out.ascii in=basin out=basinmask.asc null_value=0 --o


# Now overlay basins
# Automatically comes with b_row and b_col from the grid
v.overlay ainput=basins binput=grid output=tmp op=and --o

