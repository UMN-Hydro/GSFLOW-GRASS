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
