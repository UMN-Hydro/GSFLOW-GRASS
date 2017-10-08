g.region -p rast=DEM
g.region -p vect=grid_tmp
v.to.rast in=streams_tmp out=streams_tmp use=val val=1 --o
r.mapcalc "streams_tmp = streams_tmp * DEM" --o
g.region -p rast=DEM_coarse
r.resamp.stats in=streams_tmp out=streams_tmp_coarse method=minimum
r.resamp_stats(input=raster_input, output=raster_output, method='average', overwrite=True)
r.patch in=streams_tmp_coarse,DEM_coarse out=DEM_coarse --o
