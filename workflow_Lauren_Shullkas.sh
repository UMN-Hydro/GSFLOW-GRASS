DEM_orig=DEM_orig
DEM=DEM
DEM_coarse=DEM_coarse # MODFLOW resolution
accumulation=accumulation_tmp
streams=streams_tmp
streams_onebasin=${streams}_onebasin
basins=basins_tmp
basins_onebasin=${basins}_onebasin
segments=segments_tmp
reaches=reaches_tmp
threshold=1000000 #m2 - drainage area
grid_res=500 #150 #1000 #m2 - for MODFLOW
grid=grid_tmp
slope=slope_tmp
aspect=aspect_tmp
HRUs=HRUs_tmp
gravity_reservoirs=gravity_reservoirs_tmp
basin_mask=basin_mask_tmp
basin_mask_out=basin_mask
pour_point=pp_tmp
icalc=1 # how to compute hydraulic geometry

# Set region
g.region -p rast=$DEM_orig

# Build flow accumulation with only fully on-map flow
r.cell.area output=cellArea_meters2 units=m2 --o
r.hydrodem in=$DEM_orig out=$DEM -a --o
r.watershed elevation=$DEM flow=cellArea_meters2 accumulation=$accumulation -m --o
r.mapcalc "${accumulation}_pos = $accumulation * ($accumulation > 0)" --o
r.null map=${accumulation}_pos setnull=0
r.mapcalc "${DEM}_pos_accum = $DEM * (isnull(${accumulation}_pos) == 0)" --o
r.null map=${DEM}_pos_accum setnull=0
r.mapcalc "${accumulation}_pos = ${accumulation}_pos * ($DEM > 0)" --o
r.null map=${accumulation}_pos setnull=0
# Repeat is sometimes needed
r.mapcalc "${DEM}_pos_accum = $DEM * (isnull(${accumulation}_pos) == 0)" --o
r.null map=${DEM}_pos_accum setnull=0
r.mapcalc "${accumulation}_pos = ${accumulation}_pos * ($DEM > 0)" --o
r.null map=${accumulation}_pos setnull=0

# Build streams and sub-basins
r.stream.extract elevation=${DEM}_pos_accum accumulation=${accumulation}_pos stream_raster=$streams stream_vector=$streams threshold=$threshold direction=draindir_tmp d8cut=0 --o
r.stream.basins direction=draindir_tmp stream_rast=$streams basins=$basins --o
r.to.vect input=$basins output=$basins type=area -v --o

# Build stream network
v.stream.network map=$streams

# Restrict to a single basin
basin_outlet_cat=2638 #2840 #2485 # You must find and define this after building the stream network
v.stream.inbasin input_streams=$streams input_basins=$basins output_streams=$streams_onebasin output_basin=$basins_onebasin cat=$basin_outlet_cat output_pour_point=$pour_point --o

# GSFLOW segments: sections of stream that define subbasins
v.gsflow.segments input=$streams_onebasin output=$segments icalc=$icalc --o

# MODFLOW grid & basin mask (1s where basin exists and 0 where it doesn't)
v.gsflow.grid basin=$basins_onebasin  pour_point=$pour_point raster_input=$DEM dx=$grid_res dy=$grid_res output=$grid mask_output=$basin_mask raster_output=$DEM_coarse --o

# GSFLOW reaches: intersection of segments and grid
v.gsflow.reaches segment_input=$segments grid_input=$grid elevation=$DEM output=$reaches --o

# GSFLOW HRU parameters
r.slope.aspect elevation=$DEM slope=$slope aspect=$aspect format=percent zscale=0.01 --o
v.gsflow.hruparams input=$basins_onebasin elevation=$DEM output=$HRUs slope=$slope aspect=$aspect --o

# GSFLOW gravity reservoirs
v.gsflow.gravres hru_input=$HRUs grid_input=$grid output=$gravity_reservoirs --o

# Export DEM with MODFLOW resolution
# Also export basin mask -- 1s where basin exists and 0 where it doesn't
v.to.rast in=$basins_onebasin out=$basin_mask use=val val=1 --o
g.region vect=$grid res=$grid_res
r.resamp.stats in=$basin_mask out=$basin_mask method=sum --o
r.mapcalc "$basin_mask = $basin_mask > 0" --o
#g.region rast=$DEM
r.out.ascii input=$DEM output=$DEM.asc null_value=0 --o
r.out.ascii input=$basin_mask output=$basin_mask_out.asc null_value=0 --o
g.region rast=$DEM
#g.region vect=$basins_onebasin

# Export tables and discharge point
v.gsflow.export reaches_input=$reaches segments_input=$segments gravres_input=$gravity_reservoirs hru_input=$HRUs pour_point_input=$pour_point reaches_output=$reaches segments_output=$segments gravres_output=$gravity_reservoirs hru_output=$HRUs pour_point_output=$pour_point --o

