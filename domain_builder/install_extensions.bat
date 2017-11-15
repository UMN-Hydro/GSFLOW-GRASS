@echo off
FOR %%A in ( 
	v.gsflow.export 
	v.gsflow.gravres 
	v.gsflow.grid 
	v.gsflow.hruparams 
	v.gsflow.reaches 
	v.gsflow.segments 
	r.gsflow.hydrodem 
	v.stream.inbasin 
	v.stream.network 
	r.cell.area 
	r.stream.basins 
	r.hydrodem) DO (
	g.extension %%A
)
