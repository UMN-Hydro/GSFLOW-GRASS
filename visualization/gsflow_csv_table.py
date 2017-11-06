# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 15:33:16 2017

@author: gcng
"""

# List of gsflow-csv-file variables

# see Table 12 of GSFLOW manual

# Create cell with descriptions
varname = []
unit = []
descr = []

MODFLOW_len_unit = 'm' 
MODFLOW_time_unit = 'd' # should be day

varname.append('Date')
descr.append('Date Month, day, and year designation')
unit.append('Mon/Day/Yr')

varname.append('basinppt') # ****
descr.append('Volumetric flow rate of precipitation on modeled region')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)

varname.append('basinpervet') 
descr.append('Volumetric flow rate of evapotranspiration from pervious areas') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinimpervevap ') 
descr.append('Volumetric flow rate of evaporation from impervious areas') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinintcpevap') 
descr.append('Volumetric flow rate of evaporation of intercepted precipitation') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinsnowevap') 
descr.append('Volumetric flow rate of snowpack sublimation') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinstrmflow') # ****
descr.append('Volumetric flow rate of streamflow leaving modeled region') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinsz2gw') 
descr.append('Potential volumetric flow rate of gravity drainage from the soil zone to the unsaturated zone (before conditions of the unsaturated and saturated zones are applied)') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basingw2sz') # ****
descr.append('Volumetric flow rate of ground-water discharge from the saturated zone to the soil zone') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('gw_inout') # ****
descr.append('Volumetric flow rate to saturated zone along external boundary (negative value is flow out of modeled region)') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('stream_leakage')  # ****
descr.append('Volumetric flow rate of stream leakage to the unsaturated and saturated zones') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('uzf_recharge') 
descr.append('Volumetric flow rate of recharge from the unsaturated zone to the saturated zone') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinseepout') # **** (same as basingw2sz)
descr.append('Volumetric flow rate of ground-water discharge from the saturated zone to the soil zone') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('sat_stor') 
descr.append('Volume of water in the saturated zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('unsat_stor') 
descr.append('Volume of water in the unsaturated zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basinsoilmoist') 
descr.append('Volume of water in capillary reservoirs of the soil zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basingravstor') 
descr.append('Volume of water in gravity reservoirs of the soil zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basingwstor') 
descr.append('Volume of water in PRMS ground-water reservoirs (PRMS-only simulation)') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basinintcpstor') 
descr.append('Volume of intercepted precipitation in plant-canopy reservoirs') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basinimpervstor')
descr.append('Volume of water in impervious reservoir') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basinpweqv') 
descr.append('Volume of water in snowpack storage') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('basininterflow') # ****
descr.append('Volumetric flow rate of slow interflow to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinsroff') # ****
descr.append('Volumetric flow rate of the sum of Hortonian (Horton, 1933) and Dunnian surface runoff (Dunne and Black, 1970) to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('strm_stor')
descr.append('Volume of water in streams')
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('lake_stor')
descr.append('Volume of water in lakes')
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('obs_strmflow')
descr.append('Volumetric flow rate of streamflow measured at a gaging station')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinszreject')
descr.append('Volumetric flow rate of gravity drainage from the soil zone not accepted due to the conditions in the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('basinprefstor')
descr.append('Volume of water stores in preferential-flow reservoirs of the soil zone')
unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('uzf_et')
descr.append('Volumetric flow rate of evapotranspiration from the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('uzf_infil') # ****
descr.append('Volumetric flow rate of gravity drainage to the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('uzf_del_stor')
descr.append('Change in unsaturated-zone storage')
unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('net_sz2gw') # ****
descr.append('Net volumetric flow rate of gravity drainage from the soilzone to the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('sat_change_stor')
descr.append('Change in saturated-zone storage')
unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('streambed_loss') # ****
descr.append('Volumetric flow rate of stream leakage to unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
varname.append('sfruz_change_store')
descr.append('Change in unsaturated-zone storage under streams')
unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('gwflow2strms') # ****
descr.append('Volumetric flow rate of ground-water discharge to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
varname.append('sfruz_tot_stor')
descr.append('Volume of water in the unsaturated zone')
unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('lakebed_loss')
descr.append('Volumetric flow rate of lake leakage to the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
varname.append('lake_change_stor') 
descr.append('Change in lake storage')
unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('gwflow2lakes')
descr.append('Volumetric flow rate of ground-water discharge to lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
varname.append('basininfil') # **** (vs. uz_infil: this one is before ET?)
descr.append('Volumetric flow rate of soil infiltration including precipitation, snowmelt, and cascading Hortonian flow')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basindunnian')
descr.append('Volumetric flow rate of Dunnian runoff to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinhortonian')
descr.append('Volumetric flow rate of Hortonian runoff to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinsm2gvr')
descr.append('Volumetric flow rate of flow from capillary reservoirs to gravity reservoirs')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basingvr2sm')
descr.append('Volumetric flow rate of replenishment of capillary reservoirs from gravity reservoirs')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basininfil_tot') 
descr.append('Volumetric flow rate of soil infiltration into capillary reservoirs including precipitation, snowmelt, and cascading Hortonian and Dunnian runoff and interflow minus infiltration to preferential-flow reservoirs')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basininfil2pref')
descr.append('Volumetric flow rate of soil infiltration into preferential-flow reservoirs including precipitation, snowmelt, and cascading surface runoff')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basindnflow')
descr.append('Volumetric flow rate of cascading Dunnian runoff and interflow to HRUs')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinactet') # ****
descr.append('Volumetric flow rate of actual evapotranspiration from HRUs')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinsnowmelt')
descr.append('Volumetric flow rate of snowmelt')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinhortonianlakes')
descr.append('Volumetric flow rate of Hortonian runoff to lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinlakeinsz')
descr.append('Volumetric flow rate of Dunnian runoff and interflow to lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinlakeevap')
descr.append('Volumetric flow rate of evaporation from lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('basinlakeprecip')
descr.append('Volumetric flow rate of precipitation on lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('kkiter')
descr.append('Number of iterations for each time step')
unit.append('none')



