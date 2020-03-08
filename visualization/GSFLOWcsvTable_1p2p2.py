# -*- coding: utf-8 -*-
"""
Created on Nov 18, 2019

@author: gcng
"""

# List of gsflow-csv-file variables, for GSFLOW 1.2.2

# see Table 1-6 in "GSFLOW Input Instructions: A Supplement to Appendix 1 of 
# the GSFLOW manual (USGS TM 6-D1), 
# Version 1.2.2 GSFLOW release, February 23, 2018"

# Note: variables commented out can be selected as output variables for PRMS 
# Statistic Variables File and PRMS Animation Files (among other variables, 
# see Table 1-5)


# Create cell with descriptions
varname = []
unit = []
descr = []

MODFLOW_len_unit = 'm' 
MODFLOW_time_unit = 'd' # should be day

varname.append('Date')
descr.append('Date Month, day, and year designation')
unit.append('Mon/Day/Yr')

varname.append('Precip_Q') # ****
descr.append('Volumetric flow rate of precipitation on modeled region')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)

varname.append('CapET_Q') 
descr.append('Volumetric flow rate of evapotranspiration from pervious areas') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('ImpervEvap_Q ') 
descr.append('Volumetric flow rate of evaporation from impervious areas') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('CanopyEvap_Q') 
descr.append('Volumetric flow rate of evaporation of intercepted precipitation') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('SnowEvap_Q') 
descr.append('Volumetric flow rate of snowpack sublimation') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('StreamOut_Q') # ****
descr.append('Volumetric flow rate of streamflow leaving modeled region') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('PotGravDrn2Unsat') 
descr.append('Potential volumetric flow rate of gravity drainage from the soil zone to the unsaturated zone (before conditions of the unsaturated and saturated zones are applied)') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('Sat2Grav_Q') # ****
descr.append('Volumetric flow rate of ground-water discharge from the saturated zone to the soil zone') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('NetBoundaryFlow2Sat_Q') # ****
descr.append('Volumetric flow rate to saturated zone along external boundary (negative value is flow out of modeled region)') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('StreamExchng2Sat_Q')  # ****
descr.append('Volumetric flow rate of stream leakage to the unsaturated and saturated zones') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('RechargeUnsat2Sat_Q') 
descr.append('Volumetric flow rate of recharge from the unsaturated zone to the saturated zone') 
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
#varname.append('basinseepout') # **** (same as basingw2sz)
#descr.append('Volumetric flow rate of ground-water discharge from the saturated zone to the soil zone') 
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('Unsat_S') 
descr.append('Volume of water in the saturated zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Unsat_S') 
descr.append('Volume of water in the unsaturated zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Cap_S') 
descr.append('Volume of water in capillary reservoirs of the soil zone') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Grav_S') 
descr.append('Volume of water in gravity reservoirs of the soil zone') 
unit.append(MODFLOW_len_unit + '^3')  

## NOT IN GSFLOW 1.2.2 
#varname.append('basingwstor') 
#descr.append('Volume of water in PRMS ground-water reservoirs (PRMS-only simulation)') 
#unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Canopy_S') 
descr.append('Volume of intercepted precipitation in plant-canopy reservoirs') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Imperv_S')
descr.append('Volume of water in impervious reservoir') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('SnowPweqv_S') 
descr.append('Volume of water in snowpack storage') 
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Interflow2Stream_Q') # ****
descr.append('Volumetric flow rate of slow interflow to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
#varname.append('Sroff2Stream_Q') # ****
#descr.append('Volumetric flow rate of the sum of Hortonian (Horton, 1933) and Dunnian surface runoff (Dunne and Black, 1970) to streams')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('Stream_S')
descr.append('Volume of water in streams')
unit.append(MODFLOW_len_unit + '^3')  
 
varname.append('Lake_S')
descr.append('Volume of water in lakes')
unit.append(MODFLOW_len_unit + '^3')  
 
#varname.append('obs_strmflow')
#descr.append('Volumetric flow rate of streamflow measured at a gaging station')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
#varname.append('UnsatDrainageExcess_Q')
#descr.append('Volumetric flow rate of gravity drainage from the soil zone not accepted due to the conditions in the unsaturated and saturated zones')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
#varname.append('Pref_S')
#descr.append('Volume of water stores in preferential-flow reservoirs of the soil zone')
#unit.append(MODFLOW_len_unit + '^3') 
 
varname.append('UnsatET_Q')
descr.append('Volumetric flow rate of evapotranspiration from the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)
 
varname.append('SoilDrainage2Unsat_Q') # ****
descr.append('Volumetric flow rate of gravity drainage to the unsaturated and saturated zones')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('Unsat_dS')
#descr.append('Change in unsaturated-zone storage')
#unit.append(MODFLOW_len_unit + '^3') 
 
#varname.append('net_sz2gw') # ****
#descr.append('Net volumetric flow rate of gravity drainage from the soilzone to the unsaturated and saturated zones')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('Sat_dS')
#descr.append('Change in saturated-zone storage')
#unit.append(MODFLOW_len_unit + '^3') 
 
#varname.append('Stream2Sat_Q') # ****
#descr.append('Volumetric flow rate of stream leakage to unsaturated and saturated zones')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
#varname.append('UnsatStream_dS')
#descr.append('Change in unsaturated-zone storage under streams')
#unit.append(MODFLOW_len_unit + '^3') 
 
#varname.append('SatDisch2Stream_Q') # ****
#descr.append('Volumetric flow rate of ground-water discharge to streams')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  

varname.append('StreamExchng2Sat_Q') # ****
descr.append('Volumetric flow rate of exchange betweeen streams and the unsaturated and saturated zones (value is equal to Stream2Sat_Q minus SatDisch2Stream_Q, where a negative value indicates a net loss from streams)')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
varname.append('UnsatStream_S')
descr.append('Volume of water in the unsaturated zone')
unit.append(MODFLOW_len_unit + '^3') 
 
#varname.append('Lake2Sat_Q')
#descr.append('Volumetric flow rate of lake leakage to the unsaturated and saturated zones')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
#varname.append('Lake_dS') 
#descr.append('Change in lake storage')
#unit.append(MODFLOW_len_unit + '^3') 
 
#varname.append('SatDisch2Lake_Q')
#descr.append('Volumetric flow rate of ground-water discharge to lakes')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit)  
 
varname.append('Infil2Soil_Q') # **** (vs. uz_infil: this one is before ET?)
descr.append('Volumetric flow rate of soil infiltration including precipitation, snowmelt, and cascading Hortonian flow')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('DunnSroff2Stream_Q')
descr.append('Volumetric flow rate of Dunnian runoff to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('HortSroff2Stream_Q')
descr.append('Volumetric flow rate of Hortonian runoff to streams')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('basinsm2gvr')
#descr.append('Volumetric flow rate of flow from capillary reservoirs to gravity reservoirs')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('basingvr2sm')
#descr.append('Volumetric flow rate of replenishment of capillary reservoirs from gravity reservoirs')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('basininfil_tot') 
#descr.append('Volumetric flow rate of soil infiltration into capillary reservoirs including precipitation, snowmelt, and cascading Hortonian and Dunnian runoff and interflow minus infiltration to preferential-flow reservoirs')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('Infil2Pref_Q')
#descr.append('Volumetric flow rate of soil infiltration into preferential-flow reservoirs including precipitation, snowmelt, and cascading surface runoff')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('DunnInterflow2Cap_Q')
#descr.append('Volumetric flow rate of cascading Dunnian runoff and interflow to HRUs')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('ActualET_Q') # ****
#descr.append('Volumetric flow rate of actual evapotranspiration from HRUs')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
#varname.append('SnowMelt_Q')
#descr.append('Volumetric flow rate of snowmelt')
#unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('HortSroff2Lake_Q')
descr.append('Volumetric flow rate of Hortonian runoff to lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('DunnInterflow2Lake_Q')
descr.append('Volumetric flow rate of Dunnian runoff and interflow to lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('LakeEvap_Q')
descr.append('Volumetric flow rate of evaporation from lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('LakePrecip_Q')
descr.append('Volumetric flow rate of precipitation on lakes')
unit.append(MODFLOW_len_unit + '^3/'+ MODFLOW_time_unit) 
 
varname.append('KKITER')
descr.append('Number of iterations for each time step')
unit.append('none')



