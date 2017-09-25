# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 22:21:56 2017

Based on: GSFLOW_print_PRMSparamfile4.m

@author: gcng
"""
import numpy as np # matlab core
import scipy as sp # matlab toolboxes
import matplotlib.pyplot as plt # matlab-like plots
import os  # os functions
import pandas as pd # for data structures and reading in data from text file
from ConfigParser import SafeConfigParser

parser = SafeConfigParser()
parser.read('settings.ini')
LOCAL_DIR = parser.get('settings', 'local_dir')

GSFLOW_DIR = LOCAL_DIR + "/GSFLOW"

# GSFLOW_print_PRMSparamfile4.m
# 11/23/16
# (based on PRMS_print_paramfile3.m, 11/25/15)
# gcng
#
# This version is NOT for climate_hru (climate_hru is for pre-processed time
# series data)
#
# v4 - includes cascade

# Creates inputs files with PRMS parameters (**not GSFLOW parameters***)
# leaves out many of the "extra" options.  
# GSFLOW Input files:
#   - control file (generate with GSFLOW_print_controlfile4*.m)
#   - parameter files (generate with GSFLOW_print_PRMSparamfile4.m)
#   - variable / data files (generate with GSFLOW_print_ObsMet_files1.m)
# (Control file includes names of parameter and data files.  Thus, 
# parameter and variable file names must match those specified there!!)
#
# note on order of parameter values when there are 2 dimensions: 
# par_value[ii] = (ndim1, ndim2), par_value[ii] = par_value[ii](:)
#
# Search for 'CHANGE FOR SPECIFIC SITE
# Search for 'CHANGE FOR SPECIFIC SITE - CHIMBORAZO'

### Parameter file

# See Appendix 1 Table 1-1 (p.30-31) for Dimensions list, Appendix 1 Table 1-
# 3 for parameters (by module) (p.37-59), description p. 128

# general syntax - Dimensions
#   line 1: '####'
#   line 2: dimension name
#   line 3: dimension value


# general syntax - Parameters
#   line 1: '####'
#   line 2: parameter name
#   line 3: Number of dimensions
#   line 4 to 3+NumOfDims: Name of dimensions, 1 per line 
#   line 3+NumOfDims+1: Number of values 
#   line 3+NumOfDims+2: data type -> 1=int, 2=single prec, 3=double prec, 4=char str
#   line 3+NumOfDims+3 to end: parameter values, 1 per line

# 2-dimensional arrays:
#   read in column-by-column, i.e. [nhru x nmonth] has all Jan first, then 
#   all Feb, etc.

# *** CUSTOMIZE TO YOUR COMPUTER! *****************************************
# NOTE: '/' is directory separator for Linux, '\' for Windows!!

# directory for GSFLOW input and output files (include slash ('/') at end)
# (creates directories if they don't exist)
PRMSinput_dir = GSFLOW_DIR + '/inputs/PRMS/'
PRMSoutput_dir = GSFLOW_DIR + '/outputs/'
# 
# PRMSinput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/ChimTest/')
# PRMSoutput_dir = '')


# directory with files to be read in to generate PRMS input files (include slash ('/') at end)
in_data_dir = GSFLOW_DIR + '/DataToReadIn/'
in_GISdata_dir = in_data_dir + 'GIS/' # specifically GIS data

# parameter file that will be written (name must match that in Control file!)
parfil_pre = 'ChimTest'
fl_veg_shift = 0 # shift veg upslope (for pt_alpha)
if fl_veg_shift == 1:
    parfil_pre = parfil_pre + 'VegShift'

# GIS-generated files read in to provide values to PRMS input file
HRUfil = in_GISdata_dir  + 'HRU.csv'
segmentfil = in_GISdata_dir  + 'segment_data_4A_INFORMATION.txt'
reachfil = in_GISdata_dir  + 'reach_data.txt' # *** NEW FOR GSFLOW, only required here for NREACH, other info is for SFR file
gvrfil = in_GISdata_dir  + 'GRAVITY_RESERVOIRS.csv' # *** NEW FOR GSFLOW
MODfil = in_GISdata_dir  + 'basinmask.asc' # *** NEW FOR GSFLOW, only required here for ngwcell (NROW*NCOL)
# *************************************************************************

#%%
# Project-specific entries ->

# - choose one:
# model_mode = 'WRITE_CLIMATE' # creates pre-processed climate_hru files
# model_mode = 'PRMS' # run only PRMS
# model_mode = 'MODFLOW' # run only MODFLOW-2005
model_mode = 'GSFLOW' # run coupled PRMS-MODFLOW

parfil = PRMSinput_dir + 'py_' + parfil_pre + '_' + model_mode + '.param'

# Load GIS-generated files with HRU, segment etc information
# (pandas work like structures in matlab?)
segmentdata = pd.read_csv(segmentfil)
HRUdata = pd.read_csv(HRUfil)
reachdata = pd.read_csv(reachfil)
gvrdata = pd.read_csv(gvrfil)

f = file(MODfil, 'r')
MODdata = {}
for i in range(6):
    line = f.readline()
    line = line.rstrip() # remove newline characters
    key, value = line.split(': ')
    try:
      value = int(value)
    except:
      value = float(value)
    MODdata[key] = value
f.close()


# 2 lines available for comments
title_str1 = 'TEST'
title_str2 = 'much based on merced and sagehen examples, NOT using climate_hru'

#%%
# n_Dim, n_par, should be dynamically generated

# - initialize dimension and parameter variables
dim_name = []
dim_value = []

par_name = [] # name of parameter
par_num_dim = [] # number of dimensions
par_dim_name = [] # names of dimensions
par_num = [] # number of parameter values
par_type = [] # parameter number type 1=int, 2=single prec, 3=double prec, 4=char str
par_value = [] # array of parameter values

# ********************
# *    Dimensions    *
# ********************
# -- Fixed Dimensions (generally, do not change) --
dim_name.append('ndays') # max num days in year
dim_value.append(366)

dim_name.append('nmonths') # num months in year
dim_value.append(12)

dim_name.append('one')
dim_value.append(1)


# -- Spatial dimensions --
# ****to be read in from GIS info****
dim_name.append('nhru')  
dim_value.append(max(HRUdata.id))

# ****to be read in from GIS info****
dim_name.append('nsegment')  # num stream channel segments
# dim_value.append(max(segmentdata.data(:, strcmp(segmentdata.colheaders, 'id'))))
dim_value.append(max(segmentdata.NSEG))

dim_name.append('ngw')  # num stream channel segments
ind = np.squeeze(np.where(np.array(dim_name) == 'nhru'))
dim_value.append(dim_value[ind])  # do not change

dim_name.append('nssr')  # num subsurface reservoirs
ind = np.squeeze(np.where(np.array(dim_name) == 'nhru'))
dim_value.append(dim_value[ind])  # do not change

# new for GSFLOW:
if model_mode == 'GSFLOW':
    # ****to be read in from GIS info****
    dim_name.append('nreach')  # total num reaches (in all stream segments)
    dim_value.append(reachdata.shape[0])

# new for GSFLOW:
# ****to be read in from GIS info****
dim_name.append('ngwcell')  # total num MODFLOW GW cells 
dim_value.append(MODdata['rows'] * MODdata['cols'])

# new for GSFLOW:
# ****to be read in from GIS info****
dim_name.append('nhrucell')  # total num gravity reservoirs (intersections of hru and grid cells)
dim_value.append(max(gvrdata.cat))

# -- Time-series input data dimensions --
# (some of these data are not needed but are handy to output for calibration)
dim_name.append('nrain')  # num precip measurement stations
dim_value.append(0) # can be 0 if using climate_hru

dim_name.append('ntemp')  # num temperature measurement stations
dim_value.append(0) # can be 0 if using climate_hru

dim_name.append('nsol')  # num solar rad measurement stations
dim_value.append(0) # can be 0 if using ddsolrad (degree-day) or climate_hru

dim_name.append('nhumid')  # num humidity measurement stations
dim_value.append(0) # optional, for some  ET modules

dim_name.append('nwind')  # num wind speed measurement stations
dim_value.append(0) # optional, for some  ET modules

dim_name.append('nobs')  # num streamflow obs, replaces model when using muskingum or strmflow_in_out
dim_value.append(0)

dim_name.append('nsnow')  # num snow depth measurement stations
dim_value.append(0)  # optional


# -- Computational dimensions --
# ****to be read in from GIS info****
dim_name.append('ncascade')  # num HRU links, for cascade_flag = 1
ind = np.squeeze(np.where(np.array(dim_name) == 'nhru'))
dim_value.append(dim_value[ind])  # do not change

# (gw - can mirror surface)
dim_name.append('ncascdgw')  # num GWR links, for cascadegw_flag = 1
ind = np.squeeze(np.where(np.array(dim_name) == 'ncascade'))
dim_value.append(dim_value[ind])  # do not change

# (snow)
dim_name.append('ndepl')  # num snow-depletion curves
dim_value.append(2)

dim_name.append('ndeplval')  # num values in all snow-depletion curves
ind = np.squeeze(np.where(np.array(dim_name) == 'ndepl'))
dim_value.append(dim_value[ind] * 11) # do not change

NumDims = len(dim_name)
#%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ********************
# *   Parameters    *
# ********************

# -- Basic Computational Attributes --
par_name.append('elev_units')
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(1)  # 0: ft, 1: m) do not change

# ****to be read in from GIS info****
par_name.append('hru_area')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(HRUdata[par_name[-1]])  # [acre]

# # new to GSFLOW
# par_name.append('basin_area')
# par_dim_name.append('one')
# par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
#ind = np.squeeze(np.where(np.array(par_name) == 'hru_area'))
# par_value.append(sum(par_value[ind]))  # [acre]

# ****to be read in from GIS info****
par_name.append('hru_aspect')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(HRUdata[par_name[-1]])  # [angular degrees]

# ****to be read in from GIS info****
par_name.append('hru_elev')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
hru_elev = HRUdata[par_name[-1]] # this var gets used alot
par_value.append(hru_elev)  # elev_units = m

# ****to be read in from GIS info****
par_name.append('hru_lat')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(HRUdata[par_name[-1]])  # [angular degrees]

# ****to be read in from GIS info****
par_name.append('hru_slope')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(HRUdata[par_name[-1]])  # [-]

# -- GSFLOW: PRMS-MODFLOW mapping
ind1 = np.squeeze(np.where(np.array(dim_name) == 'nhru'))
ind2 = np.squeeze(np.where(np.array(dim_name) == 'ngwcell'))
if (model_mode == 'GSFLOW' or model_mode == 'PRMS') and (dim_value[ind1] != dim_value[ind2]):
    # new for GSFLOW, also for model_mode=PRMS which nhru~=ngwcell
    par_name.append('gvr_hru_id') # HRU corresponding to gravity reservior
    par_dim_name.append('nhrucell')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(gvrdata[par_name[-1]])  # [-]

    # new for GSFLOW, also for model_mode=PRMS which nhru~=ngwcell
    par_name.append('gvr_hru_pct') # Decimal fraction of hru area associated with gravity reservior
    par_dim_name.append('nhrucell')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(gvrdata[par_name[-1]]/100)  # convert from percent to frac [-]
    
    # new for GSFLOW, also for model_mode=PRMS which nhru~=ngwcell
    par_name.append('gvr_cell_id') # HRU corresponding to gravity reservior
    par_dim_name.append('nhrucell')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(gvrdata[par_name[-1]])  # [-]

    # new for GSFLOW, also for model_mode=PRMS which nhru~=ngwcell
    par_name.append('gvr_cell_pct') # Decimal fraction of cell area associated with gravity reservior
    par_dim_name.append('nhrucell')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(gvrdata[par_name[-1]]/100)  # convert from percent to frac [-]


# -- Measured input --
par_name.append('precip_units')
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0)  # 0: inch, 1: mm) do not change

par_name.append('temp_units')
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0)  # 0: F, 1: C) do not change

par_name.append('adj_by_hru') # flag for how to adjust precip and temp 
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(1) # adjust by: 0= sub-basin, 1= hru) do not change

par_name.append('adjmix_rain') # monthly factor to adjust rain proportion into rain/snow event
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(np.ones((12,1)))  # [-]) do not change(?)  double-check what this means??

par_name.append('rain_cbh_adj') # monthly factor to adjust rain proportion into rain/snow event
par_dim_name.append(['nhru', 'nmonths'])
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_num.append(1)
for x in par_dim_name[-1]:
    ind = np.squeeze(np.where(np.array(dim_name) == x))
    par_num.append(par_num[-1] * dim_value[ind])
par_value.append(np.ones((par_num[-1], 1)))  # [-]

par_name.append('snow_cbh_adj') # monthly factor to adjust data
par_dim_name.append(['nhru', 'nmonths'])
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_num.append(1)
for x in par_dim_name[-1]:
    ind = np.squeeze(np.where(np.array(dim_name) == x))
    par_num.append(par_num[-1] * dim_value[ind])
par_value.append(np.ones((par_num[-1], 1)))  # [-]

par_name.append('tmax_cbh_adj') # monthly factor to adjust data) 
par_dim_name.append(['nhru', 'nmonths'])
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_num.append(1)
for x in par_dim_name[-1]:
    ind = np.squeeze(np.where(np.array(dim_name) == x))
    par_num.append(par_num[-1] * dim_value[ind])
par_value.append(0*np.ones((par_num[-1], 1)))  # temp_units = F

par_name.append('tmin_cbh_adj') # monthly factor to adjust data) 
par_dim_name.append(['nhru', 'nmonths'])
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_num.append(1)
for x in par_dim_name[-1]:
    ind = np.squeeze(np.where(np.array(dim_name) == x))
    par_num.append(par_num[-1] * dim_value[ind])
par_value.append(0*np.ones((par_num[-1], 1)))  # temp_units = F

par_name.append('tmax_allrain') # if tmax > tmax_allrain, then all rain
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(np.ones((12,1))*35)  # temp_units = F

par_name.append('tmax_allsnow')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(32)  # temp_units = F) do not change


# -- Solar Radiation --
# (** No parameters here when using solrad_module = climate_hru **)
par_name.append('dday_slope') # frac of PET that is sublimated from snow 
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.35* np.ones((dim_value[ind],1)))  # [-], default: 0.4, approx sagehen: 0.35

par_name.append('dday_intcp') # transmission coeff for short-wave radiation through winter canopy
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(-17* np.ones((dim_value[ind],1)))  # [-], default: 0.4, approx sagehen: 0.35

par_name.append('radadj_slope') # frac of PET that is sublimated from snow 
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.)  # [-], sagehen

par_name.append('radadj_intcp') # transmission coeff for short-wave radiation through winter canopy
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(1* np.ones((dim_value[ind],1)))  # [-], default: 0.4, approx sagehen: 0.35

par_name.append('ppt_rad_adj') # transmission coeff for short-wave radiation through winter canopy
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.2* np.ones((dim_value[ind],1)))  # [-], default: 0.4, approx sagehen: 0.35

par_name.append('tmax_index') # transmission coeff for short-wave radiation through winter canopy
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(50* np.ones((dim_value[ind],1)))  # [-], default: 0.4, approx sagehen: 0.35


# -- Potential ET distribution --
# ** No parameters here when using et_module = climate_hru **
# ** If et_module is pet_pt (Priestly-Taylor): **
# ** CHANGE FOR SPECIFIC SITE - CHIMBORAZO ***
# bare soil pt_alpha: 0.35 [Khaldi et al. 2014] to 1.04 [Flint and Childsb
#   1991, AFM]) check McMahon et al., 2013 HESS Table S8
x,hrui_melt = hru_elev.max(0),hru_elev.argmax(0)

veg_thresh = 4400 # m (Rachel's email 12/8/16 8:33am), approx veg elev line
veg_thresh_shift = veg_thresh + 500 # Morueta-Holme et al. 2015: shift >500m since 1802
par_name.append('pt_alpha') # alpha for Priestly-Taylor, ave ~1.26 but often higher in dry and windy regions (up to 2.14 or 2.47 - but this is open water)
# Crystal had 1.7 as a windy guess
par_dim_name.append(['nhru', 'nmonths'])
par_num.append(1)

for x in par_dim_name[-1]:
    ind = np.squeeze(np.where(np.array(dim_name) == x))
    par_num.append(par_num[-1] * dim_value[ind])
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
pt_alpha = 1.26*np.ones((par_num[-1],1))  # [-], default: 1.26
if fl_veg_shift == 1: 
    ind = np.squeeze(np.where(hru_elev > veg_thresh_shift))
    pt_alpha[ind] = 1.0
else:
    ind = np.squeeze(np.where(hru_elev > veg_thresh))
    pt_alpha[ind] = 1.0
pt_alpha[hrui_melt] = 0.75
par_value.append(pt_alpha)  # [-]


# -- Evapotranspiration and sublimation --
par_name.append('potet_sublim') # frac of PET that is sublimated from snow 
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.75)  # [-], default: 0.5, mercd: 0.101, sagehen: 0.75

par_name.append('rad_trncf') # transmission coeff for short-wave radiation through winter canopy
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.5 * np.ones((dim_value[ind],1)))  # default 0.5
par_value.append(0.23 * np.ones((dim_value[ind],1)))  # approx sagehen

par_name.append('soil_type') # only for ET calc, 1: sand, 2: loam, 3: clay
par_dim_name.append('nhru')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(np.ones((dim_value[ind],1)))  # merced: all 1's
par_value.append(2*np.ones((dim_value[ind],1), int))  # sagehen: mostly 2's, some 1's

# Below: Assumes transp_module=transp_tindex (phenolgy based on temp)
#   Start summing max air temp on month 'transp_beg') when sum >=
#   'transp_tmax', transp begins) transp ends on month 'transp_end' 
#   (previous month is last one with transp > 0)
# *** CHANGE FOR SPECIFIC SITE
par_name.append('transp_beg') # month to start summing max air temp
par_dim_name.append('nhru')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(1*np.ones((dim_value[ind],1), int))  # merced: all 4's, 1 for transp anytime - CHIMBORAZO

# *** CHANGE FOR SPECIFIC SITE
par_name.append('transp_end') # month to stop transp, so previous month is last with transp
par_dim_name.append('nhru')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(13*np.ones((dim_value[ind],1), int))  # max is 13, for transp anytime - CHIMBORAZO
#par_value.append(12*np.ones((dim_value[ind],1)))  # max is 12??? I think GSFLOW manual is wrong about this

# *** CHANGE FOR SPECIFIC SITE
par_name.append('transp_tmax') # max summed temp to start transp
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
# Has ET case
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(500*np.ones((dim_value[ind],1)))  # mercd: all 500 degF
par_value.append(np.zeros((dim_value[ind],1)))  # 0 degF for transp anytime - CHIMBORAZO
# No ET case
#par_value.append(5E9*ones(dim_value(strcmp(dim_name, par_dim_name{ii})),1))

# -- Interception --
# *** CHANGE FOR SPECIFIC SITE
par_name.append('cov_type') # 0=bare soil, 1=grasses, 2=shrubs, 3=trees, 4=coniferous
                            # (0 and 4 may not be options for GSFLOW)
par_dim_name.append('nhru')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(2*np.ones((dim_value[ind],1), int))  

# *** CHANGE FOR SPECIFIC SITE
par_name.append('covden_sum') # summer veg cover density [0.0 to 1.0] (for canopy interception)
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.8*np.ones((dim_value[ind],1)))  

# *** CHANGE FOR SPECIFIC SITE
par_name.append('covden_win') # winter veg cover density
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.8*np.ones((dim_value[ind],1)))  

par_name.append('snow_intcp') # snow interception storage capacity
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.07*np.ones((dim_value[ind],1)))  # [inch] merced
par_value.append(0.1*np.ones((dim_value[ind],1)))  # [inch] sagehen

par_name.append('srain_intcp') # summer rain interception storage capacity
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.04*np.ones((dim_value[ind],1)))  # [inch] merced
par_value.append(0.05*np.ones((dim_value[ind],1)))  # [inch] sagehen

par_name.append('wrain_intcp') # winter rain interception storage capacity
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.04*np.ones((dim_value[ind],1)))  # [inch] 

# -- Snow computations (based on mercd) --
par_name.append('albset_rna')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.8)  

par_name.append('albset_rnm')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.6)  

par_name.append('albset_sna')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.05)  

par_name.append('albset_snm')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.2)  

par_name.append('cecn_coef') 
par_dim_name.append('nmonths')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
# par_value.append(np.ones((12,1)) * 0.04876140267863427) # merced
par_value.append(np.ones((12,1)) * 5) # sagehen

par_name.append('den_init')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.1)  

par_name.append('den_max')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.6)  

par_name.append('emis_noppt')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.76)  

par_name.append('freeh2o_cap')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
# par_value.append(0.19) # merced
par_value.append(0.05) # sagehen

par_name.append('hru_deplcrv')
par_dim_name.append('nhru')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(np.ones((dim_value[ind],1), int))  

par_name.append('melt_force')
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(90)

par_name.append('melt_look')
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(90)

par_name.append('settle_const')
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(0.1)

# each snow depletion curve has 11 values in increasing order ([0, 1.0])) 
# should have ndeplval=(11*ndepl) values. Below snarea_curves values are 
# used in both mercd and acf examples.
par_name.append('snarea_curve')
par_dim_name.append('ndeplval')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append([0.05, 0.24, 0.4, 0.53, 0.65, 0.74, 0.82, 0.88, 0.93, 0.97, 1.0,
    0.05, 0.25, 0.4, 0.48, 0.54, 0.58, 0.61, 0.64, 0.66, 0.68, 0.7])  

# snarea_thresh is max threshold SWE below which snow depletion curve is 
# applied (using snarea_curve). Based on mercd and acf, value could be ~1
# of hru_elev.
par_name.append('snarea_thresh')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(par_name) == 'hru_elev'))
hru_elev = par_value[ind] 
par_value.append(0.01 * hru_elev)  

par_name.append('tstorm_mo') # 0: frontal, 1: convective
par_dim_name.append('nmonths')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
# par_value.append([0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0])  # merced
par_value.append(np.zeros((12,1), int))


# -- Hortonian surface runoff, infiltration, and impervious storage --
par_name.append('carea_max')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.6*np.ones((dim_value[ind],1)))  

par_name.append('hru_percent_imperv')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(np.zeros((dim_value[ind],1)))  # [-]

par_name.append('imperv_stor_max') # max impervious area retention storage
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(np.zeros((dim_value[ind],1)))  # [-]

par_name.append('snowinfil_max')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(2*np.ones((dim_value[ind],1)))  # merced
par_value.append(2.75*np.ones((dim_value[ind],1)))  # sagehen

# Assumes srunoff_module=runoff_smidx (or equivalently runoff_smidx_casc)
par_name.append('smidx_coef')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.0011346339705706754*np.ones((dim_value[ind],1)))  # merced
par_value.append(3.7e-4*np.ones((dim_value[ind],1)))  # sagehen

par_name.append('smidx_exp')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.21845105595154202*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.3*np.ones((dim_value[ind],1)))  # sagehen


# -- Soil zone storage, interflow, gravity drainage, dunnian surface runoff --

# - GSFLOW computational parameter
par_name.append('mxsziter') # max iterations for soil zone computations
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(25)  # default is 15, sagehen example uses 25

par_name.append('mnsziter') # min iterations for soil zone computations
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(7)  # sagehen example uses 7

# - GSFLOW computational parameter
par_name.append('szconverge') # convergence criteria for soil-zone
par_dim_name.append('one')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(1e-4)  # [inches] default is 1e-8, sagehen example uses 1e-4

par_name.append('fastcoef_lin')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.8034819756283689*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.4*np.ones((dim_value[ind],1)))  # sagehen

par_name.append('fastcoef_sq')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.740714266649265*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.8*np.ones((dim_value[ind],1)))  # sagehen

par_name.append('pref_flow_den')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.1*np.ones((dim_value[ind],1)))  # sagehen

par_name.append('sat_threshold')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(999*np.ones((dim_value[ind],1)))  # 999 is inifinite soil water
par_value.append(4.5*np.ones((dim_value[ind],1)))  # approx sagehen

par_name.append('slowcoef_lin')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.015*np.ones((dim_value[ind],1)))  

par_name.append('slowcoef_sq')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.1*np.ones((dim_value[ind],1)))  

par_name.append('soil_moist_init')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(1.5*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.1*np.ones((dim_value[ind],1)))  # sagehen

par_name.append('soil_moist_max')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(3*np.ones((dim_value[ind],1)))  # merced
par_value.append(3.5*np.ones((dim_value[ind],1)))  # approx sagehen

par_name.append('soil_rechr_init')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(1*np.ones((dim_value[ind],1)))  # merced
par_value.append(0*np.ones((dim_value[ind],1)))  # sagehen

par_name.append('soil_rechr_max')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(2.7*np.ones((dim_value[ind],1)))  # merced
par_value.append(2.*np.ones((dim_value[ind],1)))  # approx sagehen

par_name.append('soil2gw_max')
par_dim_name.append('nhru')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
if model_mode == 'PRMS':
    par_value.append(0.4979401249572195*np.ones((dim_value[ind],1)))  # merced
else:
    par_value.append(0.*np.ones((dim_value[ind],1)))  #

par_name.append('ssr2gw_exp')
par_dim_name.append('nssr')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(1*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.75*np.ones((dim_value[ind],1)))  # sagehen

# sagehen - PRMS: 0.0378 to 0.05571, GSFLOW: 0.189 to 0.278 
par_name.append('ssr2gw_rate') # sagehen: GSFLOW calib to 0.2
par_dim_name.append('nssr')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
#par_value.append(0.1*np.ones((dim_value[ind],1)))  # merced
par_value.append(0.2*np.ones((dim_value[ind],1)))  # approx sagehen

par_name.append('ssstor_init')
par_dim_name.append('nssr')
par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
par_value.append(0.*np.ones((dim_value[ind],1)))  


# -- Groundwater flow --
if model_mode == 'PRMS':
    par_name.append('gwflow_coef')
    par_dim_name.append('ngw')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(0.011716573948525877*np.ones((dim_value[ind],1)))  

    par_name.append('gwsink_coef')
    par_dim_name.append('ngw')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(0.*np.ones((dim_value[ind],1)))  

    par_name.append('gwstor_init')
    par_dim_name.append('ngw')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(4.287873048430077*np.ones((dim_value[ind],1)))  

    par_name.append('gwstor_min')
    par_dim_name.append('ngw')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(0.*np.ones((dim_value[ind],1)))  

# -- Streamflow and lake routing --
# below assumes strmflow_module=strmflow_in_out (should also work for
# strmflow, params should just be ignored)
# ****to be read in from GIS info****
ind = np.squeeze(np.where(np.array(dim_name) == 'ncascade'))
if dim_value[ind] == 0:
    par_name.append('hru_segment') # stream segment that hru ultimately flows to (ignored if fl_cascade=1?)
    par_dim_name.append('nhru')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(HRUdata[par_name[-1]])

if model_mode == 'PRMS': # non-cascade
    # ****to be read in from GIS info****
    par_name.append('tosegment') # stream segment flowing into 
    par_dim_name.append('nsegment')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    # par_value.append(segmentdata.data(:, strcmp(segmentdata.colheaders, par_name{ii})))
    par_value.append(segmentdata['OUTSEG'])

    par_name.append('obsin_segment') # index of measured streamflow station that replaces inflow to a segment (can all be 0)
    par_dim_name.append('nsegment')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(0.*np.ones((dim_value[ind],1), int))  

# -- Lake routing --
# ** No parameters here when not using strmflow_module = strmflow_lake **


# -- Output options --
par_name.append('print_freq') # 0=none, 1=run totals, 2=yrly, 4=monthly, 8=daily, or additive combos
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(1)

par_name.append('print_type') # 0=measured and simulated, 1=water bal table, 2=detailed
par_dim_name.append('one')
par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
par_value.append(1)


# -- Subbasin parameters --
# ** No parameters here when using subbasin_flag=0 **
# (can map HRU's to subbasins, for subbasin statistics)

# -- Mapped results parameters --
# ** No parameters here when using mapOutON_OFF=0 **
# (can map groundwater reservoirs to grid cells)

# -- Parameters for cascading-flow simulation --
ind1 = np.squeeze(np.where(np.array(dim_name) == 'ncascade'))
ind2 = np.squeeze(np.where(np.array(dim_name) == 'ncascdgw'))
if dim_value[ind1] > 0 or dim_value[ind2] > 0: 
    par_name.append('cascade_flg') # 0: allow many to many cascade links, 1: force one to one
    par_dim_name.append('one')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(0)  

    par_name.append('cascade_tol') # cascade area below which cascade link is ignored
    par_dim_name.append('one')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(5)  # acres, unsure how to set this
 
    par_name.append('circle_switch') # error check for cascade circles
    par_dim_name.append('one')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(1) # set to 0 for computational savings    

ind = np.squeeze(np.where(np.array(dim_name) == 'ncascade'))
if dim_value[ind] > 0: 
    # ****to be read in from GIS info****
    # NOTE: HRUS are all sub-basins contributing to a segment
    par_name.append('hru_up_id') # index of upslope HRU, -1 for far afield
    par_dim_name.append('ncascade')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(range(1,dim_value[ind]+1))

    # ****to be read in from GIS info****
    par_name.append('hru_down_id') # index of downslope HRU, -1 for far afield, ignored if goes to stream seg
    par_dim_name.append('ncascade')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(range(1,dim_value[ind]+1))
    
    # ****to be read in from GIS info****
    par_name.append('hru_strmseg_down_id') # index of downslope stream segment, 0 if not to stream
    par_dim_name.append('ncascade')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    par_value.append(HRUdata['hru_segment'])

    par_name.append('hru_pct_up') # frac of HRU area contributing to downslope HRU
    par_dim_name.append('ncascade')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[-1]))
    par_value.append(1.*np.ones((dim_value[ind],1)))  

ind = np.squeeze(np.where(np.array(dim_name) == 'ncascdgw'))
if model_mode == 'PRMS' and dim_value[ind] > 0: 
    par_name.append('gw_up_id') # index of upslope GWR, -1 for far afield, ignored if goes to stream seg
    par_dim_name.append('ncascdgw')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(par_name) == 'hru_up_id'))
    par_value.append(par_value[ind])

    par_name.append('gw_down_id') # index of downslope GWR, -1 for far afield, ignored if goes to stream seg
    par_dim_name.append('ncascdgw')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(par_name) == 'hru_down_id'))
    par_value.append(par_value[ind])
    
    par_name.append('gw_strmseg_down_id') # index of downslope stream segment, 0 if not to stream
    par_dim_name.append('ncascdgw')
    par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(par_name) == 'hru_strmseg_down_id'))
    par_value.append(par_value[ind])

    par_name.append('gw_pct_up') # frac of GWR area contributing to downslope HRU
    par_dim_name.append('ncascdgw')
    par_type.append(2) # 1=int, 2=single prec, 3=double prec, 4=char str
    ind = np.squeeze(np.where(np.array(par_name) == 'hru_pct_up'))
    par_value.append(par_value[ind])


# par_name.append('')
# par_dim_name.append('')
# par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
# par_value.append(1)  
# 
# par_name.append('')
# par_dim_name.append('')
# par_type.append(1) # 1=int, 2=single prec, 3=double prec, 4=char str
# par_value.append(1)  

NumPars = len(par_name)



## -----------------------------------------------------------------------
# Many more parameters that are not specified here.  Do a check to make
# sure not using any modules that require unspecified parameters:
ind = np.squeeze(np.where(np.array(dim_name) == 'nobs'))
nobs = dim_value[ind]  
if nobs > 0:
    print('Error! Script not meant for nobs > 0! Exiting... \n')
    quit()
ind = np.squeeze(np.where(np.array(dim_name) == 'nsol'))
nsol = dim_value[ind]  
if nsol > 0:
    print('Error! Script not meant for nsol > 0! Exiting... \n')
    quit()
print('Warning!  Script assumes the following... ')
print('   climate_hru for precip_module and temp_module')
print('   climate_hru for solrad_module')
print('   climate_hru or pet_pt (Priestly-Taylor) for et_module (for PotET)')
print('   transp_tindex for transp_module')
print('   runoff_smidx for runoff_module')
print('   strmflow_in_out for strmflow_module (should also be ok for strmflow)')


#%% -----------------------------------------------------------------------
# Generally, do no change below here

# - Write to parameter file

# -- other remaing variables
par_num_dim = [];
par_num = np.ones((NumPars,1), int)
for ii in range(NumPars):
    par_num_dim.append(np.array(par_dim_name[ii]).size)
    if isinstance(par_dim_name[ii], list):
        for x in par_dim_name[ii]:
#            print x
            ind = np.squeeze(np.where(np.array(dim_name) == x))
#            print ind
            par_num[ii] = par_num[ii] * dim_value[ind]
#            print dim_value[ind]
    else:
        ind = np.squeeze(np.where(np.array(dim_name) == par_dim_name[ii]))
        par_num[ii] = par_num[ii] * dim_value[ind]        

# create PRMS input directory if it does not exist:
if not os.path.isdir(PRMSinput_dir):
    os.mkdir(PRMSinput_dir)
    
# while we're at it, create PRMS output file if it does not exist:
if not os.path.isdir(PRMSoutput_dir):
    os.mkdir(PRMSoutput_dir)

# - Write to Parameter file
line1 = '####'
fobj = open(parfil, 'w+');
fobj.write(title_str1 + '\n')
fobj.write(title_str2 + '\n')

# * Dimensions
dim_line = '** Dimensions **'
fobj.write(dim_line + '\n')
for ii in range(NumDims):
    # Line 1
    fobj.write(line1 + '\n')
    # Line 2
    fobj.write(dim_name[ii] + '\n')
    # Line 3: 
    fobj.write(str(np.squeeze(dim_value[ii])) + '\n' );


par_line = '** Parameters **'
fobj.write(par_line + '\n')

for ii in range(NumPars):
    # Line 1
    fobj.write(line1 + '\n')
    # Line 2
    fobj.write(par_name[ii] + '\n')
    # Line 3: 
    fobj.write(str(np.squeeze(par_num_dim[ii])) + '\n');
    # line 4 to 3+NumOfDims: Name of dimensions, 1 per line 
    if isinstance(par_dim_name[ii], list):
        for x in par_dim_name[ii]: # loop thru list
            fobj.write(x + '\n')
    else: # write the single dimension
        fobj.write(par_dim_name[ii] + '\n')
    # line 3+NumOfDims+1: Number of values 
    fobj.write(str(np.squeeze(par_num[ii])) + '\n')
    # line 3+NumOfDims+2: data type -> 1=int, 2=single prec, 3=double prec, 4=char str
    fobj.write(str(par_type[ii]) + '\n')
    # line 3+NumOfDims+3 to end: parameter values, 1 per line
#    if par_type[ii] == 1:
#        fobj.write(str(np.squeeze(par_value[ii])) + '\n')
#    elif par_type[ii] == 2:
#        fobj.write(str(np.squeeze(par_value[ii])) + '\n');
#    elif par_type[ii] == 4:
#        if par_num[ii] == 1:
#            fobj.write(str(np.squeeze(par_value[ii])) + '\n')
#        else:
#            for x in par_value[ii]:
#                fobj.write(str(np.squeeze(x)) + '\n')
#                
#                
    if par_num[ii] == 1:
        fobj.write(str(np.squeeze(par_value[ii])) + '\n')
#        print str(np.squeeze(par_value[ii])) + '\n'
    else:
        for x in par_value[ii]:
            fobj.write(str(np.squeeze(x)) + '\n')
                        
fobj.close();

