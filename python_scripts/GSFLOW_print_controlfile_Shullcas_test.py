"""
Created on Sun Sep 10 22:06:46 2017
Converting from GSFLOW_print_controlfile4_gcng_melt30yr.m

 Creates inputs files to run GSFLOW, based on the "Sagehen" example but
 modified for Chimborazo's Gavilan Machay watershed. Crystal's AGU2016
 poster.
 
 GSFLOW Input files:
   - control file (generate with GSFLOW_print_controlfile1.m - NOT YET WRITTEN)
   - parameter files (generate with GSFLOW_print_PRMSparamfile1.m and GSFLOW_print_GSFLOWparamfile1.m - NOT YET WRITTEN)
   - variable / data files (generate with GSFLOW_print_ObsMet_files1.m - NOT YET WRITTEN)

 (Control file includes names of parameter and data files.  Thus, 
 parameter and variable file names must match those specified there!!)

 **** Meltwater, 30-yrs spin-up ****
 - precip with melt (precip*.day)
 - no veg shift (ChimTest_GSFLOW.param)
 - no plus 1C (tmin*.day, tmax*.day)
 - run 30yrs (start_date, end_date here and in MODFLOW)

@author: gcng
"""

#==============================================================================
# ## Control file
# 
# # see PRMS manual: 
# #   - list of control variables in Appendix 1 Table 1-2 (p.33-36), 
# #   - description of control file on p.126
# 
# # general syntax (various parameters listed in succession):
# #   line 1: '####'
# #   line 2: control parameter name
# #   line 3: number of parameter values specified
# #   line 4: data type -> 1=int, 2=single prec, 3=double prec, 4=char str
# #   line 5-end: parameter values, 1 per line
# #
# # *** CUSTOMIZE TO YOUR COMPUTER! *****************************************
# # NOTE: '/' is directory separator for Linux, '\' for Windows!!
#==============================================================================

import datetime as dt
import numpy as np # matlab core
import os  # os functions
import settings_test
import platform


# creates start and end date strings
def datetime_to_list(datetime):
    dt_timetuple = datetime.utctimetuple()
    return [dt_timetuple.tm_year, dt_timetuple.tm_mon,
            dt_timetuple.tm_mday, dt_timetuple.tm_hour,
            dt_timetuple.tm_min, dt_timetuple.tm_sec]


if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

# control file that will be written with this script
# (will be in control_dir with model mode suffix)

con_filname0 = settings_test.PROJ_CODE


# - choose one:
# model_mode = 'WRITE_CLIMATE'; # creates pre-processed climate_hru files
# model_mode = 'PRMS'; # run only PRMS
# model_mode = 'MODFLOW'; # run only MODFLOW-2005
model_mode = 'GSFLOW' # run coupled PRMS-MODFLOW

# parameter file that the control file will point to (generate with PRMS_print_paramfile3.m)
parfil_pre = settings_test.PROJ_CODE # will have '_', model_mode following

# MODFLOW namefile that the control file will point to (generate with write_nam_MOD.m)
namfil = '..' + slashstr + settings_test.MODFLOWinput_dir_rel + slashstr + settings_test.PROJ_CODE + '.nam'

# input directory that the control file will point to for reading in PRMS_GSFLOW inputs
indir_rel = '..' + slashstr + settings_test.PRMSinput_dir_rel + slashstr 


# output directory that the control file will point to for creating PRMS_GSFLOW output files (include slash at end!)
outdir_rel = '..' + slashstr + settings_test.PRMSoutput_dir_rel + slashstr 

# model start and end dates
start_date = dt.datetime.strptime(settings_test.START_DATE, "%Y-%m-%d")
end_date = dt.datetime.strptime(settings_test.END_DATE, "%Y-%m-%d")

ymdhms_v = np.array([datetime_to_list(start_date),
                     datetime_to_list(end_date)])                     

#First MODFLOW initial stress period (can be earlier than model start date;
# useful when using load_init_file and modflow period represents longer 
# period that started earlier).  
ymdhms_m = ymdhms_v[0]

# initial condition files 
# (see /home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen_restart
# as example for how to stitch together many restarts using a shell script)
if model_mode == 'GSFLOW':
    fl_load_init = 0 # 1 to load previously saved initial conditions
    # load initial conditions from this file
#     load_init_file = PRMSoutput_dir + 'init_cond_infile' # load initial conditions from this file
    load_init_file = outdir_rel + 'init_cond_outfile' # eventually change so user-specified?

fl_save_init = 1 # 1 to save outputs as initial conditions
save_init_file = outdir_rel + 'init_cond_outfile' # save new results as initial conditions in this file

# 1: use all pre-processed met data
fl_all_climate_hru = 0 # (could set to = False)

if fl_all_climate_hru == 0:
    # pick one:
#     precip_module = 'precip_1sta'
    precip_module = 'climate_hru'
    
    # pick one:
#     temp_module = 'temp_1sta';
    temp_module = 'climate_hru';

    et_module = 'potet_pt' # set pt_alpha(nhru, nmonth), currently all = 1.26
    solrad_module = 'ddsolrad'

# If climate_hru, use the following file names (else ignored)
# (note that WRITE_CLIMATE will produce files with the below names)
# precip_datafil = strcat(PRMSinput_dir, 'precip_rep30yr.day'); # w/o meltwater
precip_datafil = indir_rel + 'precip.day'
# tmax_datafil = strcat(PRMSinput_dir, 'tmax_rep30yr_tadj_plus1C.day'); # to be safe: set as F
# tmin_datafil = strcat(PRMSinput_dir, 'tmin_rep30yr_tadj_plus1C.day'); # to be safe: set as F
tmax_datafil = indir_rel + 'tmax.day' # to be safe: set as F
tmin_datafil = indir_rel + 'tmin.day' # to be safe: set as F
solrad_datafil = indir_rel + 'swrad.day'
pet_datafil = indir_rel + 'potet.day'
humidity_datafil = indir_rel + 'humidity.day' # for potet_pm
transp_datafil = indir_rel + 'transp.day' # may not be needed in GSFLOW? Is needed!

# data file that the control file will point to (generate with PRMS_print_climate_hru_files2.m)
datafil = indir_rel + 'empty.day'




# *************************************************************************

# Project-specific entries ->

title_str = settings_test.PROJ_NAME

# n_par_max should be dynamically generated
con_par_name = []  # name of control file parameter
con_par_num = []  # number of values for a control parameter
con_par_type = [] # 1=int, 2=single prec, 3=double prec, 4=char str
#con_par_values = np.empty((n_par_max,1), dtype=np.object) # control parameter values
con_par_values = [] # control parameter values


# First 2 blocks should be specified, rest are optional (though default 
# values exist for all variables, see last column of App 1 Table 1-2 p.33).  

# 1 - Variables pertaining to simulation execution and required input and output files 
#     (some variable values left out if default is the only choice we want)

con_par_name.append('model_mode') # typically 'PRMS', also 'FROST' or 'WRITE_CLIMATE'
con_par_type.append(4) 
con_par_values.append(model_mode) #

con_par_name.append('modflow_name')
con_par_type.append(4) 
con_par_values.append(namfil)

# no more inputs needed for MODFLOW-only runs
if model_mode != 'MODFLOW':

    # for GSFLOW
    if model_mode == 'GSFLOW':
        con_par_name.append('csv_output_file')
        con_par_type.append(4)
        con_par_values.append(outdir_rel + 'gsflow.csv')
        
        con_par_name.append('modflow_time_zero')
        con_par_type.append(1)
        con_par_values.append(ymdhms_m) # year, month, day, hour, minute, second

        con_par_name.append('gsflow_output_file')
        con_par_type.append(4)
        con_par_values.append(outdir_rel + 'gsflow.out')

        con_par_name.append('gsf_rpt') # flag to create csv output file
        con_par_type.append(1)
        con_par_values.append(1)

        con_par_name.append('rpt_days') # Frequency with which summary tables are written;default 7
        con_par_type.append(1)
        con_par_values.append(7)
    

    con_par_name.append('start_time')
    con_par_type.append(1)
    con_par_values.append(ymdhms_v[0,:])

    con_par_name.append('end_time')
    con_par_type.append(1)
    con_par_values.append(ymdhms_v[1,:]) # year, month, day, hour, minute, second

    con_par_name.append('data_file')
    con_par_type.append(4)
    con_par_values.append(datafil)

    parfil = indir_rel + parfil_pre + '_' + model_mode + '.param'
    con_par_name.append('param_file')
    con_par_type.append(4)
    con_par_values.append(parfil)

    con_par_name.append('model_output_file')
    con_par_type.append(4 )
    con_par_values.append(outdir_rel + 'prms.out')

    # new for GSFLOW
    con_par_name.append('parameter_check_flag')
    con_par_type.append(1)
    con_par_values.append(1)

    con_par_name.append('cbh_check_flag')
    con_par_type.append(1)
    con_par_values.append(1)

    # 2 - Variables pertaining to module selection and simulation options

    # - module selection:
    #    See PRMS manual: Table 2 (pp. 12-13), summary pp. 14-16, details in
    #    Appendix 1 (pp. 29-122)

    # meteorological data
    con_par_name.append('precip_module') # precip distribution method (should match temp)
    con_par_type.append(4) 
    if fl_all_climate_hru == 1:
        con_par_values.append('climate_hru') # climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    else:
        con_par_values.append(precip_module) # climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    
    if con_par_values[-1] == 'climate_hru': # index -1 for last element
        con_par_name.append('precip_day') # file with precip data for each HRU
        con_par_type.append(4) 
        con_par_values.append(precip_datafil) # file name
    

        # Below: harmless mistake in original matlab version

    con_par_name.append('humidity_day')
    con_par_type.append(4)
    con_par_values.append(humidity_datafil) # file name

    con_par_name.append('temp_module') # temperature distribution method (should match precip)
    con_par_type.append(4) 
    if fl_all_climate_hru == 1:
        con_par_values.append('climate_hru') # climate_hru, temp_1sta, temp_dist2, temp_laps, ide_dist, xyz_dist
    else:
        con_par_values.append(temp_module) # climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    
    if con_par_values[-1] == 'climate_hru':
        con_par_name.append('tmax_day') # file with precip data for each HRU
        con_par_type.append(4) 
        con_par_values.append(tmax_datafil) # file name
        con_par_name.append('tmin_day') # file with precip data for each HRU
        con_par_type.append(4) 
        con_par_values.append(tmin_datafil) # file name
    

    con_par_name.append('solrad_module') # solar rad distribution method
    con_par_type.append(4) 
    if fl_all_climate_hru == 1:
        con_par_values.append('climate_hru') # climate_hru, temp_1sta, temp_dist2, temp_laps, ide_dist, xyz_dist
    else:
        con_par_values.append(solrad_module) # climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    
    if con_par_values[-1] == 'climate_hru':
        con_par_name.append('swrad_day') # file with precip data for each HRU
        con_par_type.append(4) 
        con_par_values.append(solrad_datafil) # file name

    con_par_name.append('et_module') # method for calculating ET
    con_par_type.append(4) 
    if fl_all_climate_hru == 1:
        con_par_values.append('climate_hru') # climate_hru, temp_1sta, temp_dist2, temp_laps, ide_dist, xyz_dist
    else:
        con_par_values.append(et_module) # climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    if con_par_values[-1] == 'climate_hru':
        con_par_name.append('potet_day,') # file with precip data for each HRU
        con_par_type.append(4) 
        con_par_values.append(pet_datafil) # file name

    con_par_name.append('transp_module') # transpiration simulation method
    con_par_type.append(4) 
    con_par_values.append('transp_tindex') # climate_hru, transp_frost, or transp_tindex
    if con_par_values[-1] == 'climate_hru':
        con_par_name.append('transp_day,') # file with precip data for each HRU
        con_par_type.append(4) 
        con_par_values.append(transp_datafil) # file name
    
    # Climate-by-HRU Files
    con_par_name.append('cbh_binary_flag') # to use binary format cbh files
    con_par_type.append(1) 
    con_par_values.append(0) # 0 for no, use default

    #Read a CBH file with humidity data
    con_par_name.append('humidity_cbh_flag') # humidity cbh files (for Penman-Monteith ET (potet_pm))
    con_par_type.append(1) 
    con_par_values.append(0) # 0 for no, use default

    #Variable orad
    con_par_name.append('orad_flag')  # humidity cbh files (not needed)
    con_par_type.append(1) 
    con_par_values.append(0) # 0 for no, use default    

    con_par_name.append('srunoff_module') # surface runoff/infil calc method
    con_par_type.append(4) 
#     con_par_values.append('srunoff_smidx_casc') # runoff_carea or srunoff_smidx
    con_par_values.append('srunoff_smidx') # runoff_carea or srunoff_smidx (updated name for GSFLOW)

    # strmflow: directly routes runoff to basin outlet 
    # muskingum: moves through stream segments, change in stream segment storages is by Muskingum eq
    # strmflow_in_out: moves through stream segments, input to stream segment = output to stream segment
    # strmflow_lake: for lakes...
    ind = np.squeeze(np.where(np.array(con_par_name) == 'model_mode'))
    if con_par_values[ind] == 'PRMS':
        con_par_name.append('strmflow_module') # streamflow routing method
        con_par_type.append(4) 
        con_par_values.append('strmflow_in_out') # strmflow, muskingum, strmflow_in_out, or strmflow_lake
    
    # cascade module
    ncascade = 0
    if ncascade > 0: # default: ncascade = 0
        con_par_name.append('cascade_flag') # runoff routing between HRU's
        con_par_type.append(1) 
        con_par_values.append(1)    
    ncascadegw = 0
    if ncascadegw > 0: # default: ncascadegw = 0
        con_par_name.append('cascadegw_flag') # gw routing between HRU's
        con_par_type.append(1) 
        con_par_values.append(1)

    con_par_name.append('dprst_flag') # flag for depression storage simulations
    con_par_type.append(1) 
    con_par_values.append(0)


    # 3 - Output file: Statistic Variables (statvar) Files
    #     See list in PRMS manual Table 1-5 pp.61-74 for variables you can print
    con_par_name.append('statsON_OFF') # flag to create Statistics output variables
    con_par_type.append(1) 
    con_par_values.append(1)
#     con_par_values.append(0)

    con_par_name.append('stat_var_file') # output Statistics file location, name
    con_par_type.append(4) 
    con_par_values.append(outdir_rel + '{}.statvar'.format(settings_test.PROJ_CODE))

    con_par_name.append('statVar_names')
    con_par_type.append(4)
    con_par_values.append(
    np.array(['basin_actet', 
    'basin_cfs', 
    'basin_gwflow_cfs',
    'basin_horad',
    'basin_imperv_evap',
    'basin_imperv_stor',
    'basin_infil',
    'basin_intcp_evap',
    'basin_intcp_stor',
    'basin_perv_et',
    'basin_pk_precip',
    'basin_potet',
    'basin_potsw',
    'basin_ppt',
    'basin_pweqv',
    'basin_rain',
    'basin_snow',
    'basin_snowcov',
    'basin_snowevap',
    'basin_snowmelt',
    'basin_soil_moist',
    'basin_soil_rechr',
    'basin_soil_to_gw',
    'basin_sroff_cfs',
    'basin_ssflow_cfs',
    'basin_ssin',
    'basin_ssstor',
    'basin_tmax',
    'basin_tmin',
    'basin_slstor',
    'basin_pref_stor']))

    # index of statVar_names to be printed to Statistics Output file
    con_par_name.append('statVar_element') # ID numbers for variables in stat_Var_names  
    ind = np.squeeze(np.where(np.array(con_par_name) == 'statVar_names'))
    con_par_num_i = con_par_values[ind].size
    con_par_type.append(4) 
    ind = np.ones((con_par_num_i, 1), int).astype(str) # index of variables (can be 1 to max index of StatVar)
    # add lines here to specify different variable indices other than 1
    con_par_values.append(ind)
    
    # 4 - "Animation" output files (spatial data)
    #     See list in: (1) PRMS manual Table 1-5 pp.61-74 and (2) GSFLOW 
    #     Input Instructions manual Table A1-2 for variables you can print
    con_par_name.append('aniOutON_OFF') # flag to create Statistics output variables
    con_par_type.append(1) 
    con_par_values.append(1)
#     con_par_values.append(0)

    con_par_name.append('ani_output_file') # output Statistics file location, name
    con_par_type.append(4) 
    con_par_values.append(outdir_rel + '{}.ani'.format(settings_test.PROJ_CODE))    
    
    con_par_name.append('aniOutVar_names')
    con_par_type.append(4) 
    con_par_values.append(np.array(['hru_ppt',  # [in] Precipitation distributed to each HRU Rain
    'hru_actet',  # [in] Actual ET for each HRU
    'actet_tot_gwsz',  # [nhru] Total actual ET from each MODFLOW cell and PRMS soil zone [in]
    'sat_recharge',  # [nhru] HRU total recharge to the saturated zone 
    'streamflow_sfr']))  # [nsegment] Streamflow as computed by SFR for each segment 
        
    # 4 - For GUI (otherwise ignored during command line execution)

    # GSFLOW: ignore these
#     con_par_name.append('ndispGraphs') # number runtime graphs with GUI
#     con_par_type.append(1) 
#     con_par_values.append(2)
# 
#     con_par_name.append('dispVar_names') # variables for runtime plot
#     con_par_type.append(4) 
#     con_par_values.append( 
#     np.array(['basin_cfs', 
#     'runoff']))
# 
#     # index of dispVar_names to be displayed in runtime plots
#     con_par_name.append('dispVar_element') # variable indices for runtime plot
#     con_par_type.append(4) 
#     ind = np.squeeze(np.where(np.array(con_par_name) == 'dispVar_names'))
#     con_par_num_i = con_par_values[ind].size
#     con_par_type.append(4) 
#     ind = np.ones((con_par_num_i, 1), int).astype(str) # index of variables (can be 1 to max index of StatVar)
#     # add lines here to specify different variable indices other than 1
#     ind[1] = '2'
#     con_par_values.append(ind)
# 
#     con_par_name.append('dispGraphsBuffSize') # num timesteps (days) before updating runtime plot
#     con_par_type.append(1) 
#     con_par_values.append(1)
# 
#     con_par_name.append('initial_deltat') # initial time step length (hrs)
#     con_par_type.append(2) 
#     con_par_values.append(24.0) # 24 hrs matches model's daily time step
# 
#     # previously for PRMS, omit
#     con_par_name.append('executable_desc')
#     con_par_type.append(4) 
#     con_par_values.append('PRMS IV')
# 
#     con_par_name.append('executable_model')
#     con_par_type.append(4) 
#     con_par_values.append(PRMS_exe)

    # 5 - Initial condition file

    # (default is init_vars_from_file = 0, but still need to specify for GUI)
    con_par_name.append('init_vars_from_file') # use IC from initial cond file
    con_par_type.append(1) 
    con_par_values.append(fl_load_init) # 0 for no, use default

    if fl_load_init == 1:
        con_par_name.append('var_init_file') # use IC from initial cond file
        con_par_type.append(4) 
        con_par_values.append(load_init_file) # 0 for no, use default
    
    # (default is save_vars_to_file = 0, but still need to specify for GUI)
    con_par_name.append('save_vars_to_file') # save IC to output file
    con_par_type.append(1) 
    con_par_values.append(fl_save_init)
    
    if fl_save_init == 1:
        con_par_name.append('var_save_file') # use IC from initial cond file
        con_par_type.append(4) 
        con_par_values.append(save_init_file) # 0 for no, use default


    # 6 - Suppress printing of some execution warnings
    con_par_name.append('print_debug')
    con_par_type.append(1) 
    con_par_values.append(-1)

# % % -----------------------------------------------------------------------
# Generally, do not change below here

con_filname = con_filname0 + '_' + model_mode + '.control'
con_filname_fullpath = settings_test.control_dir + slashstr + con_filname

# - Write to control file

if model_mode != 'MODFLOW':
    ind = np.squeeze(np.where(np.array(con_par_name) == 'statVar_names'))
    if ind.size > 0:
        con_par_name.append('nstatVars') # num output vars in statVar_names (for Statistics output file)
        con_par_type.append(1) 
        con_par_values.append(con_par_values[ind].size)
    
    ind = np.squeeze(np.where(np.array(con_par_name) == 'aniOutVar_names'))
    if ind.size > 0:
        con_par_name.append('naniOutVars') # num output vars in aniOutVar_names (for animation output file)
        con_par_type.append(1) 
        con_par_values.append(con_par_values[ind].size)

nvars = len(con_par_name)

con_par_num = np.ones((nvars,1), int)
ii = 0
for x in con_par_values:
    if isinstance(x, (list, np.ndarray)): # check if list or array
        if isinstance(x, np.ndarray):
            con_par_num[ii] = x.size
        else:
            con_par_num[ii] = len(x)
    ii = ii + 1


# - Write to control file
fobj = open(con_filname_fullpath, 'w+') # w+ for write and read
fobj.write(title_str + '\n')
line1 = '####\n'
for ii in range(0, nvars):
    # Line 1
    fobj.write(line1)
    # Line 2
    fobj.write(con_par_name[ii] + '\n');
    # Line 3: 
    fobj.write(str(np.squeeze(con_par_num[ii])) + '\n');
    # Line 4: 
    fobj.write(str(np.squeeze(con_par_type[ii])) + '\n');
    # Line 5 to end:
    if con_par_num[ii] == 1:
        fobj.write(str(np.squeeze(con_par_values[ii])) + '\n')
    else:
        for x in con_par_values[ii]:
            fobj.write(str(np.squeeze(x)) + '\n');
fobj.close()

# % % ------------------------------------------------------------------------
 # Prepare for model execution

#if model_mode != 'MODFLOW':
#
#    print 'Make sure the below data files are ready: \n   {}\n'.format(datafil)
#    print '   {}\n'.format(precip_datafil)
#    print '   {}\n'.format(tmax_datafil)
#    print '   {}\n'.format(tmin_datafil)
#    print '   {}\n'.format(solrad_datafil)
#
#    print 'Make sure the below parameter file is ready: \n   {}\n'.format(parfil)


if platform.system() == 'Linux':
    cmd_str = settings_test.GSFLOW_exe + ' ' + con_filname + ' &> out.txt'
else:
    cmd_str = settings_test.GSFLOW_exe + ' ' + con_filname


#cmd_str_cmt = '#' + GSFLOW_exe_cmt + ' ' + con_filname + ' &> out.txt'
print '*** To run command-line execution --> '
print '***   Go to ' +  settings_test.control_dir
print '***   and enter at prompt: \n  {}\n'.format(cmd_str)

runscriptfil = settings_test.control_dir + slashstr + con_filname0 + '_' + model_mode + '.sh'
fobj = open(runscriptfil, 'w+') 
#fobj.write(cmd_str_cmt);
fobj.write('pwd0=`pwd` \n')
fobj.write('cd ' + settings_test.control_dir + '\n\n')
fobj.write('echo Running GSFLOW in ' + settings_test.gsflow_simdir + ': \n')
fobj.write("echo"  + " '" + cmd_str + "' \n\n");
fobj.write(cmd_str + '\n\n');
fobj.write('cd $pwd0 \n')
fobj.close()
os.chmod(runscriptfil, 0777)


