% GSFLOW_print_controlfile1.m
% 11/23/16
% (based on PRMS_print_controlfile5.m, 11/17/15)
% gcng
%
% v2 - Leila and Crystal
%
% Creates inputs files to run GSFLOW, based on the "Sagehen" example but
% modified for Chimborazo's Gavilan Machay watershed. Crystal's AGU2016
% poster.
% 
% GSFLOW Input files:
%   - control file (generate with GSFLOW_print_controlfile1.m - NOT YET WRITTEN)
%   - parameter files (generate with GSFLOW_print_PRMSparamfile1.m and GSFLOW_print_GSFLOWparamfile1.m - NOT YET WRITTEN)
%   - variable / data files (generate with GSFLOW_print_ObsMet_files1.m - NOT YET WRITTEN)
%it 
% (Control file includes names of parameter and data files.  Thus, 
% parameter and variable file names must match those specified there!!)

% **** No Meltwater, 5-yr run-up ****
% - control file name _NoMelt
% - precip WITHOUT melt (precip*.day)
% - no veg shift (ChimTest_GSFLOW.param)
% - no plus 1C (tmin*.day, tmax*.day)
% - run 5yrs (start_date, end_date here and in MODFLOW)
% - use end of 30yr run as initial cond 
%   (test2lay_res.nam, init_var_file=)

clear all, close all, fclose all;
%% Control file

% see PRMS manual: 
%   - list of control variables in Appendix 1 Table 1-2 (p.33-36), 
%   - description of control file on p.126

% general syntax (various parameters listed in succession):
%   line 1: '####'
%   line 2: control parameter name
%   line 3: number of parameter values specified
%   line 4: data type -> 1=int, 2=single prec, 3=double prec, 4=char str
%   line 5-end: parameter values, 1 per line
%
% *** CUSTOMIZE TO YOUR COMPUTER! *****************************************
% NOTE: '/' is directory separator for Linux, '\' for Windows!!

% directory for GSFLOW input and output files (include slash ('/') at end)
control_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/control/';
PRMSinput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/PRMS/';
MODFLOWinput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
PRMSoutput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/PRMS/';

% control_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/ChimTest/';
% PRMSinput_dir = '';
% MODFLOWinput_dir = '';
% PRMSoutput_dir = '';

% GSFLOWoutput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/';
% CSVoutput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/';

% directory with files to be read in to generate PRMS input files (to get
% start/end dates)
in_data_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/DataToReadIn/';
in_climatedata_dir = strcat(in_data_dir, 'climate/'); % specifically climate data

% control file that will be written with this script
% (will be in control_dir with model mode suffix)
con_filname0 = 'ChimTest_NoMelt_Run5yrs';

% command-line executable for GSFLOW (just used to print message)
GSFLOW_exe = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/bin/gsflow';  % recompiled on Ubuntu!
GSFLOW_exe_cmt = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.1Linux/bin/gsflow'; % comment out

% % command-line executable for PRMS (not really used)
% PRMS_exe = '/home/awickert/models/prms4.0.1_linux/bin/prms'; 

% % java program for GUI (optional, only for printing execution info) (not really used)
% PRMS_java_GUI = '/home/awickert/models/prms4.0.1_linux/dist/oui4.jar'; 

% - choose one:
% model_mode = 'WRITE_CLIMATE'; % creates pre-processed climate_hru files
% model_mode = 'PRMS'; % run only PRMS
% model_mode = 'MODFLOW'; % run only MODFLOW-2005
model_mode = 'GSFLOW'; % run coupled PRMS-MODFLOW

% data file that the control file will point to (generate with PRMS_print_climate_hru_files2.m)
datafil = strcat(PRMSinput_dir, 'empty_rep30yr.day');

% parameter file that the control file will point to (generate with PRMS_print_paramfile3.m)
parfil_pre = 'ChimTest'; % will have '_', model_mode following

% MODFLOW namefile that the control file will point to (generate with write_nam_MOD.m)
namfil = strcat(MODFLOWinput_dir, 'test2lay_res.nam');

% output directory that the control file will point to for creating output files (include slash at end!)
outdir = PRMSoutput_dir;

% model start and end dates
% ymdhms_v = [ 2015  6 16 0 0 0; ...
%              2016  6 24 0 0 0];
ymdhms_v = [ 2015  6 16 0 0 0; ...
             2020  6 15 0 0 0];
% ymdhms_v = [ 2015  6 16 0 0 0; ...
%              2045  6 15 0 0 0];

%First MODFLOW initial stress period (can be earlier than model start date;
% useful when using load_init_file and modflow period represents longer 
% period that started earlier).  
ymdhms_m = [2015 6 16 0 0 0];

% initial condition files 
% (see /home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen_restart
% as example for how to stitch together many restarts using a shell script)
if strcmp(model_mode, 'GSFLOW')
    fl_load_init = 1; % 1 to load previously saved initial conditions
    % load initial conditions from this file
%     load_init_file = [PRMSoutput_dir, 'init_cond_infile']; % load initial conditions from this file
    load_init_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/PRMS/init_cond_outfile_30yrspinup_melt';
end
fl_save_init = 1; % 1 to save outputs as initial conditions
save_init_file = [PRMSoutput_dir, 'init_cond_outfile']; % save new results as initial conditions in this file

% Met (temp+precip) data files (ignored if fl_climate_hru = 1)
TempPrecip_datafil = strcat(PRMSinput_dir, 'test_boca_toma_F_in.txt');

% 1: use all pre-processed met data
fl_all_climate_hru = 0; 

if ~fl_all_climate_hru
    % pick one:
%     precip_module = 'precip_1sta';
    precip_module = 'climate_hru';
    
    % pick one:
%     temp_module = 'temp_1sta';
    temp_module = 'climate_hru';

    et_module = 'potet_pt'; % set pt_alpha(nhru, nmonth), currently all = 1.26
    solrad_module = 'ddsolrad';
end

% If climate_hru, use the following file names (else ignored)
% (note that WRITE_CLIMATE will produce files with the below names)
precip_datafil = strcat(PRMSinput_dir, 'precip_rep30yr.day'); % w/o meltwater
% precip_datafil = strcat(PRMSinput_dir, 'precip_rep30yr_melt.day');
% tmax_datafil = strcat(PRMSinput_dir, 'tmax_rep30yr_tadj_plus1C.day'); % to be safe: set as F
% tmin_datafil = strcat(PRMSinput_dir, 'tmin_rep30yr_tadj_plus1C.day'); % to be safe: set as F
tmax_datafil = strcat(PRMSinput_dir, 'tmax_rep30yr_tadj.day'); % to be safe: set as F
tmin_datafil = strcat(PRMSinput_dir, 'tmin_rep30yr_tadj.day'); % to be safe: set as F
solrad_datafil = strcat(PRMSinput_dir, 'swrad_rep30yr.day');
pet_datafil = strcat(PRMSinput_dir, 'potet_rep30yr.day');
% humidity_datafil = strcat(PRMSinput_dir, 'humidity.day'); % for potet_pm
transp_datafil = strcat(PRMSinput_dir, 'transp_rep30yr.day'); % may not be needed in GSFLOW? Is needed!



% *************************************************************************

% Project-specific entries ->

title_str = 'Chimborazo, AGU2016 poster';

% n_par_max should be dynamically generated
n_par_max = 100; % there are a lot, unsure how many total...
con_par_name = cell(n_par_max,1);  % name of control file parameter
con_par_num = zeros(n_par_max,1);  % number of values for a control parameter
con_par_type = zeros(n_par_max,1); % 1=int, 2=single prec, 3=double prec, 4=char str
con_par_values = cell(n_par_max,1);  % control parameter values

ii = 0;

% First 2 blocks should be specified, rest are optional (though default 
% values exist for all variables, see last column of App 1 Table 1-2 p.33).  

% 1 - Variables pertaining to simulation execution and required input and output files 
%     (some variable values left out if default is the only choice we want)

ii = ii+1;
con_par_name{ii} = 'model_mode'; % typically 'PRMS', also 'FROST' or 'WRITE_CLIMATE'
con_par_type(ii) = 4; 
con_par_values{ii} = {model_mode}; %

ii = ii+1;
con_par_name{ii} = 'modflow_name';
con_par_type(ii) = 4; 
con_par_values{ii} = {namfil};

% no more inputs needed for MODFLOW-only runs
if ~strcmp(model_mode, 'MODFLOW')

    % for GSFLOW
    if strcmp(model_mode, 'GSFLOW')
        ii = ii+1;
        con_par_name{ii} = 'csv_output_file';
        con_par_type(ii) = 4; 
        con_par_values{ii} = {[outdir, 'gsflow.csv']};
        
        ii = ii+1;
        con_par_name{ii} = 'modflow_time_zero';
        con_par_type(ii) = 1; 
        con_par_values{ii} = ymdhms_m(1,:); % year, month, day, hour, minute, second

        ii = ii+1;
        con_par_name{ii} = 'gsflow_output_file';
        con_par_type(ii) = 4; 
        con_par_values{ii} = {[outdir, 'gsflow.out']};

        ii = ii+1;
        con_par_name{ii} = 'gsf_rpt'; % flag to create csv output file
        con_par_type(ii) = 1; 
        con_par_values{ii} = 1;

        ii = ii+1;
        con_par_name{ii} = 'rpt_days'; % Frequency with which summary tables are written;default 7
        con_par_type(ii) = 1; 
        con_par_values{ii} = 7;
    end

    ii = ii+1;
    con_par_name{ii} = 'start_time';
    con_par_type(ii) = 1; 
    con_par_values{ii} = ymdhms_v(1,:);

    ii = ii+1;
    con_par_name{ii} = 'end_time';
    con_par_type(ii) = 1; 
    con_par_values{ii} = ymdhms_v(end,:); % year, month, day, hour, minute, second

    ii = ii+1;
    con_par_name{ii} = 'data_file';
    con_par_type(ii) = 4; 
    con_par_values{ii} = {datafil};

    parfil = strcat(PRMSinput_dir, '', parfil_pre, '_', model_mode, '.param');
    ii = ii+1;
    con_par_name{ii} = 'param_file';
    con_par_type(ii) = 4; 
    con_par_values{ii} = {parfil};

    ii = ii+1;
    con_par_name{ii} = 'model_output_file';
    con_par_type(ii) = 4; 
    con_par_values{ii} = {[outdir, 'prms.out']};

    % new for GSFLOW
    ii = ii+1;
    con_par_name{ii} = 'parameter_check_flag';
    con_par_type(ii) = 1; 
    con_par_values{ii} = 1;

    ii = ii+1;
    con_par_name{ii} = 'cbh_check_flag';
    con_par_type(ii) = 1; 
    con_par_values{ii} = 1;

    % 2 - Variables pertaining to module selection and simulation options

    % - module selection:
    %    See PRMS manual: Table 2 (pp. 12-13), summary pp. 14-16, details in
    %    Appendix 1 (pp. 29-122)

    % meteorological data
    ii = ii+1;
    con_par_name{ii} = 'precip_module'; % precip distribution method (should match temp)
    con_par_type(ii) = 4; 
    if fl_all_climate_hru
        con_par_values{ii} = {'climate_hru'}; % climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    else
        con_par_values{ii} = {precip_module}; % climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    end
    if strcmp(con_par_values{ii}, 'climate_hru')
        ii = ii+1;
        con_par_name{ii} = 'precip_day'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {precip_datafil}; % file name
    end

    if strcmp(con_par_values{ii}, 'climate_hru')
        ii = ii+1;
        con_par_name{ii} = 'humidity_day'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {humidity_datafil}; % file name
    end

    ii = ii+1;
    con_par_name{ii} = 'temp_module'; % temperature distribution method (should match precip)
    con_par_type(ii) = 4; 
    if fl_all_climate_hru
        con_par_values{ii} = {'climate_hru'}; % climate_hru, temp_1sta, temp_dist2, temp_laps, ide_dist, xyz_dist
    else
        con_par_values{ii} = {temp_module}; % climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    end
    if strcmp(con_par_values{ii}, 'climate_hru')
        ii = ii+1;
        con_par_name{ii} = 'tmax_day'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {tmax_datafil}; % file name
        ii = ii+1;
        con_par_name{ii} = 'tmin_day'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {tmin_datafil}; % file name
    end

    ii = ii+1;
    con_par_name{ii} = 'solrad_module'; % solar rad distribution method
    con_par_type(ii) = 4; 
    if fl_all_climate_hru
        con_par_values{ii} = {'climate_hru'}; % climate_hru, temp_1sta, temp_dist2, temp_laps, ide_dist, xyz_dist
    else
        con_par_values{ii} = {solrad_module}; % climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    end
    if strcmp(con_par_values{ii}, 'climate_hru')
        ii = ii+1;
        con_par_name{ii} = 'swrad_day'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {solrad_datafil}; % file name
    end

    ii = ii+1;
    con_par_name{ii} = 'et_module'; % method for calculating ET
    con_par_type(ii) = 4; 
    if fl_all_climate_hru
        con_par_values{ii} = {'climate_hru'}; % climate_hru, temp_1sta, temp_dist2, temp_laps, ide_dist, xyz_dist
    else
        con_par_values{ii} = {et_module}; % climate_hru, ide_dist, precip_1sta, precip_dist2, precip_laps, or xyz_dist
    end
    if strcmp(con_par_values{ii}, 'climate_hru')
        ii = ii+1;
        con_par_name{ii} = 'potet_day,'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {pet_datafil}; % file name
    end

    ii = ii+1;
    con_par_name{ii} = 'transp_module'; % transpiration simulation method
    con_par_type(ii) = 4; 
    con_par_values{ii} = {'transp_tindex'}; % climate_hru, transp_frost, or transp_tindex
    if strcmp(con_par_values{ii}, 'climate_hru')
        ii = ii+1;
        con_par_name{ii} = 'transp_day,'; % file with precip data for each HRU
        con_par_type(ii) = 4; 
        con_par_values{ii} = {transp_datafil}; % file name
    end
    
    % Climate-by-HRU Files
    ii = ii+1;
    con_par_name{ii} = 'cbh_binary_flag'; % to use binary format cbh files
    con_par_type(ii) = 1; 
    con_par_values{ii} = 0; % 0 for no, use default

    %Read a CBH file with humidity data
    ii = ii+1;
    con_par_name{ii} = 'humidity_cbh_flag'; % humidity cbh files (for Penman-Monteith ET (potet_pm))
    con_par_type(ii) = 1; 
    con_par_values{ii} = 0; % 0 for no, use default

    %Variable orad
    ii = ii+1;
    con_par_name{ii} = 'orad_flag';  % humidity cbh files (not needed)
    con_par_type(ii) = 1; 
    con_par_values{ii} = 0; % 0 for no, use default    

    ii = ii+1;
    con_par_name{ii} = 'srunoff_module'; % surface runoff/infil calc method
    con_par_type(ii) = 4; 
%     con_par_values{ii} = {'srunoff_smidx_casc'}; % runoff_carea or srunoff_smidx
    con_par_values{ii} = {'srunoff_smidx'}; % runoff_carea or srunoff_smidx (updated name for GSFLOW)

    % strmflow: directly routes runoff to basin outlet 
    % muskingum: moves through stream segments, change in stream segment storages is by Muskingum eq
    % strmflow_in_out: moves through stream segments, input to stream segment = output to stream segment
    % strmflow_lake: for lakes...
    ind = strcmp(con_par_name, 'model_mode');
    if strcmpi(con_par_values{ind}, 'PRMS')
        ii = ii+1;
        con_par_name{ii} = 'strmflow_module'; % streamflow routing method
        con_par_type(ii) = 4; 
        con_par_values{ii} = {'strmflow_in_out'}; % strmflow, muskingum, strmflow_in_out, or strmflow_lake
    end

    % cascade module
    ncascade = 0;
    if ncascade > 0 % default: ncascade = 0
        ii = ii+1;
        con_par_name{ii} = 'cascade_flag'; % runoff routing between HRU's
        con_par_type(ii) = 1; 
        con_par_values{ii} = 1;
    end
    ncascadegw = 0;
    if ncascadegw > 0 % default: ncascadegw = 0
        ii = ii+1;
        con_par_name{ii} = 'cascadegw_flag'; % gw routing between HRU's
        con_par_type(ii) = 1; 
        con_par_values{ii} = 1;
    end

    ii = ii+1;
    con_par_name{ii} = 'dprst_flag'; % flag for depression storage simulations
    con_par_type(ii) = 1; 
    con_par_values{ii} = 0;


    % 3 - Output file: Statistic Variables (statvar) Files
    %     See list in PRMS manual Table 1-5 pp.61-74 for variables you can print
    ii = ii+1;
    con_par_name{ii} = 'statsON_OFF'; % flag to create Statistics output variables
    con_par_type(ii) = 1; 
    con_par_values{ii} = 1;
%     con_par_values{ii} = 0;

    ii = ii+1;
    con_par_name{ii} = 'stat_var_file'; % output Statistics file location, name
    con_par_type(ii) = 4; 
    con_par_values{ii} = {[outdir, 'ChimOut.statvar']};

    ii = ii+1;
    con_par_name{ii} = 'statVar_names';
    con_par_type(ii) = 4; 
    con_par_values{ii} = {...
    'basin_actet', ...
    'basin_cfs', ...
    'basin_gwflow_cfs', ...
    'basin_horad', ...
    'basin_imperv_evap', ...
    'basin_imperv_stor', ...
    'basin_infil', ...
    'basin_intcp_evap', ...
    'basin_intcp_stor', ...
    'basin_perv_et', ...
    'basin_pk_precip', ...
    'basin_potet', ...
    'basin_potsw', ...
    'basin_ppt', ...
    'basin_pweqv', ...
    'basin_rain', ...
    'basin_snow', ...
    'basin_snowcov', ...
    'basin_snowevap', ...
    'basin_snowmelt', ...
    'basin_soil_moist', ...
    'basin_soil_rechr', ...
    'basin_soil_to_gw', ...
    'basin_sroff_cfs', ...
    'basin_ssflow_cfs', ...
    'basin_ssin', ...
    'basin_ssstor', ...
    'basin_tmax', ...
    'basin_tmin', ...
    'basin_slstor', ...
    'basin_pref_stor', ...
    }; 
% below not valid for coupled GSFLOW
%     'basin_gwin', ...
%     'basin_gwsink', ...
%     'basin_gwstor', ...
%     'basin_storage', ...

    % index of statVar_names to be printed to Statistics Output file
    ii = ii+1;
    con_par_name{ii} = 'statVar_element'; % ID numbers for variables in stat_Var_names  
    con_par_num(ii) = length(con_par_values{strcmp(con_par_name, 'statVar_names')});
    con_par_type(ii) = 4; 
    ind = ones(con_par_num(ii), 1); % index of variables (can be 1 to max index of StatVar)
    % add lines here to specify different variable indices other than 1
    ind = num2str(ind);
    con_par_values{ii} = mat2cell(ind, ones(con_par_num(ii), 1), 1);
    
    % 4 - "Animation" output files (spatial data)
    %     See list in: (1) PRMS manual Table 1-5 pp.61-74 and (2) GSFLOW 
    %     Input Instructions manual Table A1-2 for variables you can print
    ii = ii+1;
    con_par_name{ii} = 'aniOutON_OFF'; % flag to create Statistics output variables
    con_par_type(ii) = 1; 
    con_par_values{ii} = 1;
%     con_par_values{ii} = 0;

    ii = ii+1;
    con_par_name{ii} = 'ani_output_file'; % output Statistics file location, name
    con_par_type(ii) = 4; 
    con_par_values{ii} = {[outdir, 'ChimOut.ani']};    
    
    ii = ii+1;
    con_par_name{ii} = 'aniOutVar_names';
    con_par_type(ii) = 4; 
    con_par_values{ii} = {...
    'hru_ppt', ... % [in] Precipitation distributed to each HRU Rain
    'hru_actet', ... % [in] Actual ET for each HRU
    'actet_tot_gwsz', ... % [nhru] Total actual ET from each MODFLOW cell and PRMS soil zone [in]
    'sat_recharge', ... % [nhru] HRU total recharge to the saturated zone 
    'streamflow_sfr', ... % [nsegment] Streamflow as computed by SFR for each segment 
    };     

    % 4 - For GUI (otherwise ignored during command line execution)

    % GSFLOW: ignore these
%     ii = ii+1;
%     con_par_name{ii} = 'ndispGraphs'; % number runtime graphs with GUI
%     con_par_type(ii) = 1; 
%     con_par_values{ii} = 2;
% 
%     ii = ii+1;
%     con_par_name{ii} = 'dispVar_names'; % variables for runtime plot
%     con_par_type(ii) = 4; 
%     con_par_values{ii} = { ...
%     'basin_cfs', ...
%     'runoff'};
% 
%     % index of dispVar_names to be displayed in runtime plots
%     ii = ii+1;
%     con_par_name{ii} = 'dispVar_element'; % variable indices for runtime plot
%     con_par_num(ii) = length(con_par_values{strcmp(con_par_name, 'dispVar_names')});
%     con_par_type(ii) = 4; 
%     ind = ones(con_par_num(ii), 1); % index of variables
%     % add lines here to specify different variable indices other than 1
%     ind = num2str(ind);
%     con_par_values{ii} = mat2cell(ind, ones(con_par_num(ii), 1), 1);
% 
%     % which plot (of ndispGraphs) to show variable (dispVar_names) on
%     ii = ii+1;
%     con_par_name{ii} = 'dispVar_plot';
%     con_par_num(ii) = length(con_par_values{strcmp(con_par_name, 'dispVar_names')});
%     con_par_type(ii) = 4; 
%     ind = ones(con_par_num(ii), 1); % index of plots
%     % add lines here to specify different plot indices other than 1
%     ind(2) = 2;
%     ind = num2str(ind);
%     con_par_values{ii} = mat2cell(ind, ones(con_par_num(ii), 1), 1);
% 
%     ii = ii+1;
%     con_par_name{ii} = 'dispGraphsBuffSize'; % num timesteps (days) before updating runtime plot
%     con_par_type(ii) = 1; 
%     con_par_values{ii} = 1;
% 
%     ii = ii+1;
%     con_par_name{ii} = 'initial_deltat'; % initial time step length (hrs)
%     con_par_type(ii) = 2; 
%     con_par_values{ii} = 24.0; % 24 hrs matches model's daily time step
% 
%     % previously for PRMS, omit
%     ii = ii+1;
%     con_par_name{ii} = 'executable_desc';
%     con_par_type(ii) = 4; 
%     con_par_values{ii} = {'PRMS IV'};
% 
%     ii = ii+1;
%     con_par_name{ii} = 'executable_model';
%     con_par_type(ii) = 4; 
%     con_par_values{ii} = {PRMS_exe};


    % 5 - Initial condition file

    % (default is init_vars_from_file = 0, but still need to specify for GUI)
    ii = ii+1;
    con_par_name{ii} = 'init_vars_from_file'; % use IC from initial cond file
    con_par_type(ii) = 1; 
    con_par_values{ii} = fl_load_init; % 0 for no, use default

    if fl_load_init
        ii = ii+1;
        con_par_name{ii} = 'var_init_file'; % use IC from initial cond file
        con_par_type(ii) = 4; 
        con_par_values{ii} = {load_init_file}; % 0 for no, use default
    end
    
    % (default is save_vars_to_file = 0, but still need to specify for GUI)
    ii = ii+1;
    con_par_name{ii} = 'save_vars_to_file'; % save IC to output file
    con_par_type(ii) = 1; 
    con_par_values{ii} = fl_save_init;
    
    if fl_save_init
        ii = ii+1;
        con_par_name{ii} = 'var_save_file'; % use IC from initial cond file
        con_par_type(ii) = 4; 
        con_par_values{ii} = {save_init_file}; % 0 for no, use default
    end


    % 6 - Suppress printing of some execution warnings

    ii = ii+1;
    con_par_name{ii} = 'print_debug';
    con_par_type(ii) = 1; 
    con_par_values{ii} = -1;
end

%% -----------------------------------------------------------------------
% Generally, do not change below here

con_filname = strcat(control_dir, '', con_filname0, '_', model_mode, '.control');

% - Write to control file

if ~strcmp(model_mode, 'MODFLOW')
    if ~isempty(find(strcmp(con_par_name, 'statVar_names'), 1))
        ii = ii+1;
        con_par_name{ii} = 'nstatVars'; % num output vars in statVar_names (for Statistics output file)
        con_par_type(ii) = 1; 
        con_par_values{ii} = length(con_par_values{strcmp(con_par_name, 'statVar_names')});
    end

    if ~isempty(find(strcmp(con_par_name, 'aniOutVar_names'), 1))
        ii = ii+1;
        con_par_name{ii} = 'naniOutVars'; % num output vars in aniOutVar_names (for animation output file)
        con_par_type(ii) = 1; 
        con_par_values{ii} = length(con_par_values{strcmp(con_par_name, 'aniOutVar_names')});
    end
end

nvars = ii;


for ii = 1: nvars
    con_par_num(ii) = length(con_par_values{ii});
end


% - Write to control file
line1 = '####';
fmt_types = {'%d', '%f', '%f', '%s'};
fid = fopen(con_filname, 'w');
fprintf(fid,'%s\n', title_str);
for ii = 1: nvars
    % Line 1
    fprintf(fid,'%s\n', line1);
    % Line 2
    fprintf(fid,'%s\n', con_par_name{ii});
    % Line 3: 
    fprintf(fid,'%d\n', con_par_num(ii));
    % Line 4: 
    fprintf(fid,'%d\n', con_par_type(ii));
    % Line 5 to end:
    if con_par_type(ii) == 4 
        for jj = 1: con_par_num(ii) 
            fprintf(fid,[fmt_types{con_par_type(ii)}, '\n'], con_par_values{ii}{jj});
        end
    else
        fprintf(fid,[fmt_types{con_par_type(ii)}, '\n'], con_par_values{ii});
    end
end

fclose(fid);

%% ------------------------------------------------------------------------
% Prepare for model execution

if ~strcmp(model_mode, 'MODFLOW')

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    fprintf('Make sure the below data files are ready: \n  %s\n', datafil);
    fprintf('  %s\n', precip_datafil);
    fprintf('  %s\n', tmax_datafil);
    fprintf('  %s\n', tmin_datafil);
    fprintf('  %s\n', solrad_datafil);

    fprintf('Make sure the below parameter file is ready: \n  %s\n', parfil);
end

cmd_str_cmt = ['#', GSFLOW_exe_cmt, ' ', con_filname, ' &> out.txt'];
cmd_str = [GSFLOW_exe, ' ', con_filname, ' &> out.txt'];
fprintf('To run command-line execution, enter at prompt: \n %s\n', cmd_str);

runscriptfil = strcat(control_dir, '', con_filname0, '_', model_mode, '.sh');
fid = fopen(runscriptfil, 'wt');
fprintf(fid, '%s', cmd_str_cmt);
fprintf(fid, '\n\n');
fprintf(fid, '%s', cmd_str);
fclose(fid);
system(['chmod 777 ', runscriptfil]);

%     gui_str1 = [cmd_str, ' -print'];
%     gui_str2 = ['java -cp ', PRMS_java_GUI, ' oui.mms.gui.Mms ', con_filname];
%     fprintf('To launch GUI, enter the following 2 lines (one at a time) at prompt: \n %s\n %s\n', gui_str1, gui_str2);
