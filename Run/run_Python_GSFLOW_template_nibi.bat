
:: full pathname for location of GSFLOWGRASS_toolkit
SET toolkit_dir=C:\Users\gcng\GSFLOW-GRASS

:: full pathname for Settings File
SET Settings_file=C:\Users\gcng\GSFLOW-GRASS\Run\settings_template_spinup_nibi_1p2p0.ini

SET python2.7=C:\ProgramData\Anaconda2\python.exe

:: ----------------------------------------------------------------------------
:: **For default implementation, user does not need to make changes to lines below**
:: Go to directory with input file builder scripts
cd %toolkit_dir%\input_file_builder  

:: If fl_create_hydcond=1 in Settings File, this creates a file with spatially distributed hydraulic conductivity.
%python2.7% Create_hydcond_array.py %Settings_file%

:: If fl_print_climate_hru=1 in Settings File, this creates climate_hru files with uniform conditions 
%python2.7% GSFLOW_print_data_climatehru_files1_metric.py %Settings_file%

:: This creates GSFLOW control file.  See Python script to make changes beyond parameters set in Settings File.
%python2.7% GSFLOW_print_controlfile_current.py %Settings_file%

:: This creates PRMS parameter file.  See Python script to make changes beyond parameters set in Settings File.
%python2.7% GSFLOW_print_PRMSparamfile_current.py %Settings_file%

:: This creates the various MODFLOW input files.  See Python script to make changes beyond parameters set in Settings File.
%python2.7% print_MODFLOW_inputs_res_NWT_current.py %Settings_file%

:: Go to directory with run script
cd %toolkit_dir%\Run  

:: This executes GSFLOW model with the input files created above
%python2.7% run_GSFLOW.py %Settings_file%

pause
