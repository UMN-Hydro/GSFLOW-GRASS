#!/bin/bash

# full pathname for location of GSFLOWGRASS_toolkit
preproc_dir="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/"

# full pathname for Settings File
Settings_file="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/settings_AW_test.ini"

# ------------------------------------------------------------------------------
# For default implementation, user does not need to make changes below

# If fl_create_hydcond=1 in Settings File, this creates a file with spatially distributed hydraulic conductivity.
python2.7 ${preproc_dir}Create_hydcond_array.py ${Settings_file}

# If fl_print_climate_hru=1 in Settings File, this creates climate_hru files with uniform conditions 
python2.7 ${preproc_dir}GSFLOW_print_data_climatehru_files1_metric.py ${Settings_file}

# This creates GSFLOW control file.  See Python script to make changes beyond parameters set in Settings File.
python2.7 ${preproc_dir}GSFLOW_print_controlfile_current.py ${Settings_file}

# This creates PRMS parameter file.  See Python script to make changes beyond parameters set in Settings File.
python2.7 ${preproc_dir}GSFLOW_print_PRMSparamfile_current.py ${Settings_file}

# This creates the various MODFLOW input files.  See Python script to make changes beyond parameters set in Settings File.
python2.7 ${preproc_dir}print_MODFLOW_inputs_res_NWT_current.py ${Settings_file}

# This executes GSFLOW model with the input files created above
python2.7 ${preproc_dir}run_GSFLOW.py ${Settings_file}


