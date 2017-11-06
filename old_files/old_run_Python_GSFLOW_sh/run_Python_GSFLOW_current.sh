#!/bin/bash

# full pathname for location of GSFLOWGRASS_toolkit
preproc_dir="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/"

# full pathname for Settings File
Settings_file="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/settings_2lay_res.ini"

# ------------------------------------------------------------------------------
# (Do not change below this line)
cp -p $Settings_file settings.ini

# if omit input_file, will use default file: settings.ini
python2.7 ${preproc_dir}settings_test.py 

# If fl_create_hydcond=1 in Settings File, this creates a file with spatially distributed hydraulic conductivity.
python2.7 ${preproc_dir}Create_hydcond_array.py

# If fl_print_climate_hru=1 in Settings File, this creates climate_hru files with uniform conditions 
python2.7 ${preproc_dir}GSFLOW_print_data_climatehru_files1.py

# This creates GSFLOW control file.  See Python script to make changes beyond parameters set in Settings File.
python2.7 ${preproc_dir}GSFLOW_print_controlfile_current.py

# This creates PRMS parameter  file.  See Python script to make changes beyond parameters set in Settings File.
python2.7 ${preproc_dir}GSFLOW_print_PRMSparamfile_current.py

# This creates the various MODFLOW input files.  See Python script to make changes beyond parameters set in Settings File.
python2.7 ${preproc_dir}print_MODFLOW_inputs_res_NWT_current.py


python2.7 ${preproc_dir}run_GSFLOW.py

echo "press Enter"

