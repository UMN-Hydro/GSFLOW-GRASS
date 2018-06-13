#!/bin/bash
set -e # stop on error in execution

# full pathname for location of GSFLOWGRASS_toolkit
toolkit_dir="/home/gcng/workspace/ProjectFiles/GSFLOW-GRASS"

# full pathname for Settings File
Settings_file="/home/gcng/workspace/ProjectFiles/GSFLOW-GRASS/Run/settings_template_spinup_mni_1p2p0.ini"

# ------------------------------------------------------------------------------
# **For default implementation, user does not need to make changes to below lines**

pwd0=`pwd`

# (go to directory with input file builder scripts) 
cd ${toolkit_dir}/input_file_builder  

# - If fl_create_hydcond=1 in Settings File, this creates a file with spatially distributed hydraulic conductivity.
python2.7 createSpatialHydCond.py ${Settings_file}

# - If fl_print_climate_hru=1 in Settings File, this creates climate_hru files with uniform conditions 
python2.7 printClimatehru.py ${Settings_file}

# - This creates GSFLOW control file.  See Python script to make changes beyond parameters set in Settings File.
python2.7 printGSFLOWControlfile.py ${Settings_file}

# - This creates PRMS parameter file.  See Python script to make changes beyond parameters set in Settings File.
python2.7 printPRMSparamfile.py ${Settings_file}

# - This creates the various MODFLOW input files.  See Python script to make changes beyond parameters set in Settings File.
python2.7 printMODFLOWInputs.py ${Settings_file}

# (go to directory with run script)
cd ${toolkit_dir}/Run  

# - This executes GSFLOW model with the input files created above
python2.7 runGSFLOW.py ${Settings_file}

cd $pwd0


