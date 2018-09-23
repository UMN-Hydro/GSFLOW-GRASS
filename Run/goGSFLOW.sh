#!/bin/bash
set -e # stop on error in execution

# USAGE: ./goGSFLOW.sh [$SETTINGS_FILE_VAR] [$GSFLOW-GRASS_DIR]

# full pathname for Settings File
if [ -z "$1" ]
then
    Settings_file="/home/awickert/models/GSFLOW-GRASS/examples/gridTestCannon.ini"
else
    Settings_file="$1"
fi
echo Settings_file = $Settings_file

# full pathname for location of GSFLOWGRASS_toolkit
if [ -z "$2" ]
then
    toolkit_dir="/home/awickert/models/GSFLOW-GRASS"
else
    toolkit_dir="$2"
fi
echo toolkit_dir = $toolkit_dir

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


