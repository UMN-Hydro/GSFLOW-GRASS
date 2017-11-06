#!/bin/bash

#gsflow_exe=$HOME"/opt/GSFLOW_1.2.1\ Linux/bin"
#gsflow_exe=$HOME"/opt/GSFLOW-1.2.0/bin"
gsflow_exe="/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/bin"
proj_name="Shullcas"

local_dir=${PWD}

settings_ini="$local_dir/settings.ini"

preproc_dir="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor"

gsflow_dir="$local_dir/GSFLOW"

control_dir="$gsflow_dir/control/"
PRMSinput_dir="$gsflow_dir/inputs/PRMS/"
MODFLOWinput_dir="$gsflow_dir/inputs/MODFLOW_NWT/"
PRMSoutput_dir="$gsflow_dir/outputs/PRMS/"

# name of DEM files (w/ .asc suffix)
DEM="DEM"

echo "[settings]" > "$settings_ini"
echo "proj_name=$proj_name" >> "$settings_ini"
echo "local_dir=$local_dir" >> "$settings_ini"
echo "gsflow_exe=$gsflow_exe" >> "$settings_ini"

echo "control_dir=$control_dir" >> "$settings_ini"
echo "PRMSinput_dir=$PRMSinput_dir" >> "$settings_ini"
echo "MODFLOWinput_dir=$MODFLOWinput_dir" >> "$settings_ini"
echo "PRMSoutput_dir=$PRMSoutput_dir" >> "$settings_ini"
echo "DEM=$DEM" >> "$settings_ini"

mkdir -p "$control_dir" "$PRMSinput_dir" "$MODFLOWinput_dir" "$PRMSoutput_dir"


mkdir -p data/GIS
#cp -R strosa_data/GIS/* data/GIS/
#cp strosa_data/dem/* data/GIS
#cp strosa_data/climate/* "$PRMSinput_dir"
#cp -R /home/gcng/workspace/ProjectFiles/AndesWaterResources/Data2/Chimb_Newgrid/* data/GIS/
#cp strosa_data/dem/* data/GIS
#cp strosa_data/climate/* "$PRMSinput_dir"
cp /media/gcng/MY\ PASSPORT/GSFLOW/Shullcas_FirstTry_Sunday_night/Data/GIS/* data/GIS/
cp /media/gcng/MY\ PASSPORT/GSFLOW/Shullcas_FirstTry_Sunday_night/inputs/PRMS/*.day "$PRMSinput_dir"

touch ${preproc_dir}/python_scripts/MODFLOW_scripts/__init__.py

python2.7 ${preproc_dir}/python_scripts/settings.py
##python2.7 process_weatherdata.py
python2.7 ${preproc_dir}/python_scripts/GSFLOW_print_controlfile_Shullcas.py
python2.7 ${preproc_dir}/python_scripts/GSFLOW_print_PRMSparamfile_Shullcas.py
python2.7 ${preproc_dir}/python_scripts/print_MODFLOW_inputs_res_NWT_Shullcas.py
./GSFLOW/control/${proj_name}_GSFLOW.sh

