#!/bin/bash

#gsflow_exe=$HOME"/opt/GSFLOW_1.2.1\ Linux/bin"
#gsflow_exe=$HOME"/opt/GSFLOW-1.2.0/bin"
gsflow_exe="/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/bin"
proj_name="ChimbNew"

local_dir=${PWD}

settings_ini="$local_dir/settings.ini"

gsflow_dir="$local_dir/GSFLOW"

control_dir="$gsflow_dir/control/"
PRMSinput_dir="$gsflow_dir/inputs/PRMS/"
MODFLOWinput_dir="$gsflow_dir/inputs/MODFLOW_NWT/"
PRMSoutput_dir="$gsflow_dir/outputs/PRMS/"

DEM="stros_srtm_3arc_utm"

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
cp -R /home/gcng/workspace/ProjectFiles/AndesWaterResources/Data2/Chimb_New_MoreGridExtent
#cp strosa_data/dem/* data/GIS
#cp strosa_data/climate/* "$PRMSinput_dir"

touch python_scripts/MODFLOW_scripts/__init__.py

python2.7 python_scripts/settings.py
##python2.7 process_weatherdata.py
python2.7 python_scripts/GSFLOW_print_controlfile_res.py
python2.7 python_scripts/GSFLOW_print_PRMSparamfile.py
python2.7 python_scripts/print_MODFLOW_inputs_res_NWT_res.py
./GSFLOW/control/ChimbNew_GSFLOW.sh

