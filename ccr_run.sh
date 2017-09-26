#!/bin/bash

#set -vx

# (3)
# this needs to be set by the user
gsflow_exe=$HOME"/opt/GSFLOW_1.2.1 Linux/bin"
spin_up="30yr"

# (2)
# in order to fix Permission denied problems after (1)
local_dir=${PWD}

echo "[settings]" > settings.ini
echo "local_dir=$local_dir" >> settings.ini
echo "gsflow_exe=$gsflow_exe" >> settings.ini
echo "spin_up=$spin_up" >> settings.ini

gsflow_dir="$local_dir/GSFLOW"

control_dir="$gsflow_dir/control/"
PRMSinput_dir="$gsflow_dir/inputs/PRMS/"
MODFLOWinput_dir="$gsflow_dir/inputs/MODFLOW_NWT/"
#MODFLOWinput_dir="$gsflow_dir/inputs/MODFLOW_2005/'
PRMSoutput_dir="$gsflow_dir/outputs/PRMS/"

# directory with files to be read in to generate PRMS input files (to get
# start/end dates)
in_data_dir="$gsflow_dir/DataToReadIn/"

mkdir -p "$control_dir" "$PRMSinput_dir" "$MODFLOWinput_dir" "$PRMSoutput_dir" "$in_data_dir"

# (4)
tar -xvzf Data_GIS.tar.gz
mv GIS "$in_data_dir/"

tar -xvzf GSFLOW_inputs_PRMS_climate.tar.gz
mv climate "$PRMSinput_dir/"



# (1)
python2.7 "python_scripts/GSFLOW_print_controlfile4_gcng_melt$spin_up.py"
# tried to run this, needed python-tk
# sudo apt-get install python-tk
# then received error
# "Permission denied: '/home/gcng'"
# see (2)

python2.7 python_scripts/GSFLOW_print_PRMSparamfile4.py
python2.7 python_scripts/GSFLOW_print_data_climatehru_files1.py
python2.7 python_scripts/MODFLOW_scripts/print_MODFLOW_inputs_res_NWT.py 

"${gsflow_exe}/gsflow" "$gsflow_dir/control/py_ChimTest_Melt_rep${spin_up}_GSFLOW.control" &> out.txt
