#!/bin/bash

#set -vx

# (3)
# this needs to be set by the user
gsflow_exe=$HOME"/opt/GSFLOW_1.2.1\ Linux/bin"

# (2)
# in order to fix Permission denied problems after (1)
local_dir=${PWD}

echo "[settings]" > settings.ini
echo "local_dir=$local_dir" >> settings.ini
echo "gsflow_exe=$gsflow_exe" >> settings.ini

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


# (1)
python2.7 python_scripts/GSFLOW_print_controlfile4_gcng_melt30yr.py
# tried to run this, needed python-tk
# sudo apt-get install python-tk
# then received error
# "Permission denied: '/home/gcng'"
# see (2)
