#!/bin/bash

preproc_dir="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor"
gsflow_simdir="/media/gcng/STORAGE3A/GSFLOW/Shullcas_CrystalLauren5_works_res"
proj_code="Shullcas_res"

python2.7 ${preproc_dir}/python_scripts/settings_test.py
python2.7 ${preproc_dir}/python_scripts/Create_hydcond_array.py
python2.7 ${preproc_dir}/python_scripts/GSFLOW_print_data_climatehru_files1.py
python2.7 ${preproc_dir}/python_scripts/GSFLOW_print_controlfile_current.py
python2.7 ${preproc_dir}/python_scripts/GSFLOW_print_PRMSparamfile_current.py
python2.7 ${preproc_dir}/python_scripts/print_MODFLOW_inputs_res_NWT_current.py

${gsflow_simdir}/control/${proj_code}_GSFLOW.sh

