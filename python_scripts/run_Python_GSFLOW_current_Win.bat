
preproc_dir="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/"
input_file="/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/python_scripts/settings_new.ini"

:: ----------------------------------------------------------------------------
:: (Do not change below this line)
copy %input_file% settings.ini

python2.7 %preproc_dir%settings_test.py 
python2.7 %preproc_dir%Create_hydcond_array.py
python2.7 %preproc_dir%GSFLOW_print_data_climatehru_files1.py
python2.7 %preproc_dir%GSFLOW_print_controlfile_current.py
python2.7 %preproc_dir%GSFLOW_print_PRMSparamfile_current.py
python2.7 %preproc_dir%print_MODFLOW_inputs_res_NWT_current.py
python2.7 %preproc_dir%run_GSFLOW.py


