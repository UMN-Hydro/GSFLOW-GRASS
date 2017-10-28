
SET preproc_dir=D:\GSFLOW\Pre_processor_codes\Shullcas_5\python_scripts\
SET input_file=D:\GSFLOW\Pre_processor_codes\Shullcas_5\python_scripts\settings_new.ini
SET python2.7=C:\Users\lsomers\Anaconda2\python.exe

:: ----------------------------------------------------------------------------
:: (Do not change below this line)
copy %input_file% settings.ini

%python2.7% %preproc_dir%settings_test.py 
%python2.7% %preproc_dir%Create_hydcond_array.py
%python2.7% %preproc_dir%GSFLOW_print_data_climatehru_files1.py
%python2.7% %preproc_dir%GSFLOW_print_controlfile_current.py
%python2.7% %preproc_dir%GSFLOW_print_PRMSparamfile_current.py
%python2.7% %preproc_dir%print_MODFLOW_inputs_res_NWT_current.py
%python2.7% %preproc_dir%run_GSFLOW.py
pause

