# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 22:06:52 2017

Based on: print_MODFLOW_inputs_res_NWT.m

@author: gcng
"""

# print_MODFLOW_inputs

from MODFLOW_scripts import MODFLOW_NWT_lib_current as mf # functions to write individual MODFLOW files
import settings_test
import platform
import datetime as dt

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'


# - directories
sw_2005_NWT = 2 # 1 for MODFLOW-2005; 2 for MODFLOW-NWT algorithm (both can be 
                # carried out with MODFLOW-NWT code) 
fl_BoundConstH = 0 # 1 for const head at high elev boundary, needed for numerical 
                    # convergence for AGU2016 poster.  Seems resolved with MODFLOW-NWT

MODFLOW_indir = settings_test.MODFLOWinput_dir + slashstr 
MODFLOW_indir_rel = '..' + slashstr + settings_test.MODFLOWinput_dir_rel + slashstr 
MODFLOW_outdir_rel = '..' + slashstr + settings_test.MODFLOWoutput_dir_rel + slashstr 

infile_pre = settings_test.PROJ_CODE

#NLAY = 1;
#DZ = [200] # [NLAYx1] [m] ***testing
## DZ = [350, 100] # [NLAYx1] [m] ***testing
NLAY = settings_test.NLAY
DZ = settings_test.DZ


# length of transient stress period (follows 1-day steady-state period, ok if longer than actual GSFLOW period) [d]
start_date = dt.datetime.strptime(settings_test.START_DATE, "%Y-%m-%d")
end_date = dt.datetime.strptime(settings_test.END_DATE, "%Y-%m-%d")
delt = end_date - start_date
perlen_tr = delt.days

GIS_indir = settings_test.GISinput_dir + slashstr

# use restart file as initial cond (empty string to not use restart file)
fil_res_in = '' # empty string to not use restart file
if settings_test.sw_1spinup_2restart == 2:
    fil_res_in = settings_test.restart_MODfil

# for various files: ba6, dis, uzf, lpf
surfz_fil = GIS_indir + settings_test.DEM + '.asc'
# surfz_fil = GIS_indir + 'SRTM_new_20161208.asc'
# for various files: ba6, uzf
mask_fil = GIS_indir + 'basin_mask.asc'
# for ba6 (pour point)
dischargept_fil = GIS_indir + 'pp_tmp.txt'

# for sfr
reach_fil = GIS_indir + 'reaches_tmp.txt'
segment_fil_all = [GIS_indir + 'segments_tmp_4A_INFORMATION.txt', 
                   GIS_indir + 'segments_tmp_4B_UPSTREAM.txt', 
                   GIS_indir + 'segments_tmp_4C_DOWNSTREAM.txt']


## 
dis_fil = mf.write_dis_MOD2_f(MODFLOW_indir, infile_pre, surfz_fil, NLAY, DZ, perlen_tr);
ba6_fil = mf.write_ba6_MOD3_2(MODFLOW_indir, infile_pre, mask_fil, dischargept_fil, dis_fil, fl_BoundConstH); # list this below write_dis_MOD2_f

# flow algorithm
if sw_2005_NWT == 1:
    lpf_fil = mf.write_lpf_MOD2_f2_2(MODFLOW_indir, infile_pre, surfz_fil, NLAY);
elif sw_2005_NWT == 2:
    # MODFLOW-NWT files
    upw_fil = mf.write_upw_MOD2_f2_2(MODFLOW_indir, infile_pre, surfz_fil, NLAY);
    nwt_fil = mf.NWT_write_file(MODFLOW_indir, infile_pre);

# unsat zone and streamflow input files
uzf_fil = mf.make_uzf3_f_2(MODFLOW_indir, infile_pre, surfz_fil, dischargept_fil, ba6_fil); # list this below write_ba6_MOD3_2
sfr_fil = mf.make_sfr2_f_Mannings(MODFLOW_indir, infile_pre, reach_fil, dis_fil, segment_fil_all); # list this below write_dis_MOD2_f

# Write PCG file (only used for MODFLOW-2005, but this function also creates OC file)
mf.write_OC_PCG_MOD_f(MODFLOW_indir, infile_pre, perlen_tr);

# Write namefile
nam_fil = mf.write_nam_MOD_f2_NWT(MODFLOW_indir, MODFLOW_indir_rel, MODFLOW_outdir_rel, infile_pre, fil_res_in, sw_2005_NWT);

