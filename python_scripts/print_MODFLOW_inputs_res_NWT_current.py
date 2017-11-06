# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 22:06:52 2017

Based on: print_MODFLOW_inputs_res_NWT.m

@author: gcng
"""

# print_MODFLOW_inputs

from MODFLOW_scripts import MODFLOW_NWT_lib_current as mf # functions to write individual MODFLOW files
import datetime as dt
from readSettings import Settings
import platform
import sys

# Set input file
if len(sys.argv) < 2:
    settings_input_file = 'settings.ini'
    print 'Using default input file: ' + settings_input_file
else:
    settings_input_file = sys.argv[1]
    print 'Using specified input file: ' + settings_input_file
Settings = Settings(settings_input_file)

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'


# - directories
sw_2005_NWT = 2 # 1 for MODFLOW-2005; 2 for MODFLOW-NWT algorithm (both can be 
                # carried out with MODFLOW-NWT code) 

MODFLOW_indir = Settings.MODFLOWinput_dir + slashstr 
MODFLOW_indir_rel = '..' + slashstr + Settings.MODFLOWinput_dir_rel + slashstr 
MODFLOW_outdir_rel = '..' + slashstr + Settings.MODFLOWoutput_dir_rel + slashstr 

infile_pre = Settings.PROJ_CODE

#NLAY = 1;
#DZ = [200] # [NLAYx1] [m] ***testing
## DZ = [350, 100] # [NLAYx1] [m] ***testing
NLAY = Settings.NLAY
DZ = Settings.DZ


# length of transient stress period (follows 1-day steady-state period, ok if longer than actual GSFLOW period) [d]
start_date = dt.datetime.strptime(Settings.START_DATE, "%Y-%m-%d")
end_date = dt.datetime.strptime(Settings.END_DATE, "%Y-%m-%d")
delt = end_date - start_date
perlen_tr = delt.days

GIS_indir = Settings.GISinput_dir + slashstr

# use restart file as initial cond (empty string to not use restart file)
fil_res_in = '' # empty string to not use restart file
if Settings.sw_1spinup_2restart == 2:
    fil_res_in = Settings.restart_MODfil

# for various files: ba6, dis, uzf, lpf
surfz_fil = GIS_indir + 'DEM.asc'
# surfz_fil = GIS_indir + 'SRTM_new_20161208.asc'
# for various files: ba6, uzf
mask_fil = GIS_indir + 'basin_mask.asc'
# for ba6 (pour point)
dischargept_fil = GIS_indir + 'pour_point.txt'

# for sfr
reach_fil = GIS_indir + 'reaches.txt'
segment_fil_all = [GIS_indir + 'segments_4A_INFORMATION.txt', 
                   GIS_indir + 'segments_4B_UPSTREAM.txt', 
                   GIS_indir + 'segments_4C_DOWNSTREAM.txt']


## 
dis_fil = mf.write_dis_MOD2_f(MODFLOW_indir, infile_pre, surfz_fil, NLAY, DZ, perlen_tr);
ba6_fil = mf.write_ba6_MOD3_2(MODFLOW_indir, infile_pre, mask_fil, dischargept_fil, dis_fil); # list this below write_dis_MOD2_f

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

