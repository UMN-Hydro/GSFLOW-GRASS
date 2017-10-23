# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 22:06:52 2017

Based on: print_MODFLOW_inputs_res_NWT.m

@author: gcng
"""

# print_MODFLOW_inputs

import numpy as np
from MODFLOW_scripts import MODFLOW_NWT_lib_Shullcas_test as mf # functions to write individual MODFLOW files
import os  # os functions
#from ConfigParser import SafeConfigParser
import settings_test
import platform

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'

#GSFLOW_DIR = settings_test.LOCAL_DIR + "/GSFLOW/"
GSFLOW_DIR = settings_test.gsflow_simdir + "/GSFLOW/"

# - directories
sw_2005_NWT = 2 # 1 for MODFLOW-2005; 2 for MODFLOW-NWT algorithm (both can be 
                # carried out with MODFLOW-NWT code) 
fl_BoundConstH = 0 # 1 for const head at high elev boundary, needed for numerical 
                    # convergence for AGU2016 poster.  Maybe resolved with MODFLOW-NWT?

#if sw_2005_NWT == 1:
#    # MODFLOW input files
#    GSFLOW_indir = GSFLOW_DIR + slashstr + 'inputs' + slashstr + 'MODFLOW_2005' + slashstr
#    # MODFLOW output files
#    GSFLOW_outdir = GSFLOW_DIR + slashstr + 'outputs' + slashstr + 'MODFLOW_2005' + slashstr
#    
#elif sw_2005_NWT == 2:
#    # MODFLOW input files
#    GSFLOW_indir = GSFLOW_DIR + slashstr + 'inputs' + slashstr + 'MODFLOW_NWT' + slashstr
#    # MODFLOW output files
#    GSFLOW_outdir = GSFLOW_DIR + slashstr + 'outputs' + slashstr + 'MODFLOW_NWT' + slashstr
MODFLOW_indir = settings_test.MODFLOWinput_dir + slashstr ### Change here for MODFLOW_indir
MODFLOW_outdir = settings_test.MODFLOWoutput_dir + slashstr ### Change here for MODFLOW_indir

infile_pre = settings_test.PROJ_CODE

NLAY = 1;
DZ = [200] # [NLAYx1] [m] ***testing
# DZ = [350, 100] # [NLAYx1] [m] ***testing

# length of transient stress period (follows 1-day steady-state period) [d]
# perlen_tr = 365; # [d], ok if too long
# perlen_tr = 365*5 + ceil(365*5/4); # [d], includes leap years; ok if too long (I think, but maybe run time is longer?)
perlen_tr =  1131 + 1 # [d], includes leap years; ok if too long (I think, but maybe run time is longer?)

GIS_indir = settings_test.GISinput_dir + slashstr

# use restart file as initial cond (empty string to not use restart file)
fil_res_in = '' # empty string to not use restart file
#fil_res_in = '/home/gcng/workspace/Pfil_res_inrojectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/test2lay_melt_30yr.out' % empty string to not use restart file

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


# create MODFLOW input directory if it does not exist:
if not os.path.isdir(MODFLOW_indir):
    os.makedirs(MODFLOW_indir)
    
# while we're at it, create MODFLOW output file if it does not exist:
if not os.path.isdir(MODFLOW_outdir):
    os.makedirs(MODFLOW_outdir)

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
nam_fil = mf.write_nam_MOD_f2_NWT(MODFLOW_indir, MODFLOW_outdir, infile_pre, fil_res_in, sw_2005_NWT);

