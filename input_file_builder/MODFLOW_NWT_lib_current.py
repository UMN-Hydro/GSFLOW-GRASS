# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 22:34:30 2017

MODFLOW_NWT_lib

This library includes all the separate matlab functions for writing the different MODFLOW input files

@author: gcng
"""

import numpy as np
import pandas as pd # for data structures and reading in data from text file
import matplotlib.pyplot as plt # matlab-like plots
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


# function for parsing ASCII grid header in GIS data files
def read_grid_file_header(fname):
    f = open(fname, 'r')
    sdata = {}
    for i in range(6):
        line = f.readline()
        line = line.rstrip() # remove newline characters
        key, value = line.split(': ')
        try:
          value = int(value)
        except:
          value = float(value)
        sdata[key] = value
    f.close()

    return sdata


#%%
# baesd on: write_nam_MOD_f2_NWT.m

def write_nam_MOD_f2_NWT(GSFLOW_indir, GSFLOW_indir_rel, GSFLOW_outdir_rel, infile_pre, fil_res_in, sw_2005_NWT):
# v2 - allows for restart option (init)
# _NWT: from Leila's email, her work from spring 2017, incorporates
# MODFLOW-NWT.
# Edits made to Leila's version, which incorrectly called both (LPF, PCG) and UPW.
# MODFLOW-NWT will default to MODFLOW-2005 if given LPF and PCG.
#
# sw_2005_NWT: 1 for MODFLOW-2005; 2 for MODFLOW-NWT algorithm (both can be 
#              carried out with MODFLOW-NWT code)

# # - directories
# # MODFLOW input filesfil_res_in
# GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW'
# # MODFLOW output files
# GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
# infile_pre = 'test2lay_py'

    # - write to this file (within indir)
    fil_nam = infile_pre + '.nam'
    
    # all assumed to be in GSFLOW_dir
    fil_ba6 = infile_pre +'.ba6'
    if sw_2005_NWT == 1:
        fil_lpf = infile_pre +'.lpf' 
        fil_pcg = infile_pre +'.pcg' 
    elif sw_2005_NWT == 2:
        fil_upw = infile_pre +'.upw'         
        fil_nwt = infile_pre +'.nwt' 
    fil_oc = infile_pre +'.oc' 
    fil_dis = infile_pre +'.dis' 
    fil_uzf = infile_pre +'.uzf' 
    fil_sfr = infile_pre +'.sfr' 
    fil_res_out = infile_pre +'.out'  # write to restart file
    
    # ------------------------------------------------------------------------
    
    # -- .nam file with full paths
    fil_nam_0 = GSFLOW_indir + slashstr + fil_nam
    fobj = open(fil_nam_0, 'w+')
    fobj.write('LIST          7 ' + GSFLOW_outdir_rel + 'test.lst \n') # MODFLOW output file
    fobj.write('BAS6          8 ' + GSFLOW_indir_rel + fil_ba6 + '\n')
    if sw_2005_NWT == 1:    
        fobj.write('LPF          11 ' + GSFLOW_indir_rel + fil_lpf + '\n')
        fobj.write('PCG          19 ' + GSFLOW_indir_rel + fil_pcg + '\n')
    elif sw_2005_NWT == 2:
        fobj.write('UPW          25 ' + GSFLOW_indir_rel + fil_upw + '\n')
        fobj.write('NWT          17 ' + GSFLOW_indir_rel + fil_nwt + '\n')
    fobj.write('OC           22 ' + GSFLOW_indir_rel + fil_oc + '\n')
    fobj.write('DIS          10 ' + GSFLOW_indir_rel + fil_dis + '\n')
    fobj.write('UZF          12 ' + GSFLOW_indir_rel + fil_uzf + '\n')
    fobj.write('SFR          13 ' + GSFLOW_indir_rel + fil_sfr + '\n')
    if len(fil_res_in) != 0:
        fobj.write('IRED         90 ' + fil_res_in + '\n');
    fobj.write('IWRT         91 ' + GSFLOW_outdir_rel + fil_res_out + '\n')
    fobj.write('DATA(BINARY) 34 ' + GSFLOW_outdir_rel + infile_pre + '.bud \n'); # MODFLOW LPF output file, make sure 34 is unit listed in lpf file!!
    fobj.write('DATA(BINARY) 51 ' + GSFLOW_outdir_rel + infile_pre + '_head.bhd \n') # MODFLOW output file
    fobj.write('DATA(BINARY) 61 ' + GSFLOW_outdir_rel + infile_pre + '_uzf.dat \n') # MODFLOW output file
    fobj.write('DATA         52 ' + GSFLOW_outdir_rel + infile_pre + '_ibound.dat \n') # MODFLOW output file    
    fobj.close()

#%%
# based on: write_dis_MOD2_f
# write_dis_MOD (for 3D domains)
# 11/17/16
#
# v1 - 11/30/16 start to include GIS data for Chimborazo's Gavilan Machay
#      watershed; topo.asc for surface elevation (fill in bottom elevation
#      based on uniform thickness of single aquifer)
def write_dis_MOD2_f(GSFLOW_indir, infile_pre, surfz_fil, NLAY, DZ, perlen_tr):

# # ==== TO RUN AS SCRIPT ===================================================
#     # - directories
#     # MODFLOW input files
#    GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#    # MODFLOW output files
#    GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
#     
#     # infile_pre = 'test1lay';
#     # NLAY = 1;
#     # DZ = 10; # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
#     
#    infile_pre = 'test2lay_py'
#    NLAY = 2
#    DZ = [100, 50] # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
#     
#    perlen_tr = 365 # ok if too long
#     
#    GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/'
#     
#     # for various files: ba6, dis, uzf, lpf
#    surfz_fil = GIS_indir + 'topo.asc'
#    # surfz_fil = GIS_indir + 'SRTM_new_20161208.asc'
#     
#     # for various files: ba6, uzf
#    mask_fil = GIS_indir + 'basinmask_dischargept.asc'
##  =========================================================================


    # - write to this file
    # GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
    dis_fil = infile_pre + '.dis'
    
    # - read in this file for surface elevation (for TOP(NROW,NCOL))
    # surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
    
    # - read in this file for elevation of layer bottoms (for BOTM(NROW,NCOL,NLAY))
    # (layer 1 is top layer)
    botmz_fil = ''
    
    # - domain dimensions, maybe already in surfz_fil and botm_fil{}?
    # NLAY = 1;
    # NROW = 1058;
    # NCOL = 1996;
    
    # # - domain boundary (UTM zone 17S, outer boundaries)
    # north = 9841200;
    # south = 9835900;
    # east = 751500;
    # west = 741500;
    
    # # - space discretization
    # DELR = (east-west)/NCOL; # width of column [m]
    # DELC = (north-south)/NROW; # height of row [m]
    # DZ = 10; # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
    # DZ = [5; 5]; # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
    
    # - time discretization
    PERLEN = [1, perlen_tr];  # 2 periods: 1-day steady-state and multi-day transient
    
    comment1 = '# test file for Gavilan Machay'
    comment2 = '# test file'
    
    # - The following will be assumed:
    LAYCBD = np.zeros((1,NLAY), int); # no confining layer below layer
    ITMUNI = 4; # [d]
    LENUNI = 2; # [m]
    NPER = 2; # 1 SS then 1 transient
    NSTP = PERLEN 
    TSMULT = 1  # must have daily time step to correspond with PRMS
    SsTr_flag = ['ss', 'tr']
    
    ## ------------------------------------------------------------------------
    # -- Read in data from files
    sdata = read_grid_file_header(surfz_fil)
        
    NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
    NROW = sdata['rows'] 
    NCOL = sdata['cols']

    # - space discretization
    DELR = (NSEW[2]-NSEW[3])/NCOL # width of column [m]
    DELC = (NSEW[0]-NSEW[1])/NROW # height of row [m]
    
    TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)
    
    BOTM = np.zeros((NROW, NCOL, NLAY), float)
    BOTM[:,:,0] = TOP-DZ[0];
    for ilay in range(1,NLAY):
        BOTM[:,:,ilay] = BOTM[:,:,ilay-1]-DZ[ilay]

    # -- Discretization file:
    dis_fil_0 = GSFLOW_indir + slashstr + dis_fil
    fobj = open(dis_fil_0, 'w+')
#    fobj = open('/home/gcng/test.txt', 'w')
    fobj.write(comment1 + '\n');
    fobj.write(comment2 + '\n');    
    _out = np.array([NLAY, NROW, NCOL, NPER, ITMUNI, LENUNI]) # for 1D arrays
    np.savetxt(fobj, _out[None], delimiter=' ', 
               fmt='%5s %5s %5s %5s %5s %5s   NLAY, NROW, NCOL, NPER, ITMUNI, LENUNI')
    np.savetxt(fobj, LAYCBD, delimiter=' ', fmt='%d')

#    fobj.write('CONSTANT %7.3f   DELR\n' % (DELR))
#    fobj.write('CONSTANT %7.3f   DELC\n' % (DELC))
    fobj.write('CONSTANT %g   DELR\n' % (DELR))
    fobj.write('CONSTANT %g   DELC\n' % (DELC))
    fobj.write('INTERNAL  1.0  (FREE)  0                          TOP ELEVATION OF LAYER 1 \n');
    np.savetxt(fobj, TOP, delimiter=' ', fmt='%10g')
    
    for ii in range(NLAY):
        fobj.write('INTERNAL  1.0  (FREE)  0                          BOTM ELEVATION OF LAYER %d \n' 
        % (ii+1))
        np.savetxt(fobj, BOTM[:,:,ii], delimiter=' ', fmt='%10g')        

    for ii in range(NPER):
        fobj.write(' %g %d %g %s        PERLEN, NSTP, TSMULT, Ss/Tr (stress period %4d)\n' 
        % (PERLEN[ii], NSTP[ii], TSMULT, SsTr_flag[ii], ii+1))

    fobj.close()
    
    # -- plot domain discretization
    TOP_to_plot = TOP.copy()
    TOP_to_plot[TOP <= 0] == np.nan
    BOTM_to_plot = BOTM.copy()
    BOTM_to_plot[BOTM <= 0] == np.nan

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(2,2,1)
    im = ax.imshow(TOP_to_plot)
    im.set_clim(3800, 6200)
    # use im.get_clim() to get the limits and generalize
    fig.colorbar(im, orientation='horizontal')
    plt.title('TOP')
    for ilay in range(NLAY):
        plt.subplot(2,2,2+ilay)
        im = plt.imshow(BOTM_to_plot[:,:,ilay])
        im.set_clim(3800, 6200)
        fig.colorbar(im, orientation='horizontal')
        plt.title('BOTM lay' + str(ilay+1));
        
    return dis_fil_0

#    
#    figure
#    for ilay = 1:NLAY
#        subplot(2,2,double(ilay))
#        if ilay == 1
#            X = TOP - BOTM(:,:,ilay);
#        else
#            X = BOTM(:,:,ilay-1)-BOTM(:,:,ilay);
#        end
#    #     m = X(X>0); m = min(m(:));
#        imagesc(X), #caxis([m*0.9, max(X(:))]), 
#        cm = colormap;
#        cm(1,:) = [1 1 1];
#        colormap(cm);
#        colorbar
#        title(['DZ', ' lay', num2str(ilay)]);
#    end
    



#%%
# based on: write_ba6_MOD3_2.m
# (had to be careful of numerical convergence problems; set constant head for 
# outer boundary to avoid these.  Later resolved with NWT by Leila)

def write_ba6_MOD3_2(GSFLOW_indir, infile_pre, mask_fil, dischargept_fil, dis_fil):

#    # ==== TO RUN AS SCRIPT ===================================================
#    # - directories
#    # MODFLOW input files
#    GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#    # MODFLOW output files
#    GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
#    slashstr = '/'
#    
#    # infile_pre = 'test1lay';
#    # NLAY = 1;
#    # DZ = 10; # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
#    
#    infile_pre = 'test2lay_py'
#    NLAY = 2
#    DZ = [100, 50] # [NLAYx1] [m] ***testing
#    
#    GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/'
#    
#    # for various files: ba6, dis, uzf, lpf
#    surfz_fil = GIS_indir + 'topo.asc'
#    # for various files: ba6, uzf
#    mask_fil = GIS_indir + 'basinmask_dischargept.asc'
#
#    fl_BoundConstH = 1 # 1 for const head at high elev boundary, needed for 
#    numerical convergence for AGU2016 poster.  Maybe resolved with MODFLOW-NWT?
#    # =========================================================================
    
    # - write to this file
    # GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
    ba6_file = infile_pre + '.ba6'
    
    # - domain dimensions, maybe already in surfz_fil and botm_fil{}?
    # NLAY = 1;
    # NROW = 50;
    # NCOL = 50;
    
    # -- IBOUND(NROW,NCOL,NLAY): <0 const head, 0 no flow, >0 variable head
    # use basin mask (set IBOUND>0 within watershed, =0 outside watershed, <0 at discharge point and 2 neighboring pixels)
    # mask_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/basinmask_dischargept.asc';
    sdata = read_grid_file_header(mask_fil)
    NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
    NROW = sdata['rows'] 
    NCOL = sdata['cols']
    IBOUND = np.genfromtxt(mask_fil, skip_header=6, skip_footer=0, delimiter=' ', dtype=float)

#    f = open(dischargept_fil, 'r')
#    last_line = f.readlines()
#    last_line = last_line[-1].rstrip()
#    f.close()    
#    value1, value2 = last_line.split(': ')
#    value2 = value2.split(' ')
#    dischargePt_rowi = int(value2[1])
#    dischargePt_coli = int(value2[3])

    f = open(dischargept_fil, 'r')
    ctr = 0
    for line in f:
        line = line.rstrip()
        value1, value2 = line.split(': ')
        value2 = value2.split(' ')
        if ctr == 0:
            dischargePt_rowi = int(value2[1])
            dischargePt_coli = int(value2[3])
        else:
            DowngradPt_rowi = int(value2[1])
            DowngradPt_coli = int(value2[3])            
        ctr = ctr + 1
    f.close()    
      
    # - force some cells to be active to correspond to stream reaches
#    print "Warning!!  Hard-coded to set some IBOUND values to be active!! Check Andy's GIS algorithm..."
    #IBOUND[14-1,33-1] = 1
    #IBOUND[11-1,35-1] = 1
    #IBOUND[12-1,34-1] = 1
    #IBOUND[7-1,43-1] = 1
    
    # find boundary cells
    IBOUNDin = IBOUND[1:-1-1+1,1:-1-1+1]
    IBOUNDu = IBOUND[0:-1-2+1,1:-1-1+1] # up
    IBOUNDd = IBOUND[2:,1:-1-1+1] # down
    IBOUNDl = IBOUND[1:-1-1+1,0:-1-2+1] # left
    IBOUNDr = IBOUND[1:-1-1+1,2:] # right    

    # - inner boundary is constant head
    ind_bound = ((IBOUNDin==1) & ((IBOUNDin-IBOUNDu==1) | (IBOUNDin-IBOUNDd==1) 
        | (IBOUNDin-IBOUNDl==1) | (IBOUNDin-IBOUNDr==1)))
    # - outer boundary is constant head 
    # ind_bound = ((IBOUNDin==0) & ((IBOUNDin-IBOUNDu==-1) | (IBOUNDin-IBOUNDd==-1) | 
    #     (IBOUNDin-IBOUNDl==-1) | (IBOUNDin-IBOUNDr==-1)))
    
    
    # -- init head: base on TOP and BOTM
#    dis_file = GSFLOW_indir + slashstr + infile_pre + '.dis'
    # dis_file = '/home/gcng/Shortcuts/AndesWaterResources/GSFLOW/inputs/MODFLOW/test2lay.dis'
    f = open(dis_fil, 'r')
    for i in range(3): # first 2 lines are comments
        line = f.readline().rstrip()
    line = line.split()
    NLAY = int(line[0]) 
    NROW = int(line[1]) 
    NCOL = int(line[2]) 
    NPER = int(line[3]) 
    ITMUNI = int(line[4]) 
    LENUNI = int(line[5])    

    print(NROW, NCOL)

    line = f.readline().rstrip()
    LAYCBD = np.array(line.split(), float)
    
    line = f.readline().rstrip()
    line = line.split()
    DELR = float(line[1])
    line = f.readline().rstrip()
    line = line.split()
    DELC = float(line[1])
    f.close()    
    
    TOP = np.genfromtxt(dis_fil, skip_header=7, max_rows=NROW, dtype=float)
    BOTM = np.zeros((NROW, NCOL, NLAY),float);
    for ii in range(NLAY):
        BOTM[:,:,ii] = np.genfromtxt(dis_fil, skip_header=7+(ii+1)*(NROW+1), \
        max_rows=NROW, dtype=float)

   
#    # - make discharge point and neighboring cells constant head (similar to Sagehen example)
#    IBOUND[dischargePt_rowi-1,dischargePt_coli-1] = -2 # downgrad of discharge pt
    # IBOUND[dischargePt_rowi-2,dischargePt_coli-1] = -1 # neighbor points
#    if (dischargePt_rowi < NROW):
#        IBOUND[dischargePt_rowi,dischargePt_coli-1] = -1
#    if (dischargePt_coli < NCOL):
#        IBOUND[dischargePt_rowi-1,dischargePt_coli] = -2 # downgrad of discharge pt
#        if (dischargePt_rowi > 1):
#            IBOUND[dischargePt_rowi-2,dischargePt_coli] = -1 # neighbor points
#        if (dischargePt_rowi < NROW):
#            IBOUND[dischargePt_rowi,dischargePt_coli] = -1

#    # *** SPECIFIC TO SHULLCAS
#    print "setting constant head downgradient of outlet - Shullcas"
#    if (dischargePt_coli > 1):
#        IBOUND[dischargePt_rowi-1,dischargePt_coli-2] = -2 # downgrad of discharge pt
#        if (dischargePt_rowi > 1):
#            IBOUND[dischargePt_rowi-2,dischargePt_coli-2] = -1 # neighbor points
#        if (dischargePt_rowi < NROW):
#            IBOUND[dischargePt_rowi,dischargePt_coli-2] = -1
#    
    IBOUND[dischargePt_rowi-1,dischargePt_coli-1] = 1 # active cell below discharge pt 
    IBOUND[DowngradPt_rowi-1,DowngradPt_coli-1] = -1 # constant head in cell downgrad of discharge pt
    # active cells below stream reaches!
    
    # *** SPECIFIC TO Santa Rosa
#    print "setting constant head downgradient of outlet - Santa Rosa"
#    if (dischargePt_rowi < NROW):
#        IBOUND[dischargePt_rowi,dischargePt_coli-1] = -1
#    if (dischargePt_coli < NCOL):
#        IBOUND[dischargePt_rowi-1,dischargePt_coli] = -1 # downgrad of discharge pt
#        if (dischargePt_rowi > 1):
#            IBOUND[dischargePt_rowi-2,dischargePt_coli] = -1 # neighbor points
#        if (dischargePt_rowi < NROW):
#            IBOUND[dischargePt_rowi,dischargePt_coli] = -1
#    IBOUND[dischargePt_rowi-1,dischargePt_coli-1] = 1 # discharge pt
#    # hard-code to have active cells below streams 10/2/17 (Santa Rosa)
#    IBOUND[11-1,51-1] = 1
    
    print "To do for IBOUND: check if neighboring cells to discharge point are stream cells"
    
    M = np.ones((NROW,NCOL,NLAY),float)
    for ii in range(NLAY):
        M[:,:,ii] = IBOUND
    IBOUND = M
    
    # - initHead(NROW,NCOL,NLAY) 
    initHead = BOTM[:,:,0] + (TOP-BOTM[:,:,0])*0.9; # within top layer
    # # (no more than 10m below top):
    # Y = nan(NROW,NCOL,2); Y(:,:,1) = initHead; Y(:,:,2) = TOP-10; 
    # initHead = max(Y,[],3);
    M = np.ones((NROW,NCOL,NLAY),float)
    for ii in range(NLAY):
        M[:,:,ii] = initHead
    initHead = M
    
    # - assumed values
    HNOFLO = -999.99
    
    
    ## ------------------------------------------------------------------------
    # -- Write ba6 file
    fil_ba6_0 = GSFLOW_indir + slashstr + ba6_file

    fobj = open(fil_ba6_0, 'w+')
    fobj.write('# basic package file --- %d layers, %d rows, %d columns\n' % (NLAY, NROW, NCOL))
    fobj.write('FREE\n')
    
    for ilay in range(NLAY):
        fobj.write('INTERNAL          1 (FREE)  3         IBOUND for layer %d \n' % (ilay+1))
        np.savetxt(fobj, IBOUND[:,:,ii], delimiter=' ', fmt='%4d')      
    
    fobj.write('    %f  HNOFLO\n' % (HNOFLO));
    for ilay in range(NLAY):
        fobj.write('INTERNAL          1 (FREE)  3         init head for layer %d \n' % (ilay+1))
        np.savetxt(fobj, initHead[:,:,ii], delimiter=' ', fmt='%7g')      

    fobj.close()
    
    return fil_ba6_0
    
#    # -- Plot basics
#    for ii = 1:2
#        if ii == 1, 
#            X0 = IBOUND; ti0 = 'IBOUND';
#        elseif ii == 2
#            X0 = initHead; ti0 = 'init head';
#        end
#        figure
#        for ilay = 1:NLAY
#            subplot(2,2,double(ilay))
#            X = X0(:,:,ilay);
#            m = X(X>0); m = min(m(:));
#            imagesc(X), #caxis([m*0.9, max(X(:))]), 
#            cm = colormap;
#    #         cm(1,:) = [1 1 1];
#            colormap(cm);
#            colorbar
#            title([ti0, ' lay', num2str(ilay)]);
#        end
#    end
    
#%%

# based on write_lpf_MOD2_f2_2.m

def write_lpf_MOD2_f2_2(GSFLOW_indir, infile_pre, surfz_fil, NLAY):
    
#    # =========== TO RUN AS SCRIPT ===========================================
#    # - directories
#    # MODFLOW input files
#    GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#    # MODFLOW output files
#    GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
#    
#    slashstr = '/';
#     
#    # infile_pre = 'test1lay';
#    # NLAY = 1;
#    # DZ = 10; # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
#     
#    infile_pre = 'test2lay_py'
#    NLAY = 2
#    DZ = [100, 50] # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)     
#    GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/'
#     
#    # for various files: ba6, dis, uzf, lpf
#    surfz_fil = GIS_indir + 'topo.asc'
#    # for various files: ba6, uzf
#    mask_fil = GIS_indir + 'basinmask_dischargept.asc'
#     
#    # for sfr
#    segment_fil_all = []
#    segment_fil_all.append(GIS_indir + 'segment_data_4A_INFORMATION.txt')
#    segment_fil_all.append(GIS_indir + 'segment_data_4B_UPSTREAM.txt')
#    segment_fil_all.append(GIS_indir + 'segment_data_4C_DOWNSTREAM.txt')
#
#    # ====================================================================

    ## codes pixels by distance from streams, for specifying hyd cond
    #strm_buffer_fil = GIS_indir + 'segments_buffer2.asc'
    #print 'specifying hyd cond based on distance from stream!'

    # - write to this file
    # GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
    # lpf_file = 'test.lpf';
    lpf_file = infile_pre + '.lpf'
    
    # - domain dimensions, maybe already in surfz_fil and botm_fil{}?
    # NLAY = 2;
    # surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
    sdata = read_grid_file_header(surfz_fil)
        
    NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
    NROW = sdata['rows'] 
    NCOL = sdata['cols']
    print((NROW, NCOL))

    # - space discretization
    DELR = (NSEW[2]-NSEW[3])/NCOL # width of column [m]
    DELC = (NSEW[0]-NSEW[1])/NROW # height of row [m]
    
    # - set TOP to surface elevation [m]
    TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)
        
    
    ## get strm_buffer info (pixels around streams, for hyd cond)
    #strm_buffer = np.genfromtxt(strm_buffer_fil, skip_header=6, delimiter=' ', dtype=float)

    
    # -- Base hydcond, Ss (all layers), and Sy (top layer only) on data from files
    # (temp place-holder)
    hydcond = np.ones((NROW,NCOL,NLAY))
#    print "hydcond0", hydcond0
    try:
       float(Settings.hydcond0[0])
       for lay_i in range(Settings.NLAY):
           hydcond[:,:,lay_i] = float(Settings.hydcond0[lay_i]) * hydcond[:,:,lay_i]
    except ValueError:
       for ii in range(NLAY):
            hydcond[:,:,ii] = np.genfromtxt(Settings.hydcond0, skip_header=1+ii*(NROW+1), \
            max_rows=NROW, dtype=float)
           
    Ss = 2e-6*np.ones((NROW,NCOL,NLAY),float) # constant 2e-6 /m for Sagehen
    Sy = 0.15*np.ones((NROW,NCOL,NLAY),float) # 0.08-0.15 in Sagehen (lower Sy under ridges for volcanic rocks)
    WETDRY = Sy # = Sy in Sagehen (lower Sy under ridges for volcanic rocks)
        
#    hydcond[:,:,0] = 0.01#K;
#    hydcond[:,:,1] = 0.01
#
#    hydcond[:,:,0] = 0.1#K;
#    hydcond[:,:,1] = 0.1

    
    # -- assumed input values
    flow_filunit = 34 # make sure this matches namefile!!
    hdry = 1e30  # head assigned to dry cells
    nplpf = 0    # number of LPF parameters (if >0, key words would follow)
    laytyp = np.zeros((NLAY,1),int) 
    laytyp[0] = 1;  # flag, top>0: "covertible", rest=0: "confined" (does not have WT)
    layave = np.zeros((NLAY,1),int)  # flag, layave=1: harmonic mean for interblock transmissivity
    chani = np.ones((NLAY,1),int);   # flag, chani=1: constant horiz anisotropy mult factor (for each layer)
    layvka = np.zeros((NLAY,1),int)  # flag, layvka=0: vka is vert K; >0 is vertK/horK ratio
    VKA = hydcond
    laywet = np.zeros((NLAY,1),int) 
    laywet[0]=1  # flag, 1: wetting on for top convertible cells, 0: off for confined
    fl_Tr = 1 # flag, 1 for at least 1 transient stress period (for Ss and Sy)
    WETFCT = 1.001 # 1.001 for Sagehen, wetting (convert dry cells to wet)
    IWETIT = 4 # number itermations for wetting 
    IHDWET = 0 # wetting scheme, 0: equation 5-32A is used: h = BOT + WETFCT (hn - BOT)
    
    ## ------------------------------------------------------------------------

    fmt1 = ''
    for ii in range(NLAY):
        fmt1 = fmt1 + '%2d '
    
    fil_lpf_0 = GSFLOW_indir + slashstr + lpf_file
    fobj = open(fil_lpf_0, 'w+')
    fobj.write('# LPF package inputs\n');
    fobj.write('%d %g %d    ILPFCB,HDRY,NPLPF\n' % (flow_filunit, hdry, nplpf))
    _out = np.squeeze(laytyp) # for 1D arrays
    fmt12 = fmt1 + ' LAYTYP'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(layave) # for 1D arrays
    fmt12 = fmt1 + ' LAYAVE'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(chani) # for 1D arrays
    fmt12 = fmt1 + ' CHANI'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(layvka) # for 1D arrays
    fmt12 = fmt1 + ' LAYVKA'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(laywet) # for 1D arrays
    fmt12 = fmt1 + ' LAYWET'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    

    if any(laywet == 1):
        fobj.write('%g %d %d       WETFCT, IWETIT, IHDWET\n' % (WETFCT, IWETIT, IHDWET))
    
    # -- Write HKSAT and Ss, Sy (if Tr) in .lpf file
    # loop thru layers (different entry for each layer)
    for ilay in range(NLAY):
        fobj.write('INTERNAL   1.000E-00 (FREE)    0            HY layer  %d\n' % (ilay+1))
        np.savetxt(fobj, hydcond[:,:,ilay], delimiter=' ', fmt=' %4.2e')      

        fobj.write('INTERNAL   1.000E-00 (FREE)    0            VKA layer  %d\n' % (ilay+1))
        np.savetxt(fobj, VKA[:,:,ilay], delimiter=' ', fmt=' %4.2e')      
    
        if fl_Tr:
            fobj.write('INTERNAL   1.000E-00 (FREE)    0            Ss layer  %d\n' % (ilay+1))
            np.savetxt(fobj, Ss[:,:,ilay], delimiter=' ', fmt=' %4.2e')      
            if laytyp[ilay] > 0: # convertible, i.e. unconfined
                fobj.write('INTERNAL   1.000E-00 (FREE)    0            Sy layer  %d\n' % (ilay+1))
                np.savetxt(fobj, Sy[:,:,ilay], delimiter=' ', fmt=' %4.2e')      
                if laywet[ilay] > 0:
                    fobj.write('INTERNAL   1.000E-00 (FREE)    0            WETDRY layer  %d\n' % (ilay+1))
                    np.savetxt(fobj, WETDRY[:,:,ilay], delimiter=' ', fmt=' %4.2f')      
    
    fobj.write('\n')
    fobj.close()
    
    return fil_lpf_0
    
#    figure
#    for ilay = 1:NLAY
#        subplot(2,2,double(ilay))
#        X = hydcond(:,:,ilay);
#        m = X(X>0); m = min(m(:));
#        imagesc(X), #caxis([m*0.9, max(X(:))]), 
#        cm = colormap;
#    #         cm(1,:) = [1 1 1];
#        caxis([0 max(X(:))])
#        colormap(cm);
#        colorbar
#        title(['hydcond', num2str(ilay)]);
#    end

#%%

# Based on write_upw_MOD2_f2_2.m
# (which is almost identical to write_lpf_MOD2_f2_2.m; only changes for upw are 
# to comment out WETFCT, IWETIT, and IHDWET; set LAYWET all to 0.  Some edits 
# made to Leila's script.)


def write_upw_MOD2_f2_2(GSFLOW_indir, infile_pre, surfz_fil, NLAY):
    
#    # =========== TO RUN AS SCRIPT ===========================================
#    # - directories
#    # MODFLOW input files
#    GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#    # MODFLOW output files
#    GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
#    
#    slashstr = '/';
#     
#    # infile_pre = 'test1lay';
#    # NLAY = 1;
#    # DZ = 10; # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
#     
#    infile_pre = 'test2lay_py'
#    NLAY = 2
#    DZ = [100, 50] # [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)     
#    GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/'
#     
#    # for various files: ba6, dis, uzf, lpf
#    surfz_fil = GIS_indir + 'topo.asc'
#    # for various files: ba6, uzf
#    mask_fil = GIS_indir + 'basinmask_dischargept.asc'
#     
#    # for sfr
#    segment_fil_all = []
#    segment_fil_all.append(GIS_indir + 'segment_data_4A_INFORMATION.txt')
#    segment_fil_all.append(GIS_indir + 'segment_data_4B_UPSTREAM.txt')
#    segment_fil_all.append(GIS_indir + 'segment_data_4C_DOWNSTREAM.txt')
#
#    # ====================================================================

    ## codes pixels by distance from streams, for specifying hyd cond
    #strm_buffer_fil = GIS_indir + '/segments_buffer2.asc'
    #print 'specifying hyd cond based on distance from stream!'

    # - write to this file
    # GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
    # upw_file = 'test.upw';
    upw_fil = infile_pre + '.upw'
    
    # - domain dimensions, maybe already in surfz_fil and botm_fil{}?
    # NLAY = 2;
    # surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
    sdata = read_grid_file_header(surfz_fil)
        
    NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
    NROW = sdata['rows'] 
    NCOL = sdata['cols']

    # - space discretization
    DELR = (NSEW[2]-NSEW[3])/NCOL # width of column [m]
    DELC = (NSEW[0]-NSEW[1])/NROW # height of row [m]
    
    # - set TOP to surface elevation [m]
    TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)
        
    
    ## get strm_buffer info (pixels around streams, for hyd cond)
    #strm_buffer = np.genfromtxt(strm_buffer_fil, skip_header=6, delimiter=' ', dtype=float)

    
    # -- Base hydcond, Ss (all layers), and Sy (top layer only) on data from files
    hydcond = np.ones((NROW,NCOL,NLAY))
#    print "hydcond0", Settings.hydcond0
    try:
       float(Settings.hydcond0[0])
       for lay_i in range(Settings.NLAY):
           hydcond[:,:,lay_i] = float(Settings.hydcond0[lay_i]) * hydcond[:,:,lay_i]
    except ValueError:
       for ii in range(NLAY):
            hydcond[:,:,ii] = np.genfromtxt(Settings.hydcond0, skip_header=1+ii*(NROW+1), \
            max_rows=NROW, dtype=float)
           
    Ss = 2e-6*np.ones((NROW,NCOL,NLAY),float) # constant 2e-6 /m for Sagehen
    Sy = 0.15*np.ones((NROW,NCOL,NLAY),float) # 0.08-0.15 in Sagehen (lower Sy under ridges for volcanic rocks)
    WETDRY = Sy # = Sy in Sagehen (lower Sy under ridges for volcanic rocks)
    
    
    # -- assumed input values
    flow_filunit = 34 # make sure this matches namefile!!
    hdry = 1e30  # head assigned to dry cells, so easy to detect occurrence
    nplpf = 0    # number of LPF parameters (if >0, key words would follow)
    iphdry = 1 # flag to set head to HDRY when it drops below 1e-4 LENUNI above cell bottom
    laytyp = np.zeros((NLAY,1),int) 
    laytyp[0] = 1;  # flag, top>0: "covertible", rest=0: "confined" (does not have WT)
    layave = np.zeros((NLAY,1),int)  # flag, layave=1: harmonic mean for interblock transmissivity
    chani = np.ones((NLAY,1),int);   # flag, chani=1: constant horiz anisotropy mult factor (for each layer)
    layvka = np.zeros((NLAY,1),int)  # flag, layvka=0: vka is vert K; >0 is vertK/horK ratio
    VKA = hydcond
    laywet = np.zeros((NLAY,1),int) 
#    laywet[0]=1  # flag, 1: wetting on for top convertible cells, 0: off for confined
    # ***Edited from Leila's script: MODFLOW-NWT manual says LAYWET should always 
    # be set to zero in the UPW Package because all layers with LAYTYP(NLAY)>0 are assumed to be wettable
#    laywet[0]=1  # flag, 1: wetting on for top convertible cells, 0: off for confined
    fl_Tr = 1 # flag, 1 for at least 1 transient stress period (for Ss and Sy)
    # Comment these three out to turn .lpf to .upw    
#    WETFCT = 1.001 # 1.001 for Sagehen, wetting (convert dry cells to wet)
#    IWETIT = 4 # number itermations for wetting 
#    IHDWET = 0 # wetting scheme, 0: equation 5-32A is used: h = BOT + WETFCT (hn - BOT)
    
    ## ------------------------------------------------------------------------

    fmt1 = ''
    for ii in range(NLAY):
        fmt1 = fmt1 + '%2d '
    
    upw_fil_0 = GSFLOW_indir + slashstr + upw_fil
    fobj = open(upw_fil_0, 'w+')
    fobj.write('# UPW package inputs\n');
    # Edited Leila's file, which omitted IPHDRY
    fobj.write('%d %g %d %d    ILPFCB,HDRY,NPLPF,IPHDRY\n' % (flow_filunit, hdry, nplpf, iphdry))
    _out = np.squeeze(laytyp) # for 1D arrays
    fmt12 = fmt1 + ' LAYTYP'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(layave) # for 1D arrays
    fmt12 = fmt1 + ' LAYAVE'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(chani) # for 1D arrays
    fmt12 = fmt1 + ' CHANI'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(layvka) # for 1D arrays
    fmt12 = fmt1 + ' LAYVKA'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    
    _out = np.squeeze(laywet) # for 1D arrays
    fmt12 = fmt1 + ' LAYWET'
    np.savetxt(fobj, _out[None], delimiter=' ', fmt=fmt12)    

    # comment out to turn lpf to upw
#    if any(laywet == 1):
#        fobj.write('%g %d %d       WETFCT, IWETIT, IHDWET\n' % (WETFCT, IWETIT, IHDWET))
    
    # -- Write HKSAT and Ss, Sy (if Tr) in .lpf file
    # loop thru layers (different entry for each layer)
    for ilay in range(NLAY):
        fobj.write('INTERNAL   1.000E-00 (FREE)    0            HY layer  %d\n' % (ilay+1))
        np.savetxt(fobj, hydcond[:,:,ilay], delimiter=' ', fmt=' %4.2e')      

        fobj.write('INTERNAL   1.000E-00 (FREE)    0            VKA layer  %d\n' % (ilay+1))
        np.savetxt(fobj, VKA[:,:,ilay], delimiter=' ', fmt=' %4.2e')      
    
        if fl_Tr:
            fobj.write('INTERNAL   1.000E-00 (FREE)    0            Ss layer  %d\n' % (ilay+1))
            np.savetxt(fobj, Ss[:,:,ilay], delimiter=' ', fmt=' %4.2e')      
            if laytyp[ilay] > 0: # convertible, i.e. unconfined
                fobj.write('INTERNAL   1.000E-00 (FREE)    0            Sy layer  %d\n' % (ilay+1))
                np.savetxt(fobj, Sy[:,:,ilay], delimiter=' ', fmt=' %4.2e')      
                # Editing Leila's file: laywet should always be 0 for UPW
#                if laywet[ilay] > 0:
#                    fobj.write('INTERNAL   1.000E-00 (FREE)    0            WETDRY layer  %d\n' % (ilay+1))
#                    np.savetxt(fobj, WETDRY[:,:,ilay], delimiter=' ', fmt=' %4.2f')      
    
    fobj.write('\n')
    fobj.close()
    
    return upw_fil_0
    
#    figure
#    for ilay = 1:NLAY
#        subplot(2,2,double(ilay))
#        X = hydcond(:,:,ilay);
#        m = X(X>0); m = min(m(:));
#        imagesc(X), #caxis([m*0.9, max(X(:))]), 
#        cm = colormap;
#    #         cm(1,:) = [1 1 1];
#        caxis([0 max(X(:))])
#        colormap(cm);
#        colorbar
#        title(['hydcond', num2str(ilay)]);
#    end

#%%

# based on write_OC_PCG_MOD_f.m
def write_OC_PCG_MOD_f(GSFLOW_indir, infile_pre, perlen_tr):

#    # =========== TO RUN AS SCRIPT ===========================================
#    # - directories
#    # MODFLOW input files
#    GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#    # MODFLOW output files
#    GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'    
#    perlen_tr = 365*30 + np.ceil(365*30/4)
#    slashstr = '/'
#    # ========================================================================

    # - write to this file
    fil_pcg = infile_pre + '.pcg'
    fil_oc = infile_pre + '.oc'
    
    # -- shoud match .dis
    NPER = 2 # 1 SS then 1 transient
    PERLEN = [1, int(perlen_tr)]  # 2 periods: 1-day steady-state and multi-day transient
    NSTP = PERLEN
    
    # -- pcg and oc files are not changed with this script
    # fil_pcg_0 = fullfile(MODtest_dir0, fil_pcg);
    fil_pcg_0 = GSFLOW_indir + slashstr + fil_pcg
    fobj = open(fil_pcg_0, 'w+')
    # fobj.write('# Preconditioned conjugate-gradient package\n');
    # fobj.write('        50        30         1      MXITER, ITER1, NPCOND\n')
    # fobj.write('  0000.001      .001        1.         2         1         1      1.00\n')
    # fobj.write(' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n')
    
    # # sagehen example:
    # fobj.write('# Preconditioned conjugate-gradient package\n')
    # fobj.write('        1000    450         1      MXITER, ITER1, NPCOND\n')
    # fobj.write('      0.001     0.08       1.0        2         1         0      -0.05  0.70\n')
    # fobj.write(' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n')
    
    # sagehen example:
    fobj.write('# Preconditioned conjugate-gradient package\n');
    fobj.write('        1000    450         1      MXITER, ITER1, NPCOND\n')
    fobj.write('      0.001     0.08       1.0        2         1         0      -0.05  0.70\n')
    fobj.write(' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n')
    
    fobj.close();
    
    # fil_oc_0 = (MODtest_dir0, fil_oc);
    # "PRINT": to listing file
    # "SAVE": to file with unit number in name file
    fil_oc_0 = GSFLOW_indir + slashstr + fil_oc
    fobj = open(fil_oc_0, 'w+')
    fobj.write('HEAD PRINT FORMAT 20\n')
    fobj.write('HEAD SAVE UNIT 51\n')
    fobj.write('COMPACT BUDGET AUX\n')
    fobj.write('IBOUND SAVE UNIT 52\n')
    for per_i in range(NPER):
        for stp_i in range(0,int(NSTP[per_i]),30): # print every 30 days
            fobj.write('PERIOD %d STEP %d\n' % (per_i+1, stp_i+1))
            if stp_i+1 == NSTP[per_i]: # only at end of stress period
                fobj.write('   PRINT HEAD\n')
                fobj.write('   SAVE IBOUND\n')
                fobj.write('   PRINT BUDGET\n')
            fobj.write('   SAVE HEAD\n')
            fobj.write('   SAVE BUDGET\n')
    fobj.close()

#%%

# based on make_sfr2_f_Mannings
#
def make_sfr2_f_Mannings(GSFLOW_indir, infile_pre, reach_fil, dis_fil, segment_fil_all):

# Note: assume .dis file already created!! (reads in TOP for setting STRTOP)

#     # ======== TO RUN AS SCRIPT ===============================================
#     GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#     infile_pre = 'test2lay_py'
#     
#     # for sfr
#     GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/'
#     reach_fil = GIS_indir + 'reach_data.txt'
#     segment_fil_all = []
#     segment_fil_all.append(GIS_indir + 'segment_data_4A_INFORMATION_Man.txt')
#     segment_fil_all.append(GIS_indir + 'segment_data_4B_UPSTREAM_Man.txt')
#     segment_fil_all.append(GIS_indir + 'segment_data_4C_DOWNSTREAM_Man.txt')     
#     # =========================================================================
    
    ##
    
    sfr_file = infile_pre + '.sfr'
    
    # -- Refer to GSFLOW manual p.202, SFR1 manual, and SFR2 manual
    
    # - Refer to Fig. 1 of SFR1 documentation for segment vs. reach numbering
    
    # You need the following inputs (with corresponding structures)
    
    # the followings are used to write item 1
    fl_nstrm = -1  # flag for stream reaches, <0: include unsaturated zone below (sagehen: >0, but uses "REACHINPUT" keyword, which is same effect as nstrm<0)
    nsfrpar = 0  #Always Zero
    nparseg = 0    #Always Zero
    const = 86400.  #Conversion factor used in calculating depth for a stream reach (86400 in sagehen example)
    dleak = 0.0001  #Tolerance level of stream depth used in computing leakage between each stream (0.0001 in sagehen example)
    istcb1 = -1    #Flag for writing stream-aquifer leakage values (>0: file unit, <0: write to listing file)
    istcb2 = 0     #Flag for writing to a seperate formatted file information on inflows&outflows
    isfropt = 3    #defines input structure; saturated or non-saturated zone (1: No UZ; 3: UZ, unsat prop at start of simulation), sagehen uses 3
    nstrail = 10    #Number of trailing-waive increments, incr for better mass balance (10-20 rec'd, sagehen uses 8)    
    isuzn = 1   #Maximum number of vertical cells used to define the unsaturated zone beneath a stream reach (for icalc=1 (Mannings for depth): use isuzn=1)
    nsfrsets = 40  #Maximum number of different sets of trailing waves used to allocate arrays.
    irtflg = 0     #Flag whether transient streamflow routing is active (using kinematic wave approx to St Venant's eq)
    
    project_name = Settings.PROJ_CODE                                             # used to name the output file (.sfr)
    
    # data_indir = '/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/MODFLOW_scripts/sfr_final/data/';
    # data_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
    
    # items 
    reach_data_all = pd.read_csv(reach_fil)       # used to write item 2: assumes 
    
    NPER = 2       # used for item 3
    
    # items 4a: # NSEG ICALC  OUTSEG  IUPSEG  IPRIOR  NSTRPTS  FLOW  RUNOFF  ETSW  PPTSW  ROUGHCH  ROUGHBK  CDPTH  FDPTH  AWDTH  BWDTH
    segment_data_4A = pd.read_csv(segment_fil_all[0]);   # used to write items 4a
    
    segment_data_4B = pd.read_csv(segment_fil_all[1]);      # used to write items 4b (ignored for ICALC=3 in 4a)
    segment_data_4C = pd.read_csv(segment_fil_all[2]);    # used to write items 4c (ignored for ICALC=3 in 4a)
    
    
    # -------------------------------------------------------------------------
    # In case the input text files (e.g. reach_data.txt) contain header lines (comments)
    
    nstrm = reach_data_all.shape[0]
    if fl_nstrm < 0: 
        nstrm = -nstrm

    # sort rows according to increasing segment numbers
    reach_data_all = reach_data_all.sort_values(by='ISEG', ascending=1)
    #nss = max(reach_data_all['ISEG'])
    nss = segment_data_4A.shape[0]


    # sort rows according to increasing reach numbers
    colhead = reach_data_all.columns.get_values()
    ind_IREACH = np.array(np.where(colhead == 'IREACH'))
    for ii in range(nss):
        ind1 = (np.array(np.where(reach_data_all['ISEG'] == ii+1)))
        if ind1.size > 1:
            ind1 = np.squeeze(ind1)
            ind2 = np.array(reach_data_all['IREACH'])[ind1].argsort()
#        a = np.array(reach_data_all.iloc[ind1[ind2]])
#        reach_data_all.iloc[ind1] = a
            reach_data_all.iloc[ind1] = np.array(reach_data_all.iloc[ind1[ind2]])
        else:
            ind1 = np.ones((1,))* ind1[0]

        # renumber IREACH to start at 1 for each segment
        reach_data_all['IREACH'].iloc[ind1] = range(1,len(ind1)+1)
        reach_data_all.iloc[ind1,ind_IREACH[0]] = range(1,len(ind1)+1)
            
    # -- make sure STRTOP is within 1st layer 
    # - read in TOP and BOTM from .dis file
#    dis_file = GSFLOW_indir + slashstr + infile_pre + '.dis'
#    dis_file = '/home/gcng/Shortcuts/AndesWaterResources/GSFLOW/inputs/MODFLOW/test2lay.dis'
    f = open(dis_fil, 'r')
    for i in range(3): # first 2 lines are comments
        line = f.readline().rstrip()
    line = line.split()
    NLAY = int(line[0]) 
    NROW = int(line[1]) 
    NCOL = int(line[2]) 
    NPER = int(line[3]) 
    ITMUNI = int(line[4]) 
    LENUNI = int(line[5])    

    line = f.readline().rstrip()
    LAYCBD = np.array(line.split(), float)
    
    line = f.readline().rstrip()
    line = line.split()
    DELR = float(line[1])
    line = f.readline().rstrip()
    line = line.split()
    DELC = float(line[1])
    f.close()    
    
    TOP = np.genfromtxt(dis_fil, skip_header=7, max_rows=NROW, dtype=float)
    BOTM = np.zeros((NROW, NCOL, NLAY),float);
    for ii in range(NLAY):
        BOTM[:,:,ii] = np.genfromtxt(dis_fil, skip_header=7+(ii+1)*(NROW+1), \
        max_rows=NROW, dtype=float)


    # TOP for cells corresponding to reaches
    TOP_RCH = []
    for ii in range(abs(nstrm)): 
        ind_i = int(reach_data_all.loc[reach_data_all.index[ii],'IRCH'])
        ind_j = int(reach_data_all.loc[reach_data_all.index[ii],'JRCH'])
        TOP_RCH.append(TOP[ind_i-1,ind_j-1])

    # BOTM for cells corresponding to reaches
    BOTM_RCH = []
    for ii in range(abs(nstrm)): 
        ind_i = int(reach_data_all.loc[reach_data_all.index[ii],'IRCH'])
        ind_j = int(reach_data_all.loc[reach_data_all.index[ii],'JRCH'])
        ind_k = int(reach_data_all.loc[reach_data_all.index[ii],'KRCH'])
        BOTM_RCH.append(BOTM[ind_i-1, ind_j-1, ind_k-1])
    
    # - change STRTOP to be just below TOP
    print 'Setting STRTOP (stream top) to be just below (2m) corresponding grid TOP'
    STRTOP = np.array(TOP_RCH) - 2 # 2 m below TOP
    if any(STRTOP-reach_data_all['STRTHICK'] < BOTM_RCH):
        print 'Error! STRTOP is below BOTM of the corresponding layer! Exiting...'
        quit()
    reach_data_all.loc[:,'STRTOP'] = np.array(STRTOP) 

#    # -- plot stream reaches
#    RCH_mask = np.copy(TOP)
#    for ii in range(abs(nstrm)): 
#        RCH_mask[IRCH[ii]-1,JRCH[ii]-1] = np.amax(TOP)*2
#    figure
#    subplot(2,2,1)
#    imagesc(TOP), colorbar, 
#    cm = colormap; 
#    cm(end,:) = [1 1 1];
#    caxis([min(TOP(:)) max(TOP(:))* 1.25]);
#    colormap(cm);
#    subplot(2,2,2)
#    imagesc(RCH_mask), colorbar, 
#    cm = colormap; 
#    cm(end,:) = [1 1 1];
#    caxis([min(TOP(:)) max(TOP(:))* 1.25]);
#    colormap(cm);
    
#    RCH_NUM = zeros(NROW,NCOL); SEG_NUM = zeros(NROW,NCOL);
#    for ii = 1:abs(nstrm), RCH_NUM(IRCH(ii),JRCH(ii)) = IREACH(ii); end
#    for ii = 1:abs(nstrm), SEG_NUM(IRCH(ii),JRCH(ii)) = ISEG(ii); end
#    figure
#    imagesc(RCH_NUM), colorbar, 
#    colormap(jet(1+max(IREACH)));
#    caxis([-0.5 max(IREACH)+0.5])
#    figure
#    imagesc(SEG_NUM), colorbar, 
#    colormap(jet(1+max(ISEG)));
#    caxis([-0.5 max(ISEG)+0.5])

#     # when running as script: to visualize segments one at a time
#     for j = 1: max(ISEG)
#         RCH_NUM = zeros(NROW,NCOL); SEG_NUM = zeros(NROW,NCOL);
#         ind = find(ISEG==j);
#         for ii = ind(:)' 
#             RCH_NUM(IRCH(ii),JRCH(ii)) = IREACH(ii); 
#             SEG_NUM(IRCH(ii),JRCH(ii)) = ISEG(ii); 
#         end                
#         
#         figure(100)
#         imagesc(RCH_NUM), colorbar, 
#         colormap(jet(1+max(RCH_NUM(:))));
#         caxis([-0.5 max(RCH_NUM(:))+0.5])
#         title(['reaches for seg ', num2str(j)]);
#         figure(101)
#         imagesc(SEG_NUM), colorbar, 
#         colormap(jet(1+max(SEG_NUM(:))));
#         caxis([-0.5 max(SEG_NUM(:))+0.5])
#         title(['seg ', num2str(j)]);
#         pause
#     end
        
    # -- threshold slope at minimum 0.001
    ind = np.where(reach_data_all['SLOPE'] < 0.001)
    reach_data_all.loc[reach_data_all.index[ind],'SLOPE'] = 0.001

#    print 'setting various streambed properties (overwriting values in ' + reach_fil + ')'
#
#    # -- set streambed thickness (Sagehen uses constant 1m)
#    reach_data_all.loc[:,'STRTHICK'] = 1 # [m]
#
#    # -- set streambed hydraulic conductivity (Sagehen example: 5 m/d)
#    reach_data_all.loc[:,'STRHC1'] = 5 # [m]
#    
#    # set streambed theta_s
#    reach_data_all.loc[:,'THTS'] = 0.35
#    
#    # set streambed initial theta
#    reach_data_all.loc[:,'THTI'] = 0.3
#    
#    # set streambed Brooks-Corey exp (sagehen example is 3.5)
#    reach_data_all.loc[:,'EPS'] = 3.5
#    
#    # set streambed unsaturated zone saturated hydraulic conductivity
#    # (sagehen example is 0.3 m/d)
#    reach_data_all.loc[:,'UHC'] = 0.3    
#    
    nss = segment_data_4A.shape[0]
    # if isstruct(stress_periods)
    #     stress_periods = stress_periods.data;
    # end
    # - specify only for 2 stress periods:
    stress_periods = np.zeros((NPER, 3)) # itmp, irdflg, iptflg (latter 2 are set to 0)
    stress_periods[0,0] = nss
    if NPER > 1: 
        stress_periods[1:,0] = -1
    
    # -------------------------------------------------------------------------
    
    # First put 4A, 4B and 4C data all together in a cell array
    # size(cell) = nitems x 1 x nperiods
    # In this case, nitems is 3 (i.e. 4A, 4B and 4C)
    
    nitems = 3
    nperiods = stress_periods.shape[0]
    
    # -------------------------------------------------------------------------
    # validate some of the input data
    if nstrm < 0:
        if not any(isfropt == np.array([1, 2, 3, 4, 5])):
            print 'Error: ISFROPT should be set to an integer of 1, 2, 3, 4 or 5.'
            quit()

    if nsfrpar != 0:
        print 'Error: nsfrpar must be 0 because parameters not supported in GSFLOW'
        quit()

    if nparseg != 0:
        print 'Error: nparseg must be 0 because parameters not supported in GSFLOW'
        quit()    
    
    # -------------------------------------------------------------------------
    
    # Ouput file
    sfr_file_0 = GSFLOW_indir + slashstr + sfr_file
    fobj = open(sfr_file_0, 'w+')    
    
    # Write header lines (item 0)
    heading = '# Streamflow-Routing (SFR7) input file.\n'
    fobj.write(heading);
    fobj.write('# %s simulation \n' % (project_name));
    
    
    # Item 1
    fobj.write('  %5d  %5d  %5d  %5d  %8.2f  %8.4f  %5d  %5d' 
    % (nstrm, nss, nsfrpar, nparseg, const, dleak, istcb1, istcb2))
    
    if isfropt >= 1:
        fobj.write('  %5d' % (isfropt))
        if isfropt == 1:
            fobj.write('  %5d\n' % (irtflg))
        elif isfropt > 1:
            fobj.write('  %5d  %5d  %5d  %5d\n' % (nstrail, isuzn, nsfrsets, irtflg))
    else:
        fobj.write('\n');    
    
    # Item 2
    if isfropt == 1:
        ncols_reach = 10
    elif isfropt == 2:
        ncols_reach = 13
    elif isfropt == 3:
        ncols_reach = 14
    else:
        ncols_reach = 6
    reach_data_copy = reach_data_all.iloc[:, 0:ncols_reach]
    
    p = ncols_reach - 5
    fmt_reach = '  %5d  %5d  %5d  %5d  %5d'
    for p_i in range(p):
        fmt_reach = fmt_reach + '  %8.3f'
    np.savetxt(fobj, np.array(reach_data_copy), fmt=fmt_reach)
        
    # Item 3 and 4
    nper = stress_periods.shape[0]
    
    for iper in range(nper):
        
        # write item 3 to the file
        np.savetxt(fobj, stress_periods[iper,:][np.newaxis], fmt='  %5d')
        
        itmp = stress_periods[iper,0]        
        if itmp > 0:
            # segment_data_4A
            
            for iitmp in range(int(itmp)):   # start loop over itmp (num_segments)
                dummy4a = segment_data_4A.iloc[iitmp,:]
                # write item 4a to the file

                fobj.write('    %5d  %5d  %5d  %5d' % (dummy4a['NSEG'], dummy4a['ICALC'], dummy4a['OUTSEG'], dummy4a['IUPSEG']))
                
                if dummy4a['IUPSEG'] > 0:
                    fobj.write('  %5d' % (dummy4a['IPRIOR']))
                
                if dummy4a['ICALC'] == 4:
                    fobj.write('  %5d' % (dummy4a['NSTRPTS']))
                
                fobj.write('  %8.3f  %8.3f  %8.3f  %8.3f' % (dummy4a['FLOW'], dummy4a['RUNOFF'], dummy4a['ETSW'], dummy4a['PPTSW']))
                
                if (dummy4a['ICALC'] == 1) or (dummy4a['ICALC'] == 2):
                    fobj.write('  %8.3f' % (dummy4a['ROUGHCH']))
               
                if dummy4a['ICALC'] == 2:
                    fobj.write('  %8.3f' % (dummy4a['ROUGHBK']))
                
                if dummy4a['ICALC'] == 3:
                    fobj.write('  %8.3f  %8.3f  %8.3f  %8.3f' % (dummy4a['CDPTH'], dummy4a['FDPTH'], dummy4a['AWDTH'], dummy4a['BWDTH']))
                
                fobj.write('\n')
                
                # write items 4b and 4c to the file
                for i in range(2):   # start loop through 4a and 4b
                    if i == 0: 
                        dummy4bc = segment_data_4B.iloc[iitmp,:]
                    else:
                        dummy4bc = segment_data_4C.iloc[iitmp,:]                     
                    
                    fl_no_4bc = 0
                    if (any(isfropt == np.array([0, 4, 5])) and (dummy4a['ICALC'] <= 0)):
                        fmt = '      %8.3f  %8.3f  %8.3f  %8.3f  %8.3f'
                        fobj.write(fmt % (dummy4bc['HCOND'+str(i+1)], dummy4bc['THICKM'+str(i+1)], 
                                                   dummy4bc['ELEVUPDN'+str(i+1)], dummy4bc['WIDTH'+str(i+1)], 
                                                            dummy4bc['DEPTH'+str(i+1)]))
                    elif (any(isfropt == np.array([0, 4, 5])) and (dummy4a['ICALC'] == 1)):
                        if i == 0:
                            fobj.write('    %8.3f' % (dummy4bc['HCOND'+str(i+1)]))
                        
                        if (iper+1 == 1):   # only for the first period
                            fmt = '  %8.3f  %8.3f  %8.3f'
                            fobj.write(fmt % (dummy4bc['THICKM'+str(i+1)], 
                                                       dummy4bc['ELEVUPDN'+str(i+1)], dummy4bc['WIDTH'+str(i+1)]))
                            
                            if ((isfropt == 4) or (isfropt == 5)):
                                fmt = '  %8.3f  %8.3f  %8.3f'
                                fobj.write(fmt % (dummy4bc['THTS'+str(i+1)], 
                                                           dummy4bc['THTI'+str(i+1)], dummy4bc['EPS'+str(i+1)]))
                            
                            if (isfropt == 5):
                                fobj.write('  %8.3f' % (dummy4bc['UHC'+str(i+1)]))
                            
                        elif ((iper+1 > 1) and (isfropt == 0)):
                            fmt = '  %8.3f  %8.3f  %8.3f'
                            fobj.write(fmt % (dummy4bc['THICKM'+str(i+1)], 
                                                       dummy4bc['ELEVUPDN'+str(i+1)], dummy4bc['WIDTH'+str(i+1)]))
                        
                    elif (any(isfropt == np.array([0, 4, 5])) and (dummy4a['ICALC'] >= 2)):
                        fobj.write('    %8.3f' % (dummy4bc['HCOND'+str(i+1)]))
                        
                        if not (any(isfropt == np.array([4, 5])) and (iper+1 > 1) 
                        and (dummy4a['ICALC'] == 2)):
                            fobj.write('  %8.3f  %8.3f' % (dummy4bc['THICKM'+str(i+1)], dummy4bc['ELEVUPDN'+str(i+1)]))
                            
                        if (any(isfropt == np.array([4, 5])) and (iper+1 == 1) 
                        and (dummy4a['ICALC'] == 2)):
                                fmt = '  %8.3f  %8.3f  %8.3f'
                                fobj.write(fmt % (dummy4bc['THTS'+str(i+1)], 
                                                           dummy4bc['THTI'+str(i+1)], dummy4bc['EPS'+str(i+1)]))
                                
                                if (isfropt == 5):
                                    fobj.write('  %8.3f' % (dummy4bc['UHC'+str(i+1)]))                        
                    elif ((isfropt == 1) and (dummy4a['ICALC'] <= 1)):
                        fobj.write('    %8.3f' % (dummy4bc['WIDTH'+str(i+1)]))
                        
                        if (icalc <= 0):
                            fobj.write('  %8.3f' % (dummy4bc['DEPTH'+str(i+1)]))
                        
                    elif (any(isfropt == np.array([2, 3])) and (dummy4a['ICALC'] <= 1)):
                        if (iper+1 == 1):
                            fobj.write('    %8.3f' % (dummy4bc['WIDTH'+str(i+1)]))
                            
                            if (dummy4a['ICALC'] <= 0):
                                fobj.write('  %8.3f' % (dummy4bc['DEPTH'+str(i+1)]))
                    else:
                        fl_no_4bc = 1
                    if fl_no_4bc == 0:
                        fobj.write('\n')
   
    fobj.close()
    
    return sfr_file_0
    
#%%

def MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment):
    if data_type == 'INT':
        CNSTNT0 = str(int(CNSTNT))
        fmt0 = '%7d'
    elif data_type == 'REAL':
        CNSTNT0 = str(float(CNSTNT))
        fmt0 = '%7e'
    if np.array(data).size == 1:
        str0 = 'CONSTANT     ' + fmt0 + ' %s \n'
        fobj.write(str0 % (data, comment))
    else:
        fobj.write('INTERNAL  %10s%20s%10s %s \n' % (CNSTNT0, '(FREE)', '-1', comment))
        np.savetxt(fobj, data, delimiter=' ', fmt=fmt0)

def make_uzf3_f_2(GSFLOW_indir, infile_pre, surfz_fil, dischargept_fil, ba6_fil):

#    print 'UZF: Had to play around alot with finf (infiltration) to get convergence!!'


#     ######### TO RUN AS SCRIPT ##############################################
#     # - directories
#     # MODFLOW input files
#     GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
#     # MODFLOW output files
#     GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
#         
#     infile_pre = 'test2lay_py'
#     
#     NLAY = 2
#     DZ = [100, 50] # [NLAYx1] [m]
#     
#     # length of transient stress period (follows 1-day steady-state period) [d]
#     # perlen_tr = 365 # ok if too long
#     perlen_tr = 365*30 + np.ceil(365*30/4)
#     
#     GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/'
#     
#     # use restart file as initial cond (empty string to not use restart file)
#     fil_res_in = '' # empty string to not use restart file
##     fil_res_in = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/spinup30yr_constH/outputs/MODFLOW/test2lay.out' # empty string to not use restart file
#     
#     # for various files: ba6, dis, uzf, lpf
#     surfz_fil = GIS_indir + 'topo.asc'
#     # for various files: ba6, uzf
#     mask_fil = GIS_indir + 'basinmask_dischargept.asc'
#    # ########################################################################

    # -------------------------------------------------------------------------
    # You need the following inputs
    
    # - write to this file
    # GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
    uz_file = infile_pre + '.uzf'
    
    # surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
    sdata = read_grid_file_header(surfz_fil)
        
    NSEW = [sdata['north'], sdata['south'], sdata['east'], sdata['west']]
    NROW = sdata['rows'] 
    NCOL = sdata['cols']
    
    # - set TOP to surface elevation [m]
    TOP = np.genfromtxt(surfz_fil, skip_header=6, delimiter=' ', dtype=float)
    
    IBOUND = np.genfromtxt(ba6_fil, skip_header=3, max_rows=NROW, dtype=float)

    
    NPER = 2
    # **** ASSUMES PER 1 IS SS, PER 2 IS TR ****
    
    #Item 1:
    #NUZTOP: Which cell (layer) has recharge&discharge simulated in it; 
    #   1:top layer has recharge/discharge and land-surface (at top)
    #   2:layer specificed in UZFBND, top of that layer is land-surface
    #   3:highest active layer, top of UZFBND is land-surface (*sagehen example, with watershed mask of UZFBND=1)
    NUZTOP = 3     
    IUZFOPT = 1    # 1: Vertical hydraulic conductivity is specified using VKS, 2: use VKA from LPF file
    IRUNFLG = 0    # >0: Groundwater discharged to land-surface is routed to streams or lakes (using IRUNBND) (*sagehen=0, ONLY avaliable in GSFLOW for SS)
    IETFLG = 0     # ~=0: Evaporation will be simulated. (*sagehen=0; GSFLOW: =0 for PRMS ET only, =1 for additional ET from below soil zone)
    IUZFCB1 = 61    #Flag for writing rates of groundwater recharge, ET&groundwater discharge in UBUDSV format.0:wont be written, >0: file unit number
    IUZFCB2 = 0   #Writing groundwater recharge, discharge&ET in UBUDSV3 format; >0: file unit number
#    NTRAIL2 = 25   #Number of trailing waves to define theta profile, 10 to 20 usually ok, higher for better mass balance
    NTRAIL2 = 15   #Number of trailing waves to define theta profile, 10 to 20 usually ok, higher for better mass balance
    NSETS2 = 100    #Number of wave sets to simulate multiple infiltration periods, 20 usually ok.
    NUZGAG = 0     # number of cells for which to print detailed info on UZ water budget and theta (see uzgag below)
    SURFDEP = 1.0  # average undulation depth within finite diff cell (?)
    
    #Item 2-7:
    project_name = Settings.PROJ_CODE   # used to name the output file (.uzf)
    iuzfbnd = np.copy(IBOUND) # [NROW,NCOL] layer w/ top as land-surface and/or w/ discharge/recharge (see NUZTOP), default: mask with 1's
    iuzfbnd[iuzfbnd<0] = 0
    iuzfbnd[iuzfbnd>0] = 1
    if IRUNFLG > 0:
        print 'Error!  Input scripts only set up for IRUNFLG = 0!'
        quit()        
#        irunbnd = importdata('./data/irunbnd.dat'); # [NROW,NCOL] only for IRUNFLG>0, stream seg to which gw discharge is routed
    # vks = importdata('./data/vks.dat'); # [NROW,NCOL] saturated K, no needed if using value in LPF (IUZFOPT=2), [m/d]
    vks = np.copy(iuzfbnd) * 4.
    # Ok to have following parameters as SCALAR (constant for all gridcells) or as ARRAY (NROWxNCOL)
    eps = 3.5  #Brooks-Corey epsilon of the unsaturated zone.
    thts = np.copy(iuzfbnd) * 0.35    #Saturated water content of the unsaturated zone
#    eps = 4  #Brooks-Corey epsilon of the unsaturated zone.
#    thts = np.copy(iuzfbnd) * 0.18    #Saturated water content of the unsaturated zone
    if NUZGAG > 0:
        print 'Input scripts only set up for UZGAG = 0! Exiting...'
        quit()        
#        uzgag = importdata('./data/uzgag.dat'); # only for NUZGAG>0; row, col, file unit for output, IUZOPT (flag for what data to output)    
    # - infiltration (in general, not needed  bc provided by PRMS, but option to apply for initial SS period for MODFLOW)
    NUZF1 = np.array([[1], 
             -1*np.ones((NPER-1,1))]) # infiltration can ONLY specified for (initial) SS stress periods (o.w. calculated by PRMS)
    # A = importdata('./data/finf.dat'); # infiltration rate [NROW,NCOL,1], **max ONLY for 1 stress period: first SS stress period! [m/d]
    # finf = A';
    # B = reshape(A', NCOL, NROW, NPER);
    # for i=1:size(B, 3)
    #     finf(:, :, i) = transpose(B(:, :, i));
    # end
    # finf = ones(NROW,NCOL);
    # finf(:,:,1) = ones(NROW,NCOL)*8.8e-5; # m/d (8.8e-4 m/d typical in Sagehen)
    
    # - set infiltration: FINF (for S.S. initialization) according to TOP (high elev only)
    f = open(dischargept_fil, 'r')
    last_line = f.readlines()
    last_line = last_line[-1].rstrip()
    f.close()    
    value1, value2 = last_line.split(': ')
    value2 = value2.split(' ')
    dischargePt_rowi = int(value2[1])
    dischargePt_coli = int(value2[3])
 
    #value1, value2 = last_line.split(': ')
    #value2 = value2.split(' ')
    #dischargePt_rowi = int(value2[1])
    #dischargePt_coli = int(value2[3])
    TOP_mask = TOP 
    TOP_mask[IBOUND==0] = 0
    a = np.where(TOP_mask.flatten()>0)
    NZ = np.array(a).size
    z_sort = np.sort(TOP_mask.flatten())[::-1] # descending order
    
    # this works: (for no sfr: converges, but mostly all dry)
    ind = TOP_mask > z_sort[int(round(NZ/10))]
#    finf = np.zeros((NROW,NCOL))
#    finf[ind] = 4e-4; # m/d (8.8e-4 m/d typical in Sagehen)
#    
#    # testing: (for no sfr: kind of works, S.S. does not converge but transient mostly
#    # does, starts wet then dries out)
#    finf = np.zeros((NROW,NCOL))
#    ind = TOP_mask > z_sort[int(round(0.25*NZ))]
#    finf[IBOUND!=0] = 4e-4 / 10 # m/d (8.8e-4 m/d typical in Sagehen)
#    finf[ind] = 4e-4; # m/d (8.8e-4 m/d typical in Sagehen)
#    finf = np.zeros((NROW,NCOL)) + 0.001 # from Lauren's test 10/1/17
    
#    # scale FINF by elev
#    # (# Sagehen max 3*8.e-4, min 8.e-4
#    maxFINF = 3* 8.e-4 
#    minFINF = 8e-4
#    midZ_factor = (2150-1928)/(2649-1928) # based on Sagehen, where FINF levels off; sagehen: 0.31
#    midZ = midZ_factor * (np.max(TOP_mask.flatten())-np.min(TOP_mask[TOP_mask>0])) + np.min(TOP_mask[TOP_mask>0])
#    m = (maxFINF-minFINF)/(np.max(TOP_mask.flatten())-np.min(TOP_mask[TOP_mask>0]))
#    finf = m*(TOP_mask-np.min(TOP_mask[TOP_mask>0])) + minFINF
#    finf[TOP_mask<=midZ] = 4e-4
#    finf[IBOUND==0] = 0
#
#    finf[:] = 8.8e-6    
#    finf[:] = 8.8e-4    

#    finf = np.ones((NROW,NCOL))
    finf = np.copy(iuzfbnd)
    
    try:
       float(Settings.finf0)
       finf = float(Settings.finf0) * finf
    except ValueError:
       for ii in range(NLAY):
            finf[:,:,ii] = np.genfromtxt(Settings.finf0, skip_header=1+ii*(NROW+1), \
            max_rows=NROW, dtype=float)
    
    
    # # testing: 
    # finf = np.zeros((NROW,NCOL))
    # ind = TOP_mask > z_sort[int(round(0.25*NZ))]
    # finf[IBOUND!=0] = 4e-2 / 10 # m/d (8.8e-4 m/d typical in Sagehen)
    # finf[ind] = 4e-2 # m/d (8.8e-4 m/d typical in Sagehen)
    # finf[:] = 0;
    # finf[IBOUND!=0] = 4e-4 / 100 # m/d (8.8e-4 m/d typical in Sagehen)
    # finf[ind] = 4e-5 # m/d (8.8e-4 m/d typical in Sagehen)
    
    # finf[ind] = 4e-4; # m/d (8.8e-4 m/d typical in Sagehen) - just this works
#    figure, subplot(2,2,1), 
#    X = finf;
#    m = X(X>0); m = min(m(:));
#    imagesc(X(:,:)), 
#    caxis([m*0.9, max(X(:))]), 
#    cm = colormap;
#    cm(1,:) = [1 1 1];
#    colormap(cm);
#    colorbar
#    title('FINF (S.S. infil) [m/d]');
    
    
    # - ET (in general, not needed bc provided by PRMS, but option to apply excess ET below PRMS' soil-zone, i.e. apply to MODFLOW's UZ)
    # pet = 5.0E-08;  #array of ET demands rate (L/T), PET not used in GSFLOW
    # for IETFLG>0, specify how ET can be drawn from below soil-zone base;
    #   assume same for all stress periods
    if IETFLG>0:
        NUZF2 = -1*ones((NPER,1)) # use ET from below soil-zone
        NUZF3 = np.array([[1], -1*np.ones((NPER-1,1))]) # only specify extdp for first stress periods
        extdp = 15.0*np.ones((NROW,NCOL))   #array of ET extiction zone~altitude of the soil-zone base;specified at least for 1st stress period; only for IETFLG>0
        NUZF4 = np.array([[1], -1*np.ones((NPER-1,1))]) # only specify extwc for first stress periods
#        extwc = thts*0.9*np.ones((NROW,NCOL)) #array of Extinction water content; EXTWC must be between (THTS-Sy) and THTS; only for IETFLG>0 and 
        extwc = thts*0.9 #array of Extinction water content; EXTWC must be between (THTS-Sy) and THTS; only for IETFLG>0 and 
    # -------------------------------------------------------------------------
    
    # Ouput file
    uz_file_0 = GSFLOW_indir + slashstr + uz_file
    fobj = open(uz_file_0, 'w+')  
    
    # Write header lines (item 0)
    heading = '# Unsaturated-Zone Flow (UZF) input file.\n';
    fobj.write(heading);
    fobj.write('# %s simulation \n' % (project_name)) 
    
    # Write item 1
    fmtarr_int = ''
    for ii in range(9):
        fmtarr_int = fmtarr_int + '  %6d'
    fmtarr_int = fmtarr_int + '  %10.6E'
    fobj.write(fmtarr_int % (NUZTOP, IUZFOPT, IRUNFLG, IETFLG, IUZFCB1, IUZFCB2, NTRAIL2, NSETS2, NUZGAG, SURFDEP))
    
    comment = '     NUZTOP  IUZFOPT  IRUNFLG  IETFLG  IUZFCB1  IUZFCB2  NTRAIL2  NSETS2  NUZGAG  SURFDEP\n'
    fobj.write(comment)

    # Generally use these settings
    LOCAT = 'INTERNAL'
    CNSTNT = 1
    IPRN = -1
    
    # Write item2 [IUZFBND (NCOL, NROW)] - U2DINT
    # INTERNAL  1.0  (FREE)  -1
    # "INTERNAL" b/c data follows in same file,
    # 1.0 because all values scaled by 1.0
    # -1 is IPRN print flag (<0 means value NOT printed to list file)
    comment = '#IUZFBND--AREAL EXTENT OF THE ACTIVE MODEL'
    data = iuzfbnd
    data_type = 'INT'
    MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)
        
    # write item 3 [IRUNBND (NCOL, NROW)] - U2DINT
    if (IRUNFLG > 0):
        comment = '#IRUNBND--STREAM SEGMENTS OR LAKE NUMBERS'
        data = irunbnd
        data_type = 'INT'
        MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)
    
    # write item 4 [VKS (NCOL, NROW)] - U2DREL
    if (IUZFOPT == 1):
        comment = '#VKS--VERTICAL HYDRAULIC CONDUCTIVITY OF THE UNSATURATED ZONE'
        data = vks
        data_type = 'REAL'
        MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)
    
    # write items 5, 6, 7
    comment = '#EPS--BROOKS/COREY EPSILON'
    data = eps
    data_type = 'REAL'
    MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)
    
    comment = '#THTS--SATURATED WATER CONTENT'
    data = thts
    data_type = 'REAL'
    MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)
    
#    comment = '#THTI--INITIAL WATER CONTENT'
#    data = thti
#    data_type = 'REAL'
#    MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)    
    
    # write item 8
    if (NUZGAG > 0):
        comment = '#IUZROW IUZCOL IFTUNIT IUZOPT--UNSATURATED FLOW OUTPUT';
        for i in range(NUZGAG):
            dummyrow = uzgag[i, :]
            fobj.write('  %6d  %6d  %6d  %6d      %s\n' % (dummyrow, comment))
    
    # write items 9-16
    for iper in range(NPER):
        
        # write item 9
        comment = '#NUZF1 FOR STRESS PERIOD ' + str(int(iper+1)) + '\n'
        fobj.write('  %6d      %s' % (NUZF1[iper], comment))
        
        # write item 10
        if (NUZF1[iper] > 0):
    #         ent = get_file_entry(transpose(finf(:, :, iper)), 'U2DREL', ...
    #             1.0, sprintf('#FINF--STRESS PERIOD #d', iper));
            comment = '#FINF--STRESS PERIOD ' + str(int(iper+1))
            data = finf
            data_type = 'REAL'
            MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)    
        
        # ---------------------------------------------------------------------
        
        # items 11-16
        if (IETFLG > 0):
            
            # write items 11
            comment = '#NUZF2 FOR STRESS PERIOD (GSFLOW DOES NOT USE PET) ' + str(int(iper+1)) + '\n'
            fobj.write('  %6d      %s' % (NUZF2[iper], comment))
            
            # write item 12 (GSFLOW DOES NOT USE PET)
#            if (NUZF2(iper) > 0):
            
            # -----------------------------------------------------------------        
            
            # write items 13
            comment = '#NUZF3 FOR STRESS PERIOD ' + str(int(iper+1)) + '\n'
            fobj.write('  %6d      %s' % (NUZF3[iper], comment))
            
            # write item 14
            if (NUZF3[iper] > 0):
                comment = '#EXTDP FOR STRESS PERIOD ' + str(int(iper+1))
                data = extdp
                data_type = 'REAL'
                MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)    
            
            # -----------------------------------------------------------------
            
            # write items 15
            comment = '#NUZF4 FOR STRESS PERIOD ' + str(int(iper+1)) + '\n'
            fobj.write('  %6d      %s' % (NUZF4[iper], comment))

            
            # write item 16
            if (NUZF4[iper] > 0):
                comment = '#EXTWC FOR STRESS PERIOD ' + str(int(iper+1))
                data = extwc
                data_type = 'REAL'
                MOD_data_write2file(fobj, LOCAT, CNSTNT, IPRN, data_type, data, comment)    
    
    fobj.close()
    
    return uz_file_0
    
#%%
# Based on Leila's script NWT_write.m
# Adapted into function

def NWT_write_file(GSFLOW_indir, infile_pre):

# -------------------------------------------------------------------------
# input variables
# see Table 2., page 12 of MODFLOW-NWT Techniques and Methods

#     ######### TO RUN AS SCRIPT ##############################################
#     # - directories
#     # MODFLOW input files
#     GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/'
##     # MODFLOW output files
##     GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/'
#         
#     infile_pre = 'test2lay_py'
#     #########################################################################


    # Iteration Control
    headtol = 1e-4;
    fluxtol = 500;
    maxiterout = 100;
    
    # Dry Cell Tolerance
    thickfact = 0.00001;
    
    # NWT Options
    linmeth = 1;
    iprnwt = 0;
    ibotav = 0;
    options = 'SPECIFIED';
    
    # Under Relaxation Input
    dbdtheta = 0.7;
    dbdkappa = 0.0001;
    dbdgamma = 0.0;
    momfact = 0.1;
    
    #Residual Control
    backflag = 0;
    maxbackiter = 50;
    backtol = 1.2;
    backreduce = 0.75;
    
    # Linear Solution Control and Options for GMRES
    maxitinner = 50;
    ilumethod = 2;
    levfill = 1;
    stoptol = 1e-10;
    msdr = 10;
    
    # Linear Solution Control and Options for xMD
    iacl = 2;
    norder = 1;
    level = 1;
    north = 2;
    iredsys = 0;
    rrctols = 0.0;
    idroptol = 1;
    epsrn = 1e-3;
    hclosexmd = 1e-4;
    mxiterxmd = 50;
    
    nwt_file_0 = GSFLOW_indir + slashstr + infile_pre + '.nwt'
    
    headings = ['NWT Input File', 'Test Problem 3 for MODFLOW-NWT']
    
    # -------------------------------------------------------------------------
    
    fobj = open(nwt_file_0, 'w+');
    
    # item 0 -------
    for head0 in headings:
        fobj.write('# ' + head0 + '\n')
    
    # item 1 -------
    fobj.write('%10.3e  %10.3e  %4d  %10.3e  %d  %d  %d  ' %
        (headtol, fluxtol, maxiterout, thickfact, linmeth, iprnwt, ibotav))
    
    fobj.write('%s  ' % (options))
    
    if options == 'SPECIFIED':
        fobj.write('%10.3g  %10.3g  %10.3g  %10.3g  %d  ' %
            (dbdtheta, dbdkappa, dbdgamma, momfact, backflag))
        
        if backflag > 0:
            fobj.write('%d  %10.3g  %10.3' % (maxbackiter, backtol, backreduce))
        
        fobj.write('\n');
        
        # item 2a -------
        if linmeth == 1:
            fobj.write('%4d  %d  %d  %10.3g  %2d' %
                (maxitinner, ilumethod, levfill, stoptol, msdr))
        # item 2b -------
        elif linmeth == 2:
            fobj.write('%d  %d  %2d  %2d  %d  %10.3g  %d  %10.3g  %10.3g  %4d' %
                (iacl, norder, level, north, iredsys, rrctols, idroptol, 
                epsrn, hclosexmd, mxiterxmd))
        
    fobj.write('\n')
    
    fobj.close()
    
    return nwt_file_0
    # -------------------------------------------------------------------------
    # End of the script
