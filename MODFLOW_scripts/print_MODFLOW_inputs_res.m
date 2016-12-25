% print_MODFLOW_inputs.m
clear all, close all, fclose all;

addpath('/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/MODFLOW_scripts/uzf/')

% - directories
% MODFLOW input files
GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
% MODFLOW output files
GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/';

% % MODFLOW input files
% GSFLOW_indir = './';
% % MODFLOW output files
% GSFLOW_outdir = '';


% infile_pre = 'test1lay';
% NLAY = 1;
% DZ = 10; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)

infile_pre = 'test2lay_res';
% infile_nam_pre = 'test'; % shorter name bc of PRMS control file string length limit

NLAY = 2;
DZ = [100; 50]; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
% DZ = [350; 100]; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)

% length of transient stress period (follows 1-day steady-state period) [d]
% perlen_tr = 365; % ok if too long
% perlen_tr = 365*5 + ceil(365*5/4); % ok if too long (I think, but maybe run time is longer?)
perlen_tr = 365*30 + ceil(365*30/4); % ok if too long (I think, but maybe run time is longer?)

GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';

% use restart file as initial cond (empty string to not use restart file)
% fil_res_in = ''; % empty string to not use restart file
fil_res_in = '/home/gcng/workspace/Pfil_res_inrojectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/test2lay_melt_30yr.out'; % empty string to not use restart file

% for various files: ba6, dis, uzf, lpf
surfz_fil = [GIS_indir, 'topo.asc'];
% surfz_fil = [GIS_indir, 'SRTM_new_20161208.asc'];
% for various files: ba6, uzf
mask_fil = [GIS_indir, 'basinmask_dischargept.asc'];

% for sfr
reach_fil = [GIS_indir, 'reach_data.txt'];
segment_fil_all = cell(3,1);
% segment_fil_all{1} = [GIS_indir, 'segment_data_4A_INFORMATION.txt'];
% segment_fil_all{2} = [GIS_indir, 'segment_data_4B_UPSTREAM.txt'];
% segment_fil_all{3} = [GIS_indir, 'segment_data_4C_DOWNSTREAM.txt'];
segment_fil_all{1} = [GIS_indir, 'segment_data_4A_INFORMATION_Man.csv'];
segment_fil_all{2} = [GIS_indir, 'segment_data_4B_UPSTREAM_Man.csv'];
segment_fil_all{3} = [GIS_indir, 'segment_data_4C_DOWNSTREAM_Man.csv'];



%% 
write_dis_MOD2_f(GSFLOW_indir, infile_pre, surfz_fil, NLAY, DZ, perlen_tr);
% write_ba6_MOD2(GSFLOW_indir, infile_pre, surfz_fil, mask_fil, NLAY, DZ); % 
write_ba6_MOD3_2(GSFLOW_indir, infile_pre, mask_fil); % list this below write_dis_MOD2_f
write_lpf_MOD2_f2_2(GSFLOW_indir, infile_pre, surfz_fil, NLAY);

make_uzf3_f_2(GSFLOW_indir, infile_pre, surfz_fil, mask_fil);
% make_sfr2_f(GSFLOW_indir, infile_pre, reach_fil, segment_fil_all); % list this below write_dis_MOD2_f
make_sfr2_f_Mannings(GSFLOW_indir, infile_pre, reach_fil, segment_fil_all); % list this below write_dis_MOD2_f


write_OC_PCG_MOD_f(GSFLOW_indir, infile_pre, perlen_tr);

% write_nam_MOD_f(GSFLOW_indir, GSFLOW_outdir, infile_pre);
write_nam_MOD_f2(GSFLOW_indir, GSFLOW_outdir, infile_pre, fil_res_in);
