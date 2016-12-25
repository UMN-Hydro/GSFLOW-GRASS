% write_dis_MOD (for 3D domains)
% 11/17/16
%
% v1 - 11/30/16 start to include GIS data for Chimborazo's Gavilan Machay
%      watershed; topo.asc for surface elevation (fill in bottom elevation
%      based on uniform thickness of single aquifer)
function write_dis_MOD2_f(GSFLOW_indir, infile_pre, surfz_fil, NLAY, DZ)

% % ==== TO RUN AS SCRIPT ===================================================
% % - directories
% % MODFLOW input files
% GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
% % MODFLOW output files
% GSFLOW_outdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/';
% 
% % infile_pre = 'test1lay';
% % NLAY = 1;
% % DZ = 10; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
% 
% infile_pre = 'test2lay';
% NLAY = 2;
% DZ = [50; 50]; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
% 
% GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
% 
% % for various files: ba6, dis, uzf, lpf
% surfz_fil = [GIS_indir, 'topo.asc'];
% % for various files: ba6, uzf
% mask_fil = [GIS_indir, 'basinmask_dischargept.asc'];
% % =========================================================================
% %%

% clear all, close all, fclose all;

% - write to this file
% GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
dis_file = [infile_pre, '.dis'];

% - read in this file for surface elevation (for TOP(NROW,NCOL))
% surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';

% - read in this file for elevation of layer bottoms (for BOTM(NROW,NCOL,NLAY))
% (layer 1 is top layer)
botmz_fil = '';

% - domain dimensions, maybe already in surfz_fil and botm_fil{}?
% NLAY = 1;
% NROW = 1058;
% NCOL = 1996;

% % - domain boundary (UTM zone 17S, outer boundaries)
% north = 9841200;
% south = 9835900;
% east = 751500;
% west = 741500;

% % - space discretization
% DELR = (east-west)/NCOL; % width of column [m]
% DELC = (north-south)/NROW; % height of row [m]
% DZ = 10; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)
% DZ = [5; 5]; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)

% - time discretization
PERLEN = [1; 365];  % 2 periods: 1-day steady-state and multi-day transient

comment1 = '# test file for Gavilan Machay';
comment2 = '# test file';

% - The following will be assumed:
LAYCBD = zeros(NLAY,1); % no confining layer below layer
ITMUNI = 4; % [d]
LENUNI = 2; % [m]
NPER = 2; % 1 SS then 1 transient
NSTP = PERLEN; TSMULT = 1;  % must have daily time step to correspond with PRMS
SsTr_flag = ['ss'; 'tr'];

%% ------------------------------------------------------------------------

% -- Read in data from files
fid = fopen(surfz_fil, 'r');
D = textscan(fid, '%s %f', 6); 
NSEW = D{2}(1:4);
NROW = D{2}(5);
NCOL = D{2}(6);

% - space discretization
DELR = (NSEW(3)-NSEW(4))/NCOL; % width of column [m]
DELC = (NSEW(1)-NSEW(2))/NROW; % height of row [m]

% - set TOP to surface elevation [m]
D = textscan(fid, '%f'); 
fclose(fid);
fprintf('Done reading...\n');
TOP = reshape(D{1}, NCOL, NROW)'; % NROW x NCOL

BOTM = zeros(NROW, NCOL, NLAY);
BOTM(:,:,1) = TOP-DZ(1);
for ilay = 2:NLAY
    BOTM(:,:,ilay) = BOTM(:,:,ilay-1)-DZ(ilay);
end

% -- Discretization file:
fmt1 = [repmat('%4d ', 1, NCOL), '\n']; % 
fmt2 = [repmat('%10g ', 1, NCOL), '\n']; %

fid = fopen([GSFLOW_indir, '/', dis_file], 'wt');
fmt3 = [repmat(' %d', 1, NLAY), '\n']; % for LAYCBD
fprintf(fid, '%s\n', comment1);
fprintf(fid, '%s\n', comment2);
fprintf(fid, ' %d  %d  %d  %d  %d  %d ', NLAY, NROW, NCOL, NPER, ITMUNI, LENUNI);
fprintf(fid, '  NLAY,NROW,NCOL,NPER,ITMUNI,LENUNI\n');
fprintf(fid, fmt3, LAYCBD);
fprintf(fid, 'CONSTANT %14g   DELR\n', DELR);
fprintf(fid, 'CONSTANT %14g   DELC\n', DELC);
fprintf(fid, 'INTERNAL  1.0  (FREE)  0                          TOP ELEVATION OF LAYER 1 \n');
fprintf(fid, fmt2, TOP');
for ii = 1: NLAY
    fprintf(fid, 'INTERNAL  1.0  (FREE)  0                          BOTM ELEVATION OF LAYER %d \n', ii);
    fprintf(fid, fmt2, BOTM(:,:,ii)');
end
for ii = 1: NPER
    fprintf(fid, ' %g %d %g %c%c        PERLEN, NSTP, TSMULT, Ss/Tr (stress period %4d)\n', ...
        PERLEN(ii), NSTP(ii), TSMULT, SsTr_flag(ii,:), ii);
end
fclose(fid);

% -- plot domain discretization
figure
subplot(2,2,1)
imagesc(TOP), 
m = TOP(TOP>0); m = min(m(:));
caxis([m*0.9, max(TOP(:))]), 
cm = colormap;
cm(1,:) = [1 1 1];
colormap(cm);
colorbar
title('TOP')
for ilay = 1:NLAY
    subplot(2,2,1+double(ilay))
    m = BOTM(BOTM>0); m = min(m(:));
    imagesc(BOTM(:,:,ilay)), caxis([m*0.9, max(BOTM(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['BOTM', ' lay', num2str(ilay)]);
end

figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    if ilay == 1
        X = TOP - BOTM(:,:,ilay);
    else
        X = BOTM(:,:,ilay-1)-BOTM(:,:,ilay);
    end
%     m = X(X>0); m = min(m(:));
    imagesc(X), %caxis([m*0.9, max(X(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['DZ', ' lay', num2str(ilay)]);
end


