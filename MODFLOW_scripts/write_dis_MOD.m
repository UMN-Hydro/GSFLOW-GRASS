% write_dis_MOD (for 3D domains)
% 11/17/16

clear all, close all, fclose all;

% - write to this file
GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/';
dis_file = 'test.dis';

% - read in this file for surface elevation (for TOP(NROW,NCOL))
surfz_fil = '';

% - read in this file for elevation of layer bottoms (for BOTM(NROW,NCOL,NLAY))
% (layer 1 is top layer)
botmz_fil = '';

% - domain dimensions, maybe already in surfz_fil and botm_fil{}?
NLAY = 1;
NROW = 50;
NCOL = 50;

% - space discretization
DELR = 30; % width of column [m]
DELC = 30; % height of row [m]

% - time discretization
PERLEN = [1; 365];  % 2 periods: 1-day steady-state and multi-day transient

comment1 = '# test file';
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
% (temp place-holder)
TOP = zeros(NROW,NCOL)+5;
BOTM = zeros(NROW, NCOL, NLAY)+1;

% -- Discretization file:
fmt1 = [repmat('%4d ', 1, NCOL), '\n']; % 
fmt2 = [repmat('%10g ', 1, NCOL), '\n']; %

fid = fopen([GSFLOW_dir, dis_file], 'wt');
fmt3 = [repmat(' %d', 1, NLAY), '\n']; % for LAYCBD
fprintf(fid, '%s\n', comment1);
fprintf(fid, '%s\n', comment2);
fprintf(fid, ' %d  %d  %d  %d  %d  %d ', NLAY, NROW, NCOL, NPER, ITMUNI, LENUNI);
fprintf(fid, '  NLAY,NROW,NCOL,NPER,ITMUNI,LENUNI\n');
fprintf(fid, fmt3, LAYCBD);
fprintf(fid, 'CONSTANT %14g   DELR\n', DELR);
fprintf(fid, 'CONSTANT %14g   DELC\n', DELC);
fprintf(fid, 'INTERNAL  1.0  (FREE)  0                          TOP ELEVATION OF LAYER 1 \n');
for irow = 1: NROW
    fprintf(fid, fmt2, TOP(irow,:));
end
for ii = 1: NLAY
    fprintf(fid, 'INTERNAL  1.0  (FREE)  0                          BOTM ELEVATION OF LAYER %d \n', ii);
    for irow = 1: NROW
        fprintf(fid, fmt2, BOTM(irow,:,ii));
    end
end
for ii = 1: NPER
    fprintf(fid, ' %g %d %g %c%c        PERLEN, NSTP, TSMULT, Ss/Tr (stress period %4d)\n', ...
        PERLEN(ii), NSTP(ii), TSMULT, SsTr_flag(ii,:), ii);
end
fclose(fid);

% -- plot domain discretization
figure(1)
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

figure(2)
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    if ilay == 1
        X = TOP - BOTM(:,:,ilay);
    else
        X = BOTM(:,:,ilay-1)-BOTM(:,:,ilay);
    end
    m = X(X>0); m = min(m(:));
    imagesc(X), caxis([m*0.9, max(X(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['DZ', ' lay', num2str(ilay)]);
end


