% read *.dis
clear all, 
% close all, 
fclose all;

% dis_file = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen/input/modflow/sagehen.dis';

dis_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/test2lay.dis';
GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
mask_fil = [GIS_indir, 'basinmask_dischargept.asc']; % only for Chimb

fid = fopen(dis_file);

for ii = 1:2
    cmt = fgets(fid);
end

line0 = fgets(fid);
D = textscan(line0, '%d', 6);
NLAY = D{1}(1);
NROW = D{1}(2); 
NCOL = D{1}(3);
NPER = D{1}(4);
ITMUNI = D{1}(5);
LENUNI = D{1}(6);

line0 = fgets(fid);
D = textscan(line0, '%d');
LAYCBD = D{1}; % 1xNLAY (0 if no confining layer)

line0 = fgets(fid);
D = textscan(line0, '%s %d');
DELR = D{2}; % width of column
line0 = fgets(fid);
D = textscan(line0, '%s %d');
DELC = D{2}; % height of row

TOP = nan(NROW,NCOL);
line0 = fgets(fid);
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    TOP(irow,:) = D{1}(1:NCOL);
end

BOTM = nan(NROW, NCOL, NLAY);
for ilay = 1: NLAY
    line0 = fgets(fid); 
    for irow = 1: NROW
        line0 = fgets(fid);
        D = textscan(line0, '%f');
        BOTM(irow,:,ilay) = D{1}(1:NCOL);
    end
end
fclose(fid);

if exist('mask_fil', 'var')
    id = fopen(mask_fil, 'r');
    D = textscan(fid, '%s %f', 6); 
    NSEW = D{2}(1:4);
    NROW = D{2}(5);
    NCOL = D{2}(6);
    D = textscan(fid, '%f'); 
    IBOUND = reshape(D{1}, NCOL, NROW)'; % NROW x NCOL
    D = textscan(fid, '%s %s %f %s %f'); 
    dischargePt_rowi = D{3};
    dischargePt_coli = D{5};
    fclose(fid);
    TOP(IBOUND<=0) = 0;
    for ii = 1:NLAY
        X = BOTM(:,:,ii);
        X(IBOUND<=0) = 0;
        BOTM(:,:,ii) = X;
    end
end

% calculate max downward elev gradient
zgrad = zeros(NROW,NCOL);
TOP0 = TOP; TOP0(TOP<=1) = nan; %TOP0(TOP>1) = 1;
Y = TOP0;
X = zeros(NROW-2,NCOL-2,4);
X(:,:,1) = Y(2:end-1,2:end-1)-Y(1:end-2,3:end);
X(:,:,2) = Y(2:end-1,2:end-1)-Y(3:end,2:end-1);
X(:,:,3) = Y(2:end-1,2:end-1)-Y(2:end-1,1:end-2);
X(:,:,4) = Y(2:end-1,2:end-1)-Y(2:end-1,3:end);

zgrad(2:end-1,2:end-1) = max(X,[],3);
zgrad(TOP <= 1) = nan;

%%
figure
subplot(2,2,1)
imagesc(TOP), 
m = TOP(TOP>0); m = min(m(:));
caxis([m*0.9, max(BOTM(:))]), 
cm = colormap;
cm(1,:) = [1 1 1];
colormap(cm);
colorbar
title('TOP')

for ilay = 1:NLAY
    subplot(2,2,1+double(ilay))
    m = BOTM(BOTM>0); m = min(m(:));
    imagesc(BOTM(:,:,ilay)), 
    caxis([m*0.9, max(BOTM(:))]), 
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
    m = X(X>=0); m = min(m(:));
    imagesc(X), caxis([m*0.9, max(X(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['DZ', ' lay', num2str(ilay)]);
end

for ilay = 1:NLAY
    if ilay == 1
        X = TOP - BOTM(:,:,ilay);
    else
        X = X + BOTM(:,:,ilay-1)-BOTM(:,:,ilay);
    end
end    
subplot(2,2,3)
m = X(X>=0); m = min(m(:));
imagesc(X), caxis([m*0.9, max(X(:))]), 
cm = colormap;
cm(1,:) = [1 1 1];
colormap(cm);
colorbar
title(['sum(DZ)']);

figure
subplot(2,2,1)
X = zgrad;
X(isnan(X)) = -1000;
m = X(X>=0); m = min(m(:));
imagesc(X), 
caxis([m-5, max(X(:))]), 
% caxis([1000, max(X(:))]), 
cm = colormap;
cm(1,:) = [1 1 1];
colormap(cm);
colorbar
title(['max downward elev gradient']);