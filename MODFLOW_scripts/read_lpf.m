% read *.lpf
% specifically for sagehen example
clear all, 
close all, fclose all;

% NLAY = 2;
% NROW = 73;
% NCOL = 81;
% lpf_file = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen/input/modflow/sagehen.lpf';

NLAY = 2;
NROW = 27;
NCOL = 51;
lpf_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/spinup_melt_161209a/inputs/MODFLOW/test2lay.lpf';
GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
mask_fil = [GIS_indir, 'basinmask_dischargept.asc']; % only for Chimb

fid = fopen(lpf_file);

for ii = 1:8
    line0 = fgets(fid);
end

ilay = 1; % ---------------------------------------------
HK = nan(NROW,NCOL,NLAY);
line0 = fgets(fid); 
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    HK(irow,:,ilay) = D{1}(1:NCOL);
end

VKA = nan(NROW,NCOL,NLAY);
line0 = fgets(fid); 
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    VKA(irow,:,ilay) = D{1}(1:NCOL);
end

line0 = fgets(fid);
D = textscan(line0, '%s %f %s');
if strcmp(D{1}{1}, 'INTERNAL')
    Ss = nan(NROW,NCOL,NLAY);
    for irow = 1: NROW
        line0 = fgets(fid);
        D = textscan(line0, '%f');
        Ss(irow,:,ilay) = D{1}(1:NCOL);
    end
elseif strcmp(D{1}{1}, 'CONSTANT')
    Ss = nan(NLAY,1);
    Ss(ilay) = D{2};
end

Sy = nan(NROW,NCOL,NLAY);
line0 = fgets(fid); 
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    Sy(irow,:,ilay) = D{1}(1:NCOL);
end

WETDRY = nan(NROW,NCOL,NLAY);
line0 = fgets(fid); 
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    WETDRY(irow,:,ilay) = D{1}(1:NCOL);
end

ilay = 2; % ---------------------------------------------

line0 = fgets(fid); 
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    HK(irow,:,ilay) = D{1}(1:NCOL);
end

line0 = fgets(fid); 
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    VKA(irow,:,ilay) = D{1}(1:NCOL);
end

line0 = fgets(fid);
D = textscan(line0, '%s %f %s');
if strcmp(D{1}, 'INTERNAL')
    for irow = 1: NROW
        line0 = fgets(fid);
        D = textscan(line0, '%f');
        Ss(irow,:,ilay) = D{1}(1:NCOL);
    end
elseif strcmp(D{1}, 'CONSTANT')
    Ss(ilay) = D{2};
end

fclose(fid);

% read in mask if applicable
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
    
    % find boundary cells
    IBOUNDin = IBOUND(2:end-1,2:end-1);
    IBOUNDu = IBOUND(1:end-2,2:end-1); % up
    IBOUNDd = IBOUND(3:end,2:end-1); % down
    IBOUNDl = IBOUND(2:end-1,1:end-2); % left
    IBOUNDr = IBOUND(2:end-1,3:end); % right
%     % - inner boundary is constant head
%     ind_bound = IBOUNDin==1 & (IBOUNDin-IBOUNDu==1 | IBOUNDin-IBOUNDd==1 | ...
%         IBOUNDin-IBOUNDl==1 | IBOUNDin-IBOUNDr==1);
    % - outer boundary is constant head
    ind_bound = IBOUNDin==0 & (IBOUNDin-IBOUNDu==-1 | IBOUNDin-IBOUNDd==-1 | ...
        IBOUNDin-IBOUNDl==-1 | IBOUNDin-IBOUNDr==-1);

    for ii = 1:NLAY 
        X = HK(2:end-1,2:end-1,ii);
        X(ind_bound) = 0;
        HK(2:end-1,2:end-1,ii) = X;
%         Y = HK(:,:,ii);
%         Y(IBOUND==0) = 0;
%         HK(:,:,ii) = Y;
    end
end


figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = HK(HK>0); m = min(m(:));
    imagesc(HK(:,:,ilay)), caxis([m*0.9, max(HK(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['HK', ' lay', num2str(ilay)]);
    axis equal
end
for ilay = 1:NLAY
    subplot(2,2,2+double(ilay))
    m = HK(HK>0); m = min(m(:));
    imagesc(HK(:,:,ilay)), caxis([m*0.9, max(HK(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    caxis([0 10])
    title(['zoom-in HK', ' lay', num2str(ilay)]);
end

for ilay = 1:NLAY
figure
%     subplot(2,2,double(ilay))
    m = HK(HK>0); m = min(m(:));
    imagesc(HK(:,:,ilay)), caxis([m*0.9, max(HK(:))]), 
    cm = colormap;
%     cm(1,:) = [1 1 1];
    cm(1,:) = [0 0 0];
    colormap(cm);
    colorbar
    title(['HK', ' lay', num2str(ilay)]);
    axis equal
    print('-dsvg', ['K_lay', num2str(ilay), '.svg'])
end

figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = VKA(VKA>0); m = min(m(:));
    imagesc(VKA(:,:,ilay)), caxis([m*0.9, max(VKA(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['VKA', ' lay', num2str(ilay)]);
end
for ilay = 1:NLAY
    subplot(2,2,2+double(ilay))
    m = VKA(VKA>0); m = min(m(:));
    imagesc(VKA(:,:,ilay)), caxis([m*0.9, max(VKA(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    caxis([0 10])
    title(['zoom-in VKA', ' lay', num2str(ilay)]);
end

figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = Sy(Sy>0); m = min(m(:));
    imagesc(Sy(:,:,ilay)), caxis([m*0.9, max(Sy(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['Sy', ' lay', num2str(ilay)]);
end

figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = WETDRY(WETDRY>0); m = min(m(:));
    imagesc(WETDRY(:,:,ilay)), caxis([m*0.9, max(WETDRY(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['WETDRY', ' lay', num2str(ilay)]);
end
return

%% re-write horizontal HK and vertical VKA

% HK = HK*0 + 0.01;
% HK(:,:,2) = 0.001;
HK = HK*0 + 5*0.013;
HK(:,:,2) = 3*0.0009;
% HK = HK*0 + 0.5;
% HK(:,:,2) = 0.1;
VKA = HK;

fmt2 = [repmat('%6g ', 1, NCOL), '\n']; %
fid = fopen('HKVKAWETDRY.txt', 'wt');
fprintf(fid, 'INTERNAL 0.013     (FREE)     0                     HORIZONTAL HYDRAULIC CONDUCTIVITY FOR LAYER 1\n');
fprintf(fid, fmt2, HK(:,:,1)/0.013);
fprintf(fid, ' INTERNAL 0.013     (FREE)     0                      VERTICAL HYDRAULIC CONDUCTIVITY FOR LAYER 1\n');
fprintf(fid, fmt2, VKA(:,:,1)/0.013);
fprintf(fid, 'INTERNAL 1.00      (FREE)  0                       WETDRY FOR LAYER 1\n');
fprintf(fid, fmt2, WETDRY(:,:,1));
fprintf(fid, 'INTERNAL 0.0009     (FREE)     0                     HORIZONTAL HYDRAULIC CONDUCTIVITY FOR LAYER 2\n');
fprintf(fid, fmt2, HK(:,:,2)/0.0009);
fprintf(fid, ' INTERNAL 0.0009     (FREE)     0                      VERTICAL HYDRAULIC CONDUCTIVITY FOR LAYER 2\n');
fprintf(fid, fmt2, VKA(:,:,2)/0.0009);
fclose(fid);
cmd = ['cat /home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/data/sagehen/input/modflow/sagehen.mf1.lpf HKVKAWETDRY.txt > /home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/data/sagehen/input/modflow/sagehen.mf_HKVKAWETDRY.lpf'];
system(cmd)


