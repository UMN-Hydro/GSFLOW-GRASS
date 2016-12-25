% read *.dis
clear all, 
% close all, 
fclose all;

% NLAY = 2;
% NROW = 73;
% NCOL = 81;
% ba6_file = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen/input/modflow/sagehen.bas';

NLAY = 2;
NROW = 27;
NCOL = 51;
ba6_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/test2lay.ba6';

fid = fopen(ba6_file);

for ii = 1:2
    cmt = fgets(fid);
end

IBOUND = nan(NROW, NCOL, NLAY);
for ilay = 1: NLAY
    line0 = fgets(fid); 
    for irow = 1: NROW
        line0 = fgets(fid);
        D = textscan(line0, '%f');
        IBOUND(irow,:,ilay) = D{1}(1:NCOL);
    end
end

line0 = fgets(fid); 
D = textscan(line0, '%f %s');
HNOFLO = D{1};

initHead = nan(NROW, NCOL, NLAY);
for ilay = 1: NLAY
    line0 = fgets(fid); 
    for irow = 1: NROW
        line0 = fgets(fid);
        D = textscan(line0, '%f');
        initHead(irow,:,ilay) = D{1}(1:NCOL);
    end
end


fclose(fid);

% make plots

X = IBOUND;
figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = X(X>0); m = min(m(:));
    imagesc(X(:,:,ilay)), 
%     caxis([m*0.9, max(X(:))]), 
    cm = colormap;
%     cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['IBOUND', ' lay', num2str(ilay)]);
end

initHead(initHead<10) = 0; % some outer cells with initHead=1
X = initHead;
figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = X(X>0); m = min(m(:));
    imagesc(X(:,:,ilay)), 
    caxis([m*0.9, max(X(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['init head', ' lay', num2str(ilay)]);
end

initHead2 = initHead;
initHead2(IBOUND==0) = 0; % some outer cells with initHead=1
X = initHead2;
figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    m = X(X>0); m = min(m(:));
    imagesc(X(:,:,ilay)), 
    caxis([m*0.9, max(X(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['init head', ' lay', num2str(ilay)]);
end


% check: initHead only specified for IBOUND~=0? (IBOUND>0 is variable, <0
% is constnat head)
a = xor(IBOUND, initHead); % xor founds where exactly one is non-zero
X = a;
figure
for ilay = 1:NLAY
    subplot(2,2,double(ilay))
    %m = X(X>0); 
    m = X; m = min(m(:));
    imagesc(X(:,:,ilay)), 
    caxis([m*0.9, max(X(:))]), 
    cm = colormap;
    cm(1,:) = [1 1 1];
    colormap(cm);
    colorbar
    title(['initHead and ibound differ', ' lay', num2str(ilay)]);
end

