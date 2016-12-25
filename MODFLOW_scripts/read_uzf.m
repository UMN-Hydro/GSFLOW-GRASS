% read_uzf.m
% 11/21/16
%
clear all, 
% close all, 
fclose all;

NLAY = 2;
NROW = 73;
NCOL = 81;

uzf_fil = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/data/sagehen/input/modflow/sagehen.uzf';

fid = fopen(uzf_fil);

for ii = 1:2
    line0 = fgets(fid);
end

line0 = fgets(fid);
D = textscan(line0, '%d %d  %d  %d  %d  %d  %d  %d  %d %f');

NUZTOP = D{1};
IUZFOPT = D{2};
IRUNFLG = D{3};
IETFLG = D{4};
IUZFCB1  = D{5};
IUZFCB2  = D{6};
NTRAIL  = D{7};
NSETS  = D{8};
NUZGAGES  = D{9};
SURFDEP = D{10};

% UZFBND
UZFBND = nan(NROW,NCOL);
line0 = fgets(fid);
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    UZFBND(irow,:) = D{1}(1:NCOL);
end

% VKA
X = nan(NROW,NCOL);
line0 = fgets(fid);
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    X(irow,:) = D{1}(1:NCOL);
end
VKA = X;

line0 = fgets(fid);

% THTS
X = nan(NROW,NCOL);
line0 = fgets(fid);
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    X(irow,:) = D{1}(1:NCOL);
end
THTS = X;


for ii = 1:5
    fgets(fid);
end

% FINF
FINF = nan(NROW,NCOL);
line0 = fgets(fid);
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    FINF(irow,:) = D{1}(1:NCOL);
end

fclose(fid);

figure
subplot(2,2,1)
imagesc(UZFBND),
m = UZFBND(UZFBND>0); m = min(m(:));
caxis([m*0.9, max(UZFBND(:))]), 
cm = colormap;
cm(1,:) = [1 1 1];
colormap(cm);
colorbar
title('UZFBND')

subplot(2,2,2)
X = FINF;
imagesc(X),
m = X(X>0); m = min(m(:));
caxis([m*0.9, max(X(:))]), 
cm = colormap;
cm(1,:) = [1 1 1];
colormap(cm);
colorbar
title('FINF')

%%

% Rewrite FINF, multiplied by 0.0008
FINF2 = FINF*0 + 1.1;
fmt2 = [repmat('%6g ', 1, NCOL), '\n']; %
fid = fopen('FINF.txt', 'wt');
fprintf(fid, fmt2, FINF2);
fprintf(fid, '-1                                                 NUZF1 FOR STRESS PERIOD 2\n');
fclose(fid);
cmd = ['cat /home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/data/sagehen/input/modflow/sagehen1.uzf FINF.txt > /home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/data/sagehen/input/modflow/sagehen_finf.uzf'];
system(cmd)


