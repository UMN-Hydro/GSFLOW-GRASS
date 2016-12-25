% write_ba6_MOD
% 11/17/16

clear all, close all, fclose all;

% - write to this file
GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/';
ba6_file = 'test.ba6';
slashstr = '/';

% - domain dimensions, maybe already in surfz_fil and botm_fil{}?
NLAY = 1;
% NROW = 50;
% NCOL = 50;

% -- Base IBOUND and init head on data from files
% (***TEMPORARY place-holder)
surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
fid = fopen(surfz_fil, 'r');
D = textscan(fid, '%s %f', 6); 
NSEW = D{2}(1:4);
NROW = D{2}(5);
NCOL = D{2}(6);

% - space discretization
DELR = (NSEW(3)-NSEW(4))/NCOL; % width of column [m]
DELC = (NSEW(1)-NSEW(2))/NROW; % height of row [m]
DZ = 10; % [NLAYx1] ***temporary: constant 10m thick single aquifer (consider 2-layer?)

% - set TOP to surface elevation [m]
D = textscan(fid, '%f'); 
fclose(fid);
TOP = reshape(D{1}, NCOL, NROW)'; % NROW x NCOL

BOTM = zeros(NROW, NCOL, NLAY);
for ilay = 1:NLAY
    BOTM(:,:,ilay) = TOP-DZ(ilay);
end

% - IBOUND(NROW,NCOL,NLAY): <0 const head, 0 no flow, >0 variable head
%  (set IBOUND>0 within watershed, =0 outside watershed)
IBOUND = ones(NROW,NCOL,NLAY); 
% - initHead(NROW,NCOL,NLAY)
initHead = BOTM(:,:,1) + (TOP-BOTM(:,:,1))*0.7; % within top layer

% - assumed values
HNOFLO = 999.99;


%% ------------------------------------------------------------------------
% -- Write ba6 file
fil_ba6_0 = [GSFLOW_dir, slashstr, ba6_file];
fmt1 = [repmat('%4d ', 1, NCOL), '\n']; % for IBOUND 
fmt2 = [repmat('%7g ', 1, NCOL), '\n']; % for initHead

fid = fopen(fil_ba6_0, 'wt');
fprintf(fid, '# basic package file --- %d layers, %d rows, %d columns\n', NLAY, NROW, NCOL);
fprintf(fid, 'FREE\n');
for ilay = 1: NLAY
    fprintf(fid, 'INTERNAL          1 (FREE)  3         IBOUND for layer %d \n', ilay); % 1: CNSTNT multiplier, 3: IPRN>0 to print input to list file
    for irow = 1:NROW
        fprintf(fid, fmt1, IBOUND(irow,:,ilay));
    end
end
fprintf(fid, '    %f  HNOFLO\n', HNOFLO);
for ilay = 1: NLAY
    fprintf(fid, 'INTERNAL          1 (FREE)  3         init head for layer %d \n', ilay); % 1: CNSTNT multiplier, 3: IPRN>0 to print input to list file
    for irow = 1: NROW
        fprintf(fid, fmt2, initHead(irow,:,ilay));
    end
end
fclose(fid);

% -- Plot basics
for ii = 1:2
    if ii == 1, 
        X0 = IBOUND; ti0 = 'IBOUND';
    elseif ii == 2
        X0 = initHead; ti0 = 'init head';
    end
    figure
    for ilay = 1:NLAY
        subplot(2,2,double(ilay))
        X = X0(:,:,ilay);
        m = X(X>0); m = min(m(:));
        imagesc(X), caxis([m*0.9, max(X(:))]), 
        cm = colormap;
        cm(1,:) = [1 1 1];
        colormap(cm);
        colorbar
        title([ti0, ' lay', num2str(ilay)]);
    end
end