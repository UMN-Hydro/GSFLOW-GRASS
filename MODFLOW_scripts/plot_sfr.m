% make_sfr.m
% 1/8/16
% Leila Saberi
%
% 2 - gcng
clear all
fclose all;
close all

% the followings are used to write item 1
fl_nstrm = -1;  % flag for stream reaches, <0: include unsaturated zone below (sagehen: >0)
nsfrpar = 0;  %Always Zero
nparseg = 0;    %Always Zero
const = 86400.;  %Conversion factor used in calculating depth for a stream reach (86400 in sagehen example)
dleak = 0.0001;  %Tolerance level of stream depth used in computing leakage between each stream (0.0001 in sagehen example)
istcb1 = -1;    %Flag for writing stream-aquifer leakage values (>0: file unit, <0: write to listing file)
istcb2 = 0;     %Flag for writing to a seperate formatted file information on inflows&outflows
isfropt = 3;    %defines input structure; saturated or non-saturated zone (1: No UZ; 3: UZ, unsat prop at start of simulation), sagehen uses 3
nstrail = 10;    %Number of trailing-waive increments, incr for better mass balance (10-20 rec'd, sagehen uses 8)    
isuzn = 1;   %Maximum number of vertical cells used to define the unsaturated zone beneath a stream reach (for icalc=1 (Mannings for depth): use isuzn=1)
nsfrsets = 40;  %Maximum number of different sets of trailing waves used to allocate arrays.
irtflg = 0;     %Flag whether transient streamflow routing is active

project_name = 'TestProject';                                             % used to name the output file (.sfr)

% data_indir = '/home/gcng/workspace/matlab_files/GSFLOW_pre-processor/MODFLOW_scripts/sfr_final/data/';
% data_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';

% Note: assume .dis file already created!! (reads in TOP for setting STRTOP)

% ======== TO RUN AS SCRIPT ===============================================
% clear all, close all, fclose all;
GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
infile_pre = 'test2lay';

% for sfr
GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
reach_fil = [GIS_indir, 'reach_data.txt'];
segment_fil_all = cell(3,1);
segment_fil_all{1} = [GIS_indir, 'segment_data_4A_INFORMATION.txt'];
segment_fil_all{2} = [GIS_indir, 'segment_data_4B_UPSTREAM.txt'];
segment_fil_all{3} = [GIS_indir, 'segment_data_4C_DOWNSTREAM.txt'];

GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
mask_fil = [GIS_indir, 'basinmask_dischargept.asc']; % only for Chimb


% =========================================================================

%%

% -- Refer to GSFLOW manual p.202, SFR1 manual, and SFR2 manual

% - Refer to Fig. 1 of SFR1 documentation for segment vs. reach numbering

% You need the following inputs (with corresponding structures)

% items 
reach_data_all = importdata(reach_fil);       % used to write item 2: assumes 

NPER = 2;       % used for item 3

% items 4a: # NSEG ICALC  OUTSEG  IUPSEG  IPRIOR  NSTRPTS  FLOW  RUNOFF  ETSW  PPTSW  ROUGHCH  ROUGHBK  CDPTH  FDPTH  AWDTH  BWDTH
segment_data_4A = importdata(segment_fil_all{1});   % used to write items 4a

segment_data_4B = importdata(segment_fil_all{2});      % used to write items 4b (ignored for ICALC=3 in 4a)
segment_data_4C = importdata(segment_fil_all{3});    % used to write items 4c (ignored for ICALC=3 in 4a)


% -------------------------------------------------------------------------
% In case the input text files (e.g. reach_data.txt) contain header lines (comments)

if isstruct(reach_data_all)
    reach_data_all0 = reach_data_all.data;
    nstrm = size(reach_data_all0,1);
    if fl_nstrm < 0, nstrm = -nstrm; end 

    % sort rows according to increasing segment numbers
    [~, ind] = sort(reach_data_all0(:,strcmp(reach_data_all.colheaders, 'ISEG')), 'ascend');
    reach_data_all0 = reach_data_all0(ind,:);
    nss = max(reach_data_all0(:,strcmp(reach_data_all.colheaders, 'ISEG')));

    % sort rows according to increasing reach numbers
    for ii = 1: nss
        ind1 = find(reach_data_all0(:,strcmp(reach_data_all.colheaders, 'ISEG')) == ii);
        [~, ind2] = sort(reach_data_all0(ind1,strcmp(reach_data_all.colheaders, 'IREACH')), 'ascend');
        reach_data_all0(ind1,:) = reach_data_all0(ind1(ind2),:);

        % renumber IREACH to start at 1 for each segment
        reach_data_all0(ind1,strcmp(reach_data_all.colheaders, 'IREACH')) = [1:length(ind1)];
    end
    X = mat2cell(reach_data_all0, abs(nstrm), ones(size(reach_data_all0,2),1));
    [KRCH,IRCH,JRCH,ISEG,IREACH,RCHLEN,STRTOP,SLOPE,STRTHICK,STRHC1,THTS,THTI,EPS,UHC] = X{:};
    
    % -- make sure STRTOP is within 1st layer 
    % - read in TOP and BOTM from .dis file
    dis_file = [GSFLOW_indir, '/', infile_pre, '.dis'];
    fid = fopen(dis_file);
    for ii = 1:2, cmt = fgets(fid); end
    line0 = fgets(fid);
    D = textscan(line0, '%d', 6);
    NLAY = D{1}(1); NROW = D{1}(2); NCOL = D{1}(3);
    NPER = D{1}(4); ITMUNI = D{1}(5); LENUNI = D{1}(6);
    line0 = fgets(fid);
    D = textscan(line0, '%d');
    LAYCBD = D{1}; % 1xNLAY (0 if no confining layer)

    line0 = fgets(fid);
    D = textscan(line0, '%s %d');    DELR = D{2}; % width of column
    line0 = fgets(fid);
    D = textscan(line0, '%s %d');    DELC = D{2}; % height of row

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

    % TOP for cells corresponding to reaches
    TOP_RCH = nan(abs(nstrm),1);
    for ii = 1:abs(nstrm), TOP_RCH(ii) = TOP(IRCH(ii),JRCH(ii)); end

    % BOTM for cells corresponding to reaches
    BOTM_RCH = nan(abs(nstrm),1);
    for ii = 1:abs(nstrm), BOTM_RCH(ii) = BOTM(IRCH(ii),JRCH(ii),KRCH(ii)); end
    
    % - set STRTOP to be just below TOP
    STRTOP = TOP_RCH - 2;
    if ~isempty(find(STRTOP-STRTHICK < BOTM_RCH,1))
        fprintf('Error! STRTOP is below BOTM of the corresponding layer! Exiting...\n');
    end
    reach_data_all0(:,strcmp(reach_data_all.colheaders, 'STRTOP')) = STRTOP; 

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
% 
%         for ii = 1:NLAY 
%             X = HK(2:end-1,2:end-1,ii);
%             X(ind_bound) = 0;
%             HK(2:end-1,2:end-1,ii) = X;
%     %         Y = HK(:,:,ii);
%     %         Y(IBOUND==0) = 0;
%     %         HK(:,:,ii) = Y;
%         end
    end    
    
    % -- plot stream reaches
    RCH_mask = TOP;
    for ii = 1:abs(nstrm), RCH_mask(IRCH(ii),JRCH(ii)) = max(TOP(:))*2; end
    figure
    subplot(2,2,1)
    imagesc(TOP), colorbar, 
    cm = colormap; 
    cm(end,:) = [1 1 1];
    caxis([min(TOP(:)) max(TOP(:))* 1.25]);
    colormap(cm);
    subplot(2,2,2)
    imagesc(RCH_mask), colorbar, 
    cm = colormap; 
    cm(end,:) = [1 1 1];
    caxis([min(TOP(:)) max(TOP(:))* 1.25]);
    colormap(cm);
    
    RCH_NUM = zeros(NROW,NCOL); SEG_NUM = zeros(NROW,NCOL);
    for ii = 1:abs(nstrm), RCH_NUM(IRCH(ii),JRCH(ii)) = IREACH(ii); end
    for ii = 1:abs(nstrm), SEG_NUM(IRCH(ii),JRCH(ii)) = ISEG(ii); end

    % make boundary black
    X = SEG_NUM(2:end-1,2:end-1);
    X(ind_bound) = -1;
    SEG_NUM(2:end-1,2:end-1) = X;
    
    figure
    imagesc(RCH_NUM), colorbar, 
    colormap(jet(1+max(IREACH)));
    caxis([-0.5 max(IREACH)+0.5])
    figure
    imagesc(SEG_NUM), colorbar, 
    cm = colormap(jet(2+max(ISEG)));
    cm(2,:) = [1 1 1]*.8;
    cm(1,:) = [0 0 0];
    colormap(cm)
    caxis([-1.5 max(ISEG)+0.5])
    axis equal
    print('-dsvg', 'segmap.svg');
end

%     % when running as script: to visualize segments one at a time
%     for j = 1: max(ISEG)
%         RCH_NUM = zeros(NROW,NCOL); SEG_NUM = zeros(NROW,NCOL);
%         ind = find(ISEG==j);
%         for ii = ind(:)' 
%             RCH_NUM(IRCH(ii),JRCH(ii)) = IREACH(ii); 
%             SEG_NUM(IRCH(ii),JRCH(ii)) = ISEG(ii); 
%         end                
%         
%         figure(100)
%         imagesc(RCH_NUM), colorbar, 
%         colormap(jet(1+max(RCH_NUM(:))));
%         caxis([-0.5 max(RCH_NUM(:))+0.5])
%         title(['reaches for seg ', num2str(j)]);
%         figure(101)
%         imagesc(SEG_NUM), colorbar, 
%         colormap(jet(1+max(SEG_NUM(:))));
%         caxis([-0.5 max(SEG_NUM(:))+0.5])
%         title(['seg ', num2str(j)]);
%         pause
%     end
