% make_sfr.m
% 1/8/16
% Leila Saberi
%
% 2 - gcng
function make_sfr2_f(GSFLOW_indir, infile_pre, reach_fil, segment_fil_all)

% Note: assume .dis file already created!! (reads in TOP for setting STRTOP)

% % ======== TO RUN AS SCRIPT ===============================================
% clear all, close all, fclose all;
% GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
% infile_pre = 'test2lay';
% 
% % for sfr
% GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
% reach_fil = [GIS_indir, 'reach_data.txt'];
% segment_fil_all = cell(3,1);
% segment_fil_all{1} = [GIS_indir, 'segment_data_4A_INFORMATION.txt'];
% segment_fil_all{2} = [GIS_indir, 'segment_data_4B_UPSTREAM.txt'];
% segment_fil_all{3} = [GIS_indir, 'segment_data_4C_DOWNSTREAM.txt'];
% % =========================================================================

%%

sfr_file = [infile_pre, '.sfr'];

% -- Refer to GSFLOW manual p.202, SFR1 manual, and SFR2 manual

% - Refer to Fig. 1 of SFR1 documentation for segment vs. reach numbering

% You need the following inputs (with corresponding structures)

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
data_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';

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
    figure
    imagesc(RCH_NUM), colorbar, 
    colormap(jet(1+max(IREACH)));
    caxis([-0.5 max(IREACH)+0.5])
    figure
    imagesc(SEG_NUM), colorbar, 
    colormap(jet(1+max(ISEG)));
    caxis([-0.5 max(ISEG)+0.5])

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
        
    % -- threshold slope at minimum 0.001
    ind = reach_data_all0(:,strcmp(reach_data_all.colheaders, 'SLOPE')) < 0.001;
    reach_data_all0(ind,strcmp(reach_data_all.colheaders, 'SLOPE')) = 0.001;

    % -- set streambed thickness (Sagehen uses constant 1m)
    ind = strcmp(reach_data_all.colheaders, 'STRTHICK');
    reach_data_all0(:,ind) = 1; % [m]

    % -- set streambed hydraulic conductivity (Sagehen example: 5 m/d)
    ind = strcmp(reach_data_all.colheaders, 'STRHC1');
    reach_data_all0(:,ind) = 5;
    
    % set streambed theta_s
    ind = strcmp(reach_data_all.colheaders, 'THTS');
    reach_data_all0(:,ind) = 0.35;
    
    % set streambed initial theta
    ind = strcmp(reach_data_all.colheaders, 'THTI');
    reach_data_all0(:,ind) = 0.3;
    
    % set streambed Brooks-Corey exp (sagehen example is 3.5)
    ind = strcmp(reach_data_all.colheaders, 'EPS');
    reach_data_all0(:,ind) = 3.5;
    
    % set streambed unsaturated zone saturated hydraulic conductivity
    % (sagehen example is 0.3 m/d)
    ind = strcmp(reach_data_all.colheaders, 'UHC');
    reach_data_all0(:,ind) = 0.3;

    reach_data_all = reach_data_all0;
end

if isstruct(segment_data_4A)
    segment_data_4A = segment_data_4A.data;
    nss = size(segment_data_4A,1);   %Number of stream segments
end
if isstruct(segment_data_4B)
    segment_data_4B = segment_data_4B.data;
end
if isstruct(segment_data_4C)
    segment_data_4C = segment_data_4C.data;
end

% if isstruct(stress_periods)
%     stress_periods = stress_periods.data;
% end
% - specify only for 2 stress periods:
stress_periods = zeros(NPER, 3); % itmp, irdflg, iptflg (latter 2 are set to 0)
stress_periods(1,1) = nss;
if NPER > 1, stress_periods(2:end,1) = -1; end

% -------------------------------------------------------------------------

% First put 4A, 4B and 4C data all together in a cell array
% size(cell) = nitems x 1 x nperiods
% In this case, nitems is 3 (i.e. 4A, 4B and 4C)

nitems = 3;
nperiods = size(stress_periods, 1);
segment_data_all = cell(nitems, 1, nperiods);

segment_data_all{1, 1, 1} = segment_data_4A;
segment_data_all{2, 1, 1} = segment_data_4B;
segment_data_all{3, 1, 1} = segment_data_4C;


% -------------------------------------------------------------------------
% validate some of the input data

msg_invalidISFROPT = ['Error: ISFROPT should be set to an integer of', ...
    '1, 2, 3, 4 or 5.'];

if (nstrm < 0)
    if ~ismember(isfropt, [1, 2, 3, 4, 5])
        error(msg_invalidISFROPT);
    end
end


msg_notSupport = ['Error: %s: this variable must be zero because ', ...
    'parameters are not supported in GSFLOW.'];

if (nsfrpar ~= 0)
    error(msg_notSupport, 'NSFRPAR');
end

if (nparseg ~= 0)
    error(msg_notSupport, 'NPARSEG');
end


% -------------------------------------------------------------------------

% Ouput file
fid = fopen([GSFLOW_indir, '/', sfr_file], 'wt');


% Write header lines (item 0)
heading = '# Streamflow-Routing (SFR7) input file.\n';
fprintf(fid, heading);
fprintf(fid, '# %s simulation -- created on %s.\n', upper(project_name), date);


% Item 1
fprintf(fid, '  %5d  %5d  %5d  %5d  %8.2f  %8.4f  %5d  %5d', ...
    nstrm, nss, nsfrpar, nparseg, const, dleak, istcb1, istcb2);

if (isfropt >= 1)
    fprintf(fid, '  %5d', isfropt);
    if (isfropt == 1)
        fprintf(fid, '  %5d\n', irtflg);
    elseif (isfropt > 1)
        fprintf(fid, '  %5d  %5d  %5d  %5d\n', nstrail, isuzn, nsfrsets, irtflg);
    end
else
    fprintf(fid, '\n');
end


% Item 2
if (isfropt == 1)
    ncols_reach = 10;
elseif (isfropt == 2)
    ncols_reach = 13;
elseif (isfropt == 3)
    ncols_reach = 14;
else
    ncols_reach = 6;
end
reach_data_copy = reach_data_all(:, 1:ncols_reach);

p = ncols_reach - 5;
fmt_reach = [repmat('  %5d', 1, 5), repmat('  %8.3f', 1, p), '\n'];

for istrm=1:abs(nstrm)
    dummy = reach_data_copy(istrm, :);
    fprintf(fid, fmt_reach, dummy);
end


% Item 3 and 4
nper = size(stress_periods, 1);
for iper=1:nper
    
    % write item 3 to the file
    dummy3 = num2cell(stress_periods(iper, :));
    [itmp, irdflg, iptflg] = dummy3{:};
    fprintf(fid, '  %5d  %5d  %5d\n', itmp, irdflg, iptflg);
    
    if (itmp > 0)
        seg_inf_4a = segment_data_all{1, 1, iper};
        seg_inf_4b = segment_data_all{2, 1, iper};
        seg_inf_4c = segment_data_all{3, 1, iper};
        
        for iitmp=1:itmp   % start loop over itmp (num_segments)
            
            % write item 4a to the file
            dummy4a = num2cell(seg_inf_4a(iitmp, :));
            [nseg, icalc, outseg, iupseg, iprior, nstrpts, ...
                flow, runoff, etsw, pptsw, roughch, roughbk, ...
                cdpth, fdpth, awdth, bwdth] = dummy4a{:};
            
            fmt = [' ', repmat('  %5d', 1, 4)];
            fprintf(fid, fmt, nseg, icalc, outseg, iupseg);
            
            if (iupseg > 0)
                fprintf(fid, '  %5d', iprior);
            end
            
            if (icalc == 4)
                fprintf(fid, '  %5d', nstrpts);
            end
            
            fmt = repmat('  %8.3f', 1, 4);
            fprintf(fid, fmt, flow, runoff, etsw, pptsw);
            
            if ((icalc == 1) || (icalc == 2))
                fprintf(fid, '  %8.3f', roughch);
            end
            
            if (icalc == 2)
                fprintf(fid, '  %8.3f', roughbk);
            end
            
            if (icalc == 3)
                fmt = repmat('  %8.3f', 1, 4);
                fprintf(fid, fmt, cdpth, fdpth, awdth, bwdth);
            end
            
            fprintf(fid, '\n');
            
            % write items 4b and 4c to the file
            suffixes = 'bc';
            for i=1:2   % start loop over i
                suffix = suffixes(i);
                var_str = ['seg_inf_4', suffix];
                var = eval(var_str);
                
                dummy4bc = num2cell(var(iitmp, :));
                [hcond, thickm, elevupdn, width, ...
                    depth, thts, thti, eps, uhc] = dummy4bc{:};
                
                fl_no_4bc = 0;
                if (ismember(isfropt, [0, 4, 5]) && (icalc <= 0))
                    fmt = ['    ', repmat('  %8.3f', 1, 5)];
                    fprintf(fid, fmt, hcond, thickm, elevupdn, width, depth);
                    
                elseif (ismember(isfropt, [0, 4, 5]) && (icalc == 1))
                    fprintf(fid, '    %8.3f', hcond);
                    
                    if (iper == 1)   % only for the first period
                        fmt = repmat('  %8.3f', 1, 3);
                        fprintf(fid, fmt, thickm, elevupdn, width);
                        
                        if ((isfropt == 4) || (isfropt == 5))
                            fmt = repmat('  %8.3f', 1, 3);
                            fprintf(fid, fmt, thts, thti, eps);
                        end
                        
                        if (isfropt == 5)
                            fprintf(fid, '  %8.3f', uhc);
                        end
                        
                    elseif ((iper > 1) && (isfropt == 0))
                        fmt = repmat('  %8.3f', 1, 3);
                        fprintf(fid, fmt, thickm, elevupdn, width);
                    end
                    
                elseif (ismember(isfropt, [0, 4, 5]) && (icalc >= 2))
                    fprintf(fid, '    %8.3f', hcond);
                    
                    if ~(ismember(isfropt, [4, 5]) && (iper > 1) && (icalc == 2))
                        fprintf(fid, '  %8.3f  %8.3f', thickm, elevupdn);
                        
                        if (ismember(isfropt, [4, 5]) && (iper == 1) && (icalc == 2))
                            fmt = repmat('  %8.3f', 1, 3);
                            fprintf(fid, fmt, thts, thti, eps);
                            
                            if (isfropt == 5)
                                fprintf(fid, '  %8.3f', uhc);
                            end
                        end
                    end
                    
                elseif ((isfropt == 1) && (icalc <= 1))
                    fprintf(fid, '    %8.3f', width);
                    
                    if (icalc <= 0)
                        fprintf(fid, '  %8.3f', depth);
                    end
                    
                elseif (ismember(isfropt, [2, 3]) && (icalc <= 1))
                    if (iper == 1)
                        fprintf(fid, '    %8.3f', width);
                        
                        if (icalc <= 0)
                            fprintf(fid, '  %8.3f', depth);
                        end
                    end
                else
                    fl_no_4bc = 1;
                end
                if ~fl_no_4bc
                    fprintf(fid, '\n');
                end
            end   % terminate loop over i (4b and 4c, respectively)

        end   % terminate loop over itmp (num_segments)
        
    end   % enf if (itmp > 0)
    
end   % terminate loop over iper (num_periods)

fclose(fid);

% -------------------------------------------------------------------------
% End of the script
