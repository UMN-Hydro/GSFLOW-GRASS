% make_sfr.m
% 1/8/16
% Leila Saberi
%
% 2 - gcng
% function make_sfr2(GSFLOW_indir, infile_pre)

% -------------------------------------------------------------------------
GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
infile_pre = 'test2lay';

% for sfr
GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
reach_fil = [GIS_indir, 'reach_data.txt'];
segment_fil_all = cell(3,1);
segment_fil_all{1} = [data_indir, 'segment_data_4A_INFORMATION.txt'];
segment_fil_all{2} = [data_indir, 'segment_data_4B_UPSTREAM.txt'];
segment_fil_all{3} = [data_indir, 'segment_data_4C_DOWNSTREAM.txt'];

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
    reach_data_all = reach_data_all.data;
    nstrm = size(reach_data_all,1);
    if fl_nstrm < 0, nstrm = -nstrm; end 
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
fmt_reach = [repmat('  %5d', 1, 5), repmat('  %8.2f', 1, p), '\n'];

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
            
            fmt = repmat('  %8.2f', 1, 4);
            fprintf(fid, fmt, flow, runoff, etsw, pptsw);
            
            if ((icalc == 1) || (icalc == 2))
                fprintf(fid, '  %8.2f', roughch);
            end
            
            if (icalc == 2)
                fprintf(fid, '  %8.2f', roughbk);
            end
            
            if (icalc == 3)
                fmt = repmat('  %8.2f', 1, 4);
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
                    fmt = ['    ', repmat('  %8.2f', 1, 5)];
                    fprintf(fid, fmt, hcond, thickm, elevupdn, width, depth);
                    
                elseif (ismember(isfropt, [0, 4, 5]) && (icalc == 1))
                    fprintf(fid, '    %8.2f', hcond);
                    
                    if (iper == 1)   % only for the first period
                        fmt = repmat('  %8.2f', 1, 3);
                        fprintf(fid, fmt, thickm, elevupdn, width);
                        
                        if ((isfropt == 4) || (isfropt == 5))
                            fmt = repmat('  %8.2f', 1, 3);
                            fprintf(fid, fmt, thts, thti, eps);
                        end
                        
                        if (isfropt == 5)
                            fprintf(fid, '  %8.2f', uhc);
                        end
                        
                    elseif ((iper > 1) && (isfropt == 0))
                        fmt = repmat('  %8.2f', 1, 3);
                        fprintf(fid, fmt, thickm, elevupdn, width);
                    end
                    
                elseif (ismember(isfropt, [0, 4, 5]) && (icalc >= 2))
                    fprintf(fid, '    %8.2f', hcond);
                    
                    if ~(ismember(isfropt, [4, 5]) && (iper > 1) && (icalc == 2))
                        fprintf(fid, '  %8.2f  %8.2f', thickm, elevupdn);
                        
                        if (ismember(isfropt, [4, 5]) && (iper == 1) && (icalc == 2))
                            fmt = repmat('  %8.2f', 1, 3);
                            fprintf(fid, fmt, thts, thti, eps);
                            
                            if (isfropt == 5)
                                fprintf(fid, '  %8.2f', uhc);
                            end
                        end
                    end
                    
                elseif ((isfropt == 1) && (icalc <= 1))
                    fprintf(fid, '    %8.2f', width);
                    
                    if (icalc <= 0)
                        fprintf(fid, '  %8.2f', depth);
                    end
                    
                elseif (ismember(isfropt, [2, 3]) && (icalc <= 1))
                    if (iper == 1)
                        fprintf(fid, '    %8.2f', width);
                        
                        if (icalc <= 0)
                            fprintf(fid, '  %8.2f', depth);
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
