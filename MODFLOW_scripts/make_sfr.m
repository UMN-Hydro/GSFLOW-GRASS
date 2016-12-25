
% -------------------------------------------------------------------------

% You need the following inputs (with corresponding structures)

% the followings are used to write item 1
nstrm = -201;  %Number of stream reaches, negative integer to be used as a flag
nss = 15;   %Number of stream segments
nsfrpar = 0;  %Always Zero
nparseg = 0;    %Always Zero
const = 86400.;  %Conversion factor used in calculating depth for a stream reach (See P.204 GSFLOW Tutorial)
dleak = 0.0001;  %Tolerance level of stream depth used in computing leakage between each stream
istcb1 = -1;    %Flag for writing stream-aquifer leakage values (See P.204 GSFLOW Tutorial)
istcb2 = 0;     %Flag for writing to a seperate formatted file information on inflows&outflows
isfropt = 3;    %defines input structure; saturated or non-saturated zone (See P.205 GSFLOW Tutorial)
nstrail = 8;    %Number of trailing-waive incerements represent a decrease in the surface infiltration    
isuzn = 1;   %Maximum number of vertical cells used to define the unsaturated zone beneath a stream reach
nsfrsets = 40;  %Maximum number of different sets of trailing waves used to allocate arrays.
irtflg = 0;     %Flag whether transient streamflow routing is active

project_name = 'TestProject';                                             % used to name the output file (.sfr)
reach_data_all = importdata('./data/reach_data.txt');                     % used to write item 2
stress_periods = importdata('./data/stress_periods.txt');                 % used to write item 3
segment_data_4A = importdata('./data/segment_data_4A_INFORMATION.txt');   % used to write items 4a
segment_data_4B = importdata('./data/segment_data_4B_UPSTREAM.txt');      % used to write items 4b
segment_data_4C = importdata('./data/segment_data_4C_DOWNSTREAM.txt');    % used to write items 4b


% -------------------------------------------------------------------------
% In case the input text files (e.g. reach_data.txt) contain header lines (comments)

if isstruct(reach_data_all)
    reach_data_all = reach_data_all.data;
end


if isstruct(stress_periods)
    stress_periods = stress_periods.data;
end


if isstruct(segment_data_4A)
    segment_data_4A = segment_data_4A.data;
end


if isstruct(segment_data_4B)
    segment_data_4B = segment_data_4B.data;
end


if isstruct(segment_data_4C)
    segment_data_4C = segment_data_4C.data;
end


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
fname = [project_name, '.sfr'];
fid = fopen(fname, 'w');


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
                    
                end
                
                fprintf(fid, '\n');
            
            end   % terminate loop over i (4b and 4c, respectively)

        end   % terminate loop over itmp (num_segments)
        
    end   % enf if (itmp > 0)
    
end   % terminate loop over iper (num_periods)

fclose(fid);

% -------------------------------------------------------------------------
% End of the script
