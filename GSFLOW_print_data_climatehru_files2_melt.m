% GSFLOW_print_climate_hru_files1.m
% 11/24/16
%
% based on: PRMS_print_climate_hru_files2.m
%
% v2 - replicates record for very long, continuous time series
% _melt - includes extra high precip for highest elev hru

% Creates climate data inputs files to run PRMS for climate_hru mode using 
% by uniformly applying data from a single weather station:
%   Specify precip, tmax, tmin, and solrad for each HRU for each
%   simulation day.  Can have start day earlier and end day later than
%   simulation period.
% 


% *** CUSTOMIZE TO YOUR COMPUTER! *****************************************
% NOTE: '/' is directory separator for Linux, '\' for Windows!!
inname = ''; % test name

% directory for PRMS input files (include slash ('/') at end) - cbh files
% outputed to here
PRMSinput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/PRMS/';

% PRMSinput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/ChimTest';

% directory with files to be read in to generate PRMS input files in_GISdata_dir
in_data_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/PRMS/';
in_climatedata_dir = strcat(in_data_dir, 'climate/'); % specifically climate data

% GIS data (for nhru)
in_data_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/';
in_GISdata_dir = strcat(in_data_dir, 'GIS/'); % specifically GIS data
HRUfil = strcat(in_GISdata_dir, 'HRU.csv');

% *************************************************************************

% Project-specific entries ->

% -- Number of HRU's over which to apply data uniformly
% nhru should be generated dynamically (from GIS data)
HRUdata = importdata(HRUfil,',',1);
nhru = max(HRUdata.data(:, strcmp(HRUdata.colheaders, 'id')));

% -- Set any description of data here (example: location, units)
descr_str = 'BocaToma';

% -- Read in daily values from 1 station for:
%   - precip (check 'precip_units')
%   - tmax (check 'temp_units')
%   - tmin (check 'temp_units') 
%   [- swrad [langleys for F temp_units] ]-currently unavailable at Chim
%   - ymdhms_v
in_datafil = strcat(in_climatedata_dir, 'test_boca_toma_F_in.txt');
fid = fopen(in_datafil, 'r');
while(1)
    line0 = fgets(fid);
    if strncmp(line0, '####', 4)
        break
    end
end
fmt = ['%d %d %d %d %d %d %f %f %f'];
D = textscan(fid, fmt);
fclose(fid);
ymdhms_v = [D{1} D{2} D{3} D{4} D{5} D{6}];
tmax = D{7}; % must be F
tmin = D{8}; % must be F
precip = D{9}; % must be inch
% swrad = D{end}*1e6/41840; % MJ/m2 -> Langeley (1 Langely = 41840 J/m2)

% - how many times to replicate (only loop full years, starting from first record)
% (if nrep = 1 and <full year, then just print record)
nrep = 30; 

% - get hru elevation info
in_data_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/';
in_GISdata_dir = strcat(in_data_dir, 'GIS/'); % specifically GIS data
HRUfil = strcat(in_GISdata_dir, 'HRU.csv');
% load elevations of hru's
HRUdata = importdata(HRUfil,',',1);
hru_elev = HRUdata.data(:, strcmp(HRUdata.colheaders, 'hru_elev'));  % elev_units = m

% - meltwater inputs (to be added to precip at hrui_melt)
% BT 6/2015-6/2016: mean precip 0.11 in/d 
[~,hrui_melt] = max(hru_elev); % 
melt_precipequiv = 0.1; % [in/day]
melt_precipequiv = 1*0.1; % [in/day]

% temp change w/ elev, apply additive adjustment to temp, BT -> R  
% R is ~4600m, BT is ~3900m
temp_thresh = 4500; % [m]
hrui_tempadjust = find(hru_elev > temp_thresh);
tmax_adj = 5/9*-9.44; % see stefan's email 12/2/16
tmin_adj = 5/9*-2.4; % see stefan's email 12/2/16

% - warmer conditions (=0 to skip)
% tmax_add_warmer = 0; % add this many degF 
% tmin_add_warmer = 0; % add this many degF
% (see Bradley et al. 2006, 0.11degC/decade)
tmax_add_warmer = 1*9/5; % add this many degF 
tmin_add_warmer = 1*9/5; % add this many degF


%% ------------------------------------------------------------------------
% Generally, do no change below here

datenum_v = datenum(double(ymdhms_v));

% -- Replicate complete year(s) of data data nrep times (just ignore leap
% year issues)
if nrep > 1
    ymdhms_prev = datevec(datenum_v(1)-1);
    
    % get complete years record only:
    ind_end = find(ymdhms_v(:,2)==ymdhms_prev(2) & ymdhms_v(:,3)==ymdhms_prev(3), 1, 'last');
    if isempty(ind_end), fprintf('Warning! Not replicating bc there is < 1 full year!\n'); end 
    ymdhms_v = ymdhms_v(1:ind_end,:);
    tmax = tmax(1:ind_end); tmin = tmin(1:ind_end); precip = precip(1:ind_end); 
    
    % replicate nrep times
    datenum_v = datenum_v(1) + [0: nrep*size(ymdhms_v,1)-1];

    ymdhms_v = datevec(datenum_v);
    tmax = repmat(tmax, nrep, 1); 
    tmin = repmat(tmin, nrep, 1);
    precip = repmat(precip, nrep, 1);
end


% -- These files will be generated (names must match those in Control file!)
%    (generally don't change this)
empty_datafil = strcat(PRMSinput_dir, inname, '/empty_rep30yr.day');
precip_datafil = strcat(PRMSinput_dir, inname, '/precip_rep30yr_melt.day');
tmax_datafil = strcat(PRMSinput_dir, inname, '/tmax_rep30yr_tadj_plus1C.day');
tmin_datafil = strcat(PRMSinput_dir, inname, '/tmin_rep30yr_tadj_plus1C.day');
hum_datafil = strcat(PRMSinput_dir, inname, '/humidity_rep30yr.day'); % fake! 
solrad_datafil = strcat(PRMSinput_dir, inname, '/swrad_rep30yr.day');

% - how many of above met variables to process (0 thru N)
N = 4;

% - Write to data variable files
print_fmt1 = ['%4d ', repmat('%2d ', 1, 5), ' \n'];
print_fmt2 = ['%4d ', repmat('%2d ', 1, 5), repmat('%6.2f ', 1, nhru), ' \n'];
if exist('swrad', 'var'), N = 4; end 
for ii = 0: N
    switch ii
        case 0, 
            outdatafil = empty_datafil;
            data = [];
            label = {'precip 0', 'tmax 0', 'tmin 0'};
        case 1, 
            outdatafil = precip_datafil;
            data = precip;
            label = {['precip ', num2str(nhru)]};
        case 2, 
            outdatafil = tmax_datafil;
            data = tmax + tmax_add_warmer;
            label = {['tmaxf ', num2str(nhru)]};
        case 3, 
            outdatafil = tmin_datafil;
            data = tmin + tmin_add_warmer;
            label = {['tminf ', num2str(nhru)]};
        case 4, 
            outdatafil = hum_datafil;
            data = tmax/100;
            label = {['humidity_hru ', num2str(nhru)]};
        case 5, 
            outdatafil = solrad_datafil;
            data = swrad;
            label = {['swrad ', num2str(nhru)]};
            % ***NOTE: if you get errors related to swrad or orad, use instead:
            % 'orad'
            % and/or try setting orad_flag=0 in control file
    end    
    
    data0 = repmat(data, 1, nhru);
    if strncmp(label, 'precip', 6);
        data0(:,hrui_melt) = data0(:,hrui_melt) + melt_precipequiv;
    elseif strncmp(label, 'tmaxf', 5);
        data0(:,hrui_tempadjust) = data0(:,hrui_tempadjust) + tmax_adj;
    elseif strncmp(label, 'tminf', 5);
        data0(:,hrui_tempadjust) = data0(:,hrui_tempadjust) + tmin_adj;
    end
    data_all = [double(ymdhms_v), data0];

    
    fid = fopen(outdatafil, 'wt');
    fprintf(fid, '%s \n', descr_str);
    for ll = 1: length(label), fprintf(fid, '%s \n', label{ll}); end
    fprintf(fid, '########## \n');  % divider required
    
    if ii == 0
        fprintf(fid, print_fmt1, ymdhms_v');
    else
        fprintf(fid, print_fmt2, data_all');
    end
end

