% GSFLOW_print_climate_hru_files1.m
% 11/24/16
%
% based on: PRMS_print_climate_hru_files2.m

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

%% ------------------------------------------------------------------------
% Generally, do no change below here

% -- These files will be generated (names must match those in Control file!)
%    (generally don't change this)
empty_datafil = strcat(PRMSinput_dir, inname, '/empty.day');
precip_datafil = strcat(PRMSinput_dir, inname, '/precip.day');
tmax_datafil = strcat(PRMSinput_dir, inname, '/tmax.day');
tmin_datafil = strcat(PRMSinput_dir, inname, '/tmin.day');
hum_datafil = strcat(PRMSinput_dir, inname, '/humidity.day'); % fake! 
solrad_datafil = strcat(PRMSinput_dir, inname, '/swrad.day');

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
            data = tmax;
            label = {['tmaxf ', num2str(nhru)]};
        case 3, 
            outdatafil = tmin_datafil;
            data = tmin;
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
    
    data0 = [double(ymdhms_v), repmat(data, 1, nhru)];
    
    fid = fopen(outdatafil, 'wt');
    fprintf(fid, '%s \n', descr_str);
    for ll = 1: length(label), fprintf(fid, '%s \n', label{ll}); end
    fprintf(fid, '########## \n');  % divider required
    
    if ii == 0
        fprintf(fid, print_fmt1, ymdhms_v');
    else
        fprintf(fid, print_fmt2, data0');
    end
end

