% PRMS_print_climate_hru_files.m
% (12/8/15)

% Creates climate data inputs files to run PRMS for climate_hru mode using 
% by uniformly applying data from a single weather station:
%   Specify precip, tmax, tmin, and solrad for each HRU for each
%   simulation day.  Can have start day earlier and end day later than
%   simulation period.
% 

% -- Number of HRU's over which to apply data uniformly
nhru = 100;

% -- These files will be created
empty_datafil = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/inputs/empty.dat';
precip_datafil = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/inputs/precip.dat';
tmax_datafil = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/inputs/tmax.dat';
tmin_datafil = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/inputs/tmin.dat';
solrad_datafil = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/inputs/solrad.dat';

% -- Set any description of data here (example: location, units)
descr_str = 'Some station in Salta';

% -- *** FILL IN *** Read in daily values from 1 station for:
%   - precip [mm] (check 'precip_units') - check!!!
%   - tmax [degC] (check 'temp_units') - check!!!
%   - tmin [degC] (check 'temp_units') - check!!!
%   - swrad [langleys] (check 'temp_units') - check!!!
%   - ymdhms_v
in_datafil = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/inputs/SaltaTestMetInputs.csv';
fid = fopen(in_datafil, 'r');
fmt = ['%s %f %f %f %f %f %f %f'];
D = textscan(fid, fmt, 'HeaderLines', 1, 'Delimiter', ',');
fclose(fid);
ymdhms_v = datevec(D{1}, 'dd/mm/yyyy');
precip = D{3} /10 /2.54; % mm -> inch
% tmax = D{4} *(9/5)+32; % C -> F
% tmin = D{5} *(9/5)+32; % C -> F
tmax = D{4}; % C 
tmin = D{5}; % C
swrad = D{end}*1e6/41840; % MJ/m2 -> Langeley (1 Langely = 41840 J/m2)

%% ------------------------------------------------------------------------

print_fmt1 = ['%4d ', repmat('%2d ', 1, 5), ' \n'];
print_fmt2 = ['%4d ', repmat('%2d ', 1, 5), repmat('%6.2f ', 1, nhru), ' \n'];
for ii = 0: 4
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
            label = {['tmax ', num2str(nhru)]};
        case 3, 
            outdatafil = tmin_datafil;
            data = tmin;
            label = {['tmin ', num2str(nhru)]};
        case 4, 
            outdatafil = solrad_datafil;
            data = swrad;
            label = {['swrad ', num2str(nhru)]};
            % ***NOTE: if you get errors related to swrad or orad, use instead:
            % 'orad'
            % and/or try setting orad_flag=0 in control file
    end    
    
    data0 = [ymdhms_v, repmat(data, 1, nhru)];
    
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

