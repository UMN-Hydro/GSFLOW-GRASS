% Plot_PRMS_timeseries.m
%
% plots variables from StatVarFile (see 

clear all, close all, fclose all;
fig0 = 0;

% Output Statistics file name
% stat_file = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/ESCI5980/merced/output/XYZ.statvar';
% stat_file = '/home/gcng/workspace/Models/PRMS/prms4.0.1_linux/projects/merced/output/IDE.statvar';
% stat_file = '/home/gcng/workspace/ModelRuns_scratch/PRMS_projects/Salta/output/ToroOut.statvar';
stat_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/PRMS/ChimOut.statvar';

% ('basin_cfs' is simulated runoff, 'runoff' is observed)
plotVar = {'basin_cfs', 'basin_rain'};
% plotVar = {'basin_cfs', 'runoff', 'basin_rain'};
% plotVar = {'basin_cfs', 'basin_rain', 'basin_storage', 'basin_pweqv', 'basin_soil_moist', 'basin_ssstor', 'basin_gwstor', 'basin_imperv_stor', 'basin_intcp_stor'};
% plotVar = {'basin_cfs', 'basin_rain', 'basin_imperv_stor', 'basin_intcp_stor', 'basin_gwstor', 'basin_ssstor', 'basin_soil_moist', 'basin_pweqv', 'basin_storage'};
% plotVar = {'basin_ppt', 'basin_actet', 'basin_imperv_stor', 'basin_intcp_stor', 'basin_gwstor', 'basin_ssstor', 'basin_soil_moist', 'basin_pweqv', 'basin_storage', 'basin_soil_to_gw'};
% plotVar = {'basin_ppt', 'basin_intcp_evap', 'basin_actet', 'basin_potet', 'basin_infil', 'basin_soil_to_gw', 'basin_sroff_cfs', 'basin_ssflow_cfs', 'basin_cfs'};
lnstr = {'-', '-', '-', '-', '-', '-', '-', '-', '--', '-'};
figID = ones(length(plotVar),1);
% lnstr = cell(length(plotVar),1); lnstr(:) = '-';
% figID = [1: length(plotVar)];
% figID(2) = 2;
% figID(3:end) = 3;
figID(end) = 2;

%% ------------------------------------------------------------------------
% Read in statvar output file
%  Final variables with output results:
%   - varNames: cell with variable names [nVarNames,1]
%   - varIndex: array with index number for variables (1 to ndim) [nVarNames,1]
%   - yr, mon, day: arrays [nDays, 1]
%   - data_array: [nDays, nVarNames]

fid = fopen(stat_file, 'r');
numVar = fscanf(fid, '%d', 1);

% variable names 
D = textscan(fid, '%s %d', numVar);
varNames = D{1};
varIndex = D{2};

% read in date and data
fmt = [repmat('%f ', 1, 7+numVar)];
D = textscan(fid, fmt);
yr = D{2}; mon = D{3}; day = D{4}; 
data_array = cell2mat(D(8:end)); % [ndays x numVar]

fclose(fid);

%% ------------------------------------------------------------------------
datevec_all = [yr, mon, day];
datenum_all = datenum(datevec_all);

% Plot specified variables
prow = 3;
pcol = 1;
meanann = nan(length(plotVar),1);
for ii = 1: length(plotVar)
    % get time series to plot
    Y = data_array(:, strcmp(plotVar{ii}, varNames));
    meanann(ii) = mean(Y)*365;  
    
    % get plot position
    fig_i = fig0+ceil(figID(ii)/(prow*pcol));
    plot_i = mod(figID(ii)-1, (prow*pcol))+1;
    figure(fig_i), orient landscape
    subplot(prow, pcol, plot_i),
    plot(datenum_all, Y, 'LineStyle', lnstr{ii});
    hold on
    
    % get tick marks
    ind_Jan1 = find(mon == 1 & day == 1);
    ind_mon1 = find(day == 1);
    
    ind0 = ind_mon1;
    incr = ceil(length(ind0)/6);
    ticks = ind0([1:incr:length(ind0)]);
    ticklab = cell(length(ticks),1);
    ctr = 0;
    for jj = ticks(:)'
        ctr = ctr+1;
        ticklab{ctr} = [num2str(mon(jj)), '/', num2str(day(jj)), '/', num2str(yr(jj))]; 
    end
    set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
    set(gca, 'XLim', datenum_all([1,end]));
    
end

% create legend
ctr = 0;
for ii = 1: fig_i
    for jj = 1: prow*pcol
        ctr = ctr + 1;
        figure(ii+fig0), subplot(prow, pcol, jj),
        ind = find(figID==ctr);
        nlines = length(ind);
        legstr = cell(nlines,1);
        for kk = 1: nlines
            legstr{kk} = [plotVar{ind(kk)}, ' ', num2str(meanann(ind(kk)))];
        end
        legend(legstr);
        if ctr == max(figID)
            return
        end                    
    end
end
