% plot_gsflow_csv.m
% 
% List of StatVarNames: see Table 12 of GSFLOW manual and 
%  create_table_gsflowcsv.m


clear all, close all, fclose all;

% -- data files
dir0 = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/';
dir1 = 'outputs/PRMS';
gsflow_csv_fil = 'gsflow.csv';

scenario_dir_all = cell(2,1);
scen_name = cell(3,1);
clr = zeros(3,3);
% Melt subdirectories (within dir0)
scenario_dir_all{1} =  'run5yr_melt_161209b';
scen_name{1} = 'Melt';
clr(1,:) = [0 0 1];
% No Melt
scenario_dir_all{2} =  'run5yr_NoMelt_161209c';
scen_name{2} = 'No Melt';
clr(2,:) = [1 0 0];
% No Melt, Veg shift, Temp incr
scenario_dir_all{3} =  'run5yr_NoMelt_VegTemp_161209d';
scen_name{3} = 'No Melt + Veg Shift';
clr(3,:) = [0 0 0];

filpre = 'scenarios';

gwdischarge = 1.2714e+04; % from gsflow.out, constant for all scenarios so just show it once [m3/d]

% -- variables to plot:
NSubplot = 100; % maximum
PlotVar = cell(NSubplot,1); % number of different subplots (can have mutiple variables on same subplot)
units = cell(NSubplot,1);
ii = 0;
ii = ii+1; PlotVar{ii} = {'basinstrmflow'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Stream Discharge';
ii = ii+1; PlotVar{ii} = {'basinsroff','basininterflow'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Surface Runoff + Interflow to Streams';
ii = ii+1; PlotVar{ii} = {'gwflow2strms'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Groundwater Discharge to Streams';
ii = ii+1; PlotVar{ii} = {'uzf_recharge' }; units{ii} = 'm^3/d'; plot_ti{ii} = 'Recharge';
% ii = ii+1; PlotVar{ii} = {'stream_leakage'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Stream Leakage';
% ii = ii+1; PlotVar{ii} = {'basinactet'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Evapotranspiration';
% ii = ii+1; PlotVar{ii} = {'gw_inout'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Net Groundwater In';

PlotVar = PlotVar(1:ii);

% -- plot configurations
n_prow = 4;
n_pcol = 1;

% ymd_lim = [2016+28 1 1; 2016+28 12 30];
ymd_lim = [2019 6 16; 2020 6 15];
% ymd_lim = ''; % if empty, plot all data
tick_incr_mon = 3; % tick every this many months


% loop thru different scenarios
for loop = 1: length(scenario_dir_all)
    scenario_dir = scenario_dir_all{loop};
    gsflow_csv_fil0 = [dir0, '/', scenario_dir, '/', dir1, '/', gsflow_csv_fil];


    %% Read in data

    % header{1,NVars}: variable name
    % data{NVars}: all data 
    fid = fopen(gsflow_csv_fil0);
    line0 = fgets(fid);
    D = textscan(line0,'%s', 'Delimiter', ',');
    header = D{1}';

    fmt = ['%s ', repmat('%f ', 1, length(header)-1), '\n'];
    data = textscan(fid, fmt, 'Delimiter', ',');
    fclose(fid);

    datenum_all = datenum(data{strcmp(header, 'Date')});
    [yr, mon, day] = datevec(datenum_all);

    %% plot data
    for ii = 1: length(PlotVar)
        x = nan(length(datenum_all), length(PlotVar{ii}));
        for jj = 1: length(PlotVar{ii})
            x(:,jj) = data{strcmp(header, PlotVar{ii}{jj})};
        end

        fig_i = ceil(ii/(n_prow*n_pcol));
        p_i = mod(ii-1,n_prow*n_pcol)+1;
        figure(fig_i), orient tall
        subplot(n_prow, n_pcol, p_i),
        hi=plot(datenum_all, sum(x,2), 'LineWidth', 2, 'Color', clr(loop,:)); hold on,
        if p_i == 1, h(loop)= hi; end
        title(plot_ti{ii});
        if loop == length(scen_name) && p_i == 1, legend(h,scen_name); end
        ylabel(['[', units{ii}, ']']);
        if strcmp(PlotVar{ii}, 'basinstrmflow')
            line(datenum_all([1 end]), gwdischarge*ones(1,2), 'LineWidth', 2, 'Color', ones(1,3)*0.5, 'LineStyle', '--');
        end

        % get tick marks
    %     ind_Jan1 = find(mon == 1 & day == 1);
        ind_mon1 = find(day == 1);
        ind0 = ind_mon1;
        ticks = ind0([1:tick_incr_mon:length(ind0)]);
        ticklab = cell(length(ticks),1);
        ctr = 0;
        for jj = ticks(:)'
            ctr = ctr+1;
    %             ticklab{ctr} = [num2str(mon(jj)), '/', num2str(day(jj)), '/', num2str(yr(jj))]; 
            % don't include yr
            ticklab{ctr} = [num2str(mon(jj)), '/', num2str(day(jj))]; 
        end
        set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
        if ~isempty(ymd_lim)
            datenum_lim = datenum(ymd_lim);
            set(gca, 'XLim', datenum_lim);
        end
    end
    clear datenum_all data
end

for ii = 1: fig_i
    fil = [filpre, num2str(ii), '.tiff'];
    print('-dtiff', ['-f', num2str(ii)], fil);
    fil = [filpre, num2str(ii), '.svg'];
    print('-dsvg', ['-f', num2str(ii)], fil);
end
