% plot_gsflow_csv.m
% 
% List of StatVarNames: see Table 12 of GSFLOW manual and 
%  create_table_gsflowcsv.m


clear all, close all, fclose all;

% -- data files
dir0 = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/';
dir1 = 'outputs/PRMS';
gsflow_csv_fil = 'gsflow.csv';

% Choose one:
% Melt subdirectories (within dir0)
scenario_dir =  'run5yr_melt_161209b';
% % No Melt
% scenario_dir =  'run5yr_NoMelt_161209c';
% No Melt, Veg shift, Temp incr
% scenario_dir =  'run5yr_NoMelt_VegTemp_161209d';

gsflow_csv_fil0 = [dir0, '/', scenario_dir, '/', dir1, '/', gsflow_csv_fil];


% -- variables to plot:
NSubplot = 100; % maximum
PlotVar = cell(NSubplot,1); % number of different subplots (can have mutiple variables on same subplot)
ii = 0;
ii = ii+1; PlotVar{ii} = {'basinppt', 'basinactet'};
ii = ii+1; PlotVar{ii} = {'basinsroff', 'basininterflow', 'gwflow2strms'};
ii = ii+1; PlotVar{ii} = {'basinstrmflow'};
ii = ii+1; PlotVar{ii} = {'gw_inout'};

ii = ii+1; PlotVar{ii} = {'uzf_recharge' };
ii = ii+1; PlotVar{ii} = {'uzf_infil'};
ii = ii+1; PlotVar{ii} = {'basingw2sz', 'net_sz2gw'};
ii = ii+1; PlotVar{ii} = {'stream_leakage', 'gwflow2strms'};

PlotVar = PlotVar(1:ii);

% -- plot configurations
n_prow = 4;
n_pcol = 1;

% ymd_lim = [2016+28 1 1; 2016+28 12 30];
ymd_lim = [2019 6 16; 2020 6 15];
% ymd_lim = ''; % if empty, plot all data
tick_incr_mon = 3; % tick every this many months

filpre = scenario_dir;

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
for ii = 1:length(PlotVar)
    x = nan(length(datenum_all), length(PlotVar{ii}));
    for jj = 1: length(PlotVar{ii})
        x(:,jj) = data{strcmp(header, PlotVar{ii}{jj})};
    end
    
    fig_i = ceil(ii/(n_prow*n_pcol));
    p_i = mod(ii-1,n_prow*n_pcol)+1;
    figure(fig_i), orient tall
    subplot(n_prow, n_pcol, p_i),
    plot(datenum_all, x);
    legend(PlotVar{ii});

    % get tick marks
%     ind_Jan1 = find(mon == 1 & day == 1);
    ind_mon1 = find(day == 1);
    ind0 = ind_mon1;
    ticks = ind0([1:tick_incr_mon:length(ind0)]);
    ticklab = cell(length(ticks),1);
    ctr = 0;
    for jj = ticks(:)'
        ctr = ctr+1;
        ticklab{ctr} = [num2str(mon(jj)), '/', num2str(day(jj)), '/', num2str(yr(jj))]; 
    end
    set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
    if ~isempty(ymd_lim)
        datenum_lim = datenum(ymd_lim);
        set(gca, 'XLim', datenum_lim);
    end
end

for ii = 1: fig_i
    fil = [filpre, num2str(ii), '.tiff'];
    print('-dtiff', ['-f', num2str(ii)], fil);
end
