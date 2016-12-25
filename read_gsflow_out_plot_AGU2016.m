% read_gsflow_out.m 
% reads water balance results from gsflow.out 
clear all
close all
fclose all;

% *****WARNING: in gsflow.out, replace D+ with E+ for exponentials 
% (:%s/D+/E+/g)

% - Reads in:
%    varnames{1,nvar}
%    MonDayYear(ntimes,3)
%    data_cumul [m3] (ntimes x nvars)
%    data_dayrate [m3/d] (ntimes x nvars) - current day, not average

% -- data files
dir0 = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/';
dir1 = 'outputs/PRMS';
gsflow_out_fil = 'gsflow_exp.out';

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

filpre = 'scenarios_gwdis';

% -- variables to plot:
NSubplot = 100; % maximum
PlotVar = cell(NSubplot,1); % number of different subplots (can have mutiple variables on same subplot)
units = cell(NSubplot,1);
ii = 0;
% ii = ii+1; PlotVar{ii} = {'IN PRECIPITATION'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Precip';
% ii = ii+1; PlotVar{ii} = {'IN GW BOUNDARY FLOW'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Groundwater Inflow';
ii = ii+1; PlotVar{ii} = {'OUT GW BOUNDARY FLOW'}; units{ii} = 'm^3/d'; plot_ti{ii} = 'Groundwater Discharge';
PlotVar = PlotVar(1:ii);

% -- plot configurations
n_prow = 4;
n_pcol = 1;

% ymd_lim = [2016+28 1 1; 2016+28 12 30];
ymd_lim = [2019 6 16; 2020 6 15];
% ymd_lim = ''; % if empty, plot all data
ticks_ymd = [2019 7 1; 2019 10 1; 2020 1 1; 2020 4 1];


% ----------
for loop = 1: length(scenario_dir_all)
    scenario_dir = scenario_dir_all{loop};
    gsflow_out_fil0 = [dir0, '/', scenario_dir, '/', dir1, '/', gsflow_out_fil];

    fid = fopen(gsflow_out_fil0, 'r');

    varnames = cell(1,14);

    ctr = 0;
    ctr2 = 0;
    alldata = zeros(0,2);
    MonDayYear = zeros(0,3);
    line0 = fgets(fid); 
    while(1)
        line0 = fgets(fid); 
        D = textscan(fid, '%s %f %f %f %s %s %s %f \n', 1);
        if isempty(D{1}), break, end
        if mod(ctr, 1000) == 0
            MonDayYear2 = zeros(size(MonDayYear,1)+1000,3);
            MonDayYear2(1:ctr,:) = MonDayYear;
            MonDayYear = MonDayYear2;
        end    
        ctr = ctr + 1;
        MonDayYear(ctr,:) = [D{2}(1), D{3}(1), D{4}(1)]; 

        D1 = textscan(fid, '%s %s %s %f %s %s %s %f %s %f \n', 1);

        for ii = 1:8, line0 = fgets(fid); end

        while(1)
            line0 = fgets(fid);
            if isfloat(line0), break, end
            if line0(1) == '1' || isempty(line0), break, end
            if length(line0)>23
                if line0(23) == '='
                    if mod(ctr2, 1000) == 0
                        alldata2 = zeros(size(alldata,1)+1000,2);
                        alldata2(1:ctr2,:) = alldata;
                        alldata = alldata2;
                    end
                    ctr2 = ctr2 + 1;
                    if ctr == 1
                        str = strtrim(line0(1:21));
                        if ctr2<= 3
                            varnames{ctr2} = ['IN ', str];
                        elseif ctr2 <= 6
                            varnames{ctr2} = ['OUT ', str];
                        else
                            varnames{ctr2} = str;
                        end
                    end
                    alldata(ctr2,1) = str2double(line0(25:41));
                    alldata(ctr2,2) = str2double(line0(68:84));
                end
            end
        end
    end
    fclose(fid);    

    MonDayYear = MonDayYear(1:ctr,:);
    alldata = alldata(1:ctr2,:);

    %    alldata(ntimes*nvar,2); (col 1 is cumulative L^3, col 2 is rate for
    %      current day time step)
    %       alldata(1:nvar,:) is for time 1

    nvar = length(varnames);
    ntimes = size(MonDayYear, 1);
    data_cumul = reshape(alldata(:,1),[nvar,ntimes])'; % [m3] (ntimes x nvars)
    data_dayrate = reshape(alldata(:,2),[nvar,ntimes]);
    
    %% -- plot results (to compare with output from gsflow.csv)
%    varnames{1,nvar}
%    MonDayYear(ntimes,3)
%    data_cumul [m3] (ntimes x nvars)
%    data_dayrate [m3/d] (ntimes x nvars) - current day, not average    
    
    datenum_all = datenum(MonDayYear(:,[3 1 2]));
    dt = datenum_all(2)-datenum_all(1);
    % average rate over each print period:
    data = [data_cumul(1,:); data_cumul(2:end,:)-data_cumul(1:end-1,:)]/dt;
    
    for ii = 1: length(PlotVar)
        x = nan(length(datenum_all), length(PlotVar{ii}));
        for jj = 1: length(PlotVar{ii})
            x(:,jj) = data(:,strcmp(varnames, PlotVar{ii}{jj}));
        end

        fig_i = ceil(ii/(n_prow*n_pcol));
        p_i = mod(ii-1,n_prow*n_pcol)+1;
        figure(fig_i), orient tall
        subplot(n_prow, n_pcol, p_i),
        plot(datenum_all, sum(x,2), 'LineWidth', 2, 'Color', clr(loop,:)); hold on,
%         plot(datenum_all, sum(x,2)); hold on,
        title(plot_ti{ii});
        if loop == length(scen_name) && p_i == 1, legend(scen_name); end
        ylabel(['[', units{ii}, ']']);
        set(gca,'XLim', datenum_all([1 end]));

        % get tick marks
        ticks = datenum(ticks_ymd);
        ticklab = cell(length(ticks),1);
        for jj = 1:length(ticks)
            % don't include yr
            ticklab{jj} = [num2str(ticks_ymd(jj,2)), '/', num2str(ticks_ymd(jj,3))]; 
        end
        
        set(gca,'XTick', ticks, 'XTickLabel', ticklab);
           
        if ~isempty(ymd_lim)
            datenum_lim = datenum(ymd_lim);
            set(gca, 'XLim', datenum_lim);
        end
        
    end    
    
end
return
for ii = 1: fig_i
    fil = [filpre, num2str(ii), '.tiff'];
    print('-dtiff', ['-f', num2str(ii)], fil);
    fil = [filpre, num2str(ii), '.svg'];
    print('-dsvg', ['-f', num2str(ii)], fil);
end