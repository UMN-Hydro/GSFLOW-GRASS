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
clear all
fclose all;
close all

% *** CUSTOMIZE TO YOUR COMPUTER! *****************************************
% directory with files to be read in to generate PRMS input files in_GISdata_dir
in_data_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/PRMS/';
in_climatedata_dir = strcat(in_data_dir, 'climate/'); % specifically climate data

% *************************************************************************
tick_incr_mon = 3; % tick every this many months
FS = 18;

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
ymdhms_v = double([D{1} D{2} D{3} D{4} D{5} D{6}]);
tmax = D{7}; % must be F
tmin = D{8}; % must be F
precip = D{9}; % must be inch
% swrad = D{end}*1e6/41840; % MJ/m2 -> Langeley (1 Langely = 41840 J/m2)

datenum_all = datenum(ymdhms_v);

% - plot time series

% daily:
figure(1)
% get tick marks
%     ind_Jan1 = find(mon == 1 & day == 1);
ind_mon1 = find(ymdhms_v(:,3) == 1);
ind0 = ind_mon1;
ticks = ind0([1:tick_incr_mon:length(ind0)]);
ticklab = cell(length(ticks),1);
ctr = 0;
for jj = ticks(:)'
    ctr = ctr+1;
    ticklab{ctr} = [num2str(ymdhms_v(jj,2)), '/', num2str(ymdhms_v(jj,3)), '/', num2str(ymdhms_v(jj,1))]; 
end

subplot(4,1,1);
plot(datenum_all, precip*25.4), ylabel('precip [mm/d]');
set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
axis tight
% if ~isempty(ymd_lim)
%     datenum_lim = datenum(ymd_lim);
%     set(gca, 'XLim', datenum_lim);
% end

subplot(4,1,2);
plot(datenum_all, ((tmin+tmax)/2-32)*5/9), ylabel('T [^\circC]');
set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
axis tight

subplot(4,1,3);
plot(datenum_all, (tmin-32)*5/9), ylabel('Tmin [^\circC]');
set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
axis tight

subplot(4,1,4);
plot(datenum_all, (tmax-32)*5/9), ylabel('Tmax [^\circC]');
set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
axis tight

% daily final:
figure(3)
% get tick marks
%     ind_Jan1 = find(mon == 1 & day == 1);
ind_mon1 = find(ymdhms_v(:,3) == 1);
ind0 = ind_mon1;
ticks = ind0([1:tick_incr_mon:length(ind0)]);
ticklab = cell(length(ticks),1);
ctr = 0;
for jj = ticks(:)'
    ctr = ctr+1;
    ticklab{ctr} = [num2str(ymdhms_v(jj,2)), '/', num2str(ymdhms_v(jj,3)), '/', num2str(ymdhms_v(jj,1))]; 
end

subplot(2,1,1);
plot(datenum_all, precip*25.4, 'LineWidth', 2), ylabel('precip [mm/d]');
set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
axis tight
% if ~isempty(ymd_lim)
%     datenum_lim = datenum(ymd_lim);
%     set(gca, 'XLim', datenum_lim);
% end

subplot(2,1,2);
plot(datenum_all, ((tmin+tmax)/2-32)*5/9, 'LineWidth', 2), ylabel('T [^\circC]');
set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
axis tight

% subplot(4,1,3);
% plot(datenum_all, (tmin-32)*5/9, 'LineWidth', 2), ylabel('Tmin [^\circC]');
% set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
% axis tight
% 
% subplot(4,1,4);
% plot(datenum_all, (tmax-32)*5/9, 'LineWidth', 2), ylabel('Tmax [^\circC]');
% set(gca,'XTick', datenum_all(ticks), 'XTickLabel', ticklab);
% axis tight
print('-dtiff', '-f3', 'd_precip_temp_BT.tiff');
print('-dsvg', '-f3', 'd_precip_temp_BT.svg');


% monthly:
precip_mon = zeros(12,1); tmin_mon = zeros(12,1); tmax_mon = zeros(12,1);
for m = 1:12
    precip_mon(m) = mean(precip(ymdhms_v(:,2) == m));
    tmin_mon(m) = mean(tmin(ymdhms_v(:,2) == m));
    tmax_mon(m) = mean(tmax(ymdhms_v(:,2) == m));
end

figure(2)
for jj = ticks(:)'
    ctr = ctr+1;
    ticklab{ctr} = [num2str(ymdhms_v(jj,2)), '/', num2str(ymdhms_v(jj,3)), '/', num2str(ymdhms_v(jj,1))]; 
end

m_order = [6:12, 1:5];
mticklab = cell(12,1);
for ii = 1:12, mticklab{ii} = num2str(m_order(ii)); end

subplot(2,1,1);
plot([1:12], precip_mon(m_order)*25.4, '-*', 'LineWidth', 2), ylabel('precip [mm/d]');
set(gca, 'XTick', [1:12], 'XLim', [0.5 12.5], 'XTickLabel', mticklab);
hold on,
plot([0.5 12.5], ones(1,2)*mean(precip)*25.4, '--')
xlabel('month')
set(gca,'FontSize', FS)

subplot(2,1,2);
plot([1:12], ((tmin_mon(m_order)+tmax_mon(m_order))/2-32)*5/9, '-*', 'LineWidth', 2), ylabel('T [^{\circ}C]');
set(gca, 'XTick', [1:12], 'XLim', [0.5 12.5], 'XTickLabel', mticklab);
hold on
plot([0.5 12.5], ones(1,2)*mean(((tmin_mon+tmax_mon)/2-32)*5/9), '--')
xlabel('month')
set(gca,'FontSize', FS)

print('-dtiff', '-f2', 'mon_precip_temp_BT.tiff');
print('-dsvg', '-f2', 'mon_precip_temp_BT.svg');

% subplot(4,1,3);
% plot([1:12], (tmin_mon-32)*5/9, '-*', 'LineWidth', 2), ylabel('Tmin C');
% set(gca, 'XTick', [1:12], 'XLim', [0.5 12.5]);
% 
% subplot(4,1,4);
% plot([1:12], (tmax_mon-32)*5/9, '-*', 'LineWidth', 2), ylabel('Tmax C');
% set(gca, 'XTick', [1:12], 'XLim', [0.5 12.5]);
