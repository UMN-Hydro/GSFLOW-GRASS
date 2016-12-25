% test_map_gwcells.m
%
% test different spatial distributions of things
fclose all;
clear all
close all

veg_thresh = 4400; % m (Rachel's email 12/8/16 8:33am), approx veg elev line
veg_thresh_shift = veg_thresh + 500; % Morueta-Holme et al. 2015: shift >500m since 1802

in_data_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/';
in_GISdata_dir = strcat(in_data_dir, 'GIS/'); % specifically GIS data
HRUfil = strcat(in_GISdata_dir, 'HRU.csv');

% load elevations of hru's
HRUdata = importdata(HRUfil,',',1);
hru_elev = HRUdata.data(:, strcmp(HRUdata.colheaders, 'hru_elev'));  % elev_units = m

ind_veg = find(hru_elev > veg_thresh); 
ind_veg_shift = find(hru_elev > veg_thresh_shift); 

[~,hrui_melt] = max(hru_elev);


