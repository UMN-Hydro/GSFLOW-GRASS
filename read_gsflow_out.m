% read_gsflow_out.m
clear all
close all
fclose all;

% - Reads in:
%    varnames{1,nvar}
%    MonDayYear(ntimes,3)
%    data_cumul [m3] (ntimes x nvars)
%    data_dayrate [m3/d] (ntimes x nvars) - current day, not average

gsflow_out_fil = '/home/gcng/Shortcuts/AndesWaterResources/GSFLOW/simdir/run5yr_melt_161209b/outputs/PRMS/gsflow.out';

fid = fopen(gsflow_out_fil, 'r');

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
