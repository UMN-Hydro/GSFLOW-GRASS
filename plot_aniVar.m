% plot_aniVar.m
% 12/7/16

clear all, close all, fclose all;

% prefix to file names, will have suffixes with dimensions (e.g., nhru)
aniVarFile_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/PRMS/'; % end with slash
aniVarFile_pre = 'ChimOut.ani';
suffixes = {'.nhru', '.nsegment'};
dim = [11 11];

%% Read in data
field_all = cell(length(suffixes),1); nspace_all = cell(length(suffixes),1);
datestr_all = cell(length(suffixes),1); data_all = cell(length(suffixes),1); 
for ii = 1: length(suffixes)
% for ii = 1: 1
    fil = [aniVarFile_dir, aniVarFile_pre, suffixes{ii}];
    fprintf('reading in %s... \n', fil);

    fid = fopen(fil);
    
    % get header info
    ctr = 0;
    field = cell(1000,1);
    nspace = nan(1000,1);
    while(1)
        line0 = fgets(fid);
        if ~strcmp(line0(1), '#')
            break,
        elseif strncmp(line0, '# Begin', 7), 
            while(1)
                line1 = fgets(fid);
                D = textscan(line1(1:end), '# %s %s %d %d', 'Delimiter', ',');
                if strncmp(D{1}, 'End', 3)
                    break
                end
                ctr = ctr + 1;
                field{ctr} = D{1};
                nspace(ctr) = D{3};
            end
        end
    end
    field = field(1:ctr);
    nspace = nspace(1:ctr);
    
    field_all{ii} = field;
    nspace_all{ii} = nspace;
    line0 = fgets(fid);

    fmt = ['%10s ', repmat('%10f ', 1, length(field)-1)];
    ctr = 0;
    datestr = cell(1000,1);
    data = nan(1000, length(field)-1);

    while(1)
        line0 = fgets(fid);
        if isempty(line0), break, end
        D = textscan(line0, fmt, 1);
        if ~(D{1}{1}(5)=='-' && D{1}{1}(8)=='-')
            line1 = line0(dim(ii)*sum(1+nspace_all{ii}(3:end))+1:end);
            if isempty(line1), break, end
            D = textscan(line1, fmt, 2);  %***** RESUME HERE!!!
        end
        if ~mod(ctr,1000) && ctr > 0
            datestr2 = cell(ctr+1000,1);
            data2 = nan(ctr+1000,length(field)-1);
            datestr2(1:ctr) = datestr;
            data2(1:ctr,:) = data;
            datestr = datestr2;
            data = data2;
        end
        ctr = ctr+1;
        datestr{ctr} = D{1}{1};
        data(ctr,:) = cell2mat(D(2:end));
    end
    fclose(fid);
    
    datestr_all{ii} = datestr; % cell(nrecs,1)
    data_all{ii} = data;  % cell(nrecs,nfields-2)
end

    
    
