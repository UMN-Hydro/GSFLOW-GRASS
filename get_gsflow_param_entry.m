function [ output_str ] = get_gvr_entry(parvalues, parname)

% Table A1-23
if strcmp(parname, 'model_mode')
    dimention_names = {'one'};
    partype = 1;
    
elseif strcmp(parname, 'modflow_name')
    dimention_names = {'one'};
    partype = 4;
    
elseif strcmp(parname, 'mxsziter')
    dimention_names = {'one'};
    partype = 1;

% Table A1-24
elseif strcmp(parname, 'gvr_cell_id')
    dimention_names = {'nhrucell'};
    partype = 1;
    
elseif strcmp(parname, 'gvr_cell_pct')
    dimention_names = {'nhrucell'};
    partype = 2;

% Table A1-25
elseif strcmp(parname, 'gvr_hru_id')
    dimention_names = {'nhrucell'};
    partype = 1;
    
elseif strcmp(parname, 'gvr_hru_pct')
    dimention_names = {'nhrucell'};
    partype = 2;
    
elseif strcmp(parname, 'hru_area')
    dimention_names = {'nhru'};
    partype = 2;
    
elseif strcmp(parname, 'hru_type')
    dimention_names = {'nhru'};
    partype = 1;
    
elseif strcmp(parname, 'lake_hru_id')
    dimention_names = {'nhru'};
    partype = 1;
    
elseif strcmp(parname, 'hru_segment')
    dimention_names = {'nhru'};
    partype = 1;
    
elseif strcmp(parname, 'local_reachid')
    dimention_names = {nreach};
    partype = 1;
    
elseif strcmp(parname, 'numreach_segment')
    dimention_names = {'nsegment'};
    partype = 1;
    
elseif strcmp(parname, 'reach_segment')
    dimention_names = {'nreach'};
    partype = 1;
    
elseif strcmp(parname, 'segment_pct_area')
    dimention_names = {'nreach'};
    partype = 2;
    
elseif strcmp(parname, 'szconverge')
    dimention_names = {'one'};
    partype = 2;

% Table A1-26
elseif strcmp(parname, 'basin_cfs_init')
    dimention_names = {'one'};
    partype = 2;

% Table A1-27
elseif strcmp(parname, 'csv_output_file')
    dimention_names = {'one'};
    partype = 4;

elseif strcmp(parname, 'gsf_rpt')
    dimention_names = {'one'};
    partype = 1;

elseif strcmp(parname, 'gsflow_output_file')
    dimention_names = {'one'};
    partype = 4;
    
elseif strcmp(parname, 'id_obsrunoff')
    dimention_names = {'one'};
    partype = 1;
    
elseif strcmp(parname, 'model_output_file')
    dimention_names = {'one'};
    partype = 4;
    
elseif strcmp(parname, 'rpt_days')
    dimention_names = {'one'};
    partype = 4;
    
elseif strcmp(parname, 'runoff_units')
    dimention_names = {'one'};
    partype = 1;
end


ndimension = length(dimention_names);

dimention_names_str = '';
for i=1:length(dimention_names)
    dimention_names_str = strcat(dimention_names_str, dimention_names{i}, '\n');
end

if (partype == 4)
    nvalues = 1;
    values = strcat(parvalues, '\n');
else
    nvalues = length(parvalues);
    
    if (partype == 1)
        fmt = '%d';
    else
        fmt = '%f';
    end
    
    values = '';
    for i=1:nvalues
        vstr = sprintf(fmt, parvalues(i));
        values = strcat(values, vstr, '\n');
    end
end


output_str = strcat('####\n', parname, '\n', num2str(ndimension), '\n', ...
    dimention_names_str, num2str(nvalues), '\n', num2str(partype), '\n', values);

end