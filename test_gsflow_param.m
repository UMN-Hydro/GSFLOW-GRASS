% Example:
%
% suppose you want to write the following data, which are actually
% "gvr_hru_pct" and "gvr_cell_id" parameters, respectively, to your output
% file:
%
% data_A = [ 0.0265904068947, 0.0265904068947, 0.0265904068947, 
%            0.0265904068947, 0.0207706224173, 0.0115253571421,
%            0.0265904068947, 0.0237808562815 ]
%
%
% data_B = [ 64, 65, 144, 145, 146, 146, 225, 226 ];
% 
%


data_A = [ 0.0265904068947, 0.0265904068947, 0.0265904068947, ...
    0.0265904068947, 0.0207706224173, 0.0115253571421, ...
    0.0265904068947, 0.0237808562815 ];

data_B = [ 64, 65, 144, 145, 146, 146, 225, 226 ];


fid = fopen('test.file', 'w');
fprintf(fid, get_gvr_entry(data_A, 'gvr_hru_pct'));
fprintf(fid, get_gvr_entry(data_B, 'gvr_cell_id'));

fclose(fid);
