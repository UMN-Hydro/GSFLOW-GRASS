
fid = fopen(filename, 'w');
fprintf(fid, 'This is empty  \n');  % line 1: comment
fprintf(fid, 'tmax 0  \n');  % line 1: comment
fprintf(fid, 'tmin 0  \n');  % line 1: comment
fprintf(fid, 'precip 0  \n');  % line 1: comment
fprintf(fid, 'runoff 0  \n');  % line 1: comment
fprintf(fid, '######## \n');  % 


daysinmon = [31 28 31 30 31 30 31 31 30 31 30 31];
for y = 2000:2010
    for m = 1:12
        for d = 1: daysinmon(m)
            datev = [y m d];
            fprintf(fid, '%d %d %d 0 0 0 \n', datev);
            
%             % for including floating point data:
%             fmt = '%d %d %d 0 0 0 %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n';
%             fprintf(fid, fmt, [datev, precip_data]);

            
        end
    end
end
