function y = get_file_entry(x, subroutine, const, comment)
% assumes x is [nrow, ncol]

if const == 0
    const = 1;
end

if (size(x) == [1, 1])
    if strcmp(subroutine, 'U2DREL')
        y = sprintf('CONSTANT   %10.6E   %s\n', x, comment);
   
        y = sprintf('CONSTANT   %6d   %s\n', x, comment);
    end
    
else
    ncol = size(x, 2);
    if strcmp(subroutine, 'U2DREL')
        fmtin = sprintf('%dE%d.%d', ncol, 15, 6);
        fmtarr = [repmat('  %13.6E', 1, ncol), '\n'];
        ctrl_line = 'INTERNAL  %13.6E  (%s)  -1  %s\n';
    elseif strcmp(subroutine, 'U2DINT')
        fmtin = sprintf('%dI%d', ncol, 8);
        fmtarr = [repmat('  %6d', 1, ncol), '\n'];
        ctrl_line = 'INTERNAL  %6d  (%s)  -1  %s\n';
    end
    
    cline = sprintf(ctrl_line, const, fmtin, comment);
    arrvals = sprintf(fmtarr, x');
    y = [cline, arrvals];
end