
% -------------------------------------------------------------------------
% You need the following inputs

NROW = 15;
NCOL = 10;
NLAY = 1;
NPER = 12;

%Item 1:
NUZTOP = 1;     %Which cell in a vertical column for which recharge&discharge is simulated.
IUZFOPT = 1;    %Vertical hydraulic conductivity is specified using VKS
IRUNFLG = 1;    %Groundwater discharged to stream or lake
IETFLG = 1;     %Evaporation will be simulated or not.
IUZFCB1 = 0;    %Flag for writing rates of groundwater recharge, ET&groundwater discharge.0:wont be written
IUZFCB2 = 61;   %Writing groundwater recharge, discharge&ET in UBUDSV3 format
NTRAIL2 = 25;   %Number of traling waves
NSETS2 = 20;    %Number of wave sets to simulate multiple infilteration periods=20.
NUZGAG = 4;
SURFDEP = 1.0;

%Item 2-7:
project_name = 'TestProject';   % used to name the output file (.uzf)
iuzfbnd = importdata('./data/iuzfbnd.dat');
irunbnd = importdata('./data/irunbnd.dat');
vks = importdata('./data/vks.dat');
eps = 3.5;  %Brooks-Corey epsilon of the unsaturated zone.
thts = 0.35;    %Sturated water content of the unsaturated zone
thti = 0.0;     %initial water content for each vertical column of cells-not specified for steady-state simulations
uzgag = importdata('./data/uzgag.dat');

finf = zeros(NROW, NCOL, NPER);
A = importdata('./data/finf.dat');
B = reshape(A', NCOL, NROW, NPER);
for i=1:size(B, 3)
    finf(:, :, i) = transpose(B(:, :, i));
end

pet = 5.0E-08;  %array of ET demands rate (L/T)
extdp = 15.0;   %array of ET extiction zone~altitude of the soil-zone base;specified at least for 1st stress period
extwc = importdata('./data/extwc.dat'); %array of Extinction water content

% -------------------------------------------------------------------------

% Ouput file
fname = [project_name, '.uzf'];
fid = fopen(fname, 'w');


% Write header lines (item 0)
heading = '# Unsaturated-Zone Flow (UZF) input file.\n';
fprintf(fid, heading);
fprintf(fid, '# %s simulation -- %s.\n', upper(project_name), date);


% Write item 1
fmtarr_int = [repmat('  %6d', 1, 9), '  %10.6E'];
fprintf(fid, fmtarr_int, NUZTOP, IUZFOPT, IRUNFLG, IETFLG, IUZFCB1, IUZFCB2, ...
    NTRAIL2, NSETS2, NUZGAG, SURFDEP);

comment = ['     NUZTOP  IUZFOPT  IRUNFLG  IETFLG  IUZFCB1  IUZFCB2  '...
    'NTRAIL2  NSETS2  NUZGAG  SURFDEP\n'];
fprintf(fid, comment);


% Write item2 [IUZFBND (NCOL, NROW)] - U2DINT
ent = get_file_entry(iuzfbnd, 'U2DINT', 1, ...
    '#IUZFBND--AREAL EXTENT OF THE ACTIVE MODEL');
fprintf(fid, ent);


% write item 3 [IRUNBND (NCOL, NROW)] - U2DINT
if (IRUNFLG > 0)
    ent = get_file_entry(irunbnd, 'U2DINT', 1, ...
        '#IRUNBND--STREAM SEGMENTS OR LAKE NUMBERS');
    fprintf(fid, ent);
end


% write item 4 [VKS (NCOL, NROW)] - U2DREL
if (IUZFOPT == 1)
    ent = get_file_entry(vks, 'U2DREL', 1.0, ...
        '#VKS--VERTICAL HYDRAULIC CONDUCTIVITY OF THE UNSATURATED ZONE');
    fprintf(fid, ent);
end


% write items 5, 6, 7
fprintf(fid, get_file_entry(eps, 'U2DREL', 1.0, ...
    '#EPS--BROOKS/COREY EPSILON'));

fprintf(fid, get_file_entry(thts, 'U2DREL', 1.0, ...
    '#THTS--SATURATED WATER CONTENT'));

fprintf(fid, get_file_entry(thti, 'U2DREL', 1.0, ...
    '#THTI--INITIAL WATER CONTENT'));


% write item 8
if (NUZGAG > 0)
    comment = '#IUZROW IUZCOL IFTUNIT IUZOPT--UNSATURATED FLOW OUTPUT';
    for i=1:NUZGAG
        dummyrow = uzgag(i, :);
        if ~isnan(dummyrow)
            fprintf(fid, '  %6d  %6d  %6d  %6d      %s\n', dummyrow, comment);
        else
            ind = ~isnan(dummyrow);
            fprintf(fid, ['  %6d', repmat(' ', 1, 30), '%s\n'], dummyrow(ind), comment);
        end
    end
end


% write items 9-16
for iper=1:NPER
    
    % write item 9
    comment = sprintf('#NUZF1 FOR STRESS PERIOD %d\n', iper);
    if (iper <= size(finf, 3))
        NUZF1 = 1;
    else
        NUZF1 = -1;
    end
    
    fprintf(fid, '  %6d      %s', NUZF1, comment);
    
    % write item 10
    if (NUZF1 > 0)
        ent = get_file_entry(transpose(finf(:, :, iper)), 'U2DREL', ...
            1.0, sprintf('#FINF--STRESS PERIOD %d', iper));
        
        fprintf(fid, ent);
    end
    
    % ---------------------------------------------------------------------
    
    % items 11-16
    if (IETFLG > 0)
        
        % write items 11
        comment = sprintf('#NUZF2 FOR STRESS PERIOD %d\n', iper);
        if (iper <= size(pet, 3))
            NUZF2 = 1;
        else
            NUZF2 = -1;
        end
        
        fprintf(fid, '  %6d      %s', NUZF2, comment);
        
        % write item 12
        if (NUZF2 > 0)
            ent = get_file_entry(transpose(pet(:, :, iper)), 'U2DREL', ...
                1.0, sprintf('#PET--STRESS PERIOD %d', iper));
            
            fprintf(fid, ent);
        end
        
        % -----------------------------------------------------------------        
        
        % write items 13
        comment = sprintf('#NUZF3 FOR STRESS PERIOD %d\n', iper);
        if (iper <= size(extdp, 3))
            NUZF3 = 1;
        else
            NUZF3 = -1;
        end
        
        fprintf(fid, '  %6d      %s', NUZF3, comment);
        
        % write item 14
        if (NUZF3 > 0)
            ent = get_file_entry(transpose(extdp(:, :, iper)), 'U2DREL', ...
                1.0, sprintf('#EXTDP--STRESS PERIOD %d', iper));
            
            fprintf(fid, ent);
        end
        
        % -----------------------------------------------------------------
        
        % write items 15
        comment = sprintf('#NUZF4 FOR STRESS PERIOD %d\n', iper);
        if (iper <= size(extwc, 3))
            NUZF4 = 1;
        else
            NUZF4 = -1;
        end
        
        fprintf(fid, '  %6d      %s', NUZF4, comment);
        
        % write item 16
        if (NUZF4 > 0)
            ent = get_file_entry(transpose(extwc(:, :, iper)), 'U2DREL', ...
                1.0, sprintf('#EXTWC--STRESS PERIOD %d', iper));
            
            fprintf(fid, ent);
        end
        
    end
    
end

fclose(fid);

% -------------------------------------------------------------------------
% End of the script
