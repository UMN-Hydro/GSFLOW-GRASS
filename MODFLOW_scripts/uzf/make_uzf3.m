clear all, close all, fclose all;

% -------------------------------------------------------------------------
% You need the following inputs

% - write to this file
GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
uz_file = 'test.uzf';
slashstr = '/';

surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
fid = fopen(surfz_fil, 'r');
D = textscan(fid, '%s %f', 6); 
NSEW = D{2}(1:4);
NROW = D{2}(5);
NCOL = D{2}(6);
fclose(fid);
NLAY = 1;
NPER = 2;
% **** ASSUMES PER 1 IS SS, PER 2 IS TR ****

%Item 1:
%NUZTOP: Which cell (layer) has recharge&discharge simulated in it; 
%   1:top layer has recharge/discharge and land-surface (at top)
%   2:layer specificed in UZFBND, top of that layer is land-surface
%   3:highest active layer, top of UZFBND is land-surface (*sagehen example, with watershed mask of UZFBND=1)
NUZTOP = 3;     
IUZFOPT = 1;    % 1: Vertical hydraulic conductivity is specified using VKS, 2: use VKA from LPF file
IRUNFLG = 0;    % >0: Groundwater discharged to land-surface is routed to streams or lakes (using IRUNBND) (*sagehen=0, ONLY avaliable in GSFLOW for SS)
IETFLG = 0;     % ~=0: Evaporation will be simulated. (*sagehen=0; GSFLOW: =0 for PRMS ET only, =1 for additional ET from below soil zone)
IUZFCB1 = 0;    %Flag for writing rates of groundwater recharge, ET&groundwater discharge in UBUDSV format.0:wont be written, >0: file unit number
IUZFCB2 = 61;   %Writing groundwater recharge, discharge&ET in UBUDSV3 format; >0: file unit number
NTRAIL2 = 25;   %Number of trailing waves to define theta profile, 10 to 20 usually ok, higher for better mass balance
NSETS2 = 20;    %Number of wave sets to simulate multiple infiltration periods, 20 usually ok.
NUZGAG = 0;     % number of cells for which to print detailed info on UZ water budget and theta (see uzgag below)
SURFDEP = 1.0;  % average undulation depth within finite diff cell (?)

%Item 2-7:
project_name = 'TestProject';   % used to name the output file (.uzf)
iuzfbnd = ones(NROW,NCOL); % [NROW,NCOL] layer w/ top as land-surface and/or w/ discharge/recharge (see NUZTOP), default: mask with 1's
if IRUNFLG > 0
    irunbnd = importdata('./data/irunbnd.dat'); % [NROW,NCOL] only for IRUNFLG>0, stream seg to which gw discharge is routed
end
% vks = importdata('./data/vks.dat'); % [NROW,NCOL] saturated K, no needed if using value in LPF (IUZFOPT=2), [m/d]
vks = 4*ones(NROW,NCOL); 
eps = 3.5;  %Brooks-Corey epsilon of the unsaturated zone.
thts = 0.35;    %Saturated water content of the unsaturated zone
thti = 0.0;     %initial water content for each vertical column of cells-not specified for steady-state simulations
if NUZGAG > 0
    uzgag = importdata('./data/uzgag.dat'); % only for NUZGAG>0; row, col, file unit for output, IUZOPT (flag for what data to output)
end

% - infiltration (in general, not needed  bc provided by PRMS, but option to apply for initial SS period for MODFLOW)
NUZF1 = [1; -1*ones(NPER-1,1)]; % infiltration can ONLY specified for (initial) SS stress periods (o.w. calculated by PRMS)
% A = importdata('./data/finf.dat'); % infiltration rate [NROW,NCOL,1], **max ONLY for 1 stress period: first SS stress period! [m/d]
% finf = A';
% B = reshape(A', NCOL, NROW, NPER);
% for i=1:size(B, 3)
%     finf(:, :, i) = transpose(B(:, :, i));
% end
finf = ones(NROW,NCOL);
finf(:,:,1) = ones(NROW,NCOL)*0.00; % m/d (1.1m/d typical in Sagehen)

% - ET (in general, not needed bc provided by PRMS, but option to apply excess ET below PRMS' soil-zone, i.e. apply to MODFLOW's UZ)
% pet = 5.0E-08;  %array of ET demands rate (L/T), PET not used in GSFLOW
% for IETFLG>0, specify how ET can be drawn from below soil-zone base;
%   assume same for all stress periods
if IETFLG>0
    NUZF2 = [-1*ones(NPER,1)]; % use ET from below soil-zone
    NUZF3 = [1; -1*ones(NPER-1,1)]; % only specify extdp for first stress periods
    extdp = 15.0*ones(NROW,NCOL);   %array of ET extiction zone~altitude of the soil-zone base;specified at least for 1st stress period; only for IETFLG>0
    NUZF4 = [1; -1*ones(NPER-1,1)]; % only specify extwc for first stress periods
    extwc = thts*0.9*ones(NROW,NCOL); %array of Extinction water content; EXTWC must be between (THTS-Sy) and THTS; only for IETFLG>0 and 
end
% -------------------------------------------------------------------------

% Ouput file
fname = [GSFLOW_dir, '/', uz_file];
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

% fprintf(fid, get_file_entry(thti, 'U2DREL', 1.0, ...
%     '#THTI--INITIAL WATER CONTENT'));


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
%         ent = get_file_entry(transpose(finf(:, :, iper)), 'U2DREL', ...
%             1.0, sprintf('#FINF--STRESS PERIOD %d', iper));
        ent = get_file_entry(finf(:, :, iper), 'U2DREL', ...
            1.0, sprintf('#FINF--STRESS PERIOD %d', iper));
        
        fprintf(fid, ent);
    end
    
    % ---------------------------------------------------------------------
    
    % items 11-16
    if (IETFLG > 0)
        
        % write items 11
        comment = sprintf('#NUZF2 FOR STRESS PERIOD (GSFLOW DOES NOT USE PET) %d\n', iper);
        fprintf(fid, '  %6d      %s', NUZF2(iper), comment);
        
        % write item 12 (GSFLOW DOES NOT USE PET)
        if (NUZF2(iper) > 0)
        end
        
        % -----------------------------------------------------------------        
        
        % write items 13
        comment = sprintf('#NUZF3 FOR STRESS PERIOD %d\n', iper);
        fprintf(fid, '  %6d      %s', NUZF3(iper), comment);
        
        % write item 14
        format0 = [repmat(' %4.2f ', 1, NCOL), '\n'];
        if (NUZF3(iper) > 0)
            fprintf(fid, 'INTERNAL   1. (FREE)    0           EXTDP FOR STRESS PERIOD %d\n', iper); % horizontal hyd cond
            for irow = 1:NROW
                fprintf(fid, format0, extdp(irow,:));
            end
        end
        
        % -----------------------------------------------------------------
        
        % write items 15
        comment = sprintf('#NUZF4 FOR STRESS PERIOD %d\n', iper);
        fprintf(fid, '  %6d      %s', NUZF4(iper), comment);
        
        % write item 16
        if (NUZF4(iper) > 0)
            fprintf(fid, 'INTERNAL   1. (FREE)    0           EXTWC FOR STRESS PERIOD %d\n', iper); % horizontal hyd cond
            for irow = 1:NROW
                fprintf(fid, format0, extwc(irow,:));
            end            
        end
        
    end
    
end

fclose(fid);

% -------------------------------------------------------------------------
% End of the script
