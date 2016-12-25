% write_lpf_MOD
% 11/17/16

function write_lpf_MOD2_f(GSFLOW_dir, infile_pre, surfz_fil, NLAY)

% clear all, close all, fclose all;

% - write to this file
% GSFLOW_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
% lpf_file = 'test.lpf';
lpf_file = [infile_pre, '.lpf'];
slashstr = '/';

% - domain dimensions, maybe already in surfz_fil and botm_fil{}?
% NLAY = 2;
% surfz_fil = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/topo.asc';
fid = fopen(surfz_fil, 'r');
D = textscan(fid, '%s %f', 6); 
NSEW = D{2}(1:4);
NROW = D{2}(5);
NCOL = D{2}(6);
fclose(fid);

% -- Base hydcond, Ss (all layers), and Sy (top layer only) on data from files
% (temp place-holder)
hydcond = ones(NROW,NCOL,NLAY)*4; % m/d (Sagehen: 0.026 to 0.39 m/d, lower K under ridges for volcanic rocks)
Ss = ones(NROW,NCOL,NLAY)* 2e-6; % constant 2e-6 /m for Sagehen
Sy = ones(NROW,NCOL,NLAY)*0.15; % 0.08-0.15 in Sagehen (lower Sy under ridges for volcanic rocks)
WETDRY = Sy; % = Sy in Sagehen (lower Sy under ridges for volcanic rocks)

% -- assumed input values
flow_filunit = 34; % make sure this matches namefile!!
hdry = 1e30;  % head assigned to dry cells
nplpf = 0;    % number of LPF parameters (if >0, key words would follow)
laytyp = zeros(NLAY,1); laytyp(1) = 1;  % flag, top>0: "covertible", rest=0: "confined"
layave = zeros(NLAY,1);  % flag, layave=1: harmonic mean for interblock transmissivity
chani = ones(NLAY,1);   % flag, chani=1: constant horiz anisotropy mult factor (for each layer)
layvka = zeros(NLAY,1);  % flag, layvka=0: vka is vert K; >0 is vertK/horK ratio
VKA = hydcond;
laywet = zeros(NLAY,1); laywet(1)=1;  % flag, 1: wetting on for top convertible cells, 0: off for confined
fl_Tr = 1; % flag, 1 for at least 1 transient stress period (for Ss and Sy)
WETFCT = 1.001; % 1.001 for Sagehen, wetting (convert dry cells to wet)
IWETIT = 4; % number itermations for wetting 
IHDWET = 0; % wetting scheme, 0: equation 5-32A is used: h = BOT + WETFCT (hn - BOT)

%% ------------------------------------------------------------------------


fmt1 = repmat('%2d ', 1, NLAY);

fil_lpf_0 = [GSFLOW_dir, slashstr, lpf_file];
fid = fopen(fil_lpf_0, 'wt');
fprintf(fid, '# LPF package inputs\n');
fprintf(fid, '%d %g %d    ILPFCB,HDRY,NPLPF\n', flow_filunit, hdry, nplpf);
fprintf(fid, [fmt1, '     LAYTYP\n'], laytyp);
fprintf(fid, [fmt1, '     LAYAVE\n'], layave);
fprintf(fid, [fmt1, '     CHANI \n'], chani);
fprintf(fid, [fmt1, '     LAYVKA\n'], layvka);
fprintf(fid, [fmt1, '     LAYWET\n'], laywet);
if ~isempty(find(laywet,1))
    fprintf(fid, '%g %d %d       WETFCT, IWETIT, IHDWET\n', WETFCT, IWETIT, IHDWET);
end

% -- Write HKSAT and Ss, Sy (if Tr) in .lpf file
format0 = [repmat(' %4.2f ', 1, NCOL), '\n'];
format1 = [repmat(' %4.2e ', 1, NCOL), '\n'];
% loop thru layers (different entry for each layer)
for lay = 1: NLAY
    fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            HY layer  %d\n', lay); % horizontal hyd cond
    fprintf(fid, format0, hydcond(:,:,lay)');

    fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            VKA layer  %d\n', lay); % vertical hyd cond
    fprintf(fid, format0, VKA(:,:,lay)');
    
    if fl_Tr
        fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            Ss layer  %d\n', lay);
        fprintf(fid, format1, Ss(:,:,lay)');
        if laytyp(lay) > 0 % convertible, i.e. unconfined
            fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            Sy layer  %d\n', lay);
            fprintf(fid, format1, Sy(:,:,lay)');
            if laywet(lay) > 0
                fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            WETDRY layer  %d\n', lay);
                fprintf(fid, format0, WETDRY(:,:,lay)');
            end
        end
    end
end
fprintf(fid, '\n');
fclose(fid);

