% write_OC_PCG_MOD.m
% 11/20/16
function write_OC_PCG_MOD_f(GSFLOW_indir, infile_pre, perlen_tr)

% clear all, close all, fclose all;

% - write to this file
% GSFLOW_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/';
fil_pcg = [infile_pre, '.pcg'];
fil_oc = [infile_pre, '.oc'];
slashstr = '/';

% -- shoud match .dis
NPER = 2; % 1 SS then 1 transient
PERLEN = [1; perlen_tr];  % 2 periods: 1-day steady-state and multi-day transient
NSTP = PERLEN;

% -- pcg and oc files are not changed with this script
% fil_pcg_0 = fullfile(MODtest_dir0, fil_pcg);
fil_pcg_0 = [GSFLOW_indir, slashstr, fil_pcg];
fid = fopen(fil_pcg_0, 'wt');
% fprintf(fid, '# Preconditioned conjugate-gradient package\n');
% fprintf(fid, '        50        30         1      MXITER, ITER1, NPCOND\n');
% fprintf(fid, '  0000.001      .001        1.         2         1         1      1.00\n');
% fprintf(fid, ' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n');

% % sagehen example:
% fprintf(fid, '# Preconditioned conjugate-gradient package\n');
% fprintf(fid, '        1000    450         1      MXITER, ITER1, NPCOND\n');
% fprintf(fid, '      0.001     0.08       1.0        2         1         0      -0.05  0.70\n');
% fprintf(fid, ' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n');

% sagehen example:
fprintf(fid, '# Preconditioned conjugate-gradient package\n');
fprintf(fid, '        1000    450         1      MXITER, ITER1, NPCOND\n');
fprintf(fid, '      0.001     0.08       1.0        2         1         0      -0.05  0.70\n');
fprintf(fid, ' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n');

fclose(fid);

% fil_oc_0 = (MODtest_dir0, fil_oc);
% "PRINT": to listing file
% "SAVE": to file with unit number in name file
fil_oc_0 = [GSFLOW_indir, slashstr, fil_oc];
fid = fopen(fil_oc_0, 'wt');
fprintf(fid, 'HEAD PRINT FORMAT 20\n');
fprintf(fid, 'HEAD SAVE UNIT 51\n');
fprintf(fid, 'COMPACT BUDGET AUX\n');
fprintf(fid, 'IBOUND SAVE UNIT 52\n');
for per_i = 1:NPER
    for stp_i = 1:30:NSTP(per_i)
        fprintf(fid, 'PERIOD %d STEP %d\n', per_i, stp_i);
        if stp_i == NSTP(per_i) % only at end of stress period
            fprintf(fid, '   PRINT HEAD\n');
            fprintf(fid, '   SAVE IBOUND\n');
            fprintf(fid, '   PRINT BUDGET\n');
        end
        fprintf(fid, '   SAVE HEAD\n');
        fprintf(fid, '   SAVE BUDGET\n');
    end
end
fclose(fid);
