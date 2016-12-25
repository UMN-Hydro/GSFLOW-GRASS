% Plot_MODFLOW_head_flow_1D_general_quickstart.m
%
% 10/28/15 (can handle transient conditions)
%
% plots MODFLOW results: head image plots for each layer
%
% **** IN PROGRESS: PLOT WTD ***

clear all, fclose all;
close all,

% --------------------------------------------------------

% -- directory with simulation results, file names
% MOD_simdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/outputs/MODFLOW/';
dir0 = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/';
dir1 = 'outputs/MODFLOW';

scenario_dir_all = cell(2,1);
scen_name = cell(3,1);
% Melt subdirectories (within dir0)
scenario_dir_all{1} =  'run5yr_melt_161209b';
scen_name{1} = 'Melt';
% No Melt
scenario_dir_all{2} =  'run5yr_NoMelt_161209c';
scen_name{2} = 'No Melt';
% No Melt, Veg shift, Temp incr
scenario_dir_all{3} =  'run5yr_NoMelt_VegTemp_161209d';
scen_name{3} = 'No Melt + Veg Shift';

scen_i = 1;
MOD_simdir = [dir0, '/', scenario_dir_all{scen_i}, '/', dir1, '/'];
uzf_file = 'uzf.dat'; % head data

% MOD_simdir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/simdir/spinup30yr_constH/outputs/MODFLOW/';
% uzf_file = 'testhead.dat'; % head data

% MOD_simdir = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0_gcng/data/sagehen/output/modflow/';
% uzf_file = 'head_sagehen.out'; % hea5d data

% ibound_file = 'ibound.dat'; % needed to get active cell info for reading water budget file
% bud_file = 'test.bud'; % water budget file (for flow lines)
% lpf_file = 'test2lay.lpf'; % read in for K inputs

% to mask out 
GIS_indir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/Data/GIS/';
mask_fil = [GIS_indir, 'basinmask_dischargept.asc']; % only for Chimb

% Only ONE can be 1, others 0
fl_MOD_PC = 0;
fl_MOD_USGS_RedHat = 0; 
fl_MOD_UBUNTU = 1;
if fl_MOD_PC + fl_MOD_USGS_RedHat + fl_MOD_UBUNTU ~= 1
    fprintf('ERROR!! Only one of following can be equal to 1! \n');
    fprintf('fl_MOD_PC. fl_MOD_USGS_RedHat, fl_MOD_UBUNTU \n');
    fprintf('Exiting...\n');
    return
end



% =========================================================================
if ispc 
    slashstr = '\';
elseif isunix
    slashstr = '/';
end


fl_binary = 1;  % 1 for binary
fl_dble = 0;  % 1 for dble prec, 0 for single prec

dry_cell = 1e30;
inactive_cell = -999.99;

% NOTE: for some reason, int w/ bit-length info is trailed by 0!!! (USGS linux)
% windows no bit-length info 
if fl_MOD_USGS_RedHat
    nbitread = 2;
elseif fl_MOD_PC
    nbitread = 0;
elseif fl_MOD_UBUNTU % ubuntu
    nbitread = 1;
end

prec = 4; prec_str = 'real*4';
if fl_dble, prec = 8; prec_str = 'real*8'; end

% -- Get head data and plot it as contour image plots
fid = fopen([MOD_simdir, slashstr, uzf_file], 'r');

ctr = 0; tctr = 0;
% return
nvar = 5;
all_label = cell(nvar,1);
while(1)
    % NOTE: for some reason, int w/ bit-length info is trailed by 0!!!
    a_info = fread(fid, nbitread, 'integer*4');
    kstp = fread(fid, 1, 'integer*4');
    if isempty(kstp)
        break
    end
    kper = fread(fid, 1, 'integer*4');
%     pertim = fread(fid, 1, 'real*4');
%     totim = fread(fid, 1, 'real*4');
    label = fread(fid, 16, 'char');  label = char(label)'; %*** ORIG
    if ctr+1 <= 5, all_label{ctr+1} = label; end
    var_i = mod(ctr, nvar)+1;
%     pause
    % label = fread(fid, 15, 'char');  label = char(label);
    ncol = fread(fid, 1, 'integer*4');
    nrow = fread(fid, 1, 'integer*4');
    ilay = fread(fid, 1, 'integer*4');
    a_info = fread(fid, nbitread, 'integer*4');

    a_data = fread(fid, nbitread, 'integer*4');
    if nbitread == 0, 
        nn = ncol*nrow;
    else
        nn = a_data(1)/prec;
    end
    data = fread(fid, nn, prec_str);
    a_data = fread(fid, nbitread, 'integer*4');

    if nbitread == 0, 
        all_data = reshape(data,ncol,nrow)';
    else
        all_data = reshape(data,ncol,nrow,ilay);
        all_data = permute(all_data, [2 1 3]);
    end
    
    % NOTE: head = -999.99 for no flow cells; ignore these
    all_data(all_data < -999.) = 0;
    all_data(round(all_data) > 0.9e30) = 0;

    % grow array if needed
    if ctr == 0
        if nbitread == 0, 
            all_data_all = zeros(nrow,ncol,nvar,100);
        elseif nbitread == 1
            all_data_all = zeros(nrow,ncol,ilay,nvar,100);
        end
        time_info = zeros(2,100); % kstp, kper, pertim, totim
    elseif mod(ctr,100) == 0
        if nbitread == 0, 
            all_data_all2 = zeros(nrow,ncol,nvar,ctr+100);
            all_data_all2(:,:,:,1:ctr) = all_data_all;
            all_data_all = all_data_all2;
        elseif nbitread == 1
            all_data_all2 = zeros(nrow,ncol,ilay,nvar,ctr+100);
            all_data_all2(:,:,:,:,1:ctr) = all_data_all;
            all_data_all = all_data_all2;            
        end
        time_info2 = zeros(2,ctr+100); % kstp, kper, pertim, totim
        time_info2(:,1:ctr) = time_info;
        time_info = time_info2;
    end
    ctr = ctr + 1;
    if mod(ctr,nvar) == 1, tctr = tctr + 1; end
    if nbitread == 0
        all_data_all(:,:,var_i,tctr) = all_data;
    elseif nbitread == 1
        all_data_all(:,:,:,var_i,tctr) = all_data;
    end
%     time_info(:,ctr) = [kstp, kper, pertim, totim];
    time_info(:,tctr) = [kstp, kper];
end
fclose(fid);

if nbitread == 0
    all_data_all = all_data_all(:,:,1:tctr);
else
    all_data_all = all_data_all(:,:,:,:,1:tctr);
end
time_info = time_info(:,1:tctr);

% read in mask if applicable
if exist('mask_fil', 'var')
    id = fopen(mask_fil, 'r');
    D = textscan(fid, '%s %f', 6); 
    NSEW = D{2}(1:4);
    NROW = D{2}(5);
    NCOL = D{2}(6);
    D = textscan(fid, '%f'); 
    IBOUND = reshape(D{1}, NCOL, NROW)'; % NROW x NCOL
    D = textscan(fid, '%s %s %f %s %f'); 
    dischargePt_rowi = D{3};
    dischargePt_coli = D{5};
    fclose(fid);

    % find boundary cells
    IBOUNDin = IBOUND(2:end-1,2:end-1);
    IBOUNDu = IBOUND(1:end-2,2:end-1); % up
    IBOUNDd = IBOUND(3:end,2:end-1); % down
    IBOUNDl = IBOUND(2:end-1,1:end-2); % left
    IBOUNDr = IBOUND(2:end-1,3:end); % right
%     % - inner boundary is constant head
%     ind_bound = IBOUNDin==1 & (IBOUNDin-IBOUNDu==1 | IBOUNDin-IBOUNDd==1 | ...
%         IBOUNDin-IBOUNDl==1 | IBOUNDin-IBOUNDr==1);
    % - outer boundary is constant head
    ind_bound = IBOUNDin==0 & (IBOUNDin-IBOUNDu==-1 | IBOUNDin-IBOUNDd==-1 | ...
        IBOUNDin-IBOUNDl==-1 | IBOUNDin-IBOUNDr==-1);
% 
%         for ii = 1:NLAY 
%             X = HK(2:end-1,2:end-1,ii);
%             X(ind_bound) = 0;
%             HK(2:end-1,2:end-1,ii) = X;
%     %         Y = HK(:,:,ii);
%     %         Y(IBOUND==0) = 0;
%     %         HK(:,:,ii) = Y;
%         end
end    


% - my interpretation (can't find documentation):
%  var_i = 1 uzf recharge (should be recharge to sat zone?)
%  var_i = 4 surface leakage (should be leakage from surface to unsat zone?)
% for tctr_i = 1:tctr
all_tctr = [14, 15, 39];
caxis0 = zeros(length(all_tctr),2);
caxis0(1,:) = [-1, 4];
caxis0(2,:) = [-10, 165];
caxis0(3,:) = [-10, 165];
for tctr_i = all_tctr(:)'
%     for var_i = 1:nvar
    for var_i = [1] % 
        figure
        for lay_i = 1: 1
%         for lay_i = 1: ilay
%             subplot(2,2,lay_i),
            X = all_data_all(:,:,lay_i,var_i,tctr_i);
            X(IBOUND==0) = 0;
            X2 = X(2:end-1,2:end-1);
            X2(ind_bound) = -999;
            X(2:end-1,2:end-1) = X2;
            
            X = X(2:25,6:44);

            h = imagesc(X); 
            axis equal
            ytickv = get(gca,'YTick');
            yticklab = cell(size(ytickv));

            for ii = 1: length(ytickv)
                yticklab{ii} = num2str(-ytickv(ii));
            end
%             set(gca, 'YTick', ytickv, 'YTickLabel', yticklab);
            m = X(X>0);
            %             m = all_data_all;
            m1 = min(m(:)); m2 = max(m(:));
            
%             if isempty(m)
%                 caxis([0 10])
%             else
%                 caxis([m1*.9 m2])
%                 caxis([-10 m2])
%             end
%             
            caxis(caxis0(tctr_i==all_tctr,:))
            
            cm = colormap;
            cm(1,:) = [0 0 0];
            colormap(cm)

            colorbar;
%             xlabel('x [m]'), ylabel('elev [m]')
            %             time_info(:,tctr_i) % kstp, kper, pertim, totim
            title([all_label{var_i}, ' lay', num2str(lay_i), ', kper ', num2str(time_info(2,tctr_i)), ', kstp ', num2str(time_info(1,tctr_i))]);
            drawnow
        end
    end
%     pause
end

% 7: lo111
% 1111: 
% 14 to 15: transition from dry to wet

nfig = get(gcf, 'Number');
nfig = 3;
for ii = 1:nfig
    print('-dtiff', ['-f', num2str(ii)], [scen_name{scen_i}, num2str(ii), '.tiff']);
    print('-dsvg', ['-f', num2str(ii)], [scen_name{scen_i}, num2str(ii), '.svg']);
end
