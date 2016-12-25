% read *.sfr
clear all, close all, fclose all;

% NLAY = 2;
% NROW = 73;
% NCOL = 81;

% sfr_file = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen/input/modflow/sagehen.sfr';
% dis_file = '/home/gcng/workspace/Models/GSFLOW/GSFLOW_1.2.0/data/sagehen/input/modflow/sagehen.dis';
sfr_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/test2lay.sfr';
dis_file = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/inputs/MODFLOW/test2lay.dis';

fid = fopen(sfr_file);

while(1)
    line0 = fgets(fid);
    if ~(strcmp(line0(1), '#') || strncmp(line0, 'REACHINPUT', 10))
        break
    end
end

% line0 = fgets(fid);
D = textscan(line0, '%f');
NSTRM = D{1}(1);

fmt = [repmat('%f ', 1, 14), '\n'];
D = textscan(fid, fmt,abs(NSTRM));
[KRCH,IRCH,JRCH,ISEG,IREACH,RCHLEN,STRTOP,SLOPE,STRTHICK,STRHC1,THTS,THTI,EPS,UHC] = D{:};

fclose(fid);

% -- read in TOP and BOTM from .dis file
fid = fopen(dis_file);

for ii = 1:2
    cmt = fgets(fid);
end

line0 = fgets(fid);
D = textscan(line0, '%d', 6);
NLAY = D{1}(1);
NROW = D{1}(2); 
NCOL = D{1}(3);
NPER = D{1}(4);
ITMUNI = D{1}(5);
LENUNI = D{1}(6);

line0 = fgets(fid);
D = textscan(line0, '%d');
LAYCBD = D{1}; % 1xNLAY (0 if no confining layer)

line0 = fgets(fid);
D = textscan(line0, '%s %d');
DELR = D{2}; % width of column
line0 = fgets(fid);
D = textscan(line0, '%s %d');
DELC = D{2}; % height of row

TOP = nan(NROW,NCOL);
line0 = fgets(fid);
for irow = 1: NROW
    line0 = fgets(fid);
    D = textscan(line0, '%f');
    TOP(irow,:) = D{1}(1:NCOL);
end

BOTM = nan(NROW, NCOL, NLAY);
for ilay = 1: NLAY
    line0 = fgets(fid); 
    for irow = 1: NROW
        line0 = fgets(fid);
        D = textscan(line0, '%f');
        BOTM(irow,:,ilay) = D{1}(1:NCOL);
    end
end
fclose(fid);

% TOP for cells corresponding to reaches
TOP_RCH = nan(abs(NSTRM),1);
for ii = 1:abs(NSTRM), TOP_RCH(ii) = TOP(IRCH(ii),JRCH(ii)); end

% BOTM for cells corresponding to reaches
BOTM_RCH = nan(abs(NSTRM),1);
for ii = 1:abs(NSTRM), BOTM_RCH(ii) = BOTM(IRCH(ii),JRCH(ii),KRCH(ii)); end

% plot
figure
subplot(3,1,1)
plot(TOP_RCH, '-*b'), hold on
plot(BOTM_RCH, '-*k'), hold on
plot(STRTOP, 'r')
legend('top', 'bottom', 'strtop');

subplot(3,1,2)
plot(TOP_RCH - STRTOP, '-*b'), hold on
% plot(BOTM_RCH, '-*k'), hold on
% plot(STRTOP, 'r')
legend('top-strtop');

subplot(3,1,3)
plot(TOP_RCH - BOTM_RCH, '-*b'), hold on
% plot(BOTM_RCH, '-*k'), hold on
% plot(STRTOP, 'r')
legend('top-botm');

% 
% 
% IBOUND = nan(NROW, NCOL, NLAY);
% for ilay = 1: NLAY
%     line0 = fgets(fid); 
%     for irow = 1: NROW
%         line0 = fgets(fid);
%         D = textscan(line0, '%f');
%         IBOUND(irow,:,ilay) = D{1}(1:NCOL);
%     end
% end
% 
% line0 = fgets(fid); 
% D = textscan(line0, '%f %s');
% HNOFLO = D{1};
% 
% initHead = nan(NROW, NCOL, NLAY);
% for ilay = 1: NLAY
%     line0 = fgets(fid); 
%     for irow = 1: NROW
%         line0 = fgets(fid);
%         D = textscan(line0, '%f');
%         initHead(irow,:,ilay) = D{1}(1:NCOL);
%     end
% end
% 
% 
% fclose(fid);
% 
% % make plots
% 
% X = IBOUND;
% figure
% for ilay = 1:NLAY
%     subplot(2,2,double(ilay))
%     m = X(X>0); m = min(m(:));
%     imagesc(X(:,:,ilay)), 
% %     caxis([m*0.9, max(X(:))]), 
%     cm = colormap;
% %     cm(1,:) = [1 1 1];
%     colormap(cm);
%     colorbar
%     title(['IBOUND', ' lay', num2str(ilay)]);
% end
% 
% initHead(initHead<10) = 0; % some outer cells with initHead=1
% X = initHead;
% figure
% for ilay = 1:NLAY
%     subplot(2,2,double(ilay))
%     m = X(X>0); m = min(m(:));
%     imagesc(X(:,:,ilay)), 
%     caxis([m*0.9, max(X(:))]), 
%     cm = colormap;
%     cm(1,:) = [1 1 1];
%     colormap(cm);
%     colorbar
%     title(['init head', ' lay', num2str(ilay)]);
% end
% 
% % check: initHead only specified for IBOUND~=0? (IBOUND>0 is variable, <0
% % is constnat head)
% a = xor(IBOUND, initHead); % xor founds where exactly one is non-zero
% X = a;
% figure
% for ilay = 1:NLAY
%     subplot(2,2,double(ilay))
%     %m = X(X>0); 
%     m = X; m = min(m(:));
%     imagesc(X(:,:,ilay)), 
%     caxis([m*0.9, max(X(:))]), 
%     cm = colormap;
%     cm(1,:) = [1 1 1];
%     colormap(cm);
%     colorbar
%     title(['initHead and ibound differ', ' lay', num2str(ilay)]);
% end
% 
