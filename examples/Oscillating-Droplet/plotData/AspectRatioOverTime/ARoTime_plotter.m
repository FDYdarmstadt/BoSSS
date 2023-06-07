clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% 1: m = 2 Oh = 0.1, eta0 = 0.4;   
path = 'm2/m2_Oh01_eta04_';
datapath = strcat(path, 'Basaran.txt');
import = importdata(datapath);
dat1Ba = import.data;
datapath = strcat(path, 'Becker.txt');
import = importdata(datapath);
dat1Be = import.data;
datapath = strcat(path, 'Meradji.txt');
import = importdata(datapath);
dat1Me = import.data;
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat1W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat1Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat1Bthird = import.data;

% 2: m = 2 Oh = 0.1, eta0 = 0.2;
path = 'm2/m2_Oh01_eta02_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat2W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat2Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat2Bthird = import.data;

% 3: m = 2 Oh = 0.1, eta0 = 0.1;
path = 'm2/m2_Oh01_eta01_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat3W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat3Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat3Bthird = import.data;

% 4: m = 3 Oh = 0.1, eta0 = 0.4; 
path = 'm3/m3_Oh01_eta04_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat4W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat4Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat4Bthird = import.data;

% % 5: m = 3 Oh = 0.1, eta0 = 0.3; 
% path = 'm3/m3_Oh01_eta03_';
% datapath = strcat(path, 'WNLT.txt');
% import = importdata(datapath);
% dat5W = import.data;
% datapath = strcat(path, 'BoSSS-J686.txt');
% import = importdata(datapath);
% dat5B = import.data;

% 6: m = 3 Oh = 0.1, eta0 = 0.15; 
path = 'm3/m3_Oh01_eta015_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat6W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat6Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat6Bthird = import.data;

% 7: m = 4 Oh = 0.1, eta0 = 0.4; 
path = 'm4/m4_Oh01_eta04_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat7W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat7Bzero = import.data;
% datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
% import = importdata(datapath);
% dat7Bthird = import.data;

% 8: m = 4 Oh = 0.1, eta0 = 0.1; 
path = 'm4/m4_Oh01_eta01_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat8W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat8Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat8Bthird = import.data;

% 9: m = 4 Oh = 0.56, eta0 = 0.05; 
path = 'm4/m4_Oh056_eta005_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat9W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat9Bzero = import.data;
datapath = strcat(path, 'BoSSS-J686-3OrdInit.txt');
import = importdata(datapath);
dat9Bthird = import.data;
datapath = strcat(path, 'BoSSS-J432.txt');
import = importdata(datapath);
dat9Bzero2 = import.data;
datapath = strcat(path, 'BoSSS-J432-3OrdInit.txt');
import = importdata(datapath);
dat9Bthird2 = import.data;


%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;
mInt = 100; % marker interval for plotting


%%

% case 1

dat1Bthird_temp = zeros(length(dat1Bthird),2);
cnt = 1;
for i=1:length(dat1Bthird)
    if dat1Bthird(i,2) > 0.2 
        dat1Bthird_temp(cnt,1) = dat1Bthird(i,1);
        dat1Bthird_temp(cnt,2) = dat1Bthird(i,2);
        cnt = cnt + 1;
    end 
end
dat1Bthird = dat1Bthird_temp(1:cnt-1,:);

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
% plot(dat1Me(:,1), dat1Me(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'Analytic Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
plot(dat1Ba(:,1), dat1Ba(:,2), 'm:', dat1Be(:,1), dat1Be(:,2), 'g:', dat1Me(:,1), dat1Me(:,2), ':', dat1W(:,1), dat1W(:,2), 'k-', dat1Bzero(:,1), dat1Bzero(:,2), 'b-.', dat1Bthird(:,1), dat1Bthird(:,2), 'b--', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Basaran (1992)', 'Becker et al. (1994)', 'Meradji et al. (2001)', 'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m2_Oh01_eta04.tif -r256
%saveas(gcf, 'ARoverTime_m2_Oh01_eta04', 'epsc')


%%

% case 2
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat2W(:,1), dat2W(:,2), 'k-', dat2Bzero(:,1), dat2Bzero(:,2), 'b-.', dat2Bthird(:,1), dat2Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = 0.2$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m2_Oh01_eta02.tif -r256
%saveas(gcf, 'ARoverTime_m2_Oh01_eta02', 'epsc')


%%

% case 3
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat3W(:,1), dat3W(:,2), 'k-', dat3Bzero(:,1), dat3Bzero(:,2), 'b-.', dat3Bthird(:,1), dat3Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = 0.1$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m2_Oh01_eta01.tif -r256
%saveas(gcf, 'ARoverTime_m2_Oh01_eta01', 'epsc')


%%

% case 4
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat4W(:,1), dat4W(:,2), 'k', dat4Bzero(:,1), dat4Bzero(:,2), 'b', dat4Bthird(:,1), dat4Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=3$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m3_Oh01_eta04.tif -r256
%saveas(gcf, 'ARoverTime_m3_Oh01_eta04', 'epsc')


%%

% % case 5
% figure
% set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
% plot(dat5W(:,1), dat5W(:,2), 'k', dat5Bzero(:,1), dat5Bzero(:,2), 'b', dat5Bthird(:,1), dat5Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% title('$m=3$ ${\rm Oh}=0.1$ $\eta_0 = 0.3$', 'Interpreter', 'latex', 'Fontsize', size_legend)
% legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4], 'FontSize', size_tick) 
% %set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
% set(gcf, 'color', 'w')
% grid on
% 
% export_fig ARoveTime_m3_Oh01_eta03.png  %.tif -r256
% saveas(gcf, 'ARoverTime_m3_Oh01_eta03', 'epsc')


%%

% case 6
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat6W(:,1), dat6W(:,2), 'k-', dat6Bzero(:,1), dat6Bzero(:,2), 'b-.', dat6Bthird(:,1), dat6Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=3$ ${\rm Oh}=0.1$ $\eta_0 = 0.15$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m3_Oh01_eta015.tif -r256
%saveas(gcf, 'ARoverTime_m3_Oh01_eta015', 'epsc')


%%

% case 7
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat7W(:,1), dat7W(:,2), 'k-', dat7Bzero(:,1), dat7Bzero(:,2), 'b-.', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m4_Oh01_eta04.tif -r256
%saveas(gcf, 'ARoverTime_m4_Oh01_eta04', 'epsc')


%%

% case 8
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat8W(:,1), dat8W(:,2), 'k-', dat8Bzero(:,1), dat8Bzero(:,2), 'b-.', dat8Bthird(:,1), dat8Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$ $\eta_0 = 0.1$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m4_Oh01_eta01.tif -r256
%saveas(gcf, 'ARoverTime_m4_Oh01_eta01', 'epsc')


%%

% case 9
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat9W(:,1), dat9W(:,2), 'k-', dat9Bzero(:,1), dat9Bzero(:,2), 'b-.', dat9Bthird(:,1), dat9Bthird(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.56$ $\eta_0 = 0.05$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022)', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig ARoveTime_m4_Oh056_eta005.tif -r256
%saveas(gcf, 'ARoverTime_m4_Oh056_eta005', 'epsc')
