clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder  
fileExtension = ' (4, 0) J686.txt';
          
% 1: m = 2 Oh = 0.1, eta0 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};     
datapath = strcat(path, 'm2/m2_Oh01_eta01', fileExtension);
import = importdata(datapath);
dat21 = import.data;
datapath = strcat(path, 'm2/m2_Oh01_eta02', fileExtension);
import = importdata(datapath);
dat22 = import.data;
datapath = strcat(path, 'm2/m2_Oh01_eta03', fileExtension);
import = importdata(datapath);
dat23 = import.data;
datapath = strcat(path, 'm2/m2_Oh01_eta04', fileExtension);
import = importdata(datapath);
dat24 = import.data;
datapath = strcat(path, 'm2/m2_Oh01_eta05', fileExtension);
import = importdata(datapath);
dat25 = import.data;
datapath = strcat(path, 'm2/m2_Oh01_eta06', fileExtension);
import = importdata(datapath);
dat26 = import.data;
datapath = strcat(path, 'm2/m2_Oh01_eta07', fileExtension);
import = importdata(datapath);
dat27 = import.data;


% % 2: m = 3 Oh = 0.1, eta0 = {0.15, 0.3, 0.4, 0.5, 0.6, 0.7};     
% datapath = strcat(path, 'm3/m3_Oh01_eta015', fileExtension);
% import = importdata(datapath);
% dat31 = import.data;
% datapath = strcat(path, 'm3/m3_Oh01_eta03', fileExtension);
% import = importdata(datapath);
% dat32 = import.data;
% datapath = strcat(path, 'm3/m3_Oh01_eta04', fileExtension);
% import = importdata(datapath);
% dat33 = import.data;
% % datapath = strcat(path, 'm3/m3_Oh01_eta05', fileExtension);
% % import = importdata(datapath);
% % dat34 = import.data;
% % datapath = strcat(path, 'm3/m3_Oh01_eta06', fileExtension);
% % import = importdata(datapath);
% % dat35 = import.data;
% % datapath = strcat(path, 'm3/m3_Oh01_eta07', fileExtension);
% % import = importdata(datapath);
% % dat36 = import.data;


% 3: m = 4 Oh = 0.1, eta0 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};     
datapath = strcat(path, 'm4/m4_Oh01_eta01', fileExtension);
import = importdata(datapath);
dat41 = import.data;
datapath = strcat(path, 'm4/m4_Oh01_eta02', fileExtension);
import = importdata(datapath);
dat42 = import.data;
datapath = strcat(path, 'm4/m4_Oh01_eta03', fileExtension);
import = importdata(datapath);
dat43 = import.data;
datapath = strcat(path, 'm4/m4_Oh01_eta04', fileExtension);
import = importdata(datapath);
dat44 = import.data;
datapath = strcat(path, 'm4/m4_Oh01_eta05', fileExtension);
import = importdata(datapath);
dat45 = import.data;
datapath = strcat(path, 'm4/m4_Oh01_eta06', fileExtension);
import = importdata(datapath);
dat46 = import.data;
datapath = strcat(path, 'm4/m4_Oh01_eta07', fileExtension);
import = importdata(datapath);
dat47 = import.data;


%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;
mInt = 100; % marker interval for plotting

%%

% mode (0,0)
figure
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat21(:,1), dat21(:,2), dat22(:,1), dat22(:,2), dat23(:,1), dat23(:,2), dat24(:,1), dat24(:,2), dat25(:,1), dat25(:,2), dat26(:,1), dat26(:,2), dat27(:,1), dat27(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'$\eta = 0.1$', '$\eta = 0.2$', '$\eta = 0.3$', '$\eta = 0.4$', '$\eta = 0.5$', '$\eta = 0.6$', '$\eta = 0.7$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(4,0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig modeDecomp4_m2_Oh01_FirstPeriodStudy.tif -r256
saveas(gcf, 'modeDecomp4_m2_Oh01_FirstPeriodStudy', 'epsc')


%%

% % mode (0,0)
% figure
% %set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
% plot(dat31(:,1), dat31(:,2), dat32(:,1), dat32(:,2), dat33(:,1), dat33(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
% title('$m=3$ ${\rm Oh}=0.1$', 'Interpreter', 'latex', 'Fontsize', size_legend)
% legend({'$\eta = 0.15$', '$\eta = 0.3$', '$\eta = 0.4$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(0,0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
% %set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
% set(gcf, 'color', 'w')
% grid on
% 
% export_fig modeDecomp0_m3_Oh01_FirstPeriodStudy.png  %.tif -r256


%%

% mode (0,0)
figure
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat41(:,1), dat41(:,2), dat42(:,1), dat42(:,2), dat43(:,1), dat43(:,2), dat44(:,1), dat44(:,2), dat45(:,1), dat45(:,2), dat46(:,1), dat46(:,2), dat47(:,1), dat47(:,2),  'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'$\eta = 0.1$', '$\eta = 0.2$', '$\eta = 0.3$', '$\eta = 0.4$', '$\eta = 0.5$', '$\eta = 0.6$', '$\eta = 0.7$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,1.5])
ylabel('$(4,0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.5 1 1.5 2 3], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig modeDecomp4_m4_Oh01_FirstPeriodStudy.tif -r256
saveas(gcf, 'modeDecomp4_m4_Oh01_FirstPeriodStudy', 'epsc')


