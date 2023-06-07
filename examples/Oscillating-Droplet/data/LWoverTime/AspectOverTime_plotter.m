%clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% m = 2 Oh = 0.1, eta0 = 0.4;   
datapath = strcat(path, 'm2\AspRatioInTime_m2_0h01_eta04.txt');
import = importdata(datapath);
dat1 = import.data;

% m = 2 Oh = 0.1, eta0 = 0.2;   
datapath = strcat(path, 'm2\AspRatioInTime_m2_0h01_eta02.txt');
import = importdata(datapath);
dat2 = import.data;

% m = 2 Oh = 0.1, eta0 = 0.1;   
datapath = strcat(path, 'm2\AspRatioInTime_m2_0h01_eta01.txt');
import = importdata(datapath);
dat3 = import.data;

% m = 3 Oh = 0.1, eta0 = 0.4;   
datapath = strcat(path, 'm3\AspRatioInTime_m3_0h01_eta04.txt');
import = importdata(datapath);
dat4 = import.data;

% m = 3 Oh = 0.1, eta0 = 0.15;   
datapath = strcat(path, 'm3\AspRatioInTime_m3_0h01_eta015.txt');
import = importdata(datapath);
dat5 = import.data;

% m = 4 Oh = 0.1, eta0 = 0.4;   
datapath = strcat(path, 'm4\AspRatioInTime_m4_0h01_eta04.txt');
import = importdata(datapath);
dat6 = import.data;

% m = 4 Oh = 0.1, eta0 = 0.1;   
datapath = strcat(path, 'm4\AspRatioInTime_m4_0h01_eta01.txt');
import = importdata(datapath);
dat7 = import.data;


%% plot

size_legend = 14;
size_label = 16;
size_tick = 16;
size_marker = 6;

%%

figure
title('m = 2')
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1(:,1), dat1(:,2), dat2(:,1), dat2(:,2), dat3(:,1), dat3(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'eta=0.4', 'eta=0.2', 'eta=0.1'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%%

figure
title('eta_0 = 0.4')
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1(:,1), dat1(:,2), dat4(:,1), dat4(:,2), dat6(:,1), dat6(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'m=2', 'm=3', 'm=4'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on
mInt = 100; % marker interval for plotting