clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% studyName = 'case1';    % 1: m = 2 Oh = 0.1, eta0 = 0.4;   
datapath = strcat(path, 'case1\Basaran (1992).txt');
import = importdata(datapath);
dat1Ba = import.data;
datapath = strcat(path, 'case1\Becker et. al. (1994).txt');
import = importdata(datapath);
dat1Be = import.data;
datapath = strcat(path, 'case1\Meradji et. al. (2001).txt');
import = importdata(datapath);
dat1M = import.data;
datapath = strcat(path, 'case1\WNLT.txt');
import = importdata(datapath);
dat1W = import.data;
datapath = strcat(path, 'case1\BoSSS--Ainit.txt');
import = importdata(datapath);
dat1B = import.data;

% studyName = 'case2';    % 2: m = 3 Oh = 0.1, eta0 = 0.4;
datapath = strcat(path, 'case2\WNLT.txt');
import = importdata(datapath);
dat2W = import.data;
datapath = strcat(path, 'case2\BoSSS--Ainit.txt');
import = importdata(datapath);
dat2B = import.data;

% studyName = 'case3';    % 3: m = 4 Oh = 0.1, eta0 = 0.4;    

% study2Name = 'case4';     % 4: m = 2 Oh = 0.1, eta0 = 0.2; 
datapath = strcat(path, 'case4\WNLT.txt');
import = importdata(datapath);
dat4W = import.data;
datapath = strcat(path, 'case4\BoSSS--Ainit.txt');
import = importdata(datapath);
dat4B = import.data;

% study2Name = 'case5';     % 5: m = 4 Oh = 0.56, eta0 = 0.05; 
datapath = strcat(path, 'case5\WNLT.txt');
import = importdata(datapath);
dat5W = import.data;
datapath = strcat(path, 'case5\BoSSS--Ainit.txt');
import = importdata(datapath);
dat5B = import.data;


%% plot

size_legend = 14;
size_label = 16;
size_tick = 16;
size_marker = 6;
mInt = 100; % marker interval for plotting

%%

% case 1
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
% plot(dat1W(:,1), dat1W(:,2), dat1B(:,1), dat1B(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'Analytic Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
plot(dat1Ba(:,1), dat1Ba(:,2), 'm', dat1Be(:,1), dat1Be(:,2), 'g', dat1M(:,1), dat1M(:,2), dat1W(:,1), dat1W(:,2), 'k', dat1B(:,1), dat1B(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'Basaran (1992)', 'Becker et. al. (1994)', 'Meradji et. al. (2001)', 'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig ARoveTime_case1.tif -r256

%%

% case 2
figure
%set(gcf,'DefaultAxesColorOrder',[0.0 0.0 0.0])
plot(dat2W(:,1), dat2W(:,2), 'k', dat2B(:,1), dat2B(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig ARoveTime_case2.tif -r256

%%

% case 4
figure
%set(gcf,'DefaultAxesColorOrder',[0.0 0.0 0.0])
plot(dat4W(:,1), dat4W(:,2), 'k', dat4B(:,1), dat4B(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig ARoveTime_case4.tif -r256

%%

% case 5
figure
%set(gcf,'DefaultAxesColorOrder',[0.0 0.0 0.0])
plot(dat5W(:,1), dat5W(:,2), 'k', dat5B(:,1), dat5B(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig ARoveTime_case5.tif -r256
