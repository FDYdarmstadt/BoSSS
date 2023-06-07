clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% studyName = 'case1';    % 1: m = 2 Oh = 0.1, eta0 = 0.4;   
datapath = strcat(path, 'case1\StartUpStudy\WNLT.txt');
import = importdata(datapath);
dat1W = import.data;
datapath = strcat(path, 'case1\StartUpStudy\BoSSS-wallBC.txt');
import = importdata(datapath);
dat10 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\BoSSS-wallBC-Ainit_.txt');
import = importdata(datapath);
dat11 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\BoSSS-wallBC-Ainit_initCase0.txt');
import = importdata(datapath);
dat12 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\BoSSS-wallBC-Ainit_initCase1.txt');
import = importdata(datapath);
dat13 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\BoSSS-wallBC-Ainit_initCase2.txt');
import = importdata(datapath);
dat14 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\BoSSS-wallBC-Ainit_switchedComponents.txt');
import = importdata(datapath);
dat15 = import.data;

% studyName = 'case2';    % 2: m = 3 Oh = 0.1, eta0 = 0.4;  

% studyName = 'case3';    % 3: m = 4 Oh = 0.1, eta0 = 0.4;    

% study2Name = 'case4';     % 4: m = 2 Oh = 0.1, eta0 = 0.2; 
datapath = strcat(path, 'case4\StartUpStudy\WNLT.txt');
import = importdata(datapath);
dat4W = import.data;
datapath = strcat(path, 'case4\StartUpStudy\BoSSS-wallBC.txt');
import = importdata(datapath);
dat40 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\BoSSS-wallBC-Ainit_.txt');
import = importdata(datapath);
dat41 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\BoSSS-wallBC-Ainit_initCase0.txt');
import = importdata(datapath);
dat42 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\BoSSS-wallBC-Ainit_initCase1.txt');
import = importdata(datapath);
dat43 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\BoSSS-wallBC-Ainit_initCase2.txt');
import = importdata(datapath);
dat44 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\BoSSS-wallBC-Ainit_switchedComponents.txt');
import = importdata(datapath);
dat45 = import.data;

% study2Name = 'case5';     % 5: m = 4 Oh = 0.56, eta0 = 0.05; 
datapath = strcat(path, 'case5\StartUpStudy\WNLT.txt');
import = importdata(datapath);
dat5W = import.data;
datapath = strcat(path, 'case5\StartUpStudy\BoSSS-wallBC.txt');
import = importdata(datapath);
dat50 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\BoSSS-wallBC-Ainit_.txt');
import = importdata(datapath);
dat51 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\BoSSS-wallBC-Ainit_initCase0.txt');
import = importdata(datapath);
dat52 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\BoSSS-wallBC-Ainit_initCase1.txt');
import = importdata(datapath);
dat53 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\BoSSS-wallBC-Ainit_initCase2.txt');
import = importdata(datapath);
dat54 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\BoSSS-wallBC-Ainit_switchedComponents.txt');
import = importdata(datapath);
dat55 = import.data;


%% plot

size_legend = 14;
size_label = 16;
size_tick = 16;
size_marker = 6;
mInt = 100; % marker interval for plotting

%%

% case 1
figure
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1W(:,1), dat1W(:,2), 'k', dat10(:,1), dat10(:,2), 'r', dat11(:,1), dat11(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'zero velocity field', 'Analytic Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig StartUpStudy_case1.tif -r256

% complete
figure
set(gcf,'DefaultAxesColorOrder', [0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330])
plot(dat1W(:,1), dat1W(:,2), 'k', dat10(:,1), dat10(:,2), 'r', dat11(:,1), dat11(:,2), 'b', dat12(:,1), dat12(:,2), 'g', dat13(:,1), dat13(:,2), 'm', dat14(:,1), dat14(:,2), dat15(:,1), dat15(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'zero velocity field', 'Analytic Init', 'Analytic Init (neg pol vel)', 'Analytic Init (neg rad vel)', 'Analytic Init (neg vel vec)', 'Analytic Init (switched comp)'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 700, 500])

export_fig StartUpStudy_case1_complete.tif -r256

%%

% case 4
figure
%set(gcf,'DefaultAxesColorOrder',[0.0 0.0 0.0])
plot(dat4W(:,1), dat4W(:,2), 'k', dat40(:,1), dat40(:,2), 'r', dat41(:,1), dat41(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'zero velocity field', 'Analytic Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig StartUpStudy_case4.tif -r256

% complete
figure
set(gcf,'DefaultAxesColorOrder', [0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330])
plot(dat4W(:,1), dat4W(:,2), 'k', dat40(:,1), dat40(:,2), 'r', dat41(:,1), dat41(:,2), 'b', dat42(:,1), dat42(:,2), 'g', dat43(:,1), dat43(:,2), 'm', dat44(:,1), dat44(:,2), dat45(:,1), dat45(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'zero velocity field', 'Analytic Init', 'Analytic Init (neg pol vel)', 'Analytic Init (neg rad vel)', 'Analytic Init (neg vel vec)', 'Analytic Init (switched comp)'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 700, 500])

%export_fig StartUpStudy_case4_complete.tif -r256

%%

% case 5
figure
%set(gcf,'DefaultAxesColorOrder',[0.0 0.0 0.0])
plot(dat5W(:,1), dat5W(:,2), 'k', dat50(:,1), dat50(:,2), 'r', dat51(:,1), dat51(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'zero velocity field', 'Analytic Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%export_fig StartUpStudy_case5.tif -r256

% complete
figure
set(gcf,'DefaultAxesColorOrder', [0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330])
plot(dat5W(:,1), dat5W(:,2), 'k', dat50(:,1), dat50(:,2), 'r', dat51(:,1), dat51(:,2), 'b', dat52(:,1), dat52(:,2), 'g', dat53(:,1), dat53(:,2), 'm', dat54(:,1), dat54(:,2), dat55(:,1), dat55(:,2), 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'zero velocity field', 'Analytic Init', 'Analytic Init (neg pol vel)', 'Analytic Init (neg rad vel)', 'Analytic Init (neg vel vec)', 'Analytic Init (switched comp)'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 700, 500])

export_fig StartUpStudy_case5_complete.tif -r256
