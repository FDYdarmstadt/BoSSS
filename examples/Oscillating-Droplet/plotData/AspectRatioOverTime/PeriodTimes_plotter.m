clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% 1: m = 2 Oh = 0.1, eta0 = 0.4;   
path = 'm2/m2_Oh01_eta04_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat1W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat1Bzero = import.data;

datapath = 'm2/periodTimes/periodTimes_m2_Oh01_eta04.txt';
import = importdata(datapath);
dat1Wpt = import;

datapath = 'm2/dampingRates/dampingRates_m2_Oh01_eta04.txt';
import = importdata(datapath);
dat1Wdr = import;


% 2: m = 2 Oh = 0.1, eta0 = 0.2;   
path = 'm2/m2_Oh01_eta02_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat2W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat2Bzero = import.data;

datapath = 'm2/periodTimes/periodTimes_m2_Oh01_eta02.txt';
import = importdata(datapath);
dat2Wpt = import;

datapath = 'm2/dampingRates/dampingRates_m2_Oh01_eta02.txt';
import = importdata(datapath);
dat2Wdr = import;


% 3: m = 2 Oh = 0.1, eta0 = 0.1;   
path = 'm2/m2_Oh01_eta01_';
datapath = strcat(path, 'WNLT.txt');
import = importdata(datapath);
dat3W = import.data;
datapath = strcat(path, 'BoSSS-J686.txt');
import = importdata(datapath);
dat3Bzero = import.data;



%% evaluate period times

times1 = dat1Bzero(:,1);
values1 = dat1Bzero(:,2);
[periodTimes1, dampingRates1] = GetPeriodTimesAndDampingRates(times1, values1, 2);

freq0 = sqrt(8-(0.5)^2);
freqChange1 = ((2*pi / periodTimes1(1,2)) - freq0) / freq0;
LW1 = dat1Bzero(1,2);


times2 = dat2Bzero(:,1);
values2 = dat2Bzero(:,2);
[periodTimes2, dampingRates2] = GetPeriodTimesAndDampingRates(times2, values2, 2);

freqChange2 = ((2*pi / periodTimes2(1,2)) - freq0) / freq0;
LW2 = dat2Bzero(1,2);


times3 = dat3Bzero(:,1);
values3 = dat3Bzero(:,2);
[periodTimes3, dampingRates3] = GetPeriodTimesAndDampingRates(times3, values3, 2);

freqChange3 = ((2*pi / periodTimes3(1,2)) - freq0) / freq0;
LW3 = dat3Bzero(1,2);


%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;
mInt = 100; % marker interval for plotting


%% aspect ratio over time

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1W(:,1), dat1W(:,2), 'k--', dat1Bzero(:,1), dat1Bzero(:,2), 'b--', dat2W(:,1), dat2W(:,2), 'k-.', dat2Bzero(:,1), dat2Bzero(:,2), 'b-.', dat3W(:,1), dat3W(:,2), 'k:', dat3Bzero(:,1), dat3Bzero(:,2), 'b:', 'Linewidth', 1.2, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'Zrnic et al. (2022) - $\eta_0 = 0.4$', 'BoSSS - zero Init - $\eta_0 = 0.4$', 'Zrnic et al. (2022) - $\eta_0 = 0.2$', 'BoSSS - zero Init - $\eta_0 = 0.2$', 'Zrnic et al. (2022) - $\eta_0 = 0.1$', 'BoSSS - zero Init - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$L/W$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

%%

figure
scatter(periodTimes1(:,1), periodTimes1(:,2), 'b', 'filled')
hold on
scatter(dat1Wpt(:,1), dat1Wpt(:,2), 'k')
%figure
hold on
scatter(periodTimes2(:,1), periodTimes2(:,2), 'bs', 'filled')
hold on
scatter(dat2Wpt(:,1), dat2Wpt(:,2), 'ks')
%figure
hold on
scatter(periodTimes3(:,1), periodTimes3(:,2), 'bd', 'filled')

%%

figure
scatter(dampingRates1(:, 1), dampingRates1(:, 2), 'b', 'filled')
hold on
scatter(dat1Wdr(:,1), dat1Wdr(:,2), 'k')
%figure
hold on
scatter(dampingRates2(:, 1), dampingRates2(:, 2), 'bs', 'filled')
hold on
scatter(dat2Wdr(:,1), dat2Wdr(:,2), 'ks')
%figure
hold on
scatter(dampingRates3(:, 1), dampingRates3(:, 2), 'bd', 'filled')


