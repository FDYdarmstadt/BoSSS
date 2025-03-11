clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% studyName = 'case1';    % 1: m = 2 Oh = 0.1, eta0 = 0.4;   
datapath = strcat(path, 'case1\(0, 0)-WNLT.txt');
import = importdata(datapath);
dat1W0 = import.data;
datapath = strcat(path, 'case1\(2, 0)-WNLT.txt');
import = importdata(datapath);
dat1W2 = import.data;
datapath = strcat(path, 'case1\(4, 0)-WNLT.txt');
import = importdata(datapath);
dat1W4 = import.data;
datapath = strcat(path, 'case1\(6, 0)-WNLT.txt');
import = importdata(datapath);
dat1W6 = import.data;

datapath = strcat(path, 'case1\(0, 0)-Ainit.txt');
import = importdata(datapath);
dat1BA0 = import.data;
datapath = strcat(path, 'case1\(2, 0)-Ainit.txt');
import = importdata(datapath);
dat1BA2 = import.data;
datapath = strcat(path, 'case1\(4, 0)-Ainit.txt');
import = importdata(datapath);
dat1BA4 = import.data;
datapath = strcat(path, 'case1\(6, 0)-Ainit.txt');
import = importdata(datapath);
dat1BA6 = import.data;

datapath = strcat(path, 'case1\StartUpStudy\(0, 0).txt');
import = importdata(datapath);
dat1B0 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\(2, 0).txt');
import = importdata(datapath);
dat1B2 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\(4, 0).txt');
import = importdata(datapath);
dat1B4 = import.data;
datapath = strcat(path, 'case1\StartUpStudy\(6, 0).txt');
import = importdata(datapath);
dat1B6 = import.data;


% studyName = 'case2';    % 2: m = 3 Oh = 0.1, eta0 = 0.4;
datapath = strcat(path, 'case2\(0, 0)-WNLT.txt');
import = importdata(datapath);
dat2W0 = import.data;
datapath = strcat(path, 'case2\(1, 0)-WNLT.txt');
import = importdata(datapath);
dat2W1 = import.data;
datapath = strcat(path, 'case2\(2, 0)-WNLT.txt');
import = importdata(datapath);
dat2W2 = import.data;
datapath = strcat(path, 'case2\(3, 0)-WNLT.txt');
import = importdata(datapath);
dat2W3 = import.data;
datapath = strcat(path, 'case2\(4, 0)-WNLT.txt');
import = importdata(datapath);
dat2W4 = import.data;
datapath = strcat(path, 'case2\(5, 0)-WNLT.txt');
import = importdata(datapath);
dat2W5 = import.data;
datapath = strcat(path, 'case2\(6, 0)-WNLT.txt');
import = importdata(datapath);
dat2W6 = import.data;
datapath = strcat(path, 'case2\(7, 0)-WNLT.txt');
import = importdata(datapath);
dat2W7 = import.data;

datapath = strcat(path, 'case2\(0, 0)-Ainit.txt');
import = importdata(datapath);
dat2B0 = import.data;
datapath = strcat(path, 'case2\(1, 0)-Ainit.txt');
import = importdata(datapath);
dat2B1 = import.data;
datapath = strcat(path, 'case2\(2, 0)-Ainit.txt');
import = importdata(datapath);
dat2B2 = import.data;
datapath = strcat(path, 'case2\(3, 0)-Ainit.txt');
import = importdata(datapath);
dat2B3 = import.data;
datapath = strcat(path, 'case2\(4, 0)-Ainit.txt');
import = importdata(datapath);
dat2B4 = import.data;
datapath = strcat(path, 'case2\(5, 0)-Ainit.txt');
import = importdata(datapath);
dat2B5 = import.data;
datapath = strcat(path, 'case2\(6, 0)-Ainit.txt');
import = importdata(datapath);
dat2B6 = import.data;
datapath = strcat(path, 'case2\(7, 0)-Ainit.txt');
import = importdata(datapath);
dat2B7 = import.data;

% studyName = 'case3';    % 3: m = 4 Oh = 0.1, eta0 = 0.4;    


% study2Name = 'case4';     % 4: m = 2 Oh = 0.1, eta0 = 0.2; 
datapath = strcat(path, 'case4\(0, 0)-WNLT.txt');
import = importdata(datapath);
dat4W0 = import.data;
datapath = strcat(path, 'case4\(2, 0)-WNLT.txt');
import = importdata(datapath);
dat4W2 = import.data;
datapath = strcat(path, 'case4\(4, 0)-WNLT.txt');
import = importdata(datapath);
dat4W4 = import.data;
datapath = strcat(path, 'case4\(6, 0)-WNLT.txt');
import = importdata(datapath);
dat4W6 = import.data;

datapath = strcat(path, 'case4\(0, 0)-Ainit.txt');
import = importdata(datapath);
dat4BA0 = import.data;
datapath = strcat(path, 'case4\(2, 0)-Ainit.txt');
import = importdata(datapath);
dat4BA2 = import.data;
datapath = strcat(path, 'case4\(4, 0)-Ainit.txt');
import = importdata(datapath);
dat4BA4 = import.data;
datapath = strcat(path, 'case4\(6, 0)-Ainit.txt');
import = importdata(datapath);
dat4BA6 = import.data;

datapath = strcat(path, 'case4\StartUpStudy\(0, 0).txt');
import = importdata(datapath);
dat4B0 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\(2, 0).txt');
import = importdata(datapath);
dat4B2 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\(4, 0).txt');
import = importdata(datapath);
dat4B4 = import.data;
datapath = strcat(path, 'case4\StartUpStudy\(6, 0).txt');
import = importdata(datapath);
dat4B6 = import.data;


% study2Name = 'case5';     % 5: m = 4 Oh = 0.56, eta0 = 0.05; 
datapath = strcat(path, 'case5\(0, 0)-WNLT.txt');
import = importdata(datapath);
dat5W0 = import.data;
datapath = strcat(path, 'case5\(2, 0)-WNLT.txt');
import = importdata(datapath);
dat5W2 = import.data;
datapath = strcat(path, 'case5\(4, 0)-WNLT.txt');
import = importdata(datapath);
dat5W4 = import.data;
datapath = strcat(path, 'case5\(6, 0)-WNLT.txt');
import = importdata(datapath);
dat5W6 = import.data;

datapath = strcat(path, 'case5\(0, 0)-Ainit.txt');
import = importdata(datapath);
dat5BA0 = import.data;
datapath = strcat(path, 'case5\(2, 0)-Ainit.txt');
import = importdata(datapath);
dat5BA2 = import.data;
datapath = strcat(path, 'case5\(4, 0)-Ainit.txt');
import = importdata(datapath);
dat5BA4 = import.data;
datapath = strcat(path, 'case5\(6, 0)-Ainit.txt');
import = importdata(datapath);
dat5BA6 = import.data;

datapath = strcat(path, 'case5\StartUpStudy\(0, 0).txt');
import = importdata(datapath);
dat5B0 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\(2, 0).txt');
import = importdata(datapath);
dat5B2 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\(4, 0).txt');
import = importdata(datapath);
dat5B4 = import.data;
datapath = strcat(path, 'case5\StartUpStudy\(6, 0).txt');
import = importdata(datapath);
dat5B6 = import.data;



%% plot

size_legend = 14;
size_label = 16;
size_tick = 16;
size_marker = 6;
mInt = 100; % marker interval for plotting

%%

% case 1
figure
subplot(2,2,1)
plot(dat1W0(:,1), dat1W0(:,2), 'k', dat1BA0(:,1), dat1BA0(:,2), 'b', dat1B0(:,1), dat1B0(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
grid on

subplot(2,2,2)
plot(dat1W2(:,1), dat1W2(:,2), 'k', dat1BA2(:,1), dat1BA2(:,2), 'b', dat1B2(:,1), dat1B2(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5])
grid on

subplot(2,2,3)
plot(dat1W4(:,1), dat1W4(:,2), 'k', dat1BA4(:,1), dat1BA4(:,2), 'b', dat1B4(:,1), dat1B4(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
grid on

subplot(2,2,4)
plot(dat1W6(:,1), dat1W6(:,2), 'k', dat1BA6(:,1), dat1BA6(:,2), 'b', dat1B4(:,1), dat1B4(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5])
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])

%export_fig ModeDecomp_case1_StartUp.tif -r256

%%

% case 2

% equal modes
figure
subplot(2,2,1)
plot(dat2W0(:,1), dat2W0(:,2), 'k', dat2B0(:,1), dat2B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

subplot(2,2,2)
plot(dat2W2(:,1), dat2W2(:,2), 'k', dat2B2(:,1), dat2B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4])
grid on

subplot(2,2,3)
plot(dat2W4(:,1), dat2W4(:,2), 'k', dat2B4(:,1), dat2B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

subplot(2,2,4)
plot(dat2W6(:,1), dat2W6(:,2), 'k', dat2B6(:,1), dat2B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_case2equal_StartUp.tif -r256


% unequal modes
figure
subplot(2,2,1)
plot(dat2W1(:,1), dat2W1(:,2), 'k', dat2B1(:,1), dat2B1(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(1, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

subplot(2,2,2)
plot(dat2W3(:,1), dat2W3(:,2), 'k', dat2B3(:,1), dat2B3(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(3, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4])
grid on

subplot(2,2,3)
plot(dat2W5(:,1), dat2W5(:,2), 'k', dat2B5(:,1), dat2B5(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(5, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

subplot(2,2,4)
plot(dat2W7(:,1), dat2W7(:,2), 'k', dat2B7(:,1), dat2B7(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(7, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_case2unequal_StartUp.tif -r256


%%

% case 4
figure
subplot(2,2,1)
plot(dat4W0(:,1), dat4W0(:,2), 'k', dat4BA0(:,1), dat4BA0(:,2), 'b', dat4B0(:,1), dat4B0(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
grid on

subplot(2,2,2)
plot(dat4W2(:,1), dat4W2(:,2), 'k', dat4BA2(:,1), dat4BA2(:,2), 'b', dat4B2(:,1), dat4B2(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5])
grid on

subplot(2,2,3)
plot(dat4W4(:,1), dat4W4(:,2), 'k', dat4BA4(:,1), dat4BA4(:,2), 'b', dat4B4(:,1), dat4B4(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
grid on

subplot(2,2,4)
plot(dat4W6(:,1), dat4W6(:,2), 'k', dat4BA6(:,1), dat4BA6(:,2), 'b', dat4B6(:,1), dat4B6(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5])
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_case4_StartUp.tif -r256


%%

% case 5
figure
subplot(2,2,1)
plot(dat5W0(:,1), dat5W0(:,2), 'k', dat5BA0(:,1), dat5BA0(:,2), 'b', dat5B0(:,1), dat5B0(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
grid on

subplot(2,2,2)
plot(dat5W2(:,1), dat5W2(:,2), 'k', dat5BA2(:,1), dat5BA2(:,2), 'b', dat5B2(:,1), dat5B2(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southwest')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5])
grid on

subplot(2,2,3)
plot(dat5W4(:,1), dat5W4(:,2), 'k', dat5BA4(:,1), dat5BA4(:,2), 'b', dat5B4(:,1), dat5B4(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5]) 
grid on

subplot(2,2,4)
plot(dat5W6(:,1), dat5W6(:,2), 'k', dat5BA6(:,1), dat5BA6(:,2), 'b', dat5B6(:,1), dat5B6(:,2), 'r', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS-Ainit', 'BoSSS-zeroInit'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,0.5])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5])
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


export_fig ModeDecomp_case5_StartUp.tif -r256
