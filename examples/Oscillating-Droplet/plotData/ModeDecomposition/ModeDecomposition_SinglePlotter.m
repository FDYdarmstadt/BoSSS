clear all

%% import data

path = '';  %'plotData\AspectRatioOverTime'; % path to data folder     
          
% 1: m = 2 Oh = 0.1, eta0 = 0.4;   
path = 'm2/m2_Oh01_eta04';
datapath = strcat(path, ' (0, 0)-WNLT.txt');
import = importdata(datapath);
dat1W0 = import.data;
datapath = strcat(path, ' (2, 0)-WNLT.txt');
import = importdata(datapath);
dat1W2 = import.data;
datapath = strcat(path, ' (4, 0)-WNLT.txt');
import = importdata(datapath);
dat1W4 = import.data;
datapath = strcat(path, ' (6, 0)-WNLT.txt');
import = importdata(datapath);
dat1W6 = import.data;

datapath = strcat(path, ' (0, 0) J686.txt');
import = importdata(datapath);
dat1B0 = import.data;
datapath = strcat(path, ' (2, 0) J686.txt');
import = importdata(datapath);
dat1B2 = import.data;
datapath = strcat(path, ' (4, 0) J686.txt');
import = importdata(datapath);
dat1B4 = import.data;
datapath = strcat(path, ' (6, 0) J686.txt');
import = importdata(datapath);
dat1B6 = import.data;
  

% 2: m = 2 Oh = 0.1, eta0 = 0.2;   
path = 'm2/m2_Oh01_eta02';
datapath = strcat(path, ' (0, 0)-WNLT.txt');
import = importdata(datapath);
dat2W0 = import.data;
datapath = strcat(path, ' (2, 0)-WNLT.txt');
import = importdata(datapath);
dat2W2 = import.data;
datapath = strcat(path, ' (4, 0)-WNLT.txt');
import = importdata(datapath);
dat2W4 = import.data;
datapath = strcat(path, ' (6, 0)-WNLT.txt');
import = importdata(datapath);
dat2W6 = import.data;

datapath = strcat(path, ' (0, 0) J686.txt');
import = importdata(datapath);
dat2B0 = import.data;
datapath = strcat(path, ' (2, 0) J686.txt');
import = importdata(datapath);
dat2B2 = import.data;
datapath = strcat(path, ' (4, 0) J686.txt');
import = importdata(datapath);
dat2B4 = import.data;
datapath = strcat(path, ' (6, 0) J686.txt');
import = importdata(datapath);
dat2B6 = import.data;

          
% 3: m = 2 Oh = 0.1, eta0 = 0.1;   
path = 'm2/m2_Oh01_eta01';
datapath = strcat(path, ' (0, 0)-WNLT.txt');
import = importdata(datapath);
dat3W0 = import.data;
datapath = strcat(path, ' (2, 0)-WNLT.txt');
import = importdata(datapath);
dat3W2 = import.data;
datapath = strcat(path, ' (4, 0)-WNLT.txt');
import = importdata(datapath);
dat3W4 = import.data;
datapath = strcat(path, ' (6, 0)-WNLT.txt');
import = importdata(datapath);
dat3W6 = import.data;

datapath = strcat(path, ' (0, 0) J686.txt');
import = importdata(datapath);
dat3B0 = import.data;
datapath = strcat(path, ' (2, 0) J686.txt');
import = importdata(datapath);
dat3B2 = import.data;
datapath = strcat(path, ' (4, 0) J686.txt');
import = importdata(datapath);
dat3B4 = import.data;
datapath = strcat(path, ' (6, 0) J686.txt');
import = importdata(datapath);
dat3B6 = import.data;


% % 4: m = 3 Oh = 0.1, eta0 = 0.4;
% datapath = strcat(path, 'case2\(0, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W0 = import.data;
% datapath = strcat(path, 'case2\(1, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W1 = import.data;
% datapath = strcat(path, 'case2\(2, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W2 = import.data;
% datapath = strcat(path, 'case2\(3, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W3 = import.data;
% datapath = strcat(path, 'case2\(4, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W4 = import.data;
% datapath = strcat(path, 'case2\(5, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W5 = import.data;
% datapath = strcat(path, 'case2\(6, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W6 = import.data;
% datapath = strcat(path, 'case2\(7, 0)-WNLT.txt');
% import = importdata(datapath);
% dat2W7 = import.data;
% 
% datapath = strcat(path, 'case2\(0, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B0 = import.data;
% datapath = strcat(path, 'case2\(1, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B1 = import.data;
% datapath = strcat(path, 'case2\(2, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B2 = import.data;
% datapath = strcat(path, 'case2\(3, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B3 = import.data;
% datapath = strcat(path, 'case2\(4, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B4 = import.data;
% datapath = strcat(path, 'case2\(5, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B5 = import.data;
% datapath = strcat(path, 'case2\(6, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B6 = import.data;
% datapath = strcat(path, 'case2\(7, 0)-Ainit.txt');
% import = importdata(datapath);
% dat2B7 = import.data;


% 7: m = 4 Oh = 0.1, eta0 = 0.4;
path = 'm4/m4_Oh01_eta04';
datapath = strcat(path, ' (0, 0)-WNLT.txt');
import = importdata(datapath);
dat7W0 = import.data;
datapath = strcat(path, ' (2, 0)-WNLT.txt');
import = importdata(datapath);
dat7W2 = import.data;
datapath = strcat(path, ' (4, 0)-WNLT.txt');
import = importdata(datapath);
dat7W4 = import.data;
datapath = strcat(path, ' (6, 0)-WNLT.txt');
import = importdata(datapath);
dat7W6 = import.data;

datapath = strcat(path, ' (0, 0) J686.txt');
import = importdata(datapath);
dat7B0 = import.data;
datapath = strcat(path, ' (2, 0) J686.txt');
import = importdata(datapath);
dat7B2 = import.data;
datapath = strcat(path, ' (4, 0) J686.txt');
import = importdata(datapath);
dat7B4 = import.data;
datapath = strcat(path, ' (6, 0) J686.txt');
import = importdata(datapath);
dat7B6 = import.data;


% 8: m = 4 Oh = 0.1, eta0 = 0.1;
path = 'm4/m4_Oh01_eta01';
datapath = strcat(path, ' (0, 0)-WNLT.txt');
import = importdata(datapath);
dat8W0 = import.data;
datapath = strcat(path, ' (2, 0)-WNLT.txt');
import = importdata(datapath);
dat8W2 = import.data;
datapath = strcat(path, ' (4, 0)-WNLT.txt');
import = importdata(datapath);
dat8W4 = import.data;
datapath = strcat(path, ' (6, 0)-WNLT.txt');
import = importdata(datapath);
dat8W6 = import.data;

datapath = strcat(path, ' (0, 0) J686.txt');
import = importdata(datapath);
dat8B0 = import.data;
datapath = strcat(path, ' (2, 0) J686.txt');
import = importdata(datapath);
dat8B2 = import.data;
datapath = strcat(path, ' (4, 0) J686.txt');
import = importdata(datapath);
dat8B4 = import.data;
datapath = strcat(path, ' (6, 0) J686.txt');
import = importdata(datapath);
dat8B6 = import.data;


% 9: m = 4 Oh = 0.56, eta0 = 0.05; 
path = 'm4/m4_Oh056_eta005';
datapath = strcat(path, ' (0, 0)-WNLT.txt');
import = importdata(datapath);
dat9W0 = import.data;
datapath = strcat(path, ' (2, 0)-WNLT.txt');
import = importdata(datapath);
dat9W2 = import.data;
datapath = strcat(path, ' (4, 0)-WNLT.txt');
import = importdata(datapath);
dat9W4 = import.data;
datapath = strcat(path, ' (6, 0)-WNLT.txt');
import = importdata(datapath);
dat9W6 = import.data;

datapath = strcat(path, ' (0, 0) J686.txt');
import = importdata(datapath);
dat9B0 = import.data;
datapath = strcat(path, ' (2, 0) J686.txt');
import = importdata(datapath);
dat9B2 = import.data;
datapath = strcat(path, ' (4, 0) J686.txt');
import = importdata(datapath);
dat9B4 = import.data;
datapath = strcat(path, ' (6, 0) J686.txt');
import = importdata(datapath);
dat9B6 = import.data;

datapath = strcat(path, ' (0, 0) J686-3OrdInit.txt');
import = importdata(datapath);
dat9B02 = import.data;
datapath = strcat(path, ' (2, 0) J686-3OrdInit.txt');
import = importdata(datapath);
dat9B22 = import.data;
datapath = strcat(path, ' (4, 0) J686-3OrdInit.txt');
import = importdata(datapath);
dat9B42 = import.data;
datapath = strcat(path, ' (6, 0) J686-3OrdInit.txt');
import = importdata(datapath);
dat9B62 = import.data;



%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;
mInt = 100; % marker interval for plotting

%%

% case 1
figure
subplot(2,2,1)
plot(dat1W0(:,1), dat1W0(:,2), 'k', dat1B0(:,1), dat1B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,2)
plot(dat1W2(:,1), dat1W2(:,2), 'k', dat1B2(:,1), dat1B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick)
grid on

subplot(2,2,3)
plot(dat1W4(:,1), dat1W4(:,2), 'k', dat1B4(:,1), dat1B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,4)
plot(dat1W6(:,1), dat1W6(:,2), 'k', dat1B6(:,1), dat1B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_m2_Oh01_eta04.tif -r256
saveas(gcf, 'ModeDecomp_m2_Oh01_eta04', 'epsc')


%%

% case 2
figure
subplot(2,2,1)
plot(dat2W0(:,1), dat2W0(:,2), 'k', dat2B0(:,1), dat2B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,2)
plot(dat2W2(:,1), dat2W2(:,2), 'k', dat2B2(:,1), dat2B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick)
grid on

subplot(2,2,3)
plot(dat2W4(:,1), dat2W4(:,2), 'k', dat2B4(:,1), dat2B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,4)
plot(dat2W6(:,1), dat2W6(:,2), 'k', dat2B6(:,1), dat2B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_m2_Oh01_eta02.tif -r256
saveas(gcf, 'ModeDecomp_m2_Oh01_eta02', 'epsc')


%%

% case 3
figure
subplot(2,2,1)
plot(dat3W0(:,1), dat3W0(:,2), 'k', dat3B0(:,1), dat3B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,2)
plot(dat3W2(:,1), dat3W2(:,2), 'k', dat3B2(:,1), dat3B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick)
grid on

subplot(2,2,3)
plot(dat3W4(:,1), dat3W4(:,2), 'k', dat3B4(:,1), dat3B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,4)
plot(dat3W6(:,1), dat3W6(:,2), 'k', dat3B6(:,1), dat3B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_m2_Oh01_eta01.tif -r256
saveas(gcf, 'ModeDecomp_m2_Oh01_eta01', 'epsc')


%%

% % case 4
% 
% % equal modes
% figure
% subplot(2,2,1)
% plot(dat2W0(:,1), dat2W0(:,2), 'k', dat2B0(:,1), dat2B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4]) 
% grid on
% 
% subplot(2,2,2)
% plot(dat2W2(:,1), dat2W2(:,2), 'k', dat2B2(:,1), dat2B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4])
% grid on
% 
% subplot(2,2,3)
% plot(dat2W4(:,1), dat2W4(:,2), 'k', dat2B4(:,1), dat2B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4]) 
% grid on
% 
% subplot(2,2,4)
% plot(dat2W6(:,1), dat2W6(:,2), 'k', dat2B6(:,1), dat2B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4]) 
% grid on
% 
% set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])
% 
% 
% %export_fig ModeDecomp_case2equal.tif -r256
% 
% 
% % unequal modes
% figure
% subplot(2,2,1)
% plot(dat2W1(:,1), dat2W1(:,2), 'k', dat2B1(:,1), dat2B1(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(1, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4]) 
% grid on
% 
% subplot(2,2,2)
% plot(dat2W3(:,1), dat2W3(:,2), 'k', dat2B3(:,1), dat2B3(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(3, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4])
% grid on
% 
% subplot(2,2,3)
% plot(dat2W5(:,1), dat2W5(:,2), 'k', dat2B5(:,1), dat2B5(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(5, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4]) 
% grid on
% 
% subplot(2,2,4)
% plot(dat2W7(:,1), dat2W7(:,2), 'k', dat2B7(:,1), dat2B7(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
% legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
% xlim([0,4])
% ylabel('$(7, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
% set(gca, 'XTick', [0 1 2 3 4]) 
% grid on
% 
% set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])
% 
% 
% %export_fig ModeDecomp_case2unequal.tif -r256


%%

% case 7
figure
subplot(2,2,1)
plot(dat7W0(:,1), dat7W0(:,2), 'k', dat7B0(:,1), dat7B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,2)
plot(dat7W2(:,1), dat7W2(:,2), 'k', dat7B2(:,1), dat7B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick)
grid on

subplot(2,2,3)
plot(dat7W4(:,1), dat7W4(:,2), 'k', dat7B4(:,1), dat7B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,4)
plot(dat7W6(:,1), dat7W6(:,2), 'k', dat7B6(:,1), dat7B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_m4_Oh01_eta04.tif -r256
saveas(gcf, 'ModeDecomp_m4_Oh01_eta04', 'epsc')


%%

% case 8
figure
subplot(2,2,1)
plot(dat8W0(:,1), dat8W0(:,2), 'k', dat8B0(:,1), dat8B0(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,2)
plot(dat8W2(:,1), dat8W2(:,2), 'k', dat8B2(:,1), dat8B2(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick)
grid on

subplot(2,2,3)
plot(dat8W4(:,1), dat8W4(:,2), 'k', dat8B4(:,1), dat8B4(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

subplot(2,2,4)
plot(dat8W6(:,1), dat8W6(:,2), 'k', dat8B6(:,1), dat8B6(:,2), 'b', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_m4_Oh01_eta01.tif -r256
saveas(gcf, 'ModeDecomp_m4_Oh01_eta01', 'epsc')


%%

% case 9
figure
subplot(2,2,1)
plot(dat9W0(:,1), dat9W0(:,2), 'k', dat9B0(:,1), dat9B0(:,2), 'b', dat9B02(:,1), dat9B02(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(0, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

subplot(2,2,2)
plot(dat9W2(:,1), dat9W2(:,2), 'k', dat9B2(:,1), dat9B2(:,2), 'b', dat9B22(:,1), dat9B22(:,2), 'b--','Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(2, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4])
grid on

subplot(2,2,3)
plot(dat9W4(:,1), dat9W4(:,2), 'k', dat9B4(:,1), dat9B4(:,2), 'b', dat9B42(:,1), dat9B42(:,2), 'b--','Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(4, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

subplot(2,2,4)
plot(dat9W6(:,1), dat9W6(:,2), 'k', dat9B6(:,1), dat9B6(:,2), 'b', dat9B62(:,1), dat9B62(:,2), 'b--', 'Linewidth', 1.2, 'MarkerSize', size_marker)
legend({'WNLT', 'BoSSS - zero Init', 'BoSSS - 3. Order Init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('$(6, 0)$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4]) 
grid on

set(gcf, 'color', 'w', 'Position', [300, 300, 1600, 900])


%export_fig ModeDecomp_m4_Oh056_eta005.tif -r256
saveas(gcf, 'ModeDecomp_m4_Oh056_eta005', 'epsc')

