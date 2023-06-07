clear all

%% import data

path = '';  %'plotData\Energies'; % path to data folder     
          
% 1: m = 2 Oh = 0.1, eta0 = 0.4;   
path = 'm2/m2_Oh01_eta04';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat1kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat1surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat1totE = import.data;


% 2: m = 2 Oh = 0.1, eta0 = 0.2;   
path = 'm2/m2_Oh01_eta02';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat2kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat2surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat2totE = import.data;


% 3: m = 2 Oh = 0.1, eta0 = 0.1;   
path = 'm2/m2_Oh01_eta01';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat3kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat3surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat3totE = import.data;


% 4: m = 3 Oh = 0.1, eta0 = 0.4;   
path = 'm3/m3_Oh01_eta04';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat4kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat4surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat4totE = import.data;


% 5: m = 3 Oh = 0.1, eta0 = 0.3;   
path = 'm3/m3_Oh01_eta03';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat5kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat5surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat5totE = import.data;


% 6: m = 3 Oh = 0.1, eta0 = 0.15;   
path = 'm3/m3_Oh01_eta015';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat6kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat6surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat6totE = import.data;


% 7: m = 4 Oh = 0.1, eta0 = 0.4;   
path = 'm4/m4_Oh01_eta04';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat7kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat7surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat7totE = import.data;


% 8: m = 4 Oh = 0.1, eta0 = 0.1;   
path = 'm4/m4_Oh01_eta01';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat8kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat8surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat8totE = import.data;


% 9: m = 2 Oh = 0.56, eta0 = 0.05;   
path = 'm4/m4_Oh056_eta005';
datapath = strcat(path, ' kineticEnergy J686.txt');
import = importdata(datapath);
dat9kinE = import.data;
datapath = strcat(path, ' surfaceEnergy J686.txt');
import = importdata(datapath);
dat9surfE = import.data;
datapath = strcat(path, ' totalEnergy J686.txt');
import = importdata(datapath);
dat9totE = import.data;

datapath = strcat(path, ' kineticEnergy J686-3OrdInit.txt');
import = importdata(datapath);
dat9kinE2 = import.data;
datapath = strcat(path, ' surfaceEnergy J686-3OrdInit.txt');
import = importdata(datapath);
dat9surfE2 = import.data;
datapath = strcat(path, ' totalEnergy J686-3OrdInit.txt');
import = importdata(datapath);
dat9totE2 = import.data;


%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;
mInt = 100; % marker interval for plotting


%%

% case 1
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1kinE(:,1), dat1kinE(:,2), 'r', dat1surfE(:,1), dat1surfE(:,2), 'b', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'kinetic energy', 'surface energy'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig kinSurfE_m2_Oh01_eta04.tif -r256

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1totE(:,1), dat1totE(:,2), 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'total energy'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig totE_m2_Oh01_eta04.tif -r256


%%

% case 1-3
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1kinE(:,1), dat1kinE(:,2), 'r-.', dat1surfE(:,1), dat1surfE(:,2), 'b-.', dat2kinE(:,1), dat2kinE(:,2), 'r--', dat2surfE(:,1), dat2surfE(:,2), 'b--', dat3kinE(:,1), dat3kinE(:,2), 'r', dat3surfE(:,1), dat3surfE(:,2), 'b', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = \{0.4, 0.2, 0.1\}$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'kinetic energy - $\eta_0 = 0.4$', 'surface energy - $\eta_0 = 0.4$', 'kinetic energy - $\eta_0 = 0.2$', 'surface energy - $\eta_0 = 0.2$', 'kinetic energy - $\eta_0 = 0.1$', 'surface energy - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig kinSurfE_m2_Oh01.png  %.tif -r256

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat1totE(:,1), dat1totE(:,2), '-.', dat2totE(:,1), dat2totE(:,2), '--', dat3totE(:,1), dat3totE(:,2), 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=2$ ${\rm Oh}=0.1$ $\eta_0 = \{0.4, 0.2, 0.1\}$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'total energy - $\eta_0 = 0.4$', 'total energy - $\eta_0 = 0.2$', 'total energy - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig totE_m2_Oh01.png  %.tif -r256


%%

% case 4-6
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat4kinE(:,1), dat4kinE(:,2), 'r-.', dat4surfE(:,1), dat4surfE(:,2), 'b-.', dat5kinE(:,1), dat5kinE(:,2), 'r--', dat5surfE(:,1), dat5surfE(:,2), 'b--', dat6kinE(:,1), dat6kinE(:,2), 'r', dat6surfE(:,1), dat6surfE(:,2), 'b', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=3$ ${\rm Oh}=0.1$ $\eta_0 = \{0.4, 0.3, 0.15\}$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'kinetic energy - $\eta_0 = 0.4$', 'surface energy - $\eta_0 = 0.4$', 'kinetic energy - $\eta_0 = 0.2$', 'surface energy - $\eta_0 = 0.2$', 'kinetic energy - $\eta_0 = 0.1$', 'surface energy - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,2])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig kinSurfE_m3_Oh01.png  %.tif -r256

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat4totE(:,1), dat4totE(:,2), '-.', dat5totE(:,1), dat5totE(:,2), '--', dat6totE(:,1), dat6totE(:,2), 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=3$ ${\rm Oh}=0.1$ $\eta_0 = \{0.4, 0.3, 0.15\}$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'total energy - $\eta_0 = 0.4$', 'total energy - $\eta_0 = 0.2$', 'total energy - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,2])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig totE_m3_Oh01.png  %.tif -r256


%%

% case 7
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat7kinE(:,1), dat7kinE(:,2), 'r', dat7surfE(:,1), dat7surfE(:,2), 'b', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'kinetic energy', 'surface energy'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig kinSurfE_m4_Oh01_eta04.tif -r256

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat7totE(:,1), dat7totE(:,2), 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$ $\eta_0 = 0.4$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'total energy'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig totE_m4_Oh01_eta04.tif -r256


%%

% case 7,8
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat7kinE(:,1), dat7kinE(:,2), 'r--', dat7surfE(:,1), dat7surfE(:,2), 'b--', dat8kinE(:,1), dat8kinE(:,2), 'r', dat8surfE(:,1), dat8surfE(:,2), 'b', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$ $\eta_0 = \{0.4, 0.1\}$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'kinetic energy - $\eta_0 = 0.4$', 'surface energy - $\eta_0 = 0.4$', 'kinetic energy - $\eta_0 = 0.1$', 'surface energy - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig kinSurfE_m4_Oh01.png  %.tif -r256

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat7totE(:,1), dat7totE(:,2), '--', dat8totE(:,1), dat8totE(:,2), 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.1$ $\eta_0 = \{0.4, 0.1\}$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'total energy - $\eta_0 = 0.4$', 'total energy - $\eta_0 = 0.1$'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,7])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig totE_m4_Oh01.png  %.tif -r256


%%

% case 9
figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat9kinE(:,1), dat9kinE(:,2), 'r', dat9surfE(:,1), dat9surfE(:,2), 'b', dat9kinE2(:,1), dat9kinE2(:,2), 'r--', dat9surfE2(:,1), dat9surfE2(:,2), 'b--','Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.56$ $\eta_0 = 0.05$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'kinetic energy - zero init', 'surface energy - zero init', 'kinetic energy - 3. order init', 'surface energy - 3. order init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig kinSurfE_m4_Oh056_eta005.png  %.tif -r256

figure
set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(dat9totE(:,1), dat9totE(:,2), dat9totE2(:,1), dat9totE2(:,2), '--', 'Linewidth', 1.5, 'MarkerSize', size_marker)
title('$m=4$ ${\rm Oh}=0.56$ $\eta_0 = 0.05$', 'Interpreter', 'latex', 'Fontsize', size_legend)
legend({'total energy - zero init', 'total energy - 3. order init'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
xlim([0,4])
ylabel('energy', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'XTick', [0 1 2 3 4], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig totE_m4_Oh056_eta005.png  %.tif -r256
