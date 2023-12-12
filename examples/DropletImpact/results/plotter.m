clear all

%% import data

data_CoM = cell(3,1);

import = importdata("DropletRebound_8x8x16AMR1_dropVelocity0_noGravity_CenterOfMassZ");
data_CoM{1,1} = import.data;
import = importdata("DropletRebound_8x8x16AMR1_dropVelocity0_noGravity_restart_CenterOfMassZ");
data_CoM{1,1} = [data_CoM{1}; import.data];

import = importdata("DropletRebound_8x8x8AMR1_dropVelocity10vH_noGravity_CenterOfMassZ");
data_CoM{2,1} = import.data;

import = importdata("DropletRebound_8x8x8AMR1_dropVelocity25vH_noGravity_CenterOfMassZ");
data_CoM{3,1} = import.data;

import = importdata("DropletRebound_8x8x8AMR1_dropVelocity50vH_noGravity_CenterOfMassZ");
data_CoM{4,1} = import.data;

import = importdata("DropletRebound_8x8x8AMR1_dropVelocity100vH_noGravity_CenterOfMassZ");
data_CoM{5,1} = import.data;

import = importdata("DropletRebound_8x8x8AMR1_dropVelocity0_CenterOfMassZ");
data_CoM{6,1} = import.data;



% data_mV = cell(3,1);
% 
% import = importdata("DropletRebound_8x8x16AMR1_dropVelocity0_noGravity_meanVelocityZ");
% data_mV{1,1} = import.data;
% import = importdata("DropletRebound_8x8x16AMR1_dropVelocity0_noGravity_restart_meanVelocityZ");
% data_mV{1,1} = [data_mV{1}; import.data];
% 
% import = importdata("DropletRebound_8x8x8AMR1_dropVelocity10vH_noGravity_meanVelocityZ");
% data_mV{2,1} = import.data;
% 
% import = importdata("DropletRebound_8x8x8AMR1_dropVelocity25vH_noGravity_meanVelocityZ");
% data_mV{3,1} = import.data;



import = importdata("Posovertime.xlsx");
data_ref = import.data;


%% 

pos0 = data_CoM{1,1}(1,2)
% idx = 0;
% for i=1:1:length(data_ref(:,1))
%     if data_ref(i,3) < pos0
%         idx = i;
%         break
%     end
% end
m = (data_ref(79,3) - data_ref(78,3))/(data_ref(79,1) - data_ref(78,1))
t0 = data_ref(78,1) + ((pos0 - data_ref(78,3)) / m);


%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;

%%
figure
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(data_CoM{1,1}(:,1), data_CoM{1,1}(:,2), 'm', data_CoM{2,1}(:,1), data_CoM{2,1}(:,2), 'c', data_CoM{3,1}(:,1), data_CoM{3,1}(:,2), 'r', data_CoM{4,1}(:,1), data_CoM{4,1}(:,2), 'g', data_CoM{5,1}(:,1), data_CoM{5,1}(:,2), 'b', data_CoM{6,1}(:,1), data_CoM{6,1}(:,2), '--m', data_ref(:,1)-t0, data_ref(:,3), 'k')
%legend({'$u_{D,I} = 0.0$', '$u_{D,I} = -0.112 \ (10\%)$', '$u_{D,I} = -0.28 \ (25\%)$', '$u_{D,I} = -0.56 \ (50\%)$', '$u_{D,I} = -1.12 \ (100\%)$', '$u_{D,I} = 0.0 \ (g \neq 0)$', 'ref data experiment'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
ylim([0.2e-3 1.2e-3])
%xlim([-1e-4 8e-4])
%set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig dropletRebound_centerOfMass.tif -r256


%%
figure
%set(gcf,'DefaultAxesColorOrder',[0.9290 0.6940 0.1250])
plot(data_CoM{1,1}(:,1), data_CoM{1,1}(:,2), 'm', data_CoM{2,1}(:,1), data_CoM{2,1}(:,2), 'c', data_CoM{3,1}(:,1), data_CoM{3,1}(:,2), 'r', data_CoM{4,1}(:,1), data_CoM{4,1}(:,2), 'g', data_CoM{5,1}(:,1), data_CoM{5,1}(:,2), 'b', data_CoM{6,1}(:,1), data_CoM{6,1}(:,2), '--m', data_ref(:,1)-t0, data_ref(:,3), 'k')
legend({'$u_{D,I} = 0.0$', '$u_{D,I} = -0.112 \ (10\%)$', '$u_{D,I} = -0.28 \ (25\%)$', '$u_{D,I} = -0.56 \ (50\%)$', '$u_{D,I} = -1.12 \ (100\%)$', '$u_{D,I} = 0.0 \ (g \neq 0)$', 'ref data experiment'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
%ylim([0.6e-3 0.8e-3])
xlim([0 1e-4])
%set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig dropletRebound_centerOfMass_CloseUp.tif -r256


%%
% figure
% plot(data_mV{1,1}(:,1), data_mV{1,1}(:,2), data_mV{2,1}(:,1), data_mV{2,1}(:,2), data_mV{3,1}(:,1), data_mV{3,1}(:,2))
% 
% figure
% plot(data_ref(:,1), data_ref(:,3))


