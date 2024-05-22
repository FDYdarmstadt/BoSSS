clear all

%% import data

data_CoM = cell(1,1);
% 
%import = importdata("DropletRebound_8x8x8AMR1_k3_dropVelocity100vH_NeKs_withReInit10_c4t4_CenterOfMassZ");
import = importdata("DropletRebound_8x8x8AMR1_k3_dropVelocity100vH_NeKs_RI1_halfdt_c4t4_restart4_CenterOfMassZ");
data_CoM{1,1} = import.data;


import = importdata("Posovertime.xlsx");
data_ref = import.data;


%% 

pos0 = data_CoM{1,1}(1,2)
idx = 0;
for i=1:1:length(data_ref(:,1))
    if data_ref(i,3) < pos0
        idx = i;
        break
    end
end
m = (data_ref(idx,3) - data_ref(idx-1,3))/(data_ref(idx,1) - data_ref(idx-1,1))
t0 = data_ref(idx-1,1) + ((pos0 - data_ref(idx-1,3)) / m);


%% plot

size_legend = 14;
size_label = 16;
size_tick = 12;
size_marker = 6;

%%
figure
plot(data_CoM{1,1}(:,1), data_CoM{1,1}(:,2), 'b', data_ref(:,1)-t0, data_ref(:,3), 'k')
legend({'$u_{D,I} = -1.12 \ (100\%)$', 'ref data experiment'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'southeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
%ylim([0.2e-3 0.8e-3])
xlim([0e-4 12e-4])
%set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig dropletRebound_centerOfMassZ.tif -r256


%%
figure
plot(data_CoM{1,1}(:,1), data_CoM{1,1}(:,2), 'b', data_ref(:,1)-t0, data_ref(:,3), 'k')
legend({'$u_{D,I} = -1.12 \ (100\%)$', 'ref data experiment'}, 'Interpreter', 'latex', 'Fontsize', size_legend, 'Location', 'northeast')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
%ylim([0.6e-3 0.8e-3])
xlim([0 4e-4])
%set(gca, 'XTick', [0 1 2 3 4 5 6 7], 'FontSize', size_tick) 
%set(gca, 'YTick', [0.02 0.03 0.04 0.05 0.06]) 
set(gcf, 'color', 'w')
grid on

export_fig dropletRebound_centerOfMass_CloseUp2.tif -r256


