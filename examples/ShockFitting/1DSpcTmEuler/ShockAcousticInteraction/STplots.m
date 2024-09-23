% Specify the file paths for the two datasets
close all
% supersonic fast sinus
file_path1 = 'AW_p3_xCells61_tCells61_tMax1.5_1sinus_Mach1.5_sP1.5_wP0_ampneg0_amppos1E-05_wLTrue_wSTrue_wL0.8_p_per.csv';
file_path2 = 'AW_p3_xCells61_tCells61_tMax1.5_1sinus_Mach1.5_sP1.5_wP0_ampneg0_amppos1E-05_wLTrue_wSTrue_wL0.8_rho_per.csv';
export_name='supersonic_fast_sinus';

% supersonic slow sinus
file_path1 = 'AW_p3_xCells61_tCells61_tMax4_1sinus_Mach1.5_sP1.5_wP0.3_ampneg1E-05_amppos0_wLTrue_wSTrue_wL0.8_p_per.csv';
file_path2 = 'AW_p3_xCells61_tCells61_tMax4_1sinus_Mach1.5_sP1.5_wP0.3_ampneg1E-05_amppos0_wLTrue_wSTrue_wL0.8_rho_per.csv';
export_name='supersonic_slow_sinus';

% subsonic slow sinus
file_path1 = 'AW_p3_xCells61_tCells61_tMax4_1sinus_Mach1.5_sP0.5_wP0.9_ampneg1E-05_amppos0_wLTrue_wSTrue_wL0.8_p_per.csv';
file_path2 = 'AW_p3_xCells61_tCells61_tMax4_1sinus_Mach1.5_sP0.5_wP0.9_ampneg1E-05_amppos0_wLTrue_wSTrue_wL0.8_rho_per.csv';
export_name='subsonic_slow_sinus';



% Read data from both datasets to determine z-axis limits
data1 = readmatrix(file_path1);
data2 = readmatrix(file_path2);

tStart=300;
tEnd=600;
dT=1;dX=1;
solVals1 = data1(tStart:dT:tEnd, 2:dX:end);
solVals2 = data2(2:dT:end, 2:dX:end);

plot(data1(1, 2:dX:end),solVals1(2,1:dX:end));
%zlim_values = [min(min(solVals1, [], 'all'), min(solVals2, [], 'all')), ...
%               max(max(solVals1, [], 'all'), max(solVals2, [], 'all'))];
zlim_values1 = [min(solVals1, [], 'all'),max(solVals1, [], 'all')];
zlim_values2 = [min(solVals2, [], 'all'),max(solVals2, [], 'all')];
% Create a common figure
%figure;

% Plot the first subplot
plotWaterfall(file_path1, 1, '', zlim_values1,true);

% Plot the second subplot
plotWaterfall(file_path2, 2, '', zlim_values2,false);

% Save the figure
%export_fig(strcat('../figures/',export_name,'.pdf'), '-pdf');


function plotWaterfall(file_path, subplot_index, subplot_title, zlim_values,is_pressure)
    % Read data from the CSV file
    data = readmatrix(file_path);
    
    % Assuming your CSV file has two columns and you want to convert them to matrices
    dX = 1; dT = 1;
    xVals = data(1, 2:dX:end); 
    tVals = data(2:dT:end, 1); 
    solVals = data(2:dT:end, 2:dX:end);

    % Create subplot
    subplot(1, 2, subplot_index);
    waterfall(xVals, tVals, solVals);
    if(is_pressure)
        zlabel('p^\prime');
    else
        zlabel('\rho^\prime');
    end
    xlabel('x');
    ylabel('t');
    title(subplot_title);
    
    % Set z-axis limits
    zlim(zlim_values);
    clim(zlim_values);


    set(gcf, 'Color', 'w');
    azimuth = 20;  % Azimuthal angle
    elevation = 20;  % Elevation angle
    view(azimuth, elevation);
end