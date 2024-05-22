close all

%load data
% supersonic fast sinus
file_path1 = 'AW_p3_xCells61_tCells61_tMax1.5_1sinus_Mach1.5_sP1.5_wP0_ampneg0_amppos1E-05_wLTrue_wSTrue_wL0.8_p_per';

% subsonic slow sinus
file_path2 = 'AW_p3_xCells61_tCells61_tMax4_1sinus_Mach1.5_sP0.5_wP0.9_ampneg1E-05_amppos0_wLTrue_wSTrue_wL0.8_p_per';


% CNS supersonic fast sinus
file_path1 = '..\\_dat\\CompShockCap\\AW_p3_xCells61_yCells3_sP1.5_pST10_wP-0.8_ampneg0_amppos1E-05_wL0.8_Mach1.5_1sinus_RK4_p_per';
% subsonic slow sinus
file_path2 = '..\\_dat\\CompShockCap\\AW_p3_xCells61_yCells3_sP0.5_pST10_wP3.2_ampneg1E-05_amppos0_wL0.8_Mach1.5_1sinus_RK4_p_per';

solver="CNS";

% Read data from both datasets to determine z-axis limits
data1 = readmatrix(strcat(file_path1,".csv"));
data2 = readmatrix(strcat(file_path2,".csv"));


% Given variables
MBaseL = 1.5; % Base Mach number on the left
gamma = 1.4; % Specific heat ratio

% Calculate right Mach number (Mach_R)
MBaseR  = sqrt((1 + ((gamma-1)/2) * MBaseL^2) / (gamma * MBaseL^2 - (gamma-1)/2));

% Upstream fast acoustic wave amplification coefficient
delta_pPrR_over_delta_pPrL = ((1+MBaseL)^2) / (1+2*MBaseR+1/MBaseL^2) * ...
    (1 - ((gamma-1)/(gamma+1)) * (1-(1/MBaseL))^2);

% Downstream slow acoustic wave amplification coefficient
delta_pPrR_over_delta_pPrMiR = -(1- 2*MBaseR +1 / MBaseL^2) / (1+ 2*MBaseR +1 / MBaseL^2);

% Display results
fprintf('Upstream fast acoustic wave amplification coefficient: %f\n', delta_pPrR_over_delta_pPrL);
fprintf('Downstream slow acoustic wave amplification coefficient: %f\n', delta_pPrR_over_delta_pPrMiR);

% file_path1 = '..\\_dat\\AcousticShockInteraction\\AW_p3_xCells61_tCells61_tMax1.5_1sinus_Mach1.5_sP1.5_wP0_ampneg0_amppos1E-05_wLTrue_wSTrue_wL0.8_p_per.csv';
% data1 = readmatrix(file_path1);
% solVals1 = data1(2:1:end, 2:1:end);
% fprintf('Upstream fast acoustic wave max amplification: %f\n', max(max(solVals1))/1e-5);
% file_path2 = '..\\_dat\\AcousticShockInteraction\\AW_p3_xCells61_tCells61_tMax4_1sinus_Mach1.5_sP0.5_wP0.9_ampneg1E-05_amppos0_wLTrue_wSTrue_wL0.8_p_per.csv';
data2 = readmatrix(file_path2);
solVals2 = data2(2:1:end, 2:1:end);
 fprintf('Downstream slow acoustic wave max amplification: %f\n', max(solVals2)/1e-5);

figure
start = 60; endn=90;
Xmin=start/100*3.0;Xmax=endn/100*3.0;
maxPoints=max(data1( 2:end,start:endn))';
X=linspace(Xmin, Xmax, length(maxPoints));
% Assuming X and maxPoints are already defined and have the same length
dataToSave = [X', maxPoints./1e-5];  % Combine X and maxPoints into a single matrix
dataToSave2 = [[Xmin;Xmax],  2.154926*ones(2,1)];
% Write the matrix to a text file
writematrix(dataToSave, strcat(file_path1,'_MaxPoints.txt'), 'Delimiter', 'tab');
writematrix(dataToSave2, strcat(file_path1,'_Optimal.txt'), 'Delimiter', 'tab');
plot(X,maxPoints'./1e-5)

figure
start = 60; endn=99; 
Xmin=start/100*3.0;Xmax=endn/100*3.0;
minPoints=min(data2(2:end, start:1:endn))';
X=linspace(Xmin, Xmax, length(minPoints));
% Assuming X and maxPoints are already defined and have the same length
dataToSave = [X', minPoints./1e-5];  % Combine X and maxPoints into a single matrix
dataToSave2 = [[Xmin;Xmax], -0.014848*ones(2,1)];
% Write the matrix to a text file
writematrix(dataToSave, strcat(file_path2,'_MinPoints.txt'), 'Delimiter', 'tab');
writematrix(dataToSave2, strcat(file_path2,'_Optimal.txt'), 'Delimiter', 'tab');
plot(X,minPoints'./1e-5)
