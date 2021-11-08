close all;clear; clc;

%% Problem geometry description
C.L = 0.02; % Domain length, m
C.cp = 1.3; % Heat capacity;
C.initialCellNumber = 100;
useBoSSSForInitialization = true;
C.chemActive = true;
C.variableKineticParameters = true;

velocityMultiplier = 11;

%% Boundary conditions
if useBoSSSForInitialization
    % Read bosss data
    path = ['C:\tmp\BoSSS_data_VariableChemParams_cpmixture_nonunityLe\', num2str(velocityMultiplier), '\'];
    bosss_x = load([path, 'FullXCoord.txt']);
    bosssVelocityX = load([path, 'FullVelocityX.txt']);
    bosssTemperature = load([path, 'FullTemperature.txt']);
    bosssY0 = load([path, 'FullMassFraction0.txt']);
    bosssY1 = load([path, 'FullMassFraction1.txt']);
    bosssY2 = load([path, 'FullMassFraction2.txt']);
    bosssY3 = load([path, 'FullMassFraction3.txt']);
    bosssY4 = ones(1, length(bosssY0)) - (bosssY0 + bosssY1 + bosssY2 + bosssY3);


    C.vLeft = bosssVelocityX(1); % Velocity left (fuel)
    C.vRight = bosssVelocityX(end); % Velocity right (oxidizer)
    C.TF0 = bosssTemperature(1); % Temperature left (fuel)
    C.TO0 = bosssTemperature(2); % Temperature right (oxidizer)
    C.fuelInletConcentration = [bosssY0(1), bosssY1(1), bosssY2(1), bosssY3(1), bosssY4(1)];
    C.oxidizerInletConcentration = [bosssY0(end), bosssY1(end), bosssY2(end), bosssY3(end), bosssY4(end)];
else
    C.vLeft = 0.0243 * velocityMultiplier; % Velocity left (fuel)
    C.vRight = -0.0243 * velocityMultiplier * 3; % Velocity right (oxidizer)
    C.TF0 = 300; % Temperature left (fuel)
    C.TO0 = 300; % Temperature right (oxidizer)
    C.fuelInletConcentration = [0.2, 0.0, 0.0, 0.0, 0.8];
    C.oxidizerInletConcentration = [0.0, 0.23, 0.0, 0.0, 0.77];

end

%% Chemistry
C.p0 = 101325; % Pressure, Pa
C.viscosity0 = 1.716e-5; % kg/( m s) ==> viscosity at T = 273.15 for air
C.Tref = 273; % Reference temperature from powerlaw
C.a = 2 / 3; % PowerLaw exponent
C.Coef_Stoic = [-1, -2, 1, 2, 0]; % Chemical reaction 1CH4 + 2O2 -> 1CO2 + 2H2O
C.MM = [16, 32, 44, 18, 28];
C.s = (C.Coef_Stoic (2) * C.MM(2)) / (C.Coef_Stoic(1) * C.MM(1));

C.phi = C.s * C.fuelInletConcentration(1) / C.oxidizerInletConcentration(2);
C.zst = 1.0 / (1.0 + C.phi);
C.Q = 50100;
C.R = 8.314 * 1000; % gas constant, kg m^2 / (s^2 K mol)
%% Calculate flame sheet (infinite reaction rate)
x = zeros(1);
newSolution = zeros(1);
lambdaOut = 100;
if C.chemActive
    [mySolution, lambdaOut] = CounterFlowFlame_MixtureFraction(C);
    fprintf('Finished calculation of MF\n');
    %         plotMixtureFractionSolution(mySolution,C);
    x = mySolution(1, :);
    myzeros = zeros(1, length(x)); %% initial values for the derivatives not known, 
    newSolution = [mySolution(2, :); mySolution(3, :); mySolution(4, :); mySolution(5, :); ...
        myzeros; mySolution(6, :); myzeros; mySolution(7, :); myzeros; mySolution(8, :); ...
        myzeros; mySolution(9, :); myzeros; mySolution(10, :); myzeros];
end

%% Full Problem, first strain
[xint, Sxint, calculatedStrain, calculatedMaxTemperature] = CounterFlowFlame_FullChemistry(C, newSolution, x, lambdaOut);

% plotFullProblem(xint,Sxint);

figure
hold on
plot(xint, Sxint(4, :));
hold on
plot(mySolution(1, :), mySolution(5, :));
legend('')

%% Compare Matlab with BoSSS Solution

if useBoSSSForInitialization
    %bosss_x = bosss_x + 0.02; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    %VelocityX
    subplot(3, 3, 1)
    plot(xint, Sxint(1, :))

    hold on;
    plot(bosss_x, bosssVelocityX)
    legend('VelocityXMatlab', 'VelocityXBoSSS')
    xlim([xint(1), xint(end)])
    ylim auto
    title('VelocityX')
    xlabel('x')
    ylabel('u')
    hold off;

    %Temperature
    subplot(3, 3, 2)
    plot(xint, Sxint(4, :))
    hold on;
    plot(bosss_x, bosssTemperature)
    legend('Matlab', 'BoSSS')
    xlim([xint(1), xint(end)])
    ylim auto
    title('Temperature')
    xlabel('x')
    ylabel('T')

    hold off;

    %MassFractions
    subplot(3, 3, 3)
    plot(xint, Sxint(6, :))
    hold on;
    plot(bosss_x, bosssY0)
    legend('Matlab', 'BoSSS')
    xlim([xint(1), xint(end)])
    ylim auto
    title('MassFraction 0')
    xlabel('x')
    ylabel('Y0')

    %MassFractions
    subplot(3, 3, 4)
    plot(xint, Sxint(8, :))
    hold on;
    plot(bosss_x, bosssY1)
    legend('Matlab', 'BoSSS')
    xlim([xint(1), xint(end)])
    ylim auto
    title('MassFraction 1')
    xlabel('x')
    ylabel('Y1')
    hold off;


    %MassFractions
    subplot(3, 3, 5)
    plot(xint, Sxint(10, :))
    hold on;
    plot(bosss_x, bosssY2)
    legend('Matlab', 'BoSSS')
    xlim([xint(1), xint(end)])
    ylim auto
    title('MassFraction CO2')
    xlabel('x')
    ylabel('Y4')
    hold off;

    %MassFractions
    subplot(3, 3, 6)
    plot(xint, Sxint(14, :))
    hold on;
    plot(bosss_x, bosssY4)
    legend('Matlab', 'BoSSS')
    xlim([xint(1), xint(end)])
    ylim auto
    title('MassFraction N2')
    xlabel('x')
    ylabel('Y4')
    hold off;


    for k = 1:length(xint)
        phi(k) = GetPhi(Sxint(6, k), Sxint(8, k));
        if phi(k) > 1.5
            phi(k) = 1.5;
        end
    end

    subplot(3, 3, 7)
    plot(xint, phi)
    legend('Matlab')
    xlim([xint(1), xint(end)])
    ylim auto
    title('Phi, local')
    xlabel('x')
    ylabel('phi')
    hold off;

    %% Reaction rate

    for k = 1:length(xint)
        Y0 = Sxint(6, k);
        Y1 = Sxint(8, k);
        Y2 = Sxint(10, k);
        Y3 = Sxint(12, k);
        Y4 = Sxint(14, k);
        T = Sxint(4, k);


        density = 101325 / (8314 * T * (Y0 / 16 + Y1 / 32 + Y2 / 44 + Y3 / 18 + Y4 / 28));
        reacRate(k) = 6.9e11 * exp(-15900/T) * density * density * Y0 * Y1 / (16 * 32);
    end

    subplot(3, 3, 8)
    plot(xint, reacRate)
    legend('Matlab')
    xlim([xint(1), xint(end)])
    ylim auto
    title('ReactionRate')
    xlabel('x')
    ylabel('Y0Y1')
    hold off;

end

%% Save matlab result to textfile

    writetable(table(xint', Sxint(1, :)'), fullfile(path, 'FullVelocityX_MATLAB.txt'), 'Delimiter', ' ')
    writetable(table(xint', Sxint(4, :)'), fullfile(path, 'FullTemperature_MATLAB.txt'), 'Delimiter', ' ')
    writetable(table(xint', Sxint(6, :)'), fullfile(path, 'FullMassFraction0_MATLAB.txt'), 'Delimiter', ' ')
    writetable(table(xint', Sxint(8, :)'), fullfile(path, 'FullMassFraction1_MATLAB.txt'), 'Delimiter', ' ')
    writetable(table(xint', Sxint(10, :)'), fullfile(path, 'FullMassFraction2_MATLAB.txt'), 'Delimiter', ' ')
    writetable(table(xint', Sxint(12, :)'), fullfile(path, 'FullMassFraction3_MATLAB.txt'), 'Delimiter', ' ')

    writetable(table(bosss_x', bosssVelocityX'), fullfile(path, 'FullVelocityX_BOSSS.txt'), 'Delimiter', ' ')
    writetable(table(bosss_x', bosssTemperature'), fullfile(path, 'FullTemperature_BOSSS.txt'), 'Delimiter', ' ')
    writetable(table(bosss_x', bosssY0'), fullfile(path, 'FullMassFraction0_BOSSS.txt'), 'Delimiter', ' ')
    writetable(table(bosss_x', bosssY1'), fullfile(path, 'FullMassFraction1_BOSSS.txt'), 'Delimiter', ' ')
    writetable(table(bosss_x', bosssY2'), fullfile(path, 'FullMassFraction2_BOSSS.txt'), 'Delimiter', ' ')
    writetable(table(bosss_x', bosssY3'), fullfile(path, 'FullMassFraction3_BOSSS.txt'), 'Delimiter', ' ')


function viscosity = mu(T)

viscosity0 = 1.716e-5; % kg/( m s) ==> viscosity at T = 273.15 for air
viscosity_constant = viscosity0;
Tref = 273;
a = 2 / 3;
viscosity_variable = viscosity0 * (T / Tref)^(a);
viscosity = viscosity_constant * (1 - multiplier) + viscosity_variable * multiplier;
end

function phi = GetPhi(YF, YO)
s = 4;
YF0 = 0.2;
YO0 = 0.23;
phi = (s * YF0 / YO0) * (s * YF - YO + YO0) / (s * (YF0 - YF) + YO);
end
