%%Initial solution is [v,U1,U2,TO0,Y1_1,Y2_1,Y3_1,Y4_1,Y5_1]
%% Solves the one-dimensional countrflow diffusion flame in cartesian coordinates
function [xint, Sxint, calculatedStrain, maxTemperature] = CounterFlowFlame_FullChemistry2(C, newSolution, xinit, lambda)

%% Solver configuration
cpconstant = true;
chemActive = C.chemActive;
variableKineticParameters = C.variableKineticParameters;

Prandtl = 0.75;


Le = [1.0,1.0,1.0,1.0,1.0];
% Prandtl = 1.0;%

%% Problem description
L = C.L; % Domain length, m
cp = C.cp; % Heat capacity;
s = C.s;
p0 = C.p0; % Pressure, Pa
R = C.R; % gas constant, kg m^2 / (s^2 K mol)
MM = C.MM;
v0 = C.vLeft;
vl = C.vRight;
Coef_Stoic = C.Coef_Stoic;
Tref = C.Tref;
B = 6.9e11; %  Pre-exponential factor 6.9e14 cm3/(mol s) =>  6.9e11m3/(kmol s)
Ta0 = 15900; %K
massHeatRelease = C.Q * MM(1);
viscosity0 = C.viscosity0;


%% Boundary Condition values

TF0 = C.TF0; % Temperature left (fuel)
TO0 = C.TO0; % Temperature right (oxidizer)
YLeft = C.fuelInletConcentration;
YRight = C.oxidizerInletConcentration;
YF0 = YLeft(1);
YO0 = YRight(2);

if (chemActive)
    sol = bvpinit(xinit, @mat4init, lambda);
    sol.y = newSolution;
else
    sol = bvpinit(linspace(0, L, 100), @mat4init, 10);
end

%% =====================================================================
% Solve system

options = bvpset('stats', 'true', 'Vectorized', 'off', 'NMax', 10000, 'AbsTol', 1e-6, 'RelTol', 1e-3);
sol = bvp4c(@mat4ode, @mat4bc, sol, options);
lambda = sol.parameters;

fprintf('Turning on variable Cp\n');
cpconstant = false;
maxCpValue = 0;
max_x = 0;
Le = [0.97, 1.11, 1.39, 0.83, 1.0];
sol = bvp4c(@mat4ode, @mat4bc, sol, options);
Sxint = deval(sol, linspace(0, L, 1000));
calculatedStrain = max(abs(Sxint(2, :)));
fprintf('Strain (biggest du/dy magnitude): %7.3f.\n', calculatedStrain);
maxTemperature = max(Sxint(4, :));
fprintf('Maximum temperature: %7.3f.\n', maxTemperature);

fprintf('Maximum cp: %7.3f.\n', maxCpValue);
fprintf('Position of maximum cp: %7.6f.\n', max_x);

xint = linspace(0, L, 200);
Sxint = deval(sol, xint);

%% -----------------------------------------------------------------------
% Note: the system to be solved is of first order on v, second order on U,
% and second order on the scalars (T, Y1,Y2...)
% Second order equations are brought to a first order by introducing a
% transformation
    function dydx = mat4ode(x, y, lambda) % equation being solved
        v = y(1);
        U1 = y(2);
        U2 = y(3);
        T1 = y(4); % T
        T2 = y(5); % dT/dy
        Y1_1 = y(6); % Y1
        Y1_2 = y(7); % dY1/dy
        Y2_1 = y(8);
        Y2_2 = y(9);
        Y3_1 = y(10);
        Y3_2 = y(11);
        Y4_1 = y(12);
        Y4_2 = y(13);
        Y5_1 = y(14);
        Y5_2 = y(15);

        Yk_1 = [Y1_1; Y2_1; Y3_1; Y4_1; Y5_1];
        Yk_2 = [Y1_2, Y2_2, Y3_2, Y4_2, Y5_2]; % Array of derivatives
        % Try to "repair" temperature and mass fractions fields by limiting their value
        T1(T1 < 300) = 300;
        Yk_1(Yk_1 < 0) = 0;
        Yk_1(Yk_1 > 1) = 1;

        %% Density and transport parameters
        rho_ = rho(T1, Yk_1);
        rho_p = drhody(T1, T2, Yk_1, Yk_2);
        mu_ = mu(T1);
        mu_p = dmu_dy(T1, T2);
         if cpconstant
            cp = C.cp;
        else
%                        cp = getComponentCp(T1, 'N2');
              cp = getMixtureCp(T1, Yk_1, ["CH4", "O2", "CO2", "H2O", "N2"]);
              if cp > maxCpValue
                  max_x = x;
                  maxCpValue = cp;
              end

        end
        k_cp = mu_ / (Prandtl); % k/cp = mu / pr
        k_cp_p = mu_p / (Prandtl); 

        rhoD1 = mu_ / (Prandtl*Le(1));
        rhoD2 = mu_/ (Prandtl*Le(2));
        rhoD3 = mu_ / (Prandtl*Le(3));
        rhoD4 = mu_ /(Prandtl* Le(4));
        rhoD5 = mu_ / (Prandtl*Le(5));
        rhoD1_p = mu_p / (Prandtl*Le(1));
        rhoD2_p = mu_p / (Prandtl*Le(2));
        rhoD3_p = mu_p /(Prandtl* Le(3));
        rhoD4_p = mu_p / (Prandtl*Le(4));
        rhoD5_p = mu_p /(Prandtl* Le(5));
      
        Q = getHeatRelease(Yk_1);
        omega = getReactionRate(T1, Yk_1);
       
        %
        dydx = [(-rho_ * U1 - rho_p * v) / rho_; ... %Conti
            U2; ... %Mom
            (lambda + rho_ * v * U2 + rho_ * U1^2 - mu_p * U2) * 1.0 / mu_; ... %Mom
            T2; ... .%Energy
            ( rho_ * vss * T2 - k_cp_p * T2 - Q * omega/cp) / (k_cp); ... %Energy
            Y1_2; ... .%MassFraction1
            (rho_ * v * Y1_2 - rhoD1_p * Y1_2 - omega * Coef_Stoic(1) * MM(1)) / (rhoD1); ... %MassFraction1
            Y2_2; ... .%MassFraction2
            (rho_ * v * Y2_2 - rhoD2_p * Y2_2 - omega * Coef_Stoic(2) * MM(2)) / (rhoD2); ... %MassFraction2
            Y3_2; ... .%MassFraction3
            (rho_ * v * Y3_2 - rhoD3_p * Y3_2 - omega * Coef_Stoic(3) * MM(3)) / (rhoD3); ... %MassFraction3
            Y4_2; ... .%MassFraction4
            (rho_ * v * Y4_2 - rhoD4_p * Y4_2 - omega * Coef_Stoic(4) * MM(4)) / (rhoD4); ... %MassFraction4
            Y5_2; ... .%MassFraction5
            (rho_ * v * Y5_2 - rhoD5_p * Y5_2 - omega * Coef_Stoic(5) * MM(5)) / (rhoD5); ... %MassFraction5
            ];

    end
% -----------------------------------------------------------------------

    function res = mat4bc(ya, yb, lambda) % boundary conditions
        va = ya(1);
        U1a = ya(2);
        TO0a = ya(4);

        Y1_1a = ya(6);
        Y2_1a = ya(8);
        Y3_1a = ya(10);
        Y4_1a = ya(12);
        Y5_1a = ya(14);

        vb = yb(1);
        U1b = yb(2);
        TO0b = yb(4);

        Y1_1b = yb(6);
        Y2_1b = yb(8);
        Y3_1b = yb(10);
        Y4_1b = yb(12);
        Y5_1b = yb(14);
        res = [va - v0; ... % v(0) = v0
            vb - vl; ... % v(L) = vl
            U1a - 0; ... % U(0) = 0
            U1b - 0; ... % U(L) = 0
            TO0a - TF0; ... % T(0) = TF0
            TO0b - TO0; ... %T(L) = TO0
            Y1_1a - YLeft(1); ... %
            Y1_1b - YRight(1); ... %
            Y2_1a - YLeft(2); ... %
            Y2_1b - YRight(2); ... %
            Y3_1a - YLeft(3); ... %
            Y3_1b - YRight(3); ... %
            Y4_1a - YLeft(4); ... %
            Y4_1b - YRight(4); ... %
            Y5_1a - YLeft(5); ... %
            Y5_1b - YRight(5); ... %
            ];


    end
%-----------------------------------------------------------------------
    function yinit = mat4init(x) % initial guess function.
        yinit = [0.0; ... %v
            0.0; ... %U1
            0.0; ... %U2
            300; ... % TO0
            0.0; ... % T2
            1.0; ... % Y1_1
            0.0; ... % Y1_2
            0.0; ... % Y2_1
            0.0; ... % Y2_2
            0.0; ... % Y3_1
            0.0; ... % Y3_2
            0.0; ... % Y4_1
            0.0; ... % Y4_2
            0.0; ... % Y5_1
            0.0; ... % Y5_2
            ];
    end

%% Definition of helper functions for density and transport parameters

    function MW = getAverageMolecularWeight(Yk)
        mult = 0;
        for c = 1:5
            mult = mult + Yk(c) / MM(c);
        end
        MW = 1.0 / mult;

    end

    function density = rho(T, Yk)
        MW = getAverageMolecularWeight(Yk);
        density = p0 * MW / (R * T); % kg / m^3

    end

    function dRho_dy = drhody(T, dTdy, Yk, dYkdy)

        dRho_dy = -p0 * dTdy / (R * T^2 * (Yk(1) / MM(1) + Yk(2) / MM(2) + Yk(3) / MM(3) + Yk(4) / MM(4) + Yk(5) / MM(5))) - ...
            p0 * (dYkdy(1) / MM(1) + dYkdy(2) / MM(2) + dYkdy(3) / MM(3) + dYkdy(4) / MM(4) + dYkdy(5) / MM(5)) ...
            / (R * T * (Yk(1) / MM(1) + Yk(2) / MM(2) + Yk(3) / MM(3) + Yk(4) / MM(4) + Yk(5) / MM(5))^2);

    end

%% PowerLaw
%     function viscosity = mu(T)
%         viscosity = viscosity0 * (T / Tref)^(a);
%
%     end
%     function dViscosity_dy = dmu_dy(T, dTdy)
%              dViscosity_dy = (viscosity0 / Tref^a) * a * (T)^(a - 1) * dTdy;
%     end
    function viscosity = mu(T)
        S = 110.4;
        viscosity = viscosity0 * (T / Tref)^(3 / 2) * (Tref + S) / (T + S);

    end
    function dViscosity_dy = dmu_dy(T, dTdy)
        mu0 = viscosity0;
        S = 110.4;
        dViscosity_dT = 0.3e1 / 0.2e1 * mu0 * sqrt(T/Tref) * (Tref + S) / (T + S) / Tref - ...
            mu0 * (T / Tref)^(0.3e1 / 0.2e1) * (Tref + S) / (T + S)^2;
        dViscosity_dy = dViscosity_dT * dTdy;

    end


%-------------------------------------------

    function reactionRate = getReactionRate(T, Yk)
        Yf = Yk(1);
        Yo = Yk(2);
        rho_ = rho(T, Yk);
        if chemActive == true
            Ta = getActivationTemperature(T, Yk);
            reactionRate = B * exp(-Ta./T) .* (rho_ .* Yf ./ MM(1)) .* (rho_ .* Yo ./ MM(2));
        else
            reactionRate = 0.0;
        end
    end

    function heatRelease = getHeatRelease(Yk)
        if variableKineticParameters
            phi = GetPhi(Yk(1), Yk(2));
            if phi > 1.5
                phi = 1.5;
            end

            alpha = 0.21;
            if phi > 1
                heatRelease = (1.0 - alpha * (phi - 1)) * massHeatRelease;
            else
                heatRelease = massHeatRelease;
            end
        else
            heatRelease = massHeatRelease;
        end
    end

    function activationTemperature = getActivationTemperature(T, Yk)
        if variableKineticParameters
            phi = GetPhi(Yk(1), Yk(2));
            if (phi >= 1.07)
                activationTemperature = (1 + 4.443 * (phi - 1.07)^2) * Ta0;
            elseif (phi <= 1.07 && phi > 0.64)
                activationTemperature = Ta0;
            else
                activationTemperature = (1 + 8.25 * (phi - 0.64)^2) * Ta0;
            end
        else
            activationTemperature = Ta0;
        end

    end

    function phi = GetPhi(YF, YO)
        phi = (s * YF0 / YO0) * (s * YF - YO + YO0) / (s * (YF0 - YF) + YO);
    end


end
