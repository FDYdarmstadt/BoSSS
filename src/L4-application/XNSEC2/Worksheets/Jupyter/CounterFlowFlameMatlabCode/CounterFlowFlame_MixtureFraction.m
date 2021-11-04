function [solutionArray, lambdaOut] = CounterFlowFlame_MixtureFraction(myconfig)
options = bvpset('stats', 'off', 'NMax', 10000, 'AbsTol', 1e-8);

%% Solver configuration
lambda = -100; %Initial estimation for dpdz
TF0 = myconfig.TF0;
TO0 = myconfig.TO0;
Q = myconfig.Q;
cp = myconfig.cp;
zst = myconfig.zst;
L = myconfig.L;
YF0 = myconfig.fuelInletConcentration(1);
YO0 = myconfig.oxidizerInletConcentration(2);
Coef_Stoic = myconfig.Coef_Stoic;
MM = myconfig.MM;
Tref = myconfig.Tref;
a = myconfig.a;
Prandtl =1;

%%
% 
% cp_f = 1;
% t_in = 1800;
% res = 1e10;
% Y1_1 = getY1FromZ(zst);
% Y2_1 = getY2FromZ(zst);
% Y3_1 = getY3FromZ(zst);
% Y4_1 = getY4FromZ(zst);
% Y5_1 = getY5FromZ(zst);
% 
% Yk_1 = [Y1_1; Y2_1; Y3_1; Y4_1; Y5_1];
%        it = 0;
% while res > 1e-6
%     
%     cpcalc = getMixtureCp(t_in, Yk_1, ["CH4", "O2", "CO2", "H2O", "N2"]);
%     
%    Tcalc = TF0 + Q*zst*YF0/cpcalc;
%    res = abs(Tcalc-t_in);
%     t_in = Tcalc;
%     it = it+1;
% end


%% Solve system

sol = bvpinit(linspace(0, L, myconfig.initialCellNumber), @mat4init, lambda);

fprintf('Solving NS+mixture fraction equations.\n');
sol = bvp4c(@mat4ode, @mat4bc, sol, options);


Sxint = deval(sol, linspace(0, L, 200));
fprintf('Strain (biggest du/dy magnitude): %7.3f.\n', max(Sxint(2, :)));


lambdaOut = sol.parameters;

Sxint = deval(sol, sol.x);

%% Recover primitive variables from Z
ZArray = Sxint(4, :);

for i = 1:length(ZArray)
    RecoveredT(i) = getTemperatureFromZ(ZArray(i));
    RecoveredY1(i) = getY1FromZ(ZArray(i));
    RecoveredY2(i) = getY2FromZ(ZArray(i));
    RecoveredY3(i) = getY3FromZ(ZArray(i));
    RecoveredY4(i) = getY4FromZ(ZArray(i));
    RecoveredY5(i) = getY5FromZ(ZArray(i));
end

solutionArray = [sol.x; Sxint(1, :); Sxint(2, :); Sxint(3, :); RecoveredT; RecoveredY1; RecoveredY2; RecoveredY3; RecoveredY4; RecoveredY5; Sxint(4, :)];

%% -----------------------------------------------------------------------
% Note: the system to be solved is of first order on v, second order on U,
% and second order on the scalars (T, Y1,Y2...)
% Second order equations are brought to a first order by introducing a
% transformation


    function dydx = mat4ode(~, y, lambda) % equation being solved
        v = y(1);
        U1 = y(2);
        U2 = y(3);
        Z1 = y(4); % Z
        Z2 = y(5); % dZ/dy
        if (Z1 > 1)
            Z1 = 1;
        end
        if (Z1 < 0)
            Z1 = 0;
        end

        rho_ = GetRhoFromZ(Z1);
        mu_ = GetMuFromZ(Z1)/Prandtl;

        rho_p = 0.0; %drhody(Z1, Z2); % It seems to be problematic, because the derivative at zSt is not defined. Is simply ignored (seems to still give an adequate starting solution )
        mu_p = 0.0; %dmu_dy(Z1, Z2); % It seems to be problematic, because the derivative at zSt is not defined. Is simply ignored (seems to still give an adequate starting solution )


        dydx = [(-rho_ * U1 - rho_p * v) / rho_; ... %Conti
            U2; ... %Mom
            (lambda + rho_ * v * U2 + rho_ * U1^2 - mu_p * U2) * 1.0 / mu_; ... %Mom
            Z2; ... .%MassFraction
            (rho_ * v * Z2 - mu_p * Z2) / mu_; ... %MassFraction
            ];
    end
% -----------------------------------------------------------------------

    function res = mat4bc(ya, yb, lambda) % boundary conditions
        va = ya(1);
        U1a = ya(2);
        Za = ya(4);

        vb = yb(1);
        U1b = yb(2);
        Zb = yb(4);

        res = [va - myconfig.vLeft; ... % v(0) = v0
            vb - myconfig.vRight; ... % v(L) = vl
            U1a - 0; ... % U(0) = 0
            U1b - 0; ... % U(L) = 0
            Za - 1.0; ... % Z(0) = 1.0
            Zb - 0.0; ... %Z(L) = TL
            ];
    end
%-----------------------------------------------------------------------
    function yinit = mat4init(x) % initial guess function.
        yinit = [0.0; ... %v
            0.0; ... %U1
            0.0; ... %U2
            0.5; ... % Z1
            0.0; ... % Z2
            ];
    end

%% Definition of helper functions for density and transport parameters
    function viscosity = mu(T)
        viscosity = myconfig.viscosity0 * (T / Tref)^(a);
    end



    function density = rho(T, Yk)

        mult = 0.0;
        for c = 1:(length(Yk))
            mult = mult + Yk(c) / MM(c);
        end
        density = myconfig.p0 / (myconfig.R * T * mult); % kg / m^3

    end


%==============================================
%Transformations from Z
%==============================================
    function Temperature_fromZ = getTemperatureFromZ(Z)

        if (Z >= zst)

            Temperature_fromZ = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * zst * (1 - Z) / (1 - zst);
        else
            Temperature_fromZ = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * Z;
        end
    end

    function Y1_fromZ = getY1FromZ(Z)
        if (Z >= zst)
            Y1_fromZ = YF0 * (Z - zst) / (1 - zst);
        else
            Y1_fromZ = 0.0;
        end
    end


    function Y2_fromZ = getY2FromZ(Z)
        if (Z >= zst)
            Y2_fromZ = 0;
        else
            Y2_fromZ = YO0 * (1 - Z / zst);
        end
    end

    function Y3_fromZ = getY3FromZ(Z)
        nu_CO2 = Coef_Stoic(3);
        nu_O2 = Coef_Stoic(2);
        nu_CH4 = Coef_Stoic(1);
        MW_CO2 = MM(3);
        MW_O2 = MM(2);
        MW_CH4 = MM(1);
        if (Z >= zst)
            Y3_fromZ = -YO0 * (nu_CO2 * MW_CO2) / (nu_O2 * MW_O2) * (1 - Z);
        else
            Y3_fromZ = -YF0 * (nu_CO2 * MW_CO2) / (nu_CH4 * MW_CH4) * Z;
        end
    end

    function Y4_fromZ = getY4FromZ(Z)
        nu_H2O = Coef_Stoic(4);
        nu_O2 = Coef_Stoic(2);
        nu_CH4 = Coef_Stoic(1);
        MW_H2O = MM(4);
        MW_O2 = MM(2);
        MW_CH4 = MM(1);

        if (Z >= zst)
            Y4_fromZ = -YO0 * (nu_H2O * MW_H2O) / (nu_O2 * MW_O2) * (1 - Z);
        else
            Y4_fromZ = -YF0 * (nu_H2O * MW_H2O) / (nu_CH4 * MW_CH4) * Z;
        end
    end

    function Y5_fromZ = getY5FromZ(Z)
        if (Z >= zst)
            YNOxi0 = 1.0 - YO0;
            YNFuel0 = 1.0 - YF0;
            Y5_fromZ = YNOxi0 * (1 - Z) + YNFuel0 * Z;
        else
            YNOxi0 = 1.0 - YO0;
            YNFuel0 = 1.0 - YF0;
            Y5_fromZ = YNOxi0 * (1 - Z) + YNFuel0 * Z;
        end
    end
    function MuFromZ = GetMuFromZ(z)
        MuFromZ = mu(getTemperatureFromZ(z));
    end
    function rhoFromZ = GetRhoFromZ(Z)
        T = getTemperatureFromZ(Z);
        Y1 = getY1FromZ(Z);
        Y2 = getY2FromZ(Z);
        Y3 = getY3FromZ(Z);
        Y4 = getY4FromZ(Z);
        Y5 = getY5FromZ(Z);

        rhoFromZ = rho(T, [Y1, Y2, Y3, Y4, Y5]);
    end


end
