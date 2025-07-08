clc
clear all

% Material values for water at 101.4 kPa
rho_v = 0.597;
rho_l = 958.4;
mu_v = 1.26e-5;
mu_l = 2.8e-4;
c_v = 2030;
c_l = 4216;
k_v = 0.025;
k_l = 0.679;
L = 2.26e6;
sigma = 0.059;
T_sat = 373.15;

% Boundary Condition
T_w = T_sat + 10;

% Calculation of growth factor
xi_fnc = @(xi) xi * exp(xi^2) * erf(xi) - (c_v * (T_w - T_sat))/(sqrt(pi)*L);
xi = fzero(xi_fnc, [0 20]);

% interface position
X = @(t) 2 * xi * sqrt(k_v/(c_v*rho_v) * t);
fplot(X, [0 0.1]);

% Temperature profile for vapor phase
T = @(x,t) T_w + (T_sat - T_w)/erf(xi) * erf(x/(2 * sqrt(k_v/(c_v*rho_v) * t)));
t = 8e-5;
fplot(@(x) T(x,t), [0 X(t)]);