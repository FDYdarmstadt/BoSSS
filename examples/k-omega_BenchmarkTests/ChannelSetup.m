clear all

format long

%% 
Re_tau = 395;
delta = 1; % channel half width
%viscosity = 1e-4; % kinematic viscosity 
u_tau = 1;
density = 1;
beta = 3/40;

viscosity  = u_tau * delta / Re_tau % kinematic viscosity 
tau_w = u_tau^2 * density;
dpwdx = - u_tau^2 * (density / delta);
dudy = tau_w * (density / viscosity);

%% corresponding laminar velocity profile (for initial state)

U_0 = (tau_w * delta) / (2 * density * viscosity);


%% boundary condition for omega tilde

yPlus = 0.3; % in a range of [0.5, 1, 2, 3];
alpha_p = [0.37, 8.21e-2, 3.57e-2, 1.99e-2, 1.27e-2]; % for orders p = {0, 1, 2, 3, 4}
omegaT_w = log(6*tau_w/(beta*density*viscosity)) + (2*log(1./alpha_p)) + (2*log(1/yPlus));