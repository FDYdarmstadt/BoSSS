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
global T_sat
T_sat = 373.15;

% Boundary Condition
global T_w
T_w = T_sat + 5;

% some constants
alpha = k_l/(c_l*rho_l);
beta = rho_v/rho_l;
C = k_l/(rho_v*L);
global B
B = alpha / (C*beta);

% try as ivp, iterate intial gradient until T at infinity is T_w
lambda = 4.013375868517260;
T_end = T_sat;
tspan = 0:10:20;
while abs(T_end - T_w)>1e-8
    solode = ode45(@(t,y) bvpfcn(t,y,lambda),tspan, [T_sat; lambda]);
    T_end = solode.y(1,end);
    disp(T_end)
    lambda = lambda - 0.1 * (T_end - T_w);

end
x = linspace(0,10,1000);
y = deval(solode,x);
plot(x, y(1,:), '-o');

TM = [x' y'];
writematrix(TM, 'SuckingProblemRef.csv');

vs = @(t) C * sqrt(1/(2*alpha*t)) * lambda;
xs = @(t) C * lambda * 1/alpha * sqrt(2*alpha*t);

t0 = (0.001/4 * 0.1 * alpha/(C*lambda))^2/(2*alpha);
xi = linspace(0,0.001 - xs(t0),1000);
eta = @(xi,t) 1/sqrt(2*alpha*t) * xi;
% plot(xi, deval(solode,eta(xi,t0),1), '-o');
% plot(xi, deval(solode,eta(xi,0.01),1), '-o');

% fplot(xs,[0 0.01]);
eta = @(xi,t) 1/sqrt(2*alpha*t) * xi;

function dydx = bvpfcn(eta,y,lambda) % equation to solve
global B
dydx = zeros(2,1);
dydx = [y(2)
       -(eta + 1/B * lambda)*y(2)];
end
%--------------------------------
function res = bcfcn(ya,yb,lambda) % boundary conditions
global T_sat
global T_w
res = [ya(1)- T_sat
       yb(1)- T_w
       ya(2)-lambda];
end
%--------------------------------
function g = guess(eta) % initial guess for y and y'
global T_w
global T_sat
g = [(T_w + (T_sat - T_w) * exp(-eta))
     -(T_sat - T_w) * exp(-eta)];
end
%--------------------------------