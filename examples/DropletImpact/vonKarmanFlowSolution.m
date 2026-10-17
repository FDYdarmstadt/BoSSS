clear all


%% nonlinear solver
roots = fsolve(@objective, [0.51 -0.62], optimset("Display","iter")); 


%%
opts = odeset('MaxStep', 0.005);
[z, g] = ode45(@(z,g) odefun(z,g), [0 20], [0 0 0 roots(1) roots(2) 0], opts);


%% plot solutions
size_legend = 16;
size_label = 18;
size_tick = 14;
size_marker = 8;


figure
subplot(2,2,1)
plot(g(:,1), z, 'Linewidth', 1.2)
xlabel('$U$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 

subplot(2,2,2)
plot(g(:,2), z, 'Linewidth', 1.2)
xlabel('$V$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 

subplot(2,2,3)
plot(g(:,3), z, 'Linewidth', 1.2)
xlabel('$W$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 

subplot(2,2,4)
plot(g(:,6), z, 'Linewidth', 1.2)
xlabel('$P$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 


set(gcf, 'color', 'w', 'Position', [100, 000, 1200, 800])


%% save
saveData = [z, g(:,1)];
save('vonKarmanFlowSolutionMatlab_VelocityU.txt', 'saveData', '-ascii')
saveData = [z, g(:,2)];
save('vonKarmanFlowSolutionMatlab_VelocityV.txt', 'saveData', '-ascii')
saveData = [z, g(:,3)];
save('vonKarmanFlowSolutionMatlab_VelocityW.txt', 'saveData', '-ascii')
saveData = [z, g(:,6)];
save('vonKarmanFlowSolutionMatlab_PressureP.txt', 'saveData', '-ascii')


%% curve fitting

% fitU = fit(z, g(:,1), 'rat55')
fitU = fit(z, g(:,1), 'cubicspline')
Uwall = fitU(0);
Uouter = fitU(20);

% fitV = fit(z, g(:,2), 'rat55')
fitV = fit(z, g(:,2), 'cubicspline')
Vwall = fitV(0);
Vouter = fitV(20);

% fitW = fit(z, g(:,3), 'rat55')
fitW = fit(z, g(:,3), 'cubicspline')
Wwall = fitW(0);
Wouter = fitW(20);

% fitP = fit(z, g(:,6), 'rat55')
fitP = fit(z, g(:,6), 'cubicspline')
Pwall = fitP(0);
Pouter = fitP(20);


%% compare fitting

figure
subplot(2,2,1)
plot(g(:,1), z, 'b', fitU(z), z, 'r--', 'Linewidth', 1.2)
xlabel('$U$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 

subplot(2,2,2)
plot(g(:,2), z, 'b', fitV(z), z, 'r--',  'Linewidth', 1.2)
xlabel('$V$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 

subplot(2,2,3)
plot(g(:,3), z, 'b', fitW(z), z, 'r--',  'Linewidth', 1.2)
xlabel('$W$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 

subplot(2,2,4)
plot(g(:,6), z, 'b', fitP(z), z, 'r--',  'Linewidth', 1.2)
xlabel('$P$', 'Interpreter', 'latex', 'FontSize', size_label)
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', size_label)
set(gca, 'FontSize', size_tick) 


set(gcf, 'color', 'w', 'Position', [100, 000, 1200, 800])


%% function definitions

% similarity problem (shooting method)
function dgdz = odefun(z,g)
dgdz = zeros(6,1);
dgdz(1) = g(4);
dgdz(2) = g(5);
dgdz(3) = -2.0*g(1);
dgdz(4) = g(1).*g(1) - (g(2)+1.0).^2 + g(4).*g(3);
dgdz(5) = 2.0.*g(1).*(g(2)+1.0) + g(5).*g(3);
dgdz(6) = 2.0*g(3).*g(1) - 2.0*g(4);
end

% objective function for the shooting variables
function residual = objective(x)
opts = odeset('MaxStep', 0.005);
[z, g] = ode45(@(z,g) odefun(z,g), [0 20], [0 0 0 x(1) x(2) 0], opts);
target1 = 0.0;
target2 = -1.0;
residual = [g(end,1)-target1 g(end,2)-target2];
end

