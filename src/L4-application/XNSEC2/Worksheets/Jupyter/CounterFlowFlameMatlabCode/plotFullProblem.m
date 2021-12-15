function plotFullProblem(xint,Sxint)
figure

subplot(2, 3, 1)
plot(xint, Sxint(1, :))
xlim([0,xint(end)]);
ylim auto
title('Velocity')
xlabel('x')
ylabel('y')
legend('u')

subplot(2, 3, 2)
plot(xint, Sxint(4, :))
xlim([0,xint(end)]); 
ylim auto
title('Temperature')
xlabel('x')
ylabel('y')
legend('T')

subplot(2, 3, 3)
plot(xint, Sxint([6, 8, 10, 12, 14], :))
xlim([0,xint(end)]);
ylim auto
title('MassFractions')
xlabel('x')
ylabel('y')
legend('CH4', 'O2', 'CO2', 'CH4', 'N2')

% subplot(2, 3, 4)
% plot(xint, Sxint(16, :))
% xlim([0,xint(end)]);
% ylim auto
% title('Density')
% xlabel('x')
% ylabel('y')
% legend('rho')
% 
% subplot(2, 3, 5)
% plot(xint, Sxint(17, :))
% xlim([0,xint(end)]);
% ylim auto
% title('Viscosity')
% xlabel('x')
% ylabel('y')
% legend('mu')


end