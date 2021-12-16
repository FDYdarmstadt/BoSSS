function plotMixtureFractionSolution(mySolution,myconfig)
%% Find point where z = zst
x = mySolution(1, :);
z = mySolution(11, :);
zst = myconfig.zst;
[minValue, indexLow] = min(abs(z - zst));
stoicPoint = x(indexLow);



subplot(2, 2, 4)
plot(mySolution(1, :), mySolution(2, :))
% xline(stoicPoint);
axis([0, 0.02, -3, 3])
ylim auto;
title('Velocity')
xlabel('x')
ylabel('y')
legend('y')

subplot(2, 2, 3)
plot(mySolution(1, :), mySolution(5, :))
% xline(stoicPoint);
axis([0, 0.02, 0, 1])
title('Temperature')
ylim auto;
xlabel('x')
ylabel('y')
legend('T')


subplot(2, 2, 1)
plot(mySolution(1, :), mySolution(6, :));
hold on
% xline(stoicPoint);
plot(mySolution(1, :), mySolution(7, :));
plot(mySolution(1, :), mySolution(8, :));
plot(mySolution(1, :), mySolution(9, :));
plot(mySolution(1, :), mySolution(10, :));
hold off
axis([0, 0.02, 0, 1])
ylim auto;
title('Mass fractions')
xlabel('x')
ylabel('y')
legend('YCH4', 'YO2', 'CO2', 'YH2O', 'YN2')

subplot(2, 2, 2)
plot(mySolution(1, :), mySolution(11, :))
% xline(stoicPoint);
axis([0, 0.02, 0, 1])
title('MixtureFraction')
ylim auto;
xlabel('x')
ylabel('y')
legend('Z')

drawnow
end