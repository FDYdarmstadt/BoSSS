clc
clear variables
close all


%% Read data
path = 'C:\BoSSS\src\private-die\L4-application\StokesHelical_Ak\bin\Debug\';

% AV on in same cells
file = 'matrixBoSSS.txt';
matrix = ReadMsr2([path file]);


%% Calculate eigenvalues
[eigenvectors, eigenvalues] = eig(-full(matrix));
eigenvalues = diag(eigenvalues);
% fullMatrix = full(matrix);


%% Plot stability region
figure;
plot1 = plot(real(eigenvalues), imag(eigenvalues), 'bx');


%% Plot format
legend([plot1], {'Eigenvalues'});
title(['Title']);
xlabel('Real(z)');
ylabel('Imag(z)');

