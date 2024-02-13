clc
clear all

files = dir('*.csv');

for i=1:length(files)
    data = readmatrix(files(i).name); % or csvread(files(i).name)
    M{i}=data;
end

%% sample for evenly spaced fft
A = M{1};
n = length(A(:,1));
A(:,4) = linspace(A(1,1),A(end,1),n);
% A(:,3) = sin(2*pi*5*A(:,4));
A(:,3) = 0.1*exp(-A(:,4).^2/0.1^2);

y = fft(A(:,3));
fs = n/(-A(1,1)+A(end,1));
f = (0:n-1)*(fs/n); 
power = abs(y).^2/n;
% y0 = fftshift(y);         % shift y values
% f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
% power0 = abs(y0).^2/n;  
plot(f,power)
xlabel('Frequency')
ylabel('Power')

%% sample for non uniform fft
A = M{1};
N = length(A(:,1));
ts = (A(end,1)-A(1,1));
fs = N/ts;
A(:,3) = sin(2*pi*5*A(:,1));
n = A(:,1)*fs;

y = nufft(A(:,3), n,(0:N-1)/N);
plot((0:N-1)/ts,abs(y)/N*2)
xlabel('Frequency f')
ylabel('abs(Y)(f)')

%% example comparison
ts=10; 
N=426;  
fs=N/ts;
%fft for uniform data
t=(0:N-1)/fs;
x2= 2*sin(2*5*pi*t) + sin(2*2*pi*t);
Y2 = fft(x2);
figure
plot((0:N-1)/ts,abs(Y2)/N*2);
%nufft
n = A(:,1)*fs;
x = 2*sin(2*5*pi*n/fs) + sin(2*2*pi*n/fs);
Y = nufft(x,n,(0:N-1)/N);
hold on;plot((0:N-1)/ts,abs(Y)/N*2)

%%
F = struct([]);

for i=1:length(M)
    A = M{i};
    N = length(A(:,1));
    ts = (A(end,1)-A(1,1));
    fs = N/ts;
    n = A(:,1)*fs;

    y = nufft(A(:,2), n,(0:N-1)/N);
    plot((0:N-1)/ts,abs(y)/N*2);
    xlabel('Frequency f')
    ylabel('abs(Y)(f)')
    axis manual 
    axis([0 40 0 0.005]);
    F = [F, getframe(gcf)];
end

%%
fig = figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 850, 800]);
movie(fig,F,1, 10)
