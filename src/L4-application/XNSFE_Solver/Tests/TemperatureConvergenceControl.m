clc
clear all
close all

%%
[R,Phi] = meshgrid(linspace(0,0.8,100), linspace(0,2*pi,1000));

figure;
fplot(@(x) boundary(x), [0 2*pi]);

%%
X = R.*cos(Phi);
Y = R.*sin(Phi);

figure;
surf(X,Y,T(R,Phi))

%%
% extract points along y == 0
r = linspace(-0.8, 0.8, 1000);
f = T(r,0);
f = f/max(f);

x = linspace(-0.8*cos(10/90 * pi/2), 0.8*cos(10/90 * pi/2), 1000);
y = 0.8*sin(10/90 * pi/2);
rr = sqrt(x.^2+y.^2);
pp = atan2(y,x);
g = T(rr,pp);
g = g/max(g);

figure;
plot(r,f); hold on
plot(x,g);
plot(r,sin(2*pi*(r+0.8)/(1.6)));

% save to csv
dat = [r', f'];
writematrix(dat, 'boundarycond90.txt', 'Delimiter', 'tab');
dat = [x', g'];
writematrix(dat, 'boundarycond80.txt', 'Delimiter', 'tab');

%%
function f = T(r,phi)
    R = 0.8*ones(size(r));
    f = 1/(2*pi) .* integral(@(psi) boundary(psi) .* (R.^2 - r.^2)./(r.^2-2*R.*r.*cos(phi-psi)+R.^2),0,2*pi, 'ArrayValued', true);
end

function f = boundary(psi)
%     f = 1.*ones(size(psi));
%     f = kernel(psi, 1, 1.25*pi, 0.25 * pi) + kernel(psi, -1, 1.75*pi, 0.25 * pi);
    f = smootherstepkernel(psi, 1, 1.25*pi, 0.1 * pi) + smootherstepkernel(psi, -1, 1.75*pi, 0.1 * pi);
end

function f = kernel(x, scale, xm, width)
    f = ((x >= xm-width) & (x <= xm+width)) .* scale .* exp(-1) .* exp(-1/(1-((x-xm)/(width)).^2));
end

function f = smoothstepkernel(x, scale, xm, width)
    y = (x - (xm - width)) / (width);
    f = ((x > xm-width) & (x <= xm)) .* scale .* (3*y.^2-2*y.^3);
    y = (x - (xm)) / (width);
    f = f + ...
        ((x > xm) & (x <= xm+width)) .* scale .* (1-(3*y.^2-2*y.^3));
end

function f = smootherstepkernel(x, scale, xm, width)
    y = (x - (xm - width)) / (width);
    f = ((x > xm-width) & (x <= xm)) .* scale .* (6*y.^5-15*y.^4+10*y.^3);
    y = (x - (xm)) / (width);
    f = f + ...
        ((x > xm) & (x <= xm+width)) .* scale .* (1-(6*y.^5-15*y.^4+10*y.^3));
end