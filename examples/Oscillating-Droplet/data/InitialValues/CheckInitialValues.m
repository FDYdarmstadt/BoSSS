clear all

%% m3

P3=@(A,eta, theta) A(theta) + eta.*(1/2).*(5.*((cos(theta)).^3)-(3.*(cos(theta))));

A1 =@(theta) 0.996786 - 0.000315584.*cos(theta); eta1 = 0.15; 
A2 =@(theta) 0.987143 - 0.00252468.*cos(theta); eta2 = 0.3; 
A3 =@(theta) 0.977143 - 0.00598442.*cos(theta); eta3 = 0.4; 

%% input data

dat1 = importdata('m3/surfaceDrop_m3_Oh01_eta015.txt');
dat2 = importdata('m3/surfaceDrop_m3_Oh01_eta03.txt');
dat3 = importdata('m3/surfaceDrop_m3_Oh01_eta04.txt');    

%% check data

errNorm1 = 0;
expr1 = zeros(length(dat1),2);
for i=1:length(dat1)
    expr1(i,1) = dat1(i,1);
    radiusDat = dat1(i,2);
    radiusExpr = P3(A1,eta1,dat1(i,1));
    expr1(i,2) = radiusExpr;
    errNorm1 = errNorm1 + (radiusExpr - radiusDat)^2;
end
sqrt(errNorm1)

errNorm2 = 0;
expr2 = zeros(length(dat2),2);
for i=1:length(dat2)
    expr2(i,1) = dat2(i,1);
    radiusDat = dat2(i,2);
    radiusExpr = P3(A2,eta2,dat2(i,1));
    expr2(i,2) = radiusExpr;
    errNorm2 = errNorm2 + (radiusExpr - radiusDat)^2;
end
sqrt(errNorm2)

errNorm3 = 0;
expr3 = zeros(length(dat3),2);
for i=1:length(dat3)
    expr3(i,1) = dat3(i,1);
    radiusDat = dat3(i,2);
    radiusExpr = P3(A3,eta3,dat3(i,1));
    expr3(i,2) = radiusExpr;
    errNorm3 = errNorm3 + (radiusExpr - radiusDat)^2;
end
sqrt(errNorm3)


%% plot data

figure
plot(dat1(:,1),dat1(:,2),expr1(:,1),expr1(:,2))
xlabel('theta')
ylabel('radius')
title('m=3 eta=0.15')

figure
plot(dat2(:,1),dat2(:,2),expr2(:,1),expr2(:,2))
xlabel('theta')
ylabel('radius')
title('m=3 eta=0.3')

figure
plot(dat3(:,1),dat3(:,2),expr3(:,1),expr3(:,2))
xlabel('theta')
ylabel('radius')
title('m=3 eta=0.4')



%% error plot

% figure
% plot(dat1(:,1),dat1(:,2) - expr1(:,2), '-') %, dat1(:,1), 0.5*(dat1(end,2) - expr1(end,2)) + 0.5*(dat1(end,2) - expr1(end,2)).*cos(dat1(:,1)+pi), '-.')
% xlabel('theta')
% ylabel('error')
% title('m=3 eta=0.15')
% 
% figure
% plot(dat2(:,1),dat2(:,2) - expr2(:,2), '-') %, dat2(:,1), 0.5*(dat2(end,2) - expr2(end,2)) + 0.5*(dat2(end,2) - expr2(end,2)).*cos(dat2(:,1)+pi), '-.')
% xlabel('theta')
% ylabel('error')
% title('m=3 eta=0.3')
% 
% figure
% plot(dat3(:,1),dat3(:,2) - expr3(:,2), '-') %, dat3(:,1), 0.5*(dat3(end,2) - expr3(end,2)) + 0.5*(dat3(end,2) - expr3(end,2)).*cos(dat3(:,1)+pi), '-.')
% xlabel('theta')
% ylabel('error')
% title('m=3 eta=0.4')


%% check center of mass position

P3_vol1 =@(theta) 0.5*P3(A1, eta1, theta).^2; 
P3_vol2 =@(theta) 0.5*P3(A2, eta2, theta).^2; 
P3_vol3 =@(theta) 0.5*P3(A3, eta3, theta).^2; 

volume1 = integral(P3_vol1, 0, pi);
volume2 = integral(P3_vol2, 0, pi);
volume3 = integral(P3_vol3, 0, pi);

P3_cy1 =@(theta) (1/3).*cos(theta).*(P3(A1, eta1, theta).^3);

centery = integral(P3_cy1, 0, pi)/volume1


