clear

polVel = dlmread('InitialValues/polarVelCase1.txt');
radVel = dlmread('InitialValues/radialVelCase1.txt');

rS = polVel(:,1);
tS = polVel(:,2);
I = length(tS);

x = zeros(I,1);
z = zeros(I,1);
Vx = zeros(I,1);
Vz = zeros(I,1);
mV = zeros(I,1);

thetaErr = 0;
for i = 1:I
    theta = tS(i);
    radius = rS(i);
    
    %vr = -abs(polVel(i,3));
    vr = radVel(i,3);
    vt = polVel(i,3);
        
    x(i) = sin(theta)*radius;
    z(i) = cos(theta)*radius;
    theta_rec(i) = atan2(x(i),z(i));
    
    thetaErr = thetaErr + abs(theta - theta_rec(i));
    
    % what I think:
    %Vx(i) = sin(theta)*vr + cos(theta)*vt;
    %Vz(i) = -cos(theta)*vr + sin(theta)*vt ;
    
    % dino zrnic
%     Vx(i) = vr*sin(theta) + vt*cos(theta);
%     Vz(i) = vr*cos(theta) - vt*sin(theta); 

    Vx(i) = vr*sin(theta) + vt*cos(theta);
    Vz(i) = vr*cos(theta) - vt*sin(theta); 

    mV(i) = sqrt(vr^2 + vt^2);
end

figure
quiver(x, z, Vx, Vz);
axis equal

max(abs(Vx))
max(abs(Vz))

%plot(1:I,tS,'+',1:I,theta_rec,'o')

figure
scatter(x, z, mV)
axis equal