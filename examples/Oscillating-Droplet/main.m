                                                           clear

polVel = dlmread('polarVelCase1.txt');
radVel = dlmread('radialVelCase1.txt');

rS = polVel(:,1);
tS = polVel(:,2);
I = length(tS);

x = zeros(I,1);
z = zeros(I,1);
Vx = zeros(I,1);
Vz = zeros(I,1);


for i = 1:I
    theta = tS(i);
    rdius = rS(i);
    
    %vr = -abs(polVel(i,3));
    vr = polVel(i,3);
    vt = radVel(i,3);
        
    x(i) = sin(theta)*rdius;
    z(i) = cos(theta)*rdius;
    theta_rec(i) = atan2(x(i),z(i));
    
    %theta - theta_rec(i)
    
    Vx(i) = cos(theta)*vr;
    Vx(i) = Vx(i) + sin(theta)*vt;
    Vz(i) = sin(theta)*vr;
    Vz(i) = Vz(i) - cos(theta)*vt;
end


quiver(x, z, Vx, Vz);
axis equal

%plot(1:I,tS,'+',1:I,theta_rec,'o')

