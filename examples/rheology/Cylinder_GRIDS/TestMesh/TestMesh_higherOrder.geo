Geometry.CopyMeshingMethod = 1 ;
Mesh.CharacteristicLengthMax = 0.25 ;

NX1  =  3 ; RX1 = 1.00 ;
NY1  =  3 ; RY1  = 1.00 ;


Z0 = 0.0 ; // plane
SB =  0.1 ; // dummy size


X0 = -1;
X1 = 1 ;

Y0 = -1 ;
Y1 = 1 ;

Point( 0) = { X0, Y0, Z0, SB};
Point( 1) = { X1, Y0, Z0, SB};
Point( 2) = { X0, Y1, Z0, SB};
Point( 3) = { X1+0.0001, Y1, Z0, SB};


Line(1) = {0,1}; Transfinite Line{1} = NX1 Using Progression RX1 ;
Line(2) = {0,2}; Transfinite Line{2} = NY1 Using Progression RY1 ;
Line(3) = {2,3}; Transfinite Line{3} = NX1 Using Progression RX1 ;
Line(4) = {1,3}; Transfinite Line{4} = NY1 Using Progression RY1 ;

Line Loop(1) = {1,4,-3,-2};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {0,1,2,3};

Recombine Surface(1);