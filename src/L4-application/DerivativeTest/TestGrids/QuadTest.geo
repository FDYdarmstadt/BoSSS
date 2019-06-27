//gmsh geometry (.geo) for testing curved quad elements using polynomials up to order 5
// higher orders than 5 (>=6) are not implemented!


//MESHING METHODS
Geometry.CopyMeshingMethod = 1 ;
Mesh.CharacteristicLengthMax = 0.25 ;

//NUMBER OF NODES
NX1  =  30 ; RX1 = 1.00 ;
NY1  =  8 ; RY1  = 1.00 ;

//MEASUREMENTS
Z0 =  0.0 ; // plane
SB =  0.1 ; // dummy size
R  =  1.0 ; // unity radius
F  =  0.5 ; // factor around the cylinder

//COORDINATES
X0 = -15.0 ;
X1 = -(1+F)* R ;
Y0 = 0.0 ;
Y2 = 0.70710678 * R * (1+F) ;

//POINTS
Point( 0) = { 0, 0, Z0, SB};
Point( 1) = { X0, Y0, Z0, SB};
Point( 2) = { X1, Y0, Z0, SB};
Point( 3) = { X0, Y2, Z0, SB};
Point(4) = {-Y2, Y2, Z0, SB};

//LINES
Line(1) = {1,2}; Transfinite Line{1} = NX1 Using Progression RX1 ;
Line( 2) = {3,4}; Transfinite Line{2} = NX1 Using Progression RX1 ;
Line(3) = {1,3}; Transfinite Line{3} = NY1 Using Progression RY1 ;
Circle(4)= {2,0,4}; Transfinite Line{4} = NY1 Using Progression RY1 ;


//LOOPS AND SURFACES
Line Loop(1) = {1,4,-2,-3};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,4,3};
Recombine Surface(1);


