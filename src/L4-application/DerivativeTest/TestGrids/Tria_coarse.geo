//standard size of element
cl__1 = 1.4; 
//Radius of cylinder
R=0.5;
R2=1.5;

// Meshing parameters
NE_C = 80; //elements around Cylinder
NE_R = 64; //elements Region around cylinder


// Box size
Factor = 20;
X0=Factor*R;
Y0=Factor*R;
X1=Factor*R+Factor*R; //behind cylinder

// Spacing parameters
SP1 = 1; //middle of cylinder
SP2 = 0.9;  // in radial direction

// Box around cylinder
Point(1) = {-X0, -Y0, 0, cl__1};
Point(2) = {-X0, Y0, 0, cl__1};
Point(3) = {X1, Y0, 0, cl__1};
Point(4) = {X1, -Y0, 0, cl__1};
// middle
Point(5) = {0, 0, 0, 1};

// Points around cylinder
Point(6) = {-0.707106781*R, -0.707106781*R, 0, 1};

// Points around cylinder
Point(7) = {-0.707106781*R2, -0.707106781*R2, 0, 1};


// Box
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cylinder
Circle(5) = {6, 5, 6};
Transfinite Line {5} = NE_C Using Bump SP1;
// Cylinder
Circle(6) = {7, 5, 7};
Transfinite Line {6} = NE_R Using Bump SP1;

Line Loop(11) = {1,2,3,4};
Line Loop(12) = {5};
Line Loop(13) = {6};

Plane Surface(1) = {11,13};
Plane Surface(2) = {13,12};

//Recombine Surface(1);
//Recombine Surface(2);
//Recombine Surface(3);
//Recombine Surface(4);
//Recombine Surface(5);
