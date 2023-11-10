// Gmsh project created on Fri Dec 09 09:51:06 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {0.5, 0, 0, 1.0};
//+
Point(4) = {1.5, 1, 0, 1.0};
//+
Point(5) = {1.5, 0.17632698070846498, 0, 1.0};
//+
Point(6) = {1.5, 0.8188966492599938, 0, 1.0};
//+
Line(1) = {2, 4};
//+
Line(2) = {4, 6};
//+
Line(3) = {6, 5};
//+
Line(4) = {5, 3};
//+
Line(5) = {3, 1};
//+
Line(6) = {1, 2};
//+
Line(7) = {3, 6};
//+
Curve Loop(1) = {7, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 1, 2, -7};
//+
Plane Surface(2) = {2};
