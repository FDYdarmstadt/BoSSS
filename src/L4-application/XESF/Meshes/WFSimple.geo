// Gmsh project created on Sat Dec 03 16:53:15 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {-0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1.5, -0, 0, 1.0};
//+
Recursive Delete {
  Point{3}; 
}
//+
Point(3) = {0.5, 0, 0, 1.0};
//+
Point(4) = {1.5, 1, 0, 1.0};
//+
Point(5) = {1.5, 0.17632698070846498, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 5};
//+
Line(4) = {4, 5};
//+
Line(5) = {2, 4};
//+
Curve Loop(1) = {1, 2, 3, -4, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 2, 3, -4, -5};
//+
Plane Surface(2) = {2};
//+
Recursive Delete {
  Surface{1}; 
}
//+
Recursive Delete {
  Surface{2}; 
}
