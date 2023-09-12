// Gmsh project created on Thu Dec 08 21:49:19 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {0.5, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {1.5, 1, 0, 1.0};
//+
Point(6) = {1.5, 0.178, 0, 1.0};
//+
Line(1) = {6, 2};
//+
Line(2) = {2, 1};
//+
Point(7) = {-0, 0.5, 0, 1.0};
//+
Point(8) = {1.5, 0.5, 0, 1.0};
//+
Point(9) = {1.5, 0.8, 0, 1.0};
//+
Line(3) = {5, 9};
//+
Line(4) = {9, 8};
//+
Line(5) = {8, 6};
//+
Line(6) = {1, 7};
//+
Line(7) = {7, 4};
//+
Line(8) = {4, 4};
//+
Line(8) = {4, 3};
//+
Line(9) = {3, 2};
//+
Line(10) = {3, 5};
//+
Line(11) = {9, 3};
//+
Line(12) = {9, 2};
//+
Line(13) = {8, 2};
//+
Line(14) = {7, 3};
//+
Line(15) = {2, 7};
//+
Curve Loop(1) = {11, 9, -12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, -12, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, 1, -13};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, 3, 11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, -8, -7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {15, 14, 9};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {2, 6, -15};
//+
Plane Surface(7) = {7};
//+
Recursive Delete {
  Surface{2}; 
}
//+
Recursive Delete {
  Surface{1}; 
}
//+
Recursive Delete {
  Curve{13}; 
}
//+
Recursive Delete {
  Surface{3}; 
}
//+
Recursive Delete {
  Curve{11}; 
}
//+
Recursive Delete {
  Curve{11}; 
}
//+
Recursive Delete {
  Surface{4}; 
}
