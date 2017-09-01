cl1 = 1;
Point(1) = {-12, -12, 0, cl1};
Point(2) = {-12, 12, 0, cl1};
Point(3) = {12, 12, 0, cl1};
Point(4) = {12, -12, 0, cl1};
Point(5) = {0, 0, 0, cl1};
Point(6) = {0.5, 0, 0, cl1};
Point(7) = {-0.5, 0, 0, cl1};
Point(8) = {0, 0.5, 0, cl1};
Point(9) = {0, -0.5, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {7, 5, 8};
Circle(6) = {8, 5, 6};
Circle(7) = {6, 5, 9};
Circle(8) = {9, 5, 7};
Line Loop(11) = {5, 6, 7, 8, -1, -4, -3, -2};
Plane Surface(11) = {11};
Delete {
  Line{2, 3, 4};
}
Delete {
  Line{3, 2, 4};
}
Delete {
  Point{3, 4};
}
Delete {
  Line{3};
}
Delete {
  Surface{11};
}
Delete {
  Line{3, 2, 4};
}
Delete {
  Point{3, 4};
}
Point(10) = {25, 12, 0, 1.0};
Point(11) = {25, -12, 0, 1.0};
Line(12) = {2, 10};
Line(13) = {10, 11};
Line(14) = {11, 1};
Line Loop(15) = {13, 14, 1, 12};
Line Loop(16) = {7, 8, 5, 6};
Plane Surface(17) = {15, 16};
Delete {
  Surface{17};
}
Plane Surface(17) = {15, 16};
Delete {
  Surface{17};
}
Delete {
  Line{13, 12, 14};
}
Delete {
  Point{10, 11};
}
Point(11) = {12, 12, 0.5, 1.0};
Point(12) = {12, -12, 0.5, 1.0};
Line(17) = {1, 12};
Line(18) = {12, 11};
Line(19) = {11, 2};
Line Loop(20) = {19, -1, 17, 18};
Plane Surface(21) = {16, 20};
