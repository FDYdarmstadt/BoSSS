// -----------------------------------------
//
//	Sharp corner micro channel section
//
// -----------------------------------------

SizeScale = 1.0e-6;
Ls = 80.0 * SizeScale;
Le = 200.0 * SizeScale;
halfWi = 21.25 * SizeScale;
halfWe = 100.0 * SizeScale;
H = 50.0 * SizeScale;

Linlet = 200.0 * SizeScale;

nSections = 1;
Lsection = Ls + Le;
Lsystem = nSections * Lsection;

Loutlet = 200.0 * SizeScale;
Lend = Lsystem + Loutlet;


// target mesh sizes
li = 10.0 * SizeScale;		// mesh size near sharp corners
le = 20.0 * SizeScale;		// mesh size near outer walls

Point(1) = {0, -halfWi, 0, li};
Point(2) = {0, halfWi, 0, li};
Point(3) = {Ls, halfWe, 0, le};
Point(4) = {Lsection, halfWe, 0, le};
Point(5) = {Lsection, -halfWe, 0, le};
Point(6) = {Ls, -halfWe, 0, le};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = {1};




