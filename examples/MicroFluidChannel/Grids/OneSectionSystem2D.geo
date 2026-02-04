Mesh.MshFileVersion = 2.2;
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

Resolution = 1;

// target mesh sizes
li = 20.0 * SizeScale / Resolution;		// mesh size near sharp corners
le = 40.0 * SizeScale / Resolution;		// mesh size near outer walls
lin = 60.0 * SizeScale / Resolution;		// mesh size near inlet
lou = 60.0 * SizeScale / Resolution;		// mesh size near outlet

// ------------------
// mesh construction
// ------------------

// inlet
Point(1) = {0, -halfWe, 0, le};
Point(2) = {-Linlet, -halfWe, 0, lin};
Point(3) = {-Linlet, halfWe, 0, lin};
Point(4) = {0, halfWe, 0, le};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4}; 

//outlet
Point(5) = {Lend, halfWe, 0, lou};
Point(6) =  {Lend, -halfWe, 0, lou};

Line(4) = {5, 6};

// sections
nP0 = 6;
nL0 = 4;
For n In {1:nSections}
	Point(nP0 + 1) = {0, halfWi, 0, li};
	Point(nP0 + 2) = {Ls, halfWe, 0, le};
	Point(nP0 + 3) = {Lsection, halfWe, 0, le};
	Point(nP0 + 4) = {Lsection, -halfWe, 0, le};
	Point(nP0 + 5) = {Ls, -halfWe, 0, le};
	Point(nP0 + 6) = {0, -halfWi, 0, li};

	Line(nL0 + 1) = {4, 7};		// connection to upper inlet
	Line(nL0 + 2) = {7, 8};
	Line(nL0 + 3) = {8, 9};  		// section upper wall
	Line(nL0 + 4) = {9, 5};		// connection to upper outlet
	Line(nL0 + 5) = {6, 10};		// sconnection to lower outlet
	Line(nL0 + 6) = {10, 11};	// section lower wall
	Line(nL0 + 7) = {11, 12};
	Line(nL0 + 8) = {12, 1};		// connection lower inlet
EndFor

Curve Loop(1) = {1 : 12};

Plane Surface(1) = {1};




