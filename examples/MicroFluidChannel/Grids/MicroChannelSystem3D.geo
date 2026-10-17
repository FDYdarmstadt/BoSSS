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
halfH = 25.0 * SizeScale;

Linlet = 400.0 * SizeScale;

nSections = 3;
Lsection = Ls + Le;
Lsystem = nSections * Lsection;

Loutlet = 200.0 * SizeScale;
Lend = Lsystem + Loutlet;

Resolution = 2;

// target mesh sizes
li = 20.0 * SizeScale / Resolution;		// mesh size near sharp corners
le = 40.0 * SizeScale / Resolution;		// mesh size near outer walls
lin = 80.0 * SizeScale / Resolution;		// mesh size near inlet
lou = 80.0 * SizeScale / Resolution;		// mesh size near outlet

// ------------------
// mesh construction
// ------------------

// inlet
Point(1) = {0, -halfWe, -halfH, le};
Point(2) = {-Linlet, -halfWe, -halfH, lin};
Point(3) = {-Linlet, halfWe, -halfH, lin};
Point(4) = {0, halfWe, -halfH, le};

Line(1) = {1, 2};
Line(2) = {2, 3};	// inlet
Line(3) = {3, 4}; 

// sections
For n In {1:nSections}
	np0 = newp;
	Point(np0) = {(n-1)*Lsection, halfWi, -halfH, li};
	Point(np0 + 1) = {(n-1)*Lsection + Ls, halfWe, -halfH, le};
	Point(np0 + 2) = {n*Lsection, halfWe, -halfH, le};
	Point(np0 + 3) = {n*Lsection, -halfWe, -halfH, le};
	Point(np0 + 4) = {(n-1)*Lsection + Ls, -halfWe, -halfH, le};
	Point(np0 + 5) = {(n-1)*Lsection, -halfWi, -halfH, li};

	nl0 = newl;
	// upper part
	If(n == 1)
		Line(nl0) = {4, np0};		// connection to upper inlet
	EndIf
	If(n > 1)
		Line(nl0) = {7 + (n-2)*6, np0};		// connection to upper previous section
	EndIf 
	Line(nl0 + 1) = {np0, np0+1};
	Line(nl0 + 2) = {np0+1, np0+2};  		// section upper wall

	// lower part
	Line(nl0 + 3) = {np0+3, np0+4};	// section lower wall
	Line(nl0 + 4) = {np0+4, np0+5};
	If(n == 1)
		Line(nl0 + 5) = {np0+5, 1};		// connection lower inlet
	EndIf
	If(n > 1)
		Line(nl0 + 5) = {np0+5, 8 + (n-2)*6};		// connection lower previous section
	EndIf

EndFor

//outlet
np0 = newp;
Point(np0) = {Lend, halfWe, -halfH, lou};
Point(np0 + 1) =  {Lend, -halfWe, -halfH, lou};

nl0 = newl;
Line(nl0) = {np0-4, np0};		// connection to upper outlet
Line(nl0 + 1) = {np0, np0+1};				// outlet
Line(nl0 + 2) = {np0+1, np0-3};		// sconnection to lower outlet


Curve Loop(1) = {1 : newl-1};

Plane Surface(1) = {1};

ResLayers = {};
numL = Resolution * 2;
For n In {1:numL}
	ResLayers += {1};
EndFor

zLayers = {};
dL = 1.0 / (Resolution * 2.0);
For n In {1:numL}
	zLayers += {n * dL};
EndFor

Extrude {0, 0, 2*halfH} {Surface{1}; Layers{ {ResLayers[]}, {zLayers[]} }; }



