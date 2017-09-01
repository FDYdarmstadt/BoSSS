using System;
using System.Collections.Generic;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using CNS;

var c = new CNSControl();

// We load the input data from some function within the Framework. 
// The class ControlExamples provides some predefined control sets for 
// different typical test cases, i.a. the isentropic vortex.
// For the isentropic vortex, you have to specify 
// - path to a database (here: EMPTY path)
// - number of cells in each direction (here: 20)
// - DG order (here: 2)
// - advection velcoity of the vortex (here: 1.0)
int noOfCellsPerDirection = 20;
c = ControlExamples.IsentropicVortex(null,noOfCellsPerDirection,2,1.0);



// Link a database and save to database
c.DbPath = @"c:\tmp\CNS_db";
c.savetodb = true;  

// Creating the grid
c.GridFunc = delegate {       
   double[] nodes = GenericBlas.Linspace(-10, 10, noOfCellsPerDirection + 1);       
   var grid = Grid2D.Cartesian2DGrid(nodes, nodes,periodicX: true, periodicY: true);       
   return grid;       
};

c;