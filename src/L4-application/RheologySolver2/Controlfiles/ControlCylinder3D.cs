using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;

           NSEControl C = new NSEControl();

           // Solver options
           C.Rey = 200;
           C.DegreeVelocity = 2;
           C.DegreePressure = 1;
           C.ViscousScaling = 1;
           C.Tolerance = 1E-6;
           C.MaxSolverIterations = 1;
           C.savetodb = false;
           // C.DbPath = @"P:\BoSSS_DBs\Cylinder3D";
           C.ProjectName = "Cylinder3D";
           C.NoOfTimesteps = 10;
           C.dtMax = 0.1;
           C.dtMin = 0.1;
           C.StreamwisePeriodicBC = false;
           C.Steady = false;

           // Create Fields
           C.FieldOptions.Add("VelocityX", new AppControl.BaseConfig.FieldOpts() { Degree = C.DegreeVelocity, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
           C.FieldOptions.Add("VelocityY", new AppControl.BaseConfig.FieldOpts() { Degree = C.DegreeVelocity, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
           C.FieldOptions.Add("VelocityZ", new AppControl.BaseConfig.FieldOpts() { Degree = C.DegreeVelocity, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
           C.FieldOptions.Add("Pressure", new AppControl.BaseConfig.FieldOpts() { Degree = C.DegreePressure, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });


           // Load Grid
           //  C.GridGuid = new Guid("fd63edae-a816-4922-a82e-93c305c0c11d");   

           // Create Grid
           C.GridFunc = delegate { \

               // x-direction
               var _xNodes1 = Grid1D.ExponentialSpaceing(-9.5, -3, 11, 0.98); \
               _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1)); \
               var _xNodes2 = Grid1D.ExponentialSpaceing(-3, -1, 9, 0.95); \
               _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1)); \
               var _xNodes3 = Grid1D.ExponentialSpaceing(-1, 0, 8, 1);
               _xNodes3 = _xNodes3.GetSubVector(0, (_xNodes3.Length - 1));
               var _xNodes4 = Grid1D.ExponentialSpaceing(0, 2, 9, 1.05);
               _xNodes4 = _xNodes4.GetSubVector(0, (_xNodes4.Length - 1));
               var _xNodes5 = Grid1D.ExponentialSpaceing(2, 8.5, 16, 1.02);
               _xNodes5 = _xNodes5.GetSubVector(0, (_xNodes5.Length - 1));
               var _xNodes6 = Grid1D.ExponentialSpaceing(8.5, 12.5, 5, 1);

               var _xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3, _xNodes4, _xNodes5, _xNodes6);
               //var _xNodes = Grid1D.Linspace(-9.5, 12.5, 16);

               // y-direction
               var _yNodes1 = Grid1D.ExponentialSpaceing(-9, -2.5, 8, 0.91);
               _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
               var _yNodes2 = Grid1D.ExponentialSpaceing(-2.5, -0.5, 8, 0.95);
               _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
               var _yNodes3 = Grid1D.ExponentialSpaceing(-0.5, 0.5, 8, 1.0);
               _yNodes3 = _yNodes3.GetSubVector(0, (_yNodes3.Length - 1));
               var _yNodes4 = Grid1D.ExponentialSpaceing(0.5, 2.5, 8, 1.05);
               _yNodes4 = _yNodes4.GetSubVector(0, (_yNodes4.Length - 1));
               var _yNodes5 = Grid1D.ExponentialSpaceing(2.5, 9, 8, 1.1);

               var _yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3, _yNodes4, _yNodes5);
               //var _yNodes = Grid1D.Linspace(-9, 9, 11);

               // z-direction
               var _zNodes = Grid1D.Linspace(-3, 3, 11);
               //var _zNodes = Grid1D.Linspace(-3, 3, 6);


               // Cut Out
               double[] CutOutPoint1 = new double[3];
               CutOutPoint1[0] = -1;
               CutOutPoint1[1] = -0.5;
               CutOutPoint1[2] = -3;

               double[] CutOutPoint2 = new double[3];
               CutOutPoint2[0] = 0;
               CutOutPoint2[1] = 0.5;
               CutOutPoint2[2] = 3;

               var CutOut = new BoundingBox(3);
               CutOut.AddPoint(CutOutPoint1);
               CutOut.AddPoint(CutOutPoint2);

               var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, true, false, CellType.Cube_Linear, CutOut);

               grd.EdgeTagNames.Add(1, "Velocity_inlet");
               grd.EdgeTagNames.Add(2, "Pressure_Outlet");
               grd.EdgeTagNames.Add(3, "Wall");

               grd.DefineEdgeTags(delegate(double[] _X) {
                   var X = _X;
                   double x = X[0];
                   double y = X[1];
                   double z = X[2];

                   if (Math.Abs(x - (-9.5)) < 1.0e-6)
                       // inlet
                       return 1;

                   if (Math.Abs(x - (12.5)) < 1.0e-6)
                       // outlet
                       return 2;

                   if (Math.Abs(z - (-3)) < 1.0e-6)
                       // left
                       return 2;

                   if (Math.Abs(z - (3)) < 1.0e-6)
                       // right
                       return 2;

                   if (Math.Abs(x - (-1)) < 1.0e-6)
                       // Cube front
                       return 3;

                   if (Math.Abs(x - (0)) < 1.0e-6)
                       // cube back
                       return 3;

                   if (Math.Abs(y - (-0.5)) < 1.0e-6)
                       // cube left 
                       return 3;

                   if (Math.Abs(y - (0.5)) < 1.0e-6)
                       // cube right
                       return 3;

                   throw new ArgumentOutOfRangeException();
               });

               Console.WriteLine("Cylinder3D");

               return grd;
           };


           // Analytic Solution

           // Boundary conditions
           C.AddBoundaryCondition("Velocity_inlet", "VelocityX", X => 1);
           C.AddBoundaryCondition("Velocity_inlet", "VelocityY", X => 0);
           C.AddBoundaryCondition("Velocity_inlet", "VelocityZ", X => 0);
           C.AddBoundaryCondition("Wall");
           C.AddBoundaryCondition("Pressure_Outlet");

           // Set Initial Conditions
           C.InitialValues.Add("VelocityX", X => 0);
           C.InitialValues.Add("VelocityY", X => 0);
           C.InitialValues.Add("VelocityZ", X => 0);
           C.InitialValues.Add("Pressure", X => 0);

           return C;
