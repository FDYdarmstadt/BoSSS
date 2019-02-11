    //two
using BoSSS.Application.XdgPoisson3;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System;
using BoSSS.Platform;
using System.Diagnostics;

System.Console.WriteLine("this is the solution to -div(mu grad(u)) = _4_  with a circular interface and dirichlet boundaries");

XdgPoisson3Control C = new XdgPoisson3Control();

C.DbPath = @"C:\BoSSS\newDB"; \
C.savetodb = true;

C.ProjectName = "XDGPoisson/Circle";
C.ProjectDescription = "XDG Poisson Circle";

//Diffusion Parameter
C.MU_A = 1.0;
C.MU_B = 1000.0;

//Funktion on the RHS of the Poisson-Equation
C.InitialValues.Add("rhs#A", X => 4.0);
C.InitialValues.Add("rhs#B", X => 4.0);


/// Problem Definition
C.FieldOptions.Add("Phi", new AppControl.BaseConfig.FieldOpts() { Degree = 3, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
double RADIUS = 0.75;
C.InitialValues.Add("Phi", X => X[0].Pow2()+X[1].Pow2() - (RADIUS).Pow2());

C.FieldOptions.Add("u", new AppControl.BaseConfig.FieldOpts() { Degree = 3, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });

    //Definition of Grid - Only Dirichlet Boundaries
int Resolution = 13;
 C.GridFunc = delegate{  \
                return Grid2D.Cartesian2DGrid(Grid1D.Linspace(-1, 1, Resolution), Grid1D.Linspace(-1, 1, Resolution)); \
   };

C.g_Diri = (inp => ((inp.X[0].Pow2()+inp.X[1].Pow2())/C.MU_B - (RADIUS).Pow2()/C.MU_B + (RADIUS).Pow2()/C.MU_A)       );

C.IsDirichlet = (inp => true);
//C.g_Neum = (X => 0.0);


/// Exact Solution
C.ExcactSolSupported = true;
C.FieldOptions.Add("uEx", new AppControl.BaseConfig.FieldOpts() { Degree = 3, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
C.InitialValues.Add("uEx#A", X => (X[0].Pow2()+X[1].Pow2())/C.MU_A );
C.InitialValues.Add("uEx#B", X => (X[0].Pow2()+X[1].Pow2())/C.MU_B - (RADIUS).Pow2()/C.MU_B + (RADIUS).Pow2()/C.MU_A  );

/// Solver Parameters
//C.solverName = "pcg+mg+schwarz";
C.solverName = "direct";

/// Discretization Parameters
C.ViscosityMode = Lapalace_Interface.Mode.SWIP;
         
C;
