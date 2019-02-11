    //two
using BoSSS.Application.XdgPoisson3;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;


XdgPoisson3Control C = new XdgPoisson3Control();

C.DbPath = @"D:\\BoSSS-db"; \
C.savetodb = false;

C.ProjectName = "XDGPoisson/Circle";
C.ProjectDescription = "XDG Poisson Circle";

//Diffusion Parameter
C.MU_A = -1.0;
C.MU_B = -1.0;

//Funktion on the RHS of the Poisson-Equation
C.InitialValues.Add("rhs#A", X => Math.Sin((X[0] + 1.0)*3));
C.InitialValues.Add("rhs#B", X => Math.Sin((X[0] + 1.0)*3));


/// Problem Definition
C.FieldOptions.Add("Phi", new AppControl.BaseConfig.FieldOpts() { Degree = 3, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
double RADIUS = 0.8;
C.InitialValues.Add("Phi", X => X[0].Pow2()+X[1].Pow2() - (RADIUS).Pow2());

C.FieldOptions.Add("u", new AppControl.BaseConfig.FieldOpts() { Degree = 3, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });

    //Definition of Grid - Only Dirichlet Boundaries
int Resolution = 13;
 C.GridFunc = delegate() {                                                                                 \
    return Grid2D.Cartesian2DGrid(Grid1D.Linspace(-1, 1, Resolution), Grid1D.Linspace(-1, 1, Resolution)); \
   };

C.g_Diri = (inp => 0.0);
C.g_Neum = delegate(CommonParamsBnd inp) {                                                  \
                if(Math.Abs(inp.X[1] - 1.0) < 1.0e-8 || Math.Abs(inp.X[1] + 1.0) < 1.0e-8)  \
                    return 0;                                                               \
                return Math.Cos(2.0*3.0)*(1.0/3.0);                                                       \
            };
C.IsDirichlet = (inp => Math.Abs(inp.X[0] + 1.0) <= 1.0e-8);


/// Exact Solution
C.ExcactSolSupported = true;
C.FieldOptions.Add("uEx", new AppControl.BaseConfig.FieldOpts() { Degree = 3, SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE });
C.InitialValues.Add("uEx#A", X => Math.Sin((X[0] + 1.0)*3)*(1.0/9.0));
C.InitialValues.Add("uEx#B", X => Math.Sin((X[0] + 1.0)*3)*(1.0/9.0));

/// Solver Parameters
//C.solverName = "pcg+mg+schwarz";
C.solverName = "direct";

/// Discretization Parameters
C.ViscosityMode = Laplace_Interface.Mode.SIP;
         
C;
