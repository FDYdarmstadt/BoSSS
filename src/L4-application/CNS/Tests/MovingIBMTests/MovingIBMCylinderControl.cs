using System;
using System.Linq;
using System.Collections.Generic;
using System.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Control;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.Exception;
using CNS.IBM;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.Solution;
using CNS.Source;

string dbPath = @"e:\bosss_db\GridOfTomorrow";
double Mach = 0.2;
int dgDegree = 2;
//double cylinderStart = 2.0; // Causes problems somehow
double cylinderStart = 2.1;
double agglomerationThreshold = 0.3;
int levelSetQuadratureOrder = 10;
Guid gridGuid = new Guid("1b60b59e-2d5a-43dc-9625-dc4c97028153");
int cutCellWeighting = 20;

IBMControl c = new IBMControl();
c.DbPath = dbPath;
c.savetodb = dbPath != null;
			
c.ProjectName = "Moving IBM cylinder, cfl 0.1 restart, small dt";
c.Tags.Add("Moving cylinder");
c.Tags.Add("IBM");
c.Tags.Add("Agglomeration");

c.DomainType = DomainTypes.MovingImmersedBoundary;
c.ActiveOperators = Operators.Convection;
c.EquationOfState = IdealGas.Air;
c.MachNumber = Mach;

c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
c.ExplicitScheme = ExplicitSchemes.RungeKutta;
c.ExplicitOrder = 1;

c.AddVariable(Variables.Density, dgDegree);
c.AddVariable(Variables.Momentum.xComponent, dgDegree);
c.AddVariable(Variables.Momentum.yComponent, dgDegree);
c.AddVariable(Variables.Energy, dgDegree);
c.AddVariable(Variables.Rank, 0);
c.AddVariable(IBMVariables.LevelSet, 2);

c.GridGuid = gridGuid;

c.GridPartType = GridPartType.ParMETIS;
c.GridPartOptions = "10";
c.CutCellWeightingFactor = cutCellWeighting;

double gamma = c.EquationOfState.HeatCapacityRatio;
c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
c.InitialValues_Evaluators.Add(Variables.Velocity[0], X => 0.0);
c.InitialValues_Evaluators.Add(Variables.Velocity[1], X => 0.0);
c.InitialValues_Evaluators.Add(Variables.Pressure, X => 1.0);

double cylinderVelocity = 1.0;
c.LevelSetFunction = delegate (double[] X, double t) { \
	double x = X[0] - cylinderStart - cylinderVelocity * t; \
	double y = X[1]; \
	return x * x + y * y - 1.0; \
};

c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => 0.0);
c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 1.0);

c.AddBoundaryCondition("symmetryPlane");

c.AddBoundaryCondition("adiabaticSlipWall", Variables.Velocity[0], (X, t) => cylinderVelocity);
c.AddBoundaryCondition("adiabaticSlipWall", Variables.Velocity[1], (X, t) => 0.0);
c.LevelSetBoundaryTag = "adiabaticSlipWall";

c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
c.SurfaceHMF_ProjectNodesToLevelSet = false;
c.SurfaceHMF_RestrictNodes = true;
c.SurfaceHMF_UseGaussNodes = false;
c.VolumeHMF_NodeCountSafetyFactor = 3.0;
c.VolumeHMF_RestrictNodes = true;
c.VolumeHMF_UseGaussNodes = false;

c.LevelSetQuadratureOrder = levelSetQuadratureOrder;
c.AgglomerationThreshold = agglomerationThreshold;

c.dtMin = 0.0;
c.dtMax = 1.0;
c.CFLFraction = 0.01;
c.Endtime = 150.0;
c.NoOfTimesteps = int.MaxValue;

c.saveperiod = 200;
c.PrintInterval = 10;

c;