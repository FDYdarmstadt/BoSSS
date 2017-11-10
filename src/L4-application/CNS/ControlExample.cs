using System;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Queries;
using CNS;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using CNS.Solution;
using ilPSP.Utils;

CNSControl c = new CNSControl();

c.DbPath = @"e:\bosss_db\GridOfTomorrow";
int noOfCellsPerDirection = 40;
int dgDegree = 2;
double advectionVelocity = 1.0;

c.savetodb = true;
c.saveperiod = 1;
c.ProjectName = "Isentropic vortex";
c.ProjectDescription = String.Format(
    "Isentropic vortex in a periodic domain with {0}x{0} cells using polynomials of order {1}",
    noOfCellsPerDirection,
    dgDegree);
c.Tags.Add("Isentropic vortex");

c.ActiveOperators = Operators.Convection;
c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
c.TimeSteppingScheme = TimeSteppingSchemes.Explicit;
c.ExplicitScheme = ExplicitSchemes.RungeKutta;
c.ExplicitOrder = 4;

c.EquationOfState = IdealGas.Air;
c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

c.AddVariable(Variables.Density, dgDegree);
c.AddVariable(Variables.Momentum.xComponent, dgDegree);
c.AddVariable(Variables.Momentum.yComponent, dgDegree);
c.AddVariable(Variables.Energy, dgDegree);
c.AddVariable(Variables.Velocity.xComponent, dgDegree);
c.AddVariable(Variables.Velocity.yComponent, dgDegree);
c.AddVariable(Variables.Pressure, dgDegree);
c.AddVariable(Variables.Entropy, dgDegree);

c.GridFunc = delegate {
    double[] nodes = GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1);
    var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: true);
    return grid;
};

double gamma = c.EquationOfState.HeatCapacityRatio;
Func<double[], double, double> x = (X, t) => X[0] - advectionVelocity * t;
Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
Func<double[], double, double> rho = (X, t) => Math.Pow( \
    1.0 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)), \
    1.0 / (gamma - 1.0));
Func<double[], double, double> uAbs = (X, t) => r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma);

c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => advectionVelocity - Math.Sin(phi(X, 0.0)) * uAbs(X, 0.0));
c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => Math.Cos(phi(X, 0.0)) * uAbs(X, 0.0));
c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));

c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, rho));
c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(Variables.Pressure, p));
c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(Variables.Entropy, (X, t) => 1.0));

c.dtMin = 0.0;
c.dtMax = 1.0;
c.CFLFraction = 0.1;
c.Endtime = 2.0;
c.NoOfTimesteps = int.MaxValue;

c;