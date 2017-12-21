/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.LoadBalancing;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.ShockCapturing;
using CNS.Source;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CNS {

    /// <summary>
    /// A set of exemplary control files for CNS. You can try the individual
    /// examples by executing
    /// <code>
    /// ./CNS.exe -c cs:CNS.ControlExamples.ExampleName(options)
    /// </code>
    /// on the console. Using the isentropic vortex example, the calling
    /// sequence could be
    /// <code>
    /// ./CNS.exe -c cs:CNS.ControlExamples.IsentropicVortex(@\"c:\bosss_db\",20,2,1.0)
    /// </code>
    /// or, for parallel runs,
    /// <code>
    /// mpiexec -n 2 ./CNS.exe -c "cs:CNS.ControlExamples.IsentropicVortex(@\"c:\bosss_db\",20,2,1.0)"
    /// </code>
    /// (note the additional quotes).
    /// </summary>
    public static class ControlExamples {

        /// <summary>
        /// Isentropic vortex in a fully periodic domain.
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <param name="dgDegree"></param>
        /// <param name="advectionVelocity"></param>
        /// <returns></returns>
        public static CNSControl IsentropicVortex(string dbPath, int noOfCellsPerDirection, int dgDegree, double advectionVelocity) {
            CNSControl c = new CNSControl();

            c.dtMin = 0.0;
            c.dtMax = 5.0e-1;
            c.CFLFraction = 0.9;
            c.Endtime = 2.0;
            c.NoOfTimesteps = int.MaxValue;

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1000;
            c.PrintInterval = 10;
            c.ProjectDescription = String.Format(
                "Isentropic vortex in a periodic domain with {0}x{0} cells using polynomials of order {1}",
                noOfCellsPerDirection,
                dgDegree);
            c.Tags.Add("Isentropic vortex");

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 0.2;

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-10, 10, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: true);
                //var grid = Grid2D.HangingNodes2D(true, true,
                //                new GridBox(-10.0, -10.0, 10.0, 10.0, noOfCellsPerDirection),
                //                new GridBox(-10.0, -2.0, 10.0, 2.0, (int)noOfCellsPerDirection * 2, (int)noOfCellsPerDirection * 8 / 20)
                //                );
                return grid;
            };

            double gamma = c.EquationOfState.HeatCapacityRatio;
            double Umax = 0.25;
            double MachScaling = gamma * c.MachNumber * c.MachNumber;
            Func<double[], double, double> x = delegate (double[] X, double t) {
                double x1 = X[0] - advectionVelocity * t;
                while (x1 < -10.0) {
                    x1 += 20.0;
                }
                return x1;
            };
            Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
            Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
            Func<double[], double, double> rho = (X, t) => Math.Pow(
                1.0 - 0.5 * MachScaling * (gamma - 1.0) / gamma * Umax * Umax * Math.Exp(1.0 - r(X, t) * r(X, t)),
                1.0 / (gamma - 1.0));
            Func<double[], double, double> uAbs = (X, t) => Umax * r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
            Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma);
            Func<double[], double, double> u = (X, t) => advectionVelocity - Math.Sin(phi(X, t)) * uAbs(X, t);
            Func<double[], double, double> v = (X, t) => Math.Cos(phi(X, t)) * uAbs(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));


            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, rho));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(Variables.Pressure, p));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(Variables.Velocity.xComponent, u));
            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(Variables.Entropy, (X, t) => 1.0));

            c.ProjectName = "IsentropivVortex_Mach=" + c.MachNumber + "_cells=" + noOfCellsPerDirection + "_p=" + dgDegree + "_flux=" + c.ConvectiveFluxType;

            return c;
        }

        /// <summary>
        /// Exemplary complete hp parameter study
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="noOfRefinements"></param>
        /// <param name="maxDegree"></param>
        /// <param name="minDegree"></param>
        /// <param name="advectionVelocity"></param>
        /// <returns></returns>
        public static CNSControl[] IsentropicVortexStudy(string dbPath, int noOfRefinements, int maxDegree, double advectionVelocity, int minDegree = 2) {
            CNSControl[] controls = new CNSControl[(maxDegree + 1 - minDegree) * noOfRefinements];
            int ii = 0;
            for (int i = minDegree; i <= maxDegree; i++) {
                for (int j = 0; j < noOfRefinements; j++) {
                    double power = 3 + (double)j;
                    int noOfCellsPerDirection = (int)Math.Pow(2, power);
                    controls[ii] = IsentropicVortex(dbPath, noOfCellsPerDirection, i, advectionVelocity);
                    controls[ii].Paramstudy_ContinueOnError = true;
                    controls[ii].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                    new Tuple<string, object>("dgDegree", i),
                    new Tuple<string, object>("noOfCellsPerDirection", noOfCellsPerDirection)
                    };
                    ii++;
                }
            }
            return controls;
        }

        /// <summary>
        /// A Gaussian pulse in uniform flow with simple non-reflecting
        /// boundary conditions
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <param name="dgDegree"></param>
        /// <param name="advectionVelocity"></param>
        /// <param name="dampingWidth"></param>
        /// <param name="dampingStrength"></param>
        /// <returns></returns>
        public static CNSControl GaussianPulse(string dbPath, int noOfCellsPerDirection, int dgDegree, double advectionVelocity, double dampingWidth = 3.0, double dampingStrength = 10.0) {
            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = true;
            c.saveperiod = 10;
            c.ProjectName = "Gaussian pulse";
            c.ProjectDescription = String.Format(
                "Gaussian pulse in a domain with non-reflecting bcs. Uses {0}x{0} cells and polynomials of order {1}",
                noOfCellsPerDirection,
                dgDegree);
            c.Tags.Add("Gaussian pulse");

            c.ActiveOperators = Operators.Convection | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            double extent = 10.0;
            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-extent, extent, noOfCellsPerDirection + 1);
                Grid2D grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                //grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.EdgeTagNames.Add(1, "subsonicOutlet");
                grid.DefineEdgeTags(x => 1);
                return grid;
            };

            StateVector refState = StateVector.FromPrimitiveQuantities(
                new Material(c),
                1.0,
                new Vector3D(advectionVelocity, 0.0, 0.0),
                8.0);

            Func<double[], double> pulse =
                X => (1 + 0.01 * Math.Exp(-(X[0] * X[0] + X[1] * X[1]) / 2 / 0.8 / 0.8));
            c.InitialValues_Evaluators.Add(Variables.Density, X => refState.Density * pulse(X));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => refState.Velocity.x);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => refState.Velocity.y);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => refState.Pressure * pulse(X));

            c.AddBoundaryCondition("subsonicOutlet", Variables.Pressure, (X, t) => refState.Pressure);

            double dampingStart = extent - dampingWidth;
            Func<double[], double> quadraticDamping = delegate (double[] x) {
                double s = (Math.Abs(x[0]) - dampingStart) / dampingWidth;
                if (s <= 0.0) {
                    return 0.0;
                }

                return dampingStrength * s * s;
            };
            c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                (x, t, state) => quadraticDamping(x) * (state.Density - refState.Density)));
            c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                (x, t, state) => quadraticDamping(x) * (state.Momentum.x - refState.Momentum.x)));
            c.CustomMomentumSources[1].Add(map => new AdHocSourceTerm(map,
                (x, t, state) => quadraticDamping(x) * (state.Momentum.y - refState.Momentum.y)));
            c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                (x, t, state) => quadraticDamping(x) * (state.Energy - refState.Energy)));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.2;
            c.Endtime = 4.0;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }

        /// <summary>
        /// Flow around a cylinder using a boundary-fitted grid
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="dgDegree"></param>
        /// <param name="grid"></param>
        /// <returns></returns>
        public static CNSControl CylinderFlow(string dbPath, int dgDegree, string grid) {
            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = true;
            c.saveperiod = 500;

            c.ProjectName = "Cylinder flow";
            c.Tags.Add("Cylinder flow");
            c.Tags.Add("Re=20");
            c.Tags.Add("Mach=0.2");

            c.ActiveOperators = Operators.Convection | Operators.Diffusion;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
            c.SIPGPenaltyScaling = 1.0;

            c.ExplicitScheme = ExplicitSchemes.LTS;
            c.ExplicitOrder = 2;
            c.NumberOfSubGrids = 3;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 136.37563;
            c.PrandtlNumber = 0.72;
            c.ViscosityLaw = new PowerLaw();

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Temperature, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);

            //Karman_Mesh[-5,10]x[-5,5]
            //bilinear, r=48elements, total=2260:
            //9a3a21dd-45b6-42c6-b028-195ac073d52b

            //KarmanMesh
            //fine bilinear
            //9cbe1f99-5d5c-4212-aefd-67b9b4b464a9
            //coarse bilinear
            //5c0bcb24-1ec9-4f33-ae88-018bd74a7e66

            //coarsefine bilinear
            //5b60f9ec-9b73-4979-ba89-aa96728eb37a

            //coarsfine2 biquadtratic

            //coarse quadratic
            //f0e17dc5-0c60-45b0-86d1-4ac69898d0de

            //Tria Mesh
            //coarse bilinear: a87955bd-f836-4cce-9d7c-0a9bef08299e

            Guid gridGuid;
            if (Guid.TryParse(grid, out gridGuid)) {
                c.GridGuid = gridGuid;
            } else {
                c.GridFunc = delegate () {
                    GridCommons _grid;
                    if (File.Exists(grid)) {
                        _grid = GridImporter.Import(grid);
                    } else {
                        throw new Exception(
                            "Could not find a grid at " + grid);
                    }

                    _grid.EdgeTagNames.Add(1, "supersonicinlet");
                    _grid.EdgeTagNames.Add(2, "subsonicinlet");
                    _grid.EdgeTagNames.Add(3, "subsonicoutlet");
                    _grid.EdgeTagNames.Add(4, "adiabaticWall");

                    _grid.DefineEdgeTags(delegate (double[] X) {
                        throw new NotImplementedException();
                    });

                    return _grid;
                };
            }

            c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Momentum[0], X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Momentum[1], X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Energy, X => 45.14285714);

            c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 17.85714286);

            c.AddBoundaryCondition("subsonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("subsonicInlet", Variables.Velocity[0], (X, t) => 1.0);
            c.AddBoundaryCondition("subsonicInlet", Variables.Velocity[1], (X, t) => 0.0);

            c.AddBoundaryCondition("subsonicOutlet", Variables.Pressure, (X, t) => 17.85714286);

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.4;
            c.Endtime = 10.0;
            c.NoOfTimesteps = int.MaxValue;
            c.ResidualInterval = 100;
            c.PrintInterval = 10;

            return c;
        }

        /// <summary>
        /// An isentropic vortex in a domain with a straight immersed boundary.
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="dgDegree"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <param name="levelSetPosition"></param>
        /// <param name="agglomerationThreshold"></param>
        /// <returns></returns>
        public static IBMControl IBMIsentropicVortex(string dbPath, int dgDegree = 2, int noOfCellsPerDirection = 20, double levelSetPosition = -0.4, double agglomerationThreshold = 0.5) {
            IBMControl c = new IBMControl();

            double advectionVelocity = 1.0;

            c.DbPath = dbPath;
            c.savetodb = true;
            c.saveperiod = 10;
            c.ProjectName = "IBM Isentropic vortex";
            c.ProjectDescription = String.Format(
                "Isentropic vortex in a periodic domain (x-direction) bounded"
                    + " by an immersed boundary at the bottom. Uses {0}x{0}"
                    + " cells using polynomials of order {1}",
                noOfCellsPerDirection,
                dgDegree);
            c.Tags.Add("Isentropic vortex");
            c.Tags.Add("IBM");

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "adiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 10;
            c.LevelSetBoundaryTag = "supersonicInlet";

            double gamma = c.EquationOfState.HeatCapacityRatio;
            Func<double[], double, double> x = (X, t) => X[0] - advectionVelocity * t;
            Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
            Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
            Func<double[], double, double> rho = (X, t) => Math.Pow(
                1.0 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)),
                1.0 / (gamma - 1.0));
            Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma);
            Func<double[], double, double> uAbs = (X, t) => r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
            Func<double[], double, double> u = (X, t) => advectionVelocity - Math.Sin(phi(X, t)) * uAbs(X, t);
            Func<double[], double, double> v = (X, t) => Math.Cos(phi(X, t)) * uAbs(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - levelSetPosition;

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], u);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], v);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, p);

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(Variables.Density, rho));
            c.Queries.Add("L2ErrorPressure", IBMQueries.L2Error(state => state.Pressure, p));
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.2;
            c.Endtime = 1;
            c.NoOfTimesteps = 10;

            return c;
        }

        /// <summary>
        /// Flow around an immersed cylinder defined by a quadratic level set
        /// function
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="Mach"></param>
        /// <param name="refinements"></param>
        /// <param name="dgDegree"></param>
        /// <param name="agglomerationThreshold"></param>
        /// <param name="levelSetQuadratureOrder"></param>
        /// <returns></returns>
        public static IBMControl IBMCylinder(string dbPath = null, double Mach = 0.2, int refinements = 0, int dgDegree = 2, double agglomerationThreshold = 0.1, int levelSetQuadratureOrder = 10) {
            IBMControl c = new IBMControl();
            c.DbPath = dbPath;
            c.savetodb = dbPath != null;

            c.ProjectName = "IBM cylinder";
            c.ProjectDescription = String.Format(
                "Flow around cylinder represented by a level set at Mach {0}" +
                    " with cell agglomeration threshold {1} and {2}th order" +
                    " HMF quadrature (classic variant)",
                Mach,
                agglomerationThreshold,
                levelSetQuadratureOrder);

            c.Tags.Add("Cylinder");
            c.Tags.Add("IBM");
            c.Tags.Add("Agglomeration");

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 2);

            c.GridFunc = delegate {
                // Refined Region
                double refinedWidthX = 20.0;
                double refinedWidthY = 20.0;

                int baseNumberOfCells = 16 << refinements;
                double spacingFactor = 2.0;
                double D = 40.0;

                double[] xLeft = Grid1D.TanhSpacing(-D, -0.5 * refinedWidthX, baseNumberOfCells / 2 + 1, spacingFactor, false);
                double[] xCenter = GenericBlas.Linspace(-0.5 * refinedWidthX, 0.5 * refinedWidthX, 2 * baseNumberOfCells + 1);
                double[] xRight = Grid1D.TanhSpacing(0.5 * refinedWidthX, D, baseNumberOfCells / 2 + 1, spacingFactor, true);
                double[] xNodes = xLeft.Take(xLeft.Length - 1).Concat(xCenter.Take(xCenter.Length - 1)).Concat(xRight).ToArray();

                double[] yCenter = GenericBlas.Linspace(0, 0.5 * refinedWidthY, baseNumberOfCells + 1);
                double[] yTop = Grid1D.TanhSpacing(0.5 * refinedWidthY, D, baseNumberOfCells / 2 + 1, spacingFactor, true);
                double[] yNodes = yCenter.Take(yCenter.Length - 1).Concat(yTop).ToArray();

                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.EdgeTagNames.Add(2, "adiabaticSlipWall");
                Func<double[], byte> func = delegate (double[] x) {
                    if (Math.Abs(x[1]) < 1e-14) {
                        return 2;
                    } else {
                        return 1;
                    }
                };
                grid.DefineEdgeTags(func);

                grid.Name = String.Format("Refined region: {0}x{1}", 2 * baseNumberOfCells, baseNumberOfCells);

                return grid;
            };

            c.GridPartType = GridPartType.ParMETIS;
            c.GridPartOptions = "10";

            double gamma = c.EquationOfState.HeatCapacityRatio;
            c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity[0], X => Mach * Math.Sqrt(gamma));
            c.InitialValues_Evaluators.Add(Variables.Velocity[1], X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 1.0);

            c.LevelSetFunction = (X, t) => X[0] * X[0] + X[1] * X[1] - 1.0 * 1.0;

            c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => Mach * Math.Sqrt(gamma));
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 1.0);

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.LevelSetBoundaryTag = "adiabaticSlipWall";

            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
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
            c.CFLFraction = 0.1;
            c.Endtime = double.MaxValue;
            c.NoOfTimesteps = int.MaxValue;

            c.saveperiod = 1000;
            c.PrintInterval = 1;

            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("query_L2ErrorEntropy_change", 1e-13);

            return c;
        }

        /// <summary>
        /// Parameter study for the flow around an immersed cylinder
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="maxRefinements"></param>
        /// <param name="maxDegree"></param>
        /// <param name="Mach"></param>
        /// <param name="agglomerationThreshold"></param>
        /// <param name="levelSetQuadratureOrder"></param>
        /// <returns></returns>
        public static IBMControl[] IBMCylinderStudy(string dbPath = @"\\fdyprime\userspace\mueller\cluster\db", double maxRefinements = 2, double maxDegree = 4, double Mach = 0.2, double agglomerationThreshold = 0.3, int levelSetQuadratureOrder = 10) {
            int parmetiesRefinements = 10;

            List<IBMControl> controls = new List<IBMControl>();
            int i = 1;
            for (int refinements = 0; refinements <= maxRefinements; refinements++) {
                for (int dgDegree = 0; dgDegree <= maxDegree; dgDegree++) {
                    IBMControl c = IBMCylinder(dbPath, Mach, refinements, dgDegree, agglomerationThreshold, levelSetQuadratureOrder);

                    c.ProjectName += String.Format(" case {0}, agg {1}", i, agglomerationThreshold);

                    c.CFLFraction = 0.2;
                    c.Endtime = 5000.0;
                    c.NoOfTimesteps = int.MaxValue;

                    c.saveperiod = 2500;
                    c.PrintInterval = 20;
                    c.ResidualInterval = 200;

                    c.Paramstudy_ContinueOnError = true;
                    c.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("dgDegree", dgDegree),
                        new Tuple<string, object>("refinements", refinements),
                    };

                    c.GridPartType = GridPartType.ParMETIS;
                    c.GridPartOptions = parmetiesRefinements.ToString();

                    controls.Add(c);
                    i++;
                }
            }

            return controls.ToArray();
        }

        public static IBMControl NavierStokes(string dbPath = null, double Mach = 0.2, int dgDegree = 2, double agglomerationThreshold = 0.1, int levelSetQuadratureOrder = 10) {
            IBMControl c = new IBMControl();
            c.DbPath = dbPath;
            c.savetodb = dbPath != null;

            c.ProjectName = "Navier-Stokes Test";
            c.ProjectDescription = String.Format(
                "Bla",
                Mach,
                agglomerationThreshold,
                levelSetQuadratureOrder);

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection | Operators.Diffusion;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 136.37563;
            c.PrandtlNumber = 0.72;
            c.ViscosityLaw = new PowerLaw();

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 2);

            c.GridFunc = delegate () {
                int numberOfCells = 64;
                double stretching = 2.0;
                double D = 40.0;

                double[] nodesA = Grid1D.TanhSpacing(-D, 0, numberOfCells / 2 + 1, stretching, false);
                double[] nodesB = Grid1D.TanhSpacing(0, D, numberOfCells / 2 + 1, stretching, true);
                double[] nodes = nodesA.Take(nodesA.Length - 1).Concat(nodesB).ToArray();
                GridCommons grid = Grid2D.Cartesian2DGrid(nodes, nodes);

                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.EdgeTagNames.Add(2, "adiabaticWall");
                Func<double[], byte> func = delegate (double[] x) {
                    return 1;
                };
                grid.DefineEdgeTags(func);

                return grid;
            };

            c.GridPartType = GridPartType.ParMETIS;
            c.GridPartOptions = "5";

            c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Momentum[0], X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Momentum[1], X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Energy, X => 45.14285714);

            c.LevelSetFunction = (X, t) => X[0] * X[0] + X[1] * X[1] - 1.0 * 1.0;

            c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 17.85714286);

            c.AddBoundaryCondition("adiabaticWall");
            c.LevelSetBoundaryTag = "adiabaticWall";

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = levelSetQuadratureOrder;
            c.AgglomerationThreshold = agglomerationThreshold;

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = double.MaxValue;
            c.NoOfTimesteps = int.MaxValue;

            c.saveperiod = 100;
            c.PrintInterval = 10;

            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("query_L2ErrorEntropy_change", 1e-12);

            return c;
        }

        public static IBMControl IBMSphere(string dbPath = null, double Mach = 0.2, int refinements = 0, int dgDegree = 2, double agglomerationThreshold = 0.1, int levelSetQuadratureOrder = 6) {
            IBMControl c = new IBMControl();
            c.DbPath = dbPath;
            c.savetodb = dbPath != null;

            c.ProjectName = "IBM sphere";
            c.ProjectDescription = String.Format(
                "Flow around a sphere represented by a level set at Mach {0}" +
                    " with cell agglomeration threshold {1} and {2}th order" +
                    " HMF quadrature (classic variant)",
                Mach,
                agglomerationThreshold,
                levelSetQuadratureOrder);

            c.Tags.Add("Sphere");
            c.Tags.Add("IBM");
            c.Tags.Add("Agglomeration");

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 2);

            c.GridFunc = delegate () {
                int numberOfCells = 32 << refinements;
                double stretching = 2.0;
                double D = 40.0;

                double[] nodesA = Grid1D.TanhSpacing(-D, 0, numberOfCells / 2 + 1, stretching, false);
                double[] nodesB = Grid1D.TanhSpacing(0, D, numberOfCells / 2 + 1, stretching, true);
                double[] nodes = nodesA.Take(nodesA.Length - 1).Concat(nodesB).ToArray();
                GridCommons grid = Grid3D.Cartesian3DGrid(nodes, nodesB, nodesB);

                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.EdgeTagNames.Add(2, "adiabaticSlipWall");
                Func<double[], byte> func = delegate (double[] x) {
                    if (Math.Abs(x[1]) < 1e-14 || Math.Abs(x[2]) < 1e-14) {
                        return 2;
                    } else {
                        return 1;
                    }
                };
                grid.DefineEdgeTags(func);

                return grid;
            };

            c.GridPartType = GridPartType.ParMETIS;
            c.GridPartOptions = "5";

            double gamma = c.EquationOfState.HeatCapacityRatio;
            c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity[0], X => Mach * Math.Sqrt(gamma));
            c.InitialValues_Evaluators.Add(Variables.Velocity[1], X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 1.0);
            c.InitialValues_Evaluators.Add(IBMVariables.LevelSet, X => X[0] * X[0] + X[1] * X[1] - 1.0 * 1.0);

            c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => Mach * Math.Sqrt(gamma));
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 1.0);

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.LevelSetBoundaryTag = "adiabaticSlipWall";

            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            c.SurfaceHMF_ProjectNodesToLevelSet = true;
            c.SurfaceHMF_RestrictNodes = true;
            c.VolumeHMF_RestrictNodes = true;
            c.VolumeHMF_NodeCountSafetyFactor = 5.0;

            c.LevelSetQuadratureOrder = levelSetQuadratureOrder;
            c.AgglomerationThreshold = agglomerationThreshold;

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = double.MaxValue;
            c.NoOfTimesteps = int.MaxValue;

            c.saveperiod = 500;
            c.PrintInterval = 1;

            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("query_L2ErrorEntropy_change", 1e-12);

            return c;
        }

        /// <summary>
        /// Isentropic Vortex in an small channel with analytical boundary condition at top and bottom and periodic in x direction.
        /// It is used to perform a convergence study for the Local Time Stepping algorithm in the immersed boundary case
        /// </summary>
        /// <returns></returns>
        public static IBMControl IBMIsentropicVortex_LTS(double dtMax, int cycles = 120, double advectionVelocity = 50) {
            IBMControl c = new IBMControl();

            int dgDegree = 4;
            c.CFLFraction = 0.25;
            int unitResolution = 8;

            double BoxSizeX = 8;
            double BoxSizeY = 2;
            int nodesInX = (int)(BoxSizeX * unitResolution);
            int nodesInY = (int)(BoxSizeY * unitResolution);

            double levelSetPosition = -1 + 1 / (double)(2 * unitResolution);

            // Time Stepping
            c.ExplicitScheme = ExplicitSchemes.LTS;
            c.ExplicitOrder = 3;
            c.dtMin = 0.0;
            c.dtMax = dtMax;

            // Solver Settings
            c.Endtime = cycles * BoxSizeX / advectionVelocity;
            c.NoOfTimesteps = int.MaxValue;

            c.PrintInterval = 100;
            c.ResidualInterval = 10000;
            c.ResidualLoggerType = ResidualLoggerTypes.Query;

            c.DbPath = @"c:\bosss_dbv2\cutVortex";
            c.savetodb = true;
            c.saveperiod = 50000;
            c.ProjectName = "[" + BoxSizeX + "x" + BoxSizeY + "]_(" + nodesInX + "x" + nodesInY + ")_vel=" + advectionVelocity + "_cyles=" + cycles + "_CFL=" + c.CFLFraction + "_p=" + dgDegree + "_" + c.ExplicitScheme + "-" + c.ExplicitOrder;
            c.ProjectDescription = String.Format(
                "Isentropic vortex in a periodic domain (x-direction) bounded"
                    + " by an immersed boundary at the bottom. Uses {0}x{1}"
                    + " cells using polynomials of order {2}",
                nodesInX,
                nodesInY,
                dgDegree);
            c.Tags.Add("Isentropic vortex");
            c.Tags.Add("IBM");


            //IBM Settings
            c.LevelSetQuadratureOrder = 12;
            c.AgglomerationThreshold = 0.3;
            c.GridPartOptions = "0";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;


            // Solver Type
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            //Material Settings
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);


            // Grid
            c.GridFunc = delegate {
                double[] nodesX = GenericBlas.Linspace(-BoxSizeX / 2.0, BoxSizeX / 2.0, nodesInX + 1);
                double[] nodesY = GenericBlas.Linspace(-BoxSizeY / 2.0, BoxSizeY / 2.0, nodesInY + 1);
                var grid = Grid2D.Cartesian2DGrid(nodesX, nodesY, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };
            c.GridPartType = GridPartType.ParMETIS;

            double gamma = c.EquationOfState.HeatCapacityRatio;
            Func<double[], double, double> x = delegate (double[] X, double t) {
                double x1 = X[0] - advectionVelocity * t;
                while (x1 < -BoxSizeX / 2.0) {
                    x1 += BoxSizeX;
                }
                return x1;
            };
            Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
            Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
            Func<double[], double, double> rho = (X, t) => Math.Pow(1.0 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)), 1.0 / (gamma - 1.0));
            Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma);
            Func<double[], double, double> uAbs = (X, t) => r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
            Func<double[], double, double> u = (X, t) => advectionVelocity - Math.Sin(phi(X, t)) * uAbs(X, t);
            Func<double[], double, double> v = (X, t) => Math.Cos(phi(X, t)) * uAbs(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));
            c.InitialValues_Evaluators.Add(IBMVariables.LevelSet, X => X[1] - levelSetPosition);

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], u);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], v);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, p);
            c.LevelSetBoundaryTag = "supersonicInlet";

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(Variables.Density, rho));
            c.Queries.Add("L2ErrorPressure", IBMQueries.L2Error(state => state.Pressure, p));
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

            return c;
        }

        public static CNSControl ShockTube(string dbPath = null, int dgDegree = 2, int numOfCellsX = 50, int numOfCellsY = 1, double sensorLimit = 1e-3, bool true1D = false, bool saveToDb = false) {

            CNSControl c = new CNSControl();

            // Load balancing
            //c.DynamicLoadBalancing_CellCostEstimatorFactory = delegate (IApplication<AppControl> app, int performanceClassCount) {
            //    if (performanceClassCount != 2) {
            //        throw new ConfigurationException();
            //    }

            //    int[] performanceClassToCostMap = new int[] { 1, 10 };
            //    return new StaticCellCostEstimator(performanceClassToCostMap);
            //};
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = 10;

            dbPath = @"c:\bosss_db\";
            //dbPath = @"e:\bosss_db\GridOfTomorrow\";
            //dbPath = @"\\fdyprime\userspace\geisenhofer\bosss_db\";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = 10;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            bool AV = true;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            }

            // Runge-Kutta schemes
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            //Adams-Bashforth
            //c.ExplicitScheme = ExplicitSchemes.AdamsBashforth;
            //c.ExplicitOrder = 3;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 3;
            //c.NumberOfSubGrids = 1;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            if (true1D == false) {
                c.AddVariable(Variables.Momentum.yComponent, dgDegree);
                c.AddVariable(Variables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(Variables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(Variables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);

                if (true1D) {
                    var grid = Grid1D.LineGrid(xNodes, periodic: false);
                    // Boundary conditions
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });
                    return grid;
                } else {
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                    // Boundary conditions
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });
                    return grid;
                }
            };

            c.AddBoundaryCondition("AdiabaticSlipWall");

            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector2D normalVector = new Vector2D(1, 0);

            // Direction vector of initial shock
            Vector2D r = new Vector2D(normalVector.y, -normalVector.x);
            r.Normalize();

            // Distance from a point X to the initial shock
            double[] p = new double[] { 0.5, 0.0 };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            Func<double, double> Jump = (x => x <= 0.5 ? 0 : 1);

            // Initial conditions
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;

            //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            if (true1D == false) {
                c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1.0e-3;
            c.CFLFraction = 0.3;
            c.Endtime = 0.25;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "Shock tube";
            if (true1D) {
                c.SessionName = String.Format("Shock tube, 1D, dgDegree = {0}, noOfCellsX = {1}, sensorLimit = {2:0.00E-00}", dgDegree, numOfCellsX, sensorLimit);
            } else {
                c.SessionName = String.Format("Shock tube, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);
            }

            return c;
        }

        public static IBMControl IBMContactDiscontinuity(string dbPath = null, int dgDegree = 2, int numOfCellsX = 20, int numOfCellsY = 5, double sensorLimit = 1e-3, bool saveToDb = false) {

            IBMControl c = new IBMControl();

            dbPath = @"c:\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            c.DomainType = DomainTypes.StaticImmersedBoundary;

            // Adjust height of cut cells such that we obtain AVCFL_cutcell = 0.5 * AVCFL
            // Here, this only depends on h_min
            double width = (xMax - xMin) / numOfCellsX;
            double height = (yMax - yMin) / numOfCellsY;
            double heightCutCell = (-2.0 * width * height) / (2.0 * height - 2.0 * Math.Sqrt(2) * width - 2.0 * Math.Sqrt(2) * height);

            double levelSetPosition = 2 * height + (height - heightCutCell);

            c.LevelSetFunction = delegate (double[] X, double t) {
                double y = X[1];
                return y - levelSetPosition;
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";

            //c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
            //c.SurfaceHMF_ProjectNodesToLevelSet = false;
            //c.SurfaceHMF_RestrictNodes = true;
            //c.SurfaceHMF_UseGaussNodes = false;
            //c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            //c.VolumeHMF_RestrictNodes = true;
            //c.VolumeHMF_UseGaussNodes = false;
            //c.LevelSetQuadratureOrder = 6;

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.2;
            c.AddVariable(IBMVariables.LevelSet, 1);

            bool AV = true;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            Variable sensorVariable = Variables.Density;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);

            if (AV) {
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            // Runge-Kutta
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 1;
            //c.NumberOfSubGrids = 4;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            c.AddVariable(Variables.ShockSensor, 0);

            if (AV) {
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            c.AddVariable(Variables.CFLArtificialViscosity, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "SubsonicInlet");
                grid.EdgeTagNames.Add(2, "SubsonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {   // bottom
                        return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                        return 3;
                    } else if (Math.Abs(X[0]) < 1e-14) {    // left
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                        return 2;
                    } else {
                        throw new System.Exception("Problem with definition of boundary conditions");
                    }
                });

                return grid;
            };

            Func<double[], double, double> DistanceToLine = delegate (double[] X, double t) {
                // direction vector
                Vector2D p1 = new Vector2D(0.5, 0.0);
                Vector2D p2 = new Vector2D(0.5, 1.0);
                Vector2D p = p2 - p1;

                // normal vector
                Vector2D n = new Vector2D(p.y, -p.x);
                n.Normalize();

                // angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = Math.PI / 2;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * (X[0]) + n.y * (X[1]);

                // shock speed
                double vs = 1;

                // distance to line
                double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);

                return distance;
            };

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double densityLeft = 100.0;
            double densityRight = 1.0;
            double pressure = 1.0;
            double velocityXLeft = 2.0;
            double velocityY = 0.0;

            c.AddBoundaryCondition("SubsonicInlet", Variables.Density, (X, t) => densityLeft);
            c.AddBoundaryCondition("SubsonicInlet", Variables.Velocity.xComponent, (X, t) => velocityXLeft);
            c.AddBoundaryCondition("SubsonicInlet", Variables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryCondition("SubsonicOutlet", Variables.Pressure, (X, t) => pressure);
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToLine(X, 0)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressure);
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => velocityXLeft);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => velocityY);

            // Time config 
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1.0e-4;
            c.CFLFraction = 0.3;
            //c.Endtime = 0.05;
            //c.NoOfTimesteps = int.MaxValue;
            c.NoOfTimesteps = 500;

            c.ProjectName = "IBM contact discontinuity";
            c.SessionName = String.Format("IBM contact discontinuity, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, CFLFraction = {3:0.00E-00}, ALTS {4}/{5}", dgDegree, numOfCellsX, numOfCellsY, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);

            return c;
        }

        public static CNSControl ContactDiscontinuity(double levelSetPosition = 0.45, string dbPath = null, int dgDegree = 2, int numOfCellsX = 20, int numOfCellsY = 5, double sensorLimit = 1, bool saveToDb = false) {

            CNSControl c = new CNSControl();

            dbPath = @"c:\bosss_db";
            //dbPath = @"\\fdyprime\userspace\geisenhofer\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            bool AV = true;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            Variable sensorVariable = Variables.Density;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);

            if (AV) {
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            // Runge-Kutta
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 1;
            //c.NumberOfSubGrids = 4;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            c.AddVariable(Variables.ShockSensor, 0);

            if (AV) {
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            c.AddVariable(Variables.CFL, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "SubsonicInlet");
                grid.EdgeTagNames.Add(2, "SubsonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {   // bottom
                        return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                        return 3;
                    } else if (Math.Abs(X[0]) < 1e-14) {    // left
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                        return 2;
                    } else {
                        throw new System.Exception("Problem with definition of boundary conditions");
                    }
                });

                return grid;
            };

            Func<double[], double, double> DistanceToLine = delegate (double[] X, double t) {
                // direction vector
                Vector2D p1 = new Vector2D(0.5, 0.0);
                Vector2D p2 = new Vector2D(0.5, 1.0);
                Vector2D p = p2 - p1;

                // normal vector
                Vector2D n = new Vector2D(p.y, -p.x);
                n.Normalize();

                // angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = Math.PI / 2;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * (X[0]) + n.y * (X[1]);

                // shock speed
                double vs = 1;

                // distance to line
                double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);

                return distance;
            };

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double densityLeft = 100.0;
            double densityRight = 1.0;
            double pressure = 1.0;
            double velocityXLeft = 2.0;
            double velocityY = 0.0;

            c.AddBoundaryCondition("SubsonicInlet", Variables.Density, (X, t) => densityLeft);
            c.AddBoundaryCondition("SubsonicInlet", Variables.Velocity.xComponent, (X, t) => velocityXLeft);
            c.AddBoundaryCondition("SubsonicInlet", Variables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryCondition("SubsonicOutlet", Variables.Pressure, (X, t) => pressure);
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToLine(X, 0)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressure);
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => velocityXLeft);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => velocityY);

            // Time config 
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.dtFixed = 1.0e-5;
            //c.CFLFraction = 0.3;
            //c.Endtime = 0.05;
            //c.NoOfTimesteps = int.MaxValue;
            c.NoOfTimesteps = 500;

            c.ProjectName = "Contact discontinuity";
            c.SessionName = String.Format("Contact discontinuity, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, CFLFraction = {3:0.00E-00}, ALTS {4}/{5}", dgDegree, numOfCellsX, numOfCellsY, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);

            return c;
        }

        public static CNSControl IBMObliqueShockTube(string dbPath = null, int dgDegree = 2, int numOfCellsX = 50, int numOfCellsY = 50, double sensorLimit = 1e-3, bool saveToDb = true) {
            IBMControl c = new IBMControl();

            dbPath = @"c:\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = 100;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1.3;
            double yMin = 0;
            double yMax = 0.65;
            double shockPosition = 0.8;

            bool AV = true;

            c.DomainType = DomainTypes.StaticImmersedBoundary;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            Variable sensorVariable = Variables.Density;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            c.AddVariable(Variables.ShockSensor, 0);

            if (AV) {
                double epsilon0 = 1.0;
                double kappa = 0.5;
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            }

            // Level set
            double beta = Math.PI / 16;
            double[] startOfLine1 = new double[] { 0.0, 0.05 };
            double distanceBetweenLines = 0.3;
            double[] startOfLine2 = new double[] { 0.0, startOfLine1[1] + distanceBetweenLines * Math.Cos(beta) };

            Func<double, double> line1 = delegate (double x) {
                return Math.Tan(beta) * (x - startOfLine1[0]) + startOfLine1[1];
            };
            Func<double, double> line2 = delegate (double x) {
                return Math.Tan(beta) * (x - startOfLine2[0]) + startOfLine2[1];
            };

            c.LevelSetFunction = delegate (double[] X, double t) {
                return -(X[1] - line1(X[0])) * (X[1] - line2(X[0]));
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.9;
            c.SaveAgglomerationPairs = true;
            c.AddVariable(IBMVariables.LevelSet, 2);

            // Runge-Kutta schemes
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 3;
            //c.NumberOfSubGrids = 3;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(Variables.ArtificialViscosity, 2);
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });

                return grid;
            };

            c.AddBoundaryCondition("AdiabaticSlipWall");

            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector2D p1 = new Vector2D(0.5, line1(0.5));
            Vector2D p2 = new Vector2D(0.6, line1(0.6));
            Vector2D normalVector = p2 - p1;

            // Direction vector of initial shock
            Vector2D r = new Vector2D(normalVector.y, -normalVector.x);
            r.Normalize();

            // Distance from a point X to the initial shock
            double[] p = new double[] { shockPosition, line1(shockPosition) };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            // Initial conditions
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;

            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.6;
            c.Endtime = 0.25;
            c.NoOfTimesteps = int.MaxValue;
            //c.NoOfTimesteps = 2;

            c.ProjectName = "IBM oblique shock tube";
            c.SessionName = String.Format("IBM oblique shock tube, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);

            return c;
        }

        public static CNSControl DoubleMachReflection(string dbPath = null, int dgDegree = 2, int numOfCellsX = 400, int numOfCellsY = 100, double xMax = 4, double sensorLimit = 1e-3) {
            CNSControl c = new CNSControl();

            //dbPath = @"e:\bosss_db\GridOfTomorrow";
            string dbPath2 = @"\\dc1\userspace\stange\HiWi_database\tests";
            //string dbPath2 = @"/work/scratch/ws35kire/work_db";
            //string dbPath2 = @"/home/ws35kire/test_db";
            c.DbPath = dbPath2;
            c.savetodb = true; //dbPath != null;
            c.saveperiod = 100;
            c.PrintInterval = 10;

            //c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            //c.ExplicitOrder = 1;
            c.ExplicitScheme = ExplicitSchemes.LTS;
            c.ExplicitOrder = 1;
            c.NumberOfSubGrids = 3;
            c.ReclusteringInterval = 5;
            c.FluxCorrection = false;

            // Add one balance constraint for each subgrid
            c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            c.DynamicLoadBalancing_Period = 5;

            bool AV = true;

            c.GridPartType = GridPartType.ParMETIS;
            //c.GridPartType = GridPartType.none;

            double xMin = 0;
            //double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 1.0;
            double lambdaMax = 20;

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax);
            }


            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Viscosity, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            c.AddVariable(Variables.Schlieren, dgDegree - 1);
            if (AV) {
                c.AddVariable(Variables.ShockSensor, dgDegree);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            // LTS variables
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {// unten
                        if (X[0] < xWall) {// unten (links)
                            return 1;
                        } else {// unten (Rest)
                            return 3;
                        }
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {// oben
                        return 1;
                    } else if (Math.Abs(X[0]) < 1e-14) { // links
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {// rechts
                        return 2;
                    } else {
                        throw new System.Exception("bla");
                    }
                });

                return grid;
            };

            Func<double[], double, double> DistanceToLine = delegate (double[] X, double t) {
                // direction vector
                Vector2D p1 = new Vector2D(xWall, 0.0);
                Vector2D p2 = new Vector2D(xWall + 1 / Math.Tan(Math.PI / 3), 1.0);
                Vector2D p = p2 - p1;

                // normal vector
                Vector2D n = new Vector2D(p.y, -p.x);
                n.Normalize();

                // angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = Math.PI / 3;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * (X[0]) + n.y * (X[1]);

                // shock speed
                double vs = 10;

                // distance to line
                double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);

                return distance;
            };

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            // Boundary conditions
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - Jump(X[0] - (0.1 + (X[1] + 20 * t) / 1.732)) * (8.0 - 1.4));
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 7.14471 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (7.14471 - 0.0));
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => -4.125 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (-4.125 - 0.0));
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (116.5 - 1.0));

            c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(DistanceToLine(X, t)) * (8.0 - 1.4));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToLine(X, t)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToLine(X, t)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - SmoothJump(DistanceToLine(X, t)) * (116.5 - 1.0));

            c.AddBoundaryCondition("SupersonicOutlet", Variables.Pressure, (X, t) => 1.0);
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            //c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (8.0 - 1.4));
            //c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => 57.157 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (57.157 - 0.0));
            //c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => -33.0 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (-33 - 0.0));
            //c.InitialValues_Evaluators.Add(Variables.Energy, X => 563.544 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (563.544 - 2.5));

            c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(DistanceToLine(X, 0)) * (8.0 - 1.4));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToLine(X, 0)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToLine(X, 0)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 116.5 - SmoothJump(DistanceToLine(X, 0)) * (116.5 - 1.0));

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            //c.dtFixed = 1.0e-6;
            //c.CFLFraction = 0.5; // altes Setting fuer Rechnungen auf Lichtenberg
            c.CFLFraction = 0.05;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "Double Mach reflection";
            c.SessionName = String.Format("DMR_withDLB, dgDegree = {0}, numOfCellsX = {1}, numOfCellsY = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}, lamdaMax = {7}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, lambdaMax);
            //c.Tags.Add("Double Mach reflection");
            //c.Tags.Add("Artificial viscosity");
            //c.Tags.Add("Adaptive local time stepping");

            return c;
        }

        public static IBMControl IBMDoubleMachReflection(string dbPath = null, int dgDegree = 2, int numOfCellsX = 50, int numOfCellsY = 50, double sensorLimit = 1e-3) {
            IBMControl c = new IBMControl();

            //dbPath = @"c:\bosss_db";
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db";
            //dbPath = @"/work/scratch/yp19ysog/bosss_db_lb_scratch";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            c.DomainType = DomainTypes.StaticImmersedBoundary;

            c.DynamicLoadBalancing_CellClassifier = new IBMCellClassifier();
            c.DynamicLoadBalancing_Period = 0;
            c.DynamicLoadBalancing_RedistributeAtStartup = true;
            c.DynamicLoadBalancing_CellCostEstimatorFactories.Add(IBMCellCostEstimator.GetStaticCostBasedEstimator());
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(IBMCellCostEstimator.GetMultiBalanceConstraintedBasedEstimators());

            double xMin = 0.0;
            double xMax = 0.5;
            double yMin = 0.0;
            double yMax = 0.5;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            // Level set
            double angleInDegree = 30;
            double beta = 2 * Math.PI / 360 * angleInDegree;   // the wall has an angle of 60 degree
            double[] startOfRamp = new double[] { xWall, 0.0 };

            Func<double, double> ramp = delegate (double x) {
                return Math.Tan(beta) * (x - startOfRamp[0]) + startOfRamp[1];
            };

            c.LevelSetFunction = delegate (double[] X, double t) {
                return X[1] - ramp(X[0]);
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.9;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, 2);

            bool AV = true;

            // Runge-Kutta
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 1;
            //c.NumberOfSubGrids = 3;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

            c.GridPartType = GridPartType.ParMETIS;
            //c.GridPartType = GridPartType.none;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            if (dgDegree >= 1) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            }

            double lambdaMax = 20;
            if (AV) {
                double epsilon0 = 1.0;
                double kappa = 1.0;
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax);
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Viscosity, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            if (dgDegree >= 1) {
                c.AddVariable(Variables.Schlieren, dgDegree - 1);
            }

            if (AV) {
                c.AddVariable(Variables.ShockSensor, dgDegree);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            // LTS variables
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {// unten
                        if (X[0] < xWall) {// unten (links)
                            return 1;
                        } else {// unten (Rest)
                            return 3;
                        }
                        //return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {// oben
                        return 1;
                    } else if (Math.Abs(X[0]) < 1e-14) { // links
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {// rechts
                        return 2;
                    } else {
                        throw new System.Exception("bla");
                    }
                });

                return grid;
            };

            // Direction vector of initial shock (vertical)
            Vector2D r = new Vector2D(0.0, 1.0);

            // Current x-position of the shock
            double shockSpeed = 10;
            Func<double, double> getShockXPosition = delegate (double time) {
                return xWall + shockSpeed * time;
            };


            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            // Boundary conditions
            c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(X[0] - getShockXPosition(t)) * (8.0 - 1.4));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 8.25 - SmoothJump(X[0] - getShockXPosition(t)) * (8.25 - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - SmoothJump(X[0] - getShockXPosition(t)) * (116.5 - 1.0));
            c.AddBoundaryCondition("SupersonicOutlet");
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(X[0] - getShockXPosition(0)) * (8.0 - 1.4));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 8.25 - SmoothJump(X[0] - getShockXPosition(0)) * (8.25 - 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 116.5 - SmoothJump(X[0] - getShockXPosition(0)) * (116.5 - 1.0));

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.05;
            //c.dtFixed = 1.0e-6;
            c.CFLFraction = 0.3;
            c.NoOfTimesteps = 2;

            c.ProjectName = "IBM double Mach reflection";
            c.SessionName = String.Format("IBM DMR LOAD BAL, dgDegree = {0}, numOfCellsX = {1}, numOfCellsY = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}, lamdaMax = {7}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, lambdaMax);

            return c;
        }

        public static CNSControl DoubleMachReflectionStudy(int gridFactor = 1, int dgDegree = 0, double cfl = 0.5) {
            //string dbPath = @"\\dc1\userspace\stange\HiWi_database\tests";
            string dbPath = @"/work/scratch/ws35kire/work_db/";
            // int gridFactor = 1; // Number of refinements, 1 equals 400x100 cells
            //int dgDegree = 0;
            //double cfl = 0.5;

            CNSControl c = new CNSControl();
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            int numOfCellsX = 400 * gridFactor;
            int numOfCellsY = 100 * gridFactor;
            int saves = Convert.ToInt32(50 * gridFactor * 0.5 / cfl);

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 50;
            c.PrintInterval = 50;

            c.GridPartType = GridPartType.ParMETIS;

            double xMin = 0;
            double xMax = 4;
            double yMin = 0;
            double yMax = 1;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            c.AddVariable(Variables.Schlieren, dgDegree);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {// unten
                        if (X[0] < xWall) {// unten (links)
                            return 1;
                        } else {// unten (Rest)
                            return 3;
                        }
                        //return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {// oben
                        return 1;
                    } else if (Math.Abs(X[0]) < 1e-14) { // links
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {// rechts
                        return 2;
                    } else {
                        throw new System.Exception("bla");
                    }
                });

                return grid;
            };

            Func<double[], double, double> DistanceToLine = delegate (double[] X, double t) {
                // direction vector
                Vector2D p1 = new Vector2D(xWall, 0.0);
                Vector2D p2 = new Vector2D(xWall + 1 / Math.Tan(Math.PI / 3), 1.0);
                Vector2D p = p2 - p1;

                // normal vector
                Vector2D n = new Vector2D(p.y, -p.x);
                n.Normalize();

                // angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = Math.PI / 3;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * (X[0]) + n.y * (X[1]);

                // shock speed
                double vs = 10;

                // distance to line
                double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);

                return distance;
            };

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(DistanceToLine(X, t)) * (8.0 - 1.4));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToLine(X, t)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToLine(X, t)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - SmoothJump(DistanceToLine(X, t)) * (116.5 - 1.0));

            c.AddBoundaryCondition("SupersonicOutlet", Variables.Pressure, (X, t) => 1.0);
            c.AddBoundaryCondition("AdiabaticSlipWall");

            c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(DistanceToLine(X, 0)) * (8.0 - 1.4));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToLine(X, 0)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToLine(X, 0)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 116.5 - SmoothJump(DistanceToLine(X, 0)) * (116.5 - 1.0));

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = cfl;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "Double Mach reflection 0th";
            c.SessionName = String.Format("DMR, dgDegree = {0}, numOfCellsX = {1}, numOfCellsY = {2}, CFLFraction = {3:0.00E-00}", dgDegree, numOfCellsX, numOfCellsY, c.CFLFraction);
            c.Tags.Add("Double Mach reflection");
            c.Tags.Add("Runge-Kutta");

            return c;
        }

        public static IBMControl MovingIBMIsentropicVortex(string dbPath = null, int dgDegree = 3, int noOfCellsPerDirection = 20, double initialLevelSetPosition = -0.9, double agglomerationThreshold = 0.2) {
            IBMControl c = new IBMControl();

            double advectionVelocity = 1.0;

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1;

            c.ProjectName = "Moving IBM Isentropic vortex";
            c.Tags.Add("Isentropic vortex");
            c.Tags.Add("Moving IBM");

            c.DomainType = DomainTypes.MovingImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "adiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.SurfaceHMF_ProjectNodesToLevelSet = false;
            c.SurfaceHMF_RestrictNodes = true;
            c.SurfaceHMF_UseGaussNodes = false;
            c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            c.VolumeHMF_RestrictNodes = true;
            c.VolumeHMF_UseGaussNodes = false;
            c.LevelSetQuadratureOrder = 10;
            c.LevelSetBoundaryTag = "supersonicInlet";

            c.AgglomerationThreshold = agglomerationThreshold;
            c.SaveAgglomerationPairs = true;

            double gamma = c.EquationOfState.HeatCapacityRatio;
            Func<double[], double, double> x = (X, t) => X[0] - advectionVelocity * t;
            Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
            Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
            Func<double[], double, double> rho = (X, t) => Math.Pow(
                1.0 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)),
                1.0 / (gamma - 1.0));
            Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma);
            Func<double[], double, double> uAbs = (X, t) => r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
            Func<double[], double, double> u = (X, t) => advectionVelocity - Math.Sin(phi(X, t)) * uAbs(X, t);
            Func<double[], double, double> v = (X, t) => Math.Cos(phi(X, t)) * uAbs(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));

            c.LevelSetFunction = delegate (double[] X, double time) {
                double amplitude = 0.3;
                double newLevelSetPosition = initialLevelSetPosition + amplitude * Math.Sin(10.0 * time);
                return X[1] - newLevelSetPosition;
            };

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], u);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], v);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, p);

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = 10;
            c.NoOfTimesteps = 1000;

            return c;
        }

        public static IBMControl MovingIBMCylinder(string dbPath = @"e:\bosss_db\GridOfTomorrow\", double Mach = 0.2, int refinements = 0, int dgDegree = 2, double agglomerationThreshold = 0.3, int levelSetQuadratureOrder = 10) {
            IBMControl c = new IBMControl();
            c.DbPath = dbPath;
            c.savetodb = dbPath != null;

            c.ProjectName = "Moving IBM cylinder";
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
            c.AddVariable(IBMVariables.LevelSet, 2);

            //c.GridFunc = delegate {
            //    double totalHeight = 40.0;
            //    int lengthToHeightRatio = 5;
            //    int heightToRefinedHeightRatio = 4;
            //    int refinementFactor = 4;

            //    int baseNumberOfCells = 4 << refinements;
            //    double refinedHeigth = totalHeight / heightToRefinedHeightRatio;

            //    double[] lowerLeft = new double[] { -2.0, 0.0 };
            //    double[] upperRight = new double[] { lengthToHeightRatio * totalHeight, totalHeight };
            //    double[] upperRightRefined = new double[] { lengthToHeightRatio * totalHeight, refinedHeigth };

            //    int noOfCellsX = lengthToHeightRatio * baseNumberOfCells;
            //    GridCommons grid = Grid2D.HangingNodes2D(
            //        false,
            //        false,
            //        new GridBox(lowerLeft, upperRight,        noOfCellsX,                    baseNumberOfCells),
            //        new GridBox(lowerLeft, upperRightRefined, refinementFactor * noOfCellsX, refinementFactor * baseNumberOfCells / refinementFactor));
            //    grid.EdgeTagNames.Add(1, "supersonicInlet");
            //    grid.EdgeTagNames.Add(2, "adiabaticSlipWall");
            //    Func<double[], byte> func = x => 1;
            //    grid.DefineEdgeTags(func);

            //    return grid;
            //};
            c.GridGuid = new Guid("f80f004a-cfe2-44cb-bf26-6b3837c1dcb2");

            c.GridPartType = GridPartType.ParMETIS;
            c.GridPartOptions = "10";

            double gamma = c.EquationOfState.HeatCapacityRatio;
            c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity[0], X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity[1], X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 1.0);

            double cylinderVelocity = 1.0;
            c.LevelSetFunction = delegate (double[] X, double t) {
                double x = X[0] - cylinderVelocity * t;
                double y = X[1];
                return x * x + y * y - 1.0;
            };

            c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 1.0);

            c.AddBoundaryCondition("adiabaticSlipWall", Variables.Velocity[0], (X, t) => cylinderVelocity);
            c.AddBoundaryCondition("adiabaticSlipWall", Variables.Velocity[1], (X, t) => 0.0);
            c.LevelSetBoundaryTag = "adiabaticSlipWall";

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
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
            c.CFLFraction = 0.1;
            c.Endtime = 200.0;
            c.NoOfTimesteps = int.MaxValue;

            c.saveperiod = 100;
            c.PrintInterval = 10;

            return c;
        }

        public static CNSControl CylindricalExplosion(string dbPath = @"e:\bosss_db\cylindricalExplosion\", int dgDegree = 0, int numCells = 20, double sensorLimit = 1e-2) {
            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 25;
            c.PrintInterval = 1;
            c.GridPartType = GridPartType.none;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.ProjectName = "Cylindrical explosion";
            c.SessionName = String.Format(
                "Cylindrical explosion, degree={0}, numCells={1}, sensorLimit = {2:0.00E-00}", dgDegree, numCells, sensorLimit);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.CFL, 0);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-1.0, 1.0, numCells + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicOutlet");
                grid.DefineEdgeTags(X => 1);

                return grid;
            };

            double cellSize = 2.0 / numCells;
            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 1.5 * cellSize / Math.Max(dgDegree, 1);
                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            c.AddBoundaryCondition("SupersonicOutlet");

            double radius = 0.4;
            Func<double[], double> levelSet = X => Math.Sqrt(X[0] * X[0] + X[1] * X[1]) - radius;

            double rhoIn = 1.0;
            double rhoOut = 0.125;
            double pIn = 1.0;
            double pOut = 0.1;
            c.InitialValues_Evaluators.Add(Variables.Density, X => rhoIn - SmoothJump(levelSet(X)) * (rhoIn - rhoOut));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pIn - SmoothJump(levelSet(X)) * (pIn - pOut));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            // Shock-capturing
            {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ActiveOperators |= Operators.ArtificialViscosity;

                double epsilon0 = 1.0;
                double kappa = 1.0;

                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(
                    c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);

                c.AddVariable(Variables.ShockSensor, dgDegree);
                c.AddVariable(Variables.ArtificialViscosity, 2);

                c.Tags.Add("Artificial viscosity");

                c.PrandtlNumber = 0.71;
                c.ReynoldsNumber = 1.0;
            }

            return c;
        }

        public static CNSControl CylindricalExplosion1D(string dbPath = @"e:\bosss_db\GridOfTomorrow", int numCells = 1000) {
            int dgDegree = 0;

            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1000;
            c.PrintInterval = 10;
            c.GridPartType = GridPartType.none;

            c.ActiveOperators = Operators.Convection | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.CustomContinuitySources.Add(
                map => new AdHocSourceTerm(map, (X, t, state) =>
                    1.0 / X[0] * state.Momentum[0]));
            c.CustomMomentumSources[0].Add(
                map => new AdHocSourceTerm(map, (X, t, state) =>
                    1.0 / X[0] * state.Momentum[0] * state.Velocity[0]));
            c.CustomEnergySources.Add(
                map => new AdHocSourceTerm(map, (X, t, state) =>
                    1.0 / X[0] * state.Velocity[0] * (state.Energy + state.Pressure)));

            c.ProjectName = "Cylindrical explosion 1D";
            c.SessionName = String.Format(
                "Cylindrical explosion 1D, numCells={1}", dgDegree, numCells);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.CFL, 0);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(0.0, 1.0, numCells + 1);
                var grid = Grid1D.LineGrid(nodes);

                grid.EdgeTagNames.Add(1, "SupersonicOutlet");
                grid.DefineEdgeTags(X => 1);

                return grid;
            };

            c.AddBoundaryCondition("SupersonicOutlet");

            double R = 0.4;

            double rhoIn = 1.0;
            double rhoOut = 0.125;
            double pIn = 1.0;
            double pOut = 0.1;

            Func<double, double> Jump = (x => x < 0.0 ? 0.0 : 1.0);
            c.InitialValues_Evaluators.Add(Variables.Density, X => rhoIn - Jump(X[0] - R) * (rhoIn - rhoOut));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pIn - Jump(X[0] - R) * (pIn - pOut));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }

        public static CNSControl SphericalExplosion(string dbPath = @"e:\bosss_db\GridOfTomorrow", int dgDegree = 0, int numCells = 10, double sensorLimit = 1e-2) {
            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 25;
            c.PrintInterval = 1;
            c.GridPartType = GridPartType.none;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.ProjectName = "Spherical explosion";
            c.SessionName = String.Format(
                "Spherical explosion, degree={0}, numCells={1}, sensorLimit = {2:0.00E-00}", dgDegree, numCells, sensorLimit);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.CFL, 0);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-1.0, 1.0, numCells + 1);
                var grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes);

                grid.EdgeTagNames.Add(1, "SupersonicOutlet");
                grid.DefineEdgeTags(X => 1);

                return grid;
            };

            double cellSize = 2.0 / numCells;
            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.5 * cellSize / Math.Max(dgDegree, 1);
                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            c.AddBoundaryCondition("SupersonicOutlet");

            double radius = 0.4;
            Func<double[], double> levelSet = X => Math.Sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]) - radius;

            double rhoIn = 1.0;
            double rhoOut = 0.125;
            double pIn = 1.0;
            double pOut = 0.1;
            c.InitialValues_Evaluators.Add(Variables.Density, X => rhoIn - SmoothJump(levelSet(X)) * (rhoIn - rhoOut));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.zComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pIn - SmoothJump(levelSet(X)) * (pIn - pOut));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            // Shock-capturing
            {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ActiveOperators |= Operators.ArtificialViscosity;

                double epsilon0 = 1.0;
                double kappa = 1.0;
                //double lambdaMax = 2.0;    // determined manually in Visit
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(
                    c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);

                c.AddVariable(Variables.ShockSensor, dgDegree);
                c.AddVariable(Variables.ArtificialViscosity, 3);

                c.Tags.Add("Artificial viscosity");

                c.PrandtlNumber = 0.71;
                c.ReynoldsNumber = 1.0;
            }

            return c;
        }


        /// <summary>
        /// Taylor Green Vortex testcase, cf. 1st International Workshop on High-Order CFD Methods
        /// http://dept.ku.edu/~cfdku/hiocfd.html
        /// Testcase: C3.5
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="dgDegree"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <returns></returns>
        public static CNSControl TaylorGreenVortex(int dgDegree, int noOfCellsPerDirection = 20, string dbPath = @"c:\bosss_dbv2\HHLR\TGV_3D\") {
            CNSControl c = new CNSControl();

            //@"/home/kraemer/CNS/dbe"
            //@"/work/scratch/kraemer/TGV_3D/"
            //@"c:\bosss_dbv2\HHLR\TGV_3D\"

            c.DbPath = dbPath;
            c.savetodb = true;
            c.saveperiod = 100;
            c.PrintInterval = 1;

            c.ProjectName = "Taylor Green Vortex";
            c.Tags.Add("Re=1600");
            c.Tags.Add("Mach=0.1");

            c.ActiveOperators = Operators.Convection | Operators.Diffusion;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
            c.SIPGPenaltyScaling = 1.0;

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;
            c.NumberOfSubGrids = 3;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 0.1;
            c.ReynoldsNumber = 1600;
            c.PrandtlNumber = 0.71;
            c.ViscosityLaw = new ConstantViscosity();

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Momentum.zComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Velocity.zComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.KineticEnergy, dgDegree);


            //c.GridFunc = delegate () {
            //    double[] nodes = GenericBlas.Linspace(-Math.PI, Math.PI, noOfCellsPerDirection + 1);
            //    GridCommons grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes, true, true, true);
            //    return grid;
            //};

            String grid = "";

            switch (noOfCellsPerDirection) {
                case 16:
                    grid = "fad59f19-910b-4bad-ae28-b4c959f4b179";
                    break;
                case 32:
                    grid = "2db1c6de-3aa2-447f-b0ce-092dda48f0c5";
                    break;
                case 64:
                    grid = "36271cdf-94e3-45a9-a8c5-a305a50cf5a1";
                    break;
            }

            Guid gridGuid;
            if (Guid.TryParse(grid, out gridGuid)) {
                c.GridGuid = gridGuid;
            } else {
                throw new Exception(
                 "Could not find a grid at " + grid);
            }




            Func<double[], double, double> u = (X, t) => Math.Sin(X[0]) * Math.Cos(X[1]) * Math.Cos(X[2]);
            Func<double[], double, double> v = (X, t) => -Math.Cos(X[0]) * Math.Sin(X[1]) * Math.Cos(X[2]);
            Func<double[], double, double> p = (X, t) => 1.0 + 1.0 / 16.0 * (Math.Cos(2 * X[0]) * Math.Cos(2 * X[1])) * (Math.Cos(2 * X[2]) + 2);
            Func<double[], double, double> rho = (X, t) => p(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.zComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));

            //c.Queries.Add("IntegralKineticEnergy", QueryLibrary.Integral(Variables.KineticEnergy.Name));

            //c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;

            c.dtMin = 1.0e-3;
            c.dtMax = 1.0e-3;
            c.CFLFraction = 0.5;
            c.Endtime = 20.0;
            c.NoOfTimesteps = int.MaxValue;
            c.ResidualInterval = 10;
            c.PrintInterval = 10;

            return c;
        }
    }
}
