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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.Source;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CNS {

    /// <summary>
    /// A set of exemplary control files for the subsonic CNS IBM solver. You can try the individual
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
    public static class ControlExamples_Subsonic {

        /// <summary>
        /// Isentropic vortex in a fully periodic domain.
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <param name="dgDegree"></param>
        /// <param name="advectionVelocity"></param>
        /// <returns></returns>
        public static CNSControl IsentropicVortex(string dbPath = null, int noOfCellsPerDirection = 20, int dgDegree = 2, double advectionVelocity = 0.1) {
            CNSControl c = new CNSControl();

            c.dtMin = 0.0;
            c.dtMax = 5.0e-1;
            //c.CFLFraction = 0.9;
            c.dtFixed = 1e-2;
            c.Endtime = 2.0;
            c.NoOfTimesteps = int.MaxValue;

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1000;
            c.PrintInterval = 1;
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
            c.ExplicitOrder = 1;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => p(X, 0.0));


            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(CompressibleVariables.Density, rho));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(CNSVariables.Pressure, p));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(CNSVariables.Velocity.xComponent, u));
            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(CNSVariables.Entropy, (X, t) => 1.0));

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
                    controls[ii].Paramstudy_CaseIdentification.AddRange(new Tuple<string, object>[] {
                    new Tuple<string, object>("dgDegree", i),
                    new Tuple<string, object>("noOfCellsPerDirection", noOfCellsPerDirection)
                    });
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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            double extent = 10.0;
            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-extent, extent, noOfCellsPerDirection + 1);
                Grid2D grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                //grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.EdgeTagNames.Add(1, "subsonicOutlet");
                grid.DefineEdgeTags((Vector X) => 1);
                return grid;
            };

            StateVector refState = StateVector.FromPrimitiveQuantities(
                c.GetMaterial(),
                1.0,
                new Vector(advectionVelocity, 0.0, 0.0),
                8.0);

            Func<double[], double> pulse =
                X => (1 + 0.01 * Math.Exp(-(X[0] * X[0] + X[1] * X[1]) / 2 / 0.8 / 0.8));
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => refState.Density * pulse(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => refState.Velocity.x);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => refState.Velocity.y);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => refState.Pressure * pulse(X));

            c.AddBoundaryValue("subsonicOutlet", CNSVariables.Pressure, (X, t) => refState.Pressure);

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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Temperature, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

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

                    Func<Vector, string> etf = delegate (Vector X) {
                        throw new NotImplementedException();
                    };

                    _grid.DefineEdgeTags(etf);

                    return _grid;
                };
            }

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum[0], X => 1.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum[1], X => 0.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Energy, X => 45.14285714);

            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, (X, t) => 17.85714286);

            c.AddBoundaryValue("subsonicInlet", CompressibleVariables.Density, (X, t) => 1.0);
            c.AddBoundaryValue("subsonicInlet", CNSVariables.Velocity[0], (X, t) => 1.0);
            c.AddBoundaryValue("subsonicInlet", CNSVariables.Velocity[1], (X, t) => 0.0);

            c.AddBoundaryValue("subsonicOutlet", CNSVariables.Pressure, (X, t) => 17.85714286);

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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "adiabaticSlipWall");
                grid.DefineEdgeTags((Vector X) => 1);
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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => p(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - levelSetPosition;

            c.AddBoundaryValue("adiabaticSlipWall");
            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, rho);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], u);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], v);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, p);

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(CompressibleVariables.Density, rho));
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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
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
                Func<Vector, byte> func = delegate (Vector x) {
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
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[0], X => Mach * Math.Sqrt(gamma));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[1], X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 1.0);

            c.LevelSetFunction = (X, t) => X[0] * X[0] + X[1] * X[1] - 1.0 * 1.0;

            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], (X, t) => Mach * Math.Sqrt(gamma));
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, (X, t) => 1.0);

            c.AddBoundaryValue("adiabaticSlipWall");
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
                    c.Paramstudy_CaseIdentification.AddRange(new Tuple<string, object>[] {
                        new Tuple<string, object>("dgDegree", dgDegree),
                        new Tuple<string, object>("refinements", refinements),
                    });

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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
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
                Func<Vector, byte> func = delegate (Vector x) {
                    return 1;
                };
                grid.DefineEdgeTags(func);

                return grid;
            };

            c.GridPartType = GridPartType.ParMETIS;
            c.GridPartOptions = "5";

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum[0], X => 1.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum[1], X => 0.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Energy, X => 45.14285714);

            c.LevelSetFunction = (X, t) => X[0] * X[0] + X[1] * X[1] - 1.0 * 1.0;

            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, (X, t) => 17.85714286);

            c.AddBoundaryValue("adiabaticWall");
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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
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
                Func<Vector, byte> func = delegate (Vector x) {
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
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[0], X => Mach * Math.Sqrt(gamma));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[1], X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 1.0);
            c.InitialValues_Evaluators.Add(IBMVariables.LevelSet, X => X[0] * X[0] + X[1] * X[1] - 1.0 * 1.0);

            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], (X, t) => Mach * Math.Sqrt(gamma));
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, (X, t) => 1.0);

            c.AddBoundaryValue("adiabaticSlipWall");
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

            // Primary CNSVariables
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);


            // Grid
            c.GridFunc = delegate {
                double[] nodesX = GenericBlas.Linspace(-BoxSizeX / 2.0, BoxSizeX / 2.0, nodesInX + 1);
                double[] nodesY = GenericBlas.Linspace(-BoxSizeY / 2.0, BoxSizeY / 2.0, nodesInY + 1);
                var grid = Grid2D.Cartesian2DGrid(nodesX, nodesY, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.DefineEdgeTags((Vector X) => 1);
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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => p(X, 0.0));
            c.InitialValues_Evaluators.Add(IBMVariables.LevelSet, X => X[1] - levelSetPosition);

            c.AddBoundaryValue("adiabaticSlipWall");
            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, rho);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], u);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], v);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, p);
            c.LevelSetBoundaryTag = "supersonicInlet";

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(CompressibleVariables.Density, rho));
            c.Queries.Add("L2ErrorPressure", IBMQueries.L2Error(state => state.Pressure, p));
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "adiabaticSlipWall");
                grid.DefineEdgeTags((Vector X) => 1);
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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => p(X, 0.0));

            c.LevelSetFunction = delegate (double[] X, double time) {
                double amplitude = 0.3;
                double newLevelSetPosition = initialLevelSetPosition + amplitude * Math.Sin(10.0 * time);
                return X[1] - newLevelSetPosition;
            };

            c.AddBoundaryValue("adiabaticSlipWall");
            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, rho);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], u);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], v);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, p);

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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
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
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[0], X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[1], X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 1.0);

            double cylinderVelocity = 1.0;
            c.LevelSetFunction = delegate (double[] X, double t) {
                double x = X[0] - cylinderVelocity * t;
                double y = X[1];
                return x * x + y * y - 1.0;
            };

            c.AddBoundaryValue("supersonicInlet", CompressibleVariables.Density, (X, t) => 1.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], (X, t) => 0.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, (X, t) => 1.0);

            c.AddBoundaryValue("adiabaticSlipWall", CNSVariables.Velocity[0], (X, t) => cylinderVelocity);
            c.AddBoundaryValue("adiabaticSlipWall", CNSVariables.Velocity[1], (X, t) => 0.0);
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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.zComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.zComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.KineticEnergy, dgDegree);


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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.zComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => p(X, 0.0));

            //c.Queries.Add("IntegralKineticEnergy", QueryLibrary.Integral(CNSVariables.KineticEnergy.Name));

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

        public static IBMControl IBMGaussianBump(string dbPath = null, int savePeriod = 1000, int noOfCellsY = 16 * 4, int dgDegree = 2, int lsDegree = 8, double CFL = 0.3, double agg = 0.3, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 100, int maxNumOfSubSteps = 0, double epsilonX = 0.0, double epsilonY = 0.0) {
            IBMControl c = new IBMControl();

            // Session Settings
            //dbPath = @"c:\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = c.DbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 1000.0;
            c.CFLFraction = CFL;
            //c.dtFixed = 4.8e-3;
            c.NoOfTimesteps = 10;   // CNS Config

            // Residual logging
            c.ResidualInterval = 100;   // CNS Config
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            //c.ResidualBasedTerminationCriteria.Add("changeRate_abs_rhoE", 1E-3);
            c.ResidualBasedTerminationCriteria.Add("changeRate_abs_rhoE", 1E-8);     // CNS Config

            // Queries
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 2.8571428571428));

            // IBM Settings
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.LevelSetQuadratureOrder = 2 * dgDegree;
            c.AgglomerationThreshold = agg;

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            //c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            //c.SurfaceHMF_ProjectNodesToLevelSet = false;
            //c.SurfaceHMF_RestrictNodes = true;
            //c.SurfaceHMF_UseGaussNodes = false;
            //c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            //c.VolumeHMF_RestrictNodes = true;
            //c.VolumeHMF_UseGaussNodes = false;

            //c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("60688cbc-707d-4777-98e6-d237796ec14c"), -1);

            // Solver Type
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Time-Stepping Settings
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Material Settings
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            // Primary CNSVariables
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(IBMVariables.LevelSet, lsDegree);

            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            // Grid

            //switch (noOfCells) {
            //    case 8:
            //        c.GridGuid = new Guid("7337e273-542f-4b97-b592-895ac3422621");
            //        break;
            //    case 16:
            //        c.GridGuid = new Guid("32e5a779-2aef-4ea2-bdef-b158ae785f01");
            //        break;
            //    case 32:
            //        c.GridGuid = new Guid("e96c9f83-3486-4e45-aa3b-9a436445a059");
            //        break;
            //    case 64:
            //        c.GridGuid = new Guid("a86f1b67-4fa3-48ed-b6df-dcea370eb2c0");
            //        break;
            //    default:
            //        throw new ArgumentException("Wrong Grid Input");
            //}

            c.GridFunc = delegate {
                double xBoundary = 20.0;
                double yBoundary = 20.0;
                double yBottom = 0.0;

                double[] xnodes = GenericBlas.Linspace(-xBoundary, xBoundary, 2 * noOfCellsY + 1);

                //double ySplit = 6.0;
                //int ySplitNoOfCells = (int) (0.5*noOfCells);
                //double[] ynodes1 = GenericBlas.Linspace(yBottom, ySplit, ySplitNoOfCells + 1);
                //double[] ynodes2 = GenericBlas.Linspace(ySplit, yBoundary, noOfCells-ySplitNoOfCells + 1);
                //ynodes1 = ynodes1.GetSubVector(0, ynodes1.Length - 1);
                //double[] ynodes = ArrayTools.Cat(ynodes1, ynodes2);

                double[] ynodes = GenericBlas.Linspace(yBottom, yBoundary, noOfCellsY + 1);

                GridCommons grid = Grid2D.Cartesian2DGrid(xnodes, ynodes);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                //grid.EdgeTagNames.Add(2, "SubsonicInlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                //grid.EdgeTagNames.Add(4, "SubsonicOutlet");

                Func<double[], byte> func = delegate (double[] x) {

                    if (Math.Abs(x[0] + xBoundary) < 1e-5) { // Inflow
                        return 1;
                    } else if (Math.Abs(x[0] - xBoundary) < 1e-5) { // Outflow
                        return 1;
                    } else if (Math.Abs(x[1] - yBoundary) < 1e-5) { // Top
                        return 1;
                    } else { // Bottom
                        return 3;
                    }
                };
                grid.DefineEdgeTags(func);
                grid.Name = "IBM-[" + -xBoundary + "," + xBoundary + "]x[" + yBottom + "," + yBoundary + "]_Cells:(" + 2 * noOfCellsY + "x" + noOfCellsY + ")";
                return grid;
            };

            // Functions
            Func<double[], double, double> rho = (X, t) => 1.0;
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 0.0;
            // M_infty = 0.5, set u0 = 1.0, using u0 = M_infty * sqrt(gamma * p / rho) ==> p = rho * u0 * u0 / (M_infty^2 * gamma)
            Func<double[], double, double> pressure = (X, t) => 2.8571428571428;

            // Initial Values
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u0(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => u1(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - epsilonY - 0.01 - 0.3939 * Math.Exp(-0.5 * (X[0] - epsilonX) * (X[0] - epsilonX));

            // Supersonic boundary conditions
            c.AddBoundaryValue("AdiabaticSlipWall");
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, rho);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, u0);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, u1);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, pressure);

            // Subsonic boundary conditions
            //c.AddBoundaryValue("SubsonicInlet", CompressibleVariables.Density, rho);
            //c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.xComponent, u0);
            //c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.yComponent, u1);
            //c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, pressure);

            c.ProjectName = "IBMGaussianBump";
            c.SessionName = c.CutCellQuadratureType + "_(" + 2 * noOfCellsY + "x" + noOfCellsY + ")_CFL=" + c.CFLFraction + "_lsQuadOrder=" + c.LevelSetQuadratureOrder + "_p=" + dgDegree + "_agg=" + c.AgglomerationThreshold + "_epsX=" + epsilonX + "_epsY=" + epsilonY;

            return c;
        }

        public static IBMControl IBMGaussianBump_HHLR(int savePeriod = 100, int noOfCellsY = 16, int dgDegree = 2, int lsDegree = 8, double CFL = 0.3, double agg = 0.0) {
            // Lichtenberg
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_ibmgaussianbump";

            IBMControl c = IBMGaussianBump(dbPath, savePeriod, noOfCellsY, dgDegree, lsDegree, CFL, agg);

            c.ProjectName = string.Format("IBMGaussianBump_HHLR_agg{0}", agg);

            return c;
        }
    }
}
