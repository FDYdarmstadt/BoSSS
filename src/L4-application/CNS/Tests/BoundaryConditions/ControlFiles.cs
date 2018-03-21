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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.Source;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CNS.Tests.BoundaryConditions {

    /// <summary>
    /// A set of test cases to verify the correct implementation of different
    /// types of boundary conditions for the Euler equations. To separate the
    /// influences of different types of errors, most of them use
    /// <see cref="SupersonicInlet"/> at some edges, i.e. a non-physical
    /// Dirichlet boundary condition. This speeds up convergence and does not
    /// spoil the exact solution, as long as the Dirichlet values do not
    /// contradict the analytical solution
    /// </summary>
    public static class ControlFiles {

        private static CNSControl[] GetTemplates() {
            int minDivisions = 0;
            int maxDivisions = 2;
            int minDegree = 0;
            int maxDegree = 3;

            List<CNSControl> result = new List<CNSControl>();
            for (int dgDegree = minDegree; dgDegree <= maxDegree; dgDegree++) { 
                for (int divisions = minDivisions; divisions <= maxDivisions; divisions++) {

                    CNSControl c = new CNSControl();

                    //c.DbPath = @"c:\bosss_dbv2\exp";
                    c.savetodb = false;
                    int noOfCells = (2 << divisions) * 8;
                    
                    c.PrintInterval = 1000;
                    c.ResidualInterval = 100;
                    c.saveperiod = 100000;

                    c.ActiveOperators = Operators.Convection | Operators.CustomSource;
                    c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
                    c.ExplicitScheme = ExplicitSchemes.RungeKutta;
                    c.ExplicitOrder = 1;
                    c.EquationOfState = IdealGas.Air;

                    double gamma = c.EquationOfState.HeatCapacityRatio;

                    c.MachNumber = 1.0;
                    double MaSquared = c.MachNumber * c.MachNumber;


                    c.AddVariable(Variables.Density, dgDegree);
                    c.AddVariable(Variables.Momentum.xComponent, dgDegree);
                    c.AddVariable(Variables.Energy, dgDegree);
                    c.AddVariable(Variables.Velocity.xComponent, dgDegree);
                    c.AddVariable(Variables.Pressure, dgDegree);
                    c.AddVariable(Variables.LocalMachNumber, dgDegree);

                    Func<double[], double> exactDensity = X =>
                        (2.0 + Math.Cos(2 * X[0])) / 2.0;
                    Func<double[], double> exactMomentum = X =>
                        (7.0 - Math.Cos(4.0 * X[0])) / 16.0;
                    Func<double[], double> exactEnergy = X =>
                            (25.0 + 7.0 / 16.0 * MaSquared * (-2.0 + Math.Cos(2.0 * X[0])) * (-2.0 + Math.Cos(2.0 * X[0]))) * (2.0 + Math.Cos(2.0 * X[0])) / 20.0;

                    // Add all types of boundary conditions for convenience
                    Func<double[], double> exactVelocity = X =>
                        exactMomentum(X) / exactDensity(X);
                    Func<double[], double> exactPressure = X =>
                        1 + Math.Cos(2.0 * X[0])/2 ;

                    c.InitialValues_Evaluators.Add(Variables.Density, exactDensity);
                    c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, exactVelocity);
                    c.InitialValues_Evaluators.Add(Variables.Pressure, exactPressure);

                    c.AddBoundaryCondition("supersonicInlet", Variables.Density, exactDensity);
                    c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], exactVelocity);
                    c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, exactPressure);

                    c.AddBoundaryCondition("subsonicInlet", Variables.Density, exactDensity);
                    c.AddBoundaryCondition("subsonicInlet", Variables.Velocity[0], exactVelocity);

                    c.AddBoundaryCondition("subsonicOutlet", Variables.Pressure, exactPressure);

                    Func<double[], double> localMach = X => c.MachNumber * exactVelocity(X) / (Math.Sqrt(exactPressure(X) / exactDensity(X)));
                    // Spurk: p305, eq(9.100)
                    Func<double[], double> totalPressure = X => exactPressure(X) * Math.Pow((gamma - 1.0) / 2 * localMach(X) * localMach(X) + 1.0, gamma / (gamma - 1.0));
                    // Spurk: p305, eq(9.101)
                    Func<double[], double> totalTemperature = X => exactPressure(X) / exactDensity(X) * ((gamma - 1.0) / 2 * localMach(X) * localMach(X) + 1.0);
                    c.AddBoundaryCondition("subsonicPressureInlet", "p0", totalPressure);
                    c.AddBoundaryCondition("subsonicPressureInlet", "T0", totalTemperature);

                    c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                        (x, t, state) => -Math.Sin(4.0 * x[0]) / 4.0));
                    c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                        (x, t, state) => (160.0 - 35.0 * MaSquared - 56.0 * MaSquared * Math.Cos(2.0 * x[0]) + 21.0 * MaSquared * Math.Cos(4.0 * x[0])) * Math.Sin(2.0 * x[0]) / (224.0 * MaSquared)));


                    c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                            (x, t, state) => -(7.0 / 640.0) * Math.Sin(2 * x[0]) * ((160.0 + 3.0 * MaSquared) * Math.Cos(2 * x[0]) + MaSquared * (10.0 - 6.0 * Math.Cos(4.0 * x[0]) + Math.Cos(6.0 * x[0])))));


                    c.Queries.Add("densityError", QueryLibrary.L2Error(Variables.Density, exactDensity, 10));
                    c.Queries.Add("momentumError", QueryLibrary.L2Error(Variables.Momentum[0], exactMomentum, 10));
                    c.Queries.Add("energyError", QueryLibrary.L2Error(Variables.Energy, exactEnergy, 10));

                    c.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("divisions", divisions),
                        new Tuple<string, object>("dgDegree", dgDegree)
                    };

                    c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
                    c.ResidualBasedTerminationCriteria.Add("changeRate_abs_rhoE", 1E-9);

                    c.dtMin = 0.0;
                    c.dtMax = 1.0;
                    c.CFLFraction = 0.1;
                    c.Endtime = 10000;
                    c.NoOfTimesteps = 1000000;

                    c.ProjectName = "MMS1D_Euler_Mach=" + c.MachNumber + "_dg=" + dgDegree + "_cells=" + noOfCells;

                    result.Add(c);
                }
            }

            return result.ToArray();
        }

        /// <summary>
        /// Test using <see cref="SupersonicInlet"/> (Dirichlet) everywhere
        /// </summary>
        /// <returns></returns>
        public static CNSControl[] EulerSupersonicInlet1D() {
            CNSControl[] templates = GetTemplates();

            foreach (CNSControl c in templates) {
                int divisions = (int)c.Paramstudy_CaseIdentification.Single(t => t.Item1 == "divisions").Item2;

                int noOfCells = (2 << divisions) * 8;
                c.GridFunc = delegate {
                    GridCommons grid = Grid1D.LineGrid(
                        GenericBlas.Linspace(0.0, Math.PI / 2.0 + 0.0, noOfCells + 1));
                    grid.EdgeTagNames.Add(1, "supersonicInlet");
                    grid.DefineEdgeTags(x => 1);
                    return grid;
                };
                c.ProjectName += "_supersonicAll";
            }

            return templates.ToArray();
        }

        /// <summary>
        /// Test using <see cref="SupersonicInlet"/> (Dirichlet) at the inlet
        /// and <see cref="SubsonicOutlet"/> at the outlet
        /// </summary>
        /// <returns></returns>
        public static CNSControl[] EulerSubsonicOutlet1D() {
            CNSControl[] templates = GetTemplates();

            foreach (CNSControl c in templates) {
                int divisions = (int)c.Paramstudy_CaseIdentification.Single(t => t.Item1 == "divisions").Item2;

                int noOfCells = (2 << divisions) * 8;
                c.GridFunc = delegate {
                    GridCommons grid = Grid1D.LineGrid(
                        GenericBlas.Linspace(0.0, Math.PI / 2.0 + 0.0, noOfCells + 1));
                    grid.EdgeTagNames.Add(1, "supersonicInlet");
                    grid.EdgeTagNames.Add(2, "subsonicOutlet");
                    grid.DefineEdgeTags(x => Math.Abs(x[0]) < 1e-14 ? (byte)1 : (byte)2);
                    return grid;
                };
                c.ProjectName += "_subsonicOutlet";
            }

            return templates;
        }

        /// <summary>
        /// Test using <see cref="SubsonicInlet"/> at the inlet and
        /// <see cref="SupersonicInlet"/> (Dirichlet) at the outlet
        /// </summary>
        /// <returns></returns>
        public static CNSControl[] EulerSubsonicInlet1D() {
            CNSControl[] templates = GetTemplates();

            foreach (CNSControl c in templates) {
                int divisions = (int)c.Paramstudy_CaseIdentification.Single(t => t.Item1 == "divisions").Item2;

                int noOfCells = (2 << divisions) * 8;
                c.GridFunc = delegate {
                    GridCommons grid = Grid1D.LineGrid(
                        GenericBlas.Linspace(0.0, Math.PI / 2.0 + 0.0, noOfCells + 1));
                    grid.EdgeTagNames.Add(1, "supersonicInlet");
                    grid.EdgeTagNames.Add(2, "subsonicInlet");
                    grid.DefineEdgeTags(x => Math.Abs(x[0]) < 1e-14 ? (byte)2 : (byte)1);
                    return grid;
                };
                c.ProjectName += "_subsonicInlet2";
            }

            return templates;
        }

        /// <summary>
        /// Test using <see cref="SubsonicInlet"/> at the inlet and
        /// <see cref="SubsonicOutlet"/> at the outlet. That is, this test case
        /// uses physically correct boundary conditions
        /// </summary>
        /// <returns></returns>
        public static CNSControl[] EulerSubsonicInletAndOutlet1D() {
            CNSControl[] templates = GetTemplates();

            foreach (CNSControl c in templates) {
                int divisions = (int)c.Paramstudy_CaseIdentification.Single(t => t.Item1 == "divisions").Item2;

                int noOfCells = (2 << divisions) * 8;
                c.GridFunc = delegate {
                    GridCommons grid = Grid1D.LineGrid(
                        GenericBlas.Linspace(0.0, Math.PI / 2.0 + 0.0, noOfCells + 1));
                    grid.EdgeTagNames.Add(1, "subsonicInlet");
                    grid.EdgeTagNames.Add(2, "subsonicOutlet");
                    grid.DefineEdgeTags(x => Math.Abs(x[0]) < 1e-14 ? (byte)1 : (byte)2);
                    return grid;
                };
                c.ProjectName += "_subsonicAll";
            }

            return templates;
        }

        /// <summary>
        /// Uses <see cref="SubsonicPressureInlet"/> at the inlet and
        /// <see cref="SupersonicInlet"/> (Dirichlet) at the outlet
        /// </summary>
        /// <returns></returns>
        public static CNSControl[] EulerSubsonicPressureInletTest1D() {
            CNSControl[] templates = GetTemplates();

            foreach (CNSControl c in templates) {
                int divisions = (int)c.Paramstudy_CaseIdentification.Single(t => t.Item1 == "divisions").Item2;

                int noOfCells = (2 << divisions) * 8;
                c.GridFunc = delegate {
                    GridCommons grid = Grid1D.LineGrid(
                        GenericBlas.Linspace(0.0, Math.PI / 2.0 + 0.0, noOfCells + 1));
                    grid.EdgeTagNames.Add(1, "supersonicInlet");
                    grid.EdgeTagNames.Add(2, "subsonicPressureInlet");
                    grid.DefineEdgeTags(x => Math.Abs(x[0]) < 1e-14 ? (byte)2 : (byte)1);
                    return grid;
                };
            }

            return templates;
        }

        /// <summary>
        /// Test using <see cref="SubsonicPressureInlet"/> at the inlet and
        /// <see cref="SubsonicOutlet"/> at the outlet. That is, this test case
        /// uses physically correct boundary conditions
        /// </summary>
        /// <returns></returns>
        public static CNSControl[] EulerSubsonicPressureInletAndOutletTest1D() {
            CNSControl[] templates = GetTemplates();

            foreach (CNSControl c in templates) {
                int divisions = (int)c.Paramstudy_CaseIdentification.Single(t => t.Item1 == "divisions").Item2;

                int noOfCells = (2 << divisions) * 8;
                c.GridFunc = delegate {
                    GridCommons grid = Grid1D.LineGrid(
                        GenericBlas.Linspace(0.0, Math.PI / 2.0 + 0.0, noOfCells + 1));
                    grid.EdgeTagNames.Add(1, "subsonicPressureInlet");
                    grid.EdgeTagNames.Add(2, "subsonicOutlet");
                    grid.DefineEdgeTags(x => Math.Abs(x[0]) < 1e-14 ? (byte)1 : (byte)2);
                    return grid;
                };
            }

            return templates;
        }
    }
}
