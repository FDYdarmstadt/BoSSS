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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.Source;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace CNS.Tests.MMS {
    /// <summary>
    /// Various Manufactured Solutions with steady state solutions
    /// </summary>
    class MMS_steady {

        /// <summary>
        /// Performs a convergence study for Euler equations in 1D and uses
        /// a manufactured solution <see cref="Euler1D(int, int)"/>
        /// </summary>
        /// <param name="maxDgDegree"></param>
        /// <param name="refinements"></param>
        /// <returns></returns>
        public static CNSControl[] Euler1D_Study(int maxDgDegree, int refinements) {
            List<CNSControl> controls = new List<CNSControl>();
            int ii = 0;
            for (int i = 0; i <= maxDgDegree; i++) {
                for (int j = 0; j < refinements; j++) {
                    int noOfCells = (int)Math.Pow(2, 3 + j);
                    CNSControl c = Euler1D(i, noOfCells);
                    c.Paramstudy_ContinueOnError = true;
                    c.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("dgDegree", i),
                        new Tuple<string, object>("refinement", j),
                    };
                    controls.Add(c);
                    ii++;
                }

            }
            return controls.ToArray();
        }

        /// <summary>
        /// 1D Manufactured solution for the Euler equations with variable Mach number
        /// </summary>
        /// <param name="dgDegree"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <returns></returns>
        public static CNSControl Euler1D(int dgDegree, int noOfCellsPerDirection) {
            CNSControl c = new CNSControl();

            c.PrintInterval = 1;
            c.ResidualInterval = 1;
            c.saveperiod = 1;

            c.DbPath = @"";
            c.savetodb = false;
            c.Tags.Add("MMS");

            c.ActiveOperators = Operators.Convection | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.Rusanov;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);
                var grid = Grid1D.LineGrid(nodes);
                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            double gamma = c.EquationOfState.HeatCapacityRatio;

            c.MachNumber = 1 / Math.Sqrt(gamma);

            double a = 2.0 * Math.PI;
            double b = 0.25;
            double c1 = 2.0;
            double MachScaling = (gamma * c.MachNumber * c.MachNumber);

            Func<double[], double, double> rho = (X, t) => c1 + b * Math.Sin(a * X[0]);
            Func<double[], double, double> m0 = (X, t) => c1 + b * Math.Sin(a * X[0]);
            Func<double[], double, double> rhoE = (X, t) => (c1 + b * Math.Sin(a * X[0])) * (c1 + b * Math.Sin(a * X[0]));
            Func<double[], double, double> u0 = (X, t) => 1;
            Func<double[], double, double> p = (X, t) => (gamma - 1) * (rhoE(X, t) - 0.5 * MachScaling * rho(X, t) * u0(X, t) * u0(X, t));

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => m0(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Energy, X => rhoE(X, 0.0));

            c.AddBoundaryValue("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity.xComponent, u0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, p);

            // MMS Sources
            c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                   (X, t, state) => -(a * b * Math.Cos(a * X[0]))
                   ));
            c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                   (X, t, state) => -((a * b * Math.Cos(a * X[0])) * (1 + (gamma - 1) / MachScaling * (2 * rho(X, t) - 0.5 * MachScaling)))
                   ));
            c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                    (X, t, state) => -((a * b * Math.Cos(a * X[0])) * (2 * rho(X, t) + (gamma - 1) * (2 * rho(X, t) - 0.5 * MachScaling)))
                    ));

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, rho));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(CNSVariables.Pressure, p));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(CNSVariables.Velocity.xComponent, u0));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.4;
            c.Endtime = 10000.0;
            c.NoOfTimesteps = int.MaxValue;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("changeRate_abs_rhoE", 1E-14);

            c.ProjectName = "MMS_Euler1D_Flux=" + c.ConvectiveFluxType + "dg=" + dgDegree + "_cells=" + noOfCellsPerDirection;

            return c;
        }

        /// <summary>
        /// Performs a convergence study for Euler equations in 2D and uses
        /// a manufactured solution <see cref="Euler2D(int, int)"/>
        /// </summary>
        /// <param name="maxDgDegree"></param>
        /// <param name="maxRefinements"></param>
        /// <param name="minDgDegree"></param>
        /// <param name="minRefinements"></param>
        /// <returns></returns>
        public static CNSControl[] Euler2D_Study(int maxDgDegree, int maxRefinements, int minDgDegree = 0, int minRefinements = 0) {
            List<CNSControl> controls = new List<CNSControl>();
            int ii = 0;
            for (int i = minDgDegree; i <= maxDgDegree; i++) {
                for (int j = minRefinements; j < maxRefinements; j++) {
                    int noOfCells = (int)Math.Pow(2, 3 + j);
                    CNSControl c = Euler2D(i, noOfCells);
                    c.Paramstudy_ContinueOnError = true;
                    c.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("dgDegree", i),
                        new Tuple<string, object>("refinement", j),
                    };
                    controls.Add(c);
                    ii++;
                }

            }
            return controls.ToArray();
        }

        /// <summary>
        /// 2D manufactured solution for the Euler equations with variable Mach number
        /// </summary>
        /// <param name="dgDegree"></param>
        /// <param name="noOfCellsPerDirection"></param>
        /// <returns></returns>
        public static CNSControl Euler2D(int dgDegree, int noOfCellsPerDirection) {
            CNSControl c = new CNSControl();

            c.PrintInterval = 100;
            c.ResidualInterval = 100;
            c.saveperiod = 100000;

            c.DbPath = @"\\fdyprime\scratch\kraemer-eis\bosss_dbv2";
            c.savetodb = true;
            c.Tags.Add("MMS");

            c.ActiveOperators = Operators.Convection | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(0.0, Math.PI / 2, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            double gamma = c.EquationOfState.HeatCapacityRatio;

            c.MachNumber = 1;
            double Mach = c.MachNumber;
            double MachSq = Mach * Mach;
            double MachScaling = (gamma * c.MachNumber * c.MachNumber);

            Func<double[], double, double> rho = (X, t) => 1.0 + 0.5 * Math.Cos(2.0 * (X[0] + X[1]));
            Func<double[], double, double> u0 = (X, t) => 1.0 - 0.25 * Math.Cos(2.0 * (X[0] + X[1]));
            Func<double[], double, double> u1 = (X, t) => 1.0 + 0.5 * Math.Cos(4.0 * (X[0] + X[1]));
            Func<double[], double, double> p = (X, t) => 1.0 + 0.5 * Math.Cos(2.0 * (X[0] + X[1]));

            Func<double[], double, double> rhoE = (X, t) => p(X, t) / (gamma - 1) + 0.5 * gamma * Mach * Mach * rho(X, t) * (u0(X, t) * u0(X, t) + u1(X, t) * u1(X, t));
            Func<double[], double, double> m0 = (X, t) => rho(X, t) * u0(X, t);
            Func<double[], double, double> m1 = (X, t) => rho(X, t) * u1(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => m0(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => m1(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Energy, X => rhoE(X, 0.0));

            c.AddBoundaryValue("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity.xComponent, u0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity.yComponent, u1);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, p);

            // MMS Sources
            c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                   (X, t, state) => (7.0 * Math.Sin(2.0 * (X[0] + X[1])) + 7.0 * Math.Sin(4.0 * (X[0] + X[1])) + 3.0 * Math.Sin(6.0 * (X[0] + X[1]))) / 4.0
                   ));
            c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                   (X, t, state) => ((160.0 + 245.0 * MachSq + 504.0 * MachSq * Math.Cos(2.0 * (X[0] + X[1])) + 189.0 * MachSq * Math.Cos(4.0 * (X[0] + X[1])) - 56.0 * MachSq * Math.Cos(6.0 * (X[0] + X[1]))) * Math.Sin(2.0 * (X[0] + X[1]))) / (224.0 * MachSq)
                   ));
            c.CustomMomentumSources[1].Add(map => new AdHocSourceTerm(map,
                    (X, t, state) => ((40.0 + 259.0 * MachSq + 728.0 * MachSq * Math.Cos(2.0 * (X[0] + X[1])) + 266.0 * MachSq * Math.Cos(4.0 * (X[0] + X[1])) + 98.0 * MachSq * Math.Cos(6.0 * (X[0] + X[1])) + 35.0 * MachSq * Math.Cos(8.0 * (X[0] + X[1]))) * Math.Sin(2.0 * (X[0] + X[1]))) / (56.0 * MachSq)
       ));
            c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                    (X, t, state) => (7.0 * (1600.0 + 1019.0 * MachSq + (2240.0 + 3762.0 * MachSq) * Math.Cos(2.0 * (X[0] + X[1])) + 12.0 * (80.0 + 113.0 * MachSq) * Math.Cos(4.0 * (X[0] + X[1])) + 978.0 * MachSq * Math.Cos(6.0 * (X[0] + X[1])) + 333.0 * MachSq * Math.Cos(8.0 * (X[0] + X[1])) + 84.0 * MachSq * Math.Cos(10.0 * (X[0] + X[1])) + 28.0 * MachSq * Math.Cos(12.0 * (X[0] + X[1]))) * Math.Sin(2.0 * (X[0] + X[1]))) / 1280.0
                    ));

            c.Queries.Add("densityError", QueryLibrary.L2Error(Variables.Density, rho));
            c.Queries.Add("momentum0Error", QueryLibrary.L2Error(Variables.Momentum.xComponent, m0));
            c.Queries.Add("momentum1Error", QueryLibrary.L2Error(Variables.Momentum.xComponent, m1));
            c.Queries.Add("energyError", QueryLibrary.L2Error(Variables.Energy, rhoE));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = 10000.0;
            c.NoOfTimesteps = int.MaxValue;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("changeRate_abs_rhoE", 1E-12);

            c.ProjectName = "MMS2D_Euler_Flux=" + c.ConvectiveFluxType + "_Mach=" + Mach + "_dg=" + dgDegree + "_cells=" + noOfCellsPerDirection;

            return c;
        }
    }
}
