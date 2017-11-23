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
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.Source;
using ilPSP.Utils;
using System;
using static BoSSS.Foundation.Grid.Classic.GridCommons;

namespace CNS.Tests.MMS {

    /// <summary>
    /// Various unsteady Manufactured solutions
    /// </summary>
    public static class MMS_unsteady {

        /// <summary>
        /// Performs a convergence study for the compressible Navier-Stokes equations in 2D
        /// and uses a manufactured solution <see cref="Gassner2D_conserved(int, int, string)"/>
        /// </summary>
        /// <param name="noOfRefinements"></param>
        /// <param name="maxDegree"></param>
        /// <param name="minDegree"></param>
        /// <param name="dbPath"></param>
        /// <returns></returns>
        public static CNSControl[] Gassner2DStudy_conserved(int noOfRefinements, int maxDegree, int minDegree=1, string dbPath = @"c:\bosss_dbv2\exp") {
            CNSControl[] controls = new CNSControl[(maxDegree+1-minDegree) * noOfRefinements];
            int ii = 0;
            for (int i = minDegree; i <= maxDegree; i++) {
                for (int j = 0; j < noOfRefinements; j++) {
                    double power = 3 + (double)j;
                    int noOfCellsPerDirection = (int)Math.Pow(2, power);
                    controls[ii] = Gassner2D_conserved(noOfCellsPerDirection, i, dbPath);
                    controls[ii].savetodb = true;
                    controls[ii].Paramstudy_ContinueOnError = true;
                    controls[ii].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("divisions", i),
                        new Tuple<string, object>("dgDegree", i)
                    };
                    ii++;
                }
            }
            return controls;
        }

        /// <summary>
        /// Perfoms a complete time convergence study, i.e on the same grid several runs are done with decreasing constant timestep.
        /// </summary>
        /// <param name="timeLevel">Number of refinement levels</param>
        /// <param name="dtStart">Starting timestep, i.e coarsest timesteps</param>
        /// <param name="order">Explicite timestepping order</param>
        /// <param name="timeStepper">Explicte timestepping scheme, default: Local Timestepping (LTS)</param>
        /// <returns></returns>
        public static CNSControl[] Gassner2DStudy_time(int timeLevel, double dtStart, int order, string timeStepper="LTS") {
            CNSControl[] controls = new CNSControl[timeLevel];
            for (int i=0; i < timeLevel; i++) {
                double factor = Math.Pow(2, i);
                double dt = dtStart / factor;
                controls[i] = Gassner2D_time(dt, order, timeStepper);
                controls[i].Paramstudy_ContinueOnError = true;
                controls[i].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("dt", dt)
                };
            }
            return controls;
        }

        /// <summary>
        /// Fixed run for a convergence study of Local Timestepping (LTS) for 2 and 3 order
        /// -c cs:CNS.Tests.MMS.MMS_unsteady.Gassner2DStudy_LTS(4)
        /// </summary>
        /// <param name="timeLevel">Number of refinement levels</param>
        /// <returns></returns>
        public static CNSControl[] Gassner2DStudy_LTS(int timeLevel) {
            CNSControl[] controls = new CNSControl[2* timeLevel];
            int ii = 0;
            double dt =3.8E-4; 
            CNSControl[] OneOrderSet = Gassner2DStudy_time(timeLevel, dt, 3);
            foreach (CNSControl run in OneOrderSet) {
                controls[ii] = run;
                ii++;
            }
            OneOrderSet = Gassner2DStudy_time(timeLevel, dt, 2);
            foreach (CNSControl run in OneOrderSet) {
                controls[ii] = run;
                ii++;
            }

            return controls;
        }

        /// <summary>
        /// Does one run with a given timestep on a mesh with a refinment in the middel, i.e 96 cells with h=1/10
        /// and 16 cells with h=1/20, DGorder=11! It is a very fine spatial resolution to measure time discetizations errors
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="timeStepper">Explicte timestepping scheme</param>
        /// <param name="order">Explicite timestepping order</param>
        /// <returns></returns>
        public static CNSControl Gassner2D_time(double dt, int order, string timeStepper) {
            CNSControl c = Gassner2D_conserved(10, 9);
            c.GridFunc = delegate {
                var grid = Grid2D.HangingNodes2D(true, true,
                                new GridBox(0, 0, 1, 1, 10),
                                new GridBox(0.4, 0.4, 0.6, 0.6, 4)
                                );
                grid.Name = "[0,1]x[0,1], h=1/10, Ratio 1:2";
                return grid;
            };

            c.ReynoldsNumber = 1E7;

            c.CFLFraction = 0.2;
            c.dtMax = dt;

            c.NoOfTimesteps = int.MaxValue;
            c.NumberOfSubGrids = 2;
            c.Endtime = 0.1;

            c.savetodb = true;

            switch (timeStepper) {
                case "LTS":
                    c.ExplicitScheme = ExplicitSchemes.LTS;
                    break;
                case "AB":
                case "AdamsBashforth":
                    c.ExplicitScheme = ExplicitSchemes.AdamsBashforth;
                    break;
                case "RK":
                case "RungeKutta":
                    c.ExplicitScheme = ExplicitSchemes.RungeKutta;
                    break;
            }
            c.ResidualLoggerType = ResidualLoggerTypes.Query;
            c.ResidualInterval = 100;
            c.PrintInterval = 100;
            c.ExplicitOrder = order;
            c.ProjectName = "TimeMMS_" + c.ExplicitScheme+ c.ExplicitOrder + "_dt=" + dt;
            return c;
        }

        /// <summary>
        /// 2D Manufactured Solution for the compressible Navier-Stokes equations with variable Mach number
        /// and constant viscosity. It is based on conserved variables. Further details: Gassner, G., Lörcher, F., &amp; Munz, C. D. (2008). A discontinuous 
        /// Galerkin scheme based on a space-time expansion II. Viscous flow equations in multi dimensions. 
        /// Journal of Scientific Computing, 34(3), 260-286.
        /// </summary>
        /// <param name="noOfCellsPerDirection"></param>
        /// <param name="dgDegree"></param>
        /// <param name="dbPath"></param>
        /// <returns></returns>
        public static CNSControl Gassner2D_conserved(int noOfCellsPerDirection, int dgDegree, string dbPath = @"c:\bosss_dbv2\exp") {
            CNSControl c = new CNSControl();

            // Session Settings
            c.DbPath = dbPath;
            c.savetodb = false;
            c.saveperiod = 10000;
            c.ProjectDescription = "Manufactured Solution by Gassner et al. (2008)";
            c.Tags.Add("MMS");
            c.Tags.Add("Diffusive Test");

            // Solver Type
            c.ActiveOperators = Operators.Convection | Operators.Diffusion | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
            c.SIPGPenaltyScaling = 1.3;

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.1;
            c.CFLFraction = 0.25;
            c.NoOfTimesteps = 2000;

            c.PrintInterval = 100;
            c.ResidualInterval = 1000;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;

            // Time-Stepping Settings
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;


            //Material Settings
            c.EquationOfState = IdealGas.Air;
            c.ViscosityLaw = new ConstantViscosity();
            c.PrandtlNumber = 0.72;
            c.ReynoldsNumber = 1000.0;
            c.MachNumber = 0.5;
        
            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            // Parameters
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);


            // Grid
            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: true);
                return grid;
            };

            // Functions
            double b = 0.25;
            double c1 = 1.0;
            double k = 2.0 * Math.PI;
            double omega = 20.0 * Math.PI;
            double gamma = c.EquationOfState.HeatCapacityRatio;
            double MachScaling = gamma * c.MachNumber * c.MachNumber;

            Func<double[], double, double> rho = (X, t) => c1 + b * Math.Sin(k * (X[0] + X[1]) - omega * t);
            Func<double[], double, double> m = (X, t) => rho(X, t);
            Func<double[], double, double> rhoE = (X, t) => rho(X, t) * rho(X, t);

            Func<double[], double, double> u1 = (X,t) => 1.0;
            Func<double[], double, double> u2 = (X, t) => 1.0;
            Func<double[], double, double> pressure = (X, t) => (gamma - 1.0) * (rhoE(X, t) - 0.5 * MachScaling * rho(X, t) * (u1(X, t) * u1(X, t) + u2(X, t) * u2(X, t))); 

            //Initial Values
            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => m(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => m(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Energy, X => rhoE(X, 0.0));

            //c.InitialValues.Add(Variables.Velocity.xComponent, X => u0(X, 0.0));
            //c.InitialValues.Add(Variables.Velocity.yComponent, X => u1(X, 0.0));
            //c.InitialValues.Add(Variables.Pressure, X => pressure(X, 0.0));


            //MMS Sources

            c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega + 2.0 * k)));
            c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega + 2.0 * k - (gamma - 1.0) * k + (gamma - 1.0) /MachScaling * (2.0 * k * rho(X, t)))
                            ));
            c.CustomMomentumSources[1].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega + 2.0 * k - (gamma - 1.0)* k + (gamma - 1.0) / MachScaling * (2.0 * k * rho(X, t)))
                            ));
            c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-2.0 * omega * rho(X, t) + 4.0 * k * rho(X, t) + 2.0 * (gamma - 1.0) * (2.0 * k * rho(X, t) - MachScaling * k)) - 2.0 * gamma * b * k * k / (c.PrandtlNumber * c.ReynoldsNumber) * Math.Sin(k * (X[0] + X[1]) - omega * t)
                            ));


            // Queries
            int QueryDegree = 14;
            c.Queries.Add("densityError", QueryLibrary.L2Error(Variables.Density, rho, QueryDegree));
            c.Queries.Add("momentumError", QueryLibrary.L2Error(Variables.Momentum.xComponent, m, QueryDegree));
            c.Queries.Add("energyError", QueryLibrary.L2Error(Variables.Energy, rhoE, QueryDegree));
            c.Queries.Add("velocityError", QueryLibrary.L2Error(Variables.Velocity.xComponent, u1, QueryDegree));

            c.ProjectName = "MMS-Gassner2D_Mach=" + c.MachNumber + "_Flux=" + c.DiffusiveFluxType + "_h=1/" + noOfCellsPerDirection + "_p=" + dgDegree;


            return c;
        }

        /// <summary>
        /// 2D Manufactured Solution for the compressible Navier-Stokes equations with variable Mach number
        /// and constant viscosity. It is based on primitive variables.
        /// </summary>
        /// <param name="noOfCellsPerDirection"></param>
        /// <param name="dgDegree"></param>
        /// <param name="dbPath"></param>
        /// <returns></returns>
        public static CNSControl Gassner2D_primitive(int noOfCellsPerDirection, int dgDegree, string dbPath = @"c:\bosss_dbv2\exp") {
            CNSControl c = new CNSControl();

            // Session Settings
            c.DbPath = dbPath;
            c.savetodb = true;
            c.saveperiod = 10000;
            c.ProjectDescription = "Manufactured Solution with primitive variables. Based on Gassner et al. (2008)";
            c.Tags.Add("MMS");
            c.Tags.Add("Diffusive Test");

            // Solver Type
            c.ActiveOperators = Operators.Convection | Operators.Diffusion | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
            c.SIPGPenaltyScaling = 1.3;

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.1;
            c.CFLFraction = 0.25;
            c.NoOfTimesteps = 2000;

            c.PrintInterval = 100;
            c.ResidualInterval = 1000;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;

            // Time-Stepping Settings
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;


            //Material Settings
            c.EquationOfState = IdealGas.Air;
            c.ViscosityLaw = new ConstantViscosity();
            c.PrandtlNumber = 0.72;
            c.ReynoldsNumber = 1.0e6;
            c.MachNumber = 1/Math.Sqrt(1.4);

            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            // Parameters
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);


            // Grid
            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: true);
                return grid;
            };

            // Functions

            double b = 0.25;
            double c1 = 2.0;
            double k = 2.0 * Math.PI;
            double omega = 20.0 * Math.PI;
            double gamma = c.EquationOfState.HeatCapacityRatio;
            double MachScaling = gamma * c.MachNumber * c.MachNumber;

            double PressureScale = 1.0;

            Func<double[], double, double> rho = (X, t) => c1 + b * Math.Sin(k * (X[0] + X[1]) - omega * t);
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 1.0;
            Func<double[], double, double> pressure = (X, t) => rho(X, t)/PressureScale;

            Func<double[], double, double> m0 = (X, t) => u0(X, t) * rho(X, t);
            Func<double[], double, double> m1 = (X, t) => u1(X, t) * rho(X, t);
            Func<double[], double, double> rhoE = (X, t) => pressure(X, t) / (gamma - 1) + 0.5 * MachScaling * rho(X, t) * (u0(X, t) * u0(X, t) + u1(X, t) * u1(X, t));

            //Initial Values
            //c.InitialValues.Add(Variables.Density, X => rho(X, 0.0));
            //c.InitialValues.Add(Variables.Momentum.xComponent, X => m0(X, 0.0));
            //c.InitialValues.Add(Variables.Momentum.yComponent, X => m1(X, 0.0));
            //c.InitialValues.Add(Variables.Energy, X => rhoE(X, 0.0));

            //c.InitialValues.Add(Variables.Velocity.xComponent, X => u0(X, 0.0));
            //c.InitialValues.Add(Variables.Velocity.yComponent, X => u1(X, 0.0));
            //c.InitialValues.Add(Variables.Pressure, X => pressure(X, 0.0));


            //MMS Sources

            c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega + 2.0 * k)));
            c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega + 2.0 * k + k / (PressureScale * MachScaling))
                            ));
            c.CustomMomentumSources[1].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega + 2.0 * k + k / (PressureScale * MachScaling))
                            ));
            c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                (X, t, state) => -b * Math.Cos(k * (X[0] + X[1]) - omega * t) * (-omega / (PressureScale * (gamma - 1.0)) - MachScaling * omega + 2.0 * k / (PressureScale * (gamma - 1.0)) + 1.0 + MachScaling + 2.0 + k / PressureScale)
                            ));


            // Queries
            int QueryDegree = 10;
            c.Queries.Add("densityError", QueryLibrary.L2Error(Variables.Density, rho, QueryDegree));
            c.Queries.Add("momentum0Error", QueryLibrary.L2Error(Variables.Momentum.xComponent, m0, QueryDegree));
            c.Queries.Add("energyError", QueryLibrary.L2Error(Variables.Energy, rhoE, QueryDegree));
            c.Queries.Add("velocity0Error", QueryLibrary.L2Error(Variables.Velocity.xComponent, u0, QueryDegree));

            c.ProjectName = "MMS-Gassner2D_Mach=" + c.MachNumber + "_Flux=" + c.DiffusiveFluxType + "_h=1/" + noOfCellsPerDirection + "_p=" + dgDegree;


            return c;
        }

        public static CNSControl[] Gassner3DStudy(int noOfRefinements, int maxDegree, int minDegree = 2, string dbPath = @"c:\bosss_dbv2\MMS\Gassner3d_unsteady") {
            CNSControl[] controls = new CNSControl[(maxDegree + 1 - minDegree) * noOfRefinements];
            int ii = 0;
            for (int i = minDegree; i <= maxDegree; i++) {
                for (int j = 0; j < noOfRefinements; j++) {
                    double power = 3 + (double)j;
                    int noOfCellsPerDirection = (int)Math.Pow(2, power);
                    controls[ii] = Gassner3D_unsteady(noOfCellsPerDirection, i, dbPath);
                    controls[ii].savetodb = true;
                    controls[ii].Paramstudy_ContinueOnError = true;
                    controls[ii].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                        new Tuple<string, object>("divisions", i),
                        new Tuple<string, object>("dgDegree", i)
                    };
                    ii++;
                }
            }
            return controls;
        }

        public static CNSControl Gassner3D_unsteady(int noOfCellsPerDirection, int dgDegree, string dbPath= @"/home/kraemer/CNS/dbe") {
            CNSControl c = new CNSControl();

            // Session Settings
            c.DbPath = dbPath;
            c.savetodb = true;
            c.saveperiod = 10000;
            c.ProjectName = "MMS_Gassner3D_h=1/" + noOfCellsPerDirection + "_p=" + dgDegree;
            c.ProjectDescription = "Manufactured Solution by Gassner et al. (2008)";
            c.Tags.Add("MMS");
            c.Tags.Add("Diffusive Test");
            c.GridPartType = GridPartType.ParMETIS;

            // Solver Type
            c.ActiveOperators = Operators.Convection | Operators.Diffusion | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
            c.SIPGPenaltyScaling = 1.3;

            // Time-Stepping Settings
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;

            //Material Settings
            c.EquationOfState = IdealGas.Air;
            c.ViscosityLaw = new ConstantViscosity();
            c.PrandtlNumber = 0.72;
            c.ReynoldsNumber = 1000.0;

            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Momentum.zComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            // Parameters
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Velocity.zComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);


            // Grid
            //c.GridFunc = delegate {
            //    double[] nodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);
            //    var grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes, periodicX: true, periodicY: true, periodicZ: true, _CellType: CellType.Cube_Linear);
            //    return grid;
            //};
            String grid = "";

            switch (noOfCellsPerDirection) {
                case 16:
                    grid = "06d43b3a-c630-445b-b17f-5c0a17d7e290";
                    break;
                case 32:
                    grid = "c3276077-946a-4360-9917-eba9a87fcbfc";
                    break;
                case 64:
                    grid = "c80eac65-6e56-420d-ae65-be8d75001539";
                    break;
            }

            Guid gridGuid;
            if (Guid.TryParse(grid, out gridGuid)) {
                c.GridGuid = gridGuid;
            } else {
                throw new Exception(
                 "Could not find a grid at " + grid);
            }

            // Functions
            Func<double[], double, double> rho = (X, t) => 2.0 + 0.25 * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 0 * Math.PI * t);
            Func<double[], double, double> m = (X, t) => 2.0 + 0.25 * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 0 * Math.PI * t);
            Func<double[], double, double> rhoE = (X, t) => (2.0 + 0.25 * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 0 * Math.PI * t)) * (2.0 + 0.25 * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 0 * Math.PI * t));

            Func<double[], double, double> u = (X, t) => 1.0;
            Func<double[], double, double> pressure = (X, t) => (1.4 - 1) * Math.Pow((2.0 + 0.25 * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t)), 2)
                                                                    - 3 / 2 * (2.0 + 0.25 * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t));

            //Initial Values
            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => m(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => m(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Momentum.zComponent, X => m(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Energy, X => rhoE(X, 0.0));

            //c.InitialValues.Add(Variable.Velocity.xComponent, X => u0(X, 0.0));
            //c.InitialValues.Add(Variable.Velocity.yComponent, X => u1(X, 0.0));
            //c.InitialValues.Add(Variable.Velocity.zComponent, X => u2(X, 0.0));
            //c.InitialValues.Add(Variable.Pressure, X => pressure(X, 0.0));

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.1;
            c.CFLFraction = 0.1;
            c.NoOfTimesteps = 5000000;

            c.PrintInterval = 10;
            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;

            //MMS Sources 
            c.CustomContinuitySources.Add(map => new AdHocSourceTerm(map,
                (X, t, state) => 3.5 * Math.PI * Math.Cos(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t)));
            c.CustomMomentumSources[0].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => 3 * Math.PI * Math.Cos(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t) -
                                    0.05 * Math.PI * Math.Sin(4 * Math.PI * (X[0] + X[1] + X[2]) - 40 * Math.PI * t)));
            c.CustomMomentumSources[1].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => 3 * Math.PI * Math.Cos(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t) -
                                    0.05 * Math.PI * Math.Sin(4 * Math.PI * (X[0] + X[1] + X[2]) - 40 * Math.PI * t)));
            c.CustomMomentumSources[2].Add(map => new AdHocSourceTerm(map,
                (X, t, state) => 3 * Math.PI * Math.Cos(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t) -
                                    0.05 * Math.PI * Math.Sin(4 * Math.PI * (X[0] + X[1] + X[2]) - 40 * Math.PI * t)));
            c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
                (X, t, state) => 12.5 * Math.PI * Math.Cos(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t)
                                    + 0.725 * Math.PI * Math.Sin(4 * Math.PI * (X[0] + X[1] + X[2]) - 40 * Math.PI * t)
                                    - 0.005833333332 * Math.PI * Math.PI * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 20 * Math.PI * t)));
            //c.CustomEnergySources.Add(map => new AdHocSourceTerm(map,
            //    (X, t, state) => -0.005833333332 * Math.PI * Math.PI * Math.Sin(2 * Math.PI * (X[0] + X[1] + X[2]) - 0 * Math.PI * t)));


            // Queries
            int QueryDegree = 10;
            c.Queries.Add("densityError", QueryLibrary.L2Error(Variables.Density, rho, QueryDegree));
            c.Queries.Add("momentum0Error", QueryLibrary.L2Error(Variables.Momentum.xComponent, m, QueryDegree));
            c.Queries.Add("energyError", QueryLibrary.L2Error(Variables.Energy, rhoE, QueryDegree));
            c.Queries.Add("velocity0Error", QueryLibrary.L2Error(Variables.Velocity.xComponent, u, QueryDegree));


            return c;
        }
    }
}
