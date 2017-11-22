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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using CNS.Residual;
using ilPSP.Utils;
using System;

namespace CNS.Tests.IBMTests {

    public static class GaussianBump {


        private static double[] bumpBottom(double x, double y, double height) {
            return new double[2] { x, y + (height - y) / height * 0.3939 * Math.Exp(-0.5 * x * x) };
        }

        public static IBMControl[] IBMBumpStudy(int noOfGridLevels, int maxDgDegree, int lsDegree, double CFL, double epsX = 0.0, double epsY = 0.0, int minDgDegree = 0) {
            IBMControl[] controls = new IBMControl[noOfGridLevels * ((maxDgDegree - minDgDegree) + 1)];

            for (int i = minDgDegree; i <= maxDgDegree; i++) {
                for (int j = 0; j < noOfGridLevels; j++) {
                    int noOfCells = (int)Math.Pow(2, 3 + j);
                    controls[(i - minDgDegree) * noOfGridLevels + j] = IBMBump(noOfCells, i, lsDegree, CFL, epsX, epsY);
                    controls[(i - minDgDegree) * noOfGridLevels + j].Paramstudy_ContinueOnError = true;
                    controls[(i - minDgDegree) * noOfGridLevels + j].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                    new Tuple<string, object>("refinement+dgDegree", i*100+noOfCells)
                };
                }
            }

            return controls;
        }

        public static CNSControl[] BumpStudy(int noOfGridLevels, int maxDgDegree, double CFL, int minDgDegree = 0) {
            CNSControl[] controls = new CNSControl[noOfGridLevels * ((maxDgDegree - minDgDegree) + 1)];

            for (int i = minDgDegree; i <= maxDgDegree; i++) {
                for (int j = 0; j < noOfGridLevels; j++) {
                    int noOfCells = (int)Math.Pow(2, 3 + j);
                    controls[(i - minDgDegree) * noOfGridLevels + j] = Bump(noOfCells, i, CFL);
                    controls[(i - minDgDegree) * noOfGridLevels + j].Paramstudy_ContinueOnError = true;
                    controls[(i - minDgDegree) * noOfGridLevels + j].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                    new Tuple<string, object>("refinement+dgDegree", i*100+noOfCells)
                };
                }
            }

            return controls;
        }


        public static IBMControl IBMBump(int noOfCells, int dgDegree, int lsDegree, double CFL, double epsilonX = 0.0, double epsilonY = 0.0) {
            IBMControl c = new IBMControl();

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 1000.0;
            c.CFLFraction = CFL;
            c.NoOfTimesteps = 200000;

            c.PrintInterval = 100;
            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("changeRate_L2_abs_rhoE", 1E-8);

            //IBM Settings
            c.LevelSetBoundaryTag = "adiabaticSlipWall";
            c.LevelSetQuadratureOrder = 2 * dgDegree;
            c.AgglomerationThreshold = 0.3;

            // NEXT STEP: SET THIS BOOL TO FALSE AND JUST USE IN POSITIVE SUB_VOLUME;
            // THEN TRY BOUNDING BOX APPROACH?
            // WHY THE HELL DOES THIS CONFIGURATION FAIL!??!?!?!?
            c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.SurfaceHMF_ProjectNodesToLevelSet = false;
            c.SurfaceHMF_RestrictNodes = true;
            c.SurfaceHMF_UseGaussNodes = false;
            c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            c.VolumeHMF_RestrictNodes = true;
            c.VolumeHMF_UseGaussNodes = false;

            //Guid restart = new Guid(" 60688cbc-707d-4777-98e6-d237796ec14c");
            //c.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restart, -1);

            // Session Settings
            c.DbPath = @"\\fdyprime\userspace\kraemer-eis\FDY-Cluster\dbe_bump\";
            //c.DbPath = @"/home/kraemer/GaussianBump/dbev2/";
            c.savetodb = true;
            c.saveperiod = 20000;
            c.ProjectName = "BoxHMF=" + c.MomentFittingVariant + "_Ma=0.5_(" + 2 * noOfCells + "x" + noOfCells + ")_CFL=" + c.CFLFraction + "_lsQuadOrder=" + c.LevelSetQuadratureOrder + "_p=" + dgDegree + "_agg=" + c.AgglomerationThreshold + "_epsX=" + epsilonX + "_epsY=" + epsilonY;
            c.ProjectDescription = "GaussianBump with Ma=0.5";
            c.Tags.Add("Gaussian Bump");
            c.Tags.Add("IBM Test");

            // Solver Type
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Time-Stepping Settings
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            //Material Settings
            c.EquationOfState = IdealGas.Air;



            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(IBMVariables.LevelSet, lsDegree);



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
                double xBoundary = 12.0;
                double yBoundary = 12.0;
                double yBottom = 0.0;

                double[] xnodes = GenericBlas.Linspace(-xBoundary, xBoundary, 2 * noOfCells + 1);

                //double ySplit = 6.0;
                //int ySplitNoOfCells = (int) (0.5*noOfCells);
                //double[] ynodes1 = GenericBlas.Linspace(yBottom, ySplit, ySplitNoOfCells + 1);
                //double[] ynodes2 = GenericBlas.Linspace(ySplit, yBoundary, noOfCells-ySplitNoOfCells + 1);
                //ynodes1 = ynodes1.GetSubVector(0, ynodes1.Length - 1);
                //double[] ynodes = ArrayTools.Cat(ynodes1, ynodes2);

                double[] ynodes = GenericBlas.Linspace(yBottom, yBoundary, noOfCells + 1);

                GridCommons grid = Grid2D.Cartesian2DGrid(
                    xnodes,
                    ynodes
                    );

                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.EdgeTagNames.Add(2, "adiabaticSlipWall");

                Func<double[], byte> func = delegate (double[] x) {

                    if (Math.Abs(x[0] + xBoundary) < 1e-5) { // Inflow
                        return 1;
                    } else if (Math.Abs(x[0] - xBoundary) < 1e-5) { // Outflow
                        return 1;
                    } else if (Math.Abs(x[1] - yBoundary) < 1e-5) { // Top
                        return 1;
                    } else { // Bottom
                        return 2;
                    }
                };
                grid.DefineEdgeTags(func);
                grid.Name = "IBM-[" + -xBoundary + "," + xBoundary + "]x[" + yBottom + "," + yBoundary + "]_Cells:(" + 2 * noOfCells + "x" + noOfCells + ")";
                return grid;
            };

            // Functions
            Func<double[], double, double> rho = (X, t) => 1.0;
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 0.0;
            Func<double[], double, double> pressure = (X, t) => 2.8571428571428;

            //Initial Values
            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u0(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => u1(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressure(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - epsilonY - 0.01 - 0.3939 * Math.Exp(-0.5 * (X[0] - epsilonX) * (X[0] - epsilonX));

            //BoundaryConditions
            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity.xComponent, u0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity.yComponent, u1);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, pressure);


            // Queries

            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 2.8571428571428));
            return c;
        }

        public static CNSControl Bump(int noOfCells, int dgDegree, double CFL) {
            CNSControl c = new CNSControl();

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 1000.0;
            c.CFLFraction = CFL;
            c.NoOfTimesteps = 150000;

            c.PrintInterval = 100;
            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("changeRate_L2_abs_rhoE", 1E-8);

            // Session Settings
            c.DbPath = @"C:\bosss_dbv2\GaussianBump";
            c.savetodb = true;
            c.saveperiod = 5000;
            c.ProjectName = "GaussianBump-Lieb_Ma=0.5_(" + 2 * noOfCells + "x" + noOfCells + ")_CFL=" + c.CFLFraction + "_p=" + dgDegree;
            c.ProjectDescription = "GaussianBump with Ma=0.5";
            c.Tags.Add("Gaussian Bump");

            // Solver Type
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Time-Stepping Settings
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            //Material Settings
            c.EquationOfState = IdealGas.Air;



            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            // Parameters
            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);



            // Grid
            // Load from db
            switch (noOfCells) {
                case 8:
                    c.GridGuid = new Guid("75d75bdc-1ca6-47e2-a2ec-394cd5728e9c");
                    break;
                case 16:
                    c.GridGuid = new Guid("20baeb6f-62d9-4276-8100-d391dcbb07dd");
                    break;
                case 32:
                    c.GridGuid = new Guid("773b127e-4591-4044-8b6b-6ec46a627d5e");
                    break;
                case 64:
                    c.GridGuid = new Guid("84e5271f-637b-40b8-801b-64109a0af4d0");
                    break;
                default:
                    throw new ArgumentException("Wrong Grid Input");
            }

            //// generate new grid
            //c.GridFunc = delegate {
            //    double xBoundary = 12.0;
            //    double yBoundary = 12.0;
            //    double yBottom = 0.0;

            //    double[] xnodes = GenericBlas.Linspace(-xBoundary, xBoundary, 2 * noOfCells + 1);
            //    double[] ynodes = GenericBlas.Linspace(yBottom, yBoundary, noOfCells + 1);

            //    GridCommons grid = Grid2D.CurvedSquareGridChannel(xnodes, ynodes, CellType.Square_12, false, (x, y) => bumpBottom(x, y, yBoundary));

            //    grid.EdgeTagNames.Add(1, "supersonicinlet");
            //    grid.EdgeTagNames.Add(2, "adiabaticSlipWall");

            //    Func<double[], byte> func = delegate(double[] x) {

            //        if (Math.Abs(x[0] + xBoundary) < 1e-5) { // Inflow
            //            return 1;
            //        } else if (Math.Abs(x[0] - xBoundary) < 1e-5) { // Outflow
            //            return 1;
            //        } else if (Math.Abs(x[1] - yBoundary) < 1e-5) { // Top
            //            return 1;
            //        } else { // Bottom
            //            return 2;
            //        }
            //    };
            //    grid.DefineEdgeTags(func);
            //    grid.Name = "[" + -xBoundary + "," + xBoundary + "]x[" + yBottom + "," + yBoundary + "]_Cells:(" + 2 * noOfCells + "x" + noOfCells + ")";
            //    return grid;
            //};

            // Functions
            Func<double[], double, double> rho = (X, t) => 1.0;
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 0.0;
            Func<double[], double, double> pressure = (X, t) => 2.8571428571428;

            //Initial Values
            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u0(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => u1(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressure(X, 0.0));

            //BoundaryConditions
            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity.xComponent, u0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity.yComponent, u1);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, pressure);


            // Queries
            Func<double[], double, double> solEntropy = (X, t) => 2.8571428571428;

            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(Variables.Entropy, solEntropy, 10));
            return c;
        }

        public static IBMControl IBMBumpTest(int noOfCells, int dgDegree, int lsDegree, double CFL) {
            IBMControl c = new IBMControl();

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 1000.0;
            c.CFLFraction = CFL;
            c.NoOfTimesteps = 250000;

            c.PrintInterval = 100;
            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("changeRate_L2_abs_rhoE", 1E-8);

            //Guid restart = new Guid(" 60688cbc-707d-4777-98e6-d237796ec14c");
            //c.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restart, -1);

            // Session Settings
            c.DbPath = @"C:\bosss_dbv2\GaussianBump_hhlr";
            //c.DbPath = @"\\fdyprime\userspace\kraemer-eis\FDY-Cluster\dbe_bump\";
            //c.DbPath = @"/home/kraemer/GaussianBump/dbev2/";
            c.savetodb = true;
            c.saveperiod = 100;
            c.ProjectName = "TestsCutCells_(" + 2 * noOfCells + "x" + noOfCells + ")_CFL=" + c.CFLFraction + "_ls=" + lsDegree + "_p=" + dgDegree + "_agg=" + 0.53;
            c.ProjectDescription = "GaussianBump with Ma=0.5";
            c.Tags.Add("Gaussian Bump");
            c.Tags.Add("IBM Test");

            // Solver Type
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Time-Stepping Settings
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            //Material Settings
            c.EquationOfState = IdealGas.Air;

            //IBM Settings
            c.LevelSetBoundaryTag = "adiabaticSlipWall";
            c.LevelSetQuadratureOrder = 2 * lsDegree;
            c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.AgglomerationThreshold = 0.3;

            // Primary Variables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(IBMVariables.LevelSet, lsDegree);

            c.GridFunc = delegate {
                double xBoundary = 2.0625;
                double yBoundary = 2.0525;
                double yBottom = -0.01;

                double[] xnodes = GenericBlas.Linspace(-xBoundary, xBoundary, 2 * noOfCells + 1);
                double[] ynodes = GenericBlas.Linspace(yBottom, yBoundary, noOfCells + 1);

                GridCommons grid = Grid2D.Cartesian2DGrid(
                    xnodes,
                    ynodes
                    );

                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.EdgeTagNames.Add(2, "adiabaticSlipWall");

                Func<double[], byte> func = delegate (double[] x) {

                    if (Math.Abs(x[0] + xBoundary) < 1e-5) { // Inflow
                        return 1;
                    } else if (Math.Abs(x[0] - xBoundary) < 1e-5) { // Outflow
                        return 1;
                    } else if (Math.Abs(x[1] - yBoundary) < 1e-5) { // Top
                        return 1;
                    } else { // Bottom
                        return 2;
                    }
                };
                grid.DefineEdgeTags(func);
                grid.Name = "IBM-[" + -xBoundary + "," + xBoundary + "]x[" + yBottom + "," + yBoundary + "]_Cells:(" + 2 * noOfCells + "x" + noOfCells + ")";
                return grid;
            };

            // Functions
            Func<double[], double, double> rho = (X, t) => 1.0;
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 0.0;
            Func<double[], double, double> pressure = (X, t) => 2.8571428571428;

            //Initial Values
            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u0(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => u1(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressure(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - 0.3939 * Math.Exp(-0.5 * X[0] * X[0]);

            //BoundaryConditions
            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity.xComponent, u0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity.yComponent, u1);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, pressure);


            // Queries

            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 2.8571428571428));
            return c;
        }
    }
}
