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
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.Solution;
using ilPSP.Utils;
using System;

namespace CNS.Tests.IBMTests {

    public static class NACA0012 {

        public static IBMControl[] IBMNACA0012Study(int noOfGridLevels, int maxDgDegree, double CFL, double agglomeration, double alpha, int minDgDegree = 0) {
            IBMControl[] controls = new IBMControl[noOfGridLevels * ((maxDgDegree - minDgDegree) + 1)];

            for (int i = minDgDegree; i <= maxDgDegree; i++) {
                for (int j = 0; j < noOfGridLevels; j++) {
                    int noOfCells = (int)Math.Pow(2, 3 + j);
                    controls[(i - minDgDegree) * noOfGridLevels + j] = IBMNACA0012(noOfCells, i, CFL, agglomeration, alpha);
                    controls[(i - minDgDegree) * noOfGridLevels + j].Paramstudy_ContinueOnError = true;
                    controls[(i - minDgDegree) * noOfGridLevels + j].Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                    new Tuple<string, object>("refinement+dgDegree", i*100+noOfCells)
                };
                }
            }

            return controls;
        }


        public static IBMControl IBMNACA0012(int MeshPara, int dgDegree, double CFL, double agglomeration, double alpha) {
            IBMControl c = new IBMControl();

            c.savetodb = true;
            
            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 1000.0;
            c.CFLFraction = CFL;
            c.NoOfTimesteps = 200000;

            c.PrintInterval = 10;
            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            c.ResidualBasedTerminationCriteria.Add("changeRate_L2_abs_rhoE", 1E-8);

            //IBM Settings
            c.LevelSetBoundaryTag = "adiabaticSlipWall";
            c.LevelSetQuadratureOrder = 2*dgDegree+2;
            c.AgglomerationThreshold = agglomeration;

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

            //Guid restart = new Guid("cd061fe3-3215-483a-8790-f6fd686d6676");
            //c.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restart, -1);

            // Session Settings
            c.DbPath = @"\\fdyprime\userspace\kraemer-eis\FDY-Cluster\dbe_NACA\";
            //c.DbPath = @"C:\bosss_dbv2\NACA0012";
            c.savetodb = true;
            c.saveperiod = 10000;
            c.ProjectName = "MeshPara:" + MeshPara  + "_CFL=" + c.CFLFraction + "_p=" + dgDegree + "_agg=" + c.AgglomerationThreshold + "_alpha="+ alpha + "_HMF="+ c.MomentFittingVariant;
            c.ProjectDescription = "NACA0012 Steady Test with Ma=0.5";
            c.Tags.Add("NACA0012");
            c.Tags.Add("IBM Test");
            c.Tags.Add("steady");

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

            c.AddVariable(IBMVariables.LevelSet, 8);

            // Refined Region
            double xBegin = -0.012;
            double xEnd = 1.01;

            c.GridFunc = delegate
            {

                int chords = 100;

                int xleft = -chords;
                int xRight = chords;
                int yBottom = -chords;
                int yTop = chords;

                double spacingFactor = 3.95;

                double[] xnodes1 = Grid1D.TanhSpacing(xleft, xBegin, MeshPara + 1, spacingFactor, false);
                double[] xnodes2 = Grid1D.TanhSpacing_DoubleSided(xBegin, xEnd, MeshPara + 1, 2.0);
                double[] xnodes3 = Grid1D.TanhSpacing(xEnd, xRight, MeshPara + 1, spacingFactor, true);

                double[] xComplete = new double[xnodes1.Length + xnodes2.Length + xnodes3.Length - 2];
                for (int i = 0; i < xnodes1.Length; i++) {
                    xComplete[i] = xnodes1[i];
                }
                for (int i = 1; i < xnodes2.Length; i++) {
                    xComplete[i + xnodes1.Length - 1] = xnodes2[i];
                }
                for (int i = 1; i < xnodes3.Length; i++) {
                    xComplete[i + xnodes1.Length + xnodes2.Length - 2] = xnodes3[i];
                }

                double yrefinedTop = 0.2;
                double yrefinedBottom = -0.2;

                double[] ynodes1 = Grid1D.TanhSpacing(yBottom, yrefinedBottom, MeshPara + 1, spacingFactor, false);
                double[] ynodes2 = GenericBlas.Linspace(yrefinedBottom, yrefinedTop, (int)(0.75 * MeshPara) + 1);
                double[] ynodes3 = Grid1D.TanhSpacing(yrefinedTop, yTop, MeshPara + 1, spacingFactor, true);

                double[] yComplete = new double[ynodes1.Length + ynodes2.Length + ynodes3.Length - 2];
                for (int i = 0; i < ynodes1.Length; i++) {
                    yComplete[i] = ynodes1[i];
                }
                for (int i = 1; i < ynodes2.Length; i++) {
                    yComplete[i + ynodes1.Length - 1] = ynodes2[i];
                }
                for (int i = 1; i < ynodes3.Length; i++) {
                    yComplete[i + ynodes1.Length + ynodes2.Length - 2] = ynodes3[i];
                }


                int numOfCellsX = (xRight - xleft) * MeshPara;
                int numOfCellsY = (yTop - yBottom) * MeshPara;

                GridCommons grid = Grid2D.Cartesian2DGrid(
                    xComplete,
                    yComplete
                    );

                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.DefineEdgeTags(x => 1);
                grid.Name = "[" + xleft + "," + xRight + "]x[" + yBottom + "," + yTop + "]_Cells:(" + (xComplete.Length-1) + "x" + (yComplete.Length-1) + ")";

                return grid;
            };



            // Functions
            Func<double[], double, double> rho = (X, t) => 1.0;
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 0.0;
            Func<double[], double, double> pressure = (X, t) => 2.8571428571428;

            Func<double[], double> test = X => 1 - 0.05 * 0.4 / 1.4 * Math.Pow(Math.Exp(1 - X[0] * X[0] - X[1] * X[1]), (1 / 0.4));

            Func<double[], double, double> levelSet = delegate (double[] X, double t) {
                double value = 0.0;

                double radian = alpha * Math.PI / 180;

                double xRotated =1 + Math.Cos(radian) * (X[0]-1) - Math.Sin(radian) * (X[1]);
                double yRotated = Math.Sin(radian) * (X[0]-1) + Math.Cos(radian) * (X[1]);

                double a = 0.6;
                //double b = 0.2969;
                double c1 = 0.126;
                double d = 0.3516;
                double e = 0.2843;
                double f = 0.1036;


                if (yRotated >= 0.0 || (X[0]>0.562875 && X[1]>0)) { 
                //if (yRotated >= 0.0 ){
                    value = Math.Pow((yRotated + a * (c1 * xRotated + d * Math.Pow(xRotated,2) - e * Math.Pow(xRotated, 3) + f * Math.Pow(xRotated, 4))),2) - 0.0317338596 * xRotated;
                } else {
                    value = Math.Pow((-yRotated + a * (c1 * xRotated + d * Math.Pow(xRotated, 2) - e * Math.Pow(xRotated, 3) + f * Math.Pow(xRotated, 4))), 2) - 0.0317338596 * xRotated;
                }

                //value = yRotated - Math.Tan(radian)*xRotated;

                return value;
            };

            c.LevelSetFunction = levelSet;

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

            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 2.8571428571428));
            c.Queries.Add("IBMDragForce", IBMQueries.IBMDragForce());
            c.Queries.Add("IBMLiftForce", IBMQueries.IBMLiftForce());
            return c;
        }
    }
}
