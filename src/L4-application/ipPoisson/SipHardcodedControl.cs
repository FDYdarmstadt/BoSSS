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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Application.SipPoisson {

    /// <summary>
    /// predefined control-objects
    /// </summary>
    static public class SipHardcodedControl {

        /// <summary>
        /// Test on a curved grid.
        /// </summary>
        public static SipControl TestCurved() {
            var R = new SipControl();
            R.ProjectName = "ipPoison/curved";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 15 });
            R.InitialValues_Evaluators.Add("RHS", X => 0.0);
            R.InitialValues_Evaluators.Add("Tex", X => (Math.Log(X[0].Pow2() + X[1].Pow2()) / Math.Log(4.0)) + 1.0);
            R.ExactSolution_provided = true;

            R.GridFunc = delegate() {
                var grd = Grid2D.CurvedSquareGrid(GenericBlas.Linspace(1, 2, 3), GenericBlas.Linspace(0, 1, 11), CellType.Square_9, true);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(X => 1);
                return grd;
            };

            R.AddBoundaryCondition(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     double x = X[0], y = X[1];
                     return Math.Sqrt(x * x + y * y);
                 });

            R.solver_name = null;

            return R;
        }

        /// <summary>
        /// Creates Nodes, yes it really does!
        /// </summary>
        /// <param name="res"></param>
        /// <param name="stetch">
        /// Factor which determines how much the intervals in the output grow; 1.0 is no stretching.
        /// </param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        static double[] CreateNodes(int res, double stetch, double min, double max) {
            if(stetch == 1.0)
                return GenericBlas.Linspace(min, max, res + 1);
            else
                return Grid1D.ExponentialSpaceing(min, max, res + 1, stetch); // without proper preconditioning,
            // a stretched grid is much more expensive than
            // an equidistant grid !!!
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// </summary>
        public static SipControl TestCartesian1(int xRes = 32, double xStretch = 1.0, int yRes = 16, double yStretch = 1.01, int pDG = 2) {
            var RR = new SipControl();
            RR.ProjectName = "ipPoison/cartesian";
            RR.savetodb = false;

            RR.FieldOptions.Add("T", new FieldOpts() { Degree = pDG, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            RR.FieldOptions.Add("Tex", new FieldOpts() { Degree = pDG * 2 });
            RR.InitialValues_Evaluators.Add("RHS", X => 1.0);
            RR.InitialValues_Evaluators.Add("Tex", X => (0.5 * X[0].Pow2() - 10 * X[0]));
            RR.ExactSolution_provided = true;

            RR.GridFunc = delegate() {
                double[] xNodes = CreateNodes(xRes, xStretch, 0, 10);
                double[] yNodes = CreateNodes(yRes, yStretch, -1, +1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });
                
                return grd;
            };


            RR.AddBoundaryCondition(BoundaryType.Dirichlet.ToString());
            RR.AddBoundaryCondition(BoundaryType.Neumann.ToString());

            //RR.solver_name = "direct";
            RR.solver_name = null;

            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;


            return RR;
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// </summary>
        public static SipControl TestCartesian3D(int xRes = 32, double xStretch = 1.0, int yRes = 16, double yStretch = 1.0, int zRes = 16, double zStretch = 1.0) {
            var R = new SipControl();
            R.ProjectName = "ipPoison/cartesian";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = 6, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 6 });
            R.InitialValues_Evaluators.Add("RHS", X => 1.0);
            R.InitialValues_Evaluators.Add("Tex", X => (0.5 * X[0].Pow2() - 10 * X[0]));
            R.ExactSolution_provided = true;

            R.GridFunc = delegate() {
                double[] xNodes = CreateNodes(xRes, xStretch, 0, 10);
                double[] yNodes = CreateNodes(yRes, yStretch, -1, +1);
                double[] zNodes = CreateNodes(zRes, zStretch, -1, +1);

                var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };
            
            R.AddBoundaryCondition(BoundaryType.Dirichlet.ToString());
            R.AddBoundaryCondition(BoundaryType.Neumann.ToString());

            R.solver_name = null;

            return R;
        }


        /// <summary>
        /// Test on a Cartesian grid, with a sinusodial solution.
        /// </summary>
        /// <param name="Res">
        /// Grid resolution, 2 entries for 2D test, 3 entries for 3D test.
        /// </param>
        /// <param name="Stretch">
        /// Grid stretching factors.
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        public static SipControl TestCartesian2(int[] Res, double[] Stretch = null, string solver_name = "softpcg+schwarz+directcoarse", int deg = 3) {
            if(Res.Length != 2 && Res.Length != 3)
                throw new ArgumentOutOfRangeException();
            if(Stretch == null) {
                Stretch = new double[Res.Length];
                Stretch.SetAll(1.0);
            } else {
                if(Res.Length != Stretch.Length)
                    throw new ArgumentException();
            }

            var R = new SipControl();
            R.ProjectName = "ipPoison/cartesian";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg*2 });
            R.InitialValues_Evaluators.Add("RHS", X => -Math.Sin(X[0]));
            R.InitialValues_Evaluators.Add("Tex", X => Math.Sin(X[0]));
            R.ExactSolution_provided = true;
            R.NoOfMultigridLevels = 3;

            R.GridFunc = delegate() {
                GridCommons grd = null;
                if(Res.Length == 2) {
                    double[] xNodes = CreateNodes(Res[0], Stretch[0], 0, 10);
                    double[] yNodes = CreateNodes(Res[1], Stretch[1], -1, +1);

                    grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                } else if(Res.Length == 3) {
                    double[] xNodes = CreateNodes(Res[0], Stretch[0], 0, 10);
                    double[] yNodes = CreateNodes(Res[1], Stretch[1], -1, +1);
                    double[] zNodes = CreateNodes(Res[2], Stretch[2], -1, +1);

                    grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                } else {
                    throw new NotSupportedException();
                }
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            R.AddBoundaryCondition(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     double x = X[0], y = X[1];
                     return 0.0;
                 });

            R.AddBoundaryCondition(BoundaryType.Neumann.ToString(), "T",
                 delegate (double[] X) {
                     if(Math.Abs(X[1] - 1.0) < 1.0e-8 || Math.Abs(X[1] + 1.0) < 1.0e-8)
                         return 0;

                     if(X.Length > 2 && (Math.Abs(X[2] - 1.0) < 1.0e-8 || Math.Abs(X[2] + 1.0) < 1.0e-8))
                         return 0;

                     return Math.Cos(10.0);
                 });


            R.solver_name = solver_name;

            R.NoOfMultigridLevels = 3;

            return R;
        }


        /// <summary>
        /// Poisson Equation on a (-1,1)x(-1,1), Dirichlet everywhere
        /// </summary>
        public static SipControl Square(int xRes = 21, int yRes = 16, int deg = 5) {

            //Func<double[], double> exRhs = X => 2 * X[0] * X[0] + 2 * X[1] * X[1] - 4;
            //Func<double[], double> exSol = X => (1.0 - X[0] * X[0]) * (1.0 - X[1] * X[1]);

            //Func<double[], double> exSol = X => (1.0 - X[1]);
            //Func<double[], double> exRhs = X => 0.0;

            Func<double[], double> exSol = X => -Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5);
            Func<double[], double> exRhs = X => (Math.PI * Math.PI * 0.5 * Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5)); // == - /\ exSol


            var R = new SipControl();
            R.ProjectName = "ipPoison/square";
            R.savetodb = false;
            //R.DbPath = "D:\\BoSSS-db";

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 4 });
            R.InitialValues_Evaluators.Add("RHS", exRhs);
            R.InitialValues_Evaluators.Add("Tex", exSol);
            R.ExactSolution_provided = true;

            R.GridFunc = delegate() {
                double[] xNodes = GenericBlas.Linspace(-1,1,xRes);
                double[] yNodes = GenericBlas.Linspace(-1,1,yRes);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret = 1;
                    return ret;
                });


                return grd;
            };

            R.AddBoundaryCondition(BoundaryType.Dirichlet.ToString(), "T", exSol);

            R.solver_name = null;
            R.NoOfSolverRuns = 1;

            return R;
        }
    }
}
