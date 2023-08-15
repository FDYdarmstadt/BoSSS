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
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Connectors.Matlab;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.LevelSetTools.TestCases;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;



namespace BoSSS.Application.CahnHilliard {

    /// <summary>
    /// predefined control-objects
    /// </summary>
    static public class Examples {


        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// Non-dimenzionalized with interface thickness xi=1
        /// </summary>
        public static CahnHilliardControl TestCartesian(int xRes = 20, int yRes = 20, int pDG = 3) {
            var RR = new CahnHilliardControl();
            RR.ProjectName = "CahnHilliard/cartesian";

            RR.savetodb = false;

            RR.SetDGdegree(pDG);

            double S = 1.0;
            RR.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-S, S, xRes+1);
                double[] yNodes = GenericBlas.Linspace(-S, S, yRes+1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Wall.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Outflow.ToString());

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(Math.Abs(X[1]) - S) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(Math.Abs(X[0]) - S) <= 1.0e-8)
                        et = 1;

                    return et;
                });
                return grd;
            };

            double gridsize = 2.0 * S / xRes;
            // interface thickness in cells
            double multiplier = 4.0; // 2.0
            double thickness = multiplier * (double)4.0/pDG; //initially tried 3.0/pDG
            RR.cahn = 2 * gridsize * 1/4.164;//thickness * (1/4.164 * gridsize);

            RR.diff = 1.0; // 1.0 / RR.cahn.Pow2();// RR.cahn/(4.0 * mu);
            RR.lambda = 0.0;

            double radius = 0.3; // 0.275

            //RR.AddBoundaryValue(BoundaryType.Flow.ToString(), "c", X => Math.Tanh( (-Math.Sqrt(Math.Pow(X[1] - 1.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)));
            //RR.AddBoundaryValue(BoundaryType.Flow.ToString(), "c", X => -1.0);
            //RR.AddInitialValue("c", new Formula("X => Math.Sqrt(X[0]*X[0] + X[1]*X[1]) < 0.5 ? 1.0 : -1.0", false));
            //RR.AddInitialValue("c", new Formula("X => X[0] < 0.5 ? (X[0] > -0.5 ? (X[1] > -0.5 ? (X[1] < 0.5 ? 1.0 : -1.0) : -1.0) : -1.0) : -1.0", false));
            //RR.InitialValues_Evaluators.Add("c", X => X[0] < -0.5 ? -1.0 : Math.Tanh((-Math.Sqrt(Math.Pow(X[0] + 0.5, 2.0) + Math.Pow(X[1] - 0.0, 2.0)) + radius)/(Math.Sqrt(2) * RR.cahn)));
            //RR.InitialValues_Evaluators.Add("c", X => X[0] < -radius ? Math.Tanh((-Math.Sqrt(Math.Pow(X[0] + radius, 2.0) + Math.Pow(X[1] - 0.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)) : X[0] > radius ? Math.Tanh((-Math.Sqrt(Math.Pow(X[0] - radius, 2.0) + Math.Pow(X[1] - 0.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)) : Math.Tanh((-Math.Sqrt(Math.Pow(X[1], 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)));
            //RR.AddInitialValue("c", new Formula("X => X[0] < 0.5 ? (X[0] > -0.8 ? (X[1] > -0.5 ? (X[1] < 0.5 ? 1.0 : -1.0) : -1.0) : -1.0) : -1.0", false));
            RR.InitialValues_Evaluators.Add("c", X =>  Math.Tanh((-Math.Sqrt(Math.Pow(X[0] - 0.0, 2.0) + Math.Pow(X[1] - 0.0, 2.0)) + radius)/(Math.Sqrt(2) * RR.cahn)));
            //RR.AddInitialValue(VariableNames.VelocityX, new Formula("X => X[1] ", false));
            //RR.InitialValues_Evaluators.Add(VariableNames.VelocityX, X =>  Math.Exp(-Math.Pow(X[0],2.0)));
            //RR.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 2.0 * X[1] * X[0] * Math.Exp(-Math.Pow(X[0],2.0)));
            //RR.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => X[1] < 0.0 ? X[0] : 1.0 * X[0]);
            //RR.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => -X[1] );
            //RR.InitialValues_Evaluators.Add("c", X => X[0] < 0.0 ? -1.0 : 1.0);
            //int longation = 4;
            //RR.InitialValues_Evaluators.Add("c", X => X[0] < -longation * radius ? Math.Tanh((-Math.Sqrt(Math.Pow(X[0] + longation * radius, 2.0) + Math.Pow(X[1] - 0.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)) : X[0] > longation * radius ? Math.Tanh((-Math.Sqrt(Math.Pow(X[0] - longation * radius, 2.0) + Math.Pow(X[1] - 0.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)) : Math.Tanh((-Math.Sqrt(Math.Pow(X[1], 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)));

            RR.ModTyp = CahnHilliardControl.ModelType.modelB;
            RR.CorrectionType = CahnHilliardControl.Correction.None;
            // RR.CurvatureCorrection = false;
            // RR.UseDirectCurvature = false;
            RR.UseFDJacobian = false;

            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            RR.SuperSampling = 4;
            RR.AdaptiveMeshRefinement = false;

            RR.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            RR.dtFixed = 0.1;
            RR.NoOfTimesteps = 100;

            RR.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            RR.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            RR.NonLinearSolver.MaxSolverIterations = 25;

            RR.includeConvection = false;

            return RR;
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// Non-dimenzionalized with interface thickness xi=1
        /// </summary>
        public static CahnHilliardControl TestCartesian3D(int xRes = 20, int yRes = 20, int zRes = 20, int pDG = 3)
        {
            var RR = new CahnHilliardControl();
            RR.ProjectName = "CahnHilliard/cartesian";

            RR.savetodb = false;

            RR.SetDGdegree(pDG);

            RR.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-1.0, 1.0, xRes + 1);
                double[] yNodes = GenericBlas.Linspace(-1.0, 1.0, yRes + 1);
                double[] zNodes = GenericBlas.Linspace(-1.0, 1.0, zRes + 1);

                var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Wall.ToString());
                //grd.EdgeTagNames.Add(2, BoundaryType.Flow.ToString());

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });
                return grd;
            };

            double gridsize = 2.0 / xRes;
            // interface thickness in cells
            double thickness = 2.0 * (double)3.0 / pDG;
            RR.cahn = thickness * (0.25 * gridsize);
            RR.diff = 1.0;
            double radius = 0.5;

            //RR.AddBoundaryValue(BoundaryType.Flow.ToString(), "c", X => Math.Tanh( (-Math.Sqrt(Math.Pow(X[1] - 1.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)));
            //RR.AddBoundaryValue(BoundaryType.Flow.ToString(), "c", X => -1.0);
            //RR.AddInitialValue("c", new Formula("X => Math.Sqrt(X[0]*X[0] + X[1]*X[1]) < 0.5 ? 1.0 : -1.0", false));
            //RR.AddInitialValue("c", new Formula("X => X[0] < 0.5 ? (X[0] > -0.5 ? (X[1] > -0.5 ? (X[1] < 0.5 ? 1.0 : -1.0) : -1.0) : -1.0) : -1.0", false));
            //RR.AddInitialValue("c", new Formula("X => X[0] < 0.5 ? (X[0] > -0.8 ? (X[1] > -0.5 ? (X[1] < 0.5 ? 1.0 : -1.0) : -1.0) : -1.0) : -1.0", false));
            RR.InitialValues_Evaluators.Add("c", X => Math.Tanh((-Math.Sqrt(Math.Pow(X[0] - 0.3, 2.0) + Math.Pow(X[1] - 0.3, 2.0) + Math.Pow(X[2] - 0.0, 2.0)) + radius) / (Math.Sqrt(2) * RR.cahn)));
            //RR.AddInitialValue(VariableNames.VelocityX, new Formula("X => X[1] ", false));
            RR.AddInitialValue(VariableNames.VelocityX, new Formula("X => 0.0", false));
            RR.AddInitialValue(VariableNames.VelocityY, new Formula("X => 0.0 ", false));

            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;

            RR.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            RR.dtFixed = 0.01;
            RR.NoOfTimesteps = 1000;
            RR.SuperSampling = 2;
            //RR.cahn = 1;
            RR.AdaptiveMeshRefinement = false;
            RR.includeConvection = true;



            return RR;
        }

        public static CahnHilliardControl EllipticDroplet(int xRes = 90, int yRes = 90, int pDG = 2)
        {
            // Console.ReadLine();
            var RR = new CahnHilliardControl();
            RR.ProjectName = "CahnHilliard/cartesian";

            RR.savetodb = false;

            // RR.CurvatureCorrection = false;
            // RR.UseDirectCurvature = false;
            // RR.ModTyp = CahnHilliardControl.ModelType.modelB;

            RR.SetDGdegree(pDG);

            RR.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-15, 15, xRes + 1);
                double[] yNodes = GenericBlas.Linspace(-15, 15, yRes + 1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Wall.ToString());
                grd.DefineEdgeTags(delegate (double[] X) { return (byte)1; });

                return grd;
            };

            RR.AddBoundaryValue(BoundaryType.Wall.ToString(), "c", new Formula("X => -1"));

            //RR.AddInitialValue("c", new Formula("X => (X[0]*X[0] + X[1]*X[1]) < 0.25 ? 1.0 : -1.0", false));
            RR.AddInitialValue("c", new Formula("X => -Math.Tanh((Math.Sqrt(X[0]*X[0]*0.75 + X[1]*X[1])-5)*Math.Sqrt(2))"));
            RR.AddInitialValue(VariableNames.VelocityX, new Formula("X => 0.0 ", false));
            RR.AddInitialValue(VariableNames.VelocityY, new Formula("X => 0.0 ", false));

            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;

            RR.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            RR.dtFixed = 1e-1;
            RR.cahn = 1;
            RR.diff = 0.1;
            
            return RR;
        }

        public static CahnHilliardControl RotatingDisk(int xRes = 20, int yRes = 20, int pDG = 6)
        {
            var RR = new CahnHilliardControl();
            RR.ProjectName = "CahnHilliard/cartesian";
            RR.includeDiffusion = false;

            RR.savetodb = false;

            RR.SetDGdegree(pDG);

            RR.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-15, 15, xRes + 1);
                double[] yNodes = GenericBlas.Linspace(-15, 15, yRes + 1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Wall.ToString());
                grd.DefineEdgeTags(delegate (double[] X) { return (byte)1; });

                return grd;
            };
            
            RR.AddBoundaryValue(BoundaryType.Wall.ToString(), "c", new Formula("X => -1"));

            //RR.AddInitialValue("c", new Formula("X => (X[0]*X[0] + X[1]*X[1]) < 0.25 ? 1.0 : -1.0", false));
            RR.AddInitialValue("c", new Formula("X => -Math.Tanh((Math.Sqrt(X[0]*X[0]*0.75 + X[1]*X[1])-5)*Math.Sqrt(2))"));
            //RR.AddInitialValue("c", new Formula("X => X[0] < 1 && X[0] > -1 && X[1] > 0 ? -1 : -Math.Tanh((Math.Sqrt(X[0]*X[0] + X[1]*X[1])-5)*Math.Sqrt(2))"));
            RR.AddInitialValue(VariableNames.VelocityX, new Formula("X => -X[1] ", false));
            RR.AddInitialValue(VariableNames.VelocityY, new Formula("X => X[0] ", false));

            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;

            RR.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            RR.dtFixed = 1e-2;
            RR.NoOfTimesteps = 1000;

            RR.cahn = 1;
            RR.diff = 0.0;

            return RR;
        }

        public static CahnHilliardControl Zalesak(int xRes = 40, int yRes = 40, int pDG = 4) {
            var RR = new CahnHilliardControl();
            RR.ProjectName = "CahnHilliard/cartesian";

            RR.savetodb = false;

            double[] XCutout = { 1.0, -1.0 };
            var Disk = new ZalesaksDisk(XCutout, 1.0, 5.0);

            RR.SetDGdegree(pDG);

            RR.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-15, 15, xRes + 1);
                double[] yNodes = GenericBlas.Linspace(-15, 15, yRes + 1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Wall.ToString());
                grd.DefineEdgeTags(delegate (double[] X) { return (byte)1; });

                return grd;
            };

            RR.AddBoundaryValue(BoundaryType.Wall.ToString(), "c", new Formula("X => -1"));

            //RR.AddInitialValue("c", new Formula("X => (X[0]*X[0] + X[1]*X[1]) < 0.25 ? 1.0 : -1.0", false));
            RR.AddInitialValue("c", new Formula("X => -Math.Tanh((Math.Sqrt(X[0]*X[0]*0.75 + X[1]*X[1])-5)*Math.Sqrt(2))"));
            RR.AddInitialValue(VariableNames.VelocityX, new Formula("X => 0.0 ", false));
            RR.AddInitialValue(VariableNames.VelocityY, new Formula("X => 0.0 ", false));

            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;

            RR.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            RR.dtFixed = 1e-1;
            RR.NoOfTimesteps = 100;

            RR.cahn = 1;
            RR.diff = 0.1;

            return RR;
        }


    }
}
