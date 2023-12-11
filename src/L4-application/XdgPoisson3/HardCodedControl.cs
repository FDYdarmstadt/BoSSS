﻿/* =======================================================================
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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Application.XdgPoisson3 {

    /// <summary>
    /// Predefined example calculations.
    /// </summary>
    public static class HardCodedControl {

        /// <summary>
        /// A 45-degree interface in the 2D domain \f$ (-2,2)^2 \f$.
        /// </summary>
        public static XdgPoisson3Control Schraeg() {
            XdgPoisson3Control R = new XdgPoisson3Control();

            R.ProjectName = "XdgPoisson3/Schräg";
            R.savetodb = false;

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            R.GridFunc = delegate () {
                return Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 8), GenericBlas.Linspace(-2, 2, 8));
            };


            // set the level-set
            R.InitialValues_Evaluators.Add("Phi", X => (X[0] + X[1]) + 0.075);
            //R.InitialValues.Add("Phi", X => X[0] + +X[1]*0.001);
            R.ExcactSolSupported = true;
            R.InitialValues_Evaluators.Add("uEx#A", X => X[1]);
            R.InitialValues_Evaluators.Add("uEx#B", X => X[1]);
            R.InitialValues_Evaluators.Add("rhs#A", X => 0.0);
            R.InitialValues_Evaluators.Add("rhs#B", X => 0.0);


            R.MU_A = -1.0;
            R.MU_B = -1.0;

            R.xLaplaceBCs.g_Diri = delegate (CommonParamsBnd inp) {
                double y = inp.X[1];
                return y;
            };

            R.xLaplaceBCs.IsDirichlet = (inp => true);

            R.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            R.AgglomerationThreshold = 0.0;
            R.PrePreCond = MultigridOperator.Mode.IdMass;
            R.penalty_multiplyer = 1.1;
            R.ViscosityMode = XLaplace_Interface.Mode.SIP;

            return R;
        }

        /// <summary>
        /// A circular interface in the 2D domain \f$ (3/2,3/2)^2 \f$.
        /// </summary>
        public static XdgPoisson3Control Circle(int Resolution = 16, int p = 1, string DBPath = null, LinearSolverCode solver = LinearSolverCode.direct_pardiso) {
            //BoSSS.Application.XdgPoisson3.HardCodedControl.Circle(Resolution: 8);

            XdgPoisson3Control R = new XdgPoisson3Control();

            R.ProjectName = "XdgPoisson3/Circle";
            if (DBPath == null) {
                R.savetodb = false;
            } else {
                R.savetodb = false;
                R.DbPath = DBPath;
            }


            R.SetDGdegree(p);

            R.GridFunc = delegate () {
                return Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1.5, 1.5, Resolution + 1), GenericBlas.Linspace(-1.5, 1.5, Resolution + 1));
            };


            double RADIUS = 0.8;
            R.InitialValues_Evaluators.Add("Phi", X => -(X[0] - 0.0).Pow2() - (X[1] - 0.0).Pow2() + (RADIUS).Pow2());

            R.ExcactSolSupported = false;
            R.InitialValues_Evaluators.Add("uEx#A", X => 0.0);
            R.InitialValues_Evaluators.Add("uEx#B", X => 0.0);
            R.InitialValues_Evaluators.Add("u#A", X => 0.0);
            R.InitialValues_Evaluators.Add("u#B", X => 0.0);
            R.InitialValues_Evaluators.Add("rhs#A", X => +1.0);
            R.InitialValues_Evaluators.Add("rhs#B", X => +1.0);

            R.MU_A = -1.0;
            R.MU_B = -1000.0;

            R.xLaplaceBCs.g_Diri = (X => 0.0);
            R.xLaplaceBCs.IsDirichlet = (inp => true);

            R.LinearSolver = solver.GetConfig();

            R.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            R.AgglomerationThreshold = 0.1;
            R.PrePreCond = MultigridOperator.Mode.SymPart_DiagBlockEquilib;

            return R;
        }

        /*
        /// <summary>
        /// A parameter study over <see cref="Circle(int, int, string, string)"/>.
        /// </summary>
        public static IEnumerable<XdgPoisson3Control> CircleParameterStudy() {
            List<XdgPoisson3Control> R = new List<XdgPoisson3Control>();

            var hmfVersions = new[] {
                XQuadFactoryHelper.MomentFittingVariants.Classic,
                XQuadFactoryHelper.MomentFittingVariants.OneStepGauss,
                XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes };
            foreach (var HMFversion in hmfVersions) {
                foreach (int pOff in new int[] { -1, 0, 1, 2 }) {
                    foreach (int res in new int[] { 12, 24, 48, 96, 192 }) {

                        var C = Circle(Resolution: (res + 1));

                        C.Paramstudy_ContinueOnError = true;
                        C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Resolution", res));
                        C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("HMF", HMFversion));
                        C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("pOff", pOff));
                            

                        C.CutCellQuadratureType = HMFversion;
                        C.pOff = pOff;

                        R.Add(C);
                    }
                }
            }
            return R;
        }
        */

        /// <summary>
        /// A piecewise linear solution.
        /// </summary>
        public static XdgPoisson3Control PiecewiseLinear(double delta = 0.5) {
            XdgPoisson3Control R = new XdgPoisson3Control();

            R.ProjectName = "XdgPoisson3/PiecewiseLinear";
            R.savetodb = false;
            //R.DbPath = "C:\\BoSSS-db";

            R.SetDGdegree(1);



            R.GridFunc = delegate () {
                double[] xNodes = new double[] { -6, -3, -2, -1, 1, 2, 3, 6 };
                double[] yNodes = GenericBlas.Linspace(-3, 3, 2);
                return Grid2D.Cartesian2DGrid(xNodes, yNodes);
            };

            const double kA = -1;
            const double kB = -1;


            // set the level-set
            R.InitialValues_Evaluators.Add("Phi", X => (X[0] - delta));
            R.ExcactSolSupported = true;
            R.InitialValues_Evaluators.Add("uEx#A", X => kA * (X[0] - delta));
            R.InitialValues_Evaluators.Add("uEx#B", X => kB * (X[0] - delta));
            R.InitialValues_Evaluators.Add("rhs#A", X => 0.0);
            R.InitialValues_Evaluators.Add("rhs#B", X => 0.0);


            R.MU_A = -1.0 / kA;
            R.MU_B = -1.0 / kB;

            R.xLaplaceBCs.g_Diri = delegate (CommonParamsBnd inp) {
                double x = inp.X[0];
                if (x < 0)
                    return (x - delta) * kA;
                else
                    return (x - delta) * kB;
            };
            R.xLaplaceBCs.g_Neum = (inp => 0.0);


            R.xLaplaceBCs.IsDirichlet = delegate (CommonParamsBnd inp) {
                double x = inp.X[0];
                if (Math.Abs(x + 6) < 1.0e-6 || Math.Abs(x - 6) < 1.0e-6)
                    return true;
                else
                    return false;
            };

            R.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            R.AgglomerationThreshold = 0.0;

            return R;
        }

        /*
        /// <summary>
        /// A piecewise linear solution.
        /// </summary>
        public static XdgPoisson3Control MultigridVsAggregation(double delta = 0.1) {
            XdgPoisson3Control R = new XdgPoisson3Control();

            R.ProjectName = "XdgPoisson3/PiecewiseLinear";
            R.savetodb = false;
            //R.DbPath = "C:\\BoSSS-db";

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("uEx", new FieldOpts() {
                Degree = 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            //var cutout = new double[2, 2] { { 0, 0 }, { 3, 3 } };
            //var tt = new BoundingBox[] { new BoundingBox(2) };


            R.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-3, 3, 10);
                double[] yNodes = GenericBlas.Linspace(-3, 3, 10);
                return Grid2D.Cartesian2DGrid(xNodes, yNodes, CutOuts: new BoundingBox[] { new BoundingBox(new double[2, 2] { { 0, 0 }, { 3, 3 } }) });
            };
            R.LinearSolver.NoOfMultigridLevels = 3;

            const double kA = -1;
            const double kB = -2;

            // set the level-set
            R.InitialValues_Evaluators.Add("Phi", X => (X[0] - delta));
            R.ExcactSolSupported = true;
            R.InitialValues_Evaluators.Add("uEx#A", X => kA * (X[0] - delta));
            R.InitialValues_Evaluators.Add("uEx#B", X => kB * (X[0] - delta));
            R.InitialValues_Evaluators.Add("rhs#A", X => 0.0);
            R.InitialValues_Evaluators.Add("rhs#B", X => 0.0);
            R.ExcactSolSupported = true;

            R.MU_A = -1.0 / kA;
            R.MU_B = -1.0 / kB;

            R.xLaplaceBCs.g_Diri = delegate (CommonParamsBnd inp) {
                double x = inp.X[0];
                if (x < 0)
                    return (x - delta) * kA;
                else
                    return (x - delta) * kB;
            };
            R.xLaplaceBCs.g_Neum = (inp => 0.0);


            R.xLaplaceBCs.IsDirichlet = delegate (CommonParamsBnd inp) {
                double x = inp.X[0];
                if (Math.Abs(x + 3) < 1.0e-6 || Math.Abs(x - 3) < 1.0e-6)
                    return true;
                else
                    return false;
            };

            R.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();//R.solverName = "direct";
            R.AgglomerationThreshold = 0.2;


            return R;
        }
        */

        /// <summary>
        /// a piecewise parabolic solution.
        /// </summary>
        public static XdgPoisson3Control PiecewiseParabola(int dgDeg = 3, double delta = 0.0, double muA = -20, double muB = -1, double agg = 0.0) {
            // --control cs:BoSSS.Application.XdgPoisson3.HardCodedControl.PiecewiseParabola();
            XdgPoisson3Control R = new XdgPoisson3Control();

            if (delta < -1.5 || delta > 1.5)
                throw new ArgumentOutOfRangeException();

            R.ProjectName = "XdgPoisson3/PiecewiseParabola";
            R.savetodb = false;

            R.SetDGdegree(dgDeg);

            R.GridFunc = delegate () {
                double[] xNodes = new double[] { -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6 };
                double[] yNodes = GenericBlas.Linspace(-3, 3, 7);
                //double[] xNodes = new double[] { -6, -4,  -2, -1, 1,  6 };
                //double[] yNodes = GenericBlas.Linspace(-3, 3, 7);
                return Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicY: true);
            };

            R.MU_A = muA;
            R.MU_B = muB;

            double A0 = 36 / (muB + muA);
            double A1 = (3 * (-muB + muA)) / (muA * (muB + muA));
            double A2 = -1 / (2 * muA);
            double B0 = 36 / (muB + muA);
            double B1 = (3 * (-muB + muA)) / ((muB + muA) * muB);
            double B2 = -1 / (2 * muB);

            Func<double[], double> uEx_A = (X => -(A2 * (X[0] - delta).Pow2() + A1 * (X[0] - delta) + A0));
            Func<double[], double> uEx_B = (X => -(B2 * (X[0] - delta).Pow2() + B1 * (X[0] - delta) + B0));

            // set the level-set
            R.InitialValues_Evaluators.Add("Phi", X => (X[0] - delta));
            R.ExcactSolSupported = true;
            R.InitialValues_Evaluators.Add("uEx#A", uEx_A);
            R.InitialValues_Evaluators.Add("uEx#B", uEx_B);
            R.InitialValues_Evaluators.Add("rhs#A", X => +1.0);
            R.InitialValues_Evaluators.Add("rhs#B", X => +1.0);

            double Diri_left = uEx_A(new double[] { -6, 0 });
            double Diri_rigt = uEx_B(new double[] { +6, 0 });


            R.xLaplaceBCs.g_Diri = (inp => ((inp.X[0] < 0) ? Diri_left : Diri_rigt));
            R.xLaplaceBCs.g_Neum = (inp => 0.0);

            R.xLaplaceBCs.IsDirichlet = delegate (CommonParamsBnd inp) {
                double x = inp.X[0];
                if (Math.Abs(x + 6) < 1.0e-6 || Math.Abs(x - 6) < 1.0e-6)
                    return true;
                else
                    return false;
            };

            R.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            R.AgglomerationThreshold = agg;
            R.PrePreCond = MultigridOperator.Mode.SymPart_DiagBlockEquilib;

            return R;
        }

        /// <summary>
        /// A parameter study over <see cref="PiecewiseParabola"/>.
        /// </summary>
        public static IEnumerable<XdgPoisson3Control> PiecewiseParabola_Parameterstudy() {
            List<XdgPoisson3Control> cases = new List<XdgPoisson3Control>();


            foreach (var vmode in new XLaplace_Interface.Mode[] { XLaplace_Interface.Mode.SIP }) {
                foreach (var ppcMode in new MultigridOperator.Mode[] { MultigridOperator.Mode.IdMass, MultigridOperator.Mode.SymPart_DiagBlockEquilib }) {
                    foreach (var _delta in new double[] { 1, 3, 6, 10, 20 }) {
                        double delta = 1.0 - 1.0 / _delta;

                        var C = PiecewiseParabola(3, delta, -1.0, -1.0);
                        C.ViscosityMode = vmode;


                        C.PrePreCond = ppcMode;

                        C.AgglomerationThreshold = 0.0;

                        C.Paramstudy_CaseIdentification.AddRange(new[] {
                                new Tuple<string, object>("Interface_Mode", vmode),
                                new Tuple<string, object>("delta", delta),
                                new Tuple<string, object>("antidelta", 1.0-delta),
                                new Tuple<string, object>("PrePre", ppcMode)
                            });


                        cases.Add(C);
                    }

                }
            }


            return cases;
        }

        /// <summary>
        /// A parameter study over <see cref="Schraeg"/>.
        /// </summary>
        public static IEnumerable<XdgPoisson3Control> Schraeg_Parameterstudy() {
            List<XdgPoisson3Control> cases = new List<XdgPoisson3Control>();

            {
                foreach (var vmode in new XLaplace_Interface.Mode[] { XLaplace_Interface.Mode.SIP }) {
                    foreach (var ppcMode in new MultigridOperator.Mode[] { MultigridOperator.Mode.IdMass, MultigridOperator.Mode.SymPart_DiagBlockEquilib }) {
                        foreach (var mul in new double[] { 0.2, 0.5, 1, 1.5, 2, 2.5 }) {

                            var C = Schraeg();

                            C.ViscosityMode = vmode;

                            C.PrePreCond = ppcMode;

                            C.AgglomerationThreshold = 1.0e-7;
                            C.penalty_multiplyer = mul;

                            C.Paramstudy_CaseIdentification.AddRange(new[] {
                                new Tuple<string, object>("Interface_Mode", vmode),
                                new Tuple<string, object>("mul", mul),
                                new Tuple<string, object>("PrePre", ppcMode)
                            });


                            cases.Add(C);
                        }

                    }
                }
            }

            return cases;
        }


        /// <summary>
        /// A spherical interface in the 3D domain \f$ (-2, 2)^3 \f$.
        /// </summary>
        public static XdgPoisson3Control Ball3D(int pDeg, int Res, LinearSolverCode solverCode = LinearSolverCode.direct_pardiso) {
            //--control "cs:BoSSS.Application.XdgPoisson3.HardCodedControl.Ball3D(2, 5)"
            XdgPoisson3Control R = new XdgPoisson3Control();

            R.ProjectName = "XdgPoisson3/Ball3D";
            R.savetodb = false;

            R.SetDGdegree(pDeg);

            R.GridFunc = delegate () {
                return Grid3D.Cartesian3DGrid(GenericBlas.Linspace(-2, 2, Res), GenericBlas.Linspace(-2, 2, Res), GenericBlas.Linspace(-2, 2, Res));
            };


            // set the level-set
            R.InitialValues_Evaluators.Add("Phi", X => (X[0].Pow2() + X[1].Pow2() + X[2].Pow2()) - 1.0.Pow2());
            //R.InitialValues.Add("Phi", X => X[0] + 0.1);
            R.ExcactSolSupported = false;
            //R.InitialValues.Add("uEx#A", X => X[1]);
            //R.InitialValues.Add("uEx#B", X => X[1]);
            R.InitialValues_Evaluators.Add("rhs#A", X => 1.0);
            R.InitialValues_Evaluators.Add("rhs#B", X => 1.0);


            R.MU_A = -1.0;
            R.MU_B = -1.0;

            R.xLaplaceBCs.g_Diri = ((CommonParamsBnd inp) => 0.0);
            R.xLaplaceBCs.IsDirichlet = (inp => true);

            R.LinearSolver = solverCode.GetConfig();
            R.AgglomerationThreshold = 0.1;
            R.PrePreCond = MultigridOperator.Mode.DiagBlockEquilib;
            R.penalty_multiplyer = 1.1;
            R.ViscosityMode = XLaplace_Interface.Mode.SIP;

            return R;
        }

        /*
        /// <summary>
        /// This is a testrun similar to the one in \public\doc\handbook\apdx-NodeSolverPerformance\XDGPoisson\Part1-Calculations.bws
        /// phase A: 3D sphere, pahse B: [-1,1]^3\sphere
        /// </summary>
        public static XdgPoisson3Control TestOrTreat(int solver = 4) {

            XdgPoisson3Control C = new XdgPoisson3Control();

            switch(solver) {
                case 0:
                C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
                break;
                case 1:
                C.LinearSolver.SolverCode = LinearSolverCode.exp_Kcycle_schwarz;
                break;
                case 2:
                C.LinearSolver.SolverCode = LinearSolverCode.exp_gmres_ILU;
                break;

                case 3:
                C.LinearSolver.SolverCode = LinearSolverCode.exp_gmres_levelpmg;
                break;

                case 4:
                throw new NotImplementedException("deactivated old solver config.");
                //C.LinearSolver.SolverCode = LinearSolverCode.exp_OrthoS_pMG;
                //break;
                default:
                throw new NotImplementedException("guess again");
            }

            C.savetodb = false;
            //C.savetodb = true;
            C.DbPath = @"D:\trash_db";
            //C.DbPath = @"D:\Xdg_Poisson_CondNum";

            int Res = 20;
            int blocksize = 10000;

            C.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-1, +1, Res + 1);
                double[] yNodes = GenericBlas.Linspace(-1, +1, Res + 1);
                //double[] zNodes = GenericBlas.Linspace(-1, +1, Res + 1);
                //int J = (xNodes.Length - 1) * (yNodes.Length - 1) * (zNodes.Length - 1);
                int J = (xNodes.Length - 1) * (yNodes.Length - 1);
                //var grid = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grid.Name = "thisisatestgrid";
                grid.EdgeTagNames.Add(1, "Dirichlet");
                grid.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });
                return grid;
            };

            C.SessionName = String.Format("XDGPoison_solver{0}_blsz{1}_Xdg2lowB", solver, blocksize);
            C.ProjectName = "PoisonTest";
            C.GridPartType = GridPartType.Hilbert;
            C.LinearSolver.TargetBlockSize = blocksize;
            C.SetDGdegree(4);

            //C.PrePreCond = MultigridOperator.Mode.Eye;
            C.LinearSolver.NoOfMultigridLevels = 3;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MaxKrylovDim = 50;
            //C.LinearSolver.pMaxOfCoarseSolver = 1;
            //C.LinearSolver.TargetBlockSize = 79;
            C.ExcactSolSupported = false;
            double radius = 0.7;
            //C.InitialValues_Evaluators.Add("Phi", X => X[0].Pow2() + X[1].Pow2() + X[2].Pow2() - radius.Pow2());
            C.InitialValues_Evaluators.Add("Phi", X => X[0].Pow2() + X[1].Pow2() - radius.Pow2());
            C.MU_A = -1;
            C.MU_B = -1000;
            C.InitialValues_Evaluators.Add("rhs#A", X => 1.0);
            C.InitialValues_Evaluators.Add("rhs#B", X => 1.0);
            C.InitialValues_Evaluators.Add("u#A", X => 0);
            C.InitialValues_Evaluators.Add("u#B", X => 0);
            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.SetDefaultDiriBndCnd = true;  //which means ...
            //C.xLaplaceBCs.g_Diri = ((CommonParamsBnd inp) => 0.0);
            //C.xLaplaceBCs.IsDirichlet = (inp => true);
            // ... but stuff is not serializable, therefore this workaround.
            C.ViscosityMode = XLaplace_Interface.Mode.SIP;

            C.AgglomerationThreshold = 0.1;

            return C;
        }
        */

        /// <summary>
        /// A circular interface within the 2D domain \f$ (-1,1)^2 \f$, with Dirichlet boundary conditions at \f$ x = -1 \f$ and Neumann boundary conditions elsewhere.
        /// </summary>
        public static XdgPoisson3Control CircleNeum(int Resolution = 13) {

            XdgPoisson3Control C = new XdgPoisson3Control();

            C.DbPath = @"D:\\BoSSS-db";
            C.savetodb = false;

            C.ProjectName = "XDGPoisson/Circle";
            C.ProjectDescription = "XDG Poisson Circle";

            //Diffusion Parameter
            C.MU_A = -1.0;
            C.MU_B = -1.0;

            //Funktion on the RHS of the Poisson-Equation
            C.InitialValues_Evaluators.Add("rhs#A", X => Math.Sin((X[0] + 1.0) * 3));
            C.InitialValues_Evaluators.Add("rhs#B", X => Math.Sin((X[0] + 1.0) * 3));


            // Problem Definition
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            double RADIUS = 0.8;
            C.InitialValues_Evaluators.Add("Phi", X => X[0].Pow2() + X[1].Pow2() - (RADIUS).Pow2());

            C.FieldOptions.Add("u", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            //Definition of Grid - Only Dirichlet Boundaries
            C.GridFunc = delegate () {
                return Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, 1, Resolution), GenericBlas.Linspace(-1, 1, Resolution));
            };

            C.xLaplaceBCs.g_Diri = (inp => 0.0);
            C.xLaplaceBCs.g_Neum = delegate (CommonParamsBnd inp) {
                if (Math.Abs(inp.X[1] - 1.0) < 1.0e-8 || Math.Abs(inp.X[1] + 1.0) < 1.0e-8)
                    return 0;
                return Math.Cos(2.0 * 3.0) * (1.0 / 3.0);
            };
            C.xLaplaceBCs.IsDirichlet = (inp => Math.Abs(inp.X[0] + 1.0) <= 1.0e-8);


            // Exact Solution
            C.ExcactSolSupported = true;
            C.FieldOptions.Add("uEx", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.InitialValues_Evaluators.Add("uEx#A", X => Math.Sin((X[0] + 1.0) * 3) * (1.0 / 9.0));
            C.InitialValues_Evaluators.Add("uEx#B", X => Math.Sin((X[0] + 1.0) * 3) * (1.0 / 9.0));

            // Solver Parameters
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            // Discretization Parameters
            C.ViscosityMode = XLaplace_Interface.Mode.SIP;

            return C;

        }

        /*
       
        /// <summary>
        /// A parameter study over <see cref="CircleNeum(int)"/>.
        /// </summary>
        public static IEnumerable<XdgPoisson3Control> CircleNeumParameterStudy() {
            List<XdgPoisson3Control> R = new List<XdgPoisson3Control>();

            var hmfVersions = new[] {
                XQuadFactoryHelper.MomentFittingVariants.Classic,
                XQuadFactoryHelper.MomentFittingVariants.OneStepGauss,
                XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes };
            foreach (var HMFversion in hmfVersions) {
                foreach (int pOff in new int[] { -1, 0, 1, 2 }) {
                    foreach (int res in new int[] { 12, 24, 48, 96, 192 }) {


                        var Ri = CircleNeum(res);

                        Ri.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Resolution", res));

                        Ri.pOff = pOff;
                        Ri.CutCellQuadratureType = HMFversion;

                        R.Add(Ri);

                    }
                }

            }
            return R;
        }
        */
    }

}
