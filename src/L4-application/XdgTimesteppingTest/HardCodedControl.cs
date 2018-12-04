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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using System.Diagnostics;
using System.Linq;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Application.XdgTimesteppingTest {

    public static class HardCodedControl {

        public static XdgTimesteppingTestControl Gerade(
            double angle = 5 * Math.PI / 180.0, int degree = 0, int GridResolutionFactor = 1, double t_offset = 0.0) {
            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/Gerade";
            R.DbPath = null;
            R.savetodb = false;

            


            // DG degree
            // =========

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            R.GridFunc = delegate () {
                var grd = Grid2D.Cartesian2DGrid(
                    GenericBlas.Linspace(-7, 7, 7 * GridResolutionFactor + 1), 
                    GenericBlas.Linspace(-7, 7, 7 * GridResolutionFactor + 1)
                    );
                grd.EdgeTagNames.Add(1, "Inflow");
                grd.DefineEdgeTags(X => (byte)1);
                return grd;
            };

            // level-set over time
            // ===================

            const double S = 0.9;
            R.S = ((double[] X, double t) => S);

            double Nx = Math.Cos(angle);
            double Ny = Math.Sin(angle);
            R.Phi = ((double[] X, double t) => (X[0] - Nx * S * (t + t_offset)) * Nx + (X[1] - Ny * S * (t + t_offset)) * Ny);

            // exact solution
            // ==============

            R.uA_Ex = ((X, t) => 0.8);
            R.uB_Ex = ((X, t) => 1.66);

            // Initial Values
            // ==============
            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("Vx", X => Nx * S);
            R.InitialValues_Evaluators.Add("Vy", X => Ny * S);
            R.InitialValues_Evaluators.Add("u#A", X => R.uA_Ex(X, 0.0));
            R.InitialValues_Evaluators.Add("u#B", X => R.uB_Ex(X, 0.0));

            // Boundary values
            // ===============

            R.AddBoundaryValue("Inflow", "u", (X, t) => 0.8);

            // Timestepping config
            // ===================

            R.NoOfTimesteps = 1;
            R.Endtime = 2;
            R.dtFixed = R.Endtime / R.NoOfTimesteps;


            R.TimeSteppingScheme = TimeSteppingScheme.ExplicitEuler;
            R.InterfaceMode = InterfaceMode.MovingInterface;
            R.AgglomerationThreshold = 0.1;

            // return
            // ======

            return R;
        }


        public static XdgTimesteppingTestControl Burgers(
            double angle = 0.0,
            int degree = 2,
            int NoOfTimesteps = 80,
            InterfaceMode tsm = InterfaceMode.MovingInterface,
            int GridResolutionFactor = 2) {

            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/Burgers";
            R.savetodb = false;
            R.DbPath = null;

            // DG config
            // =========

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            R.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-2, 2, 7 * GridResolutionFactor + 1);
                double[] Ynodes = GenericBlas.Linspace(-2, 2, 7 * GridResolutionFactor + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Neumann");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if (Math.Abs(X[0] - (-1)) <= 1.0e-8 || Math.Abs(X[0] - (+1)) <= 1.0e-8)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            // exact solution
            // ==============

            double TimeOffset = 0.1;

            
            double Nx = Math.Cos(angle);
            double Ny = Math.Sin(angle);

            R.BurgersDirection = new Platform.LinAlg.Vector(Nx, Ny);

            const double S = 0.5 * (2 + 1);
            R.S = ((double[] X, double t) => S);
            R.Phi = ((double[] X, double t) => (X[0] - Nx * S * (t + TimeOffset)) * Nx + (X[1] - Ny * S * (t + TimeOffset)) * Ny);

            Func<double[], double> projCoord = X => X[0] * Nx + X[1] * Ny;
       
            Func<double, double, double> uAEx = delegate (double t, double xi) {
                return 2.0;
            };
            
            R.uA_Ex = ((X, t) => uAEx(t, projCoord(X)));
            R.uB_Ex = ((X, t) => 1.0);

            R.Eq = Equation.Burgers;

            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("u#A", X => R.uA_Ex(X, 0.0));
            R.InitialValues_Evaluators.Add("u#B", X => R.uB_Ex(X, 0.0));


            // anderes zeugs
            // =============

            R.TimeSteppingScheme = TimeSteppingScheme.RK4;
            R.InterfaceMode = tsm;

            R.Endtime = 0.2;
            R.NoOfTimesteps = NoOfTimesteps;
            R.dtFixed = R.Endtime / R.NoOfTimesteps;

            R.AgglomerationThreshold = 0.1;

            // return
            // ======

            return R;
        }

  
    }


    public static class HardCodedControl2 {

        public static XdgTimesteppingTestControl Gerade(
            double angle = 5 * Math.PI / 180.0, int degree = 0, int GridResolutionFactor = 1) {
            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/Gerade";
            R.DbPath = @"\\fdyprime\userspace\kummer\BoSSS-db-XNSE";
            R.savetodb = false;

            // DG degree
            // =========

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            R.GridFunc = delegate () {
                var grd = Grid2D.Cartesian2DGrid(
                    GenericBlas.Linspace(-7, 7, 7 * GridResolutionFactor + 1),
                    GenericBlas.Linspace(-7, 7, 7 * GridResolutionFactor + 1));
                grd.EdgeTagNames.Add(1, "Inflow");
                grd.DefineEdgeTags(X => (byte)1);
                return grd;
            };

            // level-set over time
            // ===================

            const double S = 0.9;
            R.S = ((double[] X, double t) => S);

            double Nx = Math.Cos(angle);
            double Ny = Math.Sin(angle);
            double t_offset = 0.0;
            R.Phi = ((double[] X, double t) => (X[0] - Nx * S * (t + t_offset)) * Nx + (X[1] - Ny * S * (t + t_offset)) * Ny);

            // exact solution
            // ==============

            R.uA_Ex = ((X, t) => 0.8);
            R.uB_Ex = ((X, t) => 1.66);

            // Initial Values
            // ==============
            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("Vx", X => Nx * S);
            R.InitialValues_Evaluators.Add("Vy", X => Ny * S);
            R.InitialValues_Evaluators.Add("u#A", X => R.uA_Ex(X, 0.0));
            R.InitialValues_Evaluators.Add("u#B", X => R.uB_Ex(X, 0.0));

            // Boundary values
            // ===============

            R.AddBoundaryValue("Inflow", "u", (X, t) => 0.8);

            // Timestepping config
            // ===================

            R.NoOfTimesteps = 1;
            R.Endtime = 0.4;
            R.dtFixed = R.Endtime / R.NoOfTimesteps;


            R.TimeSteppingScheme = TimeSteppingScheme.ExplicitEuler;
            R.InterfaceMode = InterfaceMode.MovingInterface;
            R.AgglomerationThreshold = 0.1;

            // return
            // ======

            return R;
        }

        public static XdgTimesteppingTestControl[] GeradePStudy() {

            List<XdgTimesteppingTestControl> All = new List<XdgTimesteppingTestControl>();

            InterfaceMode[] mode = new[] { InterfaceMode.MovingInterface };
            TimeSteppingScheme[] schemes = new[] { TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.BDF3 };
            int[] DGdegree = new[] { 2 };
            int[] NoOfTimeSteps = new int[] { 8, 16, 32, 64, 128, 256, 512, 1024 };
            //int[] NoOfTimeSteps = new int[] { 1024, 2048, 4096 };
            //int[] GridResolutionFactorS = new int[] { 1, 2, 4, 8 };
            int[] GridResolutionFactorS = new int[] { 1, 2 };
            double[] AgglomFactors = new double[] { 0.1 };

            foreach (var tm in mode) {
                foreach (int p in DGdegree) {
                    foreach (TimeSteppingScheme scheme in schemes) {
                        foreach (int N in NoOfTimeSteps) {
                            foreach (double alpha in AgglomFactors) {
                                foreach (int grf in GridResolutionFactorS) {


                                    var C = Gerade(degree: p, GridResolutionFactor: grf, angle: 45 * Math.PI / 180);
                                    C.TimeSteppingScheme = scheme;
                                    C.NoOfTimesteps = N;
                                    C.dtFixed = C.Endtime / N;
                                    C.InterfaceMode = tm;
                                    C.savetodb = true;
                                    C.AgglomerationThreshold = alpha;

                                    C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string, object>("TimesteppingMode", tm),
                                    new Tuple<string, object>("scheme", scheme),
                                    new Tuple<string, object>("NoOfTimesteps", N),
                                    new Tuple<string, object>("DGdegree", p),
                                    new Tuple<string, object>("AgglomerationThreshold", C.AgglomerationThreshold),
                                    new Tuple<string, object>("GridResolutionFactor", grf),
                                    new Tuple<string, object>("MoveOrNotMove", true)
                                    };

                                    All.Add(C);
                                }
                            }
                        }
                    }
                }
            }

            /*

            // vergleich: Statisches Gitter

            foreach (int p in DGdegree) {
                foreach (TimeSteppingScheme scheme in schemes) {
                    foreach (int N in NoOfTimeSteps) {
                        foreach (int grf in GridResolutionFactorS) {
                            if (grf >= 8) {
                                if (N <= 4)
                                    continue;
                            }
                            if (grf >= 4) {
                                if (N <= 2)
                                    continue;
                            }

                            var C = Gerade(degree: p, GridResolutionFactor: grf);
                            C.TimeSteppingScheme = scheme;
                            C.NoOfTimesteps = N;
                            C.dtFixed = C.Endtime / N;
                            C.savetodb = true;
                            C.TimeSteppingMode = TimesteppingMode.Splitting;

                            var orgPhi = C.Phi;
                            C.Phi = (X, t) => orgPhi(X, 0.0);

                            C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string, object>("TimesteppingMode", TimesteppingMode.Splitting),
                                    new Tuple<string, object>("scheme", scheme),
                                    new Tuple<string, object>("NoOfTimesteps", N),
                                    new Tuple<string, object>("DGdegree", p),
                                    new Tuple<string, object>("AgglomerationThreshold", C.AgglomerationThreshold),
                                    new Tuple<string, object>("GridResolutionFactor", grf),
                                    new Tuple<string, object>("MoveOrNotMove", false)
                                    };

                            All.Add(C);
                        }
                    }
                }
            }

            */



            return All.ToArray();
        }

        public static XdgTimesteppingTestControl[] GeradePStudyAgglom() {

            List<XdgTimesteppingTestControl> All = new List<XdgTimesteppingTestControl>();

            InterfaceMode[] mode = new[] { InterfaceMode.MovingInterface };
            TimeSteppingScheme[] schemes = new[] { TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3, TimeSteppingScheme.BDF4 };
            int[] DGdegree = new[] { 2, 3, 4 };
            int[] NoOfTimeSteps = new int[] { 2048, 4096 };
            int[] GridResolutionFactorS = new int[] { 1, 2 };
            double[] AgglomFactors = new double[] { 0.1 };

            foreach (var tm in mode) {
                foreach (int p in DGdegree) {
                    foreach (TimeSteppingScheme scheme in schemes) {
                        foreach (int N in NoOfTimeSteps) {
                            foreach (double alpha in AgglomFactors) {
                                foreach (int grf in GridResolutionFactorS) {


                                    var C = Gerade(degree: p, GridResolutionFactor: grf, angle: 45 * Math.PI / 180);
                                    C.TimeSteppingScheme = scheme;
                                    C.NoOfTimesteps = N;
                                    C.dtFixed = C.Endtime / N;
                                    C.InterfaceMode = tm;
                                    C.savetodb = true;
                                    C.AgglomerationThreshold = alpha;
                                    C.ProjectName += "/AgglomStudy";

                                    C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string, object>("TimesteppingMode", tm),
                                    new Tuple<string, object>("scheme", scheme),
                                    new Tuple<string, object>("NoOfTimesteps", N),
                                    new Tuple<string, object>("DGdegree", p),
                                    new Tuple<string, object>("AgglomerationThreshold", C.AgglomerationThreshold),
                                    new Tuple<string, object>("GridResolutionFactor", grf),
                                    new Tuple<string, object>("MoveOrNotMove", true)
                                    };

                                    All.Add(C);
                                }
                            }
                        }
                    }
                }
            }


            return All.ToArray();
        }

        public static XdgTimesteppingTestControl AggDbg() {
            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/AgglomDebuch";
            R.savetodb = false;

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 3,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            R.GridFunc = delegate () {
                var grd = Grid2D.Cartesian2DGrid(new double[] { -4, -2, 0, 2, 4 }, GenericBlas.Linspace(-1, 1, 2));
                grd.AddPredefinedPartitioning("For2", (double[] X) => (Math.Sign(X[0]) + 1) / 2);
                return grd;
            };

            R.GridPartType = GridPartType.Predefined;
            R.GridPartOptions = "For2";

            const double S = 0.3;
            R.S = ((double[] X, double t) => S);
            R.Phi = ((double[] X, double t) => X[0] - S * t);

            // set the level-set
            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("Vx", X => 0.3);
            R.InitialValues_Evaluators.Add("Vy", X => 0.0);
            R.InitialValues_Evaluators.Add("u#A", X => 0.8);
            R.InitialValues_Evaluators.Add("u#B", X => 1.66);

            //R.InitialValues_Evaluators.Add("u#A", X => X[0] < 0 ? 0.8 : 1.66);
            //R.InitialValues_Evaluators.Add("u#B", X => X[0] < 0 ? 0.8 : 1.66);


            R.dtFixed = 0.18;


            return R;
        }

        public static XdgTimesteppingTestControl Rarefaction(int degree = 2,
            int NoOfTimesteps = 20,
            InterfaceMode tsm = InterfaceMode.MovingInterface,
            int GridResolutionFactor = 2) {

            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/Rarefaction";
            R.savetodb = false;
            //R.DbPath = @"\\fdyprime\userspace\kummer\BoSSS-db-XNSE";
            R.DbPath = null;

            // DG config
            // =========

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            R.GridFunc = delegate () {
                double[] nodes = GenericBlas.Linspace(-2.7, 2.7, 18 * GridResolutionFactor + 1);
                BoundingBox cutOut = new BoundingBox(new double[] { -0.3, -0.3 }, new double[] { 0.3, 0.3 });
                var grd = Grid2D.Cartesian2DGrid(nodes, nodes, CutOuts: cutOut);

                grd.EdgeTagNames.Add(1, "Inflow");
                grd.EdgeTagNames.Add(2, "Outflow");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret = 2;
                    if (Math.Abs(X[0]) <= 0.31 && Math.Abs(X[1]) <= 0.31)
                        ret = 1;
                    return ret;
                });

                return grd;
            };

            // exact solution
            // ==============
            R.uA_Ex = ((X, t) => 3.0 / Math.Sqrt(X[0].Pow2() + X[1].Pow2()));
            R.uB_Ex = ((X, t) => 1.0 / Math.Sqrt(X[0].Pow2() + X[1].Pow2()));

            // boundary condition
            // ==================

            R.AddBoundaryValue("Inflow", "u", R.uA_Ex);
            R.AddBoundaryValue("Outflow");

            // Initial values
            // ==============

            const double S = 1;
            R.S = ((double[] X, double t) => S);
            R.Phi = ((double[] X, double t) => X[0].Pow2() + X[1].Pow2() - (1.0 + t).Pow2());

            R.CircleRadius = t => (1.0 + t);
            R.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.ExactCircle;

            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("Vx", X => X[0] / Math.Sqrt(X[0].Pow2() + X[1].Pow2()));
            R.InitialValues_Evaluators.Add("Vy", X => X[1] / Math.Sqrt(X[0].Pow2() + X[1].Pow2()));
            R.InitialValues_Evaluators.Add("u#A", X => R.uA_Ex(X, 0.0));
            R.InitialValues_Evaluators.Add("u#B", X => R.uB_Ex(X, 0.0));

            // restart
            // =======

            //R.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(new Guid("aff36e92-1546-4fdf-a7bc-fbeff1e67f49"), 58); 
            //R.InitialValues_Evaluators.Clear();

            // anderes zeugs
            // =============

            R.TimeSteppingScheme = TimeSteppingScheme.BDF4;
            R.InterfaceMode = tsm;

            //R.Endtime = 0.05;
            R.Endtime = 0.4;
            R.NoOfTimesteps = NoOfTimesteps;
            R.dtFixed = R.Endtime / R.NoOfTimesteps;

            R.AgglomerationThreshold = 0.1;

            // return
            // ======

            return R;
        }

        public static XdgTimesteppingTestControl[] RarefactionPStudy() {

            List<XdgTimesteppingTestControl> All = new List<XdgTimesteppingTestControl>();

            InterfaceMode[] mode = new[] { InterfaceMode.Splitting, InterfaceMode.MovingInterface };
            TimeSteppingScheme[] schemes = new[] { TimeSteppingScheme.ImplicitEuler };
            int[] DGdegree = new[] { 2 };
            int[] NoOfTimeSteps = new int[] { 1, 2, 4, 8, 16, 32, 64, 128, 256 };
            int[] GridResolutionFactorS = new int[] { 16 };


            foreach (var tm in mode) {
                foreach (int p in DGdegree) {
                    foreach (TimeSteppingScheme scheme in schemes) {
                        foreach (int grf in GridResolutionFactorS) {
                            foreach (int N in NoOfTimeSteps) {

                                //if (p <= 2 && N <= 128)
                                //    continue;

                                var C = Rarefaction(degree: p, NoOfTimesteps: N, tsm: tm, GridResolutionFactor: grf);
                                C.TimeSteppingScheme = scheme;
                                C.savetodb = true;

                                C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                        new Tuple<string, object>("TimesteppingMode", tm),
                                        new Tuple<string, object>("scheme", scheme),
                                        new Tuple<string, object>("NoOfTimesteps", N),
                                        new Tuple<string, object>("DGdegree", p),
                                        new Tuple<string, object>("AgglomerationThreshold",C.AgglomerationThreshold),
                                        new Tuple<string, object>("GridResolutionFactor", grf),
                                        new Tuple<string, object>("MoveOrNotMove", true)
                                        };

                                All.Add(C);
                            }
                        }
                    }
                }
            }
            return All.ToArray();
        }

        public static XdgTimesteppingTestControl Heat1D(int degree = 1,
            int NoOfTimesteps = 8,
            InterfaceMode tsm = InterfaceMode.MovingInterface,
            int GridResolutionFactor = 1) {

            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/Heat1D";
            R.savetodb = false;
            //R.DbPath = @"\\fdyprime\userspace\kummer\BoSSS-db-XNSE";
            R.DbPath = null;

            // DG config
            // =========

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            R.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-1, 1, 18 * GridResolutionFactor + 1);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, 3);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Neumann");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if (Math.Abs(X[0] - (-1)) <= 1.0e-8 || Math.Abs(X[0] - (+1)) <= 1.0e-8)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            // exact solution
            // ==============

            double muA = 1;
            double muB = 10;
            Func<double, double> hA = t => 1 - t - t.Pow2();
            Func<double, double> hB = t => -.3125000000 * Math.Sin(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) / Math.Sin(.2513274123 + .2513274123 * t);
            Func<double, double> oA = t => (0.6250000000e-1 * (16.0 * Math.Cos(.7853981635 + .7853981635 * t) * t.Pow2() * Math.Sin(.2513274123 + .2513274123 * t) - 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) * t.Pow2() + 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t) * t.Pow2() + 16.0 * Math.Cos(.7853981635 + .7853981635 * t) * t * Math.Sin(.2513274123 + .2513274123 * t) - 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) * t + 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t) * t - 16.0 * Math.Cos(.7853981635 + .7853981635 * t) * Math.Sin(.2513274123 + .2513274123 * t) + 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) - 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t))) / Math.Sin(.2513274123 + .2513274123 * t);
            Func<double, double> oB = t => .2528178107 * Math.Sin(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) / Math.Sin(.2513274123 + .2513274123 * t);


            Func<double, double, double> uA_Ex = (t, x) => oA(t) + hA(t) * Math.Cos(x * Math.PI * 0.5 * 10.0 / 8.0);
            Func<double, double, double> uB_Ex = (t, x) => oB(t) + hB(t) * Math.Cos(x * Math.PI / 5);
            Func<double, double, double> rhsA = (t, x) => (0.6250000000e-1 * (-11.30973355 * Math.Sin(.7853981635 + .7853981635 * t) * t.Pow2() * Math.Sin(.2513274123 + .2513274123 * t) + 32.0 * Math.Cos(.7853981635 + .7853981635 * t) * t * Math.Sin(.2513274123 + .2513274123 * t) + 0.9424777962e-1 * Math.Cos(.7853981635 + .7853981635 * t) * t.Pow2() * Math.Cos(.2513274123 + .2513274123 * t) - 10.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) * t + 3.177002308 * Math.Cos(.7853981635 + .7853981635 * t) * t.Pow2() + 8.090169943 * Math.Sin(.7853981635 + .7853981635 * t) * t - 11.30973355 * Math.Sin(.7853981635 + .7853981635 * t) * t * Math.Sin(.2513274123 + .2513274123 * t) + 16.0 * Math.Cos(.7853981635 + .7853981635 * t) * Math.Sin(.2513274123 + .2513274123 * t) + 0.9424777962e-1 * Math.Cos(.7853981635 + .7853981635 * t) * t * Math.Cos(.2513274123 + .2513274123 * t) - 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) + 3.177002308 * Math.Cos(.7853981635 + .7853981635 * t) * t + 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t) + 11.30973355 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Sin(.2513274123 + .2513274123 * t) - 0.9424777962e-1 * Math.Cos(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) - 3.177002308 * Math.Cos(.7853981635 + .7853981635 * t))) / Math.Sin(.2513274123 + .2513274123 * t) - (0.1570796327e-1 * (16.0 * Math.Cos(.7853981635 + .7853981635 * t) * t.Pow2() * Math.Sin(.2513274123 + .2513274123 * t) - 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) * t.Pow2() + 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t) * t.Pow2() + 16.0 * Math.Cos(.7853981635 + .7853981635 * t) * t * Math.Sin(.2513274123 + .2513274123 * t) - 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) * t + 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t) * t - 16.0 * Math.Cos(.7853981635 + .7853981635 * t) * Math.Sin(.2513274123 + .2513274123 * t) + 5.0 * Math.Sin(.7853981635 + .7853981635 * t) * Math.Cos(.2513274123 + .2513274123 * t) - 4.045084972 * Math.Sin(.7853981635 + .7853981635 * t))) * Math.Cos(.2513274123 + .2513274123 * t) / Math.Sin(.2513274123 + .2513274123 * t).Pow2() + (-2.0 * t - 1.0) * Math.Cos(1.963495409 * x) + (3.855314220 * (-1.0 * t.Pow2() - 1.0 * t + 1.0)) * Math.Cos(1.963495409 * x);
            Func<double, double, double> rhsB = (t, x) => -0.6354004615e-1 * Math.Sin(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) * Math.Cos(.2513274123 + .2513274123 * t) / Math.Sin(.2513274123 + .2513274123 * t).Pow2() + .1985626442 * Math.Cos(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) / Math.Sin(.2513274123 + .2513274123 * t) + .2528178107 * Math.Sin(.7853981635 + .7853981635 * t) * (2.0 * t + 1.0) / Math.Sin(.2513274123 + .2513274123 * t) + 0.7853981635e-1 * Math.Sin(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) * Math.Cos(.6283185308 * x) * Math.Cos(.2513274123 + .2513274123 * t) / Math.Sin(.2513274123 + .2513274123 * t).Pow2() - .2454369261 * Math.Cos(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) * Math.Cos(.6283185308 * x) / Math.Sin(.2513274123 + .2513274123 * t) - .3125000000 * Math.Sin(.7853981635 + .7853981635 * t) * (2.0 * t + 1.0) * Math.Cos(.6283185308 * x) / Math.Sin(.2513274123 + .2513274123 * t) - 1.233700550 * Math.Sin(.7853981635 + .7853981635 * t) * (t.Pow2() + t - 1.0) * Math.Cos(.6283185308 * x) / Math.Sin(.2513274123 + .2513274123 * t);

            R.uA_Ex = ((X, t) => uA_Ex(t, X[0]));
            R.uB_Ex = ((X, t) => uB_Ex(t, X[0]));

            R.rhsA = ((X, t) => rhsA(t, X[0]));
            R.rhsB = ((X, t) => rhsB(t, X[0]));

            R.Eq = Equation.HeatEq;
            R.muA = muA;
            R.muB = muB;

            const double S = 0.4;
            R.S = ((double[] X, double t) => S);
            R.Phi = ((double[] X, double t) => X[0].Pow2() - (0.4 + t * S).Pow2());

            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("u#A", X => R.uA_Ex(X, 0));
            R.InitialValues_Evaluators.Add("u#B", X => R.uB_Ex(X, 0));

            // restart
            // =======

            //R.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(new Guid("aff36e92-1546-4fdf-a7bc-fbeff1e67f49"), 58); 
            //R.InitialValues_Evaluators.Clear();

            // anderes zeugs
            // =============

            R.TimeSteppingScheme = TimeSteppingScheme.BDF4;
            R.InterfaceMode = tsm;

            R.Endtime = 0.4;
            R.NoOfTimesteps = NoOfTimesteps;
            R.dtFixed = R.Endtime / R.NoOfTimesteps;

            R.AgglomerationThreshold = 0.1;

            // return
            // ======

            return R;
        }

        public static XdgTimesteppingTestControl[] Heat1DPStudy() {

            List<XdgTimesteppingTestControl> All = new List<XdgTimesteppingTestControl>();

            InterfaceMode[] mode = new[] { InterfaceMode.Splitting };
            TimeSteppingScheme[] schemes = new[] { TimeSteppingScheme.BDF3, TimeSteppingScheme.BDF4 };
            int[] DGdegree = new[] { 1 };
            //int[] NoOfTimeSteps = new int[] { 8, 16, };
            int[] NoOfTimeSteps = new int[] { 4096, 4096 * 2 };
            int[] GridResolutionFactorS = new int[] { 1, 2, 4, 8, 16, 32 };

            foreach (var tm in mode) {
                foreach (int p in DGdegree) {
                    foreach (TimeSteppingScheme scheme in schemes) {
                        foreach (int N in NoOfTimeSteps) {
                            foreach (int grf in GridResolutionFactorS) {

                                var C = Heat1D(degree: p, GridResolutionFactor: grf, tsm: tm);
                                C.TimeSteppingScheme = scheme;
                                C.NoOfTimesteps = N;
                                C.dtFixed = C.Endtime / N;
                                C.InterfaceMode = tm;
                                C.savetodb = true;

                                C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string, object>("TimesteppingMode", tm),
                                    new Tuple<string, object>("scheme", scheme),
                                    new Tuple<string, object>("NoOfTimesteps", N),
                                    new Tuple<string, object>("DGdegree", p),
                                    new Tuple<string, object>("AgglomerationThreshold", C.AgglomerationThreshold),
                                    new Tuple<string, object>("GridResolutionFactor", grf),
                                    new Tuple<string, object>("MoveOrNotMove", true)
                                    };

                                All.Add(C);
                            }
                        }
                    }
                }
            }

            return All.ToArray();
        }

        public static XdgTimesteppingTestControl Burgers(int degree = 1,
            int NoOfTimesteps = 40,
            InterfaceMode tsm = InterfaceMode.MovingInterface,
            int GridResolutionFactor = 2) {

            XdgTimesteppingTestControl R = new XdgTimesteppingTestControl();

            R.ProjectName = "XdgMassMatrixEvolution/Burgers";
            R.savetodb = false;
            R.DbPath = null;  //@"\\fdyprime\userspace\kummer\BoSSS-db-XNSE";

            // DG config
            // =========

            R.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("u", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vx", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            R.FieldOptions.Add("Vy", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            R.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-2, 2, 7 * GridResolutionFactor + 1);
                double[] Ynodes = GenericBlas.Linspace(-2, 2, 7 * GridResolutionFactor + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Neumann");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if (Math.Abs(X[0] - (-1)) <= 1.0e-8 || Math.Abs(X[0] - (+1)) <= 1.0e-8)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            // exact solution
            // ==============

            double TimeOffset = 0.1;

            double angle = 45 * Math.PI / 180.0;
            //double angle = 0.0;
            double Nx = Math.Cos(angle);
            double Ny = Math.Sin(angle);

            R.BurgersDirection = new Platform.LinAlg.Vector(Nx, Ny);

            const double S = 0.5 * (2 + 1);
            R.S = ((double[] X, double t) => S);
            R.Phi = ((double[] X, double t) => (X[0] - Nx * S * (t + TimeOffset)) * Nx + (X[1] - Ny * S * (t + TimeOffset)) * Ny);

            Func<double[], double> projCoord = X => X[0] * Nx + X[1] * Ny;

            // initial value for x < 0
            Func<double, double> uA0 = x => Math.Exp(-x.Pow2()) + 1.0;

            // We will need Newtons alg. to find the origin of a characteristic;
            // for a decent initial value of Newton, we sample the initial value at certain points 
            double[] xTrial = GenericBlas.Linspace(-10, 10, 1000);
            double[] u0Tril = xTrial.Select(x => uA0(x)).ToArray();

            
            Func<double, double, double> uAEx = delegate (double t, double xi) {
                double ret = 0;

                t += TimeOffset;

                if (t < 0)
                    throw new ArgumentOutOfRangeException();
                if (t == 0)
                    return uA0(xi);

                // find low ang high bound for Newton
                int I = xTrial.Length;
                int iMinDist_lo = -1;
                double MinDist_lo = double.MaxValue;
                int iMinDist_hi = -1;
                double MinDist_hi = double.MaxValue;

                for (int i = 0; i < I; i++) {
                    double x0 = xTrial[i];
                    double u0Atx0 = u0Tril[i];

                    double x1 = x0 + u0Atx0 * t; // x1 is position of the characteristic, which originates from x0,  at time t
                    double dist = Math.Abs(x1 - xi);
                    if (dist < MinDist_lo && x1 < xi) {
                        iMinDist_lo = i;
                        MinDist_lo = dist;
                    }
                    if (dist < MinDist_hi && x1 > xi) {
                        iMinDist_hi = i;
                        MinDist_hi = dist;
                    }
                }
                double xi_lo = xTrial[iMinDist_lo];
                double xi_hi = xTrial[iMinDist_hi];



                Func<double,double> func = (xi0 => (uA0(xi0) * t + xi0 - xi));
                Func<double, double> dfunc = (xi0 => (1.0 - 2.0 * t * Math.Exp(-(xi0.Pow2()))));
                double xi0i = rtsave(func, dfunc, xi_lo, xi_hi, 1e-12);

                // Probe
                //if (!converged)
                //    throw new ArithmeticException();

                // return value at origin of characteristic
                ret = uA0(xi0i);

                return ret;
            };

            {
                // test for uAEx
                double x0 = -0.5;
                double u_at_x0 = uA0(x0);
                double dt = -0.02;
                double x = x0 + u_at_x0 * (dt + TimeOffset); // characteristic through x0, at time dt

                double u_at = uAEx(dt, x);
                if (Math.Abs(u_at - u_at_x0) > 1.0e-8)
                    throw new ArithmeticException();
            }
            
            /*
            Func<double, double, double> uAEx = delegate (double t, double xi) {
                return 2.0;
            };
            */

            R.uA_Ex = ((X, t) => uAEx(t, projCoord(X)));
            R.uB_Ex = ((X, t) => 1.0);

            R.Eq = Equation.Burgers;

            R.InitialValues_Evaluators.Add("Phi", X => R.Phi(X, 0.0));
            R.InitialValues_Evaluators.Add("u#A", X => R.uA_Ex(X, 0.0));
            R.InitialValues_Evaluators.Add("u#B", X => R.uB_Ex(X, 0.0));


            // anderes zeugs
            // =============

            R.TimeSteppingScheme = TimeSteppingScheme.RK4;
            R.MultiStepInit = false;
            R.InterfaceMode = tsm;

            R.Endtime = 0.11;
            R.NoOfTimesteps = NoOfTimesteps;
            R.dtFixed = R.Endtime / R.NoOfTimesteps;

            R.AgglomerationThreshold = 0.1;

            // return
            // ======

            return R;
        }

        /// <summary>
        /// Safeguard-Newton, from Numerical Recipes, p460.
        /// </summary>
        static double rtsave(Func<double, double> func, Func<double, double> dfunc, double x1, double x2, double accuracy) {
            const int MAXIT = 10000;
            double fl = func(x1);
            double fh = func(x2);
            if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
                throw new ArithmeticException("Root must be bracketed!");
            if (fl == 0)
                return x1;
            if (fh == 0)
                return x2;
            double xl, xh;
            if (fl < 0) {
                xl = x1;
                xh = x2;
            } else {
                xh = x1;
                xl = x2;
            }
            double rts = 0.5 * (x1 + x2);
            double dxold = Math.Abs(x2 - x1);
            double dx = dxold;
            double f = func(rts);
            double df = dfunc(rts);
            for (int j = 0; j < MAXIT; j++) {
                if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0) || (Math.Abs(2.0 * f) > Math.Abs(dxold * df))) {
                    dxold = dx;
                    dx = 0.5 * (xh - xl);
                    rts = xl + dx;
                } else {
                    dxold = dx;
                    dx = f / df;
                    double temp = rts;
                    rts -= dx;
                    if (temp == rts)
                        return rts;
                }


                if (Math.Abs(dx) < accuracy)
                    return rts;

                f = func(rts);
                df = dfunc(rts);

                if (f < 0.0)
                    xl = rts;
                else
                    xh = rts;
            }

            throw new ArithmeticException("Max number of iter exceeded.");
        }

        public static XdgTimesteppingTestControl[] BurgersPStudy() {

            List<XdgTimesteppingTestControl> All = new List<XdgTimesteppingTestControl>();

            InterfaceMode[] mode = new[] { InterfaceMode.MovingInterface, InterfaceMode.Splitting };
            TimeSteppingScheme[] schemes = new[] { TimeSteppingScheme.RK4 };
            int[] DGdegree = new[] { 0, 1, 2 };
            int[] NoOfTimeSteps = new int[] { 5, 10, 20, 40, 80, 160, 320, 640 };
            int[] GridResolutionFactorS = new int[] { 1, 2, 4, 8 };

            foreach (var tm in mode) {
                foreach (int p in DGdegree) {
                    foreach (TimeSteppingScheme scheme in schemes) {
                        foreach (int N in NoOfTimeSteps) {
                            foreach (int grf in GridResolutionFactorS) {

                                var C = Burgers(degree: p, GridResolutionFactor: grf, tsm: tm);

                                C.TimeSteppingScheme = scheme;
                                C.NoOfTimesteps = N;
                                C.dtFixed = C.Endtime / N;
                                C.InterfaceMode = tm;
                                C.savetodb = true;

                                C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string, object>("TimesteppingMode", tm),
                                    new Tuple<string, object>("scheme", scheme),
                                    new Tuple<string, object>("NoOfTimesteps", N),
                                    new Tuple<string, object>("DGdegree", p),
                                    new Tuple<string, object>("AgglomerationThreshold", C.AgglomerationThreshold),
                                    new Tuple<string, object>("GridResolutionFactor", grf),
                                    new Tuple<string, object>("MoveOrNotMove", true)
                                    };

                                All.Add(C);
                            }
                        }
                    }
                }
            }

            return All.ToArray();
        }
    }
}

