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
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using ilPSP;
using System.Diagnostics;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Tests whether the XNSE solver (<see cref="XNSE_SolverMain"/>) also works MPI-parallel for 
    /// non-trivial cases.
    /// </summary>
    [TestFixture]
    public static class XNSE_Solver_MPItest {

        /// <summary>
        /// MPI initialization.
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool MpiInit;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);
        }

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }


        [Test]
        static public void ParallelRisingDroplet() {
            var C = RisingBubble();

            using (var solver = new XNSE_SolverMain()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }





        /// <summary>
        /// 
        /// </summary>
        static void Main(string[] args) {
            SetUp();
            ParallelRisingDroplet();
            /*
            int MPIrank, MPIsize;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MPIsize);

            MPI_Request[] rq = new MPI_Request[2];
            MPI_Status[] st = new MPI_Status[2];

            unsafe
            {
                double[] data, rcv_data = new double[10];
                int Dest, Src;
                int Tag, SrcTag;
                int RcvSz;
                if (MPIrank == 0) {
                    data = new double[] { 1, 2, 3 };
                    Dest = 1;
                    Src = 1;
                    Tag = 22;
                    SrcTag = 19;
                    RcvSz = 1;
                } else {
                    data = new double[] { 66 };
                    Dest = 0;
                    Src = 0;
                    Tag = 19;
                    SrcTag = 22;
                    RcvSz = 3;
                }

                Debugger.Launch();

                fixed (double* pData = data, pRcvdata = rcv_data) {
                    csMPI.Raw.Irecv((IntPtr)pRcvdata, RcvSz, csMPI.Raw._DATATYPE.DOUBLE, Src, SrcTag, csMPI.Raw._COMM.WORLD, out rq[0]);
                    csMPI.Raw.Issend((IntPtr)pData, data.Length, csMPI.Raw._DATATYPE.DOUBLE, Dest, Tag, csMPI.Raw._COMM.WORLD, out rq[1]);
                }


                csMPI.Raw.Waitall(2, rq, st);

                //Debugger.Launch();
            }
            */
            TestFixtureTearDown();
        }


        /// <summary>
        /// Configuration which performs three timesteps of the rising droplet;
        /// uses a pre-defined partitioning for 4 processors.
        /// </summary>
        public static XNSE_Control RisingBubble(int p = 3, int kelem = 20, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath =  @" C.NoOfTimesteps = 10;";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;//@"\\fdyprime\userspace\nietz\databases\big_test";
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Bubble";
            C.ProjectDescription = "rising bubble";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            /*C.FieldOptions.Add("FilteredVelocityX", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("FilteredVelocityY", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });*/
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Testcase 1");
            C.PhysicalParameters.rho_A = 100;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 10;
            C.PhysicalParameters.Sigma = 24.5; //2.0*




            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 2.0;


            //int kelem = 40;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1); //ohne 2*
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                //grd.EdgeTagNames.Add(3, "wall_left");
                //grd.EdgeTagNames.Add(4, "wall_right");
                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                grd.AddPredefinedPartitioning("ZwoProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if(x < 0.5)
                        rank = 0;
                    else
                        rank = 1;

                    return rank;
                });

                grd.AddPredefinedPartitioning("VierProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if(x < 0.35)
                        rank = 0;
                    else if(x < 0.5)
                        rank = 1;
                    else if(x < 0.75)
                        rank = 2;
                    else
                        rank = 3;

                    return rank;
                });


                return grd;
            };

            int MpiSize;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MpiSize);
            switch(MpiSize) {
                case 4:
                C.GridPartType = GridPartType.Predefined;
                C.GridPartOptions = "VierProcSplit";
                break;

                case 2:
                C.GridPartType = GridPartType.Predefined;
                C.GridPartOptions = "ZwoProcSplit";
                break;

                default:
                C.GridPartType = GridPartType.METIS;
                break;
            }

            #endregion



            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.5, 0.5 }; //0.5,0.5
            double radius = 0.25;


            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("cd1a6c18-2659-4405-bf56-1e461441c0a0");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");


            C.AddBoundaryCondition("wall_lower", VariableNames.LevelSet, PhiFunc);

            #endregion

            // Level-Set
            // =================
            #region Fourier

            int numSp = 1024;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for(int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                samplP[sp] = radius;
            }

            //double circum = 2.0 * Math.PI * radius;
            //double filter = (circum * 20.0) / ((double)numSp / 2.0);
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem)
            //{
            //    center = center,
            //    FourierEvolve = Fourier_Evolution.FourierPoints,
            //    //Timestepper = FourierLevelSet_Timestepper.AdamsBashforth2,
            //    //UnderRelax = underrelax
            //    centerMove = CenterMovement.Reconstructed,
            //    //curvComp_extended = false
            //};


            #endregion

            C.Option_LevelSetEvolution = LevelSetEvolution.ScalarConvection;

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpoint+levelset" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            //C.AdvancedDiscretizationOptions.surfTensionMode = SurfaceTensionMode.Curvature_Fourier;
            //C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            //C.Timestepper_MassMatrix = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;

            C.CompMode = AppControl._CompMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            //double dt = 75e-4; // (1.0 / (double)kelem) / 16.0;
            //CFL condition > capillary time-step constraint
            double dt = 0.25 * Math.Sqrt((C.PhysicalParameters.rho_A + C.PhysicalParameters.rho_B) * Math.Pow((Math.Min(xSize, ySize) / kelem), 3) / (4 * Math.PI * C.PhysicalParameters.Sigma));

            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 3;
            C.saveperiod = 1;



            #endregion

            return C;

        }


    }
}
