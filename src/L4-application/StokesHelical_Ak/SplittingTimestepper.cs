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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Voronoi;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net.Appender;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Providers.LinearAlgebra;
using NSE_SIMPLE;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Data;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;


namespace StokesHelical_Ak {
    public class SplittingTimestepper {

        double[] OpAffine;

        CoordinateMapping Umap => Uvec.Mapping;

        readonly internal CoordinateVector Uvec;

        internal CoordinateVector U_convec;

        internal CoordinateVector U_0;
        double dt;

        Action<BlockMsrMatrix, double[], UnsetteledCoordinateMapping, DGField[], Dictionary<SpeciesId, MultidimensionalArray>, double, bool> Impl;

        readonly int m_BdfOrder;

        readonly MultigridOperator.ChangeOfBasisConfig[][] m_mgconfig;

        readonly R0fix m_R0fix;
        readonly double m_rMin;
        readonly bool m_PRP;

        public SplittingTimestepper(
            Action<BlockMsrMatrix, double[], UnsetteledCoordinateMapping, DGField[], Dictionary<SpeciesId, MultidimensionalArray>, double, bool> __Impl,
            DifferentialOperator Expl, ITemporalOperator MassMatrixOp, CoordinateMapping __U,
            double __dt, int BdfOrder,
            BoSSS.Solution.AdvancedSolvers.ISolverFactory lc,
            MultigridOperator.ChangeOfBasisConfig[][] mgconfig,
            AggregationGridData[] mgseq,
            R0fix r0fix, double _rMin, bool _PRP) {
            Uvec = new CoordinateVector(__U);
            dt = __dt;
            Impl = __Impl;
            m_BdfOrder = BdfOrder;
            m_mgconfig = mgconfig;
            m_R0fix = r0fix;
            m_rMin = _rMin;
            m_PRP = _PRP;

            // generate implicit solver
            // ========================
            {
                OpAffine = new double[Umap.LocalLength];
                OpMatrix = new BlockMsrMatrix(Umap, Umap);
                Impl(OpMatrix, OpAffine, Umap, __U.Fields.ToArray(), null, 0.0, false);

                // Check if PRP is required
                CheckNecassarityOfPRP(OpMatrix, Umap, _rMin < 10e-6);
            }

            // generate explicit evaluator
            // ===========================
            {
                ExplicitEval = Expl.GetEvaluatorEx(Umap.Fields, Umap.Fields.Take(3).ToArray(), Umap);
            }

            // generate mass matrix
            // ====================
            {
                var mtxBuilder = MassMatrixOp.GetMassMatrixBuilder(Umap, null, Umap);
                mtxBuilder.time = 0.0;
                var dummy = new double[Umap.LocalLength];
                MassMatrix = new BlockMsrMatrix(Umap, Umap);
                mtxBuilder.ComputeMatrix(MassMatrix, dummy);

                if (_PRP == true) { // if(BC requires PRP){ //What BC exactly need PRP ?!?
                                    // ++++++++++++++++++++++++++++++++++++
                                    // apply pressure reference to mass matrix
                                    // (should not change anything, since the conti/pressure block is zero in the mass matrix)
                                    // ++++++++++++++++++++++++++++++++++++

                    // We don't manipulate the matrix, we only check that it fullfills
                    // the requirements imposed by the pressure reference point.
                    this.Umap.GridDat.LocatePoint(new double[] { 0.5, 0.5 }, out _, out long GlobalIndex, out _, out bool onthisProc);
                    if (onthisProc) {
                        long iRowGl = Umap.GlobalUniqueCoordinateIndex_FromGlobal(3, GlobalIndex, 0);

                        if (OpMatrix.GetNoOfOffDiagonalNonZerosPerRow(iRowGl) != 0)
                            throw new ArithmeticException("pressure ref pt is to be used, but there are off-diagonal non-zeros in operator matrix");
                        if (OpMatrix[iRowGl, iRowGl] != 1.0)
                            throw new ArithmeticException("pressure ref pt is to be used, diagonal element is expected to be 1.0, but is " + OpMatrix[iRowGl, iRowGl]);

                        if (MassMatrix.GetNoOfNonZerosPerRow(iRowGl) != 0)
                            throw new ArithmeticException("pressure ref pt is to be used, but mass matrix is occupied in respective row.");
                    }
                    //}
                }



                Debug.Assert(dummy.L2NormPow2() == 0);

                double factor;
                if (m_BdfOrder == 1) {
                    factor = 1 / dt;
                } else if (m_BdfOrder == 3) {
                    factor = Globals.beta0 / dt;
                } else {
                    throw new NotSupportedException("only supporting BDF1 and BDF3, not BDF(" + m_BdfOrder + ")");
                }


                if (_rMin < 10e-6) {
                    // ++++++++++++++++++++++++++++++++++++
                    // apply R0 fix to mass matrix
                    // ++++++++++++++++++++++++++++++++++++

                    var MassMatrixR0 = r0fix.GetManipulatedMatrix(MassMatrix);
                    OpMatrix.Acc(factor, MassMatrixR0);

                } else {
                    OpMatrix.Acc(factor, MassMatrix);
                }
            }

            // setup lin. solver
            // ====================
            {
                var MgBasis = AggregationGridBasis.CreateSequence(mgseq, Umap.BasisS.ToArray());

                //CheckMinEigenVecAndValue(Uvec, gridData, "BeforeSetupMgOperator", r0fix);

                MgOperator = new MultigridOperator(MgBasis, Umap, OpMatrix, null, mgconfig, null);
                var slv = lc.CreateInstance(MgOperator);
                slv.Init(MgOperator);
                newSolver = slv;
                if (newSolver is DirectSolver directSolver) {
                    directSolver.ActivateCaching = (iter, lvl) => true;
                } else if (newSolver is ISolverWithCallback solverCallback) {
                    solverCallback.IterationCallback = delegate (int iIter, double[] X, double[] Res, MultigridOperator mgop) {
                        double ResNorm = Res.MPI_L2Norm();
                        Console.WriteLine($"{iIter} {ResNorm:0.##e-00}");
                    };
                };
            }
            // Set dummy values for BDF history
            // ================================
            // Ich setze hier meinn Zeug an, um zu checken, ob ich wirklich die Exakte Lösung berechnen kann!
            // Ersetze den Pfad mit dem tatsächlichen Pfad deiner Datei


            // Lese den Vektor aus der Datei und speichere ihn als double[]
            Un = Uvec.ToArray();
            Unn = Uvec.ToArray();
            Unnn = Uvec.ToArray();
            Cnnn = Compute_C(Unnn, 0);
            Cnn = Compute_C(Unn, 0);
            Cn = Compute_C(Un, 0);
        }

        /// <summary>
        /// Mass Matrix without R0 fix
        /// </summary>
        BlockMsrMatrix MassMatrix;

        BlockMsrMatrix OpMatrix;

        MultigridOperator MgOperator;

        IEvaluatorNonLin ExplicitEval;

        ISolverSmootherTemplate newSolver;

        //public double[] Un_von_Datei;

        public double[] Un;
        public double[] Unn;
        public double[] Unnn;

        public double[] Cn;
        public double[] Cnn;
        public double[] Cnnn;

        // 
        double[] Compute_C(double[] u, double time) {
            double[] C = new double[Umap.LocalLength];
            //var UvecBkup = Uvec.ToArray();
            Uvec.SetV(u);
            ExplicitEval.time = time;
            ExplicitEval.Evaluate(1.0, 0.0, C);
            //Uvec.SetV(UvecBkup);
            return C;
        }

        /// <summary>
        /// my method to compute the correct initial values for BDF3 using the exact solution at -dt and -2*dt 
        /// </summary>
        public void SetInitialExact(double dt, double phystime, CoordinateVector _Unnn, CoordinateVector _Unn) {
            Unnn.Clear();
            Unn.Clear();

            //Unnn = Un_von_Datei;
            //Unn = Un_von_Datei;

            Unnn = _Unnn.ToArray();
            Unn = _Unn.ToArray();
            Cnn.Clear();
            Cn.Clear();
            //====================================================================================================
            // Attention. I have to use Impl to update the Boundary Values with the respective phystime for the convective terms!
            // Important for Higer Order Time Stepping Scemes! Since Values for Earlier Time Steps are necessary
            Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, phystime - 2 * dt, true);

            Cnnn = Compute_C(Unnn, phystime - 2 * dt);

            Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, phystime - dt, true);

            Cnn = Compute_C(Unn, phystime - dt);

            Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, phystime, true);

            Cn = Compute_C(Un, phystime);
            //====================================================================================================
        }

        static public int BackupTimestep = -1;

        internal static double[] Backup_Cnnn;
        internal static double[] Backup_Cnn;
        internal static double[] Backup_Cn;
        internal static double[] Check_Cnnn;
        internal static double[] Check_Cnn;
        internal static double[] Check_Cn;
        internal static double[] Backup_Unnn;
        internal static double[] Backup_Unn;
        internal static double[] Backup_Un;
        internal static double BackupTime;

        public void Solve(double time, R0fix m_R0fix, int restartTimeStep, int currentTimestep, bool restart, double[] Unnn_, double[] Unn_, double[] Un_, int bdfOrder) {
            //Shift timesteps
            // ===============
            if (restart && restartTimeStep + 1 == currentTimestep && bdfOrder == 3) {
                Unnn = Unnn_;
                Unn = Unn_;
                ComputeOldValuesFromRestart(Unnn, Unn, time);
            } else {
                if (currentTimestep != 1) {
                    Unnn = Unn;
                    Unn = Un;
                    Cnnn = Cnn;
                    Cnn = Cn;
                }
                // ComputeDifferences_Conv_Output(time, Umap.GridDat, "_Conv_Check", m_R0fix);
                //Console.WriteLine("ZEILE WIEDER AUSKOMMENTIEREN:");
                //Console.WriteLine("SPLITTING TIMESTEPPER ZEILE 290");
                //Console.WriteLine("Convektive Terme werden geplottet");

                ComputeDifferences(time, currentTimestep);
            }

            // Build RHS
            // =========
            int L = Umap.LocalLength;
            var RHS = new double[L];
            OpAffine.Clear();
            Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, time, true);
            ApplyBDFOrder(RHS, currentTimestep, restart, time, restartTimeStep); // = (1/dt)*MM*u0, 
            RHS = m_R0fix.GetRestrictedRHS(RHS); // = R0Restriction * (1/dt)*MM*u0, 

            RHS.AccV(-1.0, OpAffine); // = R0Restriction * (1/dt)*MM*u0 + RHS_SS, note that `OpAffine` already "contains" the R0fix, so we must NOT apply it to `OpAffine`

            ExecuteSolver(RHS, currentTimestep);

            // prolongate
            var sol = m_R0fix.GetExtendedSolution(this.Uvec);
            this.Uvec.SetV(sol);
        }

        public void CheckMinEigenVecAndValue(CoordinateVector solution, IGridData gridData, string name, R0fix _R0fix) {
            (var lambdaMin, var EigSolMin) = this.OpMatrix.MinimalEigen();
            Console.WriteLine("Minimal Eigenvalue is: " + lambdaMin);
            Console.WriteLine("Eigenvector for resp. mininaml Eigenvalue is: " + EigSolMin);
            solution.SetV(EigSolMin);
            solution.SetV(_R0fix.GetExtendedSolution(solution));
            Tecplot c = new Tecplot(gridData, 1);
            c.PlotFields("ur_" + name, 0.0, solution.Fields[0]);
            c.PlotFields("uxi_" + name, 0.0, solution.Fields[1]);
            c.PlotFields("ueta_" + name, 0.0, solution.Fields[2]);
            c.PlotFields("psi_" + name, 0.0, solution.Fields[3]);
        }

        public void Solve(double time, int restartTimeStep, int currentTimestep, bool restart, double[] Unnn_, double[] Unn_, double[] Un_, int bdfOrder) {
            // Shift timesteps
            // ===============
            if (restart && restartTimeStep + 1 == currentTimestep && bdfOrder == 3) { // For restart
                Unnn = Unnn_;
                Unn = Unn_;
                ComputeOldValuesFromRestart(Unnn, Unn, time);
            } else {  // Normal
                if (currentTimestep != 1) {
                    Unnn = Unn;
                    Unn = Un;
                    Cnnn = Cnn;
                    Cnn = Cn;
                }
                //Console.WriteLine("ZEILE WIEDER AUSKOMMENTIEREN:");
                //Console.WriteLine("SPLITTING TIMESTEPPER ZEILE 350");
                //Console.WriteLine("Convektive Terme werden geplottet");
                //ComputeDifferences_Conv_Output(time, Umap.GridDat, "_Conv_Check");
                ComputeDifferences(time, currentTimestep);
            }
            // Build RHS
            // =========
            int L = Umap.LocalLength;
            var RHS = new double[L];
            OpAffine.Clear();
            Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, time, true); // Update der RB?!?! 
            RHS.AccV(-1.0, OpAffine);
            ApplyBDFOrder(RHS, currentTimestep, restart, time, restartTimeStep);
            VerifyRHS(RHS);
            ExecuteSolver(RHS, currentTimestep);
        }


        void VerifyRHS(double[] RHS) {


            if (m_PRP) {
                // +++++++++++++++++++++++++++++++++++++++++++++
                // does it fulfill the pressure reference point?
                // +++++++++++++++++++++++++++++++++++++++++++++

                this.Umap.GridDat.LocatePoint(new double[] { 0.5, 0.5 }, out _, out long GlobalIndex, out _, out bool onthisProc);
                if (onthisProc) {
                    int jLoc = this.Umap.GridDat.CellPartitioning.Global2Local(GlobalIndex);
                    int iRowLoc = this.Umap.LocalUniqueCoordinateIndex(3, jLoc, 0);

                    if (RHS[iRowLoc] != 0)
                        throw new ArithmeticException("RHS is nonzero for the pressure reference point.");
                }
            }
        }

        private void ComputeOldValuesFromRestart(double[] Unnn_, double[] Unn_, double time) {
            Un = Uvec.ToArray(); // creates a copy !
                                 // explicit evaluation 
                                 // ===================
                                 // Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, time - dt * 2, true);
            Cnnn = Compute_C(Unnn, time - dt * 2);
            var Check_Cnnn = Cnnn.CloneAs();

            Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, time - dt, true);
            Cnn = Compute_C(Unn, time - dt);
            var Check_Cnn = Cnn.CloneAs();
            //Impl(null, OpAffine, Umap, Uvec.Mapping.Fields.ToArray(), null, time, true);
            Cn = Compute_C(Un, time);
            var Check_Cn = Cn.CloneAs();

        }

        private void ComputeDifferences(double time, int currentTimestep) {

            if (currentTimestep > 2) {
                //var UN_zwischen = Uvec.ToArray();
                var diff = Un.L2Distance(Unnn);
                var denNom = Un.L2Norm();
                Console.WriteLine($"{diff / denNom}");

                double typee = currentTimestep / 20.0;

                // Überprüfen, ob 'typee' eine ganze Zahl ist
                if (typee % 1 == 0) {
                    Console.WriteLine("");
                }

                if (diff / denNom <= 10E-6) {
                    Un.SaveToTextFile($"Solution_StokesSystem_with_dt={this.dt}");
                    //RHS.SaveToTextFile($"RHS_StokesSystem_with_dt={this.dt}");
                    //Assert.That(diff / denNom >= 10E-3, "SteadyStateReached");
                }
            }

            Un = Uvec.ToArray(); // creates a copy !


            // explicit evaluation
            // ===================
            Cn = Compute_C(Un, time);
        }

        private void ComputeDifferences_Conv_Output(double time, IGridData gridData, string name, R0fix _R0fix) {

            Un = Uvec.ToArray(); // creates a copy !
                                 // explicit evaluation
                                 // ===================
            U_0 = new CoordinateVector(Umap).CloneAs();
            U_0.SetV(Un);


            //Tecplot c_Initial = new Tecplot(gridData, 1);
            //c_Initial.PlotFields("u0_r", 0.0, U_0.Fields[0]);
            //c_Initial.PlotFields("u0_xi", 0.0, U_0.Fields[1]);
            //c_Initial.PlotFields("u0_eta", 0.0, U_0.Fields[2]);
            //c_Initial.PlotFields("p0", 0.0, U_0.Fields[3]);

            // calculate Convective Terms before R0 fix
            Cn = Compute_C(Un, time);
            U_convec = U_0.CloneAs(); U_convec.RenameFields(idx => "Mom" + (idx + 1));
            U_convec.SetV(Cn);
            //MassMatrix.Solve_Direct(U_convec, Cn);


            var U_0P0 = U_0.CloneAs(); U_0P0.RenameFields("urP0", "uxiP0", "uEtaP0", "pP0");
            foreach (var u in U_0P0.Fields) {
                int Np = u.Coordinates.NoOfCols;
                for (int n = 1; n < Np; n++)
                    u.Coordinates.ClearCol(n);
            }
            var U_convecP0 = U_0.CloneAs(); U_convecP0.RenameFields(idx => "Mom" + (idx + 1) + "P0");

            var Cn0 = Compute_C(U_0P0.ToArray(), 0.0);
            U_convecP0.SetV(Cn0);



            // calculate Convective Terms after R0 fix
            var U_convecExt = U_0.CloneAs(); U_convecExt.RenameFields(idx => "Mom" + (idx + 1) + "Ext");
            U_convecExt.SetV(_R0fix.GetExtendedSolution(U_convec));


            Tecplot.PlotFields(U_0.Fields.Cat(U_convec.Fields, U_0P0.Fields, U_convecP0.Fields, U_convecExt.Fields), "Convective", 0.0, 2);
            Console.WriteLine("Remove me");
        }


        private void ComputeDifferences_Conv_Output(double time, IGridData gridData, string name) {

            Un = Uvec.ToArray(); // creates a copy !
                                 // explicit evaluation
                                 // ===================
            U_0 = new CoordinateVector(Umap).CloneAs();


            Un.FillRandom(0); // FIll in Random
            // random_With_Pres_Offset.SetV(random); // random_With_Pres_Offset <- random // Creat a real Copy not Reference type!


            U_0.SetV(Un);


            //Tecplot c_Initial = new Tecplot(gridData, 1);
            //c_Initial.PlotFields("u0_r", 0.0, U_0.Fields[0]);
            //c_Initial.PlotFields("u0_xi", 0.0, U_0.Fields[1]);
            //c_Initial.PlotFields("u0_eta", 0.0, U_0.Fields[2]);
            //c_Initial.PlotFields("p0", 0.0, U_0.Fields[3]);

            // calculate Convective Terms before R0 fix
            Cn = Compute_C(Un, time);
            U_convec = U_0.CloneAs();
            U_convec.SetV(Cn);
            U_convec.RenameFields(idx => "Mom" + (idx + 1));
            //MassMatrix.Solve_Direct(U_convec, Cn);


            var U_0P0 = U_0.CloneAs();
            U_0P0.RenameFields("urP0", "uxiP0", "uEtaP0", "pP0");
            foreach (var u in U_0P0.Fields) {
                int Np = u.Coordinates.NoOfCols;
                for (int n = 1; n < Np; n++) {
                    u.Coordinates.ClearCol(n);
                }
                CellMask Domain = CellMask.GetFullMask(gridData);
                Console.WriteLine("Mean_of_" + u.ToString() + "_=_" + u.GetMeanValueTotal(Domain));
            }
            var U_convecP0 = U_0.CloneAs();
            U_convecP0.RenameFields(idx => "Mom" + (idx + 1) + "P0");

            var Cn0 = Compute_C(U_0P0.ToArray(), 0.0);
            U_convecP0.SetV(Cn0);

            foreach (var u in U_convecP0.Fields) {
                CellMask Domain = CellMask.GetFullMask(gridData);
                Console.WriteLine("Mean_of_" + u.ToString() + "_=_" + u.GetMeanValueTotal(Domain));
            }


            Tecplot.PlotFields(U_0.Fields.Cat(U_convec.Fields, U_0P0.Fields, U_convecP0.Fields), "Convective", 0.0, 2);
            Console.WriteLine("Remove me");
        }


        private void ApplyBDFOrder(double[] RHS, int currentTimestep, bool restart, double time, int restartTimeStep) {
            // Implementation for applying BDF order (either 1 or 3)
            if (m_BdfOrder == 3) {
                ApplyBDF3(RHS, currentTimestep, restart, time, restartTimeStep);
            } else if (m_BdfOrder == 1) {
                ApplyBDF1(RHS);
            } else {
                throw new NotSupportedException($"BDF order {m_BdfOrder} not supported");
            }
        }

        private void ApplyBDF3(double[] RHS, int currentTimestep, bool restart, double time, int restartTimeStep) {


            if (BackupTimestep + 1 == currentTimestep && restart == false) {

                Backup_Cn = Cn.CloneAs();
                Backup_Cnn = Cnn.CloneAs();
                Backup_Cnnn = Cnnn.CloneAs();

                Backup_Un = Un.ToArray();
                Backup_Unn = Unn.CloneAs();
                Backup_Unnn = Unnn.CloneAs();


                BackupTime = time;
            } else if (restart && restartTimeStep + 1 == currentTimestep) {

                var dist_nnn = Backup_Unnn.L2Distance(Unnn);
                var dist_nn = Backup_Unn.L2Distance(Unn);
                var dist_n = Backup_Un.L2Distance(Un);

                var cnnn_dist = Cnnn.L2Distance(Backup_Cnnn);
                var cnn_dist = Cnn.L2Distance(Backup_Cnn);
                var cn_dist = Cn.L2Distance(Backup_Cn);


                var time_dist = time - BackupTime;


                Console.WriteLine("");

                Console.Write("");
            }
            // MassMatrix
            MassMatrix.SpMV(Globals.beta1 / dt, Un, 1.0, RHS);
            MassMatrix.SpMV(Globals.beta2 / dt, Unn, 1.0, RHS);
            MassMatrix.SpMV(Globals.beta3 / dt, Unnn, 1.0, RHS);
            // RHS
            RHS.AccV(-1.0 * Globals.g1, Cn);
            RHS.AccV(-1.0 * Globals.g2, Cnn);
            RHS.AccV(-1.0 * Globals.g3, Cnnn);
        }

        private void ApplyBDF1(double[] RHS) {

            //double errUn_l2 = Un.MPI_L2Dist(Uinfty);
            //Console.WriteLine("   errUn_l2 = " + errUn_l2);
            //Un.SetV(Uinfty);

            // MassMatrix
            MassMatrix.SpMV(1.0 / dt, Un, 1.0, RHS);
            // RHS
            RHS.AccV(-1.0, Cn);

            //var Conv = this.Uvec.CloneAs();
            //Conv.RenameFields(i => "Conv" + i);
            //Conv.SetV(Cn);
            //Tecplot.PlotFields(Conv.Fields, "convective", 0.0, 3);
        }
        private void CheckNecassarityOfPRP(BlockMsrMatrix oPmatrix, UnsetteledCoordinateMapping map, bool containsR0fix) {



            // Definition
            double[] result1 = new double[map.LocalLength];
            double[] result2 = new double[map.LocalLength];
            //double[] random = new double[map.GlobalCount];
            var random = new CoordinateVector(map.BasisS.Select(basis => new SinglePhaseField(basis)));
            var random_With_Pres_Offset = new CoordinateVector(map.BasisS.Select(basis => new SinglePhaseField(basis)));

            double[] diff = new double[map.LocalLength];
            // Fill in Randoms
            random.FillRandom(0);
            random_With_Pres_Offset.SetV(random); // random_With_Pres_Offset <- random // Creat a real Copy not Reference type!

            var pres = random_With_Pres_Offset.Fields[3]; // Fields[3] for pressure of course
            if (containsR0fix) {
                var R0mask = R0fix.GetR0BndyCells(map.GridDat);
                pres.AccConstant(1.0, R0mask.Complement());

                if (map.MpiRank == 0) {
                    pres.SetMeanValue(0, pres.GetMeanValue(0) + 1.0);
                }
            } else {
                // ++++++++
                // r0 >> 0
                // ++++++++
                pres.AccConstant(1.0);
            }
            // Matrix Vector Multiplication
            oPmatrix.SpMV(1.0, random, 0.0, result1);
            oPmatrix.SpMV(1.0, random_With_Pres_Offset, 0.0, result2);

            diff.SetV(result2);
            diff.AccV(-1.0, result1);

            Assert.That(diff.MPI_L2Norm() > 1E-5, "Pressure is not FIXED! Residual should change at least by {0} but only changes {1} when constant pressure is added", 1E-5, diff.Max() - diff.Min());
        }

        private void ExecuteSolver(double[] RHS, int currentTimestep) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;

                //RHS.Clear();
                //RHS.AccV(1.0, RHS_SS);
                //RHS.AccV(1.0, RHS_instat);

                MgOperator.UseSolver(newSolver, Uvec, RHS);
                tr.Info($"Solver {newSolver} converged? {newSolver.Converged}, no of iterations: {newSolver.ThisLevelIterations}");



                //  RHS.AccV(1.0, RHS_instat);

                //double errSol_l2 = Uvec.MPI_L2Dist(Uinfty);
                //Console.WriteLine("Solution distance: " + errSol_l2);
                // Console.WriteLine("Solution distance: " + errSol_l2);

                /*

                newSolver.ResetStat();
                var Resi = RHS.CloneAs();
                double[] RHSpc = new double[RHS.Length];
                MgOperator.TransformRhsInto(RHS, RHSpc, false);
                double[] Upc = new double[Uvec.Length];
                MgOperator.TransformSolInto(Uinfty, Upc);


                Resi.SetV(RHSpc);
                MgOperator.OperatorMatrix.SpMV(-1.0, Upc, 1.0, Resi);
                double L2Resi = Resi.L2Norm();
                Console.WriteLine(" L2Resi = " + L2Resi);


                var ResiVec = new CoordinateVector(Uvec.Fields.Select(f => f.CloneAs()));
                MgOperator.TransformRhsFrom(ResiVec, Resi);
                Tecplot.PlotFields(ResiVec.Fields, "Resi", 0.0, 2);



                MgOperator.OperatorMatrix.SaveToTextFileSparse($"OpMatrix_StokesSystem_with_dt={dt}.txt");
                RHSpc.SaveToTextFile($"RHS_StokesSystem_with_dt={this.dt}.txt");
                Upc.SaveToTextFile($"Solution_StokesSystem_with_dt={this.dt}.txt");
                */
            }
        }




        /// <summary>
        /// Returns a collection of local and global condition numbers in order to assess the operators stability,
        /// <see cref="IApplication.OperatorAnalysis"/>.
        /// </summary>
        internal IDictionary<string, double> OperatorAnalysisAk(IEnumerable<int[]> VarGroups = null,
            bool calculateGlobals = true,
            bool calclulateStencils = true,
            bool plotStencilCondNumViz = false,
            string nameOfStencil = null) {

            long J = this.Umap.GridDat.CellPartitioning.TotalLength;
            var StencilCondNoVizS = new List<DGField>();
            var Ret = new Dictionary<string, double>();
            foreach (int[] varGroup in VarGroups) {
                var ana = new BoSSS.Solution.AdvancedSolvers.Testing.OpAnalysisBase(this.OpMatrix, this.OpAffine, this.Umap, m_mgconfig, null);

                ana.CalculateGlobals = calculateGlobals;
                ana.CalculateStencils = calclulateStencils;

                ana.VarGroup = varGroup;
                var Table = ana.GetNamedProperties();

                foreach (var kv in Table) {
                    if (!Ret.ContainsKey(kv.Key)) {
                        Ret.Add(kv.Key, kv.Value);
                    }
                }

                if (plotStencilCondNumViz) {
                    StencilCondNoVizS.Add(ana.StencilCondNumbersV());
                }


            }

            if (StencilCondNoVizS.Count > 0) {
                Tecplot.PlotFields(StencilCondNoVizS, nameOfStencil, 0.0, 1);
            }

            return Ret;
        }

    }
}
