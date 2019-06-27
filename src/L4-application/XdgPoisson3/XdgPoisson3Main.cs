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
using System.Linq;
using System.Collections.Generic;
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using MPI.Wrappers;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution.Tecplot;
using System.Globalization;
using System.IO;
using ilPSP.Utils;
using ilPSP.Connectors.Matlab;
using System.Text;
using ilPSP;
using BoSSS.Solution.AdvancedSolvers;
using NUnit.Framework;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;

namespace BoSSS.Application.XdgPoisson3 {

    public class XdgPoisson3Main : BoSSS.Solution.Application<XdgPoisson3Control> {



        static void Main(string[] args) {
            BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;
            //BatchmodeConnector.MatlabExecuteable = "D:\\cygwin\\bin\\bash.exe";

            BoSSS.Solution.Application<XdgPoisson3Control>._Main(args, false, delegate () {
                return new XdgPoisson3Main();
            });
        }

#pragma warning disable 649

        [LevelSetTracker("-:A +:B", 1)]
        LevelSetTracker LevelSetTracker;

        [InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        LevelSet Phi;

        [InstantiateFromControlFile("u", "u", IOListOption.ControlFileDetermined)]
        XDGField u;

        /// <summary>
        /// exact solution (if known)
        /// </summary>
        [InstantiateFromControlFile("uEx", "u", IOListOption.ControlFileDetermined)]
        XDGField uEx;

        /// <summary>
        /// error of numerical solution (if exact solution known)
        /// </summary>
        [InstantiateFromControlFile("uErr", "u", IOListOption.ControlFileDetermined)]
        XDGField uErr;

        [InstantiateFromControlFile("rhs", "u", IOListOption.ControlFileDetermined)]
        XDGField rhs;

        [InstantiateFromControlFile("residual", "u", IOListOption.ControlFileDetermined)]
        XDGField residual;

        [InstantiateFromControlFile(
            new string[] { "du_dx", "du_dy", "du_dz" },
            new string[] { "u", "u", "u" },
            true, true,
            IOListOption.ControlFileDetermined)]
        VectorField<XDGField> GradientU;
#pragma warning restore 649

        BlockMsrMatrix Op_Matrix;
        double[] Op_Affine;
        MultiphaseCellAgglomerator Op_Agglomeration;
        MassMatrixFactory Op_mass;

        static void MyHandler(object sender, UnhandledExceptionEventArgs args) {
            Exception e = (Exception)args.ExceptionObject;
            Console.WriteLine("MyHandler caught : " + e.Message);
            Console.WriteLine("Runtime terminating: {0}", args.IsTerminating);
            System.Environment.Exit(-1234);
        }

        protected override void SetInitial() {
            //this will suppress exception prompts
            //Workaround to prevent distrubance while executing batchclient
            if (this.Control.SuppressExceptionPrompt) {
                AppDomain currentDomain = AppDomain.CurrentDomain;
                currentDomain.UnhandledException += new UnhandledExceptionEventHandler(MyHandler);
            }

            base.SetInitial();
            this.LsTrk.UpdateTracker();
            base.SetInitial();
            this.LsTrk.UpdateTracker();

            this.MGColoring = new SinglePhaseField[base.MultigridSequence.Length];
            for (int iLevel = 0; iLevel < base.MultigridSequence.Length; iLevel++) {
                this.MGColoring[iLevel] = new SinglePhaseField(new Basis(this.GridData, 0), "MGColoring_level_" + iLevel);
                base.MultigridSequence[iLevel].ColorDGField(this.MGColoring[iLevel]);
            }


        }

        protected override void CreateFields() {
            base.CreateFields();
            base.LsTrk = LevelSetTracker;
            if (Control.CutCellQuadratureType != base.LsTrk.CutCellQuadratureType)
                throw new ApplicationException();
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            AssembleMatrix(this.Control.MU_A, this.Control.MU_B, out Op_Matrix, out Op_Affine, out Op_Agglomeration, out Op_mass);
            Console.WriteLine("Matrix norm: {0}", Op_Matrix.InfNorm());
            Console.WriteLine("Symm. diff: {0}", Op_Matrix.SymmetryDeviation());
        }

        protected void BlockTest() {
            double MU_A = 0.1;

            BlockMsrMatrix Msc;
            double[] bsc;
            MultiphaseCellAgglomerator aggsc;
            MassMatrixFactory mass_sc;
            AssembleMatrix(MU_A, 1.0, out Msc, out bsc, out aggsc, out mass_sc);

            BlockMsrMatrix Mref;
            double[] bRef;
            MultiphaseCellAgglomerator aggRef;
            MassMatrixFactory massRef;
            AssembleMatrix(1.0, 1.0, out Mref, out bRef, out aggRef, out massRef);

            SpeciesId IdA = this.LsTrk.GetSpeciesId("A");
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int N = this.u.Basis.NonX_Basis.Length;
            int[] ColIdx = null;
            for (int j = 0; j < J; j++) {

                ReducedRegionCode rrc;
                this.LsTrk.Regions.GetNoOfSpecies(j, out rrc);
                int SpAIdx = this.LsTrk.GetSpeciesIndex(rrc, IdA);
                if (SpAIdx < 0)
                    continue;

                int n0 = N * SpAIdx;
                int nE = n0 + N;
                for (int n = n0; n < nE; n++) {
                    int i = (int)u.Mapping.GlobalUniqueCoordinateIndex(0, j, n);
                    //var row = Mref.GetRowShallow(i);
                    //for (int k = 0; k < row.Length; k++)
                    //    row[k].Value *= MU_A;
                    int LL = Mref.GetOccupiedColumnIndices(i, ref ColIdx);
                    for (int ll = 0; ll < LL; ll++) {
                        Mref[i, ColIdx[ll]] *= MU_A;
                    }
                }
            }

            Mref.Acc(-1.0, Msc);
            double MdiffNorm = Mref.InfNorm();
            Console.WriteLine("matrix difference norm {0}", MdiffNorm);

            Console.WriteLine("Symm. diff: {0}", Msc.SymmetryDeviation());
        }

        private void AssembleMatrix(double MU_A, double MU_B, out BlockMsrMatrix M, out double[] b, out MultiphaseCellAgglomerator agg, out MassMatrixFactory massFact) {
            using (var tr = new FuncTrace()) {
                // create operator
                // ===============

                if (this.Control.SetDefaultDiriBndCnd) {
                    this.Control.xLaplaceBCs.g_Diri = ((CommonParamsBnd inp) => 0.0);
                    this.Control.xLaplaceBCs.IsDirichlet = (inp => true);
                }

                double D = this.GridData.SpatialDimension;
                int p = u.Basis.Degree;
                double penalty_base = (p + 1) * (p + D) / D;
                double penalty_multiplyer = base.Control.penalty_multiplyer;

                XQuadFactoryHelper.MomentFittingVariants momentFittingVariant;
                if (D == 3)
                    momentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;

                momentFittingVariant = this.Control.CutCellQuadratureType;

                int order = this.u.Basis.Degree * 2;

                XSpatialOperatorMk2 Op = new XSpatialOperatorMk2(1, 1, (A, B, C) => order, "u", "c1");
                var lengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)GridData).Cells.PenaltyLengthScales;
                var lap = new XLaplace_Bulk(this.LsTrk, penalty_multiplyer * penalty_base, "u", this.Control.xLaplaceBCs, 1.0, MU_A, MU_B, lengthScales, this.Control.ViscosityMode);
                Op.EquationComponents["c1"].Add(lap);      // Bulk form
                Op.EquationComponents["c1"].Add(new XLaplace_Interface(this.LsTrk, MU_A, MU_B, penalty_base * 2, this.Control.ViscosityMode));   // coupling form

                Op.Commit();

                // create agglomeration
                // ====================

                var map = new UnsetteledCoordinateMapping(u.Basis);

                //agg = new MultiphaseCellAgglomerator(
                //    new CutCellMetrics(momentFittingVariant,
                //        QuadOrderFunc.SumOfMaxDegrees(RoundUp: true)(map.BasisS.Select(bs => bs.Degree).ToArray(), new int[0], map.BasisS.Select(bs => bs.Degree).ToArray()),
                //        //this.HMFDegree,
                //        LsTrk, this.LsTrk.SpeciesIdS.ToArray()),
                //    this.Control.AgglomerationThreshold, false);
                agg = LsTrk.GetAgglomerator(this.LsTrk.SpeciesIdS.ToArray(), order, this.Control.AgglomerationThreshold);

                // compute matrix
                // =============
                using (new BlockTrace("XdgMatrixAssembly", tr)) {
                    M = new BlockMsrMatrix(map, map);
                    b = new double[M.RowPartitioning.LocalLength];

                    //Op.ComputeMatrixEx(LsTrk,
                    //    map, null, map,
                    //    M, b, false, 0.0, true,
                    //    agg.CellLengthScales, null, null, //out massFact,
                    //    this.LsTrk.SpeciesIdS.ToArray());
                    XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(this.LsTrk, map, null, map, this.LsTrk.SpeciesIdS.ToArray());
                    mtxBuilder.time = 0.0;
                    mtxBuilder.MPITtransceive = true;
                    mtxBuilder.ComputeMatrix(M, b);
                }
                // compare with linear evaluation
                // ==============================
                DGField[] testDomainFieldS = map.BasisS.Select(bb => new XDGField(bb as XDGBasis)).ToArray();
                CoordinateVector test = new CoordinateVector(testDomainFieldS);

                DGField[] errFieldS = map.BasisS.Select(bb => new XDGField(bb as XDGBasis)).ToArray();
                CoordinateVector Err = new CoordinateVector(errFieldS);

                var eval = Op.GetEvaluatorEx(LsTrk,
                    testDomainFieldS, null, map);

                foreach (var s in this.LsTrk.SpeciesIdS)
                    eval.SpeciesOperatorCoefficients[s].CellLengthScales = agg.CellLengthScales[s];

                eval.time = 0.0;
                int L = test.Count;
                Random r = new Random();
                for (int i = 0; i < L; i++) {
                    test[i] = r.NextDouble();
                }



                double[] R1 = new double[L];
                double[] R2 = new double[L];
                eval.Evaluate(1.0, 1.0, R1);

                R2.AccV(1.0, b);
                M.SpMV(1.0, test, 1.0, R2);

                Err.AccV(+1.0, R1);
                Err.AccV(-1.0, R2);

                double ErrDist = GenericBlas.L2DistPow2(R1, R2).MPISum().Sqrt();

                double Ref = test.L2NormPow2().MPISum().Sqrt();

                Assert.LessOrEqual(ErrDist, Ref * 1.0e-5, "Mismatch between explicit evaluation of XDG operator and matrix.");


                // agglomeration wahnsinn
                // ======================

                agg.ManipulateMatrixAndRHS(M, b, map, map);

                foreach (var S in this.LsTrk.SpeciesNames) {
                    Console.WriteLine("  Species {0}: no of agglomerated cells: {1}",
                        S, agg.GetAgglomerator(this.LsTrk.GetSpeciesId(S)).AggInfo.SourceCells.NoOfItemsLocally);
                }


                // mass matrix factory
                // ===================

                Basis maxB = map.BasisS.ElementAtMax(bss => bss.Degree);
                //massFact = new MassMatrixFactory(maxB, agg);
                massFact = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order).MassMatrixFactory;
            }
        }


        SinglePhaseField[] MGColoring;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //if (base.Control.timeDependent) {
            //    dt = base.GetFixedTimestep();
            //    Console.WriteLine("Timestep {0}, dt = {1} ...", TimestepNo, dt);
            //} else {
            base.TerminationKey = true;
            dt = 1.0;
            Console.WriteLine("Steady solve ...");
            //}

            double mintime, maxtime;
            bool converged;
            int NoOfIterations, DOFs;


            // direct solver 
            this.ReferenceSolve();
            //this.ConsistencyTest();


            // new solver framework: multigrid, blablablah ...

            ExperimentalSolver(out mintime, out maxtime, out converged, out NoOfIterations, out DOFs);
            this.Op_Agglomeration.Extrapolate(this.u.Mapping);

            //Stats:
            {
                int BlkSize_min = u.Mapping.MinTotalNoOfCoordinatesPerCell;
                int BlkSize_max = u.Mapping.MaxTotalNoOfCoordinatesPerCell;
                int NoOfMtxBlocks = 0;
                foreach (int[] Neigs in this.GridData.iLogicalCells.CellNeighbours) {
                    NoOfMtxBlocks++; //               diagonal block
                    NoOfMtxBlocks += Neigs.Length; // off-diagonal block
                }
                NoOfMtxBlocks = NoOfMtxBlocks.MPISum();

                int MtxBlockSize = BlkSize_max * BlkSize_max;
                int MtxSize = MtxBlockSize * NoOfMtxBlocks;

                double MtxStorage = MtxSize * (8.0 + 4.0) / (1024 * 1024); // 12 bytes (double+int) per entry

                Console.WriteLine("   System size:                 {0}", u.Mapping.TotalLength);
                Console.WriteLine("   No of blocks:                {0}", u.Mapping.TotalNoOfBlocks);
                Console.WriteLine("   No of blocks in matrix:      {0}", NoOfMtxBlocks);
                Console.WriteLine("   DG coordinates per cell:     {0}", BlkSize_max);
                Console.WriteLine("   Non-zeros per matrix block:  {0}", MtxBlockSize);
                Console.WriteLine("   Total non-zeros in matrix:   {0}", MtxSize);
                Console.WriteLine("   maximal matrix storage (MB): {0}", MtxStorage);

                Console.WriteLine("DOF: {0}", DOFs);
                Console.WriteLine("min Blocksize: {0}", u.Mapping.MinTotalNoOfCoordinatesPerCell);
                Console.WriteLine("max Blocksize: {0}", u.Mapping.MaxTotalNoOfCoordinatesPerCell);

                base.QueryHandler.ValueQuery("maxBlkSize", u.Mapping.MaxTotalNoOfCoordinatesPerCell, true);
                base.QueryHandler.ValueQuery("minBlkSize", u.Mapping.MinTotalNoOfCoordinatesPerCell, true);
                base.QueryHandler.ValueQuery("NumberOfMatrixBlox", NoOfMtxBlocks, true);
                base.QueryHandler.ValueQuery("DOFs", DOFs, true);
            }

            base.QueryHandler.ValueQuery("minSolRunT", mintime, true);
            base.QueryHandler.ValueQuery("maxSolRunT", maxtime, true);
            base.QueryHandler.ValueQuery("Conv", converged ? 1.0 : 0.0, true);
            base.QueryHandler.ValueQuery("NoIter", NoOfIterations, true);

            Console.WriteLine("maximal Multigridlevel: {0}", MaxMlevel);
            base.QueryHandler.ValueQuery("maxMultigridlvl", MaxMlevel, true);

            Console.WriteLine("done.");

            if (this.Control.ExcactSolSupported) {
                this.uErr.Clear();
                this.uErr.Acc(+1.0, u);
                this.uErr.Acc(-1.0, uEx);

                double L2_ERR = this.uErr.L2Norm();

                base.QueryHandler.ValueQuery("L2_ERR", L2_ERR, true);

                this.LsTrk.UpdateTracker();
                int order = Math.Max(u.Basis.Degree, Math.Max(uErr.Basis.Degree, uEx.Basis.Degree)) * 2 + 1;
                //var scheme = new XQuadSchemeHelper(this.LsTrk, this.Control.HMFversion, this.LsTrk.SpeciesIdS.ToArray());
                var scheme = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


                var A_scheme = scheme.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("A"));
                var B_scheme = scheme.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("B"));

                ICompositeQuadRule<QuadRule> A_rule = A_scheme.Compile(this.GridData, order);
                ICompositeQuadRule<QuadRule> B_rule = B_scheme.Compile(this.GridData, order);

                foreach (var _rule in new ICompositeQuadRule<QuadRule>[] { A_rule, B_rule }) {
                    foreach (var ckR in _rule) {
                        ckR.Rule.Weights.CheckForNanOrInf(true, true, true);
                        ckR.Rule.Nodes.CheckForNanOrInf(true, true, true);
                    }
                }

                double L2_ERR_HMF_A = this.u.GetSpeciesShadowField("A").LxError(this.Control.InitialValues_Evaluators["uEx#A"].Vectorize(), null, A_rule).Sqrt();
                double L2_ERR_HMF_B = this.u.GetSpeciesShadowField("B").LxError(this.Control.InitialValues_Evaluators["uEx#B"].Vectorize(), null, B_rule).Sqrt();
                double L2_ERR_HMF = (L2_ERR_HMF_A.Pow2() + L2_ERR_HMF_B.Pow2()).Sqrt();
                base.QueryHandler.ValueQuery("L2_ERR_HMF_A", L2_ERR_HMF_A, true);
                base.QueryHandler.ValueQuery("L2_ERR_HMF_B", L2_ERR_HMF_B, true);
                base.QueryHandler.ValueQuery("L2_ERR_HMF", L2_ERR_HMF, true);


                Console.WriteLine("Error norm (standard):       " + L2_ERR);
                Console.WriteLine("Error norm (HMF, Species A): " + L2_ERR_HMF_A);
                Console.WriteLine("Error norm (HMF, Species B): " + L2_ERR_HMF_B);
                Console.WriteLine("Error norm (HMF):            " + L2_ERR_HMF);
            }

            return dt;
        }


        private void OperatorAnalysis() {
            AggregationGridBasis[][] XAggB = AggregationGridBasis.CreateSequence(base.MultigridSequence, u.Mapping.BasisS);

            XAggB.UpdateXdgAggregationBasis(this.Op_Agglomeration);

            var MultigridOp = new MultigridOperator(XAggB, this.u.Mapping,
                this.Op_Matrix,
                this.Op_mass.GetMassMatrix(new UnsetteledCoordinateMapping(this.u.Basis), false),
                OpConfig);

            var PcOpMatrix = MultigridOp.OperatorMatrix;
            //PcOpMatrix.SaveToTextFileSparse("C:\\temp\\Xpoisson.txt");

            MultidimensionalArray ret = MultidimensionalArray.Create(1, 4);

            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(PcOpMatrix, "OpMtx");
                bmc.Cmd("OpMtxSym = 0.5*(OpMtx + OpMtx');");
                bmc.Cmd("condNo = condest(OpMtxSym);");
                bmc.Cmd("eigMaxi = eigs(OpMtxSym,1,'lm')");
                bmc.Cmd("eigMini = eigs(OpMtxSym,1,'sm')");
                bmc.Cmd("lasterr");
                bmc.Cmd("[V,r]=chol(OpMtxSym);");
                bmc.Cmd("ret = [condNo, eigMaxi, eigMini, r]");
                bmc.GetMatrix(ret, "ret");

                bmc.Execute(false);
            }

            double condNo = ret[0, 0];
            double eigMaxi = ret[0, 1];
            double eigMini = ret[0, 2];
            double posDef = ret[0, 3] == 0 ? 1 : 0;

            Console.WriteLine("Condition number: {0:0.####E-00}", condNo);

            if (posDef == 0.0)
                Console.WriteLine("WARNING: Operator matrix is not positive definite.");

            base.QueryHandler.ValueQuery("condNo", condNo, false);
            base.QueryHandler.ValueQuery("eigMaxi", eigMaxi, false);
            base.QueryHandler.ValueQuery("eigMini", eigMini, false);
            base.QueryHandler.ValueQuery("posDef", posDef, false);
        }

        MultigridOperator.ChangeOfBasisConfig[][] OpConfig {
            get {
                int p = this.u.Basis.Degree;
                return new MultigridOperator.ChangeOfBasisConfig[][] {
                    new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] { 0 }, mode = this.Control.PrePreCond, Degree = p }
                    }
                };
            }
        }

        protected void CustomItCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            //noch nix ...
            MaxMlevel=Mgop.LevelIndex;
        }

        private int m_maxMlevel;

        public int MaxMlevel {
            get{
                return m_maxMlevel;
            }
            set{
                if (value > m_maxMlevel)
                    m_maxMlevel = value;
            }
        }

        private void ExperimentalSolver(out double mintime, out double maxtime, out bool Converged, out int NoOfIter, out int DOFs) {
            using (var tr = new FuncTrace()) {
                mintime = double.MaxValue;
                maxtime = 0;
                Converged = false;
                NoOfIter = int.MaxValue;
                DOFs = 0;

                AggregationGridBasis[][] XAggB;
                using (new BlockTrace("Aggregation_basis_init", tr)) {
                    XAggB = AggregationGridBasis.CreateSequence(base.MultigridSequence, u.Mapping.BasisS);
                }
                XAggB.UpdateXdgAggregationBasis(this.Op_Agglomeration);
  
                var MassMatrix = this.Op_mass.GetMassMatrix(this.u.Mapping, new double[] { 1.0 }, false, this.LsTrk.SpeciesIdS.ToArray());
                double[] _RHSvec = this.GetRHS();

                

                Stopwatch stw = new Stopwatch();
                stw.Reset();
                stw.Start();

                Console.WriteLine("Setting up multigrid operator...");

                int p = this.u.Basis.Degree;
                var MultigridOp = new MultigridOperator(XAggB, this.u.Mapping,
                    this.Op_Matrix,
                    this.Op_mass.GetMassMatrix(new UnsetteledCoordinateMapping(this.u.Basis), false),
                    OpConfig);
                Assert.True(MultigridOp!=null);

                int L = MultigridOp.Mapping.LocalLength;
                DOFs = MultigridOp.Mapping.TotalLength;

                double[] RHSvec = new double[L];
                MultigridOp.TransformRhsInto(_RHSvec, RHSvec);


                ISolverSmootherTemplate exsolver;

                SolverFactory SF = new SolverFactory(this.Control.NonLinearSolver, this.Control.LinearSolver);
                Action<int, double[],double[],MultigridOperator>[] Callbacks = { CustomItCallback};
                SF.GenerateLinear(out exsolver, MultigridSequence, OpConfig, Callbacks);



                using (new BlockTrace("Solver_Init", tr)) {
                    exsolver.Init(MultigridOp);
                }
                /*
                string filename = "XdgPoisson" + this.Grid.SpatialDimension + "p" + this.u.Basis.Degree + "R" + this.Grid.CellPartitioning.TotalLength;
                MultigridOp.OperatorMatrix.SaveToTextFileSparse(filename + ".txt");
                RHSvec.SaveToTextFile(filename + "_rhs.txt");

                var uEx = this.u.CloneAs();
                Op_Agglomeration.ClearAgglomerated(uEx.Mapping);
                var CO = new ConvergenceObserver(MultigridOp, MassMatrix, uEx.CoordinateVector.ToArray());
                uEx = null;
                CO.TecplotOut = "PoissonConvergence";
                //CO.PlotDecomposition(this.u.CoordinateVector.ToArray());

                if (exsolver is ISolverWithCallback) {
                    ((ISolverWithCallback)exsolver).IterationCallback = CO.IterationCallback;
                }
                //*/
                XDGField u2 = u.CloneAs();
                using (new BlockTrace("Solver_Run", tr)) {
                    // use solver (on XDG-field 'u2').
                    u2.Clear();
                    MultigridOp.UseSolver(exsolver, u2.CoordinateVector, _RHSvec);
                    Console.WriteLine("Solver: {0}, converged? {1}, {2} iterations.", exsolver.GetType().Name, exsolver.Converged, exsolver.ThisLevelIterations);
                    this.Op_Agglomeration.Extrapolate(u2.Mapping);
                    Assert.IsTrue(exsolver.Converged, "Iterative solver did not converge.");
                }
                stw.Stop();
                mintime = Math.Min(stw.Elapsed.TotalSeconds, mintime);
                maxtime = Math.Max(stw.Elapsed.TotalSeconds, maxtime);
                Converged = exsolver.Converged;
                NoOfIter = exsolver.ThisLevelIterations;

                // compute error between reference solution and multigrid solver
                XDGField ErrField = u2.CloneAs();
                ErrField.Acc(-1.0, u);
                double ERR = ErrField.L2Norm();
                double RelERR = ERR / u.L2Norm();
                Assert.LessOrEqual(RelERR, 1.0e-6, "Result from iterative solver above threshold.");
           
            }
        }

        ISparseSolver ReferenceSolver;
        private void ReferenceSolve() {
            double[] RHSvec;
            if (ReferenceSolver == null) {

                ReferenceSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

                var EqSys = this.Op_Matrix.ToMsrMatrix();
                for (int iRow = EqSys.RowPartitioning.i0; iRow < EqSys.RowPartitioning.iE; iRow++) {
                    if (EqSys.GetNoOfNonZerosPerRow(iRow) <= 0)
                        EqSys[iRow, iRow] = 1.0;
                }
                ReferenceSolver.DefineMatrix(EqSys);
            }

            RHSvec = GetRHS();
            ReferenceSolver.Solve(this.u.CoordinateVector, RHSvec);

            this.Op_Agglomeration.Extrapolate(this.u.Mapping);

        }

        private double[] GetRHS() {
            double[] RHSvec;
            BlockMsrMatrix MassMatrix;
            RHSvec = this.Op_Affine.CloneAs();
            RHSvec.ScaleV(-1.0);
            MassMatrix = this.Op_mass.GetMassMatrix(this.u.Mapping, new double[] { 1.0 }, false, this.LsTrk.SpeciesIdS.ToArray());
            MassMatrix.SpMV(1.0, this.rhs.CoordinateVector, 1.0, RHSvec);
            return RHSvec;
        }

        private void ConsistencyTest() {

            // consistency test on the original matrix
            // -----------------------------------------------------

            this.residual.Clear();
            double[] RHSvec = this.GetRHS();
            this.residual.CoordinateVector.SetV(RHSvec, 1.0);
            this.Op_Matrix.SpMV(-1.0, this.u.CoordinateVector, 1.0, residual.CoordinateVector);

            double residual_L2Norm = this.residual.L2Norm();
            Console.WriteLine("Residual norm: " + residual_L2Norm);

            Assert.LessOrEqual(residual_L2Norm, 1.0e-8);



            // consistency test on the multigrid 
            // --------------------------------------------------------------------

            AggregationGridBasis[][] XAggB = AggregationGridBasis.CreateSequence(base.MultigridSequence, u.Mapping.BasisS);
            XAggB.UpdateXdgAggregationBasis(this.Op_Agglomeration);


            int p = this.u.Basis.Degree;
            var MultigridOp = new MultigridOperator(XAggB, this.u.Mapping,
                this.Op_Matrix,
                this.Op_mass.GetMassMatrix(new UnsetteledCoordinateMapping(this.u.Basis), false),
                new MultigridOperator.ChangeOfBasisConfig[][] {
                    new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] { 0 }, mode = MultigridOperator.Mode.Eye, Degree = u.Basis.Degree }
                    }
                });

            double[] mgSolVec = new double[MultigridOp.Mapping.LocalLength];
            double[] mgRhsVec = new double[MultigridOp.Mapping.LocalLength];
            MultigridOp.TransformSolInto(this.u.CoordinateVector, mgSolVec);
            MultigridOp.TransformRhsInto(RHSvec, mgRhsVec);

            MgConsistencyTestRec(MultigridOp, mgSolVec, mgRhsVec);

            // 
            /*
            {
                int Jagg1 = MgSeq[1].NoOfAggregateCells;
                MultigridOperator MgOp0 = MultigridOp;
                MultigridOperator MgOp1 = MultigridOp.CoarserLevel;
                MultigridMapping Map0 = MgOp0.Mapping;
                MultigridMapping Map1 = MgOp1.Mapping;


                double[] V0 = new double[MgOp0.Mapping.LocalLength];
                double[] V1 = new double[MgOp1.Mapping.LocalLength];

                for(int j = 0; j < Jagg1; j++) {
                    int idx = Map1.LocalUniqueIndex(0, j, 0);
                    V1[idx] = j;
                }

                MgOp1.Prolongate(1.0, V0, 0.0, V1);

                XDGField Marker = new XDGField(this.u.Basis, "Tracker");

                MgOp0.TransformSolFrom(Marker.CoordinatesAsVector, V0);
                this.Op_Agglomeration.Extrapolate(Marker.CoordinatesAsVector, Marker.Mapping);

                Tecplot.PlotFields(new DGField[] { Marker }, "Tracker", "Tracker", 0.0, 5);

            }
             */
        }

        private static void MgConsistencyTestRec(MultigridOperator mgOp, double[] mgSolVec, double[] mgRhsVec) {

            double[] mgResidual = mgRhsVec.CloneAs();
            mgOp.OperatorMatrix.SpMV(1.0, mgSolVec, -1.0, mgResidual);
            double scale = 1.0 / (mgRhsVec.L2Norm() + mgSolVec.L2Norm());

            int DOFs = mgOp.Mapping.TotalLength;
            Debug.Assert(DOFs == mgOp.OperatorMatrix.NoOfRows);
            Debug.Assert(DOFs == mgOp.OperatorMatrix.NoOfCols);


            double mgResidual_l2Norm = mgResidual.L2Norm();
            Console.WriteLine("Multigrid Residual norm, level {0}: {1} ({2} DOFs)", mgOp.LevelIndex, mgResidual_l2Norm, DOFs);
            Assert.LessOrEqual(scale * mgResidual_l2Norm, 1.0e-8);

            if (mgOp.CoarserLevel != null) {
                double[] mgCoarseSolVec = new double[mgOp.CoarserLevel.Mapping.LocalLength];
                double[] mgCoarseRhsVec = new double[mgOp.CoarserLevel.Mapping.LocalLength];

                mgOp.CoarserLevel.Restrict(mgSolVec, mgCoarseSolVec);
                mgOp.CoarserLevel.Restrict(mgRhsVec, mgCoarseRhsVec);

                MgConsistencyTestRec(mgOp.CoarserLevel, mgCoarseSolVec, mgCoarseRhsVec);
            }
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int susamp) {
            var Fields = ArrayTools.Cat<DGField>(this.GradientU, this.Phi, this.u, this.rhs, this.residual, this.uEx, this.uErr);
            if (this.MGColoring != null && this.MGColoring.Length > 0)
                Fields = ArrayTools.Cat<DGField>(Fields, this.MGColoring);
            Tecplot.PlotFields(Fields, "XPoisson" + timestepNo.ToString(), physTime, susamp);
            Tecplot.PlotFields(Fields, "grid" + timestepNo.ToString(), physTime, 0);
        }

    }
}
