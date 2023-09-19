using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Foundation.Grid;
using MPI.Wrappers;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Solution.AdvancedSolvers;
using NUnit.Framework;
using ilPSP.Tracing;
using System.IO;

namespace BoSSS.Solution.AdvancedSolvers.Testing {

    /// <summary>
    /// Numerical testing of coupled operators
    /// </summary>
    public class OpAnalysisBase {


        /// <summary>
        /// Constructor for DG solvers
        /// </summary>
        public OpAnalysisBase(BlockMsrMatrix Mtx, double[] RHS, UnsetteledCoordinateMapping Mapping, IEnumerable<MultigridOperator.ChangeOfBasisConfig[]> OpConfig, ISpatialOperator abstractOperator) 
            : this(null, Mtx, RHS, Mapping, null, null, OpConfig, abstractOperator) //
        {

        }

       


        /// <summary>
        /// Constructor for XDG solvers
        /// </summary>
        public OpAnalysisBase(LevelSetTracker LsTrk, BlockMsrMatrix Mtx, double[] RHS, UnsetteledCoordinateMapping Mapping, MultiphaseCellAgglomerator CurrentAgglomeration, BlockMsrMatrix _mass, IEnumerable<MultigridOperator.ChangeOfBasisConfig[]> OpConfig, ISpatialOperator abstractOperator) {


            //int RHSlen = Mapping.TotalLength;
            m_map = Mapping; // mapping
            VarGroup = Mapping.BasisS.Count.ForLoop(i => i); //default: all dependent variables are included in operator matrix

            m_LsTrk = LsTrk;

            m_OpMtx = Mtx.CloneAs();
            m_Mass = _mass?.CloneAs();
            localRHS = RHS.CloneAs();

            // create the Dummy XDG aggregation basis
            var baseGrid = Mapping.GridDat;
            var mgSeq = Foundation.Grid.Aggregation.CoarseningAlgorithms.CreateSequence(baseGrid, 1);
            m_XAggB = AggregationGridBasis.CreateSequence(mgSeq, Mapping.BasisS);

            //
            m_XAggB.UpdateXdgAggregationBasis(CurrentAgglomeration);

            // create multigrid operator
            m_MultigridOp = new MultigridOperator(m_XAggB, Mapping,
                m_OpMtx,
                m_Mass,
                OpConfig,
                abstractOperator);
        }

        AggregationGridBasis[][] m_XAggB;

        LevelSetTracker m_LsTrk;

        UnsetteledCoordinateMapping m_map;
        int[] _VarGroup;
        
        BlockMsrMatrix m_Mass;
        
        BlockMsrMatrix m_OpMtx; // un-treated operator matrix
        double[] localRHS;
        
        MultigridOperator m_MultigridOp;

        /// <summary>
        /// Les operaeur
        /// </summary>
        public MultigridOperator MultigridOp {
            get {
                return m_MultigridOp;
            }
        }

        private bool calculateStencils = false;


        /// <summary>
        /// user-defined indices of depended variables, if not the full matrix should be analyzed, e.g. 0 = u_x, 1=u_y, 2=u_z, 3=p ...
        /// </summary>
        public int[] VarGroup {
            get{
                return _VarGroup.CloneAs();
            }
            set{
                _VarGroup = value;
                if (_VarGroup.Distinct().Count() != _VarGroup.Count()){
                    throw new ArithmeticException("duplicate found!");
                }

                foreach (int k in _VarGroup){
                    if (k >= m_map.BasisS.Count){
                        throw new ArithmeticException("dependent variables out of range!!");
                    }
                    else if (k < 0){
                        throw new ArithmeticException("dependent variables out of range!!");
                    }
                }
            }
        }

        /// <summary>
        /// Computes the minimal Eigenvalue and related Eigenvector using PARDISO
        /// </summary>
        public (double lambdaMin, double[] V) MinimalEigen() {
            var Mtx = this.PrecondOpMatrix;
            int L = Mtx.RowPartitioning.LocalLength;

            // extract sub-matrix
            var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
            FullSel.SetVariableSelector(this.VarGroup);
            var mask = new BlockMask(FullSel);
            var Part = mask.GetSubBlockMatrix(Mtx, Mtx.MPI_Comm);

            var bla = Part.MinimalEigen();

            double[] Vret = new double[L];
            mask.AccSubVec(bla.V, Vret);

            return (bla.lambdaMin, Vret);
        }

        public (double lambdaMax, double[] V) MaximalEigen() {
            var Mtx = this.m_OpMtx;
            int L = Mtx.RowPartitioning.LocalLength;

            // extract sub-matrix
            var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
            FullSel.SetVariableSelector(this.VarGroup);
            var mask = new BlockMask(FullSel);
            var Part = mask.GetSubBlockMatrix(Mtx, Mtx.MPI_Comm);

            var bla = Part.MaximalEigen();

            double[] Vret = new double[L];
            mask.AccSubVec(bla.V, Vret);

            return (bla.lambdaMax, Vret);
        }

        /// <summary>
        /// fully analyses the matrix with
        /// - condition number,
        /// - symmetric positive definiteness,
        /// - maximum and minimum eigenvalues and
        /// - the existence of a unique solution
        /// </summary>
        public void Analyse() {

            Console.WriteLine("Starting matrix Analysis");
            Console.WriteLine("Calculating condition number");
            double CondNum_Write = CondNumMUMPS();
            Console.WriteLine("Doing symmetry test");
            bool[] Symmetry_Write = Symmetry();
            Console.WriteLine("Doing eigenvalues test");
            var Eigenval_Write = Eigenval();

            Console.WriteLine("");
            Console.WriteLine("==================================================================");
            Console.WriteLine("Log of the analysis");
            Console.WriteLine("==================================================================");
            Console.WriteLine("Condition number:");
            Console.WriteLine("full matrix: {0}", CondNum_Write);
            Console.WriteLine("==================================================================");
            Console.WriteLine("Symmetry and positive definiteness:");
            Console.WriteLine("is symmetric: {0}", Symmetry_Write[0]);
            Console.WriteLine("is positive definite: {0}", Symmetry_Write[1]);
            Console.WriteLine("==================================================================");
            Console.WriteLine("Eigenvalues:");
            Console.WriteLine("maximal eigenvalue: {0}", Eigenval_Write.maxEigen);
            Console.WriteLine("minimal eigenvalue: {0}", Eigenval_Write.minEigen);

            rankAnalysis(m_OpMtx, localRHS);

            Console.WriteLine("");

        }


        /// <summary>
        /// According to the Rouché-Capelli theorem, the system is inconsistent if rank(augMatrix) > rank(Matrix). 
        /// If rank(augMatrix) == rank(Matrix), the system has at least one solution.
        /// Additionally, if the rank is equal to the number of variables, the solution of the system is unique.
        /// Remark: This requires the whole RHS not only local!!!
        /// </summary>
        /// <param name="OpMatrix"></param>
        /// <param name="RHS"></param>
        public void rankAnalysis(BlockMsrMatrix OpMatrix, double[] RHS) {

            int RHSlen = this.localRHS.Length;
            Debug.Assert(RHSlen == m_OpMtx.RowPartitioning.LocalLength);

            MultidimensionalArray outputArray = MultidimensionalArray.Create(2, 1); // The two rank values

            //At this point OpMatrix and RHS are local, they are collected within bmc on proc rank==0
            using(var bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(OpMatrix, "OpMatrix");
                bmc.PutVector(RHS, "RHS");
                bmc.Cmd("output = zeros(2,1)"); // First value is rank(OpMatrix), second the rank of the augmented matrix = rank([Matrix|RHS])
                bmc.Cmd("");
                bmc.Cmd("fullMtx = full(OpMatrix);");
                bmc.Cmd("augmentedMtx = [fullMtx RHS];");
                bmc.Cmd("output(1) = rank(fullMtx)");
                bmc.Cmd("output(2) = rank(augmentedMtx)");
                bmc.GetMatrix(outputArray, "output");
                bmc.Execute(false);
            }


            double[] output = new double[2];
            output[0] = outputArray[0, 0]; //Rank matrix 
            output[1] = outputArray[1, 0]; //Rank augmented matrix ( [Matrix|RHS])

            double rnkMtx = output[0];
            double rnkAugmentedMtx = output[1];

            // Some tests
            Console.WriteLine("==================================================================");

            //Output
            {
                Console.WriteLine("Results of rank analysis:");
                if(rnkAugmentedMtx > rnkMtx) {
                    //throw new Exception("The rank of the augmented matrix shouldn't be greater than the one of the original matrix!!"); 
                    Console.WriteLine("======================================================");
                    Console.WriteLine("WARNING!!!!!!! The rank of the augmented matrix shouldn't be greater than the one of the original matrix!!");
                    Console.WriteLine("This means that the system does not have a solution!");
                }

                if(rnkAugmentedMtx == rnkMtx) {
                    Console.WriteLine("The system has at least a solution");
                }

                //RHS and OpMatrix will be collected in Bmc, so total length has to be considered for RHS: RHS.length will lead to errors in parallel execution
                if(rnkMtx < RHSlen) {
                    Console.WriteLine("The rank of the matrix is smaller than the number of variables. There are {0} free parameters", (RHS.Length - rnkMtx));
                } else if(rnkMtx == RHSlen) {
                    Console.WriteLine("The system has a unique solution :) ");
                } else {
                    throw new Exception("what? should not happen");
                }

                Console.WriteLine("Rank of the matrix : {0} \n" + "Rank of the augmented matrix : {1} \n" + "Number of variables: {2}", output[0], output[1], RHS.Length);
                Console.WriteLine("==================================================================");
            }

            Debug.Assert(output[0].MPIEquals(), "value does not match on procs");
            Debug.Assert(output[1].MPIEquals(), "value does not match on procs");
        }


        /// <summary>
        /// returns the condition number of the full matrix using MUMPS;
        /// From manual it is unclear in which norm the cond num is calculated
        /// results from parallel and single execution differ!
        /// </summary>
        public double CondNumMUMPS() {
            using(new FuncTrace()) {
                int[] DepVars = this.VarGroup;
                var grd = m_map.GridDat;
                int NoOfCells = grd.Grid.NumberOfCells;
                int NoOfBdryCells = grd.GetBoundaryCells().NoOfItemsLocally_WithExternal;


                var Mtx = m_MultigridOp.OperatorMatrix;
                
                // Blocks and selectors 
                // ====================
                var InnerCellsMask = grd.GetBoundaryCells().Complement();

                var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
                FullSel.SetVariableSelector(this.VarGroup);

                var InnerSel = new SubBlockSelector(m_MultigridOp.Mapping);
                InnerSel.SetVariableSelector(this.VarGroup);
                InnerSel.CellSelector(InnerCellsMask);


                // MUMPS condition number
                // ======================

                double condestFullMUMPS = (new BlockMask(FullSel)).GetSubBlockMatrix_MpiSelf(Mtx).Condest_MUMPS();
                //double condestInnerMUMPS = 1.0;

                //if(InnerCellsMask.NoOfItemsLocally.MPISum() > 0) {
                //    condestInnerMUMPS = (new BlockMask(InnerSel)).GetSubBlockMatrix(Mtx).Condest_MUMPS();
                //}


                //return new[] { condestFullMUMPS, condestInnerMUMPS };
                return condestFullMUMPS;
            }
        }

        /// <summary>
        /// returns the condition number of the full matrix and the inner matrix without boundary terms
        /// </summary>
        /// <returns>
        /// [ConditionNumberFullOp, ConditionNumberInnerOp]
        /// </returns>
        public (double,double) CondNumMatlab2() {

            int[] DepVars = this.VarGroup;
            var grd = m_map.GridDat;
            int NoOfCells = grd.Grid.NumberOfCells;
            int NoOfBdryCells = grd.GetBoundaryCells().NoOfItemsLocally_WithExternal;


            var Mtx = m_MultigridOp.OperatorMatrix;




            // Blocks and selectors 
            // ====================
            var InnerCellsMask = grd.GetBoundaryCells().Complement();

            var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
            FullSel.SetVariableSelector(this.VarGroup);

            var InnerSel = new SubBlockSelector(m_MultigridOp.Mapping);
            InnerSel.SetVariableSelector(this.VarGroup);
            InnerSel.CellSelector(InnerCellsMask);
            
            // Matlab
            // ======

            double[] Full_0Vars = (new BlockMask(FullSel)).GlobalIndices.Select(i => i + 1.0).ToArray();
            double[] Inner_0Vars = (new BlockMask(InnerSel)).GlobalIndices.Select(i => i + 1.0).ToArray();

            MultidimensionalArray output = MultidimensionalArray.Create(2, 1);
            //string[] names = new string[] { "Full_0Vars", "Inner_0Vars" };

            using(BatchmodeConnector bmc = new BatchmodeConnector()) {

                // if Octave should be used instead of Matlab....
                // BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                bmc.PutSparseMatrix(Mtx, "FullMatrix");

                if(Inner_0Vars.Length > 0)
                    bmc.PutVector(Inner_0Vars, "Inner_0Vars");
                bmc.PutVector(Full_0Vars, "Full_0Vars");

                bmc.Cmd("output = ones(2,1);");

                bmc.Cmd("output(1) = condest(FullMatrix(Full_0Vars,Full_0Vars));");
                if(Inner_0Vars.Length > 0)
                    bmc.Cmd("output(2) = condest(FullMatrix(Inner_0Vars,Inner_0Vars));");

                bmc.GetMatrix(output, "output");
                bmc.Execute(false);

                double condestFull = output[0, 0];
                double condestInner = output[1, 0];
                Debug.Assert(condestFull.MPIEquals(), "value does not match on procs");
                Debug.Assert(condestInner.MPIEquals(), "value does not match on procs");


                Console.WriteLine($"MATLAB condition number: {condestFull:0.###e-00}");

                return (condestFull, condestInner);
            }
            
        }

        public double CondLAPACK() {
            var Mtx = m_OpMtx.ToFullMatrixOnProc0();
            double res = double.NaN;
            if (this.m_OpMtx.RowPartitioning.MpiRank == 0) {
                res=Mtx.Cond('I');
            }
            return res;
        }


        /// <summary>
        /// Creates a new <see cref="MultigridOperator"/>, where all variables in <see cref="VarGroup"/> are preconditioned according to <paramref name="altMode"/>.
        /// (All Variables **not** in <see cref="VarGroup"/> are not preconditioned, i.e. the option <see cref="MultigridOperator.Mode.Eye"/> is chosen.
        /// </summary>
        MultigridOperator GetAlternativeMgOp(MultigridOperator.Mode altMode) {
            
            // Determine the operator configuration for the alternative operator
            // =================================================================

            MultigridOperator.ChangeOfBasisConfig[][] AlternativeOpConfig = new MultigridOperator.ChangeOfBasisConfig[1][];
            {
                int NoOfVars = m_MultigridOp.Mapping.NoOfVariables;

                bool[] VarGroupMarker = new bool[NoOfVars];
                this.VarGroup.ForEach((int iVar) => VarGroupMarker[iVar] = true);
                int NoOfOtherVariables = NoOfVars - this.VarGroup.Length;
                if(NoOfOtherVariables < 0)
                    throw new ApplicationException();

                AlternativeOpConfig[0] = new MultigridOperator.ChangeOfBasisConfig[1 + NoOfOtherVariables];

                int GetDegree(int iVar) {
                    foreach(var conf in this.m_MultigridOp.Config) {
                        int L = conf.VarIndex.Length;
                        for(int l = 0; l < L; l++) {
                            if(conf.VarIndex[l] == iVar)
                                return conf.DegreeS[iVar];
                        }
                    }
                    throw new ArgumentException($"unable to find variable {iVar} in original operators change-of-basis configuration");
                }

                // set the alternative mode for the vairable group which should be analyzed
                AlternativeOpConfig[0][0] = new MultigridOperator.ChangeOfBasisConfig() {
                    VarIndex = this.VarGroup.CloneAs(),
                    DegreeS = this.VarGroup.Select(iVar => GetDegree(iVar)).ToArray(),
                    mode = altMode
                };

                // set "Eye" for all other variables
                int cnt = 0;
                for(int iVar = 0; iVar < NoOfVars; iVar++) {
                    if(VarGroupMarker[iVar] == false) {
                        AlternativeOpConfig[0][cnt] = new MultigridOperator.ChangeOfBasisConfig() {
                            VarIndex = new int[] { iVar },
                            DegreeS = new int[] { GetDegree(iVar) },
                            mode = MultigridOperator.Mode.Eye
                        };
                        cnt++;

                    }
                }
                if(cnt != NoOfOtherVariables)
                    throw new ApplicationException("error in algorithm");
            }

            // create the alternative operator
            // ===============================

            return new MultigridOperator(m_XAggB, m_MultigridOp.BaseGridProblemMapping,
                m_OpMtx.CloneAs(),
                m_Mass,
                AlternativeOpConfig,
                m_MultigridOp.AbstractOperator);

        }


        //public static BlockMsrMatrix DbeMatrix;
        //public static BlockMsrMatrix LidMatrix;



        public double MatrixStabilityTest() {

            // setup original and alternative operator
            // =======================================
            var OpOrg = m_MultigridOp;
            var OpAlt = GetAlternativeMgOp(MultigridOperator.Mode.LeftInverse_DiagBlock);

            // setup of RHS
            // ============
            double[] RhsOrgMap, RhsAltMap;
            int Lbase, Lmap;
            {
                Lbase = OpOrg.BaseGridProblemMapping.LocalLength;
                double[] RandomRhs = new double[Lbase];
                RandomRhs.FillRandom();

                if(OpOrg.Mapping.LocalLength != OpAlt.Mapping.LocalLength)
                    throw new ApplicationException("Error in algorithm.");
                Lmap = OpOrg.Mapping.LocalLength;


                RhsOrgMap = new double[Lmap];
                OpOrg.TransformRhsInto(RandomRhs, RhsOrgMap, true);

                RhsAltMap = new double[Lmap];
                OpAlt.TransformRhsInto(RandomRhs, RhsAltMap, true);
            }


            // setup Matrices
            // ==============
            var MtxOrg = OpOrg.OperatorMatrix;
            var MtxAlt = OpAlt.OperatorMatrix;

            /*
            {
                var Diff2 = DbeMatrix.CloneAs();
                Diff2.Acc(-1.0, MtxOrg);
                double diffDbeNorm = Diff2.InfNorm();
                Console.WriteLine("Dbe diff = " + diffDbeNorm);
            }

            {
                var Diff1 = LidMatrix.CloneAs();
                Diff1.Acc(-1.0, MtxAlt);
                double diffLidNorm = Diff1.InfNorm();
                Console.WriteLine("Lid diff = " + diffLidNorm);
            }
            */
            //var diff = MtxAlt.CloneAs();
            //diff.Acc(-1.0, MtxOrg);
            //double r2123 = diff.InfNorm();
            //Console.WriteLine("Difference: " + r2123);


            // extract blocks corresponding to 'VarGroup'
            // ==========================================

            var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
            FullSel.SetVariableSelector(this.VarGroup);
            var mask = new BlockMask(FullSel);
            long[] GidxS = mask.GlobalIndices;

            var MtxOrg_Block = MtxOrg.GetSubMatrix(GidxS, GidxS);
            var MtxAlt_Block = MtxAlt.GetSubMatrix(GidxS, GidxS);

            int LBlock = MtxAlt_Block._RowPartitioning.LocalLength;
            long i0 = MtxOrg._RowPartitioning.i0;

            double[] RhsOrgMap_Block = new double[LBlock], RhsAltMap_Block = new double[LBlock];
            RhsOrgMap_Block.AccVi64(1.0, RhsOrgMap, acc_index: default(long[]), b_index: GidxS, b_index_shift: -i0);
            RhsAltMap_Block.AccVi64(1.0, RhsAltMap, acc_index: default(long[]), b_index: GidxS, b_index_shift: -i0);

            /*
            {
                var Diff2 = DbeMatrix.CloneAs();
                Diff2.Acc(-1.0, MtxOrg_Block);
                double diffDbeNorm = Diff2.InfNorm();
                Console.WriteLine("Dbe diff = " + diffDbeNorm);
            }

            {
                var Diff1 = LidMatrix.CloneAs();
                Diff1.Acc(-1.0, MtxAlt_Block);
                double diffLidNorm = Diff1.InfNorm();
                Console.WriteLine("Lid diff = " + diffLidNorm);
            }
            */

            // solve original and alternative system
            // =====================================

            double[] Xorg_Block, Xalt_Block;
            {
                Xorg_Block = new double[LBlock];
                MtxOrg_Block.Solve_Direct(Xorg_Block, RhsOrgMap_Block);
                MtxOrg_Block.SaveToTextFileSparse("OrgBlock.txt");
                //DbeMatrix.Solve_Direct(Xorg_Block, RhsOrgMap_Block);

                Xalt_Block = new double[LBlock];
                MtxAlt_Block.Solve_Direct(Xalt_Block, RhsAltMap_Block);
                MtxAlt_Block.SaveToTextFileSparse("AltBlock.txt");
                //LidMatrix.Solve_Direct(Xalt_Block, RhsAltMap_Block);
            }

            // transform back int order to compare
            // ===================================
            double[] XorgMap = new double[Lmap], Xorg = new double[Lbase];
            XorgMap.AccVi64(1.0, Xorg_Block, acc_index: GidxS, b_index: default(long[]), acc_index_shift: -i0);
            OpOrg.TransformSolFrom(Xorg, XorgMap);


            double[] XaltMap = new double[Lmap], Xalt = new double[Lbase];
            XaltMap.AccVi64(1.0, Xalt_Block, acc_index: GidxS, b_index: default(long[]), acc_index_shift: -i0);
            OpAlt.TransformSolFrom(Xalt, XaltMap);

            // finally, compare!
            // =================

            double Ref = Math.Max(Xorg.MPI_L2Norm(), Xalt.MPI_L2Norm());
            double dist = Xorg.MPI_L2Dist(Xalt);
            double DisturbanceMeasure =  dist / Ref;
            Console.WriteLine("Matrix stability test disturbance measure: " + dist);

            if(DisturbanceMeasure > 1.0e-8) {
                Console.Error.WriteLine("!!!!!!!!!!!! Probably unstable discretization !!!!!!!!!!");
            }

            return DisturbanceMeasure;
        }

        /// <summary>
        /// The operator matrix after compactification (i.e. elimination of un-used DOFs),
        /// application of reference points (<see cref="ISpatialOperator.FreeMeanValue"/>),
        /// and block-preconditioning.
        /// </summary>
        public BlockMsrMatrix PrecondOpMatrix {
            get {
                return m_MultigridOp.OperatorMatrix.CloneAs();
            }
        }


        /// <summary>
        /// returns the condition number of the full matrix
        /// </summary>
        public double CondNumMatlab() {

            int[] DepVars = this.VarGroup;
            var grd = m_map.GridDat;
            int NoOfCells = grd.Grid.NumberOfCells;
            int NoOfBdryCells = grd.GetBoundaryCells().NoOfItemsLocally_WithExternal;


            var Mtx = m_MultigridOp.OperatorMatrix;

            //bool Comparison = this.VarGroup.SetEquals(new int[] { 1 });


            // Blocks and selectors 
            // ====================
            var InnerCellsMask = grd.GetBoundaryCells().Complement();

            var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
            FullSel.SetVariableSelector(this.VarGroup);

            //long J = grd.CellPartitioning.TotalLength;

            // Matlab
            // ======

            double[] Full_0Vars = (new BlockMask(FullSel)).GlobalIndices.Select(i => i + 1.0).ToArray();
            
            MultidimensionalArray output = MultidimensionalArray.Create(4, 1);
            //string[] names = new string[] { "Full_0Vars", "Inner_0Vars" };

            using(BatchmodeConnector bmc = new BatchmodeConnector()) {
                //if(Comparison) {
                //    string compPath = $"Mtx_ipPoisson-J{J}.txt";
                //    File.Copy(compPath, Path.Combine(bmc.WorkingDirectory.FullName, "comp.txt"), false);
                //}

                bmc.PutSparseMatrix(Mtx, "FullMatrix");

                bmc.PutVector(Full_0Vars, "Full_0Vars");

                bmc.Cmd("output = ones(4,1);");


                bmc.Cmd("output(1) = condest(FullMatrix(Full_0Vars,Full_0Vars));");

                //if(Comparison) {
                //    bmc.Cmd("compMtx = ReadMsr('comp.txt');");
                //    bmc.Cmd("output(2) = norm(compMtx - FullMatrix(Full_0Vars,Full_0Vars), inf);");
                //    bmc.Cmd("output(3) = condest(compMtx);");
                //}

                bmc.GetMatrix(output, "output");
                bmc.Execute(false);

                double condestFull = output[0, 0];
                Debug.Assert(condestFull.MPIEquals(), "value does not match on procs");
                

                Console.WriteLine($"MATLAB condition number vars {VarNames}: {condestFull:0.###e-00}");
                //if(Comparison) {
                //    Console.WriteLine($"MATLAB condition number check matrix vars {VarNames}: {output[2, 0]:0.###e-00}");
                //    Console.WriteLine($"MATLAB matrix check {VarNames}: {output[1, 0]:0.###e-00}");
                //}
                
                return condestFull;
            }
            
        }

        public double Cond2Matlab() {
            var Mtx = this.m_OpMtx;
            return Mtx.cond();
        }

        /// <summary>
        /// Test if the matrix is symmetric positive definite
        /// </summary>
        /// <returns>bool array res=[symmetry, positive definit]</returns>
        public bool[] Symmetry() {

            bool[] res = new bool[2];

            //extract submatrix for selected dependent variables
            int[] DepVars = this.VarGroup;
            double[] DepVars_subvec = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();

            //MsrMatrix OpMtxMSR = m_OpMtx.ToMsrMatrix();

            long[] SubMatrixIdx_Row = m_map.GetSubvectorIndices(false, DepVars);
            long[] SubMatrixIdx_Cols = m_map.GetSubvectorIndices(false, DepVars);
            int L = SubMatrixIdx_Row.Length;

            MsrMatrix SubOpMtx = new MsrMatrix(L, L, 1, 1);

            m_OpMtx.WriteSubMatrixTo(SubOpMtx, SubMatrixIdx_Row, default(long[]), SubMatrixIdx_Cols, default(long[]));


            // symmetry test by calculation of symmetry deviation
            bool sym = false;
            double SymmDev = m_OpMtx.SymmetryDeviation();
            if(SymmDev < 1e-5)
                sym = true;

            res[0] = sym;


            // positive definite test by Cholesky decomposition
            var FullyPopulatedMatrix = m_OpMtx.ToFullMatrixOnProc0();

            bool posDef = true;
            // only proc 0 gets info so the following is executed exclusively on rank 0
            if(ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                try {
                    FullyPopulatedMatrix.Cholesky();
                } catch(ArithmeticException) {
                    posDef = false;
                }
            }
            res[1] = MPIExtensions.MPIBroadcast(posDef, 0, ilPSP.Environment.MPIEnv.Mpi_comm);

            Debug.Assert(res[0].MPIEquals(), "value does not match on procs");
            Debug.Assert(res[1].MPIEquals(), "value does not match on procs");
            return res;
        }


        /// <summary>
        /// returns the maximum and minimum eigenvalues of the matrix
        /// </summary>
        /// <returns>Array myeigs =[MaximumEig, MinimumEig]</returns>
        public (double maxEigen, double minEigen) Eigenval() {
            
            var Mtx = m_OpMtx;

            double[] eigenvalues = new double[2];
            MultidimensionalArray eigs = MultidimensionalArray.Create(1, 2);
            MultidimensionalArray output = MultidimensionalArray.Create(2, 1);

            int[] DepVars = this.VarGroup;
            double[] DepVars_subvec = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();

            using(BatchmodeConnector bmc = new BatchmodeConnector()) {

                // if Octave should be used instead of Matlab....
                //BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                bmc.PutSparseMatrix(Mtx, "FullMatrix");
                bmc.PutVector(DepVars_subvec, "DepVars_subvec");
                bmc.Cmd("output = zeros(2,1)");
                bmc.Cmd("output(1) = eigs(FullMatrix(DepVars_subvec,DepVars_subvec),1,'lm');");
                bmc.Cmd("output(2) = eigs(FullMatrix(DepVars_subvec,DepVars_subvec),1,'sm');");
                bmc.GetMatrix(output, "output");
                bmc.Execute(false);
            }

            Debug.Assert(output[0, 0].MPIEquals(), "value does not match on procs");
            Debug.Assert(output[1, 0].MPIEquals(), "value does not match on procs");
            return (output[0, 0], output[1, 0]);
        }


        /// <summary>
        /// Local condition number for the block formed by each cell.
        /// </summary>
        /// <returns>
        /// one value per cell:
        /// - index: local cell index
        /// - content: condition number (one norm) of the local stencil
        /// </returns>
        public double[] BlockCondNumbers() {
            int J = m_map.LocalNoOfBlocks;
            Debug.Assert(J == m_map.GridDat.iLogicalCells.NoOfLocalUpdatedCells);
            
            var Mtx = m_MultigridOp.OperatorMatrix;
            Debug.Assert(Mtx._ColPartitioning.LocalNoOfBlocks == J);
            Debug.Assert(Mtx._RowPartitioning.LocalNoOfBlocks == J);


            var Sel = new SubBlockSelector(m_MultigridOp.Mapping);
            Sel.SetVariableSelector(this.VarGroup);
           
            var Mask = new BlockMask(Sel);

            MultidimensionalArray[] Blocks = Mask.GetDiagonalBlocks(Mtx, ignoreSpecCoupling: false, ignoreVarCoupling: false);
            Debug.Assert(Blocks.Length == J);

            double[] BCN = new double[J];
            for(int j = 0; j < J; j++) {
#if DEBUG
                int N = this.VarGroup.Sum(iVar => m_MultigridOp.Mapping.AggBasis[iVar].GetLength(j, m_MultigridOp.Degrees[iVar]));
                Debug.Assert(Blocks[j].NoOfCols == N);
                Debug.Assert(Blocks[j].NoOfRows == N);

                long i0 = m_MultigridOp.Mapping.GlobalUniqueIndex(this.VarGroup.Min(), j, 0);
                Debug.Assert(Mtx[i0, i0] == Blocks[j][0, 0]);
#endif

                BCN[j] = Blocks[j].Cond();
            }

            return BCN;
        }


        /// <summary>
        /// Local condition number formed by the block of each cell and its neighbors.
        /// </summary>
        /// <returns>
        /// one value per cell:
        /// - index: local cell index
        /// - content: condition number (one norm) of the local stencil
        /// </returns>
        public double[] StencilCondNumbers() {
            using(new FuncTrace()) {
                int J = m_map.LocalNoOfBlocks;
                Debug.Assert(J == m_map.GridDat.iLogicalCells.NoOfLocalUpdatedCells);

                var Mtx = m_MultigridOp.OperatorMatrix;
                Debug.Assert(Mtx._ColPartitioning.LocalNoOfBlocks == J);
                Debug.Assert(Mtx._RowPartitioning.LocalNoOfBlocks == J);

                var grd = m_MultigridOp.Mapping.AggGrid;

                double[] BCN = new double[J];

                for(int j = 0; j < J; j++) {
                    try { 
                    var LocBlk = grd.GetCellNeighboursViaEdges(j).Select(t => t.Item1).ToList();
                    LocBlk.Add(j);
                    for(int i = 0; i < LocBlk.Count; i++) {
                        if(LocBlk[i] >= J) {
                            LocBlk.RemoveAt(i);
                            i--;
                        }
                    }

                    var Sel = new SubBlockSelector(m_MultigridOp.Mapping);
                    Sel.SetVariableSelector(this.VarGroup);
                    Sel.CellSelector(LocBlk, global: false);
                    var Mask = new BlockMask(Sel);

                    MultidimensionalArray[,] Blocks = Mask.GetFullSubBlocks(Mtx, ignoreSpecCoupling: false, ignoreVarCoupling: false);

                    MultidimensionalArray FullBlock = Blocks.Cat();

                    BCN[j] = FullBlock.Cond('I');
                    } catch {
                        Console.WriteLine($"Empty block detected in stencil for j={j}");
                        BCN[j] = -999;
                    }
                }
                return BCN;
            }
        }

        /// <summary>
        /// Local condition number formed by the block of each cell and its neighbors, visualized as a degree-0 DG field.
        /// </summary>
        /// <returns>
        /// a DG field of degree 0, where the average value in each cell represents the respective local condition number
        /// </returns>
        public SinglePhaseField StencilCondNumbersV() {
            var grd = m_map.GridDat;

            SinglePhaseField ret = new SinglePhaseField(new Basis(grd, 0), "StencilCondNo-" + VarNames);

            double[] BCN = StencilCondNumbers();

            for(int j = 0; j < BCN.Length; j++)
                ret.SetMeanValue(j, BCN[j]);

            return ret;
        }

        /// <summary>
        /// Turn VarGroup into Names
        /// </summary>
        public string VarNames {
            get {
                int[] vgs = this.VarGroup;
                if(vgs.Length <= 1) {
                    return "Var" + vgs[0];
                } else {
                    Array.Sort(vgs);

                    string R = "Vars";
                    for(int i = 0; i < vgs.Length; i++) {
                        R += vgs[i];
                        if(i < vgs.Length - 1)
                            R += ".";
                    }
                    return R;
                }
            }
        }

        public bool CalculateStencils { get => calculateStencils; set => calculateStencils = value; }

        /// <summary>
        /// Various condition numbers, organized in a dictionary to create a regression over multiple meshes
        /// </summary>
        public IDictionary<string, double> GetNamedProperties() {
            using(new FuncTrace()) {
                var Ret = new Dictionary<string, double>();

                var grd = m_map.GridDat;

                // matrix stability
                // ================

                //double Stab = MatrixStabilityTest();
                //Ret.Add("MatrixStabilityMeasure-" + VarNames, Stab);

                // global condition numbers
                // ========================
                var stpw = new Stopwatch();
                stpw.Start();
                //double CondNo = this.CondNumMUMPS();
                double CondNo = this.CondNumMatlab(); // matlab seems to be more reliable
                //double CondNo = this.Cond2Matlab();
                //double CondNo = this.CondLAPACK();
                Ret.Add("TotCondNo-" + VarNames, CondNo);
                stpw.Stop();
                Console.WriteLine("- Calculated in " + stpw.Elapsed.TotalSeconds + " seconds");

                //var pair = this.Eigenval();
                //Ret.Add("maxEigenvalue", pair.maxEigen);
                //Ret.Add("minEigenvalue", pair.minEigen);
                //var mineig = this.MinimalEigen();
                //var maxeig = this.MaximalEigen();
                //Ret.Add("lambdaMin", mineig.lambdaMin);
                //Ret.Add("lambdaMax", maxeig.lambdaMax);
                if (CalculateStencils == true) { 
                    // block-wise condition numbers
                    // ============================
                    double[] bcn = this.StencilCondNumbers();

                    CellMask innerUncut, innerCut, bndyUncut, bndyCut;
                    if(m_LsTrk != null) {
                        // +++++++++
                        // using XDG
                        // +++++++++
                        innerUncut = grd.GetBoundaryCells().Complement().Except(m_LsTrk.Regions.GetCutCellMask());
                        innerCut = m_LsTrk.Regions.GetCutCellMask().Except(grd.GetBoundaryCells());

                        bndyUncut = grd.GetBoundaryCells().Except(m_LsTrk.Regions.GetCutCellMask());
                        bndyCut = grd.GetBoundaryCells().Intersect(m_LsTrk.Regions.GetCutCellMask());
                    } else {
                        // ++++++++
                        // using DG 
                        // ++++++++
                        innerUncut = grd.GetBoundaryCells().Complement();
                        innerCut = CellMask.GetEmptyMask(grd);

                        bndyUncut = grd.GetBoundaryCells();
                        bndyCut = CellMask.GetEmptyMask(grd);
                    }
                    double innerUncut_MaxCondNo = innerUncut.NoOfItemsLocally > 0 ? innerUncut.ItemEnum.Max(jCell => bcn[jCell]) : 1.0;
                    double innerCut_MaxCondNo = innerCut.NoOfItemsLocally > 0 ? innerCut.ItemEnum.Max(jCell => bcn[jCell]) : 1.0;
                    double bndyUncut_MaxCondNo = bndyUncut.NoOfItemsLocally > 0 ? bndyUncut.ItemEnum.Max(jCell => bcn[jCell]) : 1.0;
                    double bndyCut_MaxCondNo = bndyCut.NoOfItemsLocally > 0 ? bndyCut.ItemEnum.Max(jCell => bcn[jCell]) : 1.0;

                    innerUncut_MaxCondNo = innerUncut_MaxCondNo.MPIMax();
                    innerCut_MaxCondNo = innerCut_MaxCondNo.MPIMax();
                    bndyUncut_MaxCondNo = bndyUncut_MaxCondNo.MPIMax();
                    bndyCut_MaxCondNo = bndyCut_MaxCondNo.MPIMax();

                    Ret.Add("StencilCondNo-innerUncut-" + VarNames, innerUncut_MaxCondNo);
                    Ret.Add("StencilCondNo-bndyUncut-" + VarNames, bndyUncut_MaxCondNo);

                    if(m_LsTrk != null) {
                        Ret.Add("StencilCondNo-innerCut-" + VarNames, innerCut_MaxCondNo);
                        Ret.Add("StencilCondNo-bndyCut-" + VarNames, bndyCut_MaxCondNo);
                    }
                }
                return Ret;
            }
        }
    }
}
