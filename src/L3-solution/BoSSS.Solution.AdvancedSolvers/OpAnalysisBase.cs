﻿using System;
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
            localRHS = RHS.CloneAs();

            // create the Dummy XDG aggregation basis
            var baseGrid = Mapping.GridDat;
            var mgSeq = Foundation.Grid.Aggregation.CoarseningAlgorithms.CreateSequence(baseGrid, 1);
            AggregationGridBasis[][] XAggB = AggregationGridBasis.CreateSequence(mgSeq, Mapping.BasisS);

            //
            XAggB.UpdateXdgAggregationBasis(CurrentAgglomeration);

            // create multigrid operator
            m_MultigridOp = new MultigridOperator(XAggB, Mapping,
                m_OpMtx,
                _mass,
                OpConfig,
                abstractOperator.DomainVar.Select(varName => abstractOperator.FreeMeanValue[varName]).ToArray());
        }


        LevelSetTracker m_LsTrk;

        UnsetteledCoordinateMapping m_map;
        int[] _VarGroup;
        
        BlockMsrMatrix m_OpMtx;
        double[] localRHS;
        
        MultigridOperator m_MultigridOp;


        /// <summary>
        /// user-defined indices of depended variables, if not the full matrix should be analyzed, e.g. 0 = u_x, 1=u_y, 2=u_z, 3=p ...
        /// </summary>
        public int[] VarGroup{
            get{
                return _VarGroup;
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
            double[] Eigenval_Write = Eigenval();

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
            Console.WriteLine("maximal eigenvalue: {0}", Eigenval_Write[0]);
            Console.WriteLine("minimal eigenvalue: {0}", Eigenval_Write[1]);

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
        public void rankAnalysis(BlockMsrMatrix OpMatrix, double[] RHS){

            int RHSlen = this.localRHS.Length;
            Debug.Assert(RHSlen == m_OpMtx.RowPartitioning.LocalLength);

            MultidimensionalArray outputArray = MultidimensionalArray.Create(2, 1); // The two rank values

            //At this point OpMatrix and RHS are local, they are collected within bmc on proc rank==0
            using (var bmc = new BatchmodeConnector()){
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
                if (rnkAugmentedMtx > rnkMtx)
                {
                    //throw new Exception("The rank of the augmented matrix shouldn't be greater than the one of the original matrix!!"); 
                    Console.WriteLine("======================================================");
                    Console.WriteLine("WARNING!!!!!!! The rank of the augmented matrix shouldn't be greater than the one of the original matrix!!");
                    Console.WriteLine("This means that the system does not have a solution!");
                }

                if (rnkAugmentedMtx == rnkMtx)
                {
                    Console.WriteLine("The system has at least a solution");
                }

                //RHS and OpMatrix will be collected in Bmc, so total length has to be considered for RHS: RHS.length will lead to errors in parallel execution
                if (rnkMtx < RHSlen)
                {
                    Console.WriteLine("The rank of the matrix is smaller than the number of variables. There are {0} free parameters", (RHS.Length - rnkMtx));
                }

                else if (rnkMtx == RHSlen)
                {
                    Console.WriteLine("The system has a unique solution :) ");
                }

                else
                {
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
        /// Note: 
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
                FullSel.VariableSelector(this.VarGroup);

                var InnerSel = new SubBlockSelector(m_MultigridOp.Mapping);
                InnerSel.VariableSelector(this.VarGroup);
                InnerSel.CellSelector(InnerCellsMask);


                // MUMPS condition number
                // ======================

                double condestFullMUMPS = (new BlockMask(FullSel)).GetSubBlockMatrix(Mtx).Condest_MUMPS();
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
            FullSel.VariableSelector(this.VarGroup);

            var InnerSel = new SubBlockSelector(m_MultigridOp.Mapping);
            InnerSel.VariableSelector(this.VarGroup);
            InnerSel.CellSelector(InnerCellsMask);
            
            // Matlab
            // ======

            double[] Full_0Vars = (new BlockMask(FullSel)).GlobalIndices.Select(i => i + 1.0).ToArray();
            double[] Inner_0Vars = (new BlockMask(InnerSel)).GlobalIndices.Select(i => i + 1.0).ToArray();

            MultidimensionalArray output = MultidimensionalArray.Create(2, 1);
            string[] names = new string[] { "Full_0Vars", "Inner_0Vars" };

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


        /// <summary>
        /// returns the condition number of the full matrix
        /// </summary>
        public double CondNumMatlab() {

            int[] DepVars = this.VarGroup;
            var grd = m_map.GridDat;
            int NoOfCells = grd.Grid.NumberOfCells;
            int NoOfBdryCells = grd.GetBoundaryCells().NoOfItemsLocally_WithExternal;


            var Mtx = m_MultigridOp.OperatorMatrix;




            // Blocks and selectors 
            // ====================
            var InnerCellsMask = grd.GetBoundaryCells().Complement();

            var FullSel = new SubBlockSelector(m_MultigridOp.Mapping);
            FullSel.VariableSelector(this.VarGroup);

            
            // Matlab
            // ======

            double[] Full_0Vars = (new BlockMask(FullSel)).GlobalIndices.Select(i => i + 1.0).ToArray();
            
            MultidimensionalArray output = MultidimensionalArray.Create(2, 1);
            string[] names = new string[] { "Full_0Vars", "Inner_0Vars" };

            using(BatchmodeConnector bmc = new BatchmodeConnector()) {

                // if Octave should be used instead of Matlab....
                // BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                bmc.PutSparseMatrix(Mtx, "FullMatrix");

                bmc.PutVector(Full_0Vars, "Full_0Vars");

                bmc.Cmd("output = ones(2,1);");

                bmc.Cmd("output(1) = condest(FullMatrix(Full_0Vars,Full_0Vars));");
               

                bmc.GetMatrix(output, "output");
                bmc.Execute(false);

                double condestFull = output[0, 0];
                Debug.Assert(condestFull.MPIEquals(), "value does not match on procs");
                

                Console.WriteLine($"MATLAB condition number: {condestFull:0.###e-00}");

                
                return condestFull;
            }
            
        }


        /// <summary>
        /// Test if the matrix is symmetric positive definite
        /// </summary>
        /// <returns>bool array res=[symmetry, positive definit]</returns>
        public bool[] Symmetry(){

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
            if (SymmDev < 1e-5)
                sym = true;

            res[0] = sym;


            // positive definite test by Cholesky decomposition
            var FullyPopulatedMatrix = m_OpMtx.ToFullMatrixOnProc0();

            bool posDef = true;
            // only proc 0 gets info so the following is executed exclusively on rank 0
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                try {
                    FullyPopulatedMatrix.Cholesky();
                } catch (ArithmeticException) {
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
        public double[] Eigenval(){

            double[] eigenvalues = new double[2];
            MultidimensionalArray eigs = MultidimensionalArray.Create(1, 2);
            MultidimensionalArray output = MultidimensionalArray.Create(2, 1);

            int[] DepVars = this.VarGroup;
            double[] DepVars_subvec = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();

            using (BatchmodeConnector bmc = new BatchmodeConnector()){

                // if Octave should be used instead of Matlab....
                //BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                bmc.PutSparseMatrix(m_OpMtx, "FullMatrix");
                bmc.PutVector(DepVars_subvec, "DepVars_subvec");
                bmc.Cmd("output = zeros(2,1)");
                bmc.Cmd("output(1) = eigs(FullMatrix(DepVars_subvec,DepVars_subvec),1,'lm');");
                bmc.Cmd("output(2) = eigs(FullMatrix(DepVars_subvec,DepVars_subvec),1,'sm');");
                bmc.GetMatrix(output, "output");
                bmc.Execute(false);     
            }

            double[] myeigs = new double[] { output[0, 0], output[1, 0] };
            Debug.Assert(output[0, 0].MPIEquals(), "value does not match on procs");
            Debug.Assert(output[1, 0].MPIEquals(), "value does not match on procs");
            return myeigs;
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
            Sel.VariableSelector(this.VarGroup);
           
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

                    var LocBlk = grd.GetCellNeighboursViaEdges(j).Select(t => t.Item1).ToList();
                    LocBlk.Add(j);
                    for(int i = 0; i < LocBlk.Count; i++) {
                        if(LocBlk[i] >= J) {
                            LocBlk.RemoveAt(i);
                            i--;
                        }
                    }

                    var Sel = new SubBlockSelector(m_MultigridOp.Mapping);
                    Sel.VariableSelector(this.VarGroup);
                    Sel.CellSelector(LocBlk, global: false);
                    var Mask = new BlockMask(Sel);

                    MultidimensionalArray[,] Blocks = Mask.GetFullSubBlocks(Mtx, ignoreSpecCoupling: false, ignoreVarCoupling: false);

                    MultidimensionalArray FullBlock = Blocks.Cat();

                    BCN[j] = FullBlock.Cond('I');

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


        string VarNames {
            get {
                int[] vgs = this.VarGroup;
                if(vgs.Length <= 1) {
                    return "Var0";
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

        /// <summary>
        /// Various condition numbers, organized in a dictionary to create a regression over multiple meshes
        /// </summary>
        public IDictionary<string, double> GetNamedProperties() {
            using(new FuncTrace()) {
                var Ret = new Dictionary<string, double>();

                var grd = m_map.GridDat;

                // global condition numbers
                // ========================
                //double CondNo = this.CondNumMUMPS();
                double CondNo = this.CondNumMatlab(); // matlab seems to be more reliable
                Ret.Add("TotCondNo-" + VarNames, CondNo);

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
                return Ret;
            }
        }
    }
}
