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
using ilPSP.LinSolvers;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using BoSSS.Foundation;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using System.Collections;
using System.IO;
using BoSSS.Foundation.Grid.Classic;
using MPI.Wrappers;

namespace BoSSS.Solution.AdvancedSolvers {
    partial class MultigridOperator {

        /// <summary>
        /// configuration item for the change-of-basis in 
        /// </summary>
        public class ChangeOfBasisConfig : ICloneable {

            /// <summary>
            /// a list of variable indices onto which this configuration item applies.
            /// </summary>
            public int[] VarIndex;

            /// <summary>
            /// DG polynomial degree of respective variables on respective multigrid level
            /// </summary>
            public int[] DegreeS;

            ///// <summary>
            ///// if true, all species in one cell 
            ///// </summary>
            //public bool SpeciesSeparately;

            /// <summary>
            /// pre-conditioner mode
            /// </summary>
            public Mode mode;

            /// <summary>
            /// 
            /// </summary>
            public object Clone() {
                var r =new ChangeOfBasisConfig();
                r.mode = this.mode;
                r.DegreeS = this.DegreeS.CloneAs();
                r.VarIndex = this.VarIndex.CloneAs();
                return r;
            }
        }

        /// <summary>
        /// The type of change-of-basis 
        /// </summary>
        public enum Mode {

            /// <summary>
            /// no change of basis; trial and test function space are 
            /// kept as they are induced by the restriction operator
            /// </summary>
            Eye,

            /// <summary>
            /// Performs a change of the (X)DG basis so that the mass matrix 
            /// with respect to this basis
            /// becomes an identity.
            /// In non-XDG cases, this should be equivalent to <see cref="Eye"/>.
            /// </summary>
            IdMass,

            /// <summary>
            /// multiplies with the inverse of the mass matrix from the left
            /// </summary>
            LeftInverse_Mass,

            /// <summary>
            /// Symmetric equilibration of the diagonal matrix block;
            /// </summary>
            DiagBlockEquilib,

            /// <summary>
            /// Like <see cref="DiagBlockEquilib"/>, but only the 
            /// symmetric part of the diagonal block is used for equilibration; i.e. this works only for un-symmetric matrices.
            /// </summary>
            SymPart_DiagBlockEquilib,

            /// <summary>
            /// multiplies from the left with the inverse of the operator matrix diagonal block.
            /// </summary>
            LeftInverse_DiagBlock,


            /// <summary>
            /// Highly experimental option, blah blah blah.
            /// </summary>
            IdMass_DropIndefinite,

            /// <summary>
            /// Highly experimental option, blah blah blah.
            /// </summary>
            SymPart_DiagBlockEquilib_DropIndefinite,

            /// <summary>
            /// Schur complement for saddle-point systems
            /// </summary>
            SchurComplement
        }


        readonly ChangeOfBasisConfig[] m_Config;

        /// <summary>
        /// Multigrid Operator configuration in current level
        /// </summary>
        public ChangeOfBasisConfig[] Config {
            get {
                return m_Config.Select(c => c.CloneAs()).ToArray();
            }
        }

        /// <summary>
        /// the DG degrees on this level
        /// </summary>
        public int[] Degrees {
            get {
                VerifyConfig();
                int[] R = new int[this.BaseGridProblemMapping.BasisS.Count()];
                foreach (var c in m_Config) {
                    if (c.DegreeS.Length != c.VarIndex.Length) {
                        throw new ArgumentException("Length of Degrees must match number of Variables.");
                    }
                    for (int i = 0; i < c.VarIndex.Length; i++) {
                        int iVar = c.VarIndex[i];
                        int deg = c.DegreeS[i];
                        R[iVar] = deg;
                    }
                }

                if (AbstractOperator != null) {
                    if (!AbstractOperator.IsValidDomainDegreeCombination(R, R)) {
                        throw new ArgumentException($"DG degree combiation [{R.ToConcatString("", ", ", "")}] is reported to be illegal for DG operator");
                    }
                }

                return R;
            }
        }


        void VerifyConfig() {
            bool[] Touch = new bool[this.BaseGridProblemMapping.BasisS.Count()];
            foreach (var c in m_Config) {
                foreach (var iVar in c.VarIndex) {
                    if (Touch[iVar] == true) {
                        throw new ArgumentException("Variable #" + iVar + " is specified (at least) twice.");
                    }
                    Touch[iVar] = true;
                }

                if(c.DegreeS.Length != c.VarIndex.Length) {
                    throw new ArgumentException("Length of Degrees must match number of Variables.");
                }
            }

            for (int iVar = 0; iVar < Touch.Length; iVar++) {
                if (Touch[iVar] == false) {
                    throw new ArgumentException("No configuration specified for variable #" + iVar + ".");
                }
            }

           
        }


        /// <summary>
        /// applies the pre-conditioning to the operator matrix
        /// (passed in the constructor)
        /// and returns the pre-conditioned matrix
        /// </summary>
        /// <param name="LeftPreCond">
        /// left pre-conditioning matrix
        /// </param>
        /// <param name="RightPreCond">
        /// right pre-conditioning matrix
        /// </param>
        /// <param name="RightPreCondInv">
        /// the inverse of <paramref name="RightPreCond"/> -- usually required to transform an initial guess.
        /// </param>
        /// <param name="LeftPreCondInv"></param>
        /// <param name="MassMatrix">
        /// on entry the mass matrix w.r.t. the XDG basis
        /// </param>
        /// <param name="OpMatrix">
        /// </param>
        /// <returns>
        /// List of indefinite row indices.
        /// </returns>
        long[] ComputeChangeOfBasis(BlockMsrMatrix OpMatrix, BlockMsrMatrix MassMatrix, out BlockMsrMatrix LeftPreCond, out BlockMsrMatrix RightPreCond, out BlockMsrMatrix LeftPreCondInv, out BlockMsrMatrix RightPreCondInv) {
            using (var tr = new FuncTrace()) {
                // test arguments
                // ==============
                VerifyConfig();
                Debug.Assert(OpMatrix.RowPartitioning.LocalLength == this.Mapping.LocalLength);
                Debug.Assert(OpMatrix.ColPartition.LocalLength == this.Mapping.LocalLength);
                Debug.Assert(MassMatrix == null || (MassMatrix.RowPartitioning.LocalLength == this.Mapping.LocalLength));
                Debug.Assert(MassMatrix == null || (MassMatrix.ColPartition.LocalLength == this.Mapping.LocalLength));

                AggregationGridBasis[] basisS = this.Mapping.AggBasis;
                int[] Degrees = this.Mapping.DgDegree;
                bool marker = false;

                List<long> IndefRows = new List<long>();


                // compute preconditioner matrices
                // ===============================
                using (var bt = new BlockTrace("compute-pc", tr)) {
                    Stopwatch stw_Data = new Stopwatch(); stw_Data.Reset();
                    Stopwatch stw_Comp = new Stopwatch(); stw_Comp.Reset();


                    LeftPreCond = new BlockMsrMatrix(OpMatrix._RowPartitioning, OpMatrix._ColPartitioning);
                    RightPreCond = new BlockMsrMatrix(OpMatrix._RowPartitioning, OpMatrix._ColPartitioning);
                    RightPreCondInv = new BlockMsrMatrix(OpMatrix._RowPartitioning, OpMatrix._ColPartitioning);
                    LeftPreCondInv = new BlockMsrMatrix(OpMatrix._RowPartitioning, OpMatrix._ColPartitioning);
                    LeftPreCond.AccEyeSp(1.0);
                    RightPreCond.AccEyeSp(1.0);
                    LeftPreCondInv.AccEyeSp(1.0);
                    RightPreCondInv.AccEyeSp(1.0);


                    int LL = this.m_Config.Length;
                    MultidimensionalArray[] MassBlock = new MultidimensionalArray[LL];
                    MultidimensionalArray[] OperatorBlock = new MultidimensionalArray[LL];
                    MultidimensionalArray[] PCleftBlock = new MultidimensionalArray[LL];
                    MultidimensionalArray[] work = new MultidimensionalArray[LL];
                    MultidimensionalArray[] PCrightBlock_inv = new MultidimensionalArray[LL];
                    MultidimensionalArray[] PCleftBlock_inv = new MultidimensionalArray[LL];
                    MultidimensionalArray[] PCrightBlock = new MultidimensionalArray[LL];
                    long[][] __i0s = new long[LL][];
                    int[][] __Lns = new int[LL][];

                    for (int i = 0; i < LL; i++) {
                        var conf = m_Config[i];
                        __i0s[i] = new long[conf.VarIndex.Length];
                        __Lns[i] = new int[conf.VarIndex.Length];
                    }

                    int J = this.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    long i0 = this.Mapping.Partitioning.i0; //Processor offset

                    //List<Chunk> blaagr = new List<Chunk>();

                    //for (int jCell = 0; jCell < J; jCell++) {
                    //    for (int i = 0; i < LL; i++) {
                    //        long[] _i0s = __i0s[i];
                    //        int[] _Lns = __Lns[i];
                    //        var conf = m_Config[i];
                    //        int E = conf.VarIndex.Length;
                    //        for (int e = 0; e < E; e++) {
                    //            int iVar = conf.VarIndex[e];
                    //            _i0s[e] = this.Mapping.LocalUniqueIndex(iVar, jCell, 0) + i0;
                    //            _Lns[e] = basisS[iVar].GetLength(jCell, Degrees[iVar]);
                    //        }
                    //        ExtractBlock(_i0s, _Lns, true, MassMatrix, ref MassBlock[i]);
                    //        ExtractBlock(_i0s, _Lns, true, OpMatrix, ref OperatorBlock[i]);

                    //        if (MassBlock[i].InfNorm() == 0) {
                    //            Console.WriteLine("zero massmtx block in " + jCell);
                    //            double[] coords = this.Mapping.AggGrid.iLogicalCells.GetCenter(jCell);
                    //            foreach (var c in coords)
                    //                Console.WriteLine(c);
                    //            blaagr.Add(Chunk.GetSingleElementChunk(jCell));
                    //        }

                    //    }
                    //}

                    //(new CellMask(GridData, blaagr)).SaveToTextFile("failCell.csv");

                    for (int jCell = 0; jCell < J; jCell++) { // loop over cells...
                                                              //ReducedRegionCode rrc;
                                                              //int NoOfSpc = LsTrk.GetNoOfSpecies(jCell, out rrc);

                        //if (this.Mapping.GetLength(jCell) == 0)
                        //    // void cell
                        //    continue;

                        for (int i = 0; i < LL; i++) { // for each configuration item...
                            var conf = m_Config[i];

                            int E = conf.VarIndex.Length; //number of variables (e.g. 3 for u_x, u_y, p)
                            long[] _i0s = __i0s[i];
                            int[] _Lns = __Lns[i];
                            //AggregationGridBasis basis = null;

                            int DOF = 0;
                            bool AnyZeroLength = false;
                            for (int e1 = 0; e1 < E; e1++) {
                                int dof_var = this.Mapping.GetLengthForVar(jCell, conf.VarIndex[e1]);
                                DOF += dof_var;
                                AnyZeroLength |= (dof_var == 0);
                            }
                            if (AnyZeroLength && DOF > 0)
                                throw new ApplicationException();

                            if (DOF == 0)
                                // void cell
                                continue;

                            for (int e = 0; e < E; e++) {
                                int iVar = conf.VarIndex[e];
                                _i0s[e] = this.Mapping.LocalUniqueIndex(iVar, jCell, 0) + i0;
                                _Lns[e] = basisS[iVar].GetLength(jCell, Degrees[iVar]);
                            }

                            // extract blocks from operator and mass matrix
                            // --------------------------------------------

                            stw_Data.Start();
                            ExtractBlock(_i0s, _Lns, true, MassMatrix, ref MassBlock[i]);
                            ExtractBlock(_i0s, _Lns, true, OpMatrix, ref OperatorBlock[i]);
                            stw_Data.Stop();
                            double MassBlkNrm = MassBlock[i].InfNorm();
                            double OperatorBlkNrm = OperatorBlock[i].InfNorm();
                            int NN = MassBlock[i].NoOfRows;

                            if (MassBlkNrm == 0) {
                                //throw new ArithmeticException("absolute zero Mass block in cell " + jCell + ".");
                                //Console.WriteLine("absolute zero Mass block in cell " + jCell + ".");

                                //if (conf.mode == Mode.IdMass_DropIndefinite || conf.mode == Mode.SymPart_DiagBlockEquilib_DropIndefinite) {
                                //    // we can deal with this ...


                                //Console.WriteLine("Error at level" + i);
                                //throw new ArithmeticException("absolute zero Mass block in cell " + jCell + ".");


                                int Length = 0;
                                foreach (int len in _Lns) {
                                    Length += len;
                                }
                                long[] ZeroIdc = Length.ForLoop(Z => _i0s[0] + Z);
                                IndefRows.AddRange(ZeroIdc);
                                if (!marker) {
                                    Console.WriteLine("Zero mass matrix blocks detected. This should not happen, but we can deal with this.");
                                    Console.WriteLine("zero blocks at cell: ");
                                    marker = true;
                                }
                                Console.WriteLine(jCell);
                            } else {

                                //if(OperatorBlkNrm == 0) {
                                //    throw new ArithmeticException("absolute zero Operator block in cell " + jCell + ".");
                                //}

                                // mem alloc
                                // ---------

                                if (PCleftBlock[i] == null || PCleftBlock[i].NoOfRows != NN) {
                                    PCleftBlock[i] = MultidimensionalArray.Create(NN, NN);
                                }
                                if (PCrightBlock[i] == null || PCrightBlock[i].NoOfRows != NN) {
                                    PCrightBlock[i] = MultidimensionalArray.Create(NN, NN);
                                }
                                if (work[i] == null || work[i].NoOfRows != NN) {
                                    work[i] = MultidimensionalArray.Create(NN, NN);
                                }

                                // compute precond
                                // ---------------


                                stw_Comp.Start();
                                int Rank;
                                PCleftBlock[i].Clear();
                                PCrightBlock[i].Clear();
                                int[] idr = ComputeChangeOfBasisBlock(_Lns, MassBlock[i], OperatorBlock[i], PCleftBlock[i], PCrightBlock[i], conf.mode, out Rank, work[i]);
                                if (Rank != NN) {
                                    IndefRows.AddRange(ConvertRowIndices(jCell, basisS, Degrees, conf, E, _i0s, idr));
                                } else {
                                    Debug.Assert(idr == null);
                                }
                                stw_Comp.Stop();

                                // write block back
                                // ----------------
                                stw_Data.Start();
                                ExtractBlock(_i0s, _Lns, false, LeftPreCond, ref PCleftBlock[i]);
                                ExtractBlock(_i0s, _Lns, false, RightPreCond, ref PCrightBlock[i]);


                                // inverse precond-matrix
                                // ----------------------
                                // right-inverse: (required for transforming solution guess)
                                if (PCrightBlock_inv[i] == null || PCrightBlock_inv[i].NoOfRows != NN) {
                                    PCrightBlock_inv[i] = MultidimensionalArray.Create(NN, NN);
                                }
                                if (Rank == NN)
                                    PCrightBlock[i].InvertTo(PCrightBlock_inv[i]);
                                else
                                     RankDefInvert(PCrightBlock[i], PCrightBlock_inv[i]);

                                ExtractBlock(_i0s, _Lns, false, RightPreCondInv, ref PCrightBlock_inv[i]);

                                // left-inverse: (required for analysis purposes, to transform residuals back onto original grid)
                                if (PCleftBlock_inv[i] == null || PCleftBlock_inv[i].NoOfRows != NN) {
                                    PCleftBlock_inv[i] = MultidimensionalArray.Create(NN, NN);
                                }
                                if (Rank == NN)
                                    PCleftBlock[i].InvertTo(PCleftBlock_inv[i]);
                                else
                                    RankDefInvert(PCleftBlock[i], PCleftBlock_inv[i]);

                                ExtractBlock(_i0s, _Lns, false, LeftPreCondInv, ref PCleftBlock_inv[i]);

                                stw_Data.Stop();
                            }
                        }
                    }

                    bt.LogDummyblock(stw_Data.Elapsed.Ticks, "Change_of_Basis_data_copy");
                    bt.LogDummyblock(stw_Comp.Elapsed.Ticks, "Change_of_Basis_compute");
                }


                return IndefRows.ToArray();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_i0s">Offset for each field/variable</param>
        /// <param name="_Lns">Lengths for degrees of freedom for each field/variable</param>
        /// <param name="Sp2Full"></param>
        /// <param name="MtxSp"></param>
        /// <param name="MtxFl"></param>
        private static void ExtractBlock(
            long[] _i0s, 
            int[] _Lns,
            bool Sp2Full,
            BlockMsrMatrix MtxSp, ref MultidimensionalArray MtxFl) //
        {
            Debug.Assert(_i0s.Length == _Lns.Length);
            int E = _i0s.Length;

            int NN = _Lns.Sum();
            if (MtxFl == null || MtxFl.NoOfRows != NN) {
                Debug.Assert(Sp2Full == true);
                MtxFl = MultidimensionalArray.Create(NN, NN);
            } else {
                if (Sp2Full) {
                    MtxFl.Clear();
                }
            }

            if(!Sp2Full) {
                Debug.Assert(MtxSp != null);
            }


            int i0Rowloc = 0;
            for (int eRow = 0; eRow < E; eRow++) { // loop over variables in configuration
                long i0Row = _i0s[eRow];
                int NRow = _Lns[eRow];

                int i0Colloc = 0;
                for (int eCol = 0; eCol < E; eCol++) { // loop over variables in configuration
                    long i0Col = _i0s[eCol];
                    int NCol = _Lns[eCol]; 

                    MultidimensionalArray MtxFl_blk;
                    if (i0Rowloc == 0 && NRow == MtxFl.GetLength(0) && i0Colloc == 0 && NCol == MtxFl.GetLength(1)) {
                        MtxFl_blk = MtxFl;
                    } else {
                        MtxFl_blk = MtxFl.ExtractSubArrayShallow(new[] { i0Rowloc, i0Colloc }, new[] { i0Rowloc + NRow - 1, i0Colloc + NCol - 1 });
                    }

                    if (Sp2Full) {
                        if (MtxSp != null) {
                            MtxSp.ReadBlock(i0Row, i0Col, MtxFl_blk);
                        } else {
                            MtxFl_blk.AccEye(1.0);
                        }

                    } else {
                        MtxSp.AccBlock(i0Row, i0Col, 1.0, MtxFl_blk, 0.0);
                    }
#if DEBUG
                    for(int n_row = 0; n_row < NRow; n_row++) { // row loop...
                        for(int n_col = 0; n_col < NCol; n_col++) { // column loop...
                            Debug.Assert(MtxFl[n_row + i0Rowloc, n_col + i0Colloc] == ((MtxSp != null) ? ( MtxSp[n_row + i0Row, n_col + i0Col]) : (n_col == n_row ? 1.0 : 0.0)));
                        }
                    }
#endif
                    i0Colloc += NCol;
                }
                i0Rowloc += NRow;
            }
        }

        /// <summary>
        /// Convert local indices to global indices
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="basisS"></param>
        /// <param name="Degrees">Degree for each variable</param>
        /// <param name="conf"></param>
        /// <param name="E">Number of variables</param>
        /// <param name="_i0s"></param>
        /// <param name="LocIdx"></param>
        /// <returns></returns>
        private static long[] ConvertRowIndices(
            int jCell,
            AggregationGridBasis[] basisS, int[] Degrees,
            ChangeOfBasisConfig conf,
            int E, long[] _i0s,
            int[] LocIdx) {

            int NN = conf.VarIndex.Sum(iVar => basisS[iVar].GetLength(jCell, Degrees[iVar]));
            long[] Loc2glob = new long[NN];

            int i0Rowloc = 0;
            for (int eRow = 0; eRow < E; eRow++) { // loop over variables in configuration
                long i0Row = _i0s[eRow];
                int iVarRow = conf.VarIndex[eRow];

                int NRow = basisS[iVarRow].GetLength(jCell, Degrees[iVarRow]);

                for (int n_row = 0; n_row < NRow; n_row++) { // row loop...
                    //n_row = LocIdx[k];

                    int iRowLoc = n_row + i0Rowloc;
                    long iRowGlb = n_row + i0Row;

                    Loc2glob[iRowLoc] = iRowGlb;


                }
                i0Rowloc += NRow;
            }

            return LocIdx.Select(i => Loc2glob[i]).ToArray();
        }


        /// <summary>
        /// Inverse cholesky factorization (L * L^T) giving the result R * L = M^(-1) where R= L^T
        /// </summary>
        /// <param name="M">Input matrix (will be reduced to identity during the call)</param>
        /// <param name="L"> Inverse of L decomposition</param>
        /// <param name="R">Inverse of R decomposition, i.e. L^T </param>
        /// <exception cref="ArithmeticException"></exception>
        static void SymmInv(MultidimensionalArray M, MultidimensionalArray L, MultidimensionalArray R) {
            L.Clear();
            L.StructureType = MatrixStructure.LowerTriangular;
            R.Clear();
            R.StructureType = MatrixStructure.UpperTriangular;
            L.AccEye(1.0);
#if DEBUG
            var Mbefore = M.CloneAs();
#endif
            int n = M.NoOfRows;
            unsafe
            {

                void RowScale(double* pS, int i, double alpha, int RowCyc)
                {
                    pS += i * RowCyc;
                    for (int nn = 0; nn < n; nn++)
                    {
                        *pS *= alpha;
                        pS++;
                    }
                }

                void ColScale(double* pS, int i, double alpha, int RowCyc)
                {
                    pS += i;
                    for (int nn = 0; nn < n; nn++)
                    {
                        *pS *= alpha;
                        pS += RowCyc;
                    }
                }

                void RowAdd(double* pS, int iSrc, int iDst, double alpha, int RowCyc) 
                {
                    double* pDest = pS + iDst* RowCyc;
                    double* pSrc = pS + iSrc* RowCyc;
                    for (int l = 0; l < n; l++)
                    {
                        *pDest += *pSrc * alpha;
                        pDest++;
                        pSrc++;
                    }
                }

                void ColAdd(double* pS, int iSrc, int iDst, double alpha, int RowCyc)
                {
                    double* pDest = pS + iDst;
                    double* pSrc = pS + iSrc;
                    for (int l = 0; l < n; l++)
                    {
                        *pDest += *pSrc * alpha;
                        pDest += RowCyc;
                        pSrc += RowCyc;
                    }
                }

                fixed (double* _pL = L.Storage, _pM = M.Storage)
                {
                    int RowCycL = n > 1 ? (L.Index(1, 0) - L.Index(0, 0)) : 0;
                    int RowCycM = n > 1 ? (M.Index(1, 0) - M.Index(0, 0)) : 0;

                    double* pL = _pL + L.Index(0, 0);
                    double* pM = _pM + M.Index(0, 0);


                    for (int i = 0; i < n; i++)
                    {
                        double M_ii = M[i, i];
                        if (M_ii == 0.0)
                            throw new ArithmeticException("Zero diagonal element at " + i + "-th row.");
                        double scl = 1.0 / Math.Sqrt(Math.Abs(M_ii));
                        //M.RowScale(i, scl);
                        //L.RowScale(i, scl);
                        //M.ColScale(i, scl);
                        RowScale(pM, i, scl, RowCycM);
                        RowScale(pL, i, scl, RowCycL);
                        ColScale(pM, i, scl, RowCycM);

                        double diagsign = Math.Sign(M[i, i]);
                        if (diagsign == 0.0)
                            throw new ArithmeticException("Zero diagonal element at " + i + "-th row.");
                        if (Math.Abs(Math.Abs(M[i, i]) - 1.0) > 1.0e-8)
                            throw new ArithmeticException("Unable to create diagonal 1.0.");

                        for (int k = i + 1; k < n; k++)
                        {
                            double M_ki = M[k, i];

                            RowAdd(pM, i, k, -M_ki * diagsign, RowCycM);
                            RowAdd(pL, i, k, -M_ki * diagsign, RowCycL);
                            ColAdd(pM, i, k, -M_ki * diagsign, RowCycM);

                            Debug.Assert(Math.Abs(M[k, i]) < 1.0e-8);
                            Debug.Assert(Math.Abs(M[i, k]) < 1.0e-8);
                        }

                        /*
                        unsafe
                        {
                            fixed (double* B_entries = Lo.Storage)
                            {

                                int UPLO = 'L', DIAG = 'N';
                                LAPACK.F77_LAPACK.DSYTRF_(ref UPLO, ref DIAG, ref n, B_entries, ref n, out info);
                            }
                        }
                        */
                    }
                }
            }
            L.TransposeTo(R);

#if DEBUG
            var Test = MultidimensionalArray.Create(M.Lengths);
            //var Q = MultidimensionalArray.Create(M.Lengths);
            //Test.AccEye(1.0);
            Test.Multiply(-1.0, L, Mbefore, R, 1.0, "ij", "ik", "kl", "lj");
            for(int i = 0; i < n; i++) {
                //Debug.Assert((Test[i, i].Abs() - 1.0).Abs() < 1.0e-8);
                //Test[i, i] -= Math.Sign(Test[i, i]);
                Test[i, i] = 0;
            }

            double TestNorm = Test.InfNorm();
            double scale = Math.Max(Mbefore.InfNorm(), R.InfNorm());
            Debug.Assert(TestNorm / scale < 1.0e-4);

            /*
            if(TestNorm / scale >= 1.0e-8) {
                var MM = Mbefore.CloneAs();
                MultidimensionalArray Lo = MultidimensionalArray.Create(n, n);

                for(int j = 0; j < n; j++) {
                    double Lo_jj = MM[j, j];
                    for(int k = 0; k < j; k++) {
                        Lo_jj -= Lo[j, k].Pow2();
                    }

                    double sig = Math.Abs(Lo_jj);
                    Lo[j, j] = Math.Sqrt(Lo_jj * sig);


                    for(int i = j; i < n; i++) {
                        double acc = MM[i, j];
                        for(int k = 0; k < j; k++) {
                            acc -= Lo[i, k] * Lo[j, k];
                        }

                        Lo[i, j] = (1 / (Lo[j, j] * sig)) * acc;
                    }
                }

                int info = 0;
                unsafe {
                    fixed(double* B_entries = Lo.Storage) {

                        int UPLO = 'L', DIAG = 'N';
                        LAPACK.F77_LAPACK.DTRTRI_(ref UPLO, ref DIAG, ref n, B_entries, ref n, out info);
                    }
                }

                MultidimensionalArray Up = MultidimensionalArray.Create(n, n);
                Lo.TransposeTo(Up);
                Test.Clear();
                Test.Multiply(-1.0, Lo, Mbefore, R, 1.0, "ij", "ik", "kl", "lj");


                for(int i = 0; i < n; i++) {
                    //Debug.Assert((Test[i, i].Abs() - 1.0).Abs() < 1.0e-8);
                    Test[i, i] -= Math.Sign(Test[i, i]);
                }

                double TestNorm1 = Test.InfNorm();
                //double scale = Math.Max(Mbefore.InfNorm(), R.InfNorm());
                Console.WriteLine(TestNorm1);

            }
             */

#endif
        }

        private static int[] ComputeChangeOfBasisBlock(
            int[] BlockLen,
            MultidimensionalArray In_MassMatrixBlock, MultidimensionalArray In_OperatorMatrixBlock, 
            MultidimensionalArray OUT_LeftPC, MultidimensionalArray OUT_rightPC, Mode PCMode, out int Rank,
            MultidimensionalArray work) //
        {
            var MMclone = In_MassMatrixBlock.CloneAs();
            var OPclone = In_OperatorMatrixBlock.CloneAs();
            try {


                Rank = In_MassMatrixBlock.NoOfCols;

                int[] IndefRows = null;

                switch (PCMode) {
                    case Mode.Eye: {
                            OUT_LeftPC.AccEye(1.0);
                            OUT_rightPC.AccEye(1.0);
                            break;
                        }

                    case Mode.DiagBlockEquilib: {
                            double symmErr = In_OperatorMatrixBlock.SymmetryError();
                            double infNorm = In_OperatorMatrixBlock.InfNorm();

                            if (symmErr / infNorm > 1.0e-8)
                                throw new ArithmeticException(string.Format("LDL_DiagBlock is not supported on unsymmetric matrices (Symm-Err: {0:0.####E-00}, Inf-Norm: {1:0.####E-00}, Quotient {2:0.####E-00}).", symmErr, infNorm, symmErr / infNorm));
                            SymmInv(In_OperatorMatrixBlock, OUT_LeftPC, OUT_rightPC);

                            break;
                        }

                    case Mode.SymPart_DiagBlockEquilib: {
                            var SymmPart = work;
                            In_OperatorMatrixBlock.TransposeTo(SymmPart);
                            SymmPart.Acc(1.0, In_OperatorMatrixBlock);
                            SymmPart.Scale(0.5);

                            SymmInv(SymmPart, OUT_LeftPC, OUT_rightPC);
                            break;
                        }

                    case Mode.IdMass: {
                            In_MassMatrixBlock.SymmetricLDLInversion(OUT_rightPC, default(double[]));
                            OUT_rightPC.TransposeTo(OUT_LeftPC);
                            break;
                        }

                    case Mode.LeftInverse_DiagBlock: {
                            In_OperatorMatrixBlock.InvertTo(OUT_LeftPC);
                            OUT_rightPC.AccEye(1.0);
                            break;
                        }

                    case Mode.LeftInverse_Mass: {
                            In_MassMatrixBlock.InvertTo(OUT_LeftPC);
                            OUT_rightPC.AccEye(1.0);
                            break;
                        }

                    case Mode.IdMass_DropIndefinite: {
                            int[] ZerosEntries = ModifiedInverseChol(In_MassMatrixBlock, OUT_rightPC, 1.0e-12, false);
                            int NoOfZeros = ZerosEntries == null ? 0 : ZerosEntries.Length;
                            IndefRows = ZerosEntries;
                            Rank = OUT_LeftPC.NoOfCols - NoOfZeros;
                            OUT_rightPC.TransposeTo(OUT_LeftPC);
                            break;
                        }

                    case Mode.SymPart_DiagBlockEquilib_DropIndefinite: {

                            (Rank, IndefRows) = SymPart_DiagBlockEquilib_DropIndefinite(In_MassMatrixBlock, In_OperatorMatrixBlock, OUT_LeftPC, OUT_rightPC, work);

                            break;
                        }

                    case Mode.SchurComplement: {

                            //throw new NotImplementedException("todo");

                            int NoVars = BlockLen.Length;
                            if (NoVars <= 1)
                                throw new NotSupportedException("The Schur complement requires at least 2 variables, moron!");
                            int N1 = BlockLen.Take(NoVars - 1).Sum();
                            int N2 = BlockLen.Last();
                            int N = N1 + N2;

                            //string TestPath = @"C:\Users\flori\OneDrive\MATLAB\Schur";
                            //In_OperatorMatrixBlock.SaveToTextFile(Path.Combine(TestPath, "Opm.txt"));

                            Debug.Assert(N == In_OperatorMatrixBlock.NoOfRows);
                            Debug.Assert(N == In_OperatorMatrixBlock.NoOfCols);

                            (MultidimensionalArray M11, MultidimensionalArray M12, MultidimensionalArray M21, MultidimensionalArray M22) GetSubblox(MultidimensionalArray Mtx) {
                                return (
                                    Mtx.ExtractSubArrayShallow(new[] { 0, 0 }, new[] { N1 - 1, N1 - 1 }),
                                    Mtx.ExtractSubArrayShallow(new[] { 0, N1 }, new[] { N1 - 1, N - 1 }),
                                    Mtx.ExtractSubArrayShallow(new[] { N1, 0 }, new[] { N - 1, N1 - 1 }),
                                    Mtx.ExtractSubArrayShallow(new[] { N1, N1 }, new[] { N - 1, N - 1 })
                                );
                            }

                            var OpMtxSub = GetSubblox(In_OperatorMatrixBlock);
                            var MamaSub = GetSubblox(In_MassMatrixBlock);
                            var workSub = GetSubblox(work);
                            var lpcSub = GetSubblox(OUT_LeftPC);
                            var rpcSub = GetSubblox(OUT_rightPC);

                            (int Rank11, int[] IndefRows11) = SymPart_DiagBlockEquilib_DropIndefinite(MamaSub.M11, OpMtxSub.M11, lpcSub.M11, rpcSub.M11, workSub.M11);

                            // compute lpcSub.M11*OpMtxSub.M11*rpcSub.M11
                            // (without additional mem-alloc, yeah!)
                            var Diag11 = workSub.M11;
                            Diag11.GEMM(1.0, lpcSub.M11, OpMtxSub.M11, 0.0);
                            OpMtxSub.M11.GEMM(1.0, Diag11, rpcSub.M11, 0.0);

                            // invert diag
                            if (Rank11 == N1) {
                                OpMtxSub.M11.InvertTo(workSub.M11);
                            } else {
                                RankDefInvert(OpMtxSub.M11, workSub.M11);
                            }


                            ////Copying the matrix for multiplication (dgemm)
                            //var OpMtxSub11 = default(MultidimensionalArray);
                            //OpMtxSub11 = OpMtxSub.M11.CloneAs();
                            OpMtxSub.M11.GEMM(1.0, rpcSub.M11, OpMtxSub.M11, 0.0); // fk, 26aug22, Note: now, GEMM should be able to handle in-place correctly.

                            workSub.M11.GEMM(1.0, OpMtxSub.M11, lpcSub.M11, 0.0);
                            var Q = workSub.M11;

                            //
                            workSub.M21.GEMM(-1.0, OpMtxSub.M21, Q, 0.0);
                            workSub.M12.GEMM(-1.0, Q, OpMtxSub.M12, 0.0);

                            //rpcSub.M12.SetMatrix(workSub.M12); rpcSub.M22.AccEye(1.0);
                            //lpcSub.M21.SetMatrix(workSub.M21); lpcSub.M22.AccEye(1.0);
                            //OUT_LeftPC.SaveToTextFile(Path.Combine(TestPath, "Lpc.txt"));
                            //OUT_rightPC.SaveToTextFile(Path.Combine(TestPath, "Rpc.txt"));


                            rpcSub.M12.GEMM(1.0, Q, OpMtxSub.M12, 0.0);
                            OpMtxSub.M22.GEMM(-1.0, OpMtxSub.M21, rpcSub.M12, 1.0);


                            (int Rank22, int[] IndefRows22) = SymPart_DiagBlockEquilib_DropIndefinite(MamaSub.M22, OpMtxSub.M22, lpcSub.M22, rpcSub.M22, workSub.M22);

                            lpcSub.M21.GEMM(1.0, lpcSub.M22, workSub.M21, 0.0);
                            rpcSub.M12.GEMM(1.0, workSub.M12, rpcSub.M22, 0.0);

                            //OUT_LeftPC.SaveToTextFile(Path.Combine(TestPath, "Lpc.txt"));
                            //OUT_rightPC.SaveToTextFile(Path.Combine(TestPath, "Rpc.txt"));

                            if (IndefRows11 != null && IndefRows22 != null)
                                IndefRows = ArrayTools.Cat(IndefRows11, IndefRows22);
                            else if (IndefRows11 != null)
                                IndefRows = IndefRows11;
                            else if (IndefRows22 != null)
                                IndefRows = IndefRows22;
                            else
                                IndefRows = null;

                            Rank = Rank11 + Rank22;
                            break;
                        }

                    /*
                    case Mode.LDL_DiagBlock_DropIndefinite: {
                            if(In_OperatorMatrixBlock.SymmetryError() / In_OperatorMatrixBlock.InfNorm() > 1.0e-8)
                                throw new NotSupportedException("LDL_DiagBlock is not supported on unsymmetric matrices");
                            int N = OUT_LeftPC.NoOfCols;

                            MultidimensionalArray PL = MultidimensionalArray.Create(N, N);
                            MultidimensionalArray PR = MultidimensionalArray.Create(N, N);

                            int zeros1 = CRM114(In_MassMatrixBlock, PR, 1.0e-12);
                            Rank = N - zeros1;
                            PR.TransposeTo(PL);

                            var OpTr = (PL * In_OperatorMatrixBlock) * PR;

                            MultidimensionalArray QL = MultidimensionalArray.Create(N, N);
                            MultidimensionalArray QR =  MultidimensionalArray.Create(N, N);
                            int zeros2 = CRM114(OpTr, QR, 1.0e-12);
                            QR.TransposeTo(QL);

                            if(zeros1 != zeros2)
                                // I want this to fire also in Release mode, therfore i don't use Debug.Assert(...)
                                throw new ApplicationException();


                            OUT_LeftPC.GEMM(1.0, QL, PL, 0.0);
                            OUT_rightPC.GEMM(1.0, PR, QR, 0.0);

                            //if (zeros1 > 0 || zeros2 > 0) {
                            //    Console.WriteLine("zeros introduced: " + zeros1 + ", " + zeros2);
                            //}

                            break;
                        }
                        */
                    default:
                        throw new NotImplementedException("Unknown option: " + PCMode);
                }

                return IndefRows;
            } catch (Exception e) {
                int MPIrank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);
                MMclone.SaveToTextFile("MassMatrixBlock_" + MPIrank + ".txt");
                OPclone.SaveToTextFile("OperatorMatrixBlock_" + MPIrank + ".txt");

                throw e;
            }
        }



        /// <summary>
        /// Partial inversion of a matrix with zero rows and columns; not considered to be performance-critical,
        /// since it is only for treating pathological cases.
        /// </summary>
        static void RankDefInvert(MultidimensionalArray MtxIn, MultidimensionalArray MtxOt) {
            if(MtxIn.Dimension != 2)
                throw new ArgumentOutOfRangeException("Expecting an 2D-array, i.e. a matrix.");
            if(MtxOt.Dimension != 2)
                throw new ArgumentOutOfRangeException("Expecting an 2D-array, i.e. a matrix.");
            if(MtxIn.GetLength(0) != MtxOt.GetLength(0))
                throw new ArgumentOutOfRangeException("input and output matrix must have same size");
            if(MtxIn.GetLength(1) != MtxOt.GetLength(1))
                throw new ArgumentOutOfRangeException("input and output matrix must have same size");
            if(MtxIn.GetLength(0) != MtxOt.GetLength(1))
                throw new ArgumentOutOfRangeException("matrix must be quadratic");

            int M = MtxIn.GetLength(0);

            List<int> NonzeroIdx = new List<int>();
            for(int i = 0; i < M; i++) {
                double[] row_i = MtxIn.GetRow(i);
                double[] col_i = MtxIn.GetColumn(i);

                // If the matrix comes from Cholesky decomposition, it would have L or U shape and only the row or columns, resp., would be zeros.
                switch (MtxIn.StructureType) {
                    case MatrixStructure.LowerTriangular:
                        if (row_i.L2NormPow2() == 0.0)
                            NonzeroIdx.Add(i);
                    break;
                    case MatrixStructure.UpperTriangular:
                        if (col_i.L2NormPow2() == 0.0)
                            NonzeroIdx.Add(i);
                    break;
                    default: //check if both row and columns are zero since it is not a triangular matrix
                        if (row_i.L2NormPow2() == 0.0)
                            if (col_i.L2NormPow2() == 0.0) {
                    NonzeroIdx.Add(i);
                            } else {
                                throw new ArgumentException("Row is zero, but column is not.");
                            }
                        break;
                }
            }

            if (NonzeroIdx.Count == M) {
                MtxIn.SaveToTextFileUnsteady("d_problemMatrix");
                throw new ArgumentException("Unable to find zero row/column.");
            }

            int Q = NonzeroIdx.Count;
            MultidimensionalArray TmpIn = MultidimensionalArray.Create(Q, Q);
            MultidimensionalArray TmpOt = MultidimensionalArray.Create(Q, Q);

            for(int i = 0; i < Q; i++)
                for(int j = 0; j < Q; j++)
                    TmpIn[i, j] = MtxIn[NonzeroIdx[i], NonzeroIdx[j]];

            TmpIn.InvertTo(TmpOt);

            MtxOt.Clear();
            for(int i = 0; i < Q; i++)
                for(int j = 0; j < Q; j++)
                    MtxOt[NonzeroIdx[i], NonzeroIdx[j]] = TmpOt[i, j];
        }

        /// <summary>
        /// Implementation of <see cref="Mode.SymPart_DiagBlockEquilib_DropIndefinite"/>
        /// </summary>
        static (int Rank, int[] IndefRows) SymPart_DiagBlockEquilib_DropIndefinite(
            MultidimensionalArray In_MassMatrixBlock, MultidimensionalArray In_OperatorMatrixBlock, 
            MultidimensionalArray OUT_LeftPC, MultidimensionalArray OUT_rightPC,
            MultidimensionalArray work ) {
            int MPIrank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);

            var Owork = work.CloneAs();

            var SymmPart = work;
            In_OperatorMatrixBlock.TransposeTo(SymmPart);
            SymmPart.Acc(1.0, In_OperatorMatrixBlock);
            SymmPart.Scale(0.5);
            var ASymmPart = SymmPart.CloneAs();
            var AOUT_rightPC = OUT_rightPC.CloneAs();

            int[] ZerosEntries = ModifiedInverseChol(In_MassMatrixBlock, OUT_rightPC, 1.0e-12, false);
            var BOUT_rightPC = OUT_rightPC.CloneAs();

            int NoOfZeros = ZerosEntries == null ? 0 : ZerosEntries.Length;
            int[] _IndefRows = ZerosEntries;
            int _Rank = OUT_LeftPC.NoOfCols - NoOfZeros;

            if (NoOfZeros == 0) {
                // normal cell -- nix indefinite
                // +++++++++++++++++++++++++++++

                SymmInv(SymmPart, OUT_LeftPC, OUT_rightPC);
            } else {
                // problem-cell
                // ++++++++++++++

                var SymmPartTest = SymmPart.CloneAs();
                var OUT_rightPCTest = OUT_rightPC.CloneAs();
                var OUT_leftPCTest = OUT_LeftPC.CloneAs();
                SymmInv(SymmPartTest, OUT_leftPCTest, OUT_rightPCTest);

                SymmPartTest.SaveToTextFileUnsteady($"d_SymmPartTest.txt");
                OUT_rightPCTest.SaveToTextFileUnsteady($"d_OUT_rightPCTest.txt");
                OUT_leftPCTest.SaveToTextFileUnsteady($"d_OUT_leftPCTest.txt");



                OUT_rightPC.TransposeTo(OUT_LeftPC);
                //SymmInv(SymmPart, OUT_LeftPC, OUT_rightPC);

                //SymmPart = IMatrixExtensions.GEMM(OUT_LeftPC, SymmPart, OUT_rightPC); //SymmPart = IMatrixExtensions.GEMM(OUT_LeftPC, SymmPart, OUT_rightPC);
                var dummyOUT_rightPC = OUT_rightPC.CloneAs();
                int[] ZerosEntries2 = ModifiedInverseChol(SymmPart, dummyOUT_rightPC, 1.0e-12, true);

                if (!ZerosEntries2.SetEquals(ZerosEntries)) {

                    foreach (int n in ZerosEntries) {
                        dummyOUT_rightPC.ColScale(n, 0.0);  //this the R part, hence we need to make columns zero

                        if (dummyOUT_rightPC.GetColumn(n).L2NormPow2() != 0.0)
                            throw new Exception("Column infinity is still non-zero");

                        SymmPart.ColScale(n, 0.0);
                        SymmPart.RowScale(n, 0.0);


                        if (SymmPart.GetRow(n).L2NormPow2() != 0.0)
                            throw new Exception("SymmPart Row infinity is still non-zero");


                        if (SymmPart.GetColumn(n).L2NormPow2() != 0.0)
                            throw new Exception("SymmPart column infinity is still non-zero");

                    }
                    dummyOUT_rightPC.TransposeTo(OUT_LeftPC);

                    In_MassMatrixBlock.SaveToTextFileUnsteady($"d_In_MassMatrixBlock.txt");
                    In_OperatorMatrixBlock.SaveToTextFileUnsteady($"d_In_OperatorMatrixBlock.txt");
                    
                    Owork.SaveToTextFileUnsteady($"d_work.txt");
                    ASymmPart.SaveToTextFileUnsteady($"d_SymmPart_A.txt");
                    SymmPart.SaveToTextFileUnsteady($"d_SymmPart_B.txt");
                    if (ZerosEntries == null)
                        ZerosEntries = new int[] { };

                    if (ZerosEntries2 == null)
                        ZerosEntries2 = new int[] { };


                    ZerosEntries.SaveToTextFileDebugUnsteady($"d_ZerosEntries", ".txt");
                    ZerosEntries2.SaveToTextFileDebugUnsteady($"d_ZerosEntries2", ".txt");

                    AOUT_rightPC.SaveToTextFileUnsteady($"d_OUT_rightPC_A.txt");
                    BOUT_rightPC.SaveToTextFileUnsteady($"d_OUT_rightPC_B.txt");
                    OUT_rightPC.SaveToTextFileUnsteady($"d_OUT_rightPC_C.txt");
                    dummyOUT_rightPC.SaveToTextFileUnsteady($"d_OUT_rightPC_D.txt");

                    //throw new ArithmeticException("Zero entries are not identical with the symmetric part");
                }
            }

            return (_Rank, _IndefRows);
        }

        /// <summary>
        /// Modified inverse Cholesky - resp.LDL - factorization, 
        /// $B =  \text{mchol}^{−1}(Q)$, resp. $B = \text{mldl}^{−1}(Q)$.
        /// In difference to the classical inverse Cholesky/LDL, the algorithm works on
        /// indefinite/singular matrices Q, but may produce zero - rows in the output <paramref name="B"/>.
        /// </summary>
        /// <param name="Q">Input: some symmetric matrix $Q \in \real^{N times N}$.</param>
        /// <param name="B">
        /// Output: an upper-diagonal matrix $B \in \real^{N \times N} , so that 
        /// $ B^T Q B = D$.
        /// If $Q$ is symmetrically positive definite, $D$ is the identity matrix and $B$ is invertible.
        /// Otherwise, $D$ is a diagonal matrix containing
        /// only entries 0 and 1 in the $  \text{mchol}^{−1}$-case resp. 0, −1 and + 1 in the $  \text{mldl}^{−1}$-case.
        /// </param>
        /// <param name="threshold"></param>
        /// <param name="ldl">
        /// Switch between modified Cholesky and LDL.
        /// </param>
        static int[] ModifiedInverseChol(MultidimensionalArray Q, MultidimensionalArray B, double threshold, bool ldl) {

            if (threshold < 0.0)
                throw new ArgumentOutOfRangeException();

            var Mtx = (Q.CloneAs());
            var LOW = B;
            LOW.Clear();
            LOW.StructureType = MatrixStructure.LowerTriangular;
            LOW.AccEye(1.0);
            int N = B.NoOfRows;

            List<int> zeros = null;
            

            for(int n = 0; n < N; n++) { // loop over columns

                /*
                // pivotisierung:
                double m_nn = Mtx[n, n];
                int imax = n;
                for(int i = n + 1; i < N; i++) {
                    double m_ii = Mtx[i, i];
                    if(m_ii > m_nn) {
                        m_nn = m_ii;
                        imax = i;
                    }
                }
                if(imax != n) {
                    Mtx.SwapRows(imax, n);
                    Mtx.SwapColumns(imax, n);
                    Debug.Assert(Mtx[n, n] == m_nn);
                    LOW.SwapRows(imax, n);
                }
                */

                double m_nn = Mtx[n, n];

                double discriminator;
                if (ldl)
                    discriminator = Math.Abs(m_nn);
                else
                    discriminator = m_nn;


                if (discriminator <= threshold) {
                    // indefinite entry:
                    //Mtx.RowMul(n, -1.0);
                    //Diag[n] *= -1.0;
                    //m_nn = Mtx[n, n];

                    LOW.RowScale(n, 0.0);
                    Mtx.ColScale(n, 0.0);
                    Mtx.RowScale(n, 0.0);

                    if (zeros == null)
                        zeros = new List<int>();
                    zeros.Add(n);
                    string str = $"n={n} discriminator={discriminator} - threshold={threshold}";
                    int MPIrank;
                    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);
                    str.SaveToTextFileDebugUnsteady("d_ZeroColumn",".txt");

                } else {
                    double sign_mnn = Math.Sign(m_nn);
                    if(sign_mnn == 0.0)
                        throw new ArithmeticException();
                    double oo_mnn = Math.Sqrt(Math.Abs(1.0 / m_nn));


                    Mtx.RowScale(n, oo_mnn);
                    LOW.RowScale(n, oo_mnn);

                    //double __M_nn = Mtx[n, n];
                    //double test = __M_nn*oo_mnn; // this should be 1.0! we don't need to do the column elimination, since it would have no effect
                    //                                on the output of this method, the LOW -- matrix.

                    for(int i = n + 1; i < N; i++) { // elimination loop: we eliminate row n+1 to N in column n
                        double oldVal = Mtx[i, n];
                        double fac = -oo_mnn * oldVal * sign_mnn;

                        Mtx.RowAdd(n, i, fac);
                        LOW.RowAdd(n, i, fac);

                        Debug.Assert(Math.Abs(Mtx[i, n]) <= Math.Max(Math.Abs(oldVal) * 1.0e-14, 1.0e-10));
                        Mtx[i, n] = 0.0;
                    }

                    Mtx[n, n] = sign_mnn;
                    for(int j = n + 1; j < N; j++) { // elimination loop: we eliminate col n+1 to N in row n 
                        Mtx[n, j] = 0.0;
                    }
                }
            }

            // we produced a lower triangular matrix, but we want upper triangular
            LOW.TransposeInPlace();

            return zeros == null ? null : zeros.ToArray();
        }

    }
}
