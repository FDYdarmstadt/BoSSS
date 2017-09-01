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

namespace BoSSS.Solution.Multigrid {
    partial class MultigridOperator {

        /// <summary>
        /// configuration item for the change-of-basis in 
        /// </summary>
        public class ChangeOfBasisConfig {

            /// <summary>
            /// a list of variable indices -- given by mapping <see cref="map"/> -- onto which this configuartion item applies.
            /// </summary>
            public int[] VarIndex;
            
            /// <summary>
            /// DG polynomial degree of respective variables on respective multigrid level
            /// </summary>
            public int Degree;

            ///// <summary>
            ///// if true, all species in one cell 
            ///// </summary>
            //public bool SpeciesSeparately;

            /// <summary>
            /// pre-conditioner mode
            /// </summary>
            public Mode mode;
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
            SymPart_DiagBlockEquilib_DropIndefinite
        }


        ChangeOfBasisConfig[] m_Config;

        /// <summary>
        /// the DG degrees on this level
        /// </summary>
        int[] Degrees {
            get {
                VerifyConfig();
                int[] R = new int[this.BaseGridProblemMapping.BasisS.Count()];
                foreach(var c in m_Config) {
                    foreach(var iVar in c.VarIndex) {
                        R[iVar] = c.Degree;
                    }
                }
                return R;
            }
        }


        void VerifyConfig() {
            bool[] Touch = new bool[this.BaseGridProblemMapping.BasisS.Count()];
            foreach(var c in m_Config) {
                foreach(var iVar in c.VarIndex) {
                    if(Touch[iVar] == true) {
                        throw new ArgumentException("Variable #" + iVar + " is specified (at least) twice.");
                    }
                    Touch[iVar] = true;
                }
            }

            for(int iVar = 0; iVar < Touch.Length; iVar++) {
                if(Touch[iVar] == false) {
                    throw new ArgumentException("No configuration specified for variable #" + iVar + ".");
                }
            }
        }


        /// <summary>
        /// applies the pre-conditioning to the operator matrix
        /// (passed in the constructor)
        /// and returns the pre-conditioned matrix
        /// </summary>
        /// <param name="PCndOpMatrix">
        /// on exit,
        /// <paramref name="LeftPreCond"/>*M*<paramref name="RightPreCond"/>,
        /// where M denotes the operator matrix passed in the constructor.
        /// </param>
        /// <param name="LeftPreCond">
        /// left pre-conditioning matrix
        /// </param>
        /// <param name="RightPreCond">
        /// right pre-conditioning matrix
        /// </param>
        /// <param name="RightPreCondInv">
        /// the inverse of <paramref name="RightPreCond"/> -- usually required to transform an initial guess.
        /// </param>
        /// <returns>
        /// List of indefinite row indices.
        /// </returns>
        int[] ComputeChangeOfBasis(BlockMsrMatrix OpMatrix, BlockMsrMatrix MassMatrix, out BlockMsrMatrix LeftPreCond, out BlockMsrMatrix RightPreCond, out BlockMsrMatrix LeftPreCondInv, out BlockMsrMatrix RightPreCondInv) {
            using(var tr = new FuncTrace()) {
                // test arguments
                // ==============
                VerifyConfig();
                Debug.Assert(OpMatrix.RowPartitioning.LocalLength == this.Mapping.LocalLength);
                Debug.Assert(OpMatrix.ColPartition.LocalLength == this.Mapping.LocalLength);
                Debug.Assert(MassMatrix == null || (MassMatrix.RowPartitioning.LocalLength == this.Mapping.LocalLength));
                Debug.Assert(MassMatrix == null || (MassMatrix.ColPartition.LocalLength == this.Mapping.LocalLength));

                AggregationGridBasis[] basisS = this.Mapping.AggBasis;
                int[] Degrees = this.Mapping.DgDegree;

                List<int> IndefRows = new List<int>();
                

                // compute preconditioner matrices
                // ===============================
                using (new BlockTrace("compute-pc", tr)) {
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
                    int[][] __i0s = new int[LL][];

                    for(int i = 0; i < LL; i++) {
                        var conf = m_Config[i];
                        __i0s[i] = new int[conf.VarIndex.Length];
                    }

                    int J = this.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    int i0 = this.Mapping.Partitioning.i0;
                    for(int jCell = 0; jCell < J; jCell++) { // loop over cells...
                        //ReducedRegionCode rrc;
                        //int NoOfSpc = LsTrk.GetNoOfSpecies(jCell, out rrc);

                        //if (this.Mapping.GetLength(jCell) == 0)
                        //    // void cell
                        //    continue;

                        for (int i = 0; i < LL; i++) { // for each configuration item...
                            var conf = m_Config[i];
                            
                            int E = conf.VarIndex.Length;
                            int[] _i0s = __i0s[i];
                            AggregationGridBasis basis = null;

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
                                _i0s[e] =  this.Mapping.LocalUniqueIndex(conf.VarIndex[e], jCell, 0) + i0;
                                if(e == 0) {
                                    basis = basisS[conf.VarIndex[e]];
                                } else {
                                    if(!object.ReferenceEquals(basis, basisS[conf.VarIndex[e]])) {
                                        throw new NotSupportedException("All variables in a configuration item must share the same basis.");
                                    }
                                }
                            }
                            
                            // extract blocks from operator and mass matrix
                            // --------------------------------------------
                                                       
                            ExtractBlock(jCell, basis, Degrees, conf, E, _i0s, true, MassMatrix, ref MassBlock[i]);
                            ExtractBlock(jCell, basis, Degrees, conf, E, _i0s, true, OpMatrix, ref OperatorBlock[i]);
                            double MassBlkNrm = MassBlock[i].InfNorm();
                            double OperatorBlkNrm = OperatorBlock[i].InfNorm();
                            int NN = MassBlock[i].NoOfRows;

                            if (MassBlkNrm == 0) {
                                //throw new ArithmeticException("absolute zero Mass block in cell " + jCell + ".");
                                //Console.WriteLine("absolute zero Mass block in cell " + jCell + ".");

                                if (conf.mode == Mode.IdMass_DropIndefinite || conf.mode == Mode.SymPart_DiagBlockEquilib_DropIndefinite) {
                                    // we can deal with this ...

                                } else {
                                    throw new ArithmeticException("absolute zero Mass block in cell " + jCell + ".");
                                }
                            }
                            //if(OperatorBlkNrm == 0) {
                            //    throw new ArithmeticException("absolute zero Operator block in cell " + jCell + ".");
                            //}

                            // mem alloc
                            // ---------

                            if(PCleftBlock[i] == null || PCleftBlock[i].NoOfRows != NN) {
                                PCleftBlock[i] = MultidimensionalArray.Create(NN, NN);
                            }
                            if(PCrightBlock[i] == null || PCrightBlock[i].NoOfRows != NN) {
                                PCrightBlock[i] = MultidimensionalArray.Create(NN, NN);
                            } if(work[i] == null || work[i].NoOfRows != NN) {
                                work[i] = MultidimensionalArray.Create(NN, NN);
                            }

                            // compute precond
                            // ---------------

                            

                            int Rank;
                            PCleftBlock[i].Clear();
                            PCrightBlock[i].Clear();
                            int[] idr = ComputeChangeOfBasisBlock(MassBlock[i], OperatorBlock[i], PCleftBlock[i], PCrightBlock[i], conf.mode, out Rank, work[i]);
                            if (Rank != NN) {
                                IndefRows.AddRange(ConvertRowIndices(jCell, basis, Degrees, conf, E, _i0s, idr));
                            } else {
                                Debug.Assert(idr == null); 
                            }
                            // write block back
                            // ----------------

                            ExtractBlock(jCell, basis, Degrees, conf, E, _i0s, false, LeftPreCond, ref PCleftBlock[i]);
                            ExtractBlock(jCell, basis, Degrees, conf, E, _i0s, false, RightPreCond, ref PCrightBlock[i]);


                            // inverse precond-matrix
                            // ----------------------

                            // right-inverse: (required for transforming solution guess)
                            if(PCrightBlock_inv[i] == null || PCrightBlock_inv[i].NoOfRows != NN) {
                                PCrightBlock_inv[i] = MultidimensionalArray.Create(NN, NN);
                            }
                            if(Rank == NN)
                                PCrightBlock[i].InvertTo(PCrightBlock_inv[i]);
                            else
                                RankDefInvert(PCrightBlock[i], PCrightBlock_inv[i]);
                            ExtractBlock(jCell, basis, Degrees, conf, E, _i0s, false, RightPreCondInv, ref PCrightBlock_inv[i]);
                                                        
                            // left-inverse: (required for analysis purposes, to transform residuals back onto original grid)
                            if(PCleftBlock_inv[i] == null || PCleftBlock_inv[i].NoOfRows != NN) {
                                PCleftBlock_inv[i] = MultidimensionalArray.Create(NN, NN);
                            }
                            if(Rank == NN)
                                PCleftBlock[i].InvertTo(PCleftBlock_inv[i]);
                            else
                                RankDefInvert(PCleftBlock[i], PCleftBlock_inv[i]);
                            ExtractBlock(jCell, basis, Degrees, conf, E, _i0s, false, LeftPreCondInv, ref PCleftBlock_inv[i]);
                        }
                    }
                }


                return IndefRows.ToArray();
            }
        }


        private static void ExtractBlock(int jCell, 
            AggregationGridBasis basis, int[] Degrees,
            ChangeOfBasisConfig conf, 
            int E, int[] _i0s, bool Sp2Full, 
            IMutableMatrixEx MtxSp, ref MultidimensionalArray MtxFl) {
            
            int NN = conf.VarIndex.Sum(iVar => basis.GetLength(jCell, Degrees[iVar]));
            if(MtxFl == null || MtxFl.NoOfRows != NN) {
                MtxFl = MultidimensionalArray.Create(NN, NN);
            }


            int i0Rowloc = 0;
            for(int eRow = 0; eRow < E; eRow++) { // loop over variables in configuration
                int i0Row = _i0s[eRow];
                int iVarRow = conf.VarIndex[eRow];
                
                int NRow = basis.GetLength(jCell, Degrees[iVarRow]);

                int i0Colloc = 0;
                for(int eCol = 0; eCol < E; eCol++) { // loop over variables in configuration

                    int i0Col = _i0s[eCol];
                    int iVarCol = conf.VarIndex[eCol];
                    
                    int NCol = basis.GetLength(jCell, Degrees[iVarCol]);

                    for(int n_row = 0; n_row < NRow; n_row++) { // row loop...
                        for(int n_col = 0; n_col < NCol; n_col++) { // column loop...
                            if(Sp2Full) {
                                // copy from sparse to full
                                MtxFl[n_row + i0Rowloc, n_col + i0Colloc] = (MtxSp != null) ? ( MtxSp[n_row + i0Row, n_col + i0Col]) : (n_col == n_row ? 1.0 : 0.0);
                            } else {
                                // the other way around.
                                MtxSp[n_row + i0Row, n_col + i0Col] = MtxFl[n_row + i0Rowloc, n_col + i0Colloc];
                            }
                        }
                    }
                    i0Colloc += NCol;
                }
                i0Rowloc += NRow;
            }
        }

        private static int[] ConvertRowIndices(
            int jCell,
            AggregationGridBasis basis, int[] Degrees,
            ChangeOfBasisConfig conf,
            int E, int[] _i0s,
            int[] LocIdx) {

            int NN = conf.VarIndex.Sum(iVar => basis.GetLength(jCell, Degrees[iVar]));
            int[] Loc2glob = new int[NN];

            int i0Rowloc = 0;
            for (int eRow = 0; eRow < E; eRow++) { // loop over variables in configuration
                int i0Row = _i0s[eRow];
                int iVarRow = conf.VarIndex[eRow];

                int NRow = basis.GetLength(jCell, Degrees[iVarRow]);

                for (int n_row = 0; n_row < NRow; n_row++) { // row loop...
                    //n_row = LocIdx[k];

                    int iRowLoc = n_row + i0Rowloc;
                    int iRowGlb = n_row + i0Row;

                    Loc2glob[iRowLoc] = iRowGlb;


                }
                i0Rowloc += NRow;
            }

            return LocIdx.Select(i => Loc2glob[i]).ToArray();
        }



        static void SymmInv(MultidimensionalArray M, MultidimensionalArray L, MultidimensionalArray R) {
            L.Clear();
            R.Clear();
            L.AccEye(1.0);
#if DEBUG
            var Mbefore = M.CloneAs();
#endif

            int n = M.NoOfRows;
            for(int i = 0; i < n; i++) {
                double M_ii = M[i, i];
                if (M_ii == 0.0)
                    throw new ArithmeticException("Zero diagonal element at " + i + "-th row.");
                double scl = 1.0 / Math.Sqrt(Math.Abs(M_ii));
                M.RowScale(i, scl);
                L.RowScale(i, scl);
                M.ColScale(i, scl);

                double diagsign = Math.Sign(M[i, i]);
                if(diagsign == 0.0)
                    throw new ArithmeticException("Zero diagonal element at " + i + "-th row.");
                if (Math.Abs(Math.Abs(M[i, i]) - 1.0) > 1.0e-8)
                    throw new ArithmeticException("Unable to create diagonal 1.0.");

                for(int k = i + 1; k < n; k++) {
                    double M_ki = M[k, i];

                    M.RowAdd(i, k, -M_ki * diagsign);
                    L.RowAdd(i, k, -M_ki * diagsign);
                    M.ColAdd(i, k, -M_ki * diagsign);

                    Debug.Assert(Math.Abs(M[k, i]) < 1.0e-8);
                    Debug.Assert(Math.Abs(M[i, k]) < 1.0e-8);
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

        private static int[] ComputeChangeOfBasisBlock(MultidimensionalArray In_MassMatrixBlock, MultidimensionalArray In_OperatorMatrixBlock, 
            MultidimensionalArray OUT_LeftPC, MultidimensionalArray OUT_rightPC, Mode PCMode, out int Rank,
            MultidimensionalArray work) {
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

                    if(symmErr / infNorm > 1.0e-8)
                        throw new NotSupportedException(string.Format("LDL_DiagBlock is not supported on unsymmetric matrices (Symm-Err: {0:0.####E-00}, Inf-Norm: {1:0.####E-00}, Quotient {2:0.####E-00}).", symmErr, infNorm, symmErr / infNorm));
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
                    var SymmPart = work;
                    In_OperatorMatrixBlock.TransposeTo(SymmPart);
                    SymmPart.Acc(1.0, In_OperatorMatrixBlock);
                    SymmPart.Scale(0.5);
                    
                    int[] ZerosEntries = ModifiedInverseChol(In_MassMatrixBlock, OUT_rightPC, 1.0e-12, false);
                    int NoOfZeros = ZerosEntries == null ? 0 : ZerosEntries.Length;
                    IndefRows = ZerosEntries;
                    Rank = OUT_LeftPC.NoOfCols - NoOfZeros;

                    if(NoOfZeros == 0) {
                        // normal cell -- nix indefinite
                        // +++++++++++++++++++++++++++++

                        SymmInv(SymmPart, OUT_LeftPC, OUT_rightPC);
                    } else {
                        // problem-cell
                        // ++++++++++++++

                        OUT_rightPC.TransposeTo(OUT_LeftPC);

                        SymmPart = IMatrixExtensions.GEMM(OUT_LeftPC, SymmPart, OUT_rightPC);

                        int[] ZerosEntries2 = ModifiedInverseChol(SymmPart, OUT_rightPC, 1.0e-12, true);
                        OUT_rightPC.TransposeTo(OUT_LeftPC);

                        if(!ZerosEntries2.IsSetEqual(ZerosEntries))
                            throw new ArithmeticException();

                        
                        break;

                    }

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
                throw new NotImplementedException();
            }

            return IndefRows;
        }


        /// <summary>
        /// Partial inversion of a matrix with zero rows and columns; not considered to be performance-critical,
        /// since it is only for treating phatological cases.
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

                if(row_i.L2NormPow2() == 0.0) {
                    double[] col_i = MtxIn.GetColumn(i);

                    if(col_i.L2NormPow2() != 0.0)
                        throw new ArgumentException("Row is zero, but column is not.");

                } else {
                    NonzeroIdx.Add(i);
                }
            }

            if(NonzeroIdx.Count == M)
                throw new ArgumentException("Unable to find zero row/column.");

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
