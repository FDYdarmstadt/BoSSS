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
using System.Diagnostics;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using BoSSS.Foundation;
using ilPSP.Tracing;

namespace BoSSS.Solution.Multigrid {


    public partial class MultigridOperator {

        /// <summary>
        /// DG coordinate mapping on the original grid/mesh.
        /// </summary>
        public UnsetteledCoordinateMapping BaseGridProblemMapping {
            get;
            private set;
        }

        /// <summary>
        /// Recursive constructor, i.e. constructs operators on all multigrid levels which are provided by <paramref name="basisSeq"/>.
        /// </summary>
        /// <param name="basisSeq"></param>
        /// <param name="_ProblemMapping">
        /// DG coordinate mapping on the original grid, see <see cref="BaseGridProblemMapping"/>.
        /// </param>
        /// <param name="OperatorMatrix"></param>
        /// <param name="MassMatrix">
        /// Mass matrix on the original grid. If null, the identity matrix is assumed. It is only required if the 
        /// </param>
        /// <param name="cobc">
        /// Configuration of the cell-wise, explicit block-preconditioning for each multigrid level.
        /// (Remark: this kind of preconditioning is mathematically equivalent to a change of the DG resp. XDG basis.)
        /// </param>
        public MultigridOperator(IEnumerable<AggregationGridBasis[]> basisSeq,
            UnsetteledCoordinateMapping _ProblemMapping, BlockMsrMatrix OperatorMatrix, BlockMsrMatrix MassMatrix,
            IEnumerable<ChangeOfBasisConfig[]> cobc)
            : this(null, basisSeq, _ProblemMapping, cobc) //
        {
            if (!OperatorMatrix.RowPartitioning.EqualsPartition(_ProblemMapping))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!OperatorMatrix.ColPartition.EqualsPartition(_ProblemMapping))
                throw new ArgumentException("Column partitioning mismatch.");


            if (MassMatrix != null) {
                if (!MassMatrix.RowPartitioning.Equals(_ProblemMapping))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!MassMatrix.ColPartition.Equals(_ProblemMapping))
                    throw new ArgumentException("Column partitioning mismatch.");
            }

            if (this.LevelIndex == 0) {

                if (this.IndexIntoProblemMapping_Local == null) {
                    this.m_RawOperatorMatrix = OperatorMatrix;

                    if (MassMatrix != null)
                        this.m_RawMassMatrix = MassMatrix;  //is BlockMsrMatrix) ? ((BlockMsrMatrix)MassMatrix) : MassMatrix.ToMsrMatrix();
                    else
                        this.m_RawMassMatrix = null;
                } else {
                    this.m_RawOperatorMatrix = new BlockMsrMatrix(this.Mapping, this.Mapping);
                    OperatorMatrix.WriteSubMatrixTo(this.m_RawOperatorMatrix, this.IndexIntoProblemMapping_Global, default(int[]), this.IndexIntoProblemMapping_Global, default(int[]));

                    if (MassMatrix != null) {
                        this.m_RawMassMatrix = new BlockMsrMatrix(this.Mapping, this.Mapping);
                        BlockMsrMatrix MMR = MassMatrix;
                        MMR.WriteSubMatrixTo(this.m_RawMassMatrix, this.IndexIntoProblemMapping_Global, default(int[]), this.IndexIntoProblemMapping_Global, default(int[]));
                    } else {
                        this.m_RawMassMatrix = null;
                    }
                }
            }
        }

        BlockMsrMatrix m_RawOperatorMatrix = null;
        BlockMsrMatrix m_RawMassMatrix = null;

        bool setupdone = false;

        void Setup() {
            using (new FuncTrace()) {
                if (setupdone)
                    return;
                setupdone = true;

                if (this.FinerLevel != null)
                    this.FinerLevel.Setup();


                // Construct intermediate 'raw' restriction and prolongation operators
                // ===================================================================

                BlockMsrMatrix RawRestriction, RawProlongation;
                if (this.FinerLevel == null) {
                    RawRestriction = null;
                    RawProlongation = null;
                } else {
//#if DEBUG
//                    var __PrlgOperator_Check = this.FinerLevel.Mapping.FromOtherLevelMatrix(this.Mapping);
//#endif
                    var __PrlgOperator = this.Mapping.GetProlongationOperator(this.FinerLevel.Mapping);
                    Debug.Assert(__PrlgOperator.RowPartitioning.LocalLength == this.FinerLevel.Mapping.LocalLength);
                    Debug.Assert(__PrlgOperator.ColPartition.LocalLength == this.Mapping.LocalLength);
//#if DEBUG
//                    var __Err = __PrlgOperator_Check.CloneAs();
//                    __Err.Acc(-1.0, __PrlgOperator);
//                    double ErrNorm = __Err.InfNorm();
//                    Console.WriteLine("Error norm: " + ErrNorm);
//                    double Ref = Math.Max(__PrlgOperator.InfNorm(), __PrlgOperator_Check.InfNorm());
//                    //Debug.Assert(ErrNorm < Ref*1.0e-8);  
//#endif


                    if (this.FinerLevel.RightChangeOfBasis_Inverse != null)
                        RawProlongation = BlockMsrMatrix.Multiply(this.FinerLevel.RightChangeOfBasis_Inverse, __PrlgOperator);
                    else
                        RawProlongation = __PrlgOperator;

                    RawRestriction = RawProlongation.Transpose();
                }
                this.m_PrologateOperator = RawProlongation;
                this.m_RestrictionOperator = RawRestriction;


                // Construct intermediate 'raw' operator and mass matrix (before change of basis)
                // ==============================================================================

                // operator matrix before change of basis
                BlockMsrMatrix RawOpMatrix;
                if (this.FinerLevel == null) {
                    RawOpMatrix = this.m_RawOperatorMatrix;
                    this.m_RawOperatorMatrix = null;
                } else {
                    BlockMsrMatrix Op = FinerLevel.OperatorMatrix;

                    RawOpMatrix = BlockMsrMatrix.Multiply(RawRestriction, BlockMsrMatrix.Multiply(Op, RawProlongation));
                }

                // mass matrix before change of basis
                BlockMsrMatrix RawMassMatrix;
                if (this.FinerLevel == null) {
                    RawMassMatrix = this.m_RawMassMatrix;
                    this.m_RawMassMatrix = null;
                } else {
                    BlockMsrMatrix MM = FinerLevel.MassMatrix;

                    RawMassMatrix = BlockMsrMatrix.Multiply(RawRestriction, BlockMsrMatrix.Multiply(MM, RawProlongation));
                }

                Debug.Assert(RawOpMatrix.RowPartitioning.LocalLength == this.Mapping.LocalLength);

                // compute change of basis
                // =======================
                int[] IndefRows = this.ComputeChangeOfBasis(RawOpMatrix, RawMassMatrix, out m_LeftChangeOfBasis, out m_RightChangeOfBasis, out m_LeftChangeOfBasis_Inverse, out m_RightChangeOfBasis_Inverse);

                // apply change of basis to operator matrix
                // ========================================
                {
                    var Lpc = this.m_LeftChangeOfBasis;
                    var Rpc = this.m_RightChangeOfBasis;

                    BlockMsrMatrix O1;
                    if (Lpc != null) {
                        O1 = BlockMsrMatrix.Multiply(Lpc, RawOpMatrix);
                    } else {
                        O1 = RawOpMatrix;
                    }

                    if (Rpc != null) {
                        this.m_OperatorMatrix = BlockMsrMatrix.Multiply(O1, Rpc);
                    } else {
                        this.m_OperatorMatrix = O1;
                    }


                    // fix zero rows 
                    // (possible result from the Mode.IdMass_DropIndefinite or Mode.SymPart_DiagBlockEquilib_DropIndefinite -- option)

                    foreach (int _iRow in IndefRows) {
                        int iRow = Math.Abs(_iRow);
                        Debug.Assert(this.m_OperatorMatrix.GetNoOfNonZerosPerRow(iRow) == 0);
                        this.m_OperatorMatrix[iRow, iRow] = 1.0;
                    }
#if DEBUG
                for (int iRow = this.m_OperatorMatrix.RowPartitioning.i0; iRow < this.m_OperatorMatrix.RowPartitioning.iE; iRow++) {
                    Debug.Assert(this.m_OperatorMatrix.GetNoOfNonZerosPerRow(iRow) > 0);
                }
#endif
                }

                // apply change of basis to mass matrix
                // ====================================
                {
                    var Lpc = this.m_LeftChangeOfBasis;
                    var Rpc = this.m_RightChangeOfBasis;

                    if (RawMassMatrix != null) {

                        BlockMsrMatrix O1;
                        if (Lpc != null) {
                            O1 = BlockMsrMatrix.Multiply(Lpc, RawMassMatrix);
                        } else {
                            O1 = RawMassMatrix;
                        }

                        if (Rpc != null) {
                            this.m_MassMatrix = BlockMsrMatrix.Multiply(O1, Rpc);
                        } else {
                            this.m_MassMatrix = O1;
                        }

                    } else {
                        if (Lpc != null && Rpc != null) {
                            this.m_MassMatrix = BlockMsrMatrix.Multiply(Lpc, Rpc);
                        } else if (Lpc != null) {
                            Debug.Assert(Rpc == null);
                            this.m_MassMatrix = Lpc;
                        } else if (Rpc != null) {
                            Debug.Assert(Lpc == null);
                            this.m_MassMatrix = Rpc;
                        } else {
                            Debug.Assert(Lpc == null);
                            Debug.Assert(Rpc == null);
                            this.m_MassMatrix = null;
                        }

                    }

                    //{
                    //    var eyeTest = this.m_MassMatrix.CloneAs();
                    //    eyeTest.AccEyeSp(-1.0);
                    //    double infnrm = eyeTest.InfNorm();
                    //    Console.WriteLine("Mass matrix level {0} is id {1} ", this.Mapping.LevelIndex, infnrm);
                    //}
                }
            }
        }

        BlockMsrMatrix m_PrologateOperator = null;
        BlockMsrMatrix m_RestrictionOperator = null;

        /// <summary>
        /// Restricts a right-hand-side or solution vector <paramref name="IN_fine"/> 
        /// from the finer multi-grid level (see <see cref="FinerLevel"/>)
        /// to this multigrid level.
        /// </summary>
        /// <param name="IN_fine">Input: vector on finer level.</param>
        /// <param name="OUT_coarse">Output: vector on this level, i.e. the coarser level.</param>
        public void Restrict<T1, T2>(T1 IN_fine, T2 OUT_coarse)
            where T1 : IList<double>
            where T2 : IList<double> {
            if (this.FinerLevel == null)
                throw new NotSupportedException("Already on finest level -- no finer level to restrict from.");

            if (IN_fine.Count != this.FinerLevel.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of fine grid vector (input).", "IN_fine");
            if (OUT_coarse.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of coarse grid vector (output).", "OUT_coarse");

            Setup();

            this.m_RestrictionOperator.SpMV(1.0, IN_fine, 0.0, OUT_coarse);

            if (this.LeftChangeOfBasis != null) {
                double[] LB = new double[OUT_coarse.Count];
                this.LeftChangeOfBasis.SpMV(1.0, OUT_coarse, 0.0, LB);
                OUT_coarse.SetV(LB);
            }


            /*
            var M = MsrMatrix.Multiply(this.RightChangeOfBasis_Inverse, this.PrologateOperator);
            var MaMa = MsrMatrix.Multiply(M.Transpose(), M);
            int JAGG = this.Mapping.AggGrid.NoOfAggregateCells;
            MsrMatrix Lpc, Rpc, RpcI;
            this.CoarserLevel.ComputeChangeOfBasis(MaMa, MaMa, out Lpc, out Rpc, out RpcI);

            //var prl = MsrMatrix.Multiply(M, Rpc);
            //var rest = prl.Transpose();
            var rest = MsrMatrix.Multiply(Lpc, M.Transpose());

            rest.SpMVpara(1.0, IN_fine, 0.0, OUT_coarse);
             */
        }



        /// <summary>
        /// Prolongates a solution vector <paramref name="IN_coarse"/> from the this multi-grid level 
        /// to the finer multigrid level (see <see cref="FinerLevel"/>).
        /// </summary>
        /// <param name="OUT_fine">Input/Output: accumulation vector on finer level</param>
        /// <param name="IN_coarse">Input: vector on this level, i.e. the coarser level.</param>
        /// <param name="beta">
        /// scaling of the accumulation vector (see remarks).
        /// </param>
        /// <param name="alpha">
        /// scaling applied to prolongated vector (see remarks).
        /// </param>
        /// <remarks>
        /// On exit,
        /// <paramref name="OUT_fine"/> = <paramref name="OUT_fine"/>*<paramref name="beta"/> + Pr(<paramref name="IN_coarse"/>)*<paramref name="alpha"/>.
        /// </remarks>
        public void Prolongate<T1, T2>(double alpha, T1 OUT_fine, double beta, T2 IN_coarse)
            where T1 : IList<double>
            where T2 : IList<double> {
            if (this.FinerLevel == null)
                throw new NotSupportedException("Already on finest level -- no finer level to prolongate to.");

            if (OUT_fine.Count != this.FinerLevel.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of fine grid vector (Output)", "OUT_fine");
            if (IN_coarse.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of coarse grid vector (Input)", "IN_coarse");

            Setup();

            if (this.RightChangeOfBasis != null) {
                double[] RX = new double[IN_coarse.Count];
                this.RightChangeOfBasis.SpMV(1.0, IN_coarse, 0.0, RX);

                this.m_PrologateOperator.SpMV(alpha, RX, beta, OUT_fine);
            } else {
                this.m_PrologateOperator.SpMV(alpha, IN_coarse, beta, OUT_fine);
            }
        }

        /// <summary>
        /// Returns the index of this multigrid level. Smaller indices correspont to finer grids.
        /// </summary>
        public int LevelIndex {
            get {
                if (this.FinerLevel == null)
                    return 0;
                else
                    return this.FinerLevel.LevelIndex + 1;
            }
        }


        /// <summary>
        /// only used on top level:
        /// mapping: 'multigrid mapping index' --> 'full problem index'
        ///  - index: compressed index at this multigrid level
        ///  - content: index into <see cref="BaseGridProblemMapping"/> (local)
        /// </summary>
        int[] IndexIntoProblemMapping_Local;

        /// <summary>
        /// only used on top level:
        /// mapping: 'multigrid mapping index' --> 'full problem index'
        ///  - index: compressed index at this multigrid level
        ///  - content: index into <see cref="BaseGridProblemMapping"/> (global)
        /// </summary>
        int[] IndexIntoProblemMapping_Global;


        private MultigridOperator(MultigridOperator __FinerLevel, IEnumerable<AggregationGridBasis[]> basisES, UnsetteledCoordinateMapping _pm, IEnumerable<ChangeOfBasisConfig[]> cobc) {
            if (basisES.Count() <= 0) {
                throw new ArgumentException("At least one multigrid level is required.");
            }

            this.BaseGridProblemMapping = _pm;
            if (cobc.Count() < 1)
                throw new ArgumentException();
            this.m_Config = cobc.First();
            this.FinerLevel = __FinerLevel;
            this.Mapping = new MultigridMapping(_pm, basisES.First(), this.Degrees);

            if (this.Mapping.LocalLength > this.BaseGridProblemMapping.LocalLength)
                throw new ApplicationException("Something wrong.");


            if (basisES.Count() > 1) {
                this.CoarserLevel = new MultigridOperator(this, basisES.Skip(1), _pm, cobc.Count() > 1 ? cobc.Skip(1) : cobc);
            }

            if (this.LevelIndex == 0 && this.Mapping.AggBasis.Any(agb => agb.ReqModeIndexTrafo)) {
                int J = this.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                Debug.Assert(J == this.BaseGridProblemMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells);
                IndexIntoProblemMapping_Local = new int[this.Mapping.LocalLength];
#if DEBUG
                IndexIntoProblemMapping_Local.SetAll(-23456);
#endif

                AggregationGridBasis[] AgBss = this.Mapping.AggBasis;
                int[] Degrees = this.Mapping.DgDegree;
                int NoFlds = Degrees.Length;

                for (int j = 0; j < J; j++) { // loop over cells

                    for (int iFld = 0; iFld < NoFlds; iFld++) {
                        int N = AgBss[iFld].GetLength(j, Degrees[iFld]); // By using the length of the aggregate grid mapping,
                        //                                            we exclude XDG agglomerated cells.

                        if (N > 0) {
                            int k0 = this.Mapping.ProblemMapping.LocalUniqueCoordinateIndex(iFld, j, 0);
                            int i0 = this.Mapping.LocalUniqueIndex(iFld, j, 0);
                            for (int n = 0; n < N; n++) {
                                //X[i0 + n] = INOUT_X[k0 + n];
                                //B[i0 + n] = IN_RHS[k0 + n];

                                int n_trf = AgBss[iFld].N_Murks(j, n, N);
                                Debug.Assert(n_trf >= 0);
#if DEBUG
                                Debug.Assert(IndexIntoProblemMapping_Local[i0 + n] < 0);

#endif
                                IndexIntoProblemMapping_Local[i0 + n] = k0 + n_trf;
                                //Debug.Assert(i0 + n == 0 || IndexIntoProblemMapping_Local[i0 + n] > IndexIntoProblemMapping_Local[i0 + n - 1]);
                                Debug.Assert(IndexIntoProblemMapping_Local[i0 + n] >= 0);
                                Debug.Assert(IndexIntoProblemMapping_Local[i0 + n] < this.Mapping.ProblemMapping.LocalLength);
                            }
                        } else {
                            // should be a cell without any species.

                            Debug.Assert(this.Mapping.AggBasis.Where(b => !(b is XdgAggregationBasis)).Count() == 0, "all must be XDG");
                            Debug.Assert(this.Mapping.AggBasis.Where(b => ((XdgAggregationBasis)b).GetNoOfSpecies(j) != 0).Count() == 0, "no species in any cell allowed");

                            //Debug.Assert();
                        }
                    }
                }

#if DEBUG
                Debug.Assert(IndexIntoProblemMapping_Local.Any(i => i < 0) == false);
#endif

                if (this.Mapping.ProblemMapping.MpiSize == 1) {
                    this.IndexIntoProblemMapping_Global = this.IndexIntoProblemMapping_Local;
                } else {
                    this.IndexIntoProblemMapping_Global = this.IndexIntoProblemMapping_Local.CloneAs();
                    int i0Proc = this.Mapping.ProblemMapping.i0;
                    int L = this.IndexIntoProblemMapping_Global.Length;
                    for (int i = 0; i < L; i++) {
                        this.IndexIntoProblemMapping_Global[i] += i0Proc;
                    }
                }
            }
#if DEBUG
            if (IndexIntoProblemMapping_Local != null) {
                int i0 = this.Mapping.ProblemMapping.i0;
                int L = this.Mapping.ProblemMapping.LocalLength;

                int[] UseCount = new int[L];
                foreach (int i in IndexIntoProblemMapping_Global) {
                    UseCount[i - i0]++;
                }

                for (int l = 0; l < L; l++) {
                    Debug.Assert(UseCount[l] <= 1);
                }
            }
#endif
        }

        /// <summary>
        /// DG coordinate mapping which corresponds to this mesh level, resp. operator matrix (<see cref="OperatorMatrix"/>, <see cref="MassMatrix"/>).
        /// </summary>
        public MultigridMapping Mapping {
            get;
            private set;
        }

        /// <summary>
        /// Grid/mesh on which this operator is defined.
        /// </summary>
        public BoSSS.Foundation.Grid.IGridData GridData {
            get {
                return Mapping.AggGrid;
            }
        }

        
        /// <summary>
        /// Pointer to operator on finer level.
        /// </summary>
        public MultigridOperator FinerLevel {
            get;
            private set;
        }
        
        /// <summary>
        /// Pointer to operator on coarser level.
        /// </summary>
        public MultigridOperator CoarserLevel {
            get;
            private set;
        }

        BlockMsrMatrix m_LeftChangeOfBasis;

        public BlockMsrMatrix LeftChangeOfBasis {
            get {
                Setup();
                return m_LeftChangeOfBasis;
            }
        }

        BlockMsrMatrix m_RightChangeOfBasis;

        public BlockMsrMatrix RightChangeOfBasis {
            get {
                Setup();
                return m_RightChangeOfBasis;
            }
        }

        BlockMsrMatrix m_RightChangeOfBasis_Inverse;

        public BlockMsrMatrix RightChangeOfBasis_Inverse {
            get {
                Setup();
                return m_RightChangeOfBasis_Inverse;
            }
        }

        BlockMsrMatrix m_LeftChangeOfBasis_Inverse;

        public BlockMsrMatrix LeftChangeOfBasis_Inverse {
            get {
                Setup();
                return m_LeftChangeOfBasis_Inverse;
            }
        }


        BlockMsrMatrix m_OperatorMatrix = null;
        BlockMsrMatrix m_MassMatrix = null;
        
        /// <summary>
        /// Returns the Operator matrix on the current multigrid level.
        /// </summary>
        public BlockMsrMatrix OperatorMatrix {
            get {
                Setup();
               
                return m_OperatorMatrix;
            }
        }

        /// <summary>
        /// Returns the mass matrix on this multigrid level;
        /// an identity matrix is encoded as null.
        /// </summary>
        public BlockMsrMatrix MassMatrix {
            get {
                Setup();
                return m_MassMatrix;
            }
        }


        public void UseSolver<T1, T2>(ISolverSmootherTemplate solver, T1 INOUT_X, T2 IN_RHS, bool UseGuess = true)
            where T1 : IList<double>
            where T2 : IList<double> 
        {
            if(this.LevelIndex != 0)
                throw new NotSupportedException("Not Inteded to be called on any multi-grid level but the finest one.");

            int I = this.Mapping.ProblemMapping.LocalLength;
            if(INOUT_X.Count != I)
                throw new ArgumentException("Vector length mismatch.", "INOUT_X");
            if(IN_RHS.Count != I)
                throw new ArgumentException("Vector length mismatch.", "IN_RHS");

            if(this.FinerLevel != null)
                throw new NotSupportedException("This method may only be called on the top level.");

            //if(this.Mapping.AggGrid.NoOfAggregateCells != this.Mapping.ProblemMapping.GridDat.Cells.NoOfCells)
            //    throw new ArgumentException();
            //int J = this.Mapping.AggGrid.NoOfAggregateCells;
            
            
            int L = this.Mapping.LocalLength;
            double[] X = new double[L];
            double[] B = new double[L];
            this.TransformRhsInto(IN_RHS, B);
            if(UseGuess)
                this.TransformSolInto(INOUT_X, X);

            solver.ResetStat();
            solver.Solve(X, B);

            this.TransformSolFrom(INOUT_X, X);
        }

        public void TransformSolInto<T1, T2>(T1 u_IN, T2 v_OUT)
            where T1 : IList<double>
            where T2 : IList<double> 
        {
            if(this.FinerLevel != null)
                throw new NotSupportedException("Only supported on finest level.");
            if(u_IN.Count != this.Mapping.ProblemMapping.LocalLength)
                throw new ArgumentException("Mismatch in length of input vector.", "u");
            if(v_OUT.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of output vector.", "v");

            int L = this.Mapping.LocalLength;
            double[] uc = new double[L];
            uc.AccV(1.0, u_IN, default(int[]), this.IndexIntoProblemMapping_Local);

            if(this.RightChangeOfBasis_Inverse != null) {
                this.RightChangeOfBasis_Inverse.SpMV(1.0, uc, 0.0, v_OUT);
            } else {
                v_OUT.SetV(uc);
            }

        }

        public void TransformRhsInto<T1, T2>(T1 u_IN, T2 v_OUT)
            where T1 : IList<double>
            where T2 : IList<double> 
        {
            if(this.FinerLevel != null)
                throw new NotSupportedException("Only supported on finest level.");
            if(u_IN.Count != this.Mapping.ProblemMapping.LocalLength)
                throw new ArgumentException("Mismatch in length of input vector.", "u");
            if(v_OUT.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of output vector.", "v");

            int L = this.Mapping.LocalLength;
            double[] uc = new double[L];
            uc.AccV(1.0, u_IN, default(int[]), this.IndexIntoProblemMapping_Local);

            if(this.LeftChangeOfBasis != null) {
                this.LeftChangeOfBasis.SpMV(1.0, uc, 0.0, v_OUT);
            } else {
                v_OUT.SetV(uc);
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="u_Out">output</param>
        /// <param name="v_In">input</param>
        public void TransformSolFrom<T1, T2>(T1 u_Out, T2 v_In)
            where T1 : IList<double>
            where T2 : IList<double> //
        {
            if(this.FinerLevel != null)
                throw new NotSupportedException("Only supported on finest level.");
            if(u_Out.Count != this.Mapping.ProblemMapping.LocalLength)
                throw new ArgumentException("Mismatch in length of output vector.", "u");
            if(v_In.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of input vector.", "v");

            int L = this.Mapping.LocalLength;
            double[] uc = new double[L];
            
            if(this.RightChangeOfBasis != null) {
                this.RightChangeOfBasis.SpMV(1.0, v_In, 0.0, uc);
            } else {
                uc.SetV(v_In);
            }
            
            u_Out.ClearEntries();
            u_Out.AccV(1.0, uc, this.IndexIntoProblemMapping_Local, default(int[]));
        }
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u_Out">output</param>
        /// <param name="v_In">input</param>
        public void TransformRhsFrom<T1, T2>(T1 u_Out, T2 v_In)
            where T1 : IList<double>
            where T2 : IList<double>  //
        {
            if(this.FinerLevel != null)
                throw new NotSupportedException("Only supported on finest level.");
            if(u_Out.Count != this.Mapping.ProblemMapping.LocalLength)
                throw new ArgumentException("Mismatch in length of output vector.", "u");
            if(v_In.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of input vector.", "v");

            int L = this.Mapping.LocalLength;
            double[] uc = new double[L];

            if(this.LeftChangeOfBasis_Inverse != null) {
                this.LeftChangeOfBasis_Inverse.SpMV(1.0, v_In, 0.0, uc);
            } else {
                uc.SetV(v_In);
            }

            u_Out.ClearEntries();
            u_Out.AccV(1.0, uc, this.IndexIntoProblemMapping_Local, default(int[]));
        }
    }
}
