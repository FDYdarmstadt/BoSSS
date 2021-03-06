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
using MPI.Wrappers;
using BoSSS.Foundation.XDG;
using System.Collections;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.AdvancedSolvers {


    public partial class MultigridOperator {

        internal static LevelSetTracker GetTracker(UnsetteledCoordinateMapping map) {
            LevelSetTracker lsTrk = null;

            foreach (Basis b in map.BasisS) {
                if (b != null && b is XDGBasis xb) {
                    if (lsTrk == null) {
                        lsTrk = xb.Tracker;
                    } else {
                        if (!object.ReferenceEquals(lsTrk, xb.Tracker))
                            throw new ArgumentException("Tracker mismatch.");
                    }
                }
            }
            return lsTrk;
        }

        static private long[][] GetGlobalCellNeigbourship(UnsetteledCoordinateMapping map) {
            var GlobalNumberOfCells = map.GridDat.CellPartitioning.TotalLength;
            long[] externalCellsGlobalIndices = map.GridDat.iParallel.GlobalIndicesExternalCells;
            long[][] globalCellNeigbourship = new long[GlobalNumberOfCells][];
            int J = map.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            for (long j = 0; j < J; j++) {
                long globalIndex = j + map.GridDat.CellPartitioning.i0;
                // we use GetCellNeighboursViaEdges(j) to also find neigbours at periodic boundaries
                Tuple<int, int, int>[] cellNeighbours = map.GridDat.GetCellNeighboursViaEdges((int)j);
                globalCellNeigbourship[globalIndex] = new long[cellNeighbours.Length];
                for (int i = 0; i < cellNeighbours.Length; i++) {
                    globalCellNeigbourship[globalIndex][i] = cellNeighbours[i].Item1;
                }

                // translate local neighbour index into global index
                for (int i = 0; i < globalCellNeigbourship[globalIndex].Length; i++) {
                    if (globalCellNeigbourship[globalIndex][i] < J)
                        globalCellNeigbourship[globalIndex][i] = globalCellNeigbourship[globalIndex][i] + map.GridDat.CellPartitioning.i0;
                    else
                        globalCellNeigbourship[globalIndex][i] = (int)externalCellsGlobalIndices[globalCellNeigbourship[globalIndex][i] - J];
                }
            }

            for (int globalIndex = 0; globalIndex < GlobalNumberOfCells; globalIndex++) {
                int processID = !globalCellNeigbourship[globalIndex].IsNullOrEmpty() ? map.GridDat.MpiRank : 0;
                processID = processID.MPIMax();
                globalCellNeigbourship[globalIndex] = globalCellNeigbourship[globalIndex].MPIBroadcast(processID);
            }

            return globalCellNeigbourship;
        }

        /// <summary>
        /// Returns the cutcells + neighbours on a global level.
        /// </summary>
        /// <param name="globalCellNeighbourship">
        /// </param>
        static private BitArray GetGlobalNearBand(BitArray levelSetCells, UnsetteledCoordinateMapping map) {
            long[][] globalCellNeighbourship = GetGlobalCellNeigbourship(map);
            BitArray globalCutCells = new BitArray(checked((int)map.GridDat.CellPartitioning.TotalLength));
            int J = map.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < J; j++) {
                if (levelSetCells[j]) {
                    int globalIndex = checked((int)(j + map.GridDat.CellPartitioning.i0));
                    globalCutCells[globalIndex] = true;
                    for (int i = 0; i < globalCellNeighbourship[globalIndex].Length; i++) {
                        globalCutCells[(int)globalCellNeighbourship[globalIndex][i]] = true;
                    }
                }
            }

            for (int j = 0; j < globalCutCells.Length; j++) {
                globalCutCells[j] = globalCutCells[j].MPIOr();
            }
            return globalCutCells;
        }

        int FindPhaseDGCoordinate(UnsetteledCoordinateMapping map, int iVar, int jCell, AggregationGridBasis[] bases) {
            LevelSetTracker lsTrk = GetTracker(map);
            var basis = bases[iVar];
            if(lsTrk != null) {
                int N = basis.GetLength(jCell, map.BasisS[iVar].Degree);
                return basis.N_Murks(jCell, 0, N);
            } else {
                return 0;
            }
        }

        static long FindReferencePointCell(UnsetteledCoordinateMapping map, AggregationGridBasis[] bases) {
            int J = map.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            LevelSetTracker lsTrk = GetTracker(map);
            BitArray Cells2avoid;
            int neighborSearchDepth = 4;
            int jFound = -1;
            bool foundACell = false;
            //Debugger.Launch();
            while (!foundACell && neighborSearchDepth >= 0) {
                if (lsTrk != null) {
                    Cells2avoid = lsTrk.Regions.GetNearFieldMask(Math.Min(1, neighborSearchDepth)).GetBitMask();
                    for (int i = 0; i < neighborSearchDepth - 2; i++) {
                        if (neighborSearchDepth - 2 > 0)
                            Cells2avoid = GetGlobalNearBand(Cells2avoid.CloneAs(), map);
                    }
                } else {
                    Cells2avoid = null;
                }
                for (int j = 0; j < J; j++) {
                    if (bases[0].GetLength(j, 0) > 0 && bases[0].GetNoOfSpecies(j) == 1 && (Cells2avoid == null || Cells2avoid[j] == false)) {
                        jFound = j;
                        break;
                    }
                }
                foundACell = jFound.MPIMax() >= 0;
                neighborSearchDepth -= 1;
            }
            

            jFound += (int)map.GridDat.CellPartitioning.i0;

            long jFoundGlob = jFound.MPIMax();
            if(jFoundGlob < 0)
                throw new ApplicationException("unable to find reference cell.");
            return jFoundGlob;
        }

        /// <summary>
        /// global index
        /// cell in which the reference point is located
        /// </summary>
        long m_ReferenceCell;

        /// <summary>
        /// Global Indices into <see cref="BaseGridProblemMapping"/>
        /// </summary>
        long[] m_ReferenceIndices;

        void DefineReferenceIndices() {
            UnsetteledCoordinateMapping map = this.BaseGridProblemMapping;
            AggregationGridBasis[] bases = this.Mapping.AggBasis;


            if (map.BasisS.Count != bases.Length)
                throw new ArgumentException();
            if(bases.Length != FreeMeanValue.Length)
                throw new ArgumentException(); 

            if(FreeMeanValue.Any(b => b) == false) {
                // none of the solution variables contains a "free mean value"
                // => no need to do anything further
                return;
            }
            var asd = bases[0].DGBasis;
            int L = bases.Length;
            m_ReferenceCell = FindReferencePointCell(map, bases);
            bool onthisProc = BaseGridProblemMapping.GridDat.CellPartitioning.IsInLocalRange(m_ReferenceCell);

            if(onthisProc) {
                int jRefLoc = BaseGridProblemMapping.GridDat.CellPartitioning.TransformIndexToLocal(m_ReferenceCell);


                m_ReferenceIndices = new long[L];
                for(int iVar = 0; iVar < L; iVar++) {
                    if(FreeMeanValue[iVar]) {
                        int kCoord = FindPhaseDGCoordinate(map, iVar, jRefLoc, bases);
                        m_ReferenceIndices[iVar] = BaseGridProblemMapping.GlobalUniqueCoordinateIndex(iVar, jRefLoc, kCoord);
                    } else {
                        m_ReferenceIndices[iVar] = int.MinValue;
                    }
                }
            } else {
                m_ReferenceIndices = null;
            }


            int originRank = BaseGridProblemMapping.GridDat.CellPartitioning.FindProcess(m_ReferenceCell); // on this rank, the 'm_ReferenceIndices' are defined
            m_ReferenceIndices = m_ReferenceIndices.MPIBroadcast(originRank);

        }


        /// <summary>
        /// modifies a right-hand-side <paramref name="rhs"/>
        /// in order to fix the pressure at some reference point
        /// </summary>
        public (int idx, double val)[] SetPressureReferencePointRHS<T>(T rhs)
            where T : IList<double> {
            using (new FuncTrace()) {
                if(this.LevelIndex != 0)
                    throw new NotSupportedException("Can only be invoked on top level.");

                if(m_ReferenceIndices == null)
                    return new (int idx, double val)[0]; // nothing to do

                if (rhs.Count != BaseGridProblemMapping.LocalLength)
                    throw new ArgumentException("vector length mismatch");

                bool onthisProc = BaseGridProblemMapping.GridDat.CellPartitioning.IsInLocalRange(m_ReferenceCell);


                // clear row
                // ---------

                var bkup = new List<(int idx, double val)>();

                if (onthisProc) {

                    for(int iVar = 0; iVar < m_ReferenceIndices.Length; iVar++) {
                        // clear RHS
                        if(m_ReferenceIndices[iVar] >= 0) {
                            int iRowLoc = BaseGridProblemMapping.TransformIndexToLocal(m_ReferenceIndices[iVar]);
                            bkup.Add((iRowLoc, rhs[iRowLoc]));
                            rhs[iRowLoc] = 0;
                        }
                    }
                }

                return bkup.ToArray();
            }
        }


        /// <summary>
        /// modifies a matrix 
        /// in order to fix the pressure at some reference point
        /// </summary>
        public (long iRow, long jCol, double Val)[] SetPressureReferencePointMTX(IMutableMatrixEx Mtx) {
            using (new FuncTrace()) {
                
                if(m_ReferenceIndices == null)
                    return new (long iRow, long jCol, double Val)[0]; // nothing to do
                bool onthisProc = BaseGridProblemMapping.GridDat.CellPartitioning.IsInLocalRange(m_ReferenceCell);


                var map = this.BaseGridProblemMapping;
                if (!Mtx.RowPartitioning.EqualsPartition(map) || !Mtx.ColPartition.EqualsPartition(map))
                    throw new ArgumentException();


                // clear row
                // ---------

                var bkup = new List<(long, long, double)>();

                if (onthisProc) {

                    for(int iVar = 0; iVar < m_ReferenceIndices.Length; iVar++) {
                        // clear RHS
                        if(m_ReferenceIndices[iVar] >= 0) {
                            long iRowGl = m_ReferenceIndices[iVar];

                            // set matrix row to identity

                            long[] ColIdx = Mtx.GetOccupiedColumnIndices(iRowGl);
                            foreach(int ci in ColIdx) {
                                bkup.Add((iRowGl, ci, Mtx[iRowGl, ci]));
                                Mtx[iRowGl, ci] = 0;
                            }
                            Mtx.SetDiagonalElement(iRowGl, 1.0);


                        }
                    }
                }

                // clear column
                // ------------
                {
                    for(int iVar = 0; iVar < m_ReferenceIndices.Length; iVar++) {
                        if(m_ReferenceIndices[iVar] >= 0) {
                            long iRowGl = m_ReferenceIndices[iVar];
                            for(long i = Mtx.RowPartitioning.i0; i < Mtx.RowPartitioning.iE; i++) {
                                
                                if(i != iRowGl) {
                                    double a = Mtx[i, iRowGl];
                                    if(a != 0.0) {
                                        bkup.Add((i, iRowGl, a));
                                        Mtx[i, iRowGl] = 0;
                                    }
                                }
                            }
                        }
                    }
                }

                return bkup.ToArray();
            }
        }



        /// <summary>
        /// DG coordinate mapping on the original grid/mesh.
        /// </summary>
        public UnsetteledCoordinateMapping BaseGridProblemMapping {
            get;
            private set;
        }


        bool[] m_FreeMeanValue;


        /// <summary>
        /// pass-through from <see cref="ISpatialOperator.FreeMeanValue"/>
        /// </summary>
        public bool[] FreeMeanValue {
            get {
                if(m_FreeMeanValue == null) {
                    m_FreeMeanValue = FinerLevel.FreeMeanValue;
                }
                return m_FreeMeanValue;
            }
        }

        /// <summary>
        /// Recursive constructor, i.e. constructs operators on all multigrid levels which are provided by <paramref name="basisSeq"/>.
        /// </summary>
        /// <param name="basisSeq"></param>
        /// <param name="_ProblemMapping">
        /// DG coordinate mapping on the original grid, see <see cref="BaseGridProblemMapping"/>.
        /// </param>
        /// <param name="OperatorMatrix">
        /// Operator Matrix, aka. Jacobian Matrix,
        /// agglomeration (if applicable) should already be applied.
        /// </param>
        /// <param name="MassMatrix">
        /// Mass matrix on the original grid. If null, the identity matrix is assumed. 
        /// </param>
        /// <param name="cobc">
        /// Configuration of the cell-wise, explicit block-preconditioning for each multigrid level.
        /// (Remark: this kind of preconditioning is mathematically equivalent to a change of the DG resp. XDG basis.)
        /// </param>
        /// <param name="FreeMeanValue">
        /// pass-through from <see cref="ISpatialOperator.FreeMeanValue"/>
        /// </param>
        public MultigridOperator(IEnumerable<AggregationGridBasis[]> basisSeq,
            UnsetteledCoordinateMapping _ProblemMapping, BlockMsrMatrix OperatorMatrix, BlockMsrMatrix MassMatrix,
            IEnumerable<ChangeOfBasisConfig[]> cobc, bool[] FreeMeanValue)
            : this(null, basisSeq, _ProblemMapping, cobc) //
        {
            if (!OperatorMatrix.RowPartitioning.EqualsPartition(_ProblemMapping))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!OperatorMatrix.ColPartition.EqualsPartition(_ProblemMapping))
                throw new ArgumentException("Column partitioning mismatch.");

            if(FreeMeanValue == null) {
                FreeMeanValue = new bool[_ProblemMapping.BasisS.Count];
            } else {
                if(FreeMeanValue.Length != _ProblemMapping.BasisS.Count)
                    throw new ArgumentException();
            }
            m_FreeMeanValue = FreeMeanValue.CloneAs();

            if (MassMatrix != null) {
                if (!MassMatrix.RowPartitioning.Equals(_ProblemMapping))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!MassMatrix.ColPartition.Equals(_ProblemMapping))
                    throw new ArgumentException("Column partitioning mismatch.");
            }
            
            DefineReferenceIndices();

            if (this.LevelIndex == 0) {

                var bkup = SetPressureReferencePointMTX(OperatorMatrix);

                if (this.IndexIntoProblemMapping_Local == null) {
                    this.m_RawOperatorMatrix = OperatorMatrix;

                    if (MassMatrix != null)
                        this.m_RawMassMatrix = MassMatrix;  //is BlockMsrMatrix) ? ((BlockMsrMatrix)MassMatrix) : MassMatrix.ToMsrMatrix();
                    else
                        this.m_RawMassMatrix = null;
                } else {
                    this.m_RawOperatorMatrix = new BlockMsrMatrix(this.Mapping, this.Mapping);
                    
                    OperatorMatrix.WriteSubMatrixTo(this.m_RawOperatorMatrix, this.IndexIntoProblemMapping_Global, default(long[]), this.IndexIntoProblemMapping_Global, default(long[]));

                    if (MassMatrix != null) {
                        this.m_RawMassMatrix = new BlockMsrMatrix(this.Mapping, this.Mapping);
                        BlockMsrMatrix MMR = MassMatrix;
                        MMR.WriteSubMatrixTo(this.m_RawMassMatrix, this.IndexIntoProblemMapping_Global, default(long[]), this.IndexIntoProblemMapping_Global, default(long[]));
                    } else {
                        this.m_RawMassMatrix = null;
                    }
                }
            }
        }

        BlockMsrMatrix m_RawOperatorMatrix = null; // forgotten after Setup()
        BlockMsrMatrix m_RawMassMatrix = null;

        bool setupdone = false;

        void Setup() {
            using (var tr = new FuncTrace()) {
                if (setupdone)
                    return;
                setupdone = true;

                using(new BlockTrace("FinerLevel", tr)) {

                    bool bCheck = (this.FinerLevel != null);
                    bool bCheckGlob = bCheck.MPIOr();
                    if(bCheck != bCheckGlob)
                        throw new ApplicationException("Inconsistent data structure among MPI processes.");

                    if(this.FinerLevel != null)
                        this.FinerLevel.Setup();
                }

                // Construct intermediate 'raw' restriction and prolongation operators
                // ===================================================================

                BlockMsrMatrix RawRestriction, RawProlongation;
                using(new BlockTrace("ComputeProlongationAndRestriction", tr)) {
                    if(this.FinerLevel == null) {
                        RawRestriction = null;
                        RawProlongation = null;
                    } else {
                        //#if DEBUG
                        //                    var __PrlgOperator_Check = this.FinerLevel.Mapping.FromOtherLevelMatrix(this.Mapping);
                        //#endif
                        var __PrlgOperator = this.Mapping.GetProlongationOperator(this.FinerLevel.Mapping);
                        //var __PrlgOperator = this.FinerLevel.Mapping.FromOtherLevelMatrix(this.Mapping);
                        Debug.Assert(__PrlgOperator.RowPartitioning.LocalLength == this.FinerLevel.Mapping.LocalLength);
                        Debug.Assert(__PrlgOperator.ColPartition.LocalLength == this.Mapping.LocalLength);


                        if(this.FinerLevel.RightChangeOfBasis_Inverse != null)
                            RawProlongation = BlockMsrMatrix.Multiply(this.FinerLevel.RightChangeOfBasis_Inverse, __PrlgOperator);
                        else
                            RawProlongation = __PrlgOperator;

                        RawRestriction = RawProlongation.Transpose();
                    }
                    this.m_PrologateOperator = RawProlongation;
                    this.m_RestrictionOperator = RawRestriction;
                }

                // Construct intermediate 'raw' operator and mass matrix (before change of basis)
                // ==============================================================================

                // operator matrix before change of basis
                BlockMsrMatrix RawOpMatrix;
                using(new BlockTrace("RawOpMatrix", tr)) {
                    if(this.FinerLevel == null) {
                        RawOpMatrix = this.m_RawOperatorMatrix;
                    } else {
                        BlockMsrMatrix Op = FinerLevel.OperatorMatrix;
                        RawOpMatrix = BlockMsrMatrix.Multiply(RawRestriction, BlockMsrMatrix.Multiply(Op, RawProlongation));
                    }
                    this.m_RawOperatorMatrix = null;
                }

                // mass matrix before change of basis
                BlockMsrMatrix RawMassMatrix;
                using(new BlockTrace("RawMassMatrix", tr)) {
                    if(this.FinerLevel == null) {
                        RawMassMatrix = this.m_RawMassMatrix;
                    } else {
                        BlockMsrMatrix MM = FinerLevel.MassMatrix;
                        RawMassMatrix = BlockMsrMatrix.Multiply(RawRestriction, BlockMsrMatrix.Multiply(MM, RawProlongation));
                    }
                    this.m_RawMassMatrix = null;
                }
                
                if(RawOpMatrix.RowPartitioning.LocalLength != this.Mapping.LocalLength) {
                    throw new ApplicationException("internal error");
                }

                // compute change of basis
                // =======================
                long[] IndefRows;
                using(new BlockTrace("ComputeChangeOfBasis", tr)) {
                    IndefRows = this.ComputeChangeOfBasis(RawOpMatrix, RawMassMatrix, out m_LeftChangeOfBasis, out m_RightChangeOfBasis, out m_LeftChangeOfBasis_Inverse, out m_RightChangeOfBasis_Inverse);
                }


                // apply change of basis to operator matrix
                // ========================================
                { 
                    var Lpc = this.m_LeftChangeOfBasis;
                    var Rpc = this.m_RightChangeOfBasis;

                    BlockMsrMatrix O1;
                    using(new BlockTrace("ApplyChangeOfBasis_left", tr)) {
                        if(Lpc != null) {
                            O1 = BlockMsrMatrix.Multiply(Lpc, RawOpMatrix);
                        } else {
                            O1 = RawOpMatrix;
                        }
                    }

                    using(new BlockTrace("ApplyChangeOfBasis_right", tr)) {
                        if(Rpc != null) {
                            this.m_OperatorMatrix = BlockMsrMatrix.Multiply(O1, Rpc);
                        } else {
                            this.m_OperatorMatrix = O1;
                        }
                    }

                    // fix zero rows 
                    // (possible result from the Mode.IdMass_DropIndefinite or Mode.SymPart_DiagBlockEquilib_DropIndefinite -- option)

                    using(new BlockTrace("ApplyChangeOfBasis_indefRowFix", tr)) {
                        foreach(int _iRow in IndefRows) {
                            int iRow = Math.Abs(_iRow);
                            Debug.Assert(this.m_OperatorMatrix.GetNoOfNonZerosPerRow(iRow) == 0);
                            this.m_OperatorMatrix[iRow, iRow] = 1.0;
                        }
#if DEBUG
                        for(long iRow = this.m_OperatorMatrix.RowPartitioning.i0; iRow < this.m_OperatorMatrix.RowPartitioning.iE; iRow++) {
                            Debug.Assert(this.m_OperatorMatrix.GetNoOfNonZerosPerRow(iRow) > 0);
                        }
#endif
                    }
                }

                // apply change of basis to mass matrix
                // ====================================
                {
                    var Lpc = this.m_LeftChangeOfBasis;
                    var Rpc = this.m_RightChangeOfBasis;

                    if (RawMassMatrix != null) {

                        BlockMsrMatrix O1;
                        using(new BlockTrace("ApplyChangeOfBasisMass_left", tr)) {
                            if(Lpc != null) {
                                O1 = BlockMsrMatrix.Multiply(Lpc, RawMassMatrix);
                            } else {
                                O1 = RawMassMatrix;
                            }
                        }

                        using(new BlockTrace("ApplyChangeOfBasisMass_right", tr)) {
                            if(Rpc != null) {
                                this.m_MassMatrix = BlockMsrMatrix.Multiply(O1, Rpc);
                            } else {
                                this.m_MassMatrix = O1;
                            }
                        }
                    } else {
                        using(new BlockTrace("MassMatrixComp", tr)) {
                            if(Lpc != null && Rpc != null) {
                                this.m_MassMatrix = BlockMsrMatrix.Multiply(Lpc, Rpc);
                            } else if(Lpc != null) {
                                Debug.Assert(Rpc == null);
                                this.m_MassMatrix = Lpc;
                            } else if(Rpc != null) {
                                Debug.Assert(Lpc == null);
                                this.m_MassMatrix = Rpc;
                            } else {
                                Debug.Assert(Lpc == null);
                                Debug.Assert(Rpc == null);
                                this.m_MassMatrix = null;
                            }
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
        /// to this multi-grid level.
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
        long[] IndexIntoProblemMapping_Global;


        private MultigridOperator(MultigridOperator __FinerLevel, IEnumerable<AggregationGridBasis[]> basisES, UnsetteledCoordinateMapping _pm, IEnumerable<ChangeOfBasisConfig[]> cobc) {
            if (basisES.Count() <= 0) {
                throw new ArgumentException("At least one multigrid level is required.");
            }
            
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
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
                if (J != BaseGridProblemMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells)
                    throw new Exception("No of local cells wrong");
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
                    this.IndexIntoProblemMapping_Global = this.IndexIntoProblemMapping_Local.Select(i => (long)i).ToArray();
                } else {
                    this.IndexIntoProblemMapping_Global = this.IndexIntoProblemMapping_Local.Select(i => (long)i).ToArray();
                    long i0Proc = this.Mapping.ProblemMapping.i0;
                    int L = this.IndexIntoProblemMapping_Global.Length;
                    for (int i = 0; i < L; i++) {
                        this.IndexIntoProblemMapping_Global[i] += i0Proc;
                    }
                }
            }
#if DEBUG
            if (IndexIntoProblemMapping_Local != null) {
                long i0 = this.Mapping.ProblemMapping.i0;
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
        /// Pointer to operator on finest level.
        /// </summary>
        public MultigridOperator FinestLevel {
            get {
                if (FinerLevel != null)
                    return FinerLevel.FinestLevel;
                else
                    return this;
            }
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

        /// <summary>
        /// Number of Bytes used
        /// </summary>
        public long UsedMemory {
            get {
                GetMemoryInfo(out long Allocated, out long Used);
                return Used;
            }
        }

       

        /// <summary>
        /// Returns the amount of allocated/reserved and actually used memory in this level and coarser levels
        /// </summary>
        public void GetMemoryInfo(out long Allocated, out long Used) {
            Allocated = 0;
            Used = 0;
                       
            var allMtx = new HashSet<BlockMsrMatrix>(new FuncEqualityComparer<BlockMsrMatrix>((a, b) => ReferenceEquals(a, b)));
            allMtx.Add(m_LeftChangeOfBasis);
            allMtx.Add(m_LeftChangeOfBasis_Inverse);
            allMtx.Add(m_RightChangeOfBasis);
            allMtx.Add(m_RightChangeOfBasis_Inverse);
            allMtx.Add(m_RestrictionOperator);
            allMtx.Add(m_PrologateOperator);
            allMtx.Add(m_RawMassMatrix);
            allMtx.Add(m_RawOperatorMatrix);
            allMtx.Add(m_MassMatrix);
            allMtx.Add(m_OperatorMatrix);
                       


            foreach(var Mtx in allMtx) {
                if (Mtx != null) {
                    Mtx.GetMemoryInfo(out long allc, out long used);
                    Allocated += allc;
                    Used += used;
                }
            }


            if (CoarserLevel != null) {
                CoarserLevel.GetMemoryInfo(out long allc, out long used);
                Allocated += allc;
                Used += used;
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



        /// <summary>
        /// 
        /// </summary>
        /// <param name="solver"></param>
        /// <param name="INOUT_X">
        /// - on exit: (hopefully) the solution to the linear problem
        /// </param>
        /// <param name="IN_RHS">
        /// input: right-hand-side of the linear problem
        /// </param>
        /// <param name="UseGuess">
        /// - true: use <paramref name="INOUT_X"/> on entry as an initial guess
        /// </param>
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
            this.TransformRhsInto(IN_RHS, B, true);
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

        public void TransformRhsInto<T1, T2>(T1 u_IN, T2 v_OUT, bool ApplyRef)
            where T1 : IList<double>
            where T2 : IList<double> 
        {
            if(this.FinerLevel != null)
                throw new NotSupportedException("Only supported on finest level.");
            if(u_IN.Count != this.Mapping.ProblemMapping.LocalLength)
                throw new ArgumentException("Mismatch in length of input vector.", "u");
            if(v_OUT.Count != this.Mapping.LocalLength)
                throw new ArgumentException("Mismatch in length of output vector.", "v");


            (int idx, double val)[] bkup = null;
            if(ApplyRef)
                bkup = this.SetPressureReferencePointRHS(u_IN);

            int L = this.Mapping.LocalLength;
            double[] uc = new double[L];
            uc.AccV(1.0, u_IN, default(int[]), this.IndexIntoProblemMapping_Local);

            if(this.LeftChangeOfBasis != null) {
                this.LeftChangeOfBasis.SpMV(1.0, uc, 0.0, v_OUT);
            } else {
                v_OUT.SetV(uc);
            }

            if(bkup != null) {
                foreach(var t in bkup) {
                    u_IN[t.idx] = t.val;
                }
            }
        }

        /// <summary>
        /// transform from this operators domain to the original operator domain (for the trial space)
        /// </summary>
        /// <param name="u_Out">output; coordinate vector in the trial space of the original operator, <see cref="BaseGridProblemMapping"/></param>
        /// <param name="v_In">input: coordinate vector in the trial space of this operator</param>
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
        /// transform from this operators domain to the original operator domain (for the test space)
        /// </summary>
        /// <param name="u_Out">output; coordinate vector in the test space of the original operator, <see cref="BaseGridProblemMapping"/></param>
        /// <param name="v_In">input: coordinate vector in the test space of this operator</param>
        public void TransformRhsFrom<T1, T2>(T1 u_Out, T2 v_In)
            where T1 : IList<double>
            where T2 : IList<double>  //
        {
            if(this.FinerLevel != null)
                throw new NotSupportedException("Only supported on finest level.");
            if(u_Out.Count != this.Mapping.ProblemMapping.LocalLength)
                throw new ArgumentException("Mismatch in length of output vector.", "u");
            Debug.Assert(this.BaseGridProblemMapping.EqualsPartition(this.Mapping.ProblemMapping));
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
