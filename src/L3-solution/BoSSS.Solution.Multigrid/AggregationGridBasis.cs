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
using ilPSP;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using ilPSP.Utils;
using BoSSS.Foundation;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Multigrid {


 
    public class AggregationGridBasis {


        static GridData GetGridData(AggregationGrid ag) {
            if(ag.ParentGrid is GridData) {
                return ((GridData)ag.ParentGrid);
            } else if(ag.ParentGrid is AggregationGrid) {
                return GetGridData((AggregationGrid)ag.ParentGrid);
            } else {
                throw new NotSupportedException();
            }
        }

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="b">DG basis on original grid</param>
        /// <param name="ag"></param>
        public AggregationGridBasis(Basis b, AggregationGrid ag) {
            using (new FuncTrace()) {
                if (!object.ReferenceEquals(b.GridDat, GetGridData(ag)))
                    throw new ArgumentException("mismatch in grid data object.");
                this.DGBasis = b;
                this.AggGrid = ag;
                int N = b.Length;

                int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
                CompositeBasis = new MultidimensionalArray[JAGG];

                for (int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over agglomerated cells...
                    var compCell = ag.iLogicalCells.AggregateCellToParts[jAgg];


                    if (compCell.Length == 1) {
                        CompositeBasis[jAgg] = MultidimensionalArray.Create(1, N, N);
                        for (int n = 0; n < N; n++) {
                            CompositeBasis[jAgg][0, n, n] = 1.0;
                        }

                    } else {
                        // compute extrapolation basis
                        // ===========================

                        int I = compCell.Length - 1;
                        int[,] CellPairs = new int[I, 2];

                        for (int i = 0; i < I; i++) {
                            CellPairs[i, 0] = compCell[0];
                            CellPairs[i, 1] = compCell[i + 1];
                        }
                        var ExpolMtx = MultidimensionalArray.Create(I + 1, N, N);
                        b.GetExtrapolationMatrices(CellPairs, ExpolMtx.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { I, N - 1, N - 1 }));
                        for (int n = 0; n < N; n++) {
                            ExpolMtx[0, n, n] = 1.0;
                        }

                        // compute mass matrix
                        // ===================

                        var MassMatrix = MultidimensionalArray.Create(N, N);
                        MassMatrix.Multiply(1.0, ExpolMtx, ExpolMtx, 0.0, "lm", "kim", "kil");


                        // change to orthonormal basis
                        // ===========================
                        MultidimensionalArray B = MultidimensionalArray.Create(N, N);
                        MassMatrix.SymmetricLDLInversion(B, default(double[]));


                        CompositeBasis[jAgg] = MultidimensionalArray.Create(ExpolMtx.Lengths);
                        CompositeBasis[jAgg].Multiply(1.0, ExpolMtx, B, 0.0, "imn", "imk", "kn");

                        // check
                        // =====
#if DEBUG
                        MassMatrix.Clear();
                        for (int k = 0; k <= I; k++) {
                            for (int l = 0; l < N; l++) { // over rows of mass matrix ...
                                for (int m = 0; m < N; m++) { // over columns of mass matrix ...

                                    double mass_lm = 0.0;

                                    for (int i = 0; i < N; i++) {
                                        mass_lm += CompositeBasis[jAgg][k, i, m] * CompositeBasis[jAgg][k, i, l];
                                    }

                                    MassMatrix[l, m] += mass_lm;
                                }
                            }
                        }

                        MassMatrix.AccEye(-1.0);
                        Debug.Assert(MassMatrix.InfNorm() < 1.0e-8);
#endif
                    }
                }
            }
        }

        /// <summary>
        /// restricts/projects a vector from the full grid (<see cref="BoSSS.Foundation.Grid.Classic.GridData"/>)
        /// to the aggregated grid (<see cref="AggGrid"/>).
        /// </summary>
        /// <param name="FullGridVector">input;</param>
        /// <param name="AggGridVector">output;</param>
        virtual public void RestictFromFullGrid<T, V>(T FullGridVector, V AggGridVector)
            where T : IList<double>
            where V : IList<double> //
        {
            var fullMapping = new UnsetteledCoordinateMapping(this.DGBasis);
            if (FullGridVector.Count != fullMapping.LocalLength)
                throw new ArgumentException("mismatch in vector length", "FullGridVector");
            int L = this.LocalDim;
            if (AggGridVector.Count != L)
                throw new ArgumentException("mismatch in vector length", "AggGridVector");

            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int N = this.DGBasis.Length;

            var Buffer = MultidimensionalArray.Create(agCls.Max(cc => cc.Length), N);
            var Aggcoords = MultidimensionalArray.Create(N);

            for (int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over all aggregated cells...
                int[] agCl = agCls[jAgg];
                int K = agCl.Length;

                MultidimensionalArray coords = (K * N == Buffer.GetLength(0)) ? Buffer : Buffer.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, N - 1 });

                for (int k = 0; k < K; k++) { // loop over the cells which form the aggregated cell...
                    int jCell = agCl[k];
                    int j0 = fullMapping.LocalUniqueCoordinateIndex(0, jCell, 0);
                    for (int n = 0; n < N; n++) {
                        coords[k, n] = FullGridVector[n + j0];
                    }
                }

                Aggcoords.Clear();
                Aggcoords.Multiply(1.0, this.CompositeBasis[jAgg], coords, 0.0, "n", "kmn", "km");

                int i0 = jAgg * N;
                for (int n = 0; n < N; n++)
                    AggGridVector[n + i0] = Aggcoords[n];
            }
        }


        /// <summary>
        /// Prolongates/injects a vector from the full grid (<see cref="BoSSS.Foundation.Grid.GridData"/>)
        /// to the aggregated grid (<see cref="AggGrid"/>).
        /// </summary>
        /// <param name="FullGridVector">output;</param>
        /// <param name="AggGridVector">input;</param>
        virtual public void ProlongateToFullGrid<T, V>(T FullGridVector, V AggGridVector)
            where T : IList<double>
            where V : IList<double> //
        {
            var fullMapping = new UnsetteledCoordinateMapping(this.DGBasis);
            if(FullGridVector.Count != fullMapping.LocalLength)
                throw new ArgumentException("mismatch in vector length", "FullGridVector");
            int L = this.LocalDim;
            if(AggGridVector.Count != L)
                throw new ArgumentException("mismatch in vector length", "AggGridVector");


            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int N = this.DGBasis.Length;

            var FulCoords = new double[N];
            var AggCoords = new double[N];

            //MultidimensionalArray Trf = MultidimensionalArray.Create(N, N);

            for(int jAgg = 0; jAgg < JAGG; jAgg++) {
                int[] agCl = agCls[jAgg];
                int K = agCl.Length;

                int i0 = jAgg * N;
                for(int n = 0; n < N; n++)
                    AggCoords[n] = AggGridVector[n + i0];

                for(int k = 0; k < K; k++) { // loop over the cells wich form the aggregated cell...
                    int jCell = agCl[k];
                    int j0 = fullMapping.LocalUniqueCoordinateIndex(0, jCell, 0);

                    //Trf.Clear();
                    //Trf.Acc(1.0, this.CompositeBasis[jAgg].ExtractSubArrayShallow(k, -1, -1));
                    var Trf = this.CompositeBasis[jAgg].ExtractSubArrayShallow(k, -1, -1);

                    //Trf.Solve(FulCoords, AggCoords);
                    Trf.gemv(1.0, AggCoords, 0.0, FulCoords);

                    for(int n = 0; n < N; n++) {
                        FullGridVector[j0 + n] = FulCoords[n];
                    }
                }
            }
        }


        virtual public void GetRestrictionMatrix(BlockMsrMatrix rest, MultigridMapping mgMap, int iFld) {
            if(!object.ReferenceEquals(mgMap.AggBasis[iFld], this))
                throw new ArgumentException();

            UnsetteledCoordinateMapping fullmap = mgMap.ProblemMapping;
            int[] degree = mgMap.DgDegree;

            foreach(var b in fullmap.BasisS) {
                if(!b.IsSubBasis(this.DGBasis))
                    throw new ArgumentException("");
                if(b.MaximalLength != b.MinimalLength)
                    throw new NotSupportedException();
            }
            if(fullmap.BasisS.Count != degree.Length)
                throw new ArgumentException();
            int NoFld = degree.Length;

            int[] FullLength = fullmap.BasisS.Select(b => b.Length).ToArray();
            int[] RestLength = NoFld.ForLoop(
                l => fullmap.BasisS[l].Polynomials[0].Where(poly => poly.AbsoluteDegree <= degree[l]).Count());
            int[] RestOffset = new int[NoFld];
            for(int f = 1; f < NoFld; f++) {
                RestOffset[f] = RestOffset[f - 1] + RestLength[f - 1];
            }

            int NR = RestLength.Sum();
            int JAGG = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;


            //Partitioning RowMap = new Partitioning(JAGG * NR);
            int MpiOffset_row = rest.RowPartitioning.i0;
            
            for(int jagg = 0; jagg < JAGG; jagg++) {
                int[] agCl = agCls[jagg];
                int K = agCl.Length;

                for (int k = 0; k < K; k++) { // loop over the cells which form the aggregated cell...
                    int jCell = agCl[k];

                    int M = RestLength[iFld];
                    int N = FullLength[iFld];
                    var Block = this.CompositeBasis[jagg].ExtractSubArrayShallow(new int[] { k, 0, 0 }, new int[] { k - 1, N - 1, M - 1 });

                    //for(int f = 0; f < NoFld; f++) { // loop over DG fields in mapping...

                    int i0Col = fullmap.GlobalUniqueCoordinateIndex(iFld, jCell, 0);
                    int i0Row = jagg * NR + RestOffset[iFld] + MpiOffset_row;

                    //rest.AccBlock(i0Row, i0Col, 1.0, Block);
                    for (int m = 0; m < M; m++) {
                        for (int n = 0; n < N; n++) {
                            rest[i0Row + m, i0Col + n] = Block[n, m];
                        }
                    }

                }
            }
        }

        
        public AggregationGrid AggGrid {
            get;
            private set;
        }

        public Basis DGBasis {
            get;
            private set;
        }

        public virtual int MaximalLength {
            get {
                return this.DGBasis.MaximalLength;
            }
        }

        public virtual int MinimalLength {
            get {
                return this.DGBasis.MinimalLength;
            }
        }


        int[] m_Lengths;

        public virtual int GetMaximalLength(int p) {
            Debug.Assert(this.DGBasis.MaximalLength == this.DGBasis.MinimalLength);
            return this.GetLength(0, p);
        }

        public virtual int GetMinimalLength(int p) {
            Debug.Assert(this.DGBasis.MaximalLength == this.DGBasis.MinimalLength);
            return this.GetLength(0, p);
        }
        
        public virtual int GetLength(int jCell, int p) {
            if(m_Lengths == null) {
                m_Lengths = new int[this.DGBasis.Degree + 1];
                for(int pp = 0; pp < m_Lengths.Length; pp++) {
                    m_Lengths[pp] = this.DGBasis.Polynomials[0].Where(poly => poly.AbsoluteDegree <= pp).Count();
                }
            }

            return m_Lengths[p];
        }

        /// <summary>
        /// Local vector-space dimension.
        /// </summary>
        virtual public int LocalDim {
            get {
                return this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells * this.DGBasis.Length;
            }
        }

        /// <summary>
        /// The projector in the L2 Norm, from the 
        /// space defined by the basis <see cref="DGBasis"/>, onto the 
        /// the DG space on the aggregate grid.
        /// - array index: aggregate cell index; 
        /// - 1st index:  
        /// - 2nd index: 
        /// - 3rd index:
        /// </summary>
        public MultidimensionalArray[] CompositeBasis;


        /// <summary>
        /// for XDG, the cell mode index <paramref name="n"/> may not be equal
        /// in the full and the aggregated grid. This method performs the transformation.
        /// </summary>
        virtual internal int N_Murks(int j, int n, int N) {
            Debug.Assert(j >= 0 && j < this.DGBasis.GridDat.iLogicalCells.NoOfCells);
            Debug.Assert(n >= 0 && n < this.DGBasis.GetLength(j));
            Debug.Assert(n < N);
            Debug.Assert(this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells == this.DGBasis.GridDat.iLogicalCells.NoOfCells);
            return n;
        }

        virtual internal bool ReqModeIndexTrafo {
            get {
                return false;
            }
        }

        int[][] m_ModeIndexForDegree;

        virtual public int[] ModeIndexForDegree(int j, int p, int Pmax) {
            if(m_ModeIndexForDegree == null) {
                int Psupermax = this.DGBasis.Degree;
                m_ModeIndexForDegree = new int[Psupermax + 1][];
                for(int _p = 0; _p <= Psupermax; _p++) {
                    m_ModeIndexForDegree[_p] = new int[0];
                }

                Debug.Assert(this.DGBasis.Polynomials.Count == 1);
                var polys = this.DGBasis.Polynomials[0];

                for(int n = 0; n < polys.Count; n++) {
                    n.AddToArray(ref m_ModeIndexForDegree[polys[n].AbsoluteDegree]);
                }
            }
            return m_ModeIndexForDegree[p];
        }
        
    }

}
