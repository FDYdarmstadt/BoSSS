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
using BoSSS.Foundation;
using BoSSS.Platform;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.XNSECommon.Operator.SurfaceTension {


    public enum PatchRecoveryMode {
        L2_restrictedDom,

        L2_unrestrictedDom,

        ChebychevInteroplation
    }


    public class L2PatchRecovery {

        /// <summary>
        /// the domain on which the patch recovery is performed
        /// </summary>
        public CellMask Domain {
            get;
            private set;
        }

        public L2PatchRecovery(Basis __bInput, Basis __bOutput, CellMask cm, bool RestrictToCellMask = true) {
            this.m_bInput = __bInput;
            this.m_bOutput = __bOutput;
            var GridDat = m_bInput.GridDat;
            
            this.AggregateBasisTrafo = new MultidimensionalArray[GridDat.iLogicalCells.NoOfCells];
            this.UpdateDomain(cm, RestrictToCellMask);
        }


        public void UpdateDomain(CellMask cm, bool RestrictToCellMask = true) {
            var GridDat = m_bInput.GridDat;
            int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            this.Domain = cm;

            int[][] newStencils = new int[J][];

            var Mask = cm.GetBitMaskWithExternal();
            foreach(int jCell in cm.ItemEnum) {
                int[] NeighCells, dummy;
                GridDat.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out NeighCells, out dummy);

                if(RestrictToCellMask == true) {
                    NeighCells = NeighCells.Where(j => Mask[j]).ToArray();
                }
                Array.Sort(NeighCells);
                int[] newStencil = ArrayTools.Cat(new int[] { jCell }, NeighCells); 

                if(this.Stencils == null
                    || this.Stencils[jCell] == null && newStencil.Length > 0
                    || !ArrayTools.AreEqual(this.Stencils[jCell], newStencil)) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // a re-computation of the aggregate cell basis is necessary
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    ComputeAggregateBasis(jCell, NeighCells);
                }
                newStencils[jCell] = newStencil;                
            }

            for(int j = 0; j < J; j++) {
                if(newStencils[j] == null && this.AggregateBasisTrafo[j] != null)
                    this.AggregateBasisTrafo[j] = null;

                Debug.Assert((newStencils[j] == null) == (this.AggregateBasisTrafo[j] == null));
                if(newStencils[j] != null)
                    Debug.Assert(newStencils[j].Length == this.AggregateBasisTrafo[j].GetLength(0));
            }


            this.Stencils = newStencils;
        }




        public void SetLimiter(int[] limDeg) {
            if(limDeg == null) {
                this.Nlim = null;
            } else {

                int J = this.Domain.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                if(limDeg.Length != J)
                    throw new ArgumentException();

                var _Nlim = new int[J];

                Basis B = this.m_bInput.Degree > this.m_bOutput.Degree ? this.m_bInput : this.m_bOutput;
                int[] DegToN = (B.Degree + 1).ForLoop(deg => B.Polynomials[0].Where(p => p.AbsoluteDegree <= deg).Count());

                int degMax = B.Degree;
                for(int j = 0; j < J; j++) {
                    if(limDeg[j] > degMax)
                        throw new ArgumentOutOfRangeException();

                    if(limDeg[j] >= 0) {
                        _Nlim[j] = DegToN[limDeg[j]];
                        Debug.Assert(_Nlim[j] > 0);
                    }
                }
                
                this.Nlim = _Nlim;
            }
        }

        public bool notchangeunlim = false;


        private int[] Nlim = null;

        void ComputeAggregateBasis(int jCell, int[] Stencil_jCells) {
            // implementation notes: Fk, persönliche Notizen, 17oct13
            // ------------------------------------------------------
           
            int N = m_bOutput.Length;          // Basis dimension of 'g' in cell 'jCell'
            int K = Stencil_jCells.Length + 1; // number of cells in stencil (including jCell itself);

            var ExPolMtx = MultidimensionalArray.Create(K, N, N);

            int[,] CellPairs = new int[K - 1, 2];
            for(int k = 0; k < (K - 1); k++) {
                CellPairs[k, 0] = jCell;
                CellPairs[k, 1] = Stencil_jCells[k];
            }
            if(K > 1)
                m_bOutput.GetExtrapolationMatrices(CellPairs, ExPolMtx.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { K - 1, N - 1, N - 1 }), null);
            for(int l = 0; l < N; l++) {
                ExPolMtx[0, l, l] = 1.0;
            }


            MultidimensionalArray MassMatrix = MultidimensionalArray.Create(N, N);
            //MassMatrix.AccEye(1.0); // Mass matrix in jCell itself

            //for(int k = 0; k < K; k++) { // over stencil members ...
            //    for(int l = 0; l < N; l++) { // over rows of mass matrix ...
            //        for(int m = 0; m < N; m++) { // over columns of mass matrix ...

            //            double mass_lm = 0.0;

            //            for(int i = 0; i < N; i++) {
            //                mass_lm += ExPolMtx[k, i, m] * ExPolMtx[k, i, l];
            //            }

            //            MassMatrix[l, m] += mass_lm;
            //        }
            //    }
            //}
            MassMatrix.Multiply(1.0, ExPolMtx, ExPolMtx, 0.0, "lm", "kim", "kil");

            MultidimensionalArray B = MultidimensionalArray.Create(N, N);
            MassMatrix.SymmetricLDLInversion(B, default(double[]));
            

            var CompositeBasis = MultidimensionalArray.Create(ExPolMtx.Lengths);
            CompositeBasis.Multiply(1.0, ExPolMtx, B, 0.0, "imn", "imk", "kn");

            
            //MassMatrix.InvertSymmetrical();
            //this.InvMassMatrix[jCell] = MassMatrix;

            this.AggregateBasisTrafo[jCell] = CompositeBasis;
        }

        int[][] Stencils;

        /// <summary>
        /// Basis transformation between the original DG basis and the aggregate DG basis.
        /// Note that the aggregate Basis is also orthonormal.
        /// </summary>
        MultidimensionalArray[] AggregateBasisTrafo;
        
        Basis m_bInput;
        Basis m_bOutput;

        public double[] diagnosis;


        /// <summary>
        /// %
        /// </summary>
        /// <param name="output"></param>
        /// <param name="input">input DG field; unchanged on </param>
        public void Perform(SinglePhaseField output, SinglePhaseField input) {
            if(!output.Basis.Equals(this.m_bOutput))
                throw new ArgumentException("output basis mismatch");
            if(!input.Basis.Equals(this.m_bInput))
                throw new ArgumentException("output basis mismatch");

            var GridDat = this.m_bOutput.GridDat;
            int N = m_bOutput.Length;
            int Nin = m_bInput.Length;
            int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            diagnosis = new double[N * J];
            
            double[] RHS = new double[N];
            double[] f2 = new double[N];
            double[] g1 = new double[N];
            MultidimensionalArray Trf = MultidimensionalArray.Create(N, N);

            input.MPIExchange();

                    
            for(int jCell = 0; jCell < J; jCell++) {
                Debug.Assert((this.AggregateBasisTrafo[jCell] != null) == (this.Stencils[jCell] != null));

                //FullMatrix invMassM = this.InvMassMatrix[jCell];
                MultidimensionalArray ExPolMtx = this.AggregateBasisTrafo[jCell];
                
                if(ExPolMtx != null) {
                    int[] Stencil_jCells = this.Stencils[jCell];
                    int K = Stencil_jCells.Length;
                    Debug.Assert(Stencil_jCells[0] == jCell);
                    
                    var coordIn = input.Coordinates;
                    //for(int n = Math.Min(Nin, N) - 1; n >= 0; n--) {
                    //    RHS[n] = coordIn[jCell, n];
                    //}
                    RHS.ClearEntries();
                    g1.ClearEntries();

                    var oldCoords = input.Coordinates.GetRow(jCell);

                    for(int k = 0; k < K; k++) {
                        int jNeigh = Stencil_jCells[k];
                        for(int n = Math.Min(input.Basis.Length, N) - 1; n >= 0; n--) {
                            f2[n] = input.Coordinates[jNeigh, n];
                        }
                        if(Nlim != null && Nlim[jNeigh] > 0) {
                            for(int n = Nlim[jNeigh]; n < RHS.Length; n++) {
                                f2[n] = 0;
                            }
                        }

                        for(int l = 0; l < N; l++) {
                            double acc = 0;
                            for(int i = 0; i < N; i++) {
                                acc += ExPolMtx[k, i, l] * f2[i];
                            }
                            RHS[l] += acc;
                        }
                    }

                    if(Nlim != null && Nlim[jCell] > 0) {
                        for(int n = Nlim[jCell]; n < RHS.Length; n++) {
                            RHS[n] = 0;
                        }
                    }

                    for(int l = 0; l < N; l++) {
                        diagnosis[jCell * N + l] = RHS[l];
                    }

                    //invMassM.gemv(1.0, RHS, 0.0, g1);

                    Trf.Clear();
                    Trf.Acc(1.0, ExPolMtx.ExtractSubArrayShallow(0, -1, -1));

                    //Trf.Solve(FulCoords, AggCoords);
                    Trf.gemv(1.0, RHS, 0.0, g1);


                    if(Nlim != null && Nlim[jCell] > 0) {
                        for(int n = Nlim[jCell]; n < RHS.Length; n++) {
                            g1[n] = 0;
                        }
                    }

                    if(Nlim == null) {
                        output.Coordinates.SetRow(jCell, g1);
                    } else {
                        if((Nlim[jCell] <= 0 && this.notchangeunlim))
                            output.Coordinates.SetRow(jCell, oldCoords);
                        else
                            output.Coordinates.SetRow(jCell, g1);
                    }

                }
            }

            output.MPIExchange();
        }

    }


    
    public class PatchRecovery {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="output">output</param>
        /// <param name="cm">mask on which the filtering should be performed</param>
        /// <param name="input">input</param>
        /// <param name="mode"></param>
        public static void FilterStencilProjection(SinglePhaseField output, CellMask cm, SinglePhaseField input, PatchRecoveryMode mode) {
            var GridDat = input.GridDat;
            if (!object.ReferenceEquals(output.GridDat, GridDat))
                throw new ArgumentException();
            if (!object.ReferenceEquals(input.GridDat, GridDat))
                throw new ArgumentException();

            if (cm == null)
                cm = CellMask.GetFullMask(GridDat);
            
            switch(mode) {
                case PatchRecoveryMode.L2_unrestrictedDom:
                case PatchRecoveryMode.L2_restrictedDom: {

                    var Mask = cm.GetBitMaskWithExternal();
                    bool RestrictToCellMask = mode == PatchRecoveryMode.L2_restrictedDom;
                    foreach (int jCell in cm.ItemEnum) {
                        int[] NeighCells, dummy;
                        GridDat.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out NeighCells, out dummy);

                        if (RestrictToCellMask == true) {
                            NeighCells = NeighCells.Where(j => Mask[j]).ToArray();
                        }

                        FilterStencilProjection(output, jCell, NeighCells, input);
                    }

                    return;
                }

                case PatchRecoveryMode.ChebychevInteroplation: {
                    int P = input.Basis.Degree*3;
            
                    double[] ChebyNodes = ChebyshevNodes(P);
                    NodeSet Nds = new NodeSet(GridDat.iGeomCells.RefElements[0], P, 1);
                    Nds.SetColumn(0, ChebyNodes);
                    Nds.LockForever();
                    Polynomial[] Polynomials = GetNodalPolynomials(ChebyNodes);

                    MultidimensionalArray PolyVals = MultidimensionalArray.Create(P, P);


                    for (int p = 0; p < P; p++) {
                        Polynomials[p].Evaluate(PolyVals.ExtractSubArrayShallow(p, -1), Nds);

                        for (int i = 0; i < P; i++)
                            Debug.Assert(Math.Abs(((i == p) ? 1.0 : 0.0) - PolyVals[p, i]) <= 1.0e-8);
                    }


                    foreach (int jCell in cm.ItemEnum) {
                        AffineTrafo T;
                        var NodalValues = NodalPatchRecovery(jCell, ChebyNodes, input, out T);

                        /*
                        MultidimensionalArray Nodes2D = MultidimensionalArray.Create(P, P, 2);
                        for (int i = 0; i < P; i++) {
                            for (int j = 0; j < P; j++) {
                                Nodes2D[i, j, 0] = ChebyNodes[i];
                                Nodes2D[i, j, 1] = ChebyNodes[j];
                            }
                        }
                        var _Nodes2D = Nodes2D.ResizeShallow(P*P, 2);
                        var xNodes = _Nodes2D.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { P*P - 1, 0 });
                        var yNodes = _Nodes2D.ExtractSubArrayShallow(new int[] { 0, 1 }, new int[] { P*P - 1, 1 });
                        int K = P*P;

                        MultidimensionalArray xPolyVals = MultidimensionalArray.Create(P, K);
                        MultidimensionalArray yPolyVals = MultidimensionalArray.Create(P, K);

                        for (int p = 0; p < P; p++) {
                            Polynomials[p].Evaluate(yPolyVals.ExtractSubArrayShallow(p, -1), yNodes);
                            Polynomials[p].Evaluate(xPolyVals.ExtractSubArrayShallow(p, -1), xNodes);

                        }


                        double ERR = 0;
                        for (int k = 0; k < K; k++) {

                            double x = xNodes[k, 0];
                            double y = yNodes[k, 0];

                            double acc = 0;
                            double BasisHits = 0.0;
                            for (int i = 0; i < P; i++) {
                                for (int j = 0; j < P; j++) {
                            
                                    bool isNode = (Math.Abs(x - ChebyNodes[i]) < 1.0e-8) && (Math.Abs(y - ChebyNodes[j]) < 1.0e-8);
                                    double BasisSoll = isNode ? 1.0 : 0.0;
                                    if (isNode)
                                        Console.WriteLine();

                                    double Basis_ij = xPolyVals[i, k]*yPolyVals[j, k];

                                    Debug.Assert((Basis_ij - BasisSoll).Abs() < 1.0e-8);

                                    BasisHits += Math.Abs(Basis_ij);

                                    acc += Basis_ij*NodalValues[i, j];
                                }
                            }
                            Debug.Assert(Math.Abs(BasisHits - 1.0) < 1.0e-8);

                    


                            double soll = yNodes[k,0];
                            Debug.Assert((soll - acc).Pow2() < 1.0e-7);
                            ERR = (soll - acc).Pow2();
                        }
                        Console.WriteLine(ERR);
                        */

                        output.ProjectField(1.0,
                            delegate(int j0, int Len, NodeSet Nodes, MultidimensionalArray result) {
                                var K = Nodes.NoOfNodes;

                                MultidimensionalArray NT = T.Transform(Nodes);
                                
                                var xNodes = new NodeSet(null, NT.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, 0 }));
                                var yNodes = new NodeSet(null, NT.ExtractSubArrayShallow(new int[] { 0, 1 }, new int[] { K - 1, 1 }));

                                MultidimensionalArray xPolyVals = MultidimensionalArray.Create(P, K);
                                MultidimensionalArray yPolyVals = MultidimensionalArray.Create(P, K);

                                for (int p = 0; p < P; p++) {
                                    Polynomials[p].Evaluate(yPolyVals.ExtractSubArrayShallow(p, -1), yNodes);
                                    Polynomials[p].Evaluate(xPolyVals.ExtractSubArrayShallow(p, -1), xNodes);

                                }


                                for (int k = 0; k < K; k++) {

                                    double acc = 0;
                                    for (int i = 0; i < P; i++) {
                                        for (int j = 0; j < P; j++) {
                                            acc += xPolyVals[i, k]*yPolyVals[j, k]*NodalValues[i, j];
                                        }
                                    }
                                    result[0, k] = acc;
                                }
                            },
                            new CellQuadratureScheme(true, new CellMask(input.GridDat, Chunk.GetSingleElementChunk(jCell))));

                    }


                    return;
                }
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="Stencil_jCell"></param>
        /// <param name="input">input</param>
        /// <param name="output">result</param>
        static void FilterStencilProjection(SinglePhaseField output, int jCell, int[] Stencil_jCells, ConventionalDGField input) {
            // implementation notes: Fk, persönliche Notizen, 17oct13
            // ------------------------------------------------------
            if (output.Basis.Degree < output.Basis.Degree)
                throw new ArgumentException();

            int N = output.Basis.Length;        // Basis dimension of 'g' in cell 'jCell'
            int K = Stencil_jCells.Length; // number of cells in stencil

            var ExPolMtx = MultidimensionalArray.Create(K, N, N);

            int[,] CellPairs = new int[K, 2];
            for (int k = 0; k < K; k++) {
                CellPairs[k, 0] = jCell;
                CellPairs[k, 1] = Stencil_jCells[k];
            }
            output.Basis.GetExtrapolationMatrices(CellPairs, ExPolMtx, null);

            MultidimensionalArray MassMatrix = MultidimensionalArray.Create(N, N);
            MassMatrix.AccEye(1.0); // Mass matrix in jCell itself

            //for (int k = 0; k < K; k++) { // over stencil members ...
            //    for (int l = 0; l < N; l++) { // over rows of mass matrix ...
            //        for (int m = 0; m < N; m++) { // over columns of mass matrix ...

            //            double mass_lm = 0.0;

            //            for (int i = 0; i < N; i++) {
            //                mass_lm += ExPolMtx[k, i, m]*ExPolMtx[k, i, l];
            //            }

            //            MassMatrix[l, m] += mass_lm;
            //        }
            //    }
            //}
            MassMatrix.Multiply(1.0, ExPolMtx, ExPolMtx, 1.0, "lm", "kim", "kil");


            double[] RHS = new double[N];
            double[] f2 = new double[N];
            for (int n = Math.Min(input.Basis.Length, N) - 1; n >= 0; n--) {
                RHS[n] = input.Coordinates[jCell, n];
            }
            for (int k = 0; k < K; k++) {
                for (int n = Math.Min(input.Basis.Length, N) - 1; n >= 0; n--) {
                    f2[n] = input.Coordinates[Stencil_jCells[k], n];
                }
                for (int l = 0; l < N; l++) {
                    double acc = 0;
                    for (int i = 0; i < N; i++) {
                        acc += ExPolMtx[k, i, l]*f2[i];
                    }
                    RHS[l] += acc;
                }
            }

            double[] g1 = new double[N];
            MassMatrix.Solve(g1, RHS);

            output.Coordinates.SetRow(jCell, g1);
        }


        static double[] ChebyshevNodes(int N) {
            double[] Ret = new double[N];
            for (int i = 0; i < N; i++) {
                //Ret[i] = Math.Cos(Math.PI*(2.0*(i + 1.0) - 1.0)/(2.0*((double)N)));
                Ret[N - i - 1] = Math.Cos((2.0*(i+1.0) - 1.0)/(2.0*N)*Math.PI);
            }
            return Ret;
        }



        static MultidimensionalArray NodalPatchRecovery(int jCell, double[] Nodes, DGField Phi, out AffineTrafo Trafilein) {
            GridData GridData = (GridData)Phi.GridDat;
            
            int[] Neighbours, Edges;
            GridData.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out Neighbours, out Edges);
            if (Neighbours.Length != 8)
                throw new NotSupportedException();

            int[] JS = Neighbours;
            jCell.AddToArray(ref JS);

            MultidimensionalArray CellCoordinates = MultidimensionalArray.Create(JS.Length, 4, 2); // 1st index: jcell, top, bottom, left, right
            var Kref = GridData.iGeomCells.RefElements[0];
            if (GridData.iGeomCells.RefElements.Length != 1 || Kref.GetType() != typeof(BoSSS.Foundation.Grid.RefElements.Square)) {
                throw new NotImplementedException();
            }
            for (int f = 0; f < JS.Length; f++) {
                GridData.TransformLocal2Global(Kref.Vertices, CellCoordinates.ExtractSubArrayShallow(f, -1, -1), JS[f]);
            }

            double xMin = CellCoordinates.ExtractSubArrayShallow(-1, -1, 0).Min();
            double xMax = CellCoordinates.ExtractSubArrayShallow(-1, -1, 0).Max();
            double yMin = CellCoordinates.ExtractSubArrayShallow(-1, -1, 1).Min();
            double yMax = CellCoordinates.ExtractSubArrayShallow(-1, -1, 1).Max();

            int N = Nodes.Length;
            double[] xNodes = new double[N];
            double[] yNodes = new double[N];
            for (int i = 0; i < N; i++) {
                xNodes[i] = (xMax - xMin)*0.5*(Nodes[i] + 1) + xMin;
                yNodes[i] = (yMax - yMin)*0.5*(Nodes[i] + 1) + yMin;
            }
            var TrChey2Glob = new AffineTrafo(2);
            TrChey2Glob.Matrix[0, 0] = (xMax - xMin)*0.5; TrChey2Glob.Affine[0] = (xMax - xMin)*0.5 + xMin;
            TrChey2Glob.Matrix[1, 1] = (yMax - yMin)*0.5; TrChey2Glob.Affine[1] = (yMax - yMin)*0.5 + yMin;


            MultidimensionalArray[] NodeSets_global = new MultidimensionalArray[JS.Length];
            NodeSet[] NodeSets = new NodeSet[JS.Length];
            List<double> xNodesCell = new List<double>();
            List<double> yNodesCell = new List<double>();
            int[][] MapTo_i = new int[JS.Length][];
            int[][] MapTo_j = new int[JS.Length][];


            for (int f = 0; f < JS.Length; f++) {
                MultidimensionalArray CellCoords = CellCoordinates.ExtractSubArrayShallow(f, -1, -1);
                double xMinCell = CellCoords.ExtractSubArrayShallow(-1, 0).Min();
                double xMaxCell = CellCoords.ExtractSubArrayShallow(-1, 0).Max();
                double yMinCell = CellCoords.ExtractSubArrayShallow(-1, 1).Min();
                double yMaxCell = CellCoords.ExtractSubArrayShallow(-1, 1).Max();

                int iMin = xNodes.FirstIndexWhere(x => x >= xMinCell);
                int iMax = xNodes.LastIndexWhere(x => x <  xMaxCell);
                int jMin = yNodes.FirstIndexWhere(y => y >= yMinCell);
                int jMax = yNodes.LastIndexWhere(y => y <  yMaxCell);
                int NxCell = iMax - iMin + 1;
                int NyCell = jMax - jMin + 1;

                int K = NxCell*NyCell;
                NodeSets_global[f] = MultidimensionalArray.Create(K, 2);

                MapTo_i[f] = new int[K];
                MapTo_j[f] = new int[K];

                int k = 0;
                for (int i = iMin; i <= iMax; i++) {
                    for (int j = jMin; j <= jMax; j++) {
                        MapTo_i[f][k] = i;
                        MapTo_j[f][k] = j;
                        NodeSets_global[f][k, 0] = xNodes[i];
                        NodeSets_global[f][k, 1] = yNodes[j];

                        Debug.Assert(GridData.Cells.IsInCell(new double[] { xNodes[i], yNodes[j] }, JS[f], null));

                        k++;
                    }
                }

                NodeSets[f] = new NodeSet(GridData.Cells.GetRefElement(jCell), NxCell*NyCell, 2);
                GridData.TransformGlobal2Local(NodeSets_global[f], NodeSets[f], JS[f], null);
                NodeSets[f].LockForever();
            }

            var TrGl2Loc = AffineTrafo.FromPoints(NodeSets_global.Last().ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { 2, 1 }), NodeSets.Last().ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { 2, 1 }));
            Trafilein = (TrGl2Loc*TrChey2Glob).Invert();



            var Return = MultidimensionalArray.Create(N, N);
            for (int f = 0; f < JS.Length; f++) {
                
                int K = NodeSets[f].GetLength(0);
                var Result = MultidimensionalArray.Create(1, K);
                Phi.Evaluate(JS[f], 1, NodeSets[f], Result);

                for (int k = 0; k < K; k++) {

                    //double x_i = xNodes[MapTo_i[f][k]];
                    //double y_j = yNodes[MapTo_j[f][k]];
                    //double val = (x_i*0.3).Pow2()  - (y_j*0.3).Pow2();



                    Debug.Assert(Return[MapTo_i[f][k], MapTo_j[f][k]] == 0.0);
                    Return[MapTo_i[f][k], MapTo_j[f][k]] = Result[0, k];

                }
            }

            // ----------------------------------

            return Return;
        }

        static Polynomial[] GetNodalPolynomials(double[] Nodes) {

            int N = Nodes.Length;
            Polynomial[] Ret = new Polynomial[N];

            MultidimensionalArray Mtx = MultidimensionalArray.Create(N, N);
            for (int n = 0; n < N; n++) {
                double x_n = Nodes[n];
                double a = 1.0;
                for (int i = 0; i < N; i++) {
                    Mtx[n, i] = a;
                    a *= x_n;
                }
            }


            for (int n = 0; n < N; n++) {
                double[] RHS = new double[N];
                RHS[n] = 1.0;

                double[] Sol = new double[N];
                Mtx.CloneAs().Solve(Sol, RHS);

                Polynomial poly = new Polynomial();
                for (int i = 0; i < N; i++) {
                    poly.AddCoeff(Sol[i], new int[] { i });
                }

                Ret[n] = poly;
            }


            return Ret;
        }


    }
}
