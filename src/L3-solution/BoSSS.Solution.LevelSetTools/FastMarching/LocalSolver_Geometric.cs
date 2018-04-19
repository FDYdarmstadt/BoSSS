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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {


    /// <summary>
    /// Solver for the Reinit-problem on one cell, based on a gemetric approach.
    /// </summary>
    class LocalSolver_Geometric {

        /// <summary>
        /// Quadrature rules used for the local geometric solver; 
        /// index: reference element 
        /// </summary>
        QuadRule[] DaRuleS;

        NodeSet EdgeNodes;

        Basis LevelSetBasis_Geometric;
        GridData GridDat;

        // Buffers for geometric solver
        MultidimensionalArray[] PhiEvalBuffer;
        MultidimensionalArray[] CellNodesGlobal;

        /// <summary>
        /// array index: reference element
        /// </summary>
        MultidimensionalArray[] QuadNodesGlobal;

        public LocalSolver_Geometric(Basis __LevelSetBasis_Geometric) {
            LevelSetBasis_Geometric = __LevelSetBasis_Geometric;
            GridDat = (GridData)(LevelSetBasis_Geometric.GridDat);

            int D = this.GridDat.SpatialDimension;
            if(D != 2)
                throw new NotImplementedException("Geometric solver currently only implemented for 2D.");
            
            if(!(this.GridDat.Edges.EdgeRefElements.Count() == 1 && this.GridDat.Edges.EdgeRefElements.First().GetType() == typeof(Foundation.Grid.RefElements.Line)))
                throw new NotImplementedException();


            double[] _EdgeNodes = GenericBlas.Linspace(-1, 1, Math.Max(10, (this.LevelSetBasis_Geometric.Degree + 1) * 2));
            int NN = _EdgeNodes.Length;
            this.EdgeNodes = new NodeSet(GridDat.Edges.EdgeRefElements[0], NN, 1);
            EdgeNodes.ExtractSubArrayShallow(-1, 0).SetVector(_EdgeNodes);
            this.EdgeNodes.LockForever();

            DaRuleS = this.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(this.LevelSetBasis_Geometric.Degree * 2)).ToArray();
            QuadNodesGlobal = DaRuleS.Select(rule => MultidimensionalArray.Create(rule.NoOfNodes, D)).ToArray();

            int NNKmax = 4;
            PhiEvalBuffer = new MultidimensionalArray[NNKmax];
            CellNodesGlobal = new MultidimensionalArray[NNKmax];
            for(int i = 0; i < NNKmax; i++) {
                PhiEvalBuffer[i] = MultidimensionalArray.Create(1, NN);
                CellNodesGlobal[i] = MultidimensionalArray.Create(NN, D);
            }
        }


        public bool LocalSolve(int jCell, BitArray AcceptedMask, SinglePhaseField Phi, double _sign, out double Min, out double Max) {
            Debug.Assert(_sign.Abs() == 1);


            // find all accepted neighbors
            // ===========================

            var NeighCells = this.GridDat.GetCellNeighboursViaEdges(jCell);
            int iKref = this.GridDat.Cells.GetRefElementIndex(jCell);

            int NN = NeighCells.Length;
            var NeighCellsK = new Tuple<int, int, int>[NN];
            int NNK = 0;
            for(int nn = 0; nn < NN; nn++) {
                if(AcceptedMask[NeighCells[nn].Item1] == true) {
                    NeighCellsK[NNK] = NeighCells[nn];
                    NNK++;
                }
            }

            // evaluate accepted neighbors
            // ============================

            Min = double.MaxValue;
            Max = double.MinValue;

            var TrafoIdx = this.GridDat.Edges.Edge2CellTrafoIndex;
            int K = this.PhiEvalBuffer[0].GetLength(1); // Nodes per edge

            for(int nnk = 0; nnk < NNK; nnk++) { // loop over accepted neighbours
                int jNC = NeighCellsK[nnk].Item1;
                int iEdg = NeighCellsK[nnk].Item2;
                int InOrOut = NeighCellsK[nnk].Item3;
                int iTrafo = TrafoIdx[iEdg, InOrOut];

                this.PhiEvalBuffer[nnk].Clear();
                Phi.Evaluate(jNC, 1, this.EdgeNodes.GetVolumeNodeSet(this.GridDat, iTrafo), this.PhiEvalBuffer[nnk]);

                this.GridDat.TransformLocal2Global(this.EdgeNodes.GetVolumeNodeSet(this.GridDat, iTrafo), this.CellNodesGlobal[nnk], jNC);

                Max = Math.Max(Max, PhiEvalBuffer[nnk].Max());
                Min = Math.Min(Min, PhiEvalBuffer[nnk].Min());
            }

            var _QuadNodesGlobal = this.QuadNodesGlobal[iKref];
            this.GridDat.TransformLocal2Global(this.DaRuleS[iKref].Nodes, _QuadNodesGlobal, jCell);

            if(_sign > 0) {
                Max += this.GridDat.Cells.h_max[jCell];
            } else {
                Min -= this.GridDat.Cells.h_max[jCell];
            }


            // perform projection of geometric reinit
            // ======================================

            // temp storage
            double[] Y = new double[2];
            double[] X = new double[2];
            double[] X1 = new double[2];
            double[] X2 = new double[2];

            // basis values at cell quadrature nodes
            var BasisValues = this.LevelSetBasis_Geometric.CellEval(this.DaRuleS[iKref].Nodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);
            int NoOfQn = BasisValues.GetLength(0); // number of quadrature nodes

            // result at quadrature nodes
            var PhiAtQuadNodes = MultidimensionalArray.Create(NoOfQn);

            for(int iQn = NoOfQn - 1; iQn >= 0; iQn--) { // loop over all quadrature nodes

                Y[0] = _QuadNodesGlobal[iQn, 0];
                Y[1] = _QuadNodesGlobal[iQn, 1];

                double _dist_min1 = double.MaxValue, _phi_min1 = double.NaN;
                int _nnk_Min1 = int.MinValue, _k_min1 = int.MinValue;
                double _dist_min2 = double.MaxValue, _phi_min2 = double.NaN;
                int _nnk_Min2 = int.MinValue, _k_min2 = int.MinValue;


                // find closest point: brute force approach
                for(int nnk = 0; nnk < NNK; nnk++) { // loop over all edges with known values

                    double dist_min1 = double.MaxValue, phi_min1 = double.NaN;
                    int nnk_Min1 = int.MinValue, k_min1 = int.MinValue;
                    double dist_min2 = double.MaxValue, phi_min2 = double.NaN;
                    int nnk_Min2 = int.MinValue, k_min2 = int.MinValue;

                    for(int k = 0; k < K; k++) { // loop over all nodes on this edge

                        X[0] = this.CellNodesGlobal[nnk][k, 0];
                        X[1] = this.CellNodesGlobal[nnk][k, 1];
                        double phi = this.PhiEvalBuffer[nnk][0, k];
                        phi *= _sign;
                        double dist = GenericBlas.L2Dist(X, Y) + phi;

                        bool NoBlock = true;
                        if(dist < dist_min1) {
                            nnk_Min2 = nnk_Min1;
                            k_min2 = k_min1;
                            dist_min2 = dist_min1;
                            phi_min2 = phi_min1;

                            dist_min1 = dist;
                            nnk_Min1 = nnk;
                            k_min1 = k;
                            phi_min1 = phi;
                            NoBlock = false;
                        }

                        if(dist >= dist_min1 && dist < dist_min2 && NoBlock) {
                            dist_min2 = dist;
                            nnk_Min2 = nnk;
                            k_min2 = k;
                            phi_min2 = phi;
                        }


                    }

                    if(dist_min1 < _dist_min1) {
                        _dist_min1 = dist_min1;
                        _k_min1 = k_min1;
                        _nnk_Min1 = nnk_Min1;
                        _phi_min1 = phi_min1;

                        _dist_min2 = dist_min2;
                        _k_min2 = k_min2;
                        _nnk_Min2 = nnk_Min2;
                        _phi_min2 = phi_min2;
                    }
                }


                {
                    Debug.Assert(_nnk_Min1 == _nnk_Min2);
                    Debug.Assert(_k_min1 != _k_min2);

                    double PhiMin1 = this.PhiEvalBuffer[_nnk_Min1][0, _k_min1];
                    double PhiMin2 = this.PhiEvalBuffer[_nnk_Min2][0, _k_min2];

                    X1[0] = this.CellNodesGlobal[_nnk_Min1][_k_min1, 0];
                    X1[1] = this.CellNodesGlobal[_nnk_Min1][_k_min1, 1];

                    X2[0] = this.CellNodesGlobal[_nnk_Min2][_k_min2, 0];
                    X2[1] = this.CellNodesGlobal[_nnk_Min2][_k_min2, 1];


                }

                _dist_min1 *= _sign;

                PhiAtQuadNodes[iQn] = _dist_min1 * this.DaRuleS[iKref].Weights[iQn];
            }

            // finalize projection & return
            // ============================

            if(this.GridDat.Cells.IsCellAffineLinear(jCell)) {
                int N = this.LevelSetBasis_Geometric.GetLength(jCell);
                int N2 = Phi.Basis.GetLength(jCell);

                MultidimensionalArray Phi_1 = MultidimensionalArray.Create(N);
                double scale = this.GridDat.Cells.JacobiDet[jCell];
                Phi_1.Multiply(scale, BasisValues, PhiAtQuadNodes, 0.0, "m", "km", "k");
                for(int n = 0; n < N; n++)
                    Phi.Coordinates[jCell, n] = Phi_1[n];
                for(int n = N; n < N2; n++) {
                    Phi.Coordinates[jCell, n] = 0;
                }
            } else {
                throw new NotImplementedException("not implemented for curved cells");
            }

            return true;
        }

    }
}
