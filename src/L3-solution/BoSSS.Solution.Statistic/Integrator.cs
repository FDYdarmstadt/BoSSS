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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System.Collections;
using System.Collections.Generic;
using MPI.Wrappers;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Statistic {

    /// <summary>
    /// integrand, see <see cref="Integrator.PerformIntegration"/>;
    /// </summary>
    /// <param name="FieldValues">
    /// values of DG fields
    /// </param>
    /// <param name="iNode">
    /// quadrature rule node index
    /// </param>
    /// <param name="jCell">
    /// local cell index 
    /// </param>
    public delegate double Integrand(double[] FieldValues, int iNode, int jCell);
        
    /// <summary>
    /// General integration (in contrast to cell-wise integration, like for <see cref="BoSSS.Foundation.Quadrature"/>)
    /// </summary>
    public class Integrator {

        static double[] Linspace(double a, double b, int n) {
            if (a >= b)
                throw new ArgumentException("minimum >= maximum");

            double[] r = new double[n];
            double dx = (b - a) / ((double)(n - 1));
            for (int i = 0; i < n; i++) {
                r[i] = a + dx * ((double)i);
            }
            return r;
        }
        
        
        /// <summary>
        /// constructs an <see cref="Integrator"/> that integrates over the 
        /// circle with radius <paramref name="Radius"/> with center (<paramref name="Cen_x"/>,<paramref name="Cen_y"/>).
        /// The quadrature rule is of order 1 with <paramref name="Res"/> quadrature nodes.
        /// </summary>
        public static Integrator Circle(GridData grdDat, int Res, double Radius, double Cen_x = 0, double Cen_y = 0) {
            if (grdDat.SpatialDimension != 2)
                throw new NotSupportedException("this operation requires a 3D - grid");

            int D = 2;
            double[,] scheissQuadNodes = new double[Res, D];
            double[] depateGwichtln = new double[scheissQuadNodes.GetLength(0)];


            double dPhi = 2.0 * Math.PI / Res;
            for (int i = 0; i < Res; i++) {
                double Phi = dPhi*i;
                scheissQuadNodes[i, 0] = Math.Sin(Phi) * Radius;
                scheissQuadNodes[i, 1] = Math.Cos(Phi) * Radius;
            }

            double Umfang = Radius * 2.0 * Math.PI;
            depateGwichtln.SetAll(Umfang / Res);

            CellLocalization cl = new CellLocalization(grdDat);
            return new Integrator(cl, scheissQuadNodes, depateGwichtln);
        }


        /// <summary>
        /// constructs an <see cref="Integrator"/> that integrates over the 
        /// rectangle between the points (<paramref name="xMin"/>,<paramref name="yValue"/>,<paramref name="zMin"/>),
        /// and (<paramref name="xMax"/>,<paramref name="yValue"/>,<paramref name="zMax"/>) with surface normal (0,1,0) (i.e. an 'X-Z-plane',
        /// or an 'y=const. -- plane';<br/>
        /// The quadrature rule is of order 1 with <paramref name="xRes"/> times <paramref name="zRes"/> quadrature nodes.
        /// </summary>
        public static Integrator EquidistantOverXZPlane(GridData grdDat, double xMin, double zMin, double xMax, double zMax, double yValue, int xRes, int zRes) {
            if (grdDat.SpatialDimension != 3)
                throw new NotSupportedException("this operation requires a 3D - grid");
            
            CellLocalization cl = new CellLocalization(grdDat);
            
            double[] xNodes = Linspace(xMin, xMax, xRes+1);
            double[] zNodes = Linspace(zMin, zMax, zRes+1);

            double delta_x = xNodes[1] - xNodes[0];
            double delta_z = zNodes[1] - zNodes[0];

            double Area = (xMax - xMin) * (zMax - zMin);

            double wgt = delta_x * delta_z / Area;

            int D = 3;
            double[,] scheissQuadNodes = new double[xRes * zRes, D];
            double[] depateGwichtln = new double[scheissQuadNodes.GetLength(0)];

            int cnt = 0;
            for (int i = 0; i < xRes; i++) {
                for (int j = 0; j < zRes; j++) {

                    scheissQuadNodes[cnt, 0] = xNodes[i] + delta_x * 0.5;
                    scheissQuadNodes[cnt, 1] = yValue;
                    scheissQuadNodes[cnt, 2] = zNodes[i] + delta_z * 0.5;

                    depateGwichtln[cnt] = wgt;

                    cnt++;
                }
            }

            Integrator ret = new Integrator(cl, scheissQuadNodes, depateGwichtln);
            return ret;
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="loc"></param>
        /// <param name="QuadNodes">
        /// quad nodes outside the grid are ignored, i.e. the integrand is assumed to be 0
        /// </param>
        /// <param name="quadweights">
        /// </param>
        public Integrator(CellLocalization loc, double[,] QuadNodes, double[] quadweights) {
            int N = QuadNodes.GetLength(0);
            if (N != quadweights.Length)
                throw new ArgumentException("length of 0-th dimension of arrays must match.");
            int D = QuadNodes.GetLength(1);
            if (D != loc.GridBB.D)
                throw new ArgumentException("mismatch in spatial dimension.");
            int J = loc.GrdDat.Cells.NoOfLocalUpdatedCells;

            GridData gdat = loc.GrdDat;
            grdDat = gdat;
            
            // filter quad nodes outside of the bounding box
            // =============================================
            MultidimensionalArray QuadNodes2;
            double[] quadweights2;
            int[] orgindex;
            {
                double[] pt = new double[D];
                BitArray inside = new BitArray(N);

                int Found = 0;
                for (int n = 0; n < N; n++) {
                    for( int d = 0; d < D; d++)
                        pt[d] = QuadNodes[n,d];
                    bool ins = loc.GridBB.Contains(pt);

                    inside[n] = ins;
                    if (ins)
                        Found++;
                }

                QuadNodes2 = MultidimensionalArray.Create(Found, D);
                quadweights2 = new double[Found];
                orgindex = new int[Found];

                int cnt = 0;
                for (int n = 0; n < N; n++) {
                    if (inside[n]) {
                        for (int d = 0; d < D; d++)
                            QuadNodes2[cnt,d] = QuadNodes[n, d];
                        quadweights2[cnt] = quadweights[cnt];
                        orgindex[cnt] = n;
                        cnt++;
                    }
                }

                QuadNodes = null;
                quadweights = null;
                N = Found;
            }

            // build tree of quad nodes
            // ========================
            int[] Perm = new int[N];
            PointLocalization qn = new PointLocalization(QuadNodes2, loc.GridBB, Perm);
            double[] quadwegtNew = new double[N];
            int[] origindexNew = new int[N];
            for (int n = 0; n < N; n++) {
                quadwegtNew[n] = quadweights2[Perm[n]];
                origindexNew[n] = orgindex[Perm[n]];
            }
            quadweights2 = null; 

            // 1st index: cell index
            // 2nd index: quad node index within cell
            // 3rd index: spatial coordinate
            List<List<double[]>> QuadNodesPerCell = new List<List<double[]>>();

            // 1st index: cell index
            // 2nd index: quad node index within cell
            List<List<double>> QuadWeightsPerCell = new List<List<double>>();

            // 1st index: cell index
            // 2nd index: quad node index within cell
            List<List<int>> OrigIndexPerCell = new List<List<int>>();

            for (int j = 0; j < J; j++) {
                QuadNodesPerCell.Add(new List<double[]>());
                QuadWeightsPerCell.Add(new List<double>());
                OrigIndexPerCell.Add(new List<int>());
            }
                        
            // try to assign the quad nodes to cells
            // =====================================
            BitArray PointsLocatedMarker = new BitArray(N); // mark every node, that we assign to a cell with 
            int NoOfUnassignedNodes = N;

            //int[] Cell4Quadnodes = new int[N];
            
            // loop over cells ...
            MultidimensionalArray vertGlobal = new MultidimensionalArray(2);
            MultidimensionalArray vertLocal = new MultidimensionalArray(3);
            for (int j = 0; j < J; j++) {

                //if (loc.CellMaxCode[j] < Locations[0])
                //    continue; // skip the cell: contains none of the searched points
                //if (loc.CellMinCode[j] > Locations[N - 1])
                //    continue; // skip the cell: contains none of the searched points

                GeomBinTreeBranchCode bbcode; int bbBits;
                {
                    BoundingBoxCode __b = loc.GetCellBoundingBoxCode(j);
                    bbcode = __b.Branch;
                    bbBits = (int) __b.SignificantBits;
                }

                int iP0, Len;
                qn.GetPointsInBranch(bbcode, bbBits, out iP0, out Len);
                if (Len <= 0)
                    continue;

                vertGlobal.Allocate(Len, D);
                vertLocal.Allocate(1, Len, D);
                for (int n = 0; n < Len; n++)
                    for (int d = 0; d < D; d++)
                        vertGlobal[n, d] = qn.Points[n + iP0, d];
                gdat.TransformGlobal2Local(vertGlobal, vertLocal, j, 1, 0);


                var splx = gdat.Cells.GetRefElement(j);
                for (int n = 0; n < Len; n++) {
                    int nPt = n + iP0;
                    if (PointsLocatedMarker[nPt])
                        continue;

                    double[] pt = new double[D];
                    for (int d = 0; d < D; d++) pt[d] = vertLocal[0, n, d];
                    if (splx.IsWithin(pt)) {
                        PointsLocatedMarker[nPt] = true;
                        NoOfUnassignedNodes--;
                        //Cell4Quadnodes[nPt] = j;

                        QuadNodesPerCell[j].Add(pt);
                        QuadWeightsPerCell[j].Add(quadwegtNew[nPt]);
                        OrigIndexPerCell[j].Add(origindexNew[nPt]);
                    }
                }
            }

            // record final data structures
            // ============================
            //m_QuadNodesPerCell = new MultidimensionalArray[J];
            m_QuadNodesPerCell = new NodeSet[J];
            m_QuadWeightsPerCell = new double[J][];
            m_OriginalQuadNodesIndex = new int[J][];
            MaxNumberOfNodes = 0;
            for (int j = 0; j < J; j++) {
                List<double[]> NodesInCellj = QuadNodesPerCell[j];
                                
                m_QuadWeightsPerCell[j] = QuadWeightsPerCell[j].ToArray();
                m_OriginalQuadNodesIndex[j] = OrigIndexPerCell[j].ToArray();

                int NJ = NodesInCellj.Count;
                if (NJ > 0) {
                    var _QuadNodesPerCell_j = new NodeSet(grdDat.Cells.GetRefElement(j), NJ, D);
                    MaxNumberOfNodes = Math.Max(MaxNumberOfNodes, NJ);
                    for (int nn = 0; nn < NJ; nn++) {
                        for (int d = 0; d < D; d++) {
                            _QuadNodesPerCell_j[nn, d] = NodesInCellj[nn][d];
                        }
                    }
                    _QuadNodesPerCell_j.LockForever();

                    m_QuadNodesPerCell[j] = _QuadNodesPerCell_j;
                } else {
                    m_QuadNodesPerCell[j] = null;
                }
            }
        }

        GridData grdDat;

        int MaxNumberOfNodes;

        /// <summary>
        ///  - 1st index: cell index <br/>
        ///  - 2nd index: quad node index within cell <br/>
        ///  - 3rd index: spatial coordinate
        /// </summary>
        NodeSet[] m_QuadNodesPerCell;

        /// <summary>
        /// for each quad node, its original index in the quad rule that was provided with the constructor.
        ///  - 1st index: cell index <br/>
        ///  - 2nd index: quad node index within cell <br/>
        /// </summary>
        int[][] m_OriginalQuadNodesIndex;


        /// <summary>
        ///  - 1st index: cell index <br/>
        ///  - 2nd index: quad node index within cell
        /// </summary>
        double[][] m_QuadWeightsPerCell;
        

        /// <summary>
        /// Performs an integration of the integrand <paramref name="I"/>;<br/>
        /// (the domain of the integration is an issue of the quadrature rule specified in the constructor).
        /// </summary>
        /// <param name="I">
        /// integrand; for each quadrature node <em>x</em>, the values of the DG fields <paramref name="fields"/>
        /// at point <em>x</em> are passed to <paramref name="I"/>;
        /// </param>
        /// <param name="fields">
        /// DG fields which are the input of integrand <paramref name="I"/>;
        /// </param>
        /// <returns>
        /// the result of the integration/quadrature (equal on all MPI processors)
        /// </returns>
        /// <remarks>
        /// this call is MPI-collective, i.e. the returned result is the same on all MPI processors
        /// </remarks>
        public double PerformIntegration(Integrand I, params DGField[] fields) {

            //// temp vars.
            //MultidimensionalArray[] fieldValues = new MultidimensionalArray[fields.Length];
            //for (int fld = 0; fld < fieldValues.Length; fld++)
            //    fieldValues[fld] = MultidimensionalArray.Create(1, MaxNumberOfNodes);

            int FLD = fields.Length;
            int J = this.grdDat.Cells.NoOfLocalUpdatedCells;
            double[] ladygaga = new double[FLD];

            // loop over cells ...
            double ResultAcc = 0; // accumulator for result of quadrature
            for (int j = 0; j < J; j++) {

                if (m_QuadNodesPerCell[j] == null)
                    continue;

                // evaluate fields
                // ===============

                
                // temp vars.
                MultidimensionalArray[] fieldValues = new MultidimensionalArray[fields.Length];
                for (int fld = 0; fld < fieldValues.Length; fld++)
                    fieldValues[fld] = MultidimensionalArray.Create(1, m_QuadNodesPerCell[j].NoOfNodes);
                
                // eval
                for (int fld = 0; fld < FLD; fld++) {
                    fieldValues[fld].Clear();
                    fields[fld].Evaluate(j, 1, m_QuadNodesPerCell[j], fieldValues[fld]);
                }
               
                // evaluate integrand 'I' and do quadrature
                // ========================================
                double[] quadweights = m_QuadWeightsPerCell[j];
                int[] orgnodes = m_OriginalQuadNodesIndex[j];
                int N = quadweights.Length;
                for (int n = 0; n < N; n++) {
                    for (int fld = 0; fld < FLD; fld++)
                        ladygaga[fld] = fieldValues[fld][0, n];
                    
                    double val = I(ladygaga,orgnodes[n],j);

                    ResultAcc += val * quadweights[n];
                }
            }

            // reduce over all MPI processes
            unsafe {
                double GlobalAcc = 0;
                csMPI.Raw.Allreduce((IntPtr)(&ResultAcc), (IntPtr)(&GlobalAcc), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                return GlobalAcc;    
            }
        }
    }
}
