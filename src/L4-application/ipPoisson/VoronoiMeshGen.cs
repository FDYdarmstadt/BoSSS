//            _____                   _______                   _____                    _____                    _____                  
//           /\    \                 /::\    \                 /\    \                  /\    \                  /\    \                 
//          /::\    \               /::::\    \               /::\    \                /::\    \                /::\    \                
//         /::::\    \             /::::::\    \             /::::\    \              /::::\    \              /::::\    \               
//        /::::::\    \           /::::::::\    \           /::::::\    \            /::::::\    \            /::::::\    \              
//       /:::/\:::\    \         /:::/~~\:::\    \         /:::/\:::\    \          /:::/\:::\    \          /:::/\:::\    \             
//      /:::/__\:::\    \       /:::/    \:::\    \       /:::/__\:::\    \        /:::/__\:::\    \        /:::/__\:::\    \            
//     /::::\   \:::\    \     /:::/    / \:::\    \      \:::\   \:::\    \       \:::\   \:::\    \       \:::\   \:::\    \           
//    /::::::\   \:::\    \   /:::/____/   \:::\____\   ___\:::\   \:::\    \    ___\:::\   \:::\    \    ___\:::\   \:::\    \          
//   /:::/\:::\   \:::\ ___\ |:::|    |     |:::|    | /\   \:::\   \:::\    \  /\   \:::\   \:::\    \  /\   \:::\   \:::\    \         
//  /:::/__\:::\   \:::|    ||:::|____|     |:::|    |/::\   \:::\   \:::\____\/::\   \:::\   \:::\____\/::\   \:::\   \:::\____\        
//  \:::\   \:::\  /:::|____| \:::\    \   /:::/    / \:::\   \:::\   \::/    /\:::\   \:::\   \::/    /\:::\   \:::\   \::/    /        
//   \:::\   \:::\/:::/    /   \:::\    \ /:::/    /   \:::\   \:::\   \/____/  \:::\   \:::\   \/____/  \:::\   \:::\   \/____/         
//    \:::\   \::::::/    /     \:::\    /:::/    /     \:::\   \:::\    \       \:::\   \:::\    \       \:::\   \:::\    \             
//     \:::\   \::::/    /       \:::\__/:::/    /       \:::\   \:::\____\       \:::\   \:::\____\       \:::\   \:::\____\            
//      \:::\  /:::/    /         \::::::::/    /         \:::\  /:::/    /        \:::\  /:::/    /        \:::\  /:::/    /            
//       \:::\/:::/    /           \::::::/    /           \:::\/:::/    /          \:::\/:::/    /          \:::\/:::/    /             
//        \::::::/    /             \::::/    /             \::::::/    /            \::::::/    /            \::::::/    /              
//         \::::/    /               \::/____/               \::::/    /              \::::/    /              \::::/    /               
//          \::/____/                 ~~                      \::/    /                \::/    /                \::/    /                
//           ~~                                                \/____/                  \/____/                  \/____/                 
                                                                                                                                    

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using ilPSP.Connectors.Matlab;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.SipPoisson.Voronoi {
    static class VoronoiMeshGen {


        static int Mirror(ref double[] _x, ref double[] _y, AffineManifold[] bndys, Func<Vector, bool> _IsIn) {
            bool IsIn(double xp, double yp) {
                return _IsIn(new Vector(xp, yp));
            }


            if (_x.Length != _y.Length)
                throw new ArgumentException();
            var x = _x.ToList();
            var y = _y.ToList();
            int N = _x.Length;



            // filter all points that are outside of the domain
            for (int n = 0; n < N; n++) {
                if (!IsIn(x[n], y[n])) {
                    x.RemoveAt(n);
                    y.RemoveAt(n);
                    N--;
                    n--;
                }
            }
            Debug.Assert(x.Count == N);
            Debug.Assert(y.Count == N);
            for (int n = 0; n < N; n++) {
                Debug.Assert(IsIn(x[n], y[n]));
            }

            // mirror each point 
            for (int n = 0; n < N; n++) {
                double xn = x[n];
                double yn = y[n];
                for (int l = 0; l < bndys.Length; l++) {
                    var bndy_l = bndys[l];

                    double dist = bndy_l.PointDistance(xn, yn);

                    if (dist < 0) {
                        double xMirr = xn - bndy_l.Normal[0] * dist * 2;
                        double yMirr = yn - bndy_l.Normal[1] * dist * 2;

                        Debug.Assert(bndy_l.PointDistance(xMirr, yMirr) > 0);

                        if (!IsIn(xMirr, yMirr)) {
                            x.Add(xMirr);
                            y.Add(yMirr);
                        }
                    }
                }
            }

            // return
            _x = x.ToArray();
            _y = y.ToArray();
            return N;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="Nodes">
        /// Delaunay vertices
        /// - 1st index: vetrx index
        /// - 2nd index: 0,1 for x,y coordinate, respectively
        /// </param>
        /// <param name="PolygonBoundary"></param>
        /// <returns></returns>
        static public AggregationGrid FromPolygonalDomain(MultidimensionalArray Nodes, Vector[] PolygonBoundary, Func<Vector,bool> IsIn) {
            
            // check arguments
            // ===============

            if (Nodes.Dimension != 2)
                throw new ArgumentException("expecting 2D array;");
            if (Nodes.GetLength(1) != 2)
                throw new ArgumentException("only implemented for 2D");

            foreach(var V in PolygonBoundary) {
                if (V.Dim != 2)
                    throw new ArgumentException();
                if(!IsIn(V)) {
                    throw new ArgumentException();
                }
            }

            // boundaries for domain
            // =====================

            AffineManifold[] Boundaries = new AffineManifold[PolygonBoundary.Length];
            for(int i = 0; i < Boundaries.Length; i++) {
                Vector O = PolygonBoundary[i];
                Vector E = PolygonBoundary[(i + 1) % Boundaries.Length];
                Vector OE = E - O;

                Vector N = new Vector(-OE.y, OE.x);
                N.Normalize();

                Boundaries[i] = new AffineManifold(N, O);

                double eps = 1.0e-6;
                Vector Cen = O  + 0.5*OE;
                Vector IN = Cen - eps * N;
                Vector OT = Cen + eps * N;
                Debug.Assert(Boundaries[i].PointDistance(IN) < 0);
                Debug.Assert(IsIn(IN) == true);
                Debug.Assert(Boundaries[i].PointDistance(OT) > 0);
                Debug.Assert(IsIn(OT) == false);

            }

            // Vertex Mirroring
            // ================
            MultidimensionalArray NewNodes;
            {
                /*
                double[] xNodes = Nodes.GetColumn(0);
                double[] yNodes = Nodes.GetColumn(1);

                Mirror(ref xNodes, ref yNodes, Boundaries, IsIn);
                Debug.Assert(xNodes.Length == yNodes.Length);

                NewNodes = MultidimensionalArray.Create(xNodes.Length, 2);
                NewNodes.SetColumn(0, xNodes);
                NewNodes.SetColumn(1, yNodes);
                */

                NewNodes = Nodes.CloneAs();
                Nodes = null;
            }

            // Create Voronoi mesh (call Matlab)
            // =================================

            int[][] OutputVertexIndex;
            MultidimensionalArray VertexCoordinates;
            using (var Matlab = new BatchmodeConnector()) {


                Matlab.PutMatrix(NewNodes, "Nodes");

                // compute Voronoi diagramm
                Matlab.Cmd("[V, C] = voronoin(Nodes);");

                // output (export from matlab)
                OutputVertexIndex = new int[NewNodes.NoOfRows][];
                Matlab.GetStaggeredIntArray(OutputVertexIndex, "C");
                Matlab.GetMatrix(null, "V");

                // run matlab
                Matlab.Execute(false);

                // import here
                VertexCoordinates = (MultidimensionalArray)(Matlab.OutputObjects["V"]);

                // correct indices (1-based index to 0-based index)
                foreach (int[] cell in OutputVertexIndex) {
                    int K = cell.Length;
                    for (int k = 0; k < K; k++) {
                        cell[k]--;
                    }
                }
            }
                        

            // tessellation
            // ============

            List<Cell> cells = new List<Cell>();
            List<int[]> aggregation = new List<int[]>();
            for (int jV = 0; jV < OutputVertexIndex.Length; jV++) { // loop over Voronoi Cells
                
                int[] iVtxS = OutputVertexIndex[jV];
                int NV = iVtxS.Length;

                Vector[] VoronoiCell = iVtxS.Select(iVtx => VertexCoordinates.GetRowPt(iVtx)).ToArray();
                
                bool AnyIn = VoronoiCell.Any(V => IsIn(V));
                bool AnyOt = VoronoiCell.Any(V => !IsIn(V));

                if (AnyIn) {

                    FixOrientation(ref VoronoiCell, ref iVtxS);

                    int[,] iVtxTri;
                    
                    if(AnyOt) {
                        // ++++++++++++++++++++
                        // clipping required
                        // ++++++++++++++++++++

                        var VoronoiCellOld = VoronoiCell;
                        VoronoiCell = PolygonClipping.WeilerAthertonClipping(PolygonBoundary, IsIn, VoronoiCellOld);


                        iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell);
                        //iVtxTri = PolygonTesselation.TesselateConvexPolygon(VoronoiCell);

                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // all vertices inside - standard convex polygon tesselation
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        iVtxTri = PolygonTesselation.TesselateConvexPolygon(VoronoiCell);
                    }



                    List<int> Agg2Pt = new List<int>();

                    for (int iTri = 0; iTri < iVtxTri.GetLength(0); iTri++) { // loop over triangles of voronoi cell
                        int iV0 = iVtxTri[iTri, 0];
                        int iV1 = iVtxTri[iTri, 1];
                        int iV2 = iVtxTri[iTri, 2];

                        Vector V0 = VoronoiCell[iV0];
                        Vector V1 = VoronoiCell[iV1];
                        Vector V2 = VoronoiCell[iV2];

                        Vector D1 = V1 - V0;
                        Vector D2 = V2 - V0;
                        Debug.Assert(D1.CrossProduct2D(D2) > 1.0e-8);


                        Cell Cj = new Cell();
                        Cj.GlobalID = cells.Count;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                        Cj.NodeIndices = new int[] { iVtxS[iV0], iVtxS[iV1], iVtxS[iV2] };
                        Cj.TransformationParams.SetRowPt(0, V0);
                        Cj.TransformationParams.SetRowPt(1, V1);
                        Cj.TransformationParams.SetRowPt(2, V2);

                        Agg2Pt.Add(cells.Count);

                        cells.Add(Cj);
                    }

                    aggregation.Add(Agg2Pt.ToArray());
                } else {
                    // nop
                    double xMin = VoronoiCell.Min(V => V.x);
                    double xMax = VoronoiCell.Max(V => V.x);
                    double yMin = VoronoiCell.Min(V => V.y);
                    double yMax = VoronoiCell.Max(V => V.y);
                    //if (xMax < 1.1 && xMin > -1.1 && yMax < 1.1 && yMax > -1.1) {
                    //    Console.WriteLine(" Say Tschus to cell {0} : {1} // {2} : {3}", xMin, xMax, yMin, yMax);
                    //}
                }
            }

            // return grid
            // ===========

            // base grid
            GridCommons grd;
            grd = new Grid2D(Triangle.Instance);
            grd.Cells = cells.ToArray();
            grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
            //grd.Plot2DGrid();
            grd.DefineEdgeTags(X => (byte)1);

            // aggregation grid
            var agrd = new AggregationGrid(grd, aggregation.ToArray());
            return agrd;

        }


        static void FixOrientation(ref Vector[] Polygon, ref int[] iVtx) {
            int L = Polygon.Length;
            
            double[] signs = new double[L - 2];

            bool AllPos = true;
            bool AllNeg = true;

            for (int iTri = 0; iTri < L - 2; iTri++) { // loop over triangles of voronoi cell
                int iV0 = 0;
                int iV1 = iTri + 1;
                int iV2 = iTri + 2;

                Vector V0 = Polygon[iV0];
                Vector V1 = Polygon[iV1];
                Vector V2 = Polygon[iV2];

                Vector D1 = V1 - V0;
                Vector D2 = V2 - V0;

                signs[iTri] = D1.CrossProduct2D(D2);

                AllPos = AllPos && (signs[iTri] > 0);
                AllNeg = AllNeg && (signs[iTri] < 0);
            }

            if (AllNeg == AllPos)
                throw new ArithmeticException("Indefinite polygon");

            if (AllPos)
                return;

            if (AllNeg) {
                Polygon = Polygon.Reverse().ToArray();
                iVtx = iVtx.Reverse().ToArray();
                return;
            }

            throw new ArithmeticException("Indefinite polygon");
        }

    }
}
