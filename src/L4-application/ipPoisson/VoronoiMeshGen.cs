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
using ilPSP.Utils;
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


        enum VertexType {
            unspecified = 0,
            
            Inside = 1,

            Outside = 2,

            OnBoundaryplane_Inside = 3,

            //OnBoundaryplane_Outside = 4,

            OnCorner = 5,

            FarPoint = 6 // somewhere at infty
        }


        class VoronoiVertex {
            public Vector VTX;

            public VertexType type;




        }



        class VoronoiEdge {
            /// <summary>
            /// First vertex in voronoi cell
            /// </summary>
            public int iVtxA;

            /// <summary>
            /// Second vertex in voronoi cell
            /// </summary>
            public int iVtxB;


            public List<int> Cells = new List<int>();


            public override bool Equals(object obj) {
                if (iVtxA == iVtxB)
                    throw new ApplicationException();

                var E2 = obj as VoronoiEdge;
                if (iVtxA == E2.iVtxA && iVtxB == E2.iVtxB)
                    return true;
                if (iVtxA == E2.iVtxB && iVtxB == E2.iVtxA)
                    return true;

                return false;
            }


            public override int GetHashCode() {
                return iVtxA + (iVtxB << 16);
            }
        }


        static double CoordOnLine(Vector P, Vector A, Vector B) {
            Vector AP = P - A;
            Vector AB = B - A;

            if (AB.AbsSquare() <= 0.0)
                throw new ArgumentException();

            if (AP.AbsSquare() <= 0.0)
                return 0.0;

            return (AB*AP) / AB.AbsSquare();
            
        }

        static bool isFarPoint(Vector V) {
            if (double.IsInfinity(V.x))
                return true;
            if (double.IsInfinity(V.x))
                return true;

            return false;
        }


        /*
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
        static public AggregationGrid FromPolygonalDomain(MultidimensionalArray Nodes, Vector[] PolygonBoundary, Func<Vector,bool> IsIn, Func<Vector,Vector,bool> PointIdentity) {
            
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

            // check if any two boundary vertices are identical:
            for (int iBnd = 0; iBnd < PolygonBoundary.Length; iBnd++) {
                for (int iBnd2 = iBnd + 1; iBnd2 < PolygonBoundary.Length; iBnd2++) {
                    if (PointIdentity(PolygonBoundary[iBnd], PolygonBoundary[iBnd2]))
                        throw new ArgumentException("degenerate boundary");
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

                double eps = Math.Sqrt(BLAS.MachineEps)*10;
                Vector Cen = O  + 0.5*OE;
                Vector IN = Cen - eps * N;
                Vector OT = Cen + eps * N;
                Debug.Assert(Boundaries[i].PointDistance(IN) < 0);
                Debug.Assert(IsIn(IN) == true);
                Debug.Assert(Boundaries[i].PointDistance(OT) > 0);
                Debug.Assert(IsIn(OT) == false);


                Vector CenProj = Boundaries[i].ProjectPoint(Cen);
                if (PointIdentity(Cen, CenProj) == false)
                    throw new ArithmeticException("point identity does not seem to work");

                Debug.Assert(PointIdentity(Cen, Boundaries[i].ProjectPoint(Cen + N * OE.Abs()))); // tests if the 'ProjectPoint' works correctly

            }

            // Vertex Mirroring
            // ================
            MultidimensionalArray NewNodes;
            {
                
                double[] xNodes = Nodes.GetColumn(0);
                double[] yNodes = Nodes.GetColumn(1);

                Mirror(ref xNodes, ref yNodes, Boundaries, IsIn);
                Debug.Assert(xNodes.Length == yNodes.Length);

                NewNodes = MultidimensionalArray.Create(xNodes.Length, 2);
                NewNodes.SetColumn(0, xNodes);
                NewNodes.SetColumn(1, yNodes);
                

                //NewNodes = Nodes.CloneAs();
                Nodes = null;
            }

            // Create Voronoi mesh (call Matlab)
            // =================================

            int[][] VocellVertexIndex;
            List<VoronoiVertex> Verts;
            using (var Matlab = new BatchmodeConnector()) {


                Matlab.PutMatrix(NewNodes, "Nodes");

                // compute Voronoi diagramm
                Matlab.Cmd("[V, C] = voronoin(Nodes);");

                // output (export from matlab)
                VocellVertexIndex = new int[NewNodes.NoOfRows][];
                Matlab.GetStaggeredIntArray(VocellVertexIndex, "C");
                Matlab.GetMatrix(null, "V");

                // run matlab
                Matlab.Execute(false);

                // import here
                MultidimensionalArray VertexCoordinates = (MultidimensionalArray)(Matlab.OutputObjects["V"]);
                Verts = new List<VoronoiVertex>();
                for(int iVtx = 0; iVtx < VertexCoordinates.NoOfRows; iVtx++) {
                    Verts.Add(new VoronoiVertex() { VTX = VertexCoordinates.GetRowPt(iVtx) });
                }

                // correct indices (1-based index to 0-based index)
                foreach (int[] cell in VocellVertexIndex) {
                    int K = cell.Length;
                    for (int k = 0; k < K; k++) {
                        cell[k]--;
                    }
                }

                // fix voronoi cell orientation
                for (int jV = 0; jV < VocellVertexIndex.Length; jV++) {
                    int[] iVtxS = VocellVertexIndex[jV];
                    if (iVtxS.Any(i => isFarPoint(Verts[i].VTX)))
                        continue;

                    Vector[] vectors = iVtxS.Select(i => Verts[i].VTX).ToArray();
                    FixOrientation(ref vectors, ref iVtxS);
                    VocellVertexIndex[jV] = iVtxS;
                }
            }

            // vertex classification
            // =====================
            {
                for (int iVtx = 0; iVtx < Verts.Count; iVtx++) {
                    var Vert = Verts[iVtx];

                    // check if the vertex coincides with a corner
                    Debug.Assert(Vert.type == VertexType.unspecified);
                    if(isFarPoint(Vert.VTX)) {
                        Vert.type = VertexType.FarPoint;
                        continue;
                    }

                    // check if the vertex coincides with a corner
                    Debug.Assert(Vert.type == VertexType.unspecified);
                    for(int iBnd = 0; iBnd < PolygonBoundary.Length; iBnd++) {
                        Vector VB = PolygonBoundary[iBnd];

                        if(PointIdentity(VB, Vert.VTX)) {
                            if (Vert.type != VertexType.unspecified)
                                throw new ArithmeticException();
                            Vert.type = VertexType.OnCorner;
                        }
                    }

                    if (Vert.type == VertexType.OnCorner)
                        continue;

                    // check if the vertex is on one of the boundary sides
                    Debug.Assert(Vert.type == VertexType.unspecified);
                    for (int iBnd = 0; iBnd < PolygonBoundary.Length; iBnd++) {
                        AffineManifold bndy = Boundaries[iBnd];


                        var vaProj = Boundaries[iBnd].ProjectPoint(Vert.VTX);
                        if (PointIdentity(vaProj, Vert.VTX)) {
                            if (Vert.type != VertexType.unspecified)
                                throw new ArithmeticException();

                            // point is identical with its projection onto the boundary
                            // => it must be on the on one of the sides

                            Vector cornerA = PolygonBoundary[iBnd];
                            Vector cornerB = PolygonBoundary[(iBnd + 1) % PolygonBoundary.Length];

                            if(PointIdentity(vaProj, cornerA)) {
                                Vert.type = VertexType.OnCorner;
                                break;
                            }
                            if(PointIdentity(vaProj, cornerB)) {
                                Vert.type = VertexType.OnCorner;
                                break;
                            }

                            double alpha = CoordOnLine(vaProj, cornerA, cornerB);
                            Debug.Assert(PointIdentity(vaProj, cornerA * (1 - alpha) + cornerB * alpha));
                            if (alpha == 0.0)
                                throw new ArithmeticException("point identity seems inconsistent"); // should have been handled above
                            if (alpha == 1.0)
                                throw new ArithmeticException("point identity seems inconsistent"); // should have been handled above

                            if(alpha > 0 && alpha < 1.0) {
                                Vert.type = VertexType.OnBoundaryplane_Inside;
                                break;
                            }

                            //if (IsIn(Vert.VTX))
                            //    Vert.type = VertexType.OnBoundaryplane_Inside;
                            //else
                            //    Vert.type = VertexType.OnBoundaryplane_Outside;



                        }

                    }

                    if (Vert.type == VertexType.OnBoundaryplane_Inside || Vert.type == VertexType.OnCorner)
                        continue;

                    // check whether the vertex is inside or outside
                    Debug.Assert(Vert.type == VertexType.unspecified);
                    if (IsIn(Vert.VTX))
                        Vert.type = VertexType.Inside;
                    else
                        Vert.type = VertexType.Outside;
                }
            }


            // detect edges
            // ============

            {
                List<VoronoiEdge> A = new List<VoronoiEdge>();
                A.Add(new VoronoiEdge() { iVtxA = 1, iVtxB = 2 });
                A.Add(new VoronoiEdge() { iVtxA = 3, iVtxB = 2 });
                A.Add(new VoronoiEdge() { iVtxA = 4, iVtxB = 1 });
                A.Add(new VoronoiEdge() { iVtxA = 2, iVtxB = 1 });

                Debug.Assert(A.IndexOf(new VoronoiEdge() { iVtxA = 4, iVtxB = 1 }) == 2);
                Debug.Assert(A.IndexOf(new VoronoiEdge() { iVtxA = 1, iVtxB = 4 }) == 2);

            }

            List<int>[] Cell2Edge; // cell-2-edge map: 1st index voronoi cell; 2nd index: enum;
            List<VoronoiEdge> Edges; // index: edge index
            {
                Edges = new List<VoronoiEdge>();
                Cell2Edge  = VocellVertexIndex.Length.ForLoop(i => new List<int>());
                bool anyDouble = false;
                for (int jV = 0; jV < VocellVertexIndex.Length; jV++) {
                    int[] iVtxS = VocellVertexIndex[jV];

                    int I = iVtxS.Length;
                    for (int i = 0; i < I; i++) {
                        VoronoiEdge newEdge = new VoronoiEdge() {
                            iVtxA = iVtxS[i],
                            iVtxB = iVtxS[(i + 1) % I]
                        };

                        int iEdge = Edges.IndexOf(newEdge);
                        if(iEdge < 0) {
                            for(int e = 0; e < Edges.Count; e++) {
                                Debug.Assert(Edges[e].Equals(newEdge) == false);
                            }

                            Edges.Add(newEdge);
                            iEdge = Edges.Count - 1;
                        } else {
                            anyDouble = true;
                            if (Edges[iEdge].Cells.Count > 1)
                                throw new ArithmeticException("edge shared by more than two cells");
                            newEdge = Edges[iEdge];
                        }
                        Debug.Assert(newEdge.Equals(Edges[iEdge]));
                        Debug.Assert(Edges.Where(ve => ve.Equals(newEdge)).Count() == 1);

                        Cell2Edge[jV].Add(iEdge);
                        newEdge.Cells.Add(jV);
                    }
                }

                if (!anyDouble)
                    throw new ArithmeticException("Voronoi diagram seems to be completely disjoint - no edge is used by at least two cells.");
            }

            // intersections of edges and boundary
            // ===================================

            {

                // test all edges against all boundaries

                int E = Edges.Count;
                for(int e = 0; e < E; e++) { // loop over edges...
                    var Edge = Edges[e];
                    VoronoiVertex vA = Verts[Edge.iVtxA];
                    VoronoiVertex vB = Verts[Edge.iVtxB];

                    if(vA.type == VertexType.FarPoint || vB.type == VertexType.FarPoint) {
                        // cant deal with those guys at the moment
                        continue;
                    }

                    for(int iBnd = 0; iBnd < Boundaries.Length; iBnd++) { // loop over boundaries...
                        AffineManifold Bndy_i = Boundaries[iBnd];
                        Vector B1 = PolygonBoundary[iBnd];
                        Vector B2 = PolygonBoundary[(iBnd + 1)%Boundaries.Length];

                        var vaProj = Boundaries[iBnd].ProjectPoint(vA.VTX);
                        bool vAinPlane = PointIdentity(vaProj, vA.VTX);

                        var vbProj = Boundaries[iBnd].ProjectPoint(vB.VTX);
                        bool vBinPlane = PointIdentity(vbProj, vB.VTX);
                        // -----

                        if(vAinPlane && vBinPlane) {
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // special case: Voronoi edge and boundary are parallel
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                            Console.Write("");
                        } else {

                            bool cutfound = PolygonClipping.ComputeIntersection(vA.VTX, vB.VTX, B1, B2, out double alpha1, out double alpha2, out Vector I);

                            if(PointIdentity(I, vA.VTX)) {
                                alpha1 = 0;
                            }
                            if(PointIdentity(I, vB.VTX)) {
                                alpha1 = 1;
                            }

                            if(PointIdentity(I, B1)) {
                                alpha2 = 0;
                            }
                            if(PointIdentity(I, B2)) {
                                alpha2 = 1;
                            }


                            if(cutfound) {
                                // Voronoi edge and boundary are NOT parallel

                                if(alpha1 >= 0 && alpha1 <= 1 && alpha2 >= 0 && alpha2 <= 1) {
                                    // 'alpha1' is relative coordinate on the boundary
                                    // 'alpha2' is relative coordinate on the edge

                                    Verts.Add(new VoronoiVertex() {
                                        VTX = I,
                                        type = VertexType.OnBoundaryplane_Inside
                                    });
                                    int iNew = Verts.Count - 1;

                                }


                                
                                
                                
                            }
                        }
                    }
                }
            }


            //
            //
            List<int[]> ClippedPolygons;
            {
                ClippedPolygons = new List<int[]>();

                Debug.Assert(Cell2Edge.Length == VocellVertexIndex.Length);
                for(int jV = 0; jV < Cell2Edge.Length; jV++) {
                    int[] Edges_jV = Cell2Edge[jV].ToArray();
                    int[] VtxS = VocellVertexIndex[jV];

                    if (VtxS.Any(iVtx => Verts[iVtx].type == VertexType.FarPoint)) {
                        // don't know how to deal with those guys yet.
                        continue;
                    }

                    

#if DEBUG
                    {
                        var VtxAset = Edges_jV.Select(edge => Edges[edge].iVtxA);
                        var VtxBset = Edges_jV.Select(edge => Edges[edge].iVtxB);
                        var VtxSet = VtxAset.SetUnion(VtxBset);

                        Debug.Assert(VtxS.SetEquals(VtxSet));
                    }
#endif

                    
                    for(int i = 0; i < VtxS.Length; i++) {
                        int _iVtxA = VtxS[i];
                        int _iVtxB = VtxS[i % VtxS.Length];
                        int iEdge = Edges_jV.Single(e => Edges[e].Equals(new VoronoiEdge() { iVtxA = _iVtxA, iVtxB = _iVtxB }));
                        VoronoiEdge edge = Edges[iEdge];


                    }

                }

            }

            // tessellation of clipped polygons
            // ================================


            List<Cell> cells = new List<Cell>();
            List<int[]> aggregation = new List<int[]>();
            for (int jV = 0; jV < VocellVertexIndex.Length; jV++) { // loop over Voronoi Cells
                
                int[] iVtxS = VocellVertexIndex[jV];
                int NV = iVtxS.Length;

                Vector[] VoronoiCell = null;// iVtxS.Select(iVtx => VertexCoordinates.GetRowPt(iVtx)).ToArray();
                
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

    */
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
