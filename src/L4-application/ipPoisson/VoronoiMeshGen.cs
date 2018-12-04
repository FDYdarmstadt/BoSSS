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

        /*
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
        */

        class VoItem {
            //public int Index;
            public bool deleted = false;
        }

        class VoVertex : VoItem {
            public Vector VTX;

            public bool IsFar {
                get {
                    return IsFarPoint(VTX);
                }
            }
        }



        class VoEdge : VoItem {
            /// <summary>
            /// First vertex in voronoi cell
            /// </summary>
            public VoVertex VtxA;

            /// <summary>
            /// Second vertex in voronoi cell
            /// </summary>
            public VoVertex VtxB;


            public List<VoPolygon> Cells = new List<VoPolygon>();

            public bool isBoundary;


            public override bool Equals(object obj) {
                /*
                if (VtxA.Index == VtxB.Index)
                    throw new ApplicationException();

                var E2 = obj as VoEdge;
                if (VtxA.Index == E2.VtxA.Index && VtxB.Index == E2.VtxB.Index)
                    return true;
                if (VtxA.Index == E2.VtxB.Index && VtxB.Index == E2.VtxA.Index)
                    return true;
                */

                if (PointIdentity(VtxA.VTX, VtxB.VTX))
                    throw new ApplicationException();

                var E2 = obj as VoEdge;
                if (PointIdentity(VtxA.VTX, E2.VtxA.VTX) && PointIdentity(VtxB.VTX, E2.VtxB.VTX))
                    return true;
                if (PointIdentity(VtxA.VTX, E2.VtxB.VTX) && PointIdentity(VtxB.VTX, E2.VtxA.VTX))
                    return true;

                return false;
            }


            public override int GetHashCode() {
                return (int)((VtxA.VTX - VtxB.VTX).AbsSquare()*12345.0);
            }


            public AffineManifold plane {
                get {
                    return AffineManifold.FromPoints(VtxA.VTX, VtxB.VTX);
                }
            }

            public Vector Dir {
                get {
                    return VtxB.VTX - VtxA.VTX;
                }
            }

            public void Flip() {
                var t = VtxA;
                VtxA = VtxB;
                VtxB = t;
            }

            public double GetCoord(Vector V) {
                if (PointIdentity(V, VtxA.VTX))
                    return 0.0;
                if (PointIdentity(V, VtxB.VTX))
                    return 1.0;

                var prjV = plane.ProjectPoint(V);
                return CoordOnLine(prjV, VtxA.VTX, VtxB.VTX);
            }

            public Vector Interpol(double alpha) {
                return (1 - alpha) * (VtxA.VTX) + alpha * (VtxB.VTX);
            }

            public bool Intersect(VoEdge other, out double beta, out Vector I) {
                Vector Dthis = this.Dir;
                Dthis.Normalize();
                Vector Dother = other.Dir;
                Dother.Normalize();
                if(Dother*Dthis < 0) {
                    Dother.Scale(-1.0);
                }
                if(PointIdentity(Dother,Dthis)) {
                    beta = double.PositiveInfinity;
                    I = new Vector(double.PositiveInfinity, double.PositiveInfinity);
                    return false; // parallel - this method is inappropriate
                }


                bool nonParallel = PolygonClipping.ComputeIntersection(VtxA.VTX, VtxB.VTX, other.VtxA.VTX, other.VtxB.VTX, 
                    out double alpha, out beta, out I);

                if (nonParallel == false)
                    return false;

                if (PointIdentity(I, VtxA.VTX)) {
                    alpha = 0.0;
                    I = VtxA.VTX;
                }
                if (PointIdentity(I, VtxB.VTX)) {
                    alpha = 1.0;
                    I = VtxB.VTX;
                }

                if (PointIdentity(I, other.VtxA.VTX)) {
                    beta = 0.0;
                    I = other.VtxA.VTX;
                }
                if (PointIdentity(I, other.VtxB.VTX)) {
                    beta = 1.0;
                    I = other.VtxB.VTX;
                }
                
                if (alpha < 0)
                    return false;
                if (alpha > 1)
                    return false;
                if (beta < 0)
                    return false;
                if (beta > 1)
                    return false;


                return true;
            }
        }

        class VoPolygon : VoItem {
            public List<VoEdge> Edges = new List<VoEdge>();
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

        static bool IsFarPoint(Vector V) {
            if (double.IsInfinity(V.x))
                return true;
            if (double.IsInfinity(V.x))
                return true;

            return false;
        }

        static Func<Vector, Vector, bool> PointIdentity;


        static public AggregationGrid FromPolygonalDomain(MultidimensionalArray Nodes, Vector[] PolygonBoundary, Func<Vector, bool> IsIn, Func<Vector, Vector, bool> __PointIdentity) {

            // check arguments
            // ===============

            PointIdentity = __PointIdentity;

            if (Nodes.Dimension != 2)
                throw new ArgumentException("expecting 2D array;");
            if (Nodes.GetLength(1) != 2)
                throw new ArgumentException("only implemented for 2D");

            foreach (var V in PolygonBoundary) {
                if (V.Dim != 2)
                    throw new ArgumentException();
                if (!IsIn(V)) {
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
            for (int i = 0; i < Boundaries.Length; i++) {
                Vector O = PolygonBoundary[i];
                Vector E = PolygonBoundary[(i + 1) % Boundaries.Length];
                Vector OE = E - O;

                Vector N = new Vector(-OE.y, OE.x);
                N.Normalize();

                Boundaries[i] = new AffineManifold(N, O);

                double eps = Math.Sqrt(BLAS.MachineEps) * 10;
                Vector Cen = O + 0.5 * OE;
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
            List<Vector> Verts;
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
                Verts = new List<Vector>();
                for(int iVtx = 0; iVtx < VertexCoordinates.NoOfRows; iVtx++) {
                    Verts.Add(VertexCoordinates.GetRowPt(iVtx));
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
                    if (iVtxS.Any(i => IsFarPoint(Verts[i])))
                        continue;

                    Vector[] vectors = iVtxS.Select(i => Verts[i]).ToArray();
                    FixOrientation(ref vectors, ref iVtxS);
                    VocellVertexIndex[jV] = iVtxS;
                }
            }


            // build data structures
            // =====================
            List<VoVertex> verticeS = new List<VoVertex>();
            for(int i = 0; i < Verts.Count; i++) {
                verticeS.Add(new VoVertex() {
                    //Index = i,
                    VTX = Verts[i]
                });
            }

            List<VoPolygon> cellS = new List<VoPolygon>();
            List<VoEdge> edgeS = new List<VoEdge>();
            {
                

                bool anyDouble = false;
                for (int jV = 0; jV < VocellVertexIndex.Length; jV++) {
                    var cell = new VoPolygon() {
                        //Index = jV
                    };
                    
                    int[] iVtxS = VocellVertexIndex[jV];

                    int I = iVtxS.Length;
                    for (int i = 0; i < I; i++) {
                        int _iVtxA = iVtxS[i];
                        int _iVtxB = iVtxS[(i + 1) % I];
                        
                        VoVertex _VtxA = verticeS[_iVtxA];
                        VoVertex _VtxB = verticeS[_iVtxB];

                        if (_VtxA.IsFar)
                            continue;
                        if (_VtxB.IsFar)
                            continue;

                        VoEdge newEdge = new VoEdge() {
                            VtxA = _VtxA,
                            VtxB = _VtxB
                        };

                        int iEdge = edgeS.IndexOf(newEdge);
                        if(iEdge < 0) {
                            for(int e = 0; e < edgeS.Count; e++) {
                                Debug.Assert(edgeS[e].Equals(newEdge) == false);
                            }

                            edgeS.Add(newEdge);
                            iEdge = edgeS.Count - 1;
                            //newEdge.Index = iEdge;
                        } else {
                            anyDouble = true;
                            if (edgeS[iEdge].Cells.Count > 1)
                                throw new ArithmeticException("edge shared by more than two cells");
                            newEdge = edgeS[iEdge];
                        }
                        Debug.Assert(newEdge.Equals(edgeS[iEdge]));
                        Debug.Assert(edgeS.Where(ve => ve.Equals(newEdge)).Count() == 1);

                        newEdge.Cells.Add(cell);
                        cell.Edges.Add(newEdge);
                    }
                }

                if (!anyDouble)
                    throw new ArithmeticException("Voronoi diagram seems to be completely disjoint - no edge is used by at least two cells.");
            }

            // Add boundary edges 
            // ==================
            {
                List<int> idxBndyVertices = new List<int>();

                for (int i = 0; i < PolygonBoundary.Length; i++) {
                    Vector bVtx = PolygonBoundary[i];

                    int idxFound = -1;
                    for (int j = 0; j < verticeS.Count; j++) {
                        if (PointIdentity(bVtx, verticeS[j].VTX)) {
                            if (idxFound >= 0) {
                                throw new ArithmeticException();
                            }
                        }
                    }

                    if (idxFound < 0) {
                        var newVtx = new VoVertex() {
                            //Index = verticeS.Count,
                            VTX = bVtx
                        };


                        verticeS.Add(newVtx);
                        idxFound = verticeS.Count - 1;
                    }

                    idxBndyVertices.Add(idxFound);
                }


                for(int i = 0; i < PolygonBoundary.Length; i++) {
                    int i0 = i;
                    int iE = (i + 1) % PolygonBoundary.Length;

                    VoVertex V0 = verticeS[idxBndyVertices[i0]];
                    VoVertex VE = verticeS[idxBndyVertices[iE]];
                    Debug.Assert(PointIdentity(V0.VTX, PolygonBoundary[i0]));
                    Debug.Assert(PointIdentity(VE.VTX, PolygonBoundary[iE]));

                    var newEdge = new VoEdge() {
                        //Index = edgeS.Count,
                        VtxA = V0,
                        VtxB = VE,
                        isBoundary = true
                    };

                    edgeS.Add(newEdge);
                }
            }

            // compute intersections
            // =====================
            {
                List<VoEdge> bndyEdges = edgeS.Where(edge => edge.isBoundary).ToList();

                for(int iBndy = 0; iBndy < bndyEdges.Count; iBndy++) {
                    for(int iEdge = 0; iEdge < edgeS.Count; iEdge++) {
                        VoEdge bndy = bndyEdges[iBndy];
                        VoEdge edge = edgeS[iEdge];
                        if(edge.isBoundary) {
                            continue;
                        }
                        Debug.Assert(bndy.isBoundary == true);

                        // Colinear: exact overlap
                        // -----------------------

                        if(edge.Equals(bndy)) {
                            // +++++++++++++++++++++++++++++++++++++++++++
                            // special case: edge and bndy overlap exactly
                            // => remove edge
                            // +++++++++++++++++++++++++++++++++++++++++++

                            bndy.Cells.AddRange(edge.Cells);
                            foreach(var cell in edge.Cells) {
                                Debug.Assert(cell.Edges.Contains(edge));
                                cell.Edges.Remove(edge);
                                cell.Edges.Add(bndy);
                            }

                            edgeS.RemoveAt(iEdge);

                            edge.deleted = true;
                            iEdge--;
                            continue;
                        }

                        // Colinear: partial overlap
                        // -------------------------
                        var vaProj = bndy.plane.ProjectPoint(edge.VtxA.VTX);
                        bool vAinPlane = PointIdentity(vaProj, edge.VtxA.VTX);

                        var vbProj = bndy.plane.ProjectPoint(edge.VtxB.VTX);
                        bool vBinPlane = PointIdentity(vbProj, edge.VtxB.VTX);
                        {
                            if(vAinPlane && vBinPlane) {
                                // +++++++++++++++++++
                                // edges are co-linear
                                // +++++++++++++++++++

                                if (edge.Dir * bndy.Dir < 0) {
                                    edge.Flip();
                                    Vector t = vaProj;
                                    vaProj = vbProj;
                                    vbProj = t;
                                }
                                Debug.Assert(edge.Dir * bndy.Dir > 0);

                                double alphaA = bndy.GetCoord(vaProj);
                                double alphaB = bndy.GetCoord(vbProj);
                                Debug.Assert(alphaA < alphaB);
                                Debug.Assert((alphaA == 0 && alphaB == 1) == false);

                                if (alphaB <= 0.0)
                                    // no overlap
                                    continue;

                                if(alphaA >= 1.0)
                                    // no overlap
                                    continue;


                                if(alphaA == 0 && alphaB < 1) {
                                    // edge:  o---o
                                    // bndy:  o--------o
                                    //
                                    // out:   o---o----o

                                    edge.isBoundary = true;
                                    if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                        bndyEdges.Add(edge);

                                    Debug.Assert(PointIdentity(bndy.Interpol(alphaB), edge.VtxB.VTX));

                                    bndy.VtxA = edge.VtxB;

                                    iEdge--;
                                    continue;
                                }

                                if(alphaB == 0 && alphaA > 0) {
                                    // edge:       o---o
                                    // bndy:  o--------o
                                    //
                                    // out:   o----o---o

                                    edge.isBoundary = true;
                                    if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                        bndyEdges.Add(edge);

                                    Debug.Assert(PointIdentity(bndy.Interpol(alphaA), edge.VtxA.VTX));

                                    bndy.VtxB = edge.VtxA;

                                    iEdge--;
                                    continue;
                                }

                                if(alphaA == 0 && alphaB > 1) {
                                    // edge:  o------------o
                                    // bndy:  o--------o
                                    //
                                    // out:   o--------o---o

                                    edge.VtxA = bndy.VtxB;
                                    edge.isBoundary = false; // outside domain
                                    
                                    iEdge--;
                                    continue;
                                }

                                if(alphaA < 0 && alphaB == 1) {
                                    // edge:  o------------o
                                    // bndy:      o--------o
                                    //
                                    // out:   o---o--------o

                                    edge.VtxB = bndy.VtxA;
                                    edge.isBoundary = false; // outside domain
                                    
                                    iEdge--;
                                    continue;
                                }

                                Debug.Assert(alphaA != 0.0); // all special cases should treated by now
                                Debug.Assert(alphaB != 1.0); // all special cases should treated by now

                                if(alphaA < 0 && alphaB > 0) {
                                    // edge:  o-----------------o
                                    // bndy:      o--------o
                                    //
                                    // out:   o---o--------o~~~~o

                                    var t = edge.VtxB;
                                    edge.VtxB = bndy.VtxA;
                                    
                                    var newEdge = new VoEdge() {
                                        VtxA = bndy.VtxB,
                                        VtxB = t
                                    };
                                    newEdge.Cells.AddRange(edge.Cells);
                                    edgeS.Add(newEdge);

                                    bndy.Cells.AddRange(edge.Cells);

                                    iEdge--;
                                    continue;
                                }


                                if(alphaA > 0 && alphaB < 1) {
                                    // edge:      o-o
                                    // bndy:   o--------o
                                    //
                                    // out:    o--o-o~~~o
                                    Debug.Assert(alphaB > 0);

                                    var t = bndy.VtxB;
                                    bndy.VtxB = edge.VtxA;

                                    edge.isBoundary = true;
                                    bndyEdges.Add(edge);

                                    var newBndy = new VoEdge() {
                                        isBoundary = true,
                                        VtxA = edge.VtxB,
                                        VtxB = t
                                    };
                                    bndyEdges.Add(newBndy);
                                    edgeS.Add(newBndy);

                                    iEdge--;
                                    continue;
                                }

                                if(alphaA < 0 && alphaB < 1) {
                                    // edge:   o-----o
                                    // bndy:      o--------o
                                    //
                                    // out:    o~~o--o-----o

                                    var newEdge = new VoEdge() {
                                        VtxA = edge.VtxA,
                                        VtxB = bndy.VtxA
                                    };
                                    newEdge.Cells.AddRange(edge.Cells);
                                    edgeS.Add(newEdge);

                                    edge.VtxA = bndy.VtxA;
                                    edge.isBoundary = true;
                                    if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                        bndyEdges.Add(edge);

                                    bndy.VtxA = edge.VtxB;

                                    iEdge--;
                                    continue;
                                }

                                if(alphaA > 0 && alphaB > 1) {
                                    // edge:        o-----o
                                    // bndy:   o--------o
                                    //
                                    // out:    o----o---o~o

                                    var t = bndy.VtxB;
                                    bndy.VtxB = edge.VtxA;

                                    edge.isBoundary = true;
                                    var tt = edge.VtxB;
                                    edge.VtxB = t;
                                    if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                        bndyEdges.Add(edge);


                                    var newEdge = new VoEdge() {
                                        VtxA = t,
                                        VtxB = tt
                                    };
                                    newEdge.Cells.AddRange(edge.Cells);
                                    edgeS.Add(newEdge);
                                    
                                    iEdge--;
                                    continue;
                                }

                                Debug.Assert(false); // this point should never be reached
                            }
                        }

                        // regular intersections
                        // =====================

                        if (PointIdentity(bndy.VtxA.VTX, edge.VtxA.VTX)) continue; 
                        if (PointIdentity(bndy.VtxA.VTX, edge.VtxB.VTX)) continue; 
                        if (PointIdentity(bndy.VtxB.VTX, edge.VtxA.VTX)) continue; 
                        if (PointIdentity(bndy.VtxB.VTX, edge.VtxB.VTX)) continue; 


                        bool cutfound = bndy.Intersect(edge, out double alpha, out Vector I);

                        if(cutfound) {
                            if (alpha == 0.0 || alpha == 1.0) {
                                // T-junction
                                
                                // find junction vertex
                                VoVertex newVert;
                                if (alpha == 0.0) {
                                    Debug.Assert(vAinPlane);
                                    newVert = edge.VtxA;
                                } else {
                                    Debug.Assert(alpha == 1.0);
                                    Debug.Assert(vBinPlane);
                                    newVert = edge.VtxB;
                                }

                                // split bndy -- 1st part:
                                var t = bndy.VtxB;
                                bndy.VtxB = newVert;

                                // split bndy -- 2nd part:
                                var newBndy = new VoEdge() {
                                    isBoundary = true,
                                    VtxA = newVert,
                                    VtxB = t
                                };
                                newBndy.Cells.AddRange(bndy.Cells);
                                bndyEdges.Add(newBndy);
                                edgeS.Add(newBndy);

                                //
                                continue;

                            } else {
                                // X-junction
                                
                                // introduce new vertex
                                var newVert = new VoVertex() {
                                    VTX = I
                                };
                                Debug.Assert(verticeS.Where(vtx => PointIdentity(vtx.VTX, I)).Count() == 0);
                                verticeS.Add(newVert);
                                
                                // split bndy -- 1st part:
                                var t = bndy.VtxB;
                                bndy.VtxB = newVert;

                                // split bndy -- 2nd part:
                                var newBndy = new VoEdge() {
                                    isBoundary = true,
                                    VtxA = newVert,
                                    VtxB = t
                                };
                                newBndy.Cells.AddRange(bndy.Cells);
                                bndyEdges.Add(newBndy);
                                edgeS.Add(newBndy);

                                // split edge -- 1st part:
                                var tt = edge.VtxB;
                                edge.VtxB = newVert;

                                // split edge -- 1st part:
                                var newEdge = new VoEdge() {
                                    VtxA = newVert,
                                    VtxB = tt
                                };
                                newEdge.Cells.AddRange(edge.Cells);
                                edgeS.Add(newEdge);

                                //
                                iEdge--;
                                continue;
                            }
                            
                        }

                    }
                }

            }

            using (var gp = new Gnuplot()) {

                PlotFormat nrmF = new PlotFormat("-xb");
                PlotFormat bdyF = new PlotFormat("-or");

                foreach(var e in edgeS) {
                    var F = !e.isBoundary ? nrmF : bdyF;

                    double[] xS = new[] { e.VtxA.VTX.x, e.VtxB.VTX.x };
                    double[] yS = new[] { e.VtxA.VTX.y, e.VtxB.VTX.y };

                    gp.PlotXY(xS, yS, format:F);
                }



                gp.Execute();
                Console.ReadKey();
            }


            return null;
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
