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

        enum VertexType {
            unspecified = 0,
            
            Inside = 1,

            Outside = 2,

            Boundary = 3,



            FarPoint = 6 // somewhere at infty
        }

        class VoItem {
            //public int Index;
        }

        class VoVertex : VoItem {
            public Vector VTX {
                get {
                    return new Vector(x, y);
                }
            }

            double x;
            double y;

            public bool deleted = false;
            public VertexType type = VertexType.unspecified;

            public List<VoEdge> Edges;

            public bool IsFar {
                get {
                    return IsFarPoint(VTX);
                }
            }
            
            public int ID {
                get;
                private set;
            }

            static int IDcounter = 1;

            public VoVertex(Vector __VTX) {
                if (__VTX.Dim != 2)
                    throw new ArgumentException();
                x = __VTX.x;
                y = __VTX.y;

                ID = IDcounter;
                IDcounter++;
            }

            public override string ToString() {
                return (ID + ": " + VTX.ToString());
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

            public bool isBoundary {
                get {
                    return (VtxA.type == VertexType.Boundary) && (VtxB.type == VertexType.Boundary);
                }
                set {
                    if (value == false)
                        throw new NotSupportedException("unsetting not supported");
                    VtxA.type = VertexType.Boundary;
                    VtxB.type = VertexType.Boundary;
                }
            }


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

                if (PointIdentity(VtxA, VtxB))
                    throw new ApplicationException();

                var E2 = obj as VoEdge;
                if (PointIdentity(VtxA, E2.VtxA) && PointIdentity(VtxB, E2.VtxB))
                    return true;
                if (PointIdentity(VtxA, E2.VtxB) && PointIdentity(VtxB, E2.VtxA))
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
                if (PointIdentityG(V, VtxA.VTX))
                    return 0.0;
                if (PointIdentityG(V, VtxB.VTX))
                    return 1.0;

                var prjV = plane.ProjectPoint(V);
                return CoordOnLine(prjV, VtxA.VTX, VtxB.VTX);
            }

            public Vector Interpol(double alpha) {
                return (1 - alpha) * (VtxA.VTX) + alpha * (VtxB.VTX);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="other"></param>
            /// <param name="alpha">
            /// coordinate of <paramref name="I"/> on this edge
            /// </param>
            /// <param name="beta">
            /// coordinate of <paramref name="I"/> on <paramref name="other"/>
            /// </param>
            /// <param name="I"></param>
            /// <returns></returns>
            public bool Intersect(VoEdge other, out double alpha, out double beta, out Vector I) {
                Vector Dthis = this.Dir;
                Dthis.Normalize();
                Vector Dother = other.Dir;
                Dother.Normalize();
                if(Dother*Dthis < 0) {
                    Dother.Scale(-1.0);
                }
                if(PointIdentityG(Dother,Dthis)) {
                    alpha = double.PositiveInfinity;
                    beta = double.PositiveInfinity;
                    I = new Vector(double.PositiveInfinity, double.PositiveInfinity);
                    return false; // parallel - this method is inappropriate
                }


                bool nonParallel = PolygonClipping.ComputeIntersection(VtxA.VTX, VtxB.VTX, other.VtxA.VTX, other.VtxB.VTX, 
                    out alpha, out beta, out I);

                if (nonParallel == false)
                    return false;

                if (PointIdentityG(I, VtxA.VTX)) {
                    alpha = 0.0;
                    I = VtxA.VTX;
                }
                if (PointIdentityG(I, VtxB.VTX)) {
                    alpha = 1.0;
                    I = VtxB.VTX;
                }

                if (PointIdentityG(I, other.VtxA.VTX)) {
                    beta = 0.0;
                    I = other.VtxA.VTX;
                }
                if (PointIdentityG(I, other.VtxB.VTX)) {
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

            public override string ToString() {
                if (VtxA == null || VtxB == null)
                    return "undefined";
                //return VtxA.VTX.ToString() + "--" + VtxB.VTX.ToString();
                return VtxA.ID  + "--" + VtxB.ID;
            }
        }

        class VoPolygon : VoItem {
            public List<VoEdge> Edges = new List<VoEdge>();

            /// <summary>
            /// true when all vertices are <see cref="VertexType.Inside"/>
            /// </summary>
            public bool IsInside {
                get {
                    bool r = true;

                    foreach(var e in Edges) {
                        if(e.VtxA.type != VertexType.Inside) {
                            r = false;
                            break;
                        }
                        if(e.VtxB.type != VertexType.Inside) {
                            r = false;
                            break;
                        }
                    }
                    
                    return r;
                }
            }

            /// <summary>
            /// true when all vertices are <see cref="VertexType.Outside"/>
            /// </summary>
            public bool IsOutside {
                get {
                    bool r = true;

                    foreach(var e in Edges) {
                        if(e.VtxA.type != VertexType.Outside) {
                            r = false;
                            break;
                        }
                        if(e.VtxB.type != VertexType.Outside) {
                            r = false;
                            break;
                        }
                    }

                    return r;
                }
            }

            /// <summary>
            /// Removes all outside edges from <see cref="Edges"/>
            /// </summary>
            /// <param name="AlsoBoundary">
            /// - true: also edges with <see cref="VoEdge.isBoundary"/>==true will be removed
            /// - false: boundary edges remain in the polygon
            /// </param>
            public void RemoveOutsideParts(bool AlsoBoundary) {
                for(int ie = 0; ie < Edges.Count; ie++) {
                    VoEdge e = Edges[ie];

                    if (e.VtxA.type == VertexType.unspecified) {
                        throw new NotSupportedException();
                    }

                    if(   (e.VtxA.type == VertexType.Outside)
                       || (e.VtxB.type == VertexType.Outside)
                       || (e.VtxA.type == VertexType.FarPoint)
                       || (e.VtxB.type == VertexType.FarPoint)
                       || e.VtxA.IsFar
                       || e.VtxB.IsFar
                       || (AlsoBoundary && e.isBoundary)
                        ) {
                        Edges.RemoveAt(ie);
                        ie--;
                        continue;
                    }

                }
            }

            class MyTuple { // i define my own class, because i need to change values
                //             System.Tuple is read-only; System.ValueTuple is a value type (struct), but algorithm requires ref-type (class)
                public MyTuple(VoEdge e, bool s) { edge = e; used = s; }
                public VoEdge edge;
                public bool used;
            }

            public VoVertex[] GetVerticesSequence(out bool isClosed, bool includeBoundaries = false) {
                List<VoVertex> R = new List<VoVertex>();
                isClosed = false;
                
                MyTuple[] EdgesLoc = this.Edges
                    .Where(edg => includeBoundaries || edg.isBoundary == false)
                    .Select(edg => new MyTuple(edg, false))
                    .ToArray();
                
                MyTuple CurrEdge = EdgesLoc[0]; // select any edge
                CurrEdge.used = true;
                VoVertex start = CurrEdge.edge.VtxA;
                R.Add(CurrEdge.edge.VtxA);
                bool SearchReverse = false;
                while (true) {
                    VoVertex NextVtx = PointIdentity(R.Last(), CurrEdge.edge.VtxA) ? CurrEdge.edge.VtxB : CurrEdge.edge.VtxA;
                    if (PointIdentity(NextVtx, start)) {
                        isClosed = true;
                        break;
                    }
                    R.Add(NextVtx);
                    MyTuple[] NextEdgeS = EdgesLoc.Where(tt =>
                        tt.used == false
                     && (PointIdentity(tt.edge.VtxA, NextVtx) || PointIdentity(tt.edge.VtxB, NextVtx)) 
                     && (!object.ReferenceEquals(tt.edge, CurrEdge))
                    ).ToArray();

                    if (NextEdgeS.Length <= 0) {
                        SearchReverse = true;
                        break;
                    }
                    if (NextEdgeS.Length > 1)
                        throw new NotSupportedException("polygon with branches");

                    CurrEdge = NextEdgeS[0];
                    CurrEdge.used = true;
                }

                if(SearchReverse) {
                    Debug.Assert(isClosed == false);

                    while (true) {
                        VoVertex CurrVertex = R[0];
                        MyTuple[] NextEdgeS = EdgesLoc.Where(tt =>
                            tt.used == false
                         && (PointIdentity(tt.edge.VtxA, CurrVertex) || PointIdentity(tt.edge.VtxB, CurrVertex))
                         && (!object.ReferenceEquals(tt.edge, CurrEdge))
                        ).ToArray();

                        if (NextEdgeS.Length <= 0) {
                            break;
                        }
                        if (NextEdgeS.Length > 1)
                            throw new NotSupportedException("polygon with branches");

                        MyTuple NextEdge = NextEdgeS[0];
                        Debug.Assert(NextEdge.used == false);
                        VoVertex NextVertex;
                        if(PointIdentity(NextEdge.edge.VtxA, CurrVertex)) {
                            NextVertex = NextEdge.edge.VtxB;
                        } else {
                            Debug.Assert(PointIdentity(NextEdge.edge.VtxB, CurrVertex));
                            NextVertex = NextEdge.edge.VtxA;
                        }
                        NextEdge.used = true;
                        R.Insert(0, NextVertex);
                    }
                }

                                             
                return R.ToArray();
            }


            public VoVertex[] GetSegmentBetween(VoVertex Start, VoVertex End) {
                if (PointIdentity(Start, End))
                    throw new ArgumentException();

                VoVertex[] all = GetVerticesSequence(out bool isClosed, true);
                if (!isClosed)
                    throw new NotSupportedException();



                //if (!all.Contains(Start))
                //    throw new ArgumentException("start vertex is not an element of this polygon.");
                //if (!all.Contains(End))
                //    throw new ArgumentException("end vertex is not an element of this polygon.");

                // 
                List<VoVertex>[] Segments = new List<VoVertex>[2];
                double[] Length = new double[2];
                for(int iSweep = 0; iSweep < 2; iSweep++) { // sweep 1: search forward; step 2: search backward;
                    Segments[iSweep] = new List<VoVertex>();

                    int i0 = all.IndexOf(Start, PointIdentity);
                    if(i0 < 0) 
                        throw new ArgumentException("Start vertex is not an element of this polygon.", "Start");

                    bool bFound = false;
                    for(int i = 0; i < all.Length; i++) {
                        int iCnt = (i0 + i) % all.Length;
                        Segments[iSweep].Add(all[iCnt]);

                        if(i > 0) {
                            Length[iSweep] += Segments[iSweep][i].VTX.Dist(Segments[iSweep][i - 1].VTX);
                        }

                        if(PointIdentity(all[iCnt], End)) {
                            bFound = true;
                            break;
                        }
                    }

                    if(!bFound)
                        throw new ArgumentException("End vertex is not an element of this polygon.", "End");

                    // check
                    Debug.Assert(PointIdentity(Start, Segments[iSweep].First()));
                    Debug.Assert(PointIdentity(End, Segments[iSweep].Last()));
                    
                    // reverse list: 
                    all = all.Reverse().ToArray();
                }

                // return shorter segment
                if(Length[0] < Length[1]) {
                    return Segments[0].ToArray();
                } else {
                    return Segments[1].ToArray();
                }

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

        static bool IsFarPoint(Vector V) {
            if (double.IsInfinity(V.x))
                return true;
            if (double.IsInfinity(V.x))
                return true;

            return false;
        }

        static Func<Vector, Vector, bool> PointIdentityG;

        static bool PointIdentity(VoVertex V1, VoVertex V2) {
            if(object.ReferenceEquals(V1,V2)) {
                Debug.Assert(PointIdentityG(V1.VTX, V2.VTX));
                return true;
            } else {
                return PointIdentityG(V1.VTX, V2.VTX);
            }


        }


        static public AggregationGrid FromPolygonalDomain(MultidimensionalArray Nodes, Vector[] PolygonBoundary, bool mirroring, Func<Vector, bool> IsIn, Func<Vector, Vector, bool> __PointIdentity) {

            // check arguments
            // ===============

            PointIdentityG = __PointIdentity;

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
                    if (PointIdentityG(PolygonBoundary[iBnd], PolygonBoundary[iBnd2]))
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
                if (PointIdentityG(Cen, CenProj) == false)
                    throw new ArithmeticException("point identity does not seem to work");

                Debug.Assert(PointIdentityG(Cen, Boundaries[i].ProjectPoint(Cen + N * OE.Abs()))); // tests if the 'ProjectPoint' works correctly

            }


            // Vertex Mirroring
            // ================
            MultidimensionalArray NewNodes;
            if(mirroring) {
                
                double[] xNodes = Nodes.GetColumn(0);
                double[] yNodes = Nodes.GetColumn(1);

                Mirror(ref xNodes, ref yNodes, Boundaries, IsIn);
                Debug.Assert(xNodes.Length == yNodes.Length);

                NewNodes = MultidimensionalArray.Create(xNodes.Length, 2);
                NewNodes.SetColumn(0, xNodes);
                NewNodes.SetColumn(1, yNodes);
                

                
                Nodes = null;
            } else {
                NewNodes = Nodes.CloneAs();
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
                    //Vector[] _vectors = vectors;
                    int sign = CheckOrientation(vectors);
                    if(sign > 0) {
                        // nop
                    } else if(sign < 0) {
                        iVtxS = iVtxS.Reverse().ToArray();
                    } else {
                        throw new ArithmeticException("got indefinite polygon form matlab");
                    }
                    VocellVertexIndex[jV] = iVtxS;
                }
            }


            // build data structures
            // =====================
            List<VoVertex> verticeS = new List<VoVertex>();
            for(int i = 0; i < Verts.Count; i++) {
                verticeS.Add(new VoVertex(Verts[i]));
            }

            List<VoPolygon> cellS = new List<VoPolygon>();
            List<VoEdge> edgeS = new List<VoEdge>();
            {
                

                bool anyDouble = false;
                for (int jV = 0; jV < VocellVertexIndex.Length; jV++) {
                    var cell = new VoPolygon() {
                        //Index = jV
                    };
                    cellS.Add(cell);
                    
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

                    //if (cell.Edges.Count <= 1)
                    //    Console.WriteLine("warn0: degen poly");
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
                        if (PointIdentityG(bVtx, verticeS[j].VTX)) {
                            if (idxFound >= 0) {
                                throw new ArithmeticException();
                            }
                        }
                    }

                    if (idxFound < 0) {
                        var newVtx = new VoVertex(bVtx);


                        verticeS.Add(newVtx);
                        idxFound = verticeS.Count - 1;
                    }

                    idxBndyVertices.Add(idxFound);
                    verticeS[idxFound].type = VertexType.Boundary;
                }


                for(int i = 0; i < PolygonBoundary.Length; i++) {
                    int i0 = i;
                    int iE = (i + 1) % PolygonBoundary.Length;

                    VoVertex V0 = verticeS[idxBndyVertices[i0]];
                    VoVertex VE = verticeS[idxBndyVertices[iE]];
                    Debug.Assert(PointIdentityG(V0.VTX, PolygonBoundary[i0]));
                    Debug.Assert(PointIdentityG(VE.VTX, PolygonBoundary[iE]));

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
            CheckEdgeUniqueness(edgeS, true);
            {
                List<VoEdge> bndyEdges = edgeS.Where(edge => edge.isBoundary).ToList();

                int iRun = 0;


                for (int iBndy = 0; iBndy < bndyEdges.Count; iBndy++) {
                    for (int iEdge = 0; iEdge < edgeS.Count; iEdge++) {
                        iRun++;
                        VoEdge bndy = bndyEdges[iBndy];
                        VoEdge edge = edgeS[iEdge];
                        CheckEdgeUniqueness(edgeS);
                        if (edge.isBoundary) {
                            continue;
                        }
                        Debug.Assert(bndy.isBoundary == true);
                        Debug.Assert(bndy.VtxA.type == VertexType.Boundary);
                        Debug.Assert(bndy.VtxB.type == VertexType.Boundary);

                        // Colinear: exact overlap
                        // - - - - - - - - - - - - 

                        if (edge.Equals(bndy)) {
                            // +++++++++++++++++++++++++++++++++++++++++++
                            // special case: edge and bndy overlap exactly
                            // => remove edge
                            // +++++++++++++++++++++++++++++++++++++++++++

                            bndy.Cells.AddRange(edge.Cells);
                            //foreach (var cell in edge.Cells) {
                            //    Debug.Assert(cell.Edges.Contains(edge));
                            //    cell.Edges.Remove(edge);
                            //    cell.Edges.Add(bndy);
                            //}

                            edgeS.RemoveAt(iEdge);

                            //edge.deleted = true;
                            iEdge--;
                            continue;
                        }

                        // Colinear: partial overlap
                        // - - - - - - - - - - - - - 
                        var vaProj = bndy.plane.ProjectPoint(edge.VtxA.VTX);
                        bool vAinPlane = PointIdentityG(vaProj, edge.VtxA.VTX);

                        var vbProj = bndy.plane.ProjectPoint(edge.VtxB.VTX);
                        bool vBinPlane = PointIdentityG(vbProj, edge.VtxB.VTX);

                        if (vAinPlane && vBinPlane) {
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

                            if (alphaB <= 0.0) {
                                // no overlap
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            if (alphaA >= 1.0) {
                                // no overlap
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }


                            if (alphaA == 0 && alphaB < 1) {
                                // edge:  o---o
                                // bndy:  o--------o
                                //
                                // out:   o---o----o

                                // pt1 (edge)
                                edge.isBoundary = true;
                                if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                    bndyEdges.Add(edge);
                                Debug.Assert(PointIdentityG(bndy.Interpol(alphaB), edge.VtxB.VTX));
                                CheckEdgeUniqueness(edgeS);

                                // pt2 (bndy)
                                bndy.VtxA = edge.VtxB;
                                Debug.Assert(bndy.isBoundary == true);
                                CheckEdgeUniqueness(edgeS);

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            if (alphaB == 1 && alphaA > 0) {
                                // edge:       o---o
                                // bndy:  o--------o
                                //
                                // out:   o----o---o

                                // pt2 (edge)
                                edge.isBoundary = true;
                                if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                    bndyEdges.Add(edge);
                                Debug.Assert(PointIdentityG(bndy.Interpol(alphaA), edge.VtxA.VTX));
                                CheckEdgeUniqueness(edgeS);

                                // pt1 (bndy)
                                bndy.VtxB = edge.VtxA;
                                bndy.isBoundary = true;
                                CheckEdgeUniqueness(edgeS);

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            if (alphaA == 0 && alphaB > 1) {
                                // edge:  o------------o
                                // bndy:  o--------o
                                //
                                // out:   o--------o---o

                                // pt1 (bndy)
                                bndy.Cells.AddRange(edge.Cells);

                                // pt2 (edge)
                                edge.VtxA = bndy.VtxB;
                                CheckEdgeUniqueness(edgeS);
                                Debug.Assert(edge.isBoundary == false); // outside domain

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            if (alphaA < 0 && alphaB == 1) {
                                // edge:  o------------o
                                // bndy:      o--------o
                                //
                                // out:   o---o--------o

                                // pt1 (edge)
                                edge.VtxB = bndy.VtxA;
                                CheckEdgeUniqueness(edgeS);
                                Debug.Assert(edge.isBoundary == false); // outside domain

                                // pt2 (bndy) 
                                bndy.Cells.AddRange(edge.Cells);

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            Debug.Assert(alphaA != 0.0); // all special cases should treated by now
                            Debug.Assert(alphaB != 1.0); // all special cases should treated by now

                            if (alphaA < 0 && alphaB > 0) {
                                // edge:  o-----------------o
                                // bndy:      o--------o
                                //
                                // out:   o---o--------o~~~~o

                                // pt1 (edge)
                                var t = edge.VtxB;
                                edge.VtxB = bndy.VtxA;
                                Debug.Assert(edge.isBoundary == false); // outside domain

                                // pt2 (bndy)
                                bndy.Cells.AddRange(edge.Cells);
                                Debug.Assert(bndy.isBoundary == true);

                                // pt3 (new)
                                var newEdge = new VoEdge() {
                                    VtxA = bndy.VtxB,
                                    VtxB = t
                                };
                                newEdge.Cells.AddRange(edge.Cells);
                                edgeS.Add(newEdge);
                                CheckEdgeUniqueness(edgeS);
                                Debug.Assert(newEdge.isBoundary == false); // outside domain

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }


                            if (alphaA > 0 && alphaB < 1) {
                                // edge:      o-o
                                // bndy:   o--------o
                                //
                                // out:    o--o-o~~~o
                                Debug.Assert(alphaB > 0);

                                // pt1 (bndy)
                                var t = bndy.VtxB;
                                bndy.VtxB = edge.VtxA;
                                bndy.isBoundary = true;

                                // pt2 (edge)
                                edge.isBoundary = true;
                                bndyEdges.Add(edge);

                                // pt3 (new)
                                var newBndy = new VoEdge() {
                                    VtxA = edge.VtxB,
                                    VtxB = t
                                };
                                newBndy.isBoundary = true;
                                bndyEdges.Add(newBndy);
                                edgeS.Add(newBndy);
                                CheckEdgeUniqueness(edgeS);

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            if (alphaA < 0 && alphaB < 1) {
                                // edge:   o-----o
                                // bndy:      o--------o
                                //
                                // out:    o~~o--o-----o

                                // pt1 (new)
                                var newEdge = new VoEdge() {
                                    VtxA = edge.VtxA,
                                    VtxB = bndy.VtxA
                                };
                                newEdge.Cells.AddRange(edge.Cells);
                                Debug.Assert(newEdge.isBoundary == false); // outside domain
                                edgeS.Add(newEdge);
                                CheckEdgeUniqueness(edgeS);

                                // pt2 (edge)
                                edge.VtxA = bndy.VtxA;
                                edge.isBoundary = true;
                                if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                    bndyEdges.Add(edge);

                                // pt3 (bndy)
                                bndy.VtxA = edge.VtxB;
                                bndy.isBoundary = true;

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }

                            if (alphaA > 0 && alphaB > 1) {
                                // edge:        o-----o
                                // bndy:   o--------o
                                //
                                // out:    o----o---o~o

                                // pt1 (bndy)
                                var t = bndy.VtxB;
                                bndy.VtxB = edge.VtxA;
                                bndy.isBoundary = true;

                                // pt2 (edge)
                                edge.isBoundary = true;
                                var tt = edge.VtxB;
                                edge.VtxB = t;
                                if (!bndyEdges.Contains(edge, (a, b) => object.ReferenceEquals(a, b)))
                                    bndyEdges.Add(edge);

                                // pt3 (new)
                                var newEdge = new VoEdge() {
                                    VtxA = t,
                                    VtxB = tt
                                };
                                newEdge.Cells.AddRange(edge.Cells);
                                Debug.Assert(newEdge.isBoundary == false); // outside domain
                                edgeS.Add(newEdge);
                                CheckEdgeUniqueness(edgeS);

                                // cont
                                iEdge--;
                                CheckEdgeUniqueness(edgeS);
                                continue;
                            }



                        } else {
                            // +++++++++++++++++++++++++
                            // edges are NOT colinear
                            // +++++++++++++++++++++++++

                            if (PointIdentity(bndy.VtxA, edge.VtxA)) { CheckEdgeUniqueness(edgeS); continue; } // L-junction
                            if (PointIdentity(bndy.VtxA, edge.VtxB)) { CheckEdgeUniqueness(edgeS); continue; } // L-junction
                            if (PointIdentity(bndy.VtxB, edge.VtxA)) { CheckEdgeUniqueness(edgeS); continue; } // L-junction
                            if (PointIdentity(bndy.VtxB, edge.VtxB)) { CheckEdgeUniqueness(edgeS); continue; } // L-junction


                            bool cutfound = bndy.Intersect(edge, out double alpha, out double beta, out Vector I);

                            if (cutfound) {
                                //var vaProj = bndy.plane.ProjectPoint(edge.VtxA.VTX);
                                //bool vAinPlane = PointIdentity(vaProj, edge.VtxA.VTX);

                                //var vbProj = bndy.plane.ProjectPoint(edge.VtxB.VTX);
                                //bool vBinPlane = PointIdentity(vbProj, edge.VtxB.VTX);

                                bool isL = false;
                                if ((alpha == 0 || alpha == 1) && (beta == 0.0 || beta == 1.0)) {
                                    // L-junction: nothing to do

                                    isL = true;
                                    //
                                    //continue;
                                }
                                Debug.Assert(isL == false);

                                if (beta == 0.0 || beta == 1.0) {
                                    // T-junction

                                    // find junction vertex
                                    VoVertex newVert;
                                    if (beta == 0.0) {
                                        Debug.Assert(vAinPlane);
                                        newVert = edge.VtxA;
                                    } else {
                                        Debug.Assert(beta == 1.0);
                                        Debug.Assert(vBinPlane);
                                        newVert = edge.VtxB;
                                    }

                                    // split bndy -- 1st part:
                                    var t = bndy.VtxB;
                                    bndy.VtxB = newVert;
                                    bndy.isBoundary = true;
                                    CheckEdgeUniqueness(edgeS);

                                    // split bndy -- 2nd part:
                                    var newBndy = new VoEdge() {
                                        VtxA = newVert,
                                        VtxB = t
                                    };
                                    newBndy.isBoundary = true;
                                    newBndy.Cells.AddRange(bndy.Cells);
                                    bndyEdges.Add(newBndy);
                                    edgeS.Add(newBndy);
                                    CheckEdgeUniqueness(edgeS);

                                    //
                                    CheckEdgeUniqueness(edgeS);
                                    continue;

                                } else {
                                    // X-junction

                                    // introduce new vertex
                                    var newVert = new VoVertex(I);
                                    Debug.Assert(verticeS.Where(vtx => PointIdentityG(vtx.VTX, I)).Count() == 0);
                                    verticeS.Add(newVert);

                                    // split bndy -- 1st part:
                                    var t = bndy.VtxB;
                                    bndy.VtxB = newVert;
                                   
                                    // split bndy -- 2nd part:
                                    var newBndy = new VoEdge() {
                                        VtxA = newVert,
                                        VtxB = t
                                    };
                                    newBndy.isBoundary = true;
                                    newBndy.Cells.AddRange(bndy.Cells);
                                    bndyEdges.Add(newBndy);
                                    edgeS.Add(newBndy);
                                    CheckEdgeUniqueness(edgeS);

                                    // split edge -- 1st part:
                                    var tt = edge.VtxB;
                                    edge.VtxB = newVert;
                                    Debug.Assert(edge.isBoundary == false); // could fail in those cases when an edge exits the domain and enters again
                                    //                                         => store 'isBoundary' separately

                                    // split edge -- 2nd part:
                                    var newEdge = new VoEdge() {
                                        VtxA = newVert,
                                        VtxB = tt
                                    };
                                    newEdge.Cells.AddRange(edge.Cells);
                                    edgeS.Add(newEdge);
                                    CheckEdgeUniqueness(edgeS);
                                    Debug.Assert(edge.isBoundary == false); // could fail in those cases when an edge exits the domain and enters again
                                    //                                         => store 'isBoundary' separately

                                    //
                                    iEdge--;
                                    CheckEdgeUniqueness(edgeS);
                                    continue;
                                }

                            }
                        }
                    }
                }

                // pass 2: finding all T- and X-junctions
                // --------------------------------------

                //for (int iBndy = 0; iBndy < bndyEdges.Count; iBndy++) {
                //    for (int iEdge = 0; iEdge < edgeS.Count; iEdge++) {
                //        iRun++;
                //        VoEdge bndy = bndyEdges[iBndy];
                //        VoEdge edge = edgeS[iEdge];
                //        CheckEdgeUniqueness(edgeS);
                //        if (edge.isBoundary) {
                //            continue;
                //        }
                //        Debug.Assert(bndy.isBoundary == true);
                //        Debug.Assert(bndy.VtxA.type == VertexType.Boundary);
                //        Debug.Assert(bndy.VtxB.type == VertexType.Boundary);


                //    }
                //}

                Debug.Assert(bndyEdges.SetEquals(edgeS.Where(edge => edge.isBoundary)));
                bndyEdges = null; // forget


                // eliminate duplicate boundaries
                // ------------------------------

                for (int e = 0; e < edgeS.Count; e++) {
                    Debug.Assert(PointIdentity(edgeS[e].VtxA, edgeS[e].VtxB) == false);
                    //if (!edgeS[e].isBoundary)
                    //    continue;
                    for (int e2 = e + 1; e2 < edgeS.Count; e2++) {
                        Debug.Assert(object.ReferenceEquals(edgeS[e], edgeS[e2]) == false);
                        //if (!edgeS[e2].isBoundary)
                        //    continue;

                        if(edgeS[e].Equals(edgeS[e2])) {
                            // duplicate found
                            var e2Remove = edgeS[e2];
                            edgeS.RemoveAt(e2);
                            foreach (var cl in e2Remove.Cells) {
                                
                                if (!edgeS[e].Cells.ContainsRefEqual(cl))
                                    edgeS[e].Cells.Add(cl);

                                
                            }
                            e2--;
                        }
                    }
                }

                CheckEdgeUniqueness(edgeS, true);
            }

            // final point classification
            // ==========================

            foreach(var V in verticeS) {
                if(V.type == VertexType.unspecified) {
                    if(IsFarPoint(V.VTX)) {
                        V.type = VertexType.FarPoint;
                        
                    } else if(IsIn(V.VTX)) {
                        V.type = VertexType.Inside;
                    } else {
                        V.type = VertexType.Outside;
                    }
                }
            }

            // repair cell-to-edges, vertices-to-edges
            // =======================================
            {
                for(int jV = 0; jV < cellS.Count; jV++){
                    var cell = cellS[jV];
                    //if (cell.Edges.Count <= 1)
                    //    Console.WriteLine("warn: less than 2");
                    cell.Edges.Clear();
                }

                HashSet<VoPolygon> cellS4edge = new HashSet<VoPolygon>(new FuncEqualityComparer<VoPolygon>((a, b) => object.ReferenceEquals(a, b)));

                for(int ie = 0; ie < edgeS.Count; ie++) {
                    var edge = edgeS[ie];

                    cellS4edge.Clear();
                    cellS4edge.AddRange(edge.Cells);
                    Debug.Assert(cellS4edge.Count == edge.Cells.Count);
                    foreach(var cell in cellS4edge) {
                        Debug.Assert(cellS.ContainsRefEqual(cell));
                    }
                    
                    for (int ic = 0; ic < edge.Cells.Count; ic++) {
                        var cell = edge.Cells[ic];
                        Debug.Assert(cell.Edges.ContainsRefEqual(edge) == false);
                        cell.Edges.Add(edge);
                    }


                    if (edge.VtxA.Edges == null)
                        edge.VtxA.Edges = new List<VoEdge>();
                    if (edge.VtxB.Edges == null)
                        edge.VtxB.Edges = new List<VoEdge>();

                    if (!edge.VtxA.Edges.ContainsRefEqual(edge))
                        edge.VtxA.Edges.Add(edge);
                    if (!edge.VtxB.Edges.ContainsRefEqual(edge))
                        edge.VtxB.Edges.Add(edge);
                }
                

                for(int jV = 0; jV < cellS.Count; jV++){
                    var cell = cellS[jV];
                    //if (cell.Edges.Count <= 1)
                    //    Console.WriteLine("warn2: less than 2");
                }
            }

            // collect boundary polygon
            // ========================
            VoPolygon bndyPoly;
            {
                bndyPoly = new VoPolygon();
                bndyPoly.Edges = edgeS.Where(edge => edge.isBoundary == true).ToList();

                VoVertex[] bndyPolyVertices;
                bndyPolyVertices = bndyPoly.GetVerticesSequence(out bool bndyClosed, true);
                if (!bndyClosed)
                    throw new Exception();
            }


            // collect inside polygons
            // =======================
            List<VoVertex[]> Insiders = new List<VoVertex[]>();
            for(int j = 0; j < cellS.Count; j++) {
                VoPolygon Cj = cellS[j];

                if (Cj.IsOutside)
                    continue;

                int NoEdgesB4 = Cj.Edges.Count;

                bool wasInside = Cj.IsInside;
                if (!wasInside) {
                    Cj.RemoveOutsideParts(true);
                    if (Cj.Edges.Count <= 1)
                        continue; // everything culled away
                }

                int NoEdgesAf = Cj.Edges.Count;

                VoVertex[] seq = Cj.GetVerticesSequence(out bool Closed, false);
                int sign = CheckOrientation(seq);
                if (sign == 0) {
                    DebugPlot(VocellVertexIndex, Verts, Cj.Edges);
                    throw new ArithmeticException("indefinite polygon.");
                }
                if (sign < 0)
                    seq = seq.Reverse().ToArray();
                Debug.Assert(CheckOrientation(seq) > 0);
                //Console.WriteLine("{0}: {1} {2}", j, Closed, sign);

                if(!Closed) {

                    VoVertex[] closing = bndyPoly.GetSegmentBetween(seq.Last(), seq[0]);
                    if(closing.Length > 2) {
                        VoVertex[] closing_inner = closing.GetSubVector(1, closing.Length - 2);

                        seq = seq.Cat(closing_inner);
                    }

                    
                }
                Insiders.Add(seq);
            }


            //DebugPlot(VocellVertexIndex, Verts, edgeS);

            // Build BoSSS structure 
            // =====================
            List<Cell> cells = new List<Cell>();
            List<int[]> aggregation = new List<int[]>();
            {
                for (int jV = 0; jV < Insiders.Count; jV++) { // loop over Voronoi Cells

                    int[] iVtxS = Insiders[jV].Select(voVtx => voVtx.ID).ToArray();
                    int NV = iVtxS.Length;

                    Vector[] VoronoiCell = Insiders[jV].Select(voVtx => voVtx.VTX).ToArray();





                    //FixOrientation(ref VoronoiCell, ref iVtxS);

                    int[,] iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell); ;


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

                }
            }

            // return grid
            // ===========
            {
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
        }


        static void DebugPlot(int[][] VocellVertexIndex, IList<Vector> Verts, IList<VoEdge> edgeS) {
            using (var gp = new Gnuplot()) {
                if (VocellVertexIndex != null) {
                    PlotFormat orgF = new PlotFormat(":k");
                    orgF.LineWidth = 2;
                    foreach (int[] orgCell in VocellVertexIndex) {
                        List<double> xS = new List<double>();
                        List<double> yS = new List<double>();

                        int I = orgCell.Length;
                        bool containsInf = false;
                        for (int i = 0; i < I; i++) {
                            Vector vtx = Verts[orgCell[i]];
                            if (IsFarPoint(vtx)) {
                                containsInf = true;
                            } else {
                                xS.Add(vtx.x);
                                yS.Add(vtx.y);
                            }
                        }

                        if (!containsInf) {
                            xS.Add(xS[0]);
                            yS.Add(yS[0]);
                        }

                        gp.PlotXY(xS, yS, format: orgF);
                    }
                }
                if (edgeS != null) {
                    PlotFormat nrmF = new PlotFormat("-xb");
                    PlotFormat bdyF = new PlotFormat("-or");


                    foreach (var e in edgeS) {
                        var F = !e.isBoundary ? nrmF : bdyF;

                        double[] xS = new[] { e.VtxA.VTX.x, e.VtxB.VTX.x };
                        double[] yS = new[] { e.VtxA.VTX.y, e.VtxB.VTX.y };

                        gp.PlotXY(xS, yS, format: F);
                    }
                }
                gp.Execute();
                Console.WriteLine("Press any key to continue...");
                Console.ReadKey();
            }
        }



        static void CheckEdgeUniqueness(IEnumerable<VoEdge> __Allbndy, bool really = false) {
            VoEdge[] Allbndy = __Allbndy.ToArray();

            //really = true;
            if (really) {
                for (int e = 0; e < Allbndy.Length; e++) {
                    Debug.Assert(PointIdentity(Allbndy[e].VtxA, Allbndy[e].VtxB) == false);
                    for (int e2 = e + 1; e2 < Allbndy.Length; e2++) {
                        Debug.Assert(object.ReferenceEquals(Allbndy[e], Allbndy[e2]) == false);
                        Debug.Assert(Allbndy[e].Equals(Allbndy[e2]) == false);
                    }
                }
            }

            /*
            Vector a = new Vector(0.56035901324048, 0.999999999999999);
            Vector b = new Vector(1, 1);
            for (int e = 0; e < Allbndy.Length; e++) {
                var edg = Allbndy[e];

                if (PointIdentity(edg.VtxA.VTX, a) && PointIdentity(edg.VtxB.VTX, b))
                    Debug.Assert(false);
                if (PointIdentity(edg.VtxA.VTX, b) && PointIdentity(edg.VtxB.VTX, a))
                    Debug.Assert(false);
            }
            */
        }

        static int CheckOrientation(IEnumerable<VoVertex> __Polygon) {
            Vector[] Polygon = __Polygon.Select(vov => vov.VTX).ToArray();
            return CheckOrientation(Polygon);
        }


        public static int CheckOrientation(Vector[] Polygon) {
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
                return 0;
                //throw new ArithmeticException("Indefinite polygon");

            if (AllPos)
                return 1;

            if (AllNeg) {
                //Polygon = Polygon.Reverse().ToArray();
                //iVtx = iVtx.Reverse().ToArray();
                return -1;
            }

            throw new ArithmeticException("Indefinite polygon");
        }

    }
}
