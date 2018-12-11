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

            public static VoVertex Create(Vector __VTX) {
                foreach (var v in verticeS) {
                    if (PointIdentityG(v.VTX, __VTX))
                        return v;
                }
                var vv = new VoVertex(__VTX);
                verticeS.Add(vv);
                return vv;
            }

            public static List<VoVertex> verticeS = new List<VoVertex>();

            private VoVertex(Vector __VTX) {
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


            public override bool Equals(object obj) {
                if (object.ReferenceEquals(this, obj))
                    return true;

                return PointIdentity(this, (VoVertex)obj);
            }

            public override int GetHashCode() {
                return 1;
            }
        }



        class VoEdge : VoItem {
            public static VoEdge Create(VoVertex __VtxA, VoVertex __VtxB) {
                if (__VtxA.Equals(__VtxB))
                    throw new ArgumentException();

                //CheckEdgeUniqueness(true);
                foreach (var e in edgeS) {
                    if ((e.VtxA.Equals(__VtxA) && e.VtxB.Equals(__VtxB))
                       || (e.VtxA.Equals(__VtxB) && e.VtxB.Equals(__VtxA))) {
                        return e;
                    }
                }

                var ee = new VoEdge() { VtxA = __VtxA, VtxB = __VtxB };
                edgeS.Add(ee);
                //CheckEdgeUniqueness(true);
                return ee;

            }

            public int ID {
                get;
                private set;
            }

            static int IDcounter = 1;

            private VoEdge() {
                this.ID = IDcounter;
                IDcounter++;
            }

            internal void CheckCell2Edge() {
                foreach (var cl in this.Cells) {
                    Debug.Assert(cl.Edges.ContainsRefEqual(this));
                }
            }


            public void Split(VoVertex I, out VoEdge Edge1, out VoEdge Edge2) {
                Vector Ip = this.plane.ProjectPoint(I.VTX);
                Debug.Assert(PointIdentityG(I.VTX, Ip));
                CheckCell2Edge();

                bool bndy = this.isBoundary;
                var cellsClone = this.Cells.ToArray();

                // create first part (by modifying this)
                VoVertex t = this.VtxB;
                this.VtxB = I;
                Edge1 = this;
                if (bndy)
                    Edge1.isBoundary = bndy;
                Edge1.CheckCell2Edge();

                //var miselfAnd = edgeS.Where(edg => edg.Equals(this)).ToArray();
                //Debug.Assert(miselfAnd.Length == 1);
                //Debug.Assert(object.ReferenceEquals(miselfAnd[0], this));
                int w = 0;
                for (int i = 0; i < edgeS.Count; i++) {
                    var edg = edgeS[i];
                    if (object.ReferenceEquals(edg, this))
                        continue;

                    if (edg.Equals(this)) {
                        if (edg.isBoundary)
                            this.isBoundary = true;

                        VoPolygon[] otherCells = edg.Cells.ToArray();
                        foreach (var cl in otherCells) {
                            //if (!this.Cells.Contains(cl))
                            //    this.Cells.Add(cl);

                            //cl.RemoveEdge(edg);
                            //cl.AddEdge(this);
                            cl.ReplaceEdge(edg, this);

                            Debug.Assert(this.Cells.ContainsRefEqual(cl));
                        }


                        w++;
                        Debug.Assert(w <= 1);
                        edgeS.RemoveAt(i);
                        i--;
                    }
                }

                // check
                foreach(var cl in Edge1.Cells) {
                    bool s = cl.CheckOriginalDomain();
                    if(!s) {
                        DebugPlot(null, null, cl.Edges, null);
                    }
                    Debug.Assert(s);
                }
               


                // create second part
                Edge2 = Create(I, t);
                foreach (var cl in cellsClone) {
                    cl.AddEdge(Edge2);
                    if (!Edge2.Cells.Contains(cl))
                        Edge2.Cells.Add(cl);
                }
                if (bndy)
                    Edge2.isBoundary = bndy;
                Edge2.CheckCell2Edge();

                // check
                foreach(var cl in Edge1.Cells) {
                    bool s = cl.CheckOriginalDomain();
                    if(!s) {
                        DebugPlot(null, null, cl.Edges, null);
                    }
                    Debug.Assert(s);
                }
                foreach(var cl in Edge2.Cells) {
                    bool s = cl.CheckOriginalDomain();
                    if(!s) {
                        DebugPlot(null, null, cl.Edges, null);
                    }
                    Debug.Assert(s);
                }
            }



            static public List<VoEdge> edgeS = new List<VoEdge>();

            /// <summary>
            /// First vertex in Voronoi cell
            /// </summary>
            public VoVertex VtxA;

            /// <summary>
            /// Second vertex in Voronoi cell
            /// </summary>
            public VoVertex VtxB;


            public List<VoPolygon> Cells = new List<VoPolygon>();

            /// <summary>
            /// checks in DEBUG for any duplicates in <see cref="Cells"/>
            /// </summary>
            public void SanitizeCells() {
#if DEBUG
                for (int i = 0; i < Cells.Count; i++) {
                    for (int j = i + 1; j < Cells.Count; j++) {
                        var Cj = Cells[j];
                        var Ci = Cells[i];

                        if (Ci.Equals(Cj)) {
                            //Cells.RemoveAt(j);
                            Debug.Assert(false);
                            j--;
                        }
                    }
                }
#endif
            }

            bool m_IsBoundary = false;

            public bool isBoundary {
                get {
                    //return (VtxA.type == VertexType.Boundary) && (VtxB.type == VertexType.Boundary);
                    return m_IsBoundary;
                }
                set {
                    if (value == false)
                        throw new NotSupportedException("unsetting not supported");
                    VtxA.type = VertexType.Boundary;
                    VtxB.type = VertexType.Boundary;
                    m_IsBoundary = value;
                }
            }


            public override bool Equals(object obj) {
                if (object.ReferenceEquals(this, obj))
                    return true;

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
                return 0; // VtxA.ID + VtxB.ID; // (int)((VtxA.VTX - VtxB.VTX).AbsSquare()*12345.0);
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

            /// <summary>
            /// returns a 1D coordinate of <paramref name="V"/> on this edge.
            /// </summary>
            public double GetCoord(Vector V) {
                if (PointIdentityG(V, VtxA.VTX))
                    return 0.0;
                if (PointIdentityG(V, VtxB.VTX))
                    return 1.0;

                var prjV = plane.ProjectPoint(V);
                return CoordOnLine(prjV, VtxA.VTX, VtxB.VTX);
            }

            /// <summary>
            /// linear interpolation between <see cref="VtxA"/> (<paramref name="alpha"/> == 0) and <see cref="VtxB"/> (<paramref name="alpha"/> == 1).
            /// </summary>
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
                if (Dother * Dthis < 0) {
                    Dother.Scale(-1.0);
                }
                if (PointIdentityG(Dother, Dthis)) {
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
                return VtxA.ID + "--" + VtxB.ID;
            }


        }

        class VoPolygon : VoItem {

            public int ID {
                get;
                private set;
            }

            AffineManifold[] DomainCheck;

            double hmin;


            public VoPolygon(IEnumerable<VoEdge> edges, bool bConvex) {
                m_Edges.AddRange(edges);
                hmin = double.MaxValue;
                foreach (var e in this.Edges) {
                    if (!e.Cells.ContainsRefEqual(this))
                        e.Cells.Add(this);
                    hmin = Math.Min(hmin, e.Dir.Abs());
                }
                this.ID = IDcounter;
                IDcounter++;


                DomainCheck = new AffineManifold[m_Edges.Count];
                for(int ie = 0; ie < m_Edges.Count; ie++) {
                    VoVertex VA = m_Edges[ie].VtxA;
                    VoVertex VB = m_Edges[ie].VtxB;

                    DomainCheck[ie] = AffineManifold.FromPoints(VA.VTX, VB.VTX);

                    double sign = 0;
                    foreach(var vtx in this.Vertices) {
                        sign += DomainCheck[ie].PointDistance(vtx.VTX);
                    }

                    if(sign > 0) {
                        DomainCheck[ie].Normal.Scale(-1.0);
                        DomainCheck[ie].a *= -1.0;
                    }
                }
            }

            internal bool CheckOriginalDomain() {

                foreach (var v in this.Vertices) {
                    foreach (var edge in DomainCheck) {
                        if (edge.PointDistance(v.VTX) >= hmin * 1e-3) { // all inside vertices should have a *negative* distance
                            return false;
                        }
                    }
                }
                return true;
            }


            static int IDcounter = 1;

            public IReadOnlyList<VoEdge> Edges {
                get {
                    return m_Edges.AsReadOnly();
                }
            }


            List<VoEdge> m_Edges = new List<VoEdge>();

            // hack, nicht im sinne der Kapselung
            internal void AddEdge(VoEdge e) {
                if (m_Edges.Contains(e)) {
                    return;
                }
                m_Edges.Add(e);
            }

            internal void ReplaceEdge(VoEdge old, VoEdge nju) {
                if (!this.Edges.ContainsRefEqual(old))
                    throw new ApplicationException();
                if (!old.Equals(nju))
                    throw new ApplicationException();

                int i = IEnumerableExtensions.IndexOf(this.m_Edges, old);
                this.m_Edges[i] = nju;

                if (!nju.Cells.ContainsRefEqual(this))
                    nju.Cells.Add(this);
            }


            // hack, nicht im Sinne der Kapselung
            public void ClearEdges() {
                m_Edges.Clear();
            }


            /// <summary>
            /// true when all vertices are <see cref="VertexType.Inside"/>
            /// </summary>
            public bool IsInside {
                get {
                    bool r = true;

                    foreach (var e in Edges) {
                        if (e.VtxA.type != VertexType.Inside) {
                            r = false;
                            break;
                        }
                        if (e.VtxB.type != VertexType.Inside) {
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

                    foreach (var e in Edges) { // loop over edges
                        if (e.VtxA.type != VertexType.Outside) {
                            r = false;
                            break;
                        }
                        if (e.VtxB.type != VertexType.Outside) {
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
            public void RemoveOutsideParts(bool AlsoBoundary, Func<Vector, bool> IsIn) {
                for (int ie = 0; ie < m_Edges.Count; ie++) { // loop over edges
                    VoEdge e = m_Edges[ie];

                    if (e.VtxA.type == VertexType.unspecified) {
                        throw new InvalidOperationException();
                    }

                    if(   (e.VtxA.type == VertexType.Outside && e.VtxB.type == VertexType.Inside)
                       || (e.VtxA.type == VertexType.Inside && e.VtxB.type == VertexType.Outside)) {
                        throw new InvalidOperationException("found a non-intersected edge.");

                    }

                    bool bRemove = false;
                    if (e.VtxA.type == VertexType.Boundary && e.VtxB.type == VertexType.Boundary) {
                        Vector near_a = e.Interpol(0.001);
                        Vector near_b = e.Interpol(0.999);

                        bool in_near_a = IsIn(near_a);
                        bool in_near_b = IsIn(near_b);

                        if(in_near_a && in_near_b) {
                            // edge allowed to stay
                            continue;
                        }

                        if( !in_near_a && !in_near_b) {
                            bRemove = true;
                        }

                        if (in_near_b != in_near_a) {
                            /*
                            DebugPlot(null, null, this.Edges, new VoVertex[] { e.VtxA, e.VtxB });


                            throw new InvalidOperationException("found a non-intersected edge.");
                            */
                            bRemove = true;

                        }
                    }


                    if ((e.VtxA.type == VertexType.Outside)
                       || (e.VtxB.type == VertexType.Outside)
                       || (e.VtxA.type == VertexType.FarPoint)
                       || (e.VtxB.type == VertexType.FarPoint)
                       || e.VtxA.IsFar
                       || e.VtxB.IsFar
                       || (AlsoBoundary && e.isBoundary)
                       || bRemove
                        ) {

                        // remove linking of edge to this cell
                        VoEdge edge2remove = m_Edges[ie];
                        int i1 = edge2remove.Cells.IndexOf(this);
                        edge2remove.Cells.RemoveAt(i1);

                        // remove edge
                        m_Edges.RemoveAt(ie);
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


            public VoVertex[] GetIntersectionSequence(VoPolygon bndyPoly) {

                VoVertex[][] seqS = GetVerticesSequence(out bool isClosed, false);

              


                if (seqS.Length == 0)
                    return new VoVertex[0];
                if(seqS.All(seq => seq.Length == 0))
                    return new VoVertex[0];

                if (isClosed) {
                    Debug.Assert(seqS.Length == 1);
                    return seqS[0];
                }

                //if (seqS.Length > 1)
                //    throw new NotSupportedException();

                if (seqS.Length == 1) {
                    int sign = CheckOrientation(seqS[0]);
                    if (sign == 0) {
                        //DebugPlot(VocellVertexIndex, Verts, Cj.Edges);
                        //Console.WriteLine("indef polygon");
                        //NoOfIndef++;
                        //continue;
                        //throw new ArithmeticException("indefinite polygon.");
                    }
                    if (sign < 0)
                        seqS[0] = seqS[0].Reverse().ToArray();
                    //Debug.Assert(CheckOrientation(seqS[0]) > 0);


                    Debug.Assert(isClosed == false);


                    VoVertex[] closing = bndyPoly.GetSegmentBetween(seqS[0].Last(), seqS[0][0]);
                    if (closing.Length > 2) {
                        VoVertex[] closing_inner = closing.GetSubVector(1, closing.Length - 2);

                        seqS[0] = seqS[0].Cat(closing_inner);
                    }

                    return seqS[0];
                } else {

                    List<VoVertex> R = new List<VoVertex>();
                    bool[] seqUsed = new bool[seqS.Length];
                    R.AddRange(seqS[0]);
                    seqUsed[0] = true;

                    VoVertex start = R[0];
                    VoVertex current = R.Last();


                    while (seqUsed.Any(b => !b)) {

                        double minlength = double.MaxValue;
                        VoVertex[] close = null;
                        int iNext = -1;
                        for (int i = 0; i < seqS.Length; i++) {
                            if (seqUsed[i])
                                continue;
                            Debug.Assert(i != 0); // i need the sign of i

                            VoVertex[] nextTry = seqS[i];


                            VoVertex[] close1 = bndyPoly.GetSegmentBetween(current, nextTry[0]);
                            double len1 = Length(close1);
                            Debug.Assert(close1.Length > 1);
                            if (len1 < minlength) {
                                close = close1;
                                minlength = len1;
                                iNext = i;
                            }

                            VoVertex[] close2 = bndyPoly.GetSegmentBetween(current, nextTry[nextTry.Length - 1]);
                            double len2 = Length(close2);
                            Debug.Assert(close2.Length > 1);
                            if (len2 < minlength) {
                                close = close2;
                                minlength = len2;
                                iNext = -i;
                            }
                        }

                        if (close.Length > 2) {
                            R.AddRange(close.GetSubVector(1, close.Length - 2));
                        }

                        seqUsed[Math.Abs(iNext)] = true;
                        if (iNext > 0) {
                            R.AddRange(seqS[iNext]);
                        } else if (iNext < 0) {
                            R.AddRange(seqS[-iNext].Reverse());
                        } else {
                            throw new ApplicationException("error in algorithm");
                        }

                        current = R.Last();
                    }

                    if (current.Equals(start))
                        throw new ApplicationException("error in algorithm");


                    VoVertex[] finalClosing = bndyPoly.GetSegmentBetween(current, start);
                    if (finalClosing.Length > 2) {
                        VoVertex[] closing_inner = finalClosing.GetSubVector(1, finalClosing.Length - 2);

                        R.AddRange(closing_inner);
                    }
                    /*
                    using(var gp = new Gnuplot()) {
                        var x = R.Select(X => X.VTX.x).ToArray();
                        var y = R.Select(X => X.VTX.y).ToArray();
                        gp.PlotXY(x, y);
                        gp.Execute();
                        Console.ReadKey();
                    }
                    */

                    return R.ToArray();
                }
            }

            static double Length(VoVertex[] seq) {
                double a = 0;
                for (int i = 1; i < seq.Length; i++) {
                    a += seq[i].VTX.Dist(seq[i - 1].VTX);
                }

                return a;
            }


            public VoVertex[][] GetVerticesSequence(out bool isClosed, bool includeBoundaries = false) {
                List<VoVertex[]> R = new List<VoVertex[]>();
                isClosed = false;

                MyTuple[] EdgesLoc = this.Edges
                    .Where(edg => includeBoundaries || edg.isBoundary == false)
                    .Select(edg => new MyTuple(edg, false))
                    .ToArray();

                while (EdgesLoc.Any(ee => ee.used == false)) {
                    VoVertex[] rr;
                    isClosed = isClosed | VertSeq(out rr, EdgesLoc);
                    R.Add(rr);
                }

                //Debug.Assert(EdgesLoc.All(e => e.used),"öha");

                if (!EdgesLoc.All(e => e.used)) {
                    //DebugPlot(null, null, this.Edges, new VoVertex[]{s);
                    throw new ApplicationException("error in algorithm");
                }

                if (isClosed)
                    Debug.Assert(R.Count == 1);

                return R.ToArray();
            }

            private static bool VertSeq(out VoVertex[] _R, MyTuple[] EdgesLoc) {
                bool isClosed = false;
                List<VoVertex> R = new List<VoVertex>();

                MyTuple CurrEdge = EdgesLoc.First(ee => ee.used == false); // select any edge
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

                if (SearchReverse) {
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
                        if (PointIdentity(NextEdge.edge.VtxA, CurrVertex)) {
                            NextVertex = NextEdge.edge.VtxB;
                        } else {
                            Debug.Assert(PointIdentity(NextEdge.edge.VtxB, CurrVertex));
                            NextVertex = NextEdge.edge.VtxA;
                        }
                        NextEdge.used = true;
                        R.Insert(0, NextVertex);
                    }
                }

                _R = R.ToArray();
                return isClosed;
            }

            /// <summary>
            /// vertices of this polygon
            /// </summary>
            public IReadOnlyList<VoVertex> Vertices {
                get {
                    List<VoVertex> R = new List<VoVertex>();
                    foreach(var edg in this.Edges) {
                        if (!R.ContainsExactly(edg.VtxA))
                            R.Add(edg.VtxA);
                        if (!R.ContainsExactly(edg.VtxB))
                            R.Add(edg.VtxB);
                    }

                    return R.AsReadOnly();
                }
            }


            public VoVertex[] GetSegmentBetween(VoVertex Start, VoVertex End) {
                if (PointIdentity(Start, End))
                    throw new ArgumentException();

                VoVertex[][] __all = GetVerticesSequence(out bool isClosed, true);
                if (!isClosed)
                    throw new NotSupportedException();
                if (isClosed && __all.Length > 1)
                    throw new NotSupportedException();
                var all = __all[0];



                //if (!all.Contains(Start))
                //    throw new ArgumentException("start vertex is not an element of this polygon.");
                //if (!all.Contains(End))
                //    throw new ArgumentException("end vertex is not an element of this polygon.");

                // 
                List<VoVertex>[] Segments = new List<VoVertex>[2];
                double[] Length = new double[2];
                for (int iSweep = 0; iSweep < 2; iSweep++) { // sweep 1: search forward; step 2: search backward;
                    Segments[iSweep] = new List<VoVertex>();

                    int i0 = all.IndexOf(Start, PointIdentity);
                    if (i0 < 0)
                        throw new ArgumentException("Start vertex is not an element of this polygon.", "Start");

                    bool bFound = false;
                    for (int i = 0; i < all.Length; i++) {
                        int iCnt = (i0 + i) % all.Length;
                        Segments[iSweep].Add(all[iCnt]);

                        if (i > 0) {
                            Length[iSweep] += Segments[iSweep][i].VTX.Dist(Segments[iSweep][i - 1].VTX);
                        }

                        if (PointIdentity(all[iCnt], End)) {
                            bFound = true;
                            break;
                        }
                    }

                    if (!bFound) {
                        DebugPlot(null, null, this.Edges, new VoVertex[] { Start, End });


                        throw new ArgumentException("End vertex is not an element of this polygon.", "End");
                    }
                    // check
                    Debug.Assert(PointIdentity(Start, Segments[iSweep].First()));
                    Debug.Assert(PointIdentity(End, Segments[iSweep].Last()));

                    // reverse list: 
                    all = all.Reverse().ToArray();
                }

                // return shorter segment
                if (Length[0] < Length[1]) {
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

            return (AB * AP) / AB.AbsSquare();

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
            if (object.ReferenceEquals(V1, V2)) {
                Debug.Assert(PointIdentityG(V1.VTX, V2.VTX));
                return true;
            } else {
                return PointIdentityG(V1.VTX, V2.VTX);
            }


        }


        static bool Intersect(VoEdge edge, VoEdge bndy, List<VoEdge> bndyEdges) {
            // Colinear: exact overlap
            // - - - - - - - - - - - - 

            if (edge.Equals(bndy)) {
                // +++++++++++++++++++++++++++++++++++++++++++
                // special case: edge and bndy overlap exactly
                // => remove edge
                // +++++++++++++++++++++++++++++++++++++++++++

                //bndy.Cells.AddRange(edge.Cells);
                //foreach (var cell in edge.Cells) {
                //    Debug.Assert(cell.Edges.Contains(edge));
                //    cell.Edges.Remove(edge);
                //    cell.Edges.Add(bndy);
                //}

                //VoEdge.edgeS.RemoveAt(iEdge);
                //Debug.Assert(false);

                //edge.deleted = true;
                //iEdge--;
                //continue;
                throw new ApplicationException("Error in graph - duplicate edges");
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
                    CheckEdgeUniqueness();
                    return false;
                }

                if (alphaA >= 1.0) {
                    // no overlap
                    CheckEdgeUniqueness();
                    return false;
                }


                if (alphaA == 0 && alphaB < 1) {
                    // edge:  o---o
                    // bndy:  o--------o
                    //
                    // out:   o---o----o
                    
                    bndyEdges.MyRemove(bndy);
                    bndy.Split(edge.VtxB, out VoEdge pt1, out VoEdge pt2);
                    pt1.isBoundary = true;
                    pt2.isBoundary = true;
                    bndyEdges.MySetAdd(pt1);
                    bndyEdges.MySetAdd(pt2);

                    // cont
                    CheckEdgeUniqueness();
                    return true;
                }

                if (alphaB == 1 && alphaA > 0) {
                    // edge:       o---o
                    // bndy:  o--------o
                    //
                    // out:   o----o---o

                    bndyEdges.MyRemove(bndy);
                    bndy.Split(edge.VtxA, out VoEdge pt1, out VoEdge pt2);
                    pt1.isBoundary = true;
                    pt2.isBoundary = true;
                    bndyEdges.MySetAdd(pt1);
                    bndyEdges.MySetAdd(pt2);


                    // cont
                    CheckEdgeUniqueness();
                    return true;
                }

                if (alphaA == 0 && alphaB > 1) {
                    // edge:  o------------o
                    // bndy:  o--------o
                    //
                    // out:   o--------o---o

                   
                    bndyEdges.MyRemove(bndy);
                    edge.Split(bndy.VtxB, out VoEdge pt1, out VoEdge pt2);
                    bndyEdges.MySetAdd(pt1);

                    // cont
                    CheckEdgeUniqueness();
                    return true;
                }

                if (alphaA < 0 && alphaB == 1) {
                    // edge:  o------------o
                    // bndy:      o--------o
                    //
                    // out:   o---o--------o

                    bndyEdges.MyRemove(bndy);
                    edge.Split(bndy.VtxA, out VoEdge pt1, out VoEdge pt2);
                    bndyEdges.MySetAdd(pt2);

                    // cont
                    CheckEdgeUniqueness();
                    return true;
                }

                Debug.Assert(alphaA != 0.0); // all special cases should treated by now
                Debug.Assert(alphaB != 1.0); // all special cases should treated by now

                if (alphaA < 0 && alphaB > 0) {
                    // edge:  o-----------------o
                    // bndy:      o--------o
                    //
                    // out:   o---o--------o~~~~o

                                        
                    VoVertex I1 = bndy.VtxA;
                    VoVertex I2 = bndy.VtxB;
                    bndyEdges.MyRemove(bndy);
                    edge.Split(I1, out VoEdge pt1, out VoEdge temp);
                    temp.Split(I2, out VoEdge pt2, out VoEdge pt3);
                    pt2.isBoundary = true;
                    bndyEdges.Add(pt2);
                    

                    // cont
                    CheckEdgeUniqueness();
                    return true;
                }


                if (alphaA > 0 && alphaB < 1) {
                    // edge:      o-o
                    // bndy:   o--------o
                    //
                    // out:    o--o-o~~~o
                    Debug.Assert(alphaB > 0);

                   
                    VoVertex I1 = edge.VtxA;
                    VoVertex I2 = edge.VtxB;
                    bndyEdges.MyRemove(bndy);
                    bndy.Split(I1, out VoEdge pt1, out VoEdge temp);
                    temp.Split(I2, out VoEdge pt2, out VoEdge pt3);
                    pt1.isBoundary = true;
                    pt2.isBoundary = true;
                    pt3.isBoundary = true;
                    bndyEdges.MySetAdd(pt1);
                    bndyEdges.MySetAdd(pt2);
                    bndyEdges.MySetAdd(pt3);


                    // cont
                    return true;
                }

                if (alphaA < 0 && alphaB < 1) {
                    // edge:   o-----o
                    // bndy:      o--------o
                    //
                    // out:    o~~o--o-----o

                    /*
                    // pt1 (new)
                    var newEdge = VoEdge.Create(edge.VtxA, bndy.VtxA);
                    newEdge.Cells.AddRange(edge.Cells);
                    Debug.Assert(newEdge.isBoundary == false); // outside domain

                    // pt2 (edge)
                    edge.VtxA = bndy.VtxA;
                    edge.isBoundary = true;
                    bndyEdges.MySetAdd(edge);

                    // pt3 (bndy)
                    bndy.VtxA = edge.VtxB;
                    bndy.isBoundary = true;
                    */

                    VoVertex I1 = bndy.VtxA;
                    VoVertex I2 = edge.VtxB;
                    bndyEdges.MyRemove(bndy);
                    edge.Split(I1, out VoEdge pt1, out VoEdge temp);
                    bndy.Split(I2, out VoEdge pt2, out VoEdge pt3);
                    pt2.isBoundary = true;
                    pt3.isBoundary = true;
                    bndyEdges.MySetAdd(pt2);
                    bndyEdges.MySetAdd(pt3);

                    // cont
                    return true;
                }

                if (alphaA > 0 && alphaB > 1) {
                    // edge:        o-----o
                    // bndy:   o--------o
                    //
                    // out:    o----o---o~o

                    /*
                    // pt1 (bndy)
                    var t = bndy.VtxB;
                    bndy.VtxB = edge.VtxA;
                    bndy.isBoundary = true;

                    // pt2 (edge)
                    edge.isBoundary = true;
                    var tt = edge.VtxB;
                    edge.VtxB = t;
                    bndyEdges.MySetAdd(edge);

                    // pt3 (new)
                    var newEdge = VoEdge.Create(t, tt);
                    newEdge.Cells.AddRange(edge.Cells);
                    Debug.Assert(newEdge.isBoundary == false); // outside domain
                    CheckEdgeUniqueness();
                    */

                    var II = edge.VtxA;
                    var I2 = bndy.VtxB;
                    bndyEdges.MyRemove(bndy);
                    bndy.Split(II, out VoEdge pt1, out VoEdge temp);
                    edge.Split(I2, out VoEdge pt2, out VoEdge pt3);
                    pt1.isBoundary = true;
                    pt2.isBoundary = true;
                    bndyEdges.Add(pt1);
                    bndyEdges.Add(pt2);


                    // cont
                    CheckEdgeUniqueness();
                    return true;
                }



            } else {
                // +++++++++++++++++++++++++
                // edges are NOT colinear
                // +++++++++++++++++++++++++

                if (PointIdentity(bndy.VtxA, edge.VtxA)) { CheckEdgeUniqueness(); return false; } // L-junction
                if (PointIdentity(bndy.VtxA, edge.VtxB)) { CheckEdgeUniqueness(); return false; } // L-junction
                if (PointIdentity(bndy.VtxB, edge.VtxA)) { CheckEdgeUniqueness(); return false; } // L-junction
                if (PointIdentity(bndy.VtxB, edge.VtxB)) { CheckEdgeUniqueness(); return false; } // L-junction


                bool cutfound = bndy.Intersect(edge, out double alpha, out double beta, out Vector I);
                // alpha - coordinate: bndy
                // beta  - coordinate: edge

                //if ((I.x - 0.727).Abs() <= 0.02 && (I.y - 1.0).Abs() < 0.01)
                //    Console.Write("");

                if (!cutfound)
                    return false;

                if ((alpha == 0 || alpha == 1) && (beta == 0.0 || beta == 1.0)) {
                    // L-junction: nothing to do
                    // this *should* have been caught already by the L-junction detection upwards
                    throw new ApplicationException("should never reach this point - error in algorithm");
                }

                if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
                    // no intersection within the finite line segment
                    return false;
                }

                if ((beta == 0.0 || beta == 1.0) && (alpha > 0.0 && alpha < 1.0)) {
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

                    bndyEdges.MyRemove(bndy);
                    bndy.Split(newVert, out VoEdge pt1, out VoEdge pt2);
                    pt1.isBoundary = true;
                    pt2.isBoundary = true;
                    bndyEdges.MySetAdd(pt1);
                    bndyEdges.MySetAdd(pt2);


                    //
                    CheckEdgeUniqueness();
                    return true;

                }

                if ((beta > 0.0 || beta < 1.0) && (alpha == 0.0 || alpha == 1.0)) {
                    // T-junction

                    // find junction vertex
                    VoVertex newVert;
                    if (alpha == 0.0) {
                        newVert = bndy.VtxA;
                    } else {
                        Debug.Assert(alpha == 1.0);
                        newVert = bndy.VtxB;
                    }

                    edge.Split(newVert, out VoEdge ept1, out VoEdge ept2);
                    return true;

                }

                if ((beta > 0.0 || beta < 1.0) && (alpha > 0.0 || alpha < 1.0)){
                    // X-junction

                    // introduce new vertex
                    var newVert = VoVertex.Create(I);

                    bndyEdges.MyRemove(bndy);
                    bndy.Split(newVert, out VoEdge pt1, out VoEdge pt2);
                    pt1.isBoundary = true;
                    pt2.isBoundary = true;
                    bndyEdges.MySetAdd(pt1);
                    bndyEdges.MySetAdd(pt2);

                    edge.Split(newVert, out VoEdge ept1, out VoEdge ept2);


                    //
                    return true;
                }

                throw new ApplicationException("should never reach this point - error in algorithm");

            } // edges are NOT colinear

            throw new ApplicationException("should never reach this point - error in algorithm");
        }



        static public AggregationGrid FromPolygonalDomain(MultidimensionalArray Nodes, Vector[] PolygonBoundary, bool mirroring, int NoOfLyyodsIter, Func<Vector, bool> IsIn, Func<Vector, Vector, bool> __PointIdentity) {

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


           

            // Create Voronoi mesh (call Matlab)
            // =================================

            int[][] VocellVertexIndex;
            List<Vector> Verts;
            MultidimensionalArray DelaunayVertices;
            {
                VocellVertexIndex = null;
                Verts = null;
                for (int iLloyd = 0; iLloyd < NoOfLyyodsIter; iLloyd++) {

                    // Lloyds algorithm (Voronoi relaxation)
                    // -------------------------------------
                    if (iLloyd > 0) {
                        for (int jV = 0; jV < Nodes.NoOfRows; jV++) {

                            Vector oldNode = Nodes.GetRowPt(jV);

                            if (IsIn(oldNode)) {
                                // replace, for each cell, the node with the center-of-gravity

                                int[] cell = VocellVertexIndex[jV];
                                Vector COG = new Vector(0, 0);
                                for (int k = 0; k < cell.Length; k++) {
                                    COG += Verts[cell[k]];
                                }
                                COG.Scale(1.0 / cell.Length);

                                // apply gravity
                                double x = COG.x;
                                double y = COG.y;
                                double a = 2.0;
                                Vector G = new Vector(
                                    -2 * x * a * a * Math.Exp(-(x * x + y * y) * a * a),
                                    -2 * y * a * a * Math.Exp(-(x * x + y * y) * a * a)
                                    );
                                COG += G * 0.1;


                                // 
                                Nodes.SetRowPt(jV, COG);
                            }
                        }
                    }

                    // Optional Vertex Mirroring
                    // -------------------------
                    if (mirroring) {

                        double[] xNodes = Nodes.GetColumn(0);
                        double[] yNodes = Nodes.GetColumn(1);

                        Mirror(ref xNodes, ref yNodes, Boundaries, IsIn);
                        Debug.Assert(xNodes.Length == yNodes.Length);

                        DelaunayVertices = MultidimensionalArray.Create(xNodes.Length, 2);
                        DelaunayVertices.SetColumn(0, xNodes);
                        DelaunayVertices.SetColumn(1, yNodes);

                    } else {
                        DelaunayVertices = Nodes.CloneAs();
                    }

                    // Voronoi generation using matlab
                    // --------------------------------

                    using (var Matlab = new BatchmodeConnector()) {
                        Matlab.PutMatrix(DelaunayVertices, "Nodes");

                        // compute Voronoi diagramm
                        Matlab.Cmd("[V, C] = voronoin(Nodes);");

                        // output (export from matlab)
                        VocellVertexIndex = new int[DelaunayVertices.NoOfRows][];
                        Matlab.GetStaggeredIntArray(VocellVertexIndex, "C");
                        Matlab.GetMatrix(null, "V");

                        // run matlab
                        Matlab.Execute(false);

                        // import here
                        MultidimensionalArray VertexCoordinates = (MultidimensionalArray)(Matlab.OutputObjects["V"]);
                        Verts = new List<Vector>();
                        for (int iVtx = 0; iVtx < VertexCoordinates.NoOfRows; iVtx++) {
                            Verts.Add(VertexCoordinates.GetRowPt(iVtx));
                        }

                        // correct indices (1-based index to 0-based index)
                        foreach (int[] cell in VocellVertexIndex) {
                            int K = cell.Length;
                            for (int k = 0; k < K; k++) {
                                cell[k]--;
                            }
                        }
                    }
                }



                // fix Voronoi cell orientation
                // ----------------------------
                for (int jV = 0; jV < VocellVertexIndex.Length; jV++) {
                    int[] iVtxS = VocellVertexIndex[jV];
                    if (iVtxS.Any(i => IsFarPoint(Verts[i])))
                        continue;

                    Vector[] vectors = iVtxS.Select(i => Verts[i]).ToArray();
                    //Vector[] _vectors = vectors;
                    int sign = CheckOrientation(vectors);
                    if (sign > 0) {
                        // nop
                    } else if (sign < 0) {
                        iVtxS = iVtxS.Reverse().ToArray();
                    } else {
                        throw new ArithmeticException("got indefinite polygon form matlab");
                    }
                    VocellVertexIndex[jV] = iVtxS;
                }

                //
                Nodes = null;
            }


            // build data structures
            // =====================
            List<int> verticesIndices = new List<int>();
            for(int i = 0; i < Verts.Count; i++) {
                var v = VoVertex.Create(Verts[i]);
                if(v.ID != VoVertex.verticeS.Last().ID) {
                    throw new ArithmeticException("Matlab produced indistinguishable Voronoi vertices.");
                }
                verticesIndices.Add(VoVertex.verticeS.Count - 1);
            }

            List<VoPolygon> cellS = new List<VoPolygon>();
            {

                //bool anyDouble = false;
                for (int jV = 0; jV < VocellVertexIndex.Length; jV++) {

                    
                    int[] iVtxS = VocellVertexIndex[jV];

                    int I = iVtxS.Length;
                    List<VoEdge> edges_jV = new List<VoEdge>(); 
                    for (int i = 0; i < I; i++) {
                        
                        int _iVtxA = iVtxS[i];
                        int _iVtxB = iVtxS[(i + 1) % I];
                        
                        VoVertex _VtxA = VoVertex.verticeS[verticesIndices[_iVtxA]];
                        VoVertex _VtxB = VoVertex.verticeS[verticesIndices[_iVtxB]];

                        if (_VtxA.IsFar)
                            continue;
                        if (_VtxB.IsFar)
                            continue;

                        VoEdge newEdge = VoEdge.Create(_VtxA, _VtxB);

                        //int iEdge = VoEdge.edgeS.IndexOf(newEdge);
                        //if(iEdge < 0) {
                        //    for(int e = 0; e < edgeS.Count; e++) {
                        //        Debug.Assert(edgeS[e].Equals(newEdge) == false);
                        //    }

                        //    edgeS.Add(newEdge);
                        //    iEdge = edgeS.Count - 1;
                        //    //newEdge.Index = iEdge;
                        //} else {
                        //    anyDouble = true;
                        //    if (edgeS[iEdge].Cells.Count > 1)
                        //        throw new ArithmeticException("edge shared by more than two cells");
                        //    newEdge = edgeS[iEdge];
                        //}
                        //Debug.Assert(newEdge.Equals(edgeS[iEdge]));
                        Debug.Assert(VoEdge.edgeS.Where(ve => ve.Equals(newEdge)).Count() == 1);

                        edges_jV.Add(newEdge);
                    }
                    var cell = new VoPolygon(edges_jV, true);
                    cellS.Add(cell);
                    foreach (var e in cell.Edges)
                        e.CheckCell2Edge();
                }

                foreach(var e in VoEdge.edgeS)
                    e.CheckCell2Edge();

                //if (!anyDouble)
                //    throw new ArithmeticException("Voronoi diagram seems to be completely disjoint - no edge is used by at least two cells.");
            }

            //DelaunayVertices.NoOfRows

         
            // Add boundary edges 
            // ==================
            {
                List<int> idxBndyVertices = new List<int>();

                for (int i = 0; i < PolygonBoundary.Length; i++) {
                    Vector bVtx = PolygonBoundary[i];

                    //int idxFound = -1;
                    //for (int j = 0; j < VoVertex.verticeS.Count; j++) {
                    //    if (PointIdentityG(bVtx, VoVertex.verticeS[j].VTX)) {
                    //        if (idxFound >= 0) {
                    //            throw new ArithmeticException();
                    //        }
                    //    }
                    //}

                    var newVtx = VoVertex.Create(bVtx);
                    int idxFound = VoVertex.verticeS.IndexWhere(v => v.ID == newVtx.ID);
                    
                    idxBndyVertices.Add(idxFound);
                    newVtx.type = VertexType.Boundary;
                }


                for(int i = 0; i < PolygonBoundary.Length; i++) {
                    int i0 = i;
                    int iE = (i + 1) % PolygonBoundary.Length;

                    VoVertex V0 = VoVertex.verticeS[idxBndyVertices[i0]];
                    VoVertex VE = VoVertex.verticeS[idxBndyVertices[iE]];
                    Debug.Assert(PointIdentityG(V0.VTX, PolygonBoundary[i0]));
                    Debug.Assert(PointIdentityG(VE.VTX, PolygonBoundary[iE]));

                    var newEdge = VoEdge.Create(V0, VE);
                    newEdge.isBoundary = true;
                }
            }

            // compute intersections
            // =====================
            CheckEdgeUniqueness(true);
            {
                List<VoEdge> bndyEdges = VoEdge.edgeS.Where(edge => edge.isBoundary).ToList();
                foreach(var e in bndyEdges)
                    e.CheckCell2Edge();
                foreach(var e in VoEdge.edgeS)
                    e.CheckCell2Edge();

                int iRun = 0;
                bool gefinden = true;

                while (gefinden) {
                    gefinden = false;
                    for (int iBndy = 0; iBndy < bndyEdges.Count && gefinden == false; iBndy++) {
                        for (int iEdge = 0; iEdge < VoEdge.edgeS.Count && gefinden == false; iEdge++) {
                            iRun++;
                            VoEdge bndy = bndyEdges.ElementAt(iBndy);
                            VoEdge edge = VoEdge.edgeS[iEdge];
                            CheckEdgeUniqueness();
                            bool iscont = VoEdge.edgeS.Where(eee => eee.isBoundary).ToArray().ContainsExactly(bndy);

                            bndy.CheckCell2Edge();
                            edge.CheckCell2Edge();

                            if (edge.isBoundary) {
                                continue;
                            }
                            Debug.Assert(bndy.isBoundary == true);
                            Debug.Assert(bndy.VtxA.type == VertexType.Boundary);
                            Debug.Assert(bndy.VtxB.type == VertexType.Boundary);

                            //foreach (var eb in bndyEdges) {
                            //    eb.CheckCell2Edge();
                            //}

                            //if (iRun == 784) {
                            //    Debug.Assert(VoEdge.edgeS.Where(eee => eee.isBoundary).SetEquals(bndyEdges));
                            //    Debugger.Break();
                            //}

                            bool intsc = Intersect(edge, bndy, bndyEdges);


                            if (intsc) {
                                iEdge--;
                                gefinden = true; // both lists ('bndyEdges', 'edgeS') may have changed 
                                //                  the easiest approach is to restart all loops



                                Debug.Assert(VoEdge.edgeS.Where(eee => eee.isBoundary).SetEquals(bndyEdges));

                                /*
                                if (iRun == 784) {
                                    var A = bndyEdges.ToArray();
                                    var B = VoEdge.edgeS.Where(äää => äää.isBoundary).ToArray();

                                    List<VoEdge> NotInA_r = new List<VoEdge>();
                                    List<VoEdge> NotInA_e = new List<VoEdge>();
                                    foreach (var b in B) {
                                        if (!A.ContainsRefEqual(b))
                                            NotInA_r.Add(b);

                                        var bAs = A.Where(a => a.Equals(b)).ToArray();
                                        Debug.Assert(bAs.Length <= 1);
                                        if (bAs.Length <= 0)
                                            NotInA_e.Add(b);
                                    }
                                    List<VoEdge> NotInB_r = new List<VoEdge>();
                                    List<VoEdge> NotInB_e = new List<VoEdge>();
                                    foreach (var a in A) {
                                        if (!B.ContainsRefEqual(a))
                                            NotInB_r.Add(a);

                                        var aBs = B.Where(b => b.Equals(a)).ToArray();
                                        Debug.Assert(aBs.Length <= 1);
                                        if (aBs.Length <= 0)
                                            NotInB_e.Add(a);
                                    }
                                }
                                */
                                //foreach(var eb in bndyEdges) {
                                //    eb.CheckCell2Edge();
                                //}
                            }
                        }
                    }
                }
                

                {
                    /*
                    var A = bndyEdges.ToArray();
                    var B = VoEdge.edgeS.Where(edge => edge.isBoundary).ToArray();

                    List<VoEdge> NotInA_r = new List<VoEdge>();
                    List<VoEdge> NotInA_e = new List<VoEdge>();
                    foreach (var b in B) {
                        if (!A.ContainsRefEqual(b))
                            NotInA_r.Add(b);

                        var bAs = A.Where(a => a.Equals(b)).ToArray();
                        Debug.Assert(bAs.Length <= 1);
                        if (bAs.Length <= 0)
                            NotInA_e.Add(b);
                    }
                    List<VoEdge> NotInB_r = new List<VoEdge>();
                    List<VoEdge> NotInB_e = new List<VoEdge>();
                    foreach (var a in A) {
                        if (!B.ContainsRefEqual(a))
                            NotInB_r.Add(a);

                        var aBs = B.Where(b => b.Equals(a)).ToArray();
                        Debug.Assert(aBs.Length <= 1);
                        if (aBs.Length <= 0)
                            NotInB_e.Add(a);
                    }
                    */

                    Debug.Assert(VoEdge.edgeS.Where(edge => edge.isBoundary).SetEquals(bndyEdges));
                    bndyEdges = null; // forget
                }


                CheckEdgeUniqueness(true);
            }

            //DebugPlot(VocellVertexIndex, Verts, VoEdge.edgeS);

            // final point classification
            // ==========================

            foreach(var V in VoVertex.verticeS) {
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

            // test cell-to-edges, vertices-to-edges
            // =======================================
            {
#if DEBUG
                //VoEdge[][] cell2edge_check = cellS.Select(cell => cell.Edges.ToArray()).ToArray();

                for (int ie = 0; ie < VoEdge.edgeS.Count; ie++) {
                    var edge = VoEdge.edgeS[ie];
                    edge.CheckCell2Edge();
                }


                //for (int jV = 0; jV < cellS.Count; jV++){
                //    var cell = cellS[jV];
                //    //if (cell.Edges.Count <= 1)
                //    //    Console.WriteLine("warn: less than 2");
                //    cell.ClearEdges();
                //}

                List<VoEdge>[] cell2edge_check = cellS.Count.ForLoop(i => new List<VoEdge>());


                for(int ie = 0; ie < VoEdge.edgeS.Count; ie++) {
                    var edge = VoEdge.edgeS[ie];
                    edge.SanitizeCells();

                    for (int ic = 0; ic < edge.Cells.Count; ic++) {
                        var cell = edge.Cells[ic];
                        int iCell = cellS.IndexOf(cell);
                        if (!cell2edge_check[iCell].ContainsExactly(edge))
                            cell2edge_check[iCell].Add(edge);

                        //Debug.Assert(cell.Edges.ContainsRefEqual(edge) == false);
                        //cell.AddEdge(edge);

                    }

                    /*
                    if (edge.VtxA.Edges == null)
                        edge.VtxA.Edges = new List<VoEdge>();
                    if (edge.VtxB.Edges == null)
                        edge.VtxB.Edges = new List<VoEdge>();

                    if (!edge.VtxA.Edges.ContainsRefEqual(edge))
                        edge.VtxA.Edges.Add(edge);
                    if (!edge.VtxB.Edges.ContainsRefEqual(edge))
                        edge.VtxB.Edges.Add(edge);
                    */
                }
                

                for(int jV = 0; jV < cellS.Count; jV++){
                    var cell = cellS[jV];
                    //if (cell.Edges.Count <= 1)
                    //    Console.WriteLine("warn2: less than 2");
                    Debug.Assert(cell2edge_check[jV].SetEquals(cell.Edges));
                }
#endif
            }

            // collect boundary polygon
            // ========================
            VoPolygon bndyPoly;
            {
                bndyPoly = new VoPolygon(VoEdge.edgeS.Where(edge => edge.isBoundary == true), false);
                                
                bndyPoly.GetVerticesSequence(out bool bndyClosed, true); // only purpose: test if bndy polygon is closed 
                if (!bndyClosed)
                    throw new Exception();
            }

            /*
            {
                for (int j = 0; j < cellS.Count; j++) {
                    VoPolygon Cj = cellS[j];
                    var vtxS = Cj.Vertices.ToArray();
                    bool contain_Arsch = vtxS.Any(vtx => (vtx.VTX.x - (-1.0)).Abs() <= 0.001 && (vtx.VTX.y - (0.207)).Abs() <= 0.01);
                    bool contain_Brsch = vtxS.Any(vtx => (vtx.VTX.x - (-0.895)).Abs() <= 0.01 && (vtx.VTX.y - (0.0)).Abs() <= 0.001);
                    if (contain_Arsch && contain_Brsch) {

                        Console.WriteLine("poly  " + j + " is outside? " + Cj.IsOutside);
                        DebugPlot(null, null, Cj.Edges, null);
                    }
                }
            }
            */

            // collect inside polygons
            // =======================
            List<VoVertex[]> Insiders = new List<VoVertex[]>();
            for(int j = 0; j < cellS.Count; j++) {
                VoPolygon Cj = cellS[j];

                if (Cj.IsOutside)
                    continue;

                //if (j == 125)
                //    Console.Write("");

                int NoEdgesB4 = Cj.Edges.Count;

                bool wasInside = Cj.IsInside;
                if (!wasInside) {
                    Cj.RemoveOutsideParts(true, IsIn);
                    if (Cj.Edges.Count <= 0 || (Cj.Edges.Count == 1 && Cj.Edges[0].isBoundary))
                        continue; // everything culled away
                }

                
                VoVertex[] seq = Cj.GetIntersectionSequence(bndyPoly);
                if(seq.Length > 0)
                    Insiders.Add(seq);
            }
            

           

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

                    int[,] iVtxTri;
                    try {
                        iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell); ;
                    } catch(ArithmeticException ae) {
                        DebugPlot(null, null, null, Insiders[jV]);
                        iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell); ;
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

                        if (D1.CrossProduct2D(D2) < 0) {
                            int it = iV0;
                            iV0 = iV2;
                            iV2 = it;

                            Vector vt = V0;
                            V0 = V2;
                            V2 = vt;

                            D1 = V1 - V0;
                            D2 = V2 - V0;
                        }

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


        static void DebugPlot(int[][] VocellVertexIndex, IList<Vector> Verts, IEnumerable<VoEdge> edgeS, IEnumerable<VoVertex> SomePoints) {
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
                if (SomePoints != null) {
                    PlotFormat ptF = new PlotFormat(Style: Styles.LinesPoints, pointType: PointTypes.Asterisk, lineColor: LineColors.Black);



                    double[] xS = SomePoints.Select(V => V.VTX.x).ToArray();
                    double[] yS = SomePoints.Select(V => V.VTX.y).ToArray();

                    gp.PlotXY(xS, yS, format: ptF);

                }

                gp.Execute();
                Console.WriteLine("Press any key to continue...");
                Console.ReadKey();
            }
        }



        static void CheckEdgeUniqueness(bool really = false) {
            VoEdge[] Allbndy = VoEdge.edgeS.ToArray();

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

    static class MySetHelpers {

        public static void MyRemove<T, L>(this L list, T itm) where L : IList<T> {
            int found = 0;
            for (int i = 0; i < list.Count; i++) {
                if(object.ReferenceEquals(itm, list[i])) {
                    found++;
                    list.RemoveAt(i);
                    i--;
                }
            }

            if (found != 1)
                throw new InvalidOperationException();
        }


        public static void MySetAdd<T, L>(this L list, T itm) where L : IList<T> {
            int found = 0;
            for (int i = 0; i < list.Count; i++) {
                if(itm.Equals(list[i])) {
                    found++;
                }
            }

            if (found > 1)
                throw new InvalidOperationException();
            if(found == 0) {
                list.Add(itm);
            }
        }
    }
}
