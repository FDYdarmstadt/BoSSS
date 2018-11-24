using BoSSS.Platform.LinAlg;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.SipPoisson {
    static class PolygonItersection {


        // Sutherland–Hodgman for clipping polygons
        //   clip polygon: convex
        //   subject polygon: arbitrary
        //
        // List outputList = subjectPolygon;   
        // for (Edge clipEdge in clipPolygon) do
        //    List inputList = outputList;
        //    outputList.clear();
        //    Point S = inputList.last;
        //    for (Point E in inputList) do
        //       if (E inside clipEdge) then
        //          if (S not inside clipEdge) then
        //             outputList.add(ComputeIntersection(S,E,clipEdge));
        //          end if
        //          outputList.add(E);
        //       else if (S inside clipEdge) then
        //          outputList.add(ComputeIntersection(S,E,clipEdge));
        //       end if
        //       S = E;
        //    done
        // done


        public static Vector[] SutherlandHodgmanClipping(AffineManifold[] ClipPolygon, Vector[] SubjectPolygon) {
            List<Vector> outputList = new List<Vector>(SubjectPolygon);

            foreach(AffineManifold clipEdge in ClipPolygon) {
                List<Vector> InputList = new List<Vector>(outputList);
                outputList.Clear();

                Vector S = InputList.Last();

                foreach (Vector E in InputList) {
                    bool E_inside = clipEdge.PointDistance(E) <= 0;
                    bool S_inside = clipEdge.PointDistance(S) <= 0;
                    if (E_inside) {
                        if(!S_inside) {
                            outputList.Add(ComputeIntersection(S, E, clipEdge));
                        }
                        
                        outputList.Add(E);
                    } else {

                        if (!S_inside) {
                            outputList.Add(ComputeIntersection(S, E, clipEdge));
                        }
                    }
                    S = E;
                }
            }
            

            return outputList.ToArray();
        }


        static Vector ComputeIntersection(Vector S, Vector E, AffineManifold clipEdge) {
            var SE = AffineManifold.FromPoints(S, E);
            var I = AffineManifold.Intersect2D(SE, clipEdge);

            //var D = S - E;
            //D.Normalize();
            //var N1 = new Vector(D.y, -D.x);
            //double inner = N1 * clipEdge.Normal;

            Debug.Assert(!double.IsNaN(I.x));
            Debug.Assert(!double.IsInfinity(I.x));
            Debug.Assert(!double.IsNaN(I.y));
            Debug.Assert(!double.IsInfinity(I.y));

            return I;
        }



        static bool ComputeIntersection(Vector S1, Vector S2, Vector E1, Vector E2, out double alpha1, out double alpha2, out Vector I) {
            if (S1.Dim != 2)
                throw new ArgumentException("spatial dimension mismatch.");
            if (S2.Dim != 2)
                throw new ArgumentException("spatial dimension mismatch.");
            if (E1.Dim != 2)
                throw new ArgumentException("spatial dimension mismatch.");
            if (E2.Dim != 2)
                throw new ArgumentException("spatial dimension mismatch.");

            Vector S12 = S2 - S1;
            Vector E12 = E2 - E1;
            if (S12.Abs() <= 0)
                throw new ArgumentException();
            if (E12.Abs() <= 0)
                throw new ArgumentException();


            var P_S12 = AffineManifold.FromPoints(S1, S2);
            var P_E12 = AffineManifold.FromPoints(E1, E2);

            Vector NS = P_S12.Normal; NS.Normalize();
            Vector NE = P_E12.Normal; NE.Normalize();

            double parallel = NS * NE;
            if (parallel.Abs() >= 1.0) {
                alpha1 = double.PositiveInfinity;
                alpha2 = double.PositiveInfinity;
                I = new Vector(double.PositiveInfinity, double.PositiveInfinity);
                return false;
            }

            //S12.Normalize();
            //E12.Normalize();

            I = AffineManifold.Intersect2D(P_S12, P_E12);

            Vector IS1 = I - S2;
            Vector IE1 = I - E2;
            Vector IS2 = I - S1;
            Vector IE2 = I - E1;

            Vector IS;
            bool flip_1;
            if(IS1.AbsSquare() > IS2.AbsSquare()) {
                IS = IS1;
                flip_1 = true;
            } else {
                IS = IS2;
                flip_1 = false;
            }

            Vector IE;
            bool flip_2;
            if(IE1.AbsSquare() > IE2.AbsSquare()) {
                IE = IE1;
                flip_2 = true;
            } else {
                IE = IE2;
                flip_2 = false;
            }

            Debug.Assert((S12.AngleTo(IS).Abs() <= 1.0e-5) || ((S12.AngleTo(IS).Abs() - Math.PI).Abs() <= 1.0e-5));
            Debug.Assert((E12.AngleTo(IE).Abs() <= 1.0e-5) || ((E12.AngleTo(IE).Abs() - Math.PI).Abs() <= 1.0e-5));

            alpha1 = (S12 * IS)/S12.AbsSquare();
            alpha2 = (E12 * IE)/E12.AbsSquare();

            if (flip_1)
                alpha1 = 1 + alpha1;
            if (flip_2)
                alpha2 = 1 + alpha2;

            return true;
        }


        struct MyTuple {
            public MyTuple(int __iSubj_1, int __iClip_2) { iSubj_1 = __iSubj_1; iClip_2 = __iClip_2; }
            public int iSubj_1;
            public int iClip_2;
        }

        public static Vector[] WeilerAthertonClipping(Vector[] _ClipReg, Func<Vector,bool> IsInClipReg, Vector[] _SubjPoly) {
            // =========================
            // checking the input
            // =========================

            if (_ClipReg.Length < 3)
                throw new ArgumentException();
            if (_SubjPoly.Length < 3)
                throw new ArgumentException();
            foreach (Vector v in _ClipReg) {
                if (v.Dim != 2)
                    throw new ArgumentException();
                if (!IsInClipReg(v))
                    throw new ArgumentException();
            }
            foreach (Vector v in _SubjPoly) {
                if (v.Dim != 2)
                    throw new ArgumentException();
            }
            
            
            // Given polygon A as the clipping region and polygon B as the subject polygon to be clipped, 
            // the algorithm consists of the following steps:
            //
            // 1. List the vertices of the clipping-region polygon A and those of the subject polygon B.
            // 2. Label the listed vertices of subject polygon B as either inside or outside of clipping region A.
            // 3. Find all the polygon intersections and insert them into both lists, linking the lists at the intersections.
            // 4. Generate a list of "inbound" intersections – the intersections where the vector from the intersection 
            //    to the subsequent vertex of subject polygon B begins inside the clipping region.
            // 5. Follow each intersection clockwise around the linked lists until the start position is found.
            //
            // If there are no intersections then one of three conditions must be true:
            // i.   A is inside B – return A for clipping, B for merging.
            // ii.  B is inside A – return B for clipping, A for merging.
            // iii. A and B do not overlap – return None for clipping or A & B for merging.
            //
            // (from https://en.wikipedia.org/wiki/Weiler%E2%80%93Atherton_clipping_algorithm)

            // =================
            // step 1:
            // =================

            List<Vector> ClipReg = _ClipReg.ToList();
            List<Vector> SubjPoly = _SubjPoly.ToList();
            _ClipReg = null; //  avoid accidental use
            _SubjPoly = null; // avoid accidental use
            
            // =================
            // step 2:
            // =================
            List<int> SubjInside = SubjPoly.Select(X => IsInClipReg(X) ? 1 : -1).ToList();

            if (SubjInside.All(b => b > 0))
                // case (ii)
                return SubjPoly.ToArray();

            if(SubjInside.All(b => b < 0)) {
                // case (iii)
                // for the moment, we ignore case (i)

                return null;
            }

            // =================
            // step 3:
            // =================

            List<MyTuple> links = new List<MyTuple>();

            for(int iEdgeSub = 1; iEdgeSub <= SubjPoly.Count; iEdgeSub++) {
                Vector V1 = SubjPoly[iEdgeSub - 1];
                Vector V2 = SubjPoly[iEdgeSub % SubjPoly.Count];
                
                // test against all edges 
                for(int iEdgeClip = 1; iEdgeClip <= ClipReg.Count; iEdgeClip++) {
                    Vector S1 = ClipReg[iEdgeClip - 1];
                    Vector S2 = ClipReg[iEdgeClip % ClipReg.Count];

                    bool Intersect = ComputeIntersection(V1, V2, S1, S2, out double alpha1, out double alpha2, out Vector I);
                    if(Intersect) {
                        if(alpha1 >= 0.0 && alpha1 <= 1.0 && alpha2 >= 0.0 && alpha2 <= 1.0) {

                            int iEdgeSubInsert;
                            if(alpha1 <= 0.0) {

                            } else if(alpha1 >= 1.0) {

                            } else {
                                iEdgeSubInsert = iEdgeSub;
                                SubjPoly.Insert(iEdgeSub, I);
                                iEdgeSub++;
                            }

                            
                            ClipReg.Insert(iEdgeClip, I);

                            for(int i = 0; i < links.Count; i++) {
                                if (links[i].iClip_2 >= iEdgeClip) {
                                    var a = links[i];
                                    a.iClip_2++;
                                    links[i] = a;
                                }
                            }
                            Debug.Assert(links.Count <= 0 || iEdgeSubInsert > links.Max(tt => tt.iSubj_1));
                            links.Add(new MyTuple(iEdgeSubInsert, iEdgeClip));
                            SubjInside.Insert(iEdgeSub, 0);

                            
                            iEdgeClip++;
                                                       
                        } else {
                            // somewhere outside
                        }
                    } else {
                        // parallel

                    }
                }
            }

            return SubjPoly.ToArray();

            // =================
            // step 4 & 5
            // =================

            Debug.Assert(SubjPoly.Count == SubjInside.Count);
            Debug.Assert(SubjInside.Where(b => b > 0).Count() >= 1);

            int iStart = -1;
            for(int i = 0; i < SubjInside.Count; i++) {
                if(SubjInside[i] > 0) {
                    iStart = i;
                    break;
                }
            }

            List<Vector> R = new List<Vector>();

            int iCurrent = iStart;
            bool bOnSub = true;
            //int iNext = -1;
            int linksPointer = 0;
            Debug.Assert(linksPointer >= links.Count || iCurrent < links[linksPointer].iSubj_1);

            do {
                int L;
                if (bOnSub) {
                    R.Add(SubjPoly[iCurrent]);
                    L = SubjPoly.Count;
                } else {
                    R.Add(ClipReg[iCurrent]);
                    L = ClipReg.Count;
                }

                int iNext = iCurrent;
                iNext++;
                if (iNext >= L)
                    iNext = 0;

                if(bOnSub) {
                    Debug.Assert(SubjInside[iNext] >= 0);
                    Debug.Assert(linksPointer >= links.Count || iCurrent < links[linksPointer].iSubj_1);

                    if(linksPointer < links.Count && iNext == links[linksPointer].iSubj_1) {
                        // switch to clipping polygon

                        Debug.Assert(SubjInside[iNext] == 0);
                        iNext = links[linksPointer].iClip_2;
                        bOnSub = false;
                        linksPointer++;
                    }


                } else {
                    Debug.Assert(linksPointer < links.Count); // we are on clipping poly; ther must be still one link that brings us back on subject poly
                    //Debug.Assert(iCurrent < links[linksPointer].iClip_2);

                    if(iNext == links[linksPointer].iClip_2) {
                        iNext = links[linksPointer].iSubj_1;
                        bOnSub = true;
                        linksPointer++;
                    }
                }

                iCurrent = iNext;

            } while (!bOnSub || iCurrent != iStart);



            // =================
            // return
            // =================

            return R.ToArray();

        }

    }
}
