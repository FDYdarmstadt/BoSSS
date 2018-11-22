using BoSSS.Platform.LinAlg;
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

    }
}
