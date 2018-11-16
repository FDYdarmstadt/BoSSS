using System;
using System.Collections.Generic;
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
    }
}
