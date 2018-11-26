using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.SipPoisson {

    /// <summary>
    /// polygon tesselation/triangulating algorithms, 
    /// used reduce 2D Voronoi meshes to triangular meshes.
    /// </summary>
    static class PolygonTesselation {

        static double sign(Vector p1, Vector p2, Vector p3) {
            if (p1.Dim != 2)
                throw new ArgumentException();
            if (p2.Dim != 2)
                throw new ArgumentException();
            if (p3.Dim != 2)
                throw new ArgumentException();

            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        }

        static bool PointInTriangle(Vector pt, Vector v1, Vector v2, Vector v3) {


            double d1, d2, d3;
            bool has_neg, has_pos;

            d1 = sign(pt, v1, v2);
            d2 = sign(pt, v2, v3);
            d3 = sign(pt, v3, v1);

            has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

            return !(has_neg && has_pos);
        }


        static bool IsEar(List<Vector> Polygon, int iTst) {
            int L = Polygon.Count;

            if (iTst < 0)
                throw new IndexOutOfRangeException();
            if (iTst >= L)
                throw new IndexOutOfRangeException();

            int iNxt = iTst + 1;
            if (iNxt >= L)
                iNxt -= L;
            int iPrv = iTst - 1;
            if (iPrv < 0)
                iPrv += L;

            Vector vPrv = Polygon[iPrv];
            Vector vTst = Polygon[iTst];
            Vector vNxt = Polygon[iNxt];

            for(int l = 0; l < L; l++) {
                if (l == iPrv || l == iNxt || l == iTst)
                    continue;

                if (!PointInTriangle(Polygon[l], vPrv, vNxt, vTst))
                    return false;
            }

            return true;
        }
    
        public static int[,] TesselatePolygon(IEnumerable<Vector> _Polygon) {
            List<Vector> Polygon = _Polygon.ToList();
            //var R = new List<ValueTuple<int, int, int>>();
            int[,] R = new int[Polygon.Count - 2, 3];

            int I = Polygon.Count - 2;
            for(int i = 0; i < I; i--) { // Theorem: every simple polygon decomposes into I-2 triangles

                int iEar = -1;
                for (int iTst = 0; iTst < Polygon.Count; iTst++) {
                    if(IsEar(Polygon, iEar)) {
                        iTst = iEar;
                        break;
                    }
                }

                if (iEar < 0)
                    throw new ArithmeticException("unable to find ear.");


                int iNxt = iEar + 1;
                if (iNxt >= Polygon.Count)
                    iNxt -= Polygon.Count;
                int iPrv = iEar - 1;
                if (iPrv < 0)
                    iPrv += Polygon.Count;

                R[i, 0] = iPrv;
                R[i, 1] = iEar;
                R[i, 2] = iNxt;
            }

            // return
            return R;
        }


    }
}
