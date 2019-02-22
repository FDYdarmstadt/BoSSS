using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Platform
{
    public static class PolygonTesselation
    { 
        /// <summary>
        /// 
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <returns></returns>
        static double Sign(Vector p1, Vector p2, Vector p3)
        {
            if (p1.Dim != 2)
                throw new ArgumentException();
            if (p2.Dim != 2)
                throw new ArgumentException();
            if (p3.Dim != 2)
                throw new ArgumentException();

            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        }

        static bool PointInTriangle(Vector pt, Vector v1, Vector v2, Vector v3)
        {
            double d1, d2, d3;
            bool has_neg, has_pos;

            d1 = Sign(pt, v1, v2);
            d2 = Sign(pt, v2, v3);
            d3 = Sign(pt, v3, v1);

            has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

            return !(has_neg && has_pos);
        }

        static bool IsEar(List<Vector> Polygon, int iTst)
        {
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


            for (int l = 0; l < L; l++)
            {
                if (l == iPrv || l == iNxt || l == iTst)
                    continue;

                if (PointInTriangle(Polygon[l], vPrv, vNxt, vTst))
                    return false;
            }

            Vector D1 = vPrv - vTst;
            Vector D2 = vNxt - vTst;
            if (D1.AngleTo(D2) >= (179.9999999 / 180.0) * Math.PI)
            {
                //if (D1.CrossProduct2D(D2).Abs() < 1.0e-6)
                D1.AngleTo(D2);
                return false; // this is a hack to avoid very skinny triangles
            }
            return true;
        }

        public static int[,] TesselatePolygon(IEnumerable<Vector> _Polygon)
        {
            List<Vector> Polygon = _Polygon.ToList();
            List<int> orgIdx = Polygon.Count.ForLoop(i => i).ToList();

            //var R = new List<ValueTuple<int, int, int>>();
            int[,] R = new int[Polygon.Count - 2, 3];


            int I = Polygon.Count - 2;
            for (int i = 0; i < I; i++)
            { // Theorem: every simple polygon decomposes into (I-2) triangles

                int iEar = -1;
                for (int iTst = 0; iTst < Polygon.Count; iTst++)
                {
                    if (IsEar(Polygon, iTst))
                    {
                        iEar = iTst;
                        break;
                    }
                }

                if (iEar < 0)
                {
                    /*
                    using(var gp = new Gnuplot()) {
                        double[] xn = Polygon.Select(X => X.x).ToArray();
                        double[] yn = Polygon.Select(X => X.y).ToArray();
                        double[] xm = _Polygon.Select(X => X.x).ToArray();
                        double[] ym = _Polygon.Select(X => X.y).ToArray();

                        gp.PlotXY(xn, yn, "earless", new PlotFormat("-xb"));
                        gp.PlotXY(xm, ym, "orginal", new PlotFormat(":0r"));

                        gp.Execute();

                        Console.ReadKey();
                    }
                    */
                    throw new ArithmeticException("unable to find ear.");
                }

                int iNxt = iEar + 1;
                if (iNxt >= Polygon.Count)
                    iNxt -= Polygon.Count;
                int iPrv = iEar - 1;
                if (iPrv < 0)
                    iPrv += Polygon.Count;

                R[i, 0] = orgIdx[iPrv];
                R[i, 1] = orgIdx[iEar];
                R[i, 2] = orgIdx[iNxt];

                orgIdx.RemoveAt(iEar);
                Polygon.RemoveAt(iEar);
            }

            // return
            return R;
        }

        public static int[,] TesselateConvexPolygon(IEnumerable<Vector> _Polygon)
        {
            Vector[] Polygon = _Polygon.ToArray();
            int NV = _Polygon.Count();
            int[,] R = new int[NV - 2, 3];

            for (int iTri = 0; iTri < NV - 2; iTri++)
            { // loop over triangles of voronoi cell
                int iV0 = 0;
                int iV1 = iTri + 1;
                int iV2 = iTri + 2;
                R[iTri, 0] = iV0;
                R[iTri, 1] = iV1;
                R[iTri, 2] = iV2;

                /*
                Vector V0 = Polygon[iV0];
                Vector V1 = Polygon[iV1];
                Vector V2 = Polygon[iV2];

                Vector D1 = V1 - V0;
                Vector D2 = V2 - V0;
                if(!(D1.CrossProduct2D(D2) > 1.0e-8)) {
                    throw new ArithmeticException("Not a convex polygon.");
                }
                //if (D1.CrossProduct2D(D2) < 0) {
                //    Vector T = V2;
                //    int t = iV2;
                //    V2 = V1;
                //    iV2 = iV1;
                //    V1 = T;
                //    iV1 = t;
                //}
                */
            }

            return R;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_Polygon">
        /// A positively oriented polygon
        /// </param>
        /// <param name="point"></param>
        /// <returns></returns>
        public static bool PointInConvexPolygon(IEnumerable<Vector> _Polygon, Vector point)
        {
            //http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
            //Shift
            Vector[] shiftedPoly = new Vector[_Polygon.Count()];
            int i = 0;
            foreach (Vector vtx in _Polygon)
            {
                shiftedPoly[i] = vtx - point;
                ++i;
            }
            bool inside = true;
            //Check signs
            for(i = 0; i < shiftedPoly.Length - 1; ++i)
            {
                double a_i = shiftedPoly[i + 1].x * shiftedPoly[i].y - shiftedPoly[i].x * shiftedPoly[i + 1].y;
                inside &= a_i <= 0;
            }
            double a = shiftedPoly[0].x * shiftedPoly[shiftedPoly.Length-1].y 
                - shiftedPoly[shiftedPoly.Length-1].x * shiftedPoly[0].y;
            inside &= a <= 0;
            return inside;

        }
    }

    public static class PolygonClipping
    {
        /// <summary>
        /// Intersection of line <paramref name="S1"/>--<paramref name="S2"/> and <paramref name="E1"/>--<paramref name="E2"/>
        /// </summary>
        /// <param name="S1"></param>
        /// <param name="S2"></param>
        /// <param name="E1"></param>
        /// <param name="E2"></param>
        /// <param name="alpha1">
        /// coordinate of <paramref name="I"/> on the line <paramref name="S1"/>--<paramref name="S2"/>
        /// </param>
        /// <param name="alpha2">
        /// coordinate of <paramref name="I"/> on the line <paramref name="E1"/>--<paramref name="E2"/>
        /// </param>
        /// <param name="I"></param>
        /// <returns></returns>
        public static bool ComputeIntersection(Vector S1, Vector S2, Vector E1, Vector E2, out double alpha1, out double alpha2, out Vector I)
        {
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
            if (parallel.Abs() >= 1.0 - 1e-10)
            {
                alpha1 = P_S12.PointDistance(E1);
                //alpha1 = double.PositiveInfinity;
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
            if (IS1.AbsSquare() > IS2.AbsSquare())
            {
                IS = IS1;
                flip_1 = true;
            }
            else
            {
                IS = IS2;
                flip_1 = false;
            }

            Vector IE;
            bool flip_2;
            if (IE1.AbsSquare() > IE2.AbsSquare())
            {
                IE = IE1;
                flip_2 = true;
            }
            else
            {
                IE = IE2;
                flip_2 = false;
            }

            Debug.Assert((S12.AngleTo(IS).Abs() <= 1.0e-5) || ((S12.AngleTo(IS).Abs() - Math.PI).Abs() <= 1.0e-5));
            Debug.Assert((E12.AngleTo(IE).Abs() <= 1.0e-5) || ((E12.AngleTo(IE).Abs() - Math.PI).Abs() <= 1.0e-5));

            alpha1 = (S12 * IS) / S12.AbsSquare();
            alpha2 = (E12 * IE) / E12.AbsSquare();

            if (flip_1)
                alpha1 = 1 + alpha1;
            if (flip_2)
                alpha2 = 1 + alpha2;

            return true;
        }
    }
}
