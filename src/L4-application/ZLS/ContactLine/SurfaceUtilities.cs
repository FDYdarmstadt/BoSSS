using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ContactLine {
    static class SurfaceUtilities {

        public static double[,] SurfaceProjection(Vector Nsurf) {

            int D = Nsurf.Dim;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1.0 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0.0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        public static Vector Tangent(Vector Nsurf, Vector Nedge) {
            Debug.Assert(Nsurf.Dim == Nedge.Dim);

            int D = Nsurf.Dim;

            Vector tau = new Vector(D);
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            tau.NormalizeInPlace();
            return tau;
        }

        public static Vector Rotate(Vector Nsurf) {
            Debug.Assert(Nsurf.Dim == 2);

            int D = Nsurf.Dim;

            Vector tau = new Vector(D);
            tau[0] = -Nsurf[1];
            tau[1] = Nsurf[0];
            return tau;
        }
    }
}
