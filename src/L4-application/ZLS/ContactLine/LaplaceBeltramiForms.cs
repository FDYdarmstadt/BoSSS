using BoSSS.Foundation;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using NSEVariableNames = BoSSS.Solution.NSECommon.VariableNames;

namespace ZwoLevelSetSolver.ContactLine {
    class FreeSlipContactLineForm : ContactLineForm {

        int d;
        int Dim;

        public FreeSlipContactLineForm(int d, int D) : base(D) {
            this.d = d;
            this.Dim = D;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);

            Vector Tangente_IN = SurfaceUtilities.Tangent(SurfaceNormal_IN, EdgeNormal);

            double Flx_InCell = 0;
            double m_sigma = Sigma(ref cpv);

            // isotropic surface tension terms
            for(int m_d = 0; m_d < this.Dim; m_d++) {
                Flx_InCell -= m_sigma * (EdgeNormal[m_d] * Tangente_IN[m_d]) * EdgeNormal[d];
            }
            return Flx_InCell * V;
        }
    }

    class NavierSlipLinearContactLineForm : ContactLineForm {

        int d;
        int D;
        double beta;
        double theta;

        public NavierSlipLinearContactLineForm(int d, int D, double beta, double theta) : base(D) {
            this.d = d;
            this.D = D;
            this.beta = beta;
            this.theta = theta;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);
            double m_sigma = Sigma(ref cpv);
            double Flx_InCell = NavierSlip(m_sigma, U, SurfaceNormal_IN, EdgeNormal);
            return Flx_InCell * V;
        }

        double NavierSlip(double m_sigma, Vector U, Vector SurfaceNormal_IN, Vector EdgeNormal) {
            Vector PSnI = new Vector(D); // projection of surface/level-set normal onto domain boundary tangent
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    double nn = EdgeNormal[d1] * EdgeNormal[d2];
                    if(d1 == d2) {
                        PSnI[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                    } else {
                        PSnI[d1] += -nn * SurfaceNormal_IN[d2];
                    }
                }
            }
            PSnI.Normalize();


            // Young's relation (static contact angle)
            double Flx_InCell = -m_sigma * Math.Cos(theta) * PSnI[d];

            // dissipative contact line force
            // beta*(u*nL)

            for(int d1 = 0; d1 < D; d1++) {
                Flx_InCell += beta * (U[d1] * PSnI[d1]) * PSnI[d];
            }
            return Flx_InCell;
        }
    }
}
