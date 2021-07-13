using BoSSS.Foundation;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ContactLine {
    class LaplaceBeltramiEquilibriumForm : ContactLineForm {
        int d;
        int D;

        public LaplaceBeltramiEquilibriumForm(int d, int D) : base(D) {
            this.d = d;
            this.D = D;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);

            Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

            double Flx_InCell = 0;
            double m_sigma = Sigma(ref cpv);

            // isotropic surface tension terms
            Flx_InCell -= m_sigma * EdgeNormal[d] * Tangente_IN[d];
            return Flx_InCell * V;
        }
    }

    class LaplaceBeltramiSlipForm : ContactLineForm {
        int d;
        int D;

        public LaplaceBeltramiSlipForm(int d, int D) : base(D) {
            this.d = d;
            this.D = D;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);

            Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

            double Flx_InCell = 0;
            double m_sigma = Sigma(ref cpv);

            // isotropic surface tension terms
            Flx_InCell -= m_sigma * EdgeNormal[d] * (Tangente_IN * EdgeNormal) * (-1);
            return Flx_InCell * V;
        }
    }
}