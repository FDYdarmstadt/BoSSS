using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HFSISolver.ContactLine {
    class LaplaceBeltramiEquilibriumForm : ContactLineForm, ISpeciesFilter {
        int d;
        int D;
        double scale;
        string species;

        public LaplaceBeltramiEquilibriumForm(int d, int D, double scale, string species) : base(D) {
            this.d = d;
            this.D = D;
            this.scale = scale;
            this.species = species;
        }

        public string ValidSpecies => species;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);

            Vector Tangente_IN = SurfaceUtilities.Tangent(SurfaceNormal_IN, EdgeNormal);

            double m_sigma = Sigma(ref cpv);

            // isotropic surface tension terms
            double force = scale * m_sigma * Tangente_IN[d];

            return force * V;
        }
    }
}