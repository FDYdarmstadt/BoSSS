using BoSSS.Foundation;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XESF.Fluxes {
    public class RHDensityFlux : RHFlux {
        public RHDensityFlux(IBoundaryConditionMap boundaryMap, Material material) : base(boundaryMap, material) {
        }
        public override Vector GetFluxVector(int D, double[] U) {
            Vector fluxvector = new Vector(D);
            for(int d = 0; d < D; d++) {
                // state.Momentum
                double flux = U[d + 1];
                Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                fluxvector[d] = flux;
            }
            return fluxvector;
        }

    }
}
