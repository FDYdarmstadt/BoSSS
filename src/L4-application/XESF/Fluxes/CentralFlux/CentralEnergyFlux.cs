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
    internal class CentralEnergyFlux : CentralFlux {
        public CentralEnergyFlux(IBoundaryConditionMap boundaryMap, Material material) : base(boundaryMap, material) {
        }

        public override Vector GetFluxVector(int D, double[] U) {
            Vector Output = new Vector(D);
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;


            double density = U[0];
            double energy = U[D + 1];

            double momentumSquared = 0.0;
            for(int d = 0; d < D; d++) {
                momentumSquared += U[d + 1] * U[d + 1];
            }
            double pressure = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density);

            //return state.Velocity * (state.Energy + state.Pressure);
            for(int d = 0; d < D; d++) {
                double flux = U[d + 1] / density * (energy + pressure);
                //Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                Output[d] += flux;
            }
            return Output;
        }

    }
}
