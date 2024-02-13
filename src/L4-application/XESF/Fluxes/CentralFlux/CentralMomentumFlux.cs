using BoSSS.Foundation;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XESF.Fluxes {
    internal class CentralMomentumFlux : CentralFlux {
        public CentralMomentumFlux(IBoundaryConditionMap boundaryMap, Material material, int component) : base(boundaryMap, material) {
            this.component = component;
        }

        protected int component;
        public override Vector GetFluxVector(int D, double[] U) {
            Vector Output = new Vector(D);
            int L = U.Length;
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;


            double density = U[0];
            double energy = U[D + 1];

            //state.Momentum[MomentumComponent] * state.Velocity
            //    + 1/(gamma * Mach^2) * state.Pressure * ComponentVector;
            double momentumSquared = 0.0;
            for(int d = 0; d < D; d++) {
                Output[d] += U[component + 1] * U[d + 1] / density;
                momentumSquared += U[d + 1] * U[d + 1];
            }

            double flux = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density) / gammaMachSquared;
            //Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
            Output[component] += flux;
            return Output;
        }
    }
}

