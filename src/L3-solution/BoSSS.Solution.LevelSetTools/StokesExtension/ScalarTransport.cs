using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {
    
    /// <summary>
    /// Flux for a scalar transport equation
    /// </summary>
    class ScalarTransportFlux : LinearFlux {

        public ScalarTransportFlux(IncompressibleBoundaryCondMap map, int D) {
            m_map = map;
            m_spatDim = D;
        }

        IncompressibleBoundaryCondMap m_map;

        int m_spatDim;

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            var EdgeType = m_map.EdgeTag2Type[inp.EdgeTag];
            Vector n = inp.Normal;
            var vel = ((Vector)inp.Parameters_IN);

            switch(EdgeType) {
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.Wall:
                return 0.0;
                //case IncompressibleBcType.Velocity_Inlet:
                //    Vector bndVel = new Vector(m_spatDim);
                //    var velFunction = m_spatDim.ForLoop(d => m_map.bndFunction[VariableNames.Velocity_d(d)]);
                //    for(int d = 0; d < m_spatDim; d++) {
                //        bndVel[d] = velFunction[d][inp.EdgeTag](inp.X, inp.time);
                //    }
                //    return (bndVel * Uin[0]) * n;
                //case IncompressibleBcType.Dong_OutFlow:
                //    if(n * vel >= 0) {
                //        // flow from inside 
                //        return (vel * Uin[0]) * n;
                //    } else {
                //        // flow from outside into the domain
                //        return 0.0;
                //    }
                default:
                if(n * vel >= 0) {
                    // flow from inside->outside, i.e., characteristic is leaving the domain, i.e., an open boundary
                    return (vel * Uin[0]) * n;
                } else {
                    // flow from outside into the domain
                    double Uout = Uin[0];
                    return (vel * Uout) * n;
                    //return (vel * Inflow(time)) * n;
                }
            }
        }

        /// <summary>
        /// An upwind flux
        /// </summary>
        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            Vector n = inp.Normal;

            var vel = 0.5 * (((Vector)inp.Parameters_IN) + ((Vector)inp.Parameters_OUT));

            if(vel * n > 0)
                return (vel * Uin[0]) * n;
            else
                return (vel * Uout[0]) * n;

        }

        /// <summary>
        /// `$ \underline{u} \varphi `$
        /// </summary>
        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            for(int d = 0; d < inp.D; d++) {
                output[d] = inp.Parameters[d] * U[0];
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[] { "Phi" }; }
        }

        /// <summary>
        /// the transport velocity
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return Solution.NSECommon.VariableNames.VelocityVector(this.m_spatDim);
            }
        }

    }
    
}
