using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Application.AdaptiveMeshRefinementTest {

        /// <summary>
    /// Flux for an 2D scalar transport equation
    /// </summary>
    class ScalarTransportFlux : NonlinearFlux {
        
        double Inflow(double time) {
            return 0.0;
        }

        /// <summary>
        /// The predefined, div-free flow field
        /// </summary>
        Vector2D FlowField(double[] x) {
            Vector2D u;
            u.x = x[1];
            u.y = -x[0];
            return u;
        }



        protected override double BorderEdgeFlux(double time, double[] X, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            Vector2D n; n.x = normal[0]; n.y = normal[1];

            var vel = FlowField(X);

            if(n * vel >= 0) {
                // flow from inside 
                return (vel * Uin[0]) * n;
            } else {
                // flow from outside into the domain
                return (vel * Uin[0]) * n;
                //return (vel * Inflow(time)) * n;
            }
        }

        /// <summary>
        /// a new comment
        /// </summary>
        /// <param name="time"></param>
        /// <param name="X"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <returns></returns>
        protected override double InnerEdgeFlux(double time, double[] X, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            Vector2D n; n.x = normal[0]; n.y = normal[1];

            var vel = FlowField(X);

            if(vel * n > 0)
                return (vel * Uin[0]) * n;
            else
                return (vel * Uout[0]) * n;

        }

        protected override void Flux(double time, double[] X, double[] U, double[] output) {
            Vector2D o;
            o = FlowField(X) * U[0];
            output[0] = o.x;
            output[1] = o.y;
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[] { "u" }; }
        }

        /// <summary>
        /// the transport velocity
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }

}
