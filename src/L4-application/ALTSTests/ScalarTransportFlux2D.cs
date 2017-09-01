/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled or executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using BoSSS.Solution.Utils;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;

namespace ALTSTests {

    /// <summary>
    /// Flux to a 2D scalar transport equation
    /// </summary>
    class ScalarTransportFlux2D : NonlinearFlux {

        double inflow;

        public ScalarTransportFlux2D(double inflow) {
            this.inflow = inflow;
        }

        /// <summary>
        /// the predefined, div-free flow field
        /// </summary>
        Vector2D FlowField(double[] x, double[] Uin, double[] Uot) {
            Vector2D u;
            u.x = 0.5 * (Uin[1] + Uot[1]);
            u.y = 0.5 * (Uin[2] + Uot[2]);
            return u;
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            Vector2D n;
            n.x = normal[0];
            n.y = normal[1];

            var vel = FlowField(x, Uin, Uin);

            if (n * vel >= 0) {
                // flow from inside 
                return (vel * Uin[0]) * n;
            } else {
                // flow from outside into the domain
                //return (vel * Uin[0]) * n;
                return (vel * inflow) * n;
            }
        }

        /// <summary>
        /// calculating the inner edge fluxes by using a first oder upwind scheme
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <returns></returns>
        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            Vector2D n;
            n.x = normal[0];
            n.y = normal[1];

            var vel = FlowField(x, Uin, Uout);

            if (vel * n > 0)
                return (vel * Uin[0]) * n;
            else
                return (vel * Uout[0]) * n;
        }

        protected override void Flux(double time, double[] x, double[] U, double[] output) {
            Vector2D o;
            o = FlowField(x, U, U) * U[0];
            output[0] = o.x;
            output[1] = o.y;
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "c" };
            }
        }

        /// <summary>
        /// the transport velocity
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return BoSSS.Solution.NSECommon.VariableNames.VelocityVector(2);
            }
        }
    }
}