/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using BoSSS.Solution.Utils;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using System.Collections.Generic;

namespace BoSSS.Application.ScalarTransport {

    /// <summary>
    /// flux fo an 2D scalar transport equation
    /// </summary>
    class ScalarTransportFlux : NonlinearFlux {
        
        double Inflow(double time) {
            return 0.0;
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
            Vector2D n; n.x = normal[0]; n.y = normal[1];

            var vel = FlowField(x, Uin, Uin);

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
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <returns></returns>
        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            Vector2D n; n.x = normal[0]; n.y = normal[1];

            var vel = FlowField(x, Uin, Uout);

            if(vel * n > 0)
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
            get { return new string[] { "u" }; }
        }

        /// <summary>
        /// the transport velocity
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return Solution.NSECommon.VariableNames.VelocityVector(2);
            }
        }

    }
    
    
    /// <summary>
    /// flux fo an 3D scalar transport equation
    /// </summary>
    class ScalarTransportFlux3D : NonlinearFlux {



        double Inflow(double time) {
            return 0.0;
        }

        /// <summary>
        /// predefined, div-free flow field
        /// </summary>
        double[] FlowField(double[] x, double[] Uin, double[] Uout) {
            temp_u[0] = 0.5*(Uin[1] + Uout[1]);
            temp_u[1] = 0.5*(Uin[2] + Uout[2]);
            temp_u[2] = 0.5*(Uin[3] + Uout[3]);
            return temp_u;
        }

        double[] temp_u = new double[3];

        static double Dot(double[] a, double[] b) {
            return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
        }


        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            //Vector2D n; n.x = normal[0]; n.y = normal[1];
            var u = FlowField(x, Uin, Uin);

            if (Dot(u,normal) >= 0) {
                return Dot(u, normal) * Uin[0]; // (c(x) * Uin[0]) * n;
            } else {
                return Dot(u, normal) * Inflow(time); //(c(x) * Inflow(time)) * n;
            }
        }

        /// <summary>
        /// a new comment
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <returns></returns>
        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            double[] u = FlowField(x, Uin, Uout);
                        
            if (Dot(u, normal) > 0)
                return Dot(u, normal) * Uin[0]; // (c(x) * Uin[0]) * n;
            else
                return Dot(u, normal) * Uout[0]; // (c(x) * Uout[0]) * n;

        }

        protected override void Flux(double time, double[] x, double[] U, double[] output) {
            var u = FlowField(x, U, U);

            output[0] = u[0] * U[0];
            output[1] = u[1] * U[0];
            output[2] = u[2] * U[0];
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
                return Solution.NSECommon.VariableNames.VelocityVector(3);
            }
        }

    }
}
