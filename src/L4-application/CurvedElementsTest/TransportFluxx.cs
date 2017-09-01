using System;
using BoSSS.Solution.Utils;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using System.Collections.Generic;

namespace CurvedElementsTest {

    /// <summary>
    /// flux fo an 2D scalar transport equation
    /// </summary>
    class ScalarTransportFlux : NonlinearFlux {

        ///// <summary>
        ///// wind direction
        ///// </summary>
        //Vector2D c = new Vector2D(0,1);

        double Inflow(double time) {
            return 0.0;
        }

        /// <summary>
        /// the predefined, div-free flow field
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        Vector2D FlowField(double[] x) {
            Vector2D c;
            c.x = x[1]; // whirl-wind
            c.y = -x[0];
            //c.x = 1.0;
            //c.y = 1.0;
            return c;
        }



        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            Vector2D n; n.x = normal[0]; n.y = normal[1];

            if (n * FlowField(x) >= 0) {
                // flow from inside 
                return (FlowField(x) * Uin[0]) * n;
            } else {
                // flow from outside into the domain
                return (FlowField(x) * Inflow(time)) * n;
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

            if (FlowField(x) * n > 0)
                return (FlowField(x) * Uin[0]) * n;
            else
                return (FlowField(x) * Uout[0]) * n;

        }

        protected override void Flux(double time, double[] x, double[] U, double[] output) {
            Vector2D o;
            o = FlowField(x) * U[0];
            output[0] = o.x;
            output[1] = o.y;
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[] { "u" }; }
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
        /// <param name="x">input; spatial vector</param>
        /// <param name="c">output; velocity vector of the flow field at <paramref name="x"/></param>
        void FlowField(double[] x, double[] c) {
            c[0] = 1;
            c[1] = 0;
            c[2] = 0;
        }

        double[] c = new double[3];

        static double Dot(double[] a, double[] b) {
            return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
        }


        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            //Vector2D n; n.x = normal[0]; n.y = normal[1];
            FlowField(x, c);

            if (Dot(c, normal) >= 0) {
                return Dot(c, normal) * Uin[0]; // (c(x) * Uin[0]) * n;
            } else {
                return Dot(c, normal) * Inflow(time); //(c(x) * Inflow(time)) * n;
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
            FlowField(x, c);

            if (Dot(c, normal) > 0)
                return Dot(c, normal) * Uin[0]; // (c(x) * Uin[0]) * n;
            else
                return Dot(c, normal) * Uout[0]; // (c(x) * Uout[0]) * n;

        }

        protected override void Flux(double time, double[] x, double[] U, double[] output) {
            FlowField(x, c);

            output[0] = c[0] * U[0];
            output[1] = c[1] * U[0];
            output[2] = c[2] * U[0];
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[] { "u" }; }
        }


    }
}
