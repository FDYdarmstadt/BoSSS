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

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;

namespace BoSSS.Application.XdgTimesteppingTest {


    class TranportFlux_Bulk : LinearFlux {

        public Func<double[], double, double>[] Inflow;
        
        Vector2D FlowField(double[] x, double[] Uin, double[] Uot) {
            Vector2D u;
            u.x = 0.5 * (Uin[0] + Uot[0]);
            u.y = 0.5 * (Uin[1] + Uot[1]);
            return u;
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            Vector2D n; n.x = inp.Normale[0]; n.y = inp.Normale[1];

            var vel = FlowField(inp.X, inp.Parameters_IN, inp.Parameters_IN);

            if (n * vel >= 0) {
                // flow from inside 
                return (vel * Uin[0]) * n;
            } else {
                // flow from outside into the domain

                //return (vel * Uin[0]) * n;
                return (vel * Inflow[inp.EdgeTag](inp.X, inp.time)) * n;
            }
        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            Vector2D n; n.x = inp.Normale[0]; n.y = inp.Normale[1];

            var vel = FlowField(inp.X, inp.Parameters_IN, inp.Parameters_OUT);
            if (vel * n  > 0)
                return (vel * Uin[0]) * n;
            else
                return (vel * Uout[0]) * n;
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            Vector2D o;
            o = FlowField(inp.Xglobal, inp.Parameters, inp.Parameters) * U[0];
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
                return new string[] { "Vx", "Vy" };
            }
        }

    }



    class TransportFlux_Interface : ILevelSetComponent {

        LevelSetTracker m_LsTrk;

        Func<double[], double, double> m_NormalVel; 

        public TransportFlux_Interface(LevelSetTracker lstrk, Func<double[], double, double> NormalVel) {
            m_LsTrk = lstrk;
            m_NormalVel = NormalVel;
        }

        Vector2D FlowField(double[] x, double[] Uin, double[] Uot) {
            Vector2D u;
            u.x = 0.5 * (Uin[0] + Uot[0]);
            u.y = 0.5 * (Uin[1] + Uot[1]);
            return u;
        }

        public double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Vector2D V = FlowField(inp.x, inp.ParamsNeg, inp.ParamsPos);
            Vector2D N = new Vector2D(inp.n);

            double s = m_NormalVel(inp.x, inp.time);
            double RelSpeed = V * N - s;

            /*
            // static flux contribution
            // ========================
            double staticFlux;
            if (dir >= 0) { // select UP-wind!
                staticFlux = (V * N) * uA[0];
            } else {
                staticFlux = (V * N) * uB[0];
            }

            // moving-mesh-contribution
            // ========================
                        
            double movingFlux;
            if (dir > 0) { // select DOWN-wind!
                movingFlux = (-s) * uB[0];
            } else {
                movingFlux = (-s) * uA[0];
            }

            // return
            // ======

            return (staticFlux + movingFlux) * (vA - vB);
            */


            // Flux in moving frame
            // ====================

            double Flux;
            if (RelSpeed >= 0) { // UP-wind with respect to relative speed
                Flux = RelSpeed * uA[0];
            } else {
                Flux = RelSpeed * uB[0];
            }
            return Flux * (vA - vB);
        }


        

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
        }

        /// <summary>
        /// the transport velocity
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return new string[] { "Vx", "Vy" };
            }
        }

    }


}
