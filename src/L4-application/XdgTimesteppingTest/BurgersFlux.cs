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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XdgTimesteppingTest {


    class BurgersFlux_Bulk : LinearFlux {

        public Func<double[], double, double> Inflow;

        public Vector Direction;

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            Vector n = new Vector(2); n.x = inp.Normale[0]; n.y = inp.Normale[1];

            double u0In = inp.Parameters_IN[0];
            double u0Ot = Inflow(inp.X, inp.time);
            double uIn = Uin[0];
            double uOt = u0Ot;

            return UpWind(Direction, n, u0In, u0Ot, uIn, uOt);
        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            Vector n = new Vector(2); n.x = inp.Normale[0]; n.y = inp.Normale[1];

            double u0In = inp.Parameters_IN[0];
            double u0Ot = inp.Parameters_OUT[0];
            double uIn = Uin[0];
            double uOt = Uout[0];

            return UpWind(Direction, n, u0In, u0Ot, uIn, uOt);
        }

        internal static double UpWind(Vector Direction, Vector n, double u0In, double u0Ot, double uIn, double uOt) {

            // shock speed of the Riemann problem
            double s_Riemann = 0.5 * (u0In + u0Ot);

            Vector CharakterisicVel = Direction * s_Riemann;
            if (CharakterisicVel * n > 0)
                return (Direction * (u0In * uIn * 0.5)) * n;
            else
                return (Direction * (u0Ot * uOt * 0.5)) * n;
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            double u0 = inp.Parameters[0];
            output[0] = Direction.x * (u0 * U[0] * 0.5);
            output[1] = Direction.y * (u0 * U[0] * 0.5);
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
                return new string[] { "u0" };
            }
        }
    }


    public class BurgersFlux_Interface : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        Func<double[], double, double> m_NormalVel;

        Vector m_Direction;

        public BurgersFlux_Interface(LevelSetTracker lstrk, Func<double[], double, double> NormalVel, Vector Direction) {
            m_LsTrk = lstrk;
            m_NormalVel = NormalVel;
            m_Direction = Direction;
        }
               

        public double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Vector N = new Vector(inp.n);
            if (Math.Abs(N * m_Direction - 1.0) >= 1.0e-8)
                throw new ArithmeticException("Normal Vector mismatch.");

            double u0In = inp.ParamsNeg[0];
            double u0Ot = inp.ParamsPos[0];
            double uIn = uA[0];
            double uOt = uB[0];

            // interface speed
            double s = m_NormalVel(inp.x, inp.time);

            // shock speed
            double sigma = 0.5 * (u0In + u0Ot);
            sigma = 3.0 / 2.0;

            // Speed of shock relative to interface speed
            double relShockSpeed = sigma - s;

            // upwind-selection
            double uUpwnd, u0Upwnd;
            if (relShockSpeed > 0) {
                uUpwnd = uIn;
                u0Upwnd = u0In;
            } else {
                uUpwnd = uOt;
                u0Upwnd = u0Ot;
            }

            // Flux in moving frame (since level-set normal is equal to pseudo-1D-direction 'm_direction', we can drop the normal vectors.
            uUpwnd = u0Upwnd;
            double Flux = (0.5 * u0Upwnd - s) * uUpwnd;

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
                return new string[] { "u0" };
            }
        }

    }


}

