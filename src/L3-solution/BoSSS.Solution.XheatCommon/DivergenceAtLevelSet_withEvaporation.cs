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
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.XheatCommon {
    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class DivergenceAtLevelSet_withEvaporation : EvaporationAtLevelSet {

        public DivergenceAtLevelSet_withEvaporation(int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            double vorZeichen, bool RescaleConti, double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_LsTrk = lsTrk;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= rhoA;
                scaleB /= rhoB;
            }

            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.Rint = _Rint;

            this.Tsat = _Tsat;
            this.sigma = _sigma;
            this.pc = _pc;
        }


        double rhoA;
        double rhoB;

        double scaleA;
        double scaleB;



        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - GeneralizedDivergenceAtLevelSet: {0}", M);

            return M;

        }


        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.ParamsNeg, cp.ParamsPos, cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;


            double uAxN = -M * (1 / rhoA);
            double uBxN = -M * (1 / rhoB);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict;
            //uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            uAxN_fict = uBxN;

            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict;
            //uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            uBxN_fict = uAxN;


            // compute the fluxes: note that for the continuity equation, we use not a real flux,
            // but some kind of penalization, therefore the fluxes have opposite signs!
            double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side

            FlxNeg *= scaleA;
            FlxPos *= scaleB;

            double Ret = FlxNeg * vA - FlxPos * vB;

            return -Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }



        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D), VariableNames.Temperature, "Curvature", "DisjoiningPressure");
            }
        }

    }


}
