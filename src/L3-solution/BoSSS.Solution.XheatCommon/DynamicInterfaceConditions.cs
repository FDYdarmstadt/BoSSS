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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System.Collections;
using BoSSS.Solution.XheatCommon;

namespace BoSSS.Solution.XheatCommon {


    public class MassFluxAtInterface : EvaporationAtLevelSet {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        public MassFluxAtInterface(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, ThermalParameters thermParams, double _Rint, double _sigma) {
            //double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            m_LsTrk = LsTrk;
            if(_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.D = _D;
            this.m_d = _d;

            this.rhoA = _rhoA;
            this.rhoB = _rhoB;

            this.kA = thermParams.k_A;
            this.kB = thermParams.k_B;
            this.hVapA = thermParams.hVap_A;
            this.Rint = _Rint;

            this.Tsat = thermParams.T_sat;
            this.sigma = _sigma;
            this.pc = thermParams.pc;
        }

        int m_d;

        double rhoA;
        double rhoB;



        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - MassFluxAtInterface: {0}", M);

            return M;

        }


        public override double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.n;

            double M = ComputeEvaporationMass(cp.ParamsNeg, cp.ParamsPos, cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;

            double massFlux = M.Pow2() * ((1 / rhoA) - (1 / rhoB)) * Normal[m_d];
               
            double p_disp = cp.ParamsNeg[1];
            // augmented capillary pressure
            //double acp_jump = 0.0;
            //if(!double.IsNaN(p_disp))
            //    acp_jump = massFlux + p_disp;
            //else
            //    acp_jump = massFlux;


            double FlxNeg = -0.5 * massFlux;
            double FlxPos = +0.5 * massFlux;


            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            double Ret = FlxNeg * vA - FlxPos * vB;

            return Ret;
        }



        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat( new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D), VariableNames.Temperature, "Curvature", "DisjoiningPressure" ); //;
            }
        }


    }


}
