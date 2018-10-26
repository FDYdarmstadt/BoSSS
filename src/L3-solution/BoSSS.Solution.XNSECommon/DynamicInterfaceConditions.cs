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


namespace BoSSS.Solution.XNSECommon.Operator.DynamicInterfaceConditions {


    public class PrescribedMassFlux : ILevelSetForm {


        LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public PrescribedMassFlux(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _M) {
            m_LsTrk = LsTrk;
            if(_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;

            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.M = _M;
        }

        int m_D;
        int m_d;

        double rhoA;
        double rhoB;
        double M;


        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.n;

            double massFlux = -M.Pow2() * ((1/rhoA) - (1/rhoB)) * Normal[m_d];

            double FlxNeg = -0.5 * massFlux;
            double FlxPos = +0.5 * massFlux;


            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }
    }



}
