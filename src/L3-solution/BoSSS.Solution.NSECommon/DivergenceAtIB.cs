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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon.Operator.Continuity {
    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class DivergenceAtIB : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public DivergenceAtIB(int _D, LevelSetTracker lsTrk,
            double vorZeichen, Func<double, double>[] _uLevSet, Func<double, double>[] _wLevSet, double particleRadius) {
            this.D = _D;
            this.uLevSet = _uLevSet;
            this.wLevSet = _wLevSet;
            this.m_LsTrk = lsTrk;
            this.pRadius = particleRadius;
        }

        int D;
        Func<double, double>[] uLevSet,wLevSet;
        double pRadius;

        /// <summary>
        /// the penalty flux
        /// </summary>
        static double DirichletFlux(double UxN_in, double UxN_out) {
            return (UxN_in - UxN_out);
        }

        public double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            
            double uAxN = GenericBlas.InnerProd(U_Neg, cp.n);

            double[] _uLevSet = new double[D];

            _uLevSet[0] = (uLevSet[0])(cp.time)+pRadius*wLevSet[0](cp.time)*-cp.n[1];
            _uLevSet[1] = (uLevSet[1])(cp.time) + pRadius * wLevSet[0](cp.time) * cp.n[0];

            double uBxN = GenericBlas.InnerProd(_uLevSet, cp.n);
          
            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict;
            uAxN_fict = uBxN;

            double FlxNeg = -DirichletFlux(uAxN, uAxN_fict); // flux on A-side
            //double FlxPos = 0;

            return FlxNeg * v_Neg;
        }

        /*
        public override void PrimalVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos) {
            FlxNeg = 0;
            FlxPos = 0;
        }

        public override void FluxPotential(out double G, double[] U) {
            G = 0;
        }

        public override void Nu(out double NuNeg, out double NuPos, ref CommonParamsLs cp) {
            NuNeg = 1.0;
            NuPos = 1.0;
        }
        */
        
        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(this.D);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public int LevelSetIndex {
            get {
                return 0;
            }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV;
            }
        }
    }
}
