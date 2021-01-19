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
    public class DivergenceAtIB : ILevelSetForm, ISupportsJacobianComponent {

        LevelSetTracker m_LsTrk;

        public DivergenceAtIB(int _D, LevelSetTracker lsTrk, int iLevSet, string FluidSpc, string SolidSpecies) {
            this.D = _D;
            this.m_LsTrk = lsTrk;
            this.LevelSetIndex = iLevSet;
            this.PositiveSpecies = lsTrk.GetSpeciesId(SolidSpecies);
            this.NegativeSpecies = lsTrk.GetSpeciesId(FluidSpc);
        }

        int D;
        
        double pRadius;
        


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double DirichletFlux(double UxN_in, double UxN_out) {
            return (UxN_in - UxN_out);
        }

        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            
            double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);

            //var parameters_P = m_getParticleParams(cp.X, cp.time);
            //double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
            //double wLevSet = parameters_P[2];
            //pRadius = parameters_P[3];

            double[] _uLevSet = new double[D];

            //_uLevSet[0] = uLevSet[0]+pRadius*wLevSet*-cp.Normal[1];
            //_uLevSet[1] = uLevSet[1] + pRadius * wLevSet * cp.Normal[0];

            double uBxN = GenericBlas.InnerProd(_uLevSet, cp.Normal);
          
            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict;
            uAxN_fict = uBxN;

            double FlxNeg = -DirichletFlux(uAxN, uAxN_fict); // flux on A-side
            //double FlxPos = 0;

            return FlxNeg * v_Neg;
        }

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
            get;
            private set;
        }

        public SpeciesId PositiveSpecies {
            get;
            private set;
        }

        public SpeciesId NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }
}
