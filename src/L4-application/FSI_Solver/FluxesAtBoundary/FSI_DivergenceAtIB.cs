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
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using FSI_Solver;
using ilPSP;

namespace BoSSS.Solution.NSECommon.Operator.Continuity {
    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class FSI_DivergenceAtIB : ILevelSetForm {

        public FSI_DivergenceAtIB(int _D, LevelSetTracker lsTrk, Func<Vector, FSI_ParameterAtIB> coupling) {
            D = _D;
            m_LsTrk = lsTrk;
            m_Coupling = coupling;
        }

        private readonly LevelSetTracker m_LsTrk;

        private readonly int D;

        /// <summary>
        /// Describes: 0: velX, 1: velY, 2:rotVel,3:particleradius
        /// </summary>
        private readonly Func<Vector, FSI_ParameterAtIB> m_Coupling;

        public double LevelSetForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            Vector fluidVelocity = new Vector(U_Neg);
            return (m_Coupling(inp.X).VelocityAtPointOnLevelSet() - fluidVelocity) * inp.Normal * v_Neg;
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
