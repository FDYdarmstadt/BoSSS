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

using System.Collections.Generic;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Application.FSI_Solver;

namespace BoSSS.Solution.NSECommon.Operator.Continuity {
    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class FSI_DivergenceAtIB : ILevelSetForm {

        public FSI_DivergenceAtIB(int _D, LevelSetTracker lsTrk, Particle[] allParticles, double minGridLength) {
            D = _D;
            m_LsTrk = lsTrk;
            this.allParticles = allParticles;
            this.minGridLength = minGridLength;
        }

        private readonly LevelSetTracker m_LsTrk;

        private readonly int D;

        private readonly Particle[] allParticles;

        private readonly double minGridLength;

        public double InnerEdgeForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            Vector fluidVelocity = new Vector(U_Neg);

            Particle currentParticle = allParticles[0];
            for (int p = 0; p < allParticles.Length; p++) {
                bool containsParticle = allParticles[p].Contains(inp.X, 1.5 * minGridLength);
                if (containsParticle) {
                    currentParticle = allParticles[p];
                    break;
                }
            }
            Vector radialVector = currentParticle.CalculateRadialVector(inp.X);
            Vector particleVelocity = new Vector(currentParticle.Motion.GetTranslationalVelocity(0)[0] - currentParticle.Motion.GetRotationalVelocity(0) * radialVector[1],
                              currentParticle.Motion.GetTranslationalVelocity(0)[1] + currentParticle.Motion.GetRotationalVelocity(0) * radialVector[0]);

            return (particleVelocity - fluidVelocity) * inp.Normal * v_Neg;
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
