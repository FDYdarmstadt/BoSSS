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
using ilPSP.Utils;
using BoSSS.Foundation;
using FSI_Solver;
using ilPSP;
using BoSSS.Application.FSI_Solver;

namespace BoSSS.Solution.NSECommon.Operator.Convection {
    public class FSI_ConvectionAtIB : ILevelSetForm {
        public FSI_ConvectionAtIB(int currentDim, int spatialDim, IncompressibleBoundaryCondMap _bcmap, Particle[] allParticles, bool useMovingMesh, double minGridLength) {
            //m_LsTrk = LsTrk;
            m_D = spatialDim;
            m_d = currentDim;
            this.allParticles = allParticles;
            m_UseMovingMesh = useMovingMesh;
            this.minGridLength = minGridLength;
            NegFlux = new LinearizedConvection(spatialDim, _bcmap, currentDim);
        }

        //private readonly LevelSetTracker m_LsTrk;
        private readonly int m_D;
        private readonly int m_d;
        private readonly Particle[] allParticles;
        private readonly bool m_UseMovingMesh;
        private readonly double minGridLength;

        // Use Fluxes as in Bulk Convection
        private readonly LinearizedConvection NegFlux;

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Velocity_d(m_d) }; }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {

            CommonParams inp = cp;

            // Particle parameters
            // =============================
            Particle currentParticle = allParticles[0];
            for (int p = 0; p < allParticles.Length; p++) {
                bool containsParticle = allParticles[p].Contains(inp.X, 1.5 * minGridLength);
                if (containsParticle) {
                    currentParticle = allParticles[p];
                    break;
                }
            }
            double angle = currentParticle.Motion.GetAngle(0);
            Vector orientation = new Vector(Math.Cos(angle), Math.Sin(angle));
            double scaleActiveBoundary = orientation * new Vector(inp.Normal) > 0 && currentParticle.ActiveStress != 0 ? 1 : 0;

            // Level-set velocity
            // =============================
            double[] uLevSet_temp = new double[1];
            Vector radialVector = currentParticle.CalculateRadialVector(inp.X);
            Vector particleVelocity = new Vector(currentParticle.Motion.GetTranslationalVelocity(0)[0] - currentParticle.Motion.GetRotationalVelocity(0) * radialVector[1],
                              currentParticle.Motion.GetTranslationalVelocity(0)[1] + currentParticle.Motion.GetRotationalVelocity(0) * radialVector[0]);

            if (m_d == 0) {
                uLevSet_temp[0] = particleVelocity[0];
            }
            else {
                uLevSet_temp[0] = particleVelocity[1];
            }

            // Outer values for Velocity and VelocityMean
            // =============================            
            inp.Parameters_OUT = new double[inp.Parameters_IN.Length];
            inp.Parameters_OUT[0] = particleVelocity[0];
            inp.Parameters_OUT[1] = particleVelocity[1];
            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
            inp.Parameters_OUT[2] = 0;
            inp.Parameters_OUT[3] = 0;

            // Computing Flux
            // =============================
            double FlxNeg = m_UseMovingMesh
                ? 0 // Moving mesh
                : (this.NegFlux.InnerEdgeForm(ref inp, U_Neg, uLevSet_temp, null, null, v_Neg, 0, null, null)) * (1 - scaleActiveBoundary);// Splitting
            return FlxNeg;
        }
    }
}
