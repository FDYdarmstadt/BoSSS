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

namespace BoSSS.Solution.NSECommon.Operator.Convection {
    public class FSI_ConvectionAtIB : ILevelSetForm {
        public FSI_ConvectionAtIB(int currentDim, int spatialDim, LevelSetTracker LsTrk, IncompressibleBoundaryCondMap _bcmap, Func<Vector, FSI_ParameterAtIB> getParticleParams, bool useMovingMesh) {
            m_LsTrk = LsTrk;
            m_D = spatialDim;
            m_d = currentDim;
            m_getParticleParams = getParticleParams;
            m_UseMovingMesh = useMovingMesh;
            NegFlux = new LinearizedConvection(spatialDim, _bcmap, currentDim);
        }

        private readonly LevelSetTracker m_LsTrk;
        private readonly int m_D;
        private readonly int m_d;
        private readonly Func<Vector, FSI_ParameterAtIB> m_getParticleParams;
        private readonly bool m_UseMovingMesh;

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

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
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
            FSI_ParameterAtIB coupling = m_getParticleParams(inp.X);
            Vector orientation = new Vector(Math.Cos(coupling.Angle()), Math.Sin(coupling.Angle()));
            double scaleActiveBoundary = orientation * new Vector(inp.Normal) > 0 && coupling.ActiveStress() != 0 ? 1 : 0;

            // Level-set velocity
            // =============================
            double[] uLevSet_temp = new double[1];
            if (m_d == 0) {
                uLevSet_temp[0] = coupling.VelocityAtPointOnLevelSet()[0];
            }
            else {
                uLevSet_temp[0] = coupling.VelocityAtPointOnLevelSet()[1];
            }

            // Outer values for Velocity and VelocityMean
            // =============================            
            inp.Parameters_OUT = new double[inp.Parameters_IN.Length];
            inp.Parameters_OUT[0] = coupling.VelocityAtPointOnLevelSet()[0];
            inp.Parameters_OUT[1] = coupling.VelocityAtPointOnLevelSet()[1];
            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
            inp.Parameters_OUT[2] = 0;
            inp.Parameters_OUT[3] = 0;

            // Computing Flux
            // =============================
            double FlxNeg = m_UseMovingMesh == true
                ? 0 // Moving mesh
                : (this.NegFlux.InnerEdgeForm(ref inp, U_Neg, uLevSet_temp, null, null, v_Neg, 0, null, null)) * (1 - scaleActiveBoundary);// Splitting
            return FlxNeg;
        }
    }
}
