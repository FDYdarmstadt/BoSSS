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
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Solution.NSECommon.Operator.Convection {
    public class FSI_ConvectionAtIB : ILevelSetForm {
        public FSI_ConvectionAtIB(int currentDim, int spatialDim, LevelSetTracker LsTrk, IncompressibleBoundaryCondMap _bcmap, Func<Vector, double[]> getParticleParams, bool useMovingMesh) {
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

        /// <summary>
        /// Describes: 0: velX, 1: velY, 2: rotVel, 3: particle radius, 4: scaling paramete for active particles 
        /// </summary>
        private readonly Func<Vector, double[]> m_getParticleParams;
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

        public double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {

            BoSSS.Foundation.CommonParams inp;

            // Input parameters
            // =============================
            inp.Parameters_IN = cp.ParamsNeg;
            inp.Normale = cp.n;
            inp.iEdge = int.MinValue;
            inp.GridDat = this.m_LsTrk.GridDat;
            inp.X = cp.x;
            inp.time = cp.time;
            inp.Parameters_OUT = new double[inp.Parameters_IN.Length];

            // Particle parameters
            // =============================
            Vector X = new Vector(inp.X);
            double[] parameters_P = m_getParticleParams(X);
            double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
            double wLevSet = parameters_P[2];
            double[] RadialVector = new double[] { parameters_P[3], parameters_P[4] };
            double RadialLength = parameters_P[5];
            double active_stress = parameters_P[6];
            double Ang_P = parameters_P[7];
            double[] orientation = new double[] { Math.Cos(Ang_P), Math.Sin(Ang_P) };
            double scaleActiveBoundary = orientation[0] * inp.Normale[0] + orientation[1] * inp.Normale[1] > 0 && active_stress != 0 ? 1 : 0;

            // Level-set velocity
            // =============================
            double[] uLevSet_temp = new double[1];
            if (m_d == 0) {
                uLevSet_temp[0] = (uLevSet[0] - RadialLength * wLevSet * RadialVector[1]);
            }
            else {
                uLevSet_temp[0] = (uLevSet[1] + RadialLength * wLevSet * RadialVector[0]);
            }

            // Outer values for Velocity and VelocityMean
            // =============================
            inp.Parameters_OUT[0] = uLevSet[0] - RadialLength * wLevSet * RadialVector[1];
            inp.Parameters_OUT[1] = uLevSet[1] + RadialLength * wLevSet * RadialVector[0];
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
