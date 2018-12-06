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
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon.Operator.Viscosity {

    public class ActiveViscosityAtIB : BoSSS.Foundation.XDG.ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public ActiveViscosityAtIB(int _d, int _D, LevelSetTracker t, double penalty, Func<double, int, double> _PenaltyFunc, double _muA,
            Func<double[], double, double[]> getParticleParams) {

            this.m_penalty = penalty;
            this.m_PenaltyFunc = _PenaltyFunc;
            this.m_LsTrk = t;
            this.muA = _muA;
            this.component = _d;
            this.m_getParticleParams = getParticleParams;
            this.m_D = _D;
        }

        int component;
        int m_D;

        /// <summary>
        /// Describes: 0: velX, 1: velY, 2: rotVel, 3: particleradius, 4: active_stress, 5: first scaling parameter, 6: second scaling parameter, 7: current angle
        /// </summary>
        Func<double[], double, double[]> m_getParticleParams;

        /// <summary>
        /// Viscosity in species A
        /// </summary>
        double muA;

        double m_penalty;

        

        Func<double, int, double> m_PenaltyFunc;

        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;


            //Debug.Assert(!double.IsNaN(inp.PosCellLengthScale));
            //double hCutCellMin = Math.Min(inp.PosCellLengthScale, inp.NegCellLengthScale);
            //double hCutCellMin = inp.NegCellLengthScale; // for IBM, there is no positive species!
            //if (hCutCellMin <= 1.0e-10 * hCellMin)
            //    // very small cell -- clippling
            //    hCutCellMin = hCellMin;
            //double _penalty = penalty(hCutCellMin);
            double _penalty = m_PenaltyFunc(m_penalty, inp.jCell);

            int D = N.Length;


            // Particle parameters
            // ============================= 
            var parameters_P = m_getParticleParams(inp.x, inp.time);
            double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
            double wLevSet = parameters_P[2];

            double pRadius = parameters_P[3];//distance between current position and center of mass
            double active_stress = parameters_P[4];
            double scale = parameters_P[5];
            double scale_2 = parameters_P[6];
            double Ang_P = parameters_P[7];

            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);


            // Gradient of u and v 
            // ============================= 
            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[component, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
            }


            // Evaluate the complete velocity as a sum of translation and angular velocity
            // ============================= 
            double uAFict;
            double Ret = 0.0;
            double active_stress_visc;

            // 3D for IBM_Solver
            // ============================= 
            if (inp.x.Length == 3) {
                
                Ret -= Grad_uA_xN * (vA);                                     // consistency term
                Ret -= Grad_vA_xN * (uA[component] - 0) * (1 - scale);        // symmetry term
                Ret += _penalty * (uA[component] - 0) * (vA) * (1 - scale);   // penalty term

                Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
                return Ret * muA;
            }


            // 2D
            // ============================= 
            //Defining boundary conditions (no slip/slip)
            if (component == 0) {
                uAFict = (uLevSet[component] + pRadius * wLevSet * -inp.n[1]) * (1 - scale) + uA[component] * scale;
            } else {
                uAFict = (uLevSet[component] + pRadius * wLevSet * inp.n[0]) * (1 - scale) + uA[component] * scale;
            }

            //Defining active stress
            if (component == 0)
            {
                active_stress_visc = active_stress * Math.Abs(inp.n[1]) * Math.Cos(Ang_P) * scale;
            } else
            {
                active_stress_visc = active_stress * Math.Abs(inp.n[0]) * Math.Sin(Ang_P) * scale;
            }

            //Computing flux
            Ret -= Grad_uA_xN * (vA) * (1 - scale) * muA;                    // consistency term
            Ret -= Grad_vA_xN * (uA[component] - uAFict) * muA;              // symmetry term
            Ret += _penalty * (uA[component] - uAFict) * (vA) * muA;         // penalty term
            Ret += active_stress_visc * (vA);                                // active term (Neumann boundary condition)

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
    }
}