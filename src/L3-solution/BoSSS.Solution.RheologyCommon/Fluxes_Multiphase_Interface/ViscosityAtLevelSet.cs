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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Solution.RheologyCommon {
    /// <summary>
    /// Viscosity part of constitutive equations at interface for multiphase.
    /// </summary>
    public class ViscosityAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        protected double betaA;
        protected double betaB;
        protected double[] pen1;

        /// <summary>
        /// Initialize viscosity part
        /// </summary>
        public ViscosityAtLevelSet(LevelSetTracker lstrk, int Component, double beta_a, double beta_b, double[] Penalty1) {
            this.m_LsTrk = lstrk;
            this.Component = Component;
            this.betaA = 1.0 - beta_a;
            this.betaB = 1.0 - beta_b;
            this.pen1 = Penalty1;
        }

        /// <summary>
        /// Ordering of the dependencies
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                switch (Component) {
                    case 0:
                        return new string[] { VariableNames.VelocityX, VariableNames.VelocityX };
                    case 1:
                        return new string[] { VariableNames.VelocityX, VariableNames.VelocityY };
                    case 2:
                        return new string[] { VariableNames.VelocityY, VariableNames.VelocityY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// Ordering of the parameters
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParams inp,
            double[] UA, double[] UB, double[,] Grad_uA, double[,] Grad_uB,
            double VA, double VB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];

            int D = N.Length;
            Debug.Assert(this.ArgumentOrdering.Count == 2);

            double PosCellLengthScale = PosLengthScaleS[inp.jCellOut];
            double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            Debug.Assert(UA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(UB.Length == this.ArgumentOrdering.Count);

            double res = 0;

            switch (Component) {
                case 0:
                    res += 0.5 * ((UA[0] * betaA + UB[0] * betaB) * N[0] + (UA[1] * betaA + UB[1] * betaB) * N[0]); // central difference fo grad(u) and grad(u)^T
                    res += pen1[0] * ((UA[0] * betaA - UB[0] * betaB) * N[0]) + pen1[1] * ((UA[0] * betaA - UB[0] * betaB) * N[1]); // beta Penalty for grad(u)
                    res += pen1[0] * ((UA[1] * betaA - UB[1] * betaB) * N[0]) + pen1[1] * ((UA[1] * betaA - UB[1] * betaB) * N[1]); // beta penalty for grad(u)^T

                    break;
                case 1:
                    res += 0.5 * ((UA[0] * betaA + UB[0] * betaB) * N[1] + (UA[1] * betaA + UB[1] * betaB) * N[0]);
                    res += pen1[0] * ((UA[0] * betaA - UB[0] * betaB) * N[0]) + pen1[1] * ((UA[0] * betaA - UB[0] * betaB) * N[1]);
                    res += pen1[0] * ((UA[1] * betaA - UB[1] * betaB) * N[0]) + pen1[1] * ((UA[1] * betaA - UB[1] * betaB) * N[1]);

                    break;
                case 2:
                    res += 0.5 * ((UA[0] * betaA + UB[0] * betaB) * N[1] + (UA[1] * betaA + UB[1] * betaB) * N[1]);
                    res += pen1[0] * ((UA[0] * betaA - UB[0] * betaB) * N[0]) + pen1[1] * ((UA[0] * betaA - UB[0] * betaB) * N[1]);
                    res += pen1[0] * ((UA[1] * betaA - UB[1] * betaB) * N[0]) + pen1[1] * ((UA[1] * betaA - UB[1] * betaB) * N[1]);

                    break;
                default:
                    throw new NotImplementedException();
            }

            return -2 * 0.5 * res * (VA - VB);
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }
    }
}
