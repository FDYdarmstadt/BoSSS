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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XdgTimesteppingTest {

    class HeatFlux_Bulk : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public void SetParameter(string speciesName, SpeciesId SpcId, MultidimensionalArray LengthScales) {
            switch (speciesName) {
                case "A": Viscosity = m_muA; rhs = m_rhsA; complementViscosity = m_muB; break;
                case "B": Viscosity = m_muB; rhs = m_rhsB; complementViscosity = m_muA; break;
                default: throw new ArgumentException("Unknown species.");
            }
            m_LengthScales = LengthScales;
        }

        MultidimensionalArray m_LengthScales;
        public double m_muA;
        public double m_muB;
        public double m_penalty = 4.0;
        public Func<double[], double, double> m_rhsA;
        public Func<double[], double, double> m_rhsB;

        double Viscosity;
        double complementViscosity;
        Func<double[], double, double> rhs;

        public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int d = 0; d < cpv.D; d++)
                acc += GradU[0, d] * GradV[d] * Viscosity;

            acc -= rhs(cpv.Xglobal, cpv.time) * V; 

            return acc;
        }



        public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            double muA = Viscosity;
            double muB = Viscosity;


            for (int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (muA * _Grad_uA[0, d] + muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }

            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term

            return -Acc;

        }

        double g_Neu(double[] X, double[] N, int EdgeTag) {
            return 0.0;
        }

        double g_Diri(double[] X, double time, int EdgeTag) {
            return 0.0;
        }

        bool IsDiri(double[] X, double[] N, int EdgeTag) {
            if (EdgeTag == 1)
                return true;
            if (EdgeTag == 2)
                return false;

            throw new NotSupportedException();
        }

        protected double penalty(int jCellIn, int jCellOut) {

            double muFactor; // the WTF factor
            if (jCellOut >= 0)
                muFactor = 1.0;
            else
                muFactor = Math.Max(Viscosity, complementViscosity) / Viscosity;
            double penaltySizeFactor_A = 1.0 / this.m_LengthScales[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / this.m_LengthScales[jCellOut] : 0;
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);
            return this.m_penalty * penaltySizeFactor * muFactor;  // Bem.: die Viskosität wird in der swipViscosity_Term1.InnerEdgeForm(...) dazumultipliziert.

            //return (base.penalty(jCellIn, jCellOut, cj)/currentMu)*Math.Max(currentMu, complementMu);
        }

        public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);
            double muA = Viscosity;

            if (IsDiri(inp.X, inp.Normale, inp.EdgeTag)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = g_Diri(inp.X, inp.time, inp.EdgeTag);

                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normale[d];
                    Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                    Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                }

                Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;

            } else {
                // inhom. Neumann
                // ++++++++++++
                double g_N = g_Neu(inp.X, inp.Normale, inp.EdgeTag);
                Acc += muA * g_N * _vA;
            }


            return -Acc;
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV | TermActivationFlags.V;
            }
        }
    }

    public class HeatFlux_Interface : ILevelSetComponent {

        LevelSetTracker m_LsTrk;

        public HeatFlux_Interface(LevelSetTracker lstrk, Func<double[], double, double> S) {
            this.m_LsTrk = lstrk;
            this.m_D = lstrk.GridDat.SpatialDimension;
            this.m_NormalVel = S;
        }

        public double m_muA;
        public double m_muB;
        public double m_penalty = 4.0;
        int m_D;
        Func<double[], double, double> m_NormalVel;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            // symmetric interior penalty
            // ==========================

            int D = N.Length;
            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double hCutCellMin = Math.Min(inp.NegCellLengthScale, inp.PosCellLengthScale);
            Debug.Assert(!(double.IsInfinity(hCutCellMin) || double.IsNaN(hCutCellMin)));

            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            double Ret = 0.0;
            Ret -= 0.5 * (m_muA * Grad_uA_xN + m_muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret -= 0.5 * (m_muA * Grad_vA_xN + m_muB * Grad_vB_xN) * (uA[0] - uB[0]);     // symmetry term
            Ret += (m_penalty / hCutCellMin) * (uA[0] - uB[0]) * (vA - vB) * (Math.Abs(m_muA) > Math.Abs(m_muB) ? m_muA : m_muB); // penalty term


            // moving-mesh-contribution
            // ========================

            double s = m_NormalVel(inp.x, inp.time);
            double movingFlux;
            if (s > 0) { // select DOWN-wind!
                movingFlux = (-s) * uB[0];
            } else {
                movingFlux = (-s) * uA[0];
            }


            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { "u" }; }
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

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }
}
