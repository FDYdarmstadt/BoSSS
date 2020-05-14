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
using System.Collections;

namespace BoSSS.Solution.XheatCommon {

    /// <summary>
    /// 
    /// </summary>
    public class ViscosityAtLevelSet_FullySymmetric_withEvap : EvaporationAtLevelSet {


        public ViscosityAtLevelSet_FullySymmetric_withEvap(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component, 
            ThermalParameters thermParams, double _sigma) 
            : base(lstrk.GridDat.SpatialDimension, lstrk, thermParams, _sigma) {

            this.muA = _muA;
            this.muB = _muB;
            this.m_penalty_base = _penalty;
            this.component = _component;

        }

        double muA;
        double muB;

        int component;



        /// <summary>
        /// default-implementation
        /// </summary>
        public override double LevelSetForm(ref CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;
            //double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];


            //double Grad_uA_xN = 0, Grad_uB_xN = 0, 
            double Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }
            double Ret = 0.0;

            //double PosCellLengthScale = PosLengthScaleS[inp.jCellOut];
            //double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

            //double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            //if (hCutCellMin <= 1.0e-10 * hCellMin)
            //    // very small cell -- clippling
            //    hCutCellMin = hCellMin;

            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);


            double M = ComputeEvaporationMass(inp.Parameters_IN, inp.Parameters_OUT, N, inp.jCellIn);
            if (M == 0.0)
                return 0.0;

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);


            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component];     // symmetry term
            Ret -= M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component] * (vA - vB) * pnlty * muMax; // penalty term
                                                                                                               // Transpose Term
            for (int i = 0; i < m_D; i++) {
                //Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                Ret += 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * N[component] * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[i];
            }

            return -Ret;
        }



        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double Penalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];
            double penaltySizeFactor_B = 1.0 / PosLengthScaleS[jCellOut];

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;


        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        /// <param name="csA"></param>
        /// <param name="csB"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);

            double _D = m_D;
            double _p = TestDGdeg;

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.GradV | TermActivationFlags.V;
            }
        }

        //public override IList<string> ParameterOrdering {
        //    get {
        //        return ArrayTools.Cat(VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
        //    }
        //}


    }


}
