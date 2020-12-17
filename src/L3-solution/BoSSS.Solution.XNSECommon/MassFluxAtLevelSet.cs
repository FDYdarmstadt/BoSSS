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
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System.Collections;
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Solution.XheatCommon {


    public abstract class MassFluxAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        // for micro regions
        protected double m_sigma;

        //protected double rho;     // density of liquid phase 
        protected double m_rhoA;
        protected double m_rhoB;

        protected int m_D;


        public MassFluxAtLevelSet(int _D, LevelSetTracker _LsTrk, PhysicalParameters physicalParameters) {

            this.m_D = _D;
            this.m_LsTrk = _LsTrk;

            this.m_rhoA = physicalParameters.rho_A;
            this.m_rhoB = physicalParameters.rho_B;

            this.m_sigma = physicalParameters.Sigma;

        }       

        protected double MassFlux(CommonParams cp) {

            Debug.Assert(cp.Parameters_IN[0] == cp.Parameters_OUT[0], "mass flux must be continuous across interface");

            double M = 0.5 * (cp.Parameters_IN[0] + cp.Parameters_OUT[0]); // take average of massflux, should be equal in both phases anyway

            return M;
        }

        public abstract double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB);


        protected LevelSetTracker m_LsTrk;

        bool MEvapIsPrescribd = false;
        double prescrbMEvap;

        public virtual void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("prescribedMassflux")) {
                MEvapIsPrescribd = true;
                prescrbMEvap = (double)csA.UserDefinedValues["prescribedMassflux"];
            }

        }


        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return new List<string> { VariableNames.MassFluxExtension };
            }
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public virtual TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }


    }

    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class DivergenceAtLevelSet_withMassFlux : MassFluxAtLevelSet {

        public DivergenceAtLevelSet_withMassFlux(int _D, LevelSetTracker lsTrk,
            double vorZeichen, bool RescaleConti, PhysicalParameters physicalParameters)
            : base(_D, lsTrk, physicalParameters) {

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= m_rhoA;
                scaleB /= m_rhoB;
            }
        }

        double scaleA;
        double scaleB;


        public override double InnerEdgeForm(ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double M = MassFlux(cp);

            if (M == 0.0)
                return 0.0;

            double uAxN = -M * (1 / m_rhoA);
            double uBxN = -M * (1 / m_rhoB);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict;
            //uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            uAxN_fict = uBxN;

            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict;
            //uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            uBxN_fict = uAxN;


            // compute the fluxes: note that for the continuity equation, we use not a real flux,
            // but some kind of penalization, therefore the fluxes have opposite signs!
            double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side

            FlxNeg *= scaleA;
            FlxPos *= scaleB;

            double Ret = FlxNeg * vA - FlxPos * vB;

            return -Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }

    }

    /// <summary>
    /// 
    /// </summary>
    public class ViscosityAtLevelSet_FullySymmetric_withMassFlux : MassFluxAtLevelSet {


        public ViscosityAtLevelSet_FullySymmetric_withMassFlux(LevelSetTracker lstrk, double _penalty, int _component,
            PhysicalParameters physicalParameters)
            : base(lstrk.GridDat.SpatialDimension, lstrk, physicalParameters) {

            this.muA = physicalParameters.mu_A;
            this.muB = physicalParameters.mu_B;
            this.m_penalty_base = _penalty;
            this.component = _component;

        }

        double muA;
        double muB;

        int component;


        /// <summary>
        /// default-implementation
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams inp,
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

            double M = MassFlux(inp);
            if (M == 0.0)
                return 0.0;

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);


            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component];     // symmetry term
            Ret -= M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component] * (vA - vB) * pnlty * muMax; // penalty term
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

            double scaledPenalty = penaltySizeFactor * m_penalty * m_penalty_base;
            if (scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");

            return scaledPenalty;
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

    }

    public class MassFluxAtLevelSet_withMassFlux : MassFluxAtLevelSet {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        public MassFluxAtLevelSet_withMassFlux(int _d, int _D, LevelSetTracker LsTrk, PhysicalParameters physicalParameters)
            : base(_D, LsTrk, physicalParameters) {

            this.m_d = _d;
            if (m_d >= m_D)
                throw new ArgumentOutOfRangeException();
        }

        int m_d;


        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.Normal;

            double M = MassFlux(cp);
            if (M == 0.0)
                return 0.0;

            double massFlux = M.Pow2() * ((1 / m_rhoA) - (1 / m_rhoB)) * Normal[m_d];

            double FlxNeg = -0.5 * massFlux;
            double FlxPos = +0.5 * massFlux;

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            double Ret = FlxNeg * vA - FlxPos * vB;

            return Ret;
        }
    }

}
