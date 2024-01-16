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

    /// <summary>
    /// Base class for fluxes that need to incorporate the mass flux at the level set
    /// </summary>
    public abstract class MassFluxAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        // for micro regions
        protected double m_sigma;

        //protected double rho;     // density of liquid phase 
        protected double m_rhoA;
        protected double m_rhoB;

        protected int m_D;


        public MassFluxAtLevelSet(int _D, PhysicalParameters physicalParameters, string phaseA, string phaseB) {

            this.m_D = _D;
            this.m_rhoA = physicalParameters.rho_A;
            this.m_rhoB = physicalParameters.rho_B;

            this.m_sigma = physicalParameters.Sigma;

            this.NegativeSpecies = phaseA;
            this.PositiveSpecies = phaseB;

        }       

        protected double MassFlux(CommonParams cp) {

            Debug.Assert(cp.Parameters_IN[0] == cp.Parameters_OUT[0], "mass flux must be continuous across interface");

            double M = 0.5 * (cp.Parameters_IN[0] + cp.Parameters_OUT[0]); // take average of massflux, should be equal in both phases anyway

            return M;
        }

        public abstract double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB);

        bool MEvapIsPrescribd = false;
        double prescrbMEvap;

        public virtual void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("prescribedMassflux")) {
                MEvapIsPrescribd = true;
                prescrbMEvap = (double)csA.UserDefinedValues["prescribedMassflux"];
            }

        }
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
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

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }

        public virtual TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }


    }

    /// <summary>
    /// Base class for fluxes that need to incorporate the mass flux at the level set
    /// </summary>
    public abstract class MassFluxAtLevelSet_StrongCoupling : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        // for micro regions
        protected double m_sigma;

        //protected double rho;     // density of liquid phase 
        protected double m_rhoA;
        protected double m_rhoB;
        protected double m_kA;
        protected double m_kB;
        protected double m_hvap;

        protected int m_D;


        public MassFluxAtLevelSet_StrongCoupling(int _D, ThermalParameters thermalParameters, string phaseA, string phaseB) {

            this.m_D = _D;
            this.m_rhoA = thermalParameters.rho_A;
            this.m_rhoB = thermalParameters.rho_B;

            this.m_kA = thermalParameters.k_A;
            this.m_kB = thermalParameters.k_B;

            this.m_hvap = thermalParameters.hVap;
            this.NegativeSpecies = phaseA;
            this.PositiveSpecies = phaseB;
        }

        protected double MassFlux(CommonParams cp, double[,] Grad_uA, double[,] Grad_uB) {

            double M = 0.0;
            if (!MEvapIsPrescribd) {
                for (int d = 0; d < m_D; d++) {
                    M += -(m_kA * Grad_uA[0, d] - m_kB * Grad_uB[0, d]) * cp.Normal[d];
                }
                M *= 1.0 / m_hvap;
            } else {
                M = prescrbMEvap;
            }
            //Console.WriteLine("Massflux: " + M);
            return M;
        }

        public abstract double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB);

        public bool MEvapIsPrescribd = false;
        double prescrbMEvap;

        public virtual void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("PrescribedMassFlux")) {
                MEvapIsPrescribd = true;
                prescrbMEvap = (double)csA.UserDefinedValues["PrescribedMassFlux"];
            }

        }

        public virtual IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return new string[0]; //VariableNames.Velocity0MeanVector(m_D);
            }
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }


        public virtual TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.GradUxV; }

        }
    }


    /// <summary>
    /// velocity jump penalty for the divergence operator (continuity equation), on the level set
    /// </summary>
    public class DivergenceAtLevelSet_withMassFlux : MassFluxAtLevelSet {

        public DivergenceAtLevelSet_withMassFlux(int _D,
            double vorZeichen, bool RescaleConti, PhysicalParameters physicalParameters, string phaseA, string phaseB)
            : base(_D, physicalParameters, phaseA, phaseB) {

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

            return  Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }
    }

    /// <summary>
    /// velocity jump penalty for the divergence operator (continuity equation), on the level set
    /// </summary>
    public class DivergenceAtLevelSet_Evaporation_StrongCoupling : MassFluxAtLevelSet_StrongCoupling {

        public DivergenceAtLevelSet_Evaporation_StrongCoupling(int _D,
            double vorZeichen, bool RescaleConti, ThermalParameters thermalParameters, string phaseA, string phaseB)
            : base(_D, thermalParameters, phaseA, phaseB) {

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

            double M = MassFlux(cp, Grad_uA, Grad_uB);

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

            return Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }
    }


    /// <summary>
    /// velocity jump penalty of the low mach equations for the divergence operator (continuity equation), on the level set 
    /// </summary>
    public class DivergenceAtLevelSet_Evaporation_StrongCoupling_LowMach : MassFluxAtLevelSet_StrongCoupling {

        public DivergenceAtLevelSet_Evaporation_StrongCoupling_LowMach(int _D,
            double vorZeichen, bool RescaleConti, ThermalParameters thermalParameters, string phaseA, string phaseB)
            : base(_D, thermalParameters, phaseA, phaseB) {

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

            double M = MassFlux(cp, Grad_uA, Grad_uB);

            if (M == 0.0)
                return 0.0;

            //double uAxN = -M * (1 / m_rhoA);
            //double uBxN = -M * (1 / m_rhoB);

            //// transform from species B to A: we call this the "A-fictitious" value
            //double uAxN_fict;
            ////uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            //uAxN_fict = uBxN;

            //// transform from species A to B: we call this the "B-fictitious" value
            //double uBxN_fict;
            ////uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            //uBxN_fict = uAxN;


            //// compute the fluxes: note that for the continuity equation, we use not a real flux,
            //// but some kind of penalization, therefore the fluxes have opposite signs!
            //double FlxNeg = +Flux(uAxN, uAxN_fict) * m_rhoA; // flux on A-side
            //double FlxPos = -Flux(uBxN_fict, uBxN) * m_rhoB;  // flux on B-side

            //FlxNeg *= scaleA;
            //FlxPos *= scaleB;

            //double Ret = FlxNeg * vA - FlxPos * vB;






            double uAxN = GenericBlas.InnerProd(U_Neg.GetSubVector(1, m_D), cp.Normal);
            double uBxN = GenericBlas.InnerProd(U_Pos.GetSubVector(1, m_D), cp.Normal);
            double res2 = 0.5 * (m_rhoA* uAxN + m_rhoB* uBxN) * (vA - vB); // This is correct WITH mass evaporation flux

            double Ret = -(M * 0.5 * (1.0 / m_rhoA + 1.0 / m_rhoB) * (m_rhoA - m_rhoB) - (m_rhoA - m_rhoB) * 0.5 * (uAxN + uBxN)) * 0.5 * (vA + vB);
            //double Ret = -(M * (1.0 / m_rhoA) * (m_rhoA - m_rhoB) - (m_rhoA - m_rhoB) * (uAxN)) * 0.5 * (vA + vB);
            return Ret+res2;
        }

        public override TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.AllOn; }
            //get { return TermActivationFlags.GradUxV | TermActivationFlags.V; }


        }
        public override IList<string> ArgumentOrdering => base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }
    }

    /// <summary>
    /// Correction terms at level set for a non material interface, for the viscous terms.
    /// </summary>
    public class ViscosityAtLevelSet_FullySymmetric_withMassFlux : MassFluxAtLevelSet {


        public ViscosityAtLevelSet_FullySymmetric_withMassFlux(double _penalty, int _component, int D,
            PhysicalParameters physicalParameters, string phaseA, string phaseB)
            : base(D, physicalParameters, phaseA, phaseB) {

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


            //double Grad_uA_xN = 0, Grad_uB_xN = 0;
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
            for (int i = 0; i < m_D; i++) {
                //Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                Ret += 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * N[component] * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[i];
            }
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component];     // symmetry term
            Ret -= M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component] * (vA - vB) * pnlty * muMax; // penalty term

            return Ret;
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

        public class ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling_LowMach : ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling {

            public ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling_LowMach(double _penalty, int _component, int D,
            ThermalParameters thermalParameters, PhysicalParameters physicalParameters, string phaseA, string phaseB)
            : base(_penalty,_component, D, thermalParameters,physicalParameters,phaseA,phaseB) {
            }

        public override TermActivationFlags LevelSetTerms => base.LevelSetTerms | TermActivationFlags.GradUxGradV | TermActivationFlags.V | TermActivationFlags.GradV; //  TermActivationFlags.V | TermActivationFlags.GradV; terms are necesary if the massflux is prescribed

    }

        /// <summary>
        /// Correction terms at level set for a non material interface, for the viscous terms.
        /// </summary>
        public class ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling : MassFluxAtLevelSet_StrongCoupling {


        public ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling(double _penalty, int _component, int D,
            ThermalParameters thermalParameters, PhysicalParameters physicalParameters, string phaseA, string phaseB)
            : base(D, thermalParameters, phaseA, phaseB) {

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


            //double Grad_uA_xN = 0, Grad_uB_xN = 0;
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

            double M = MassFlux(inp, Grad_uA, Grad_uB);

            if (M == 0.0)
                return 0.0;

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);


            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term   
            for (int i = 0; i < m_D; i++) {
                //Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                Ret += 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * N[component] * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[i];
            }
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component];     // symmetry term
            Ret -= M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component] * (vA - vB) * pnlty * muMax; // penalty term

            return Ret;
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


        //public override TermActivationFlags LevelSetTerms => base.LevelSetTerms | TermActivationFlags.GradUxGradV /*|TermActivationFlags.V*/;
        public override TermActivationFlags LevelSetTerms => base.LevelSetTerms | TermActivationFlags.GradUxGradV;

    }

    /// <summary>
    /// Interfaceflux Extension for jump conditions in Navierstokes with Massflux
    /// </summary>
    public class MassFluxAtLevelSet_withMassFlux : MassFluxAtLevelSet {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d"></param>
        /// <param name="_D"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        /// <param name="phaseA"></param>
        /// <param name="phaseB"></param>
        public MassFluxAtLevelSet_withMassFlux(int _d, int _D, PhysicalParameters physicalParameters, bool _movingMesh, string phaseA, string phaseB)
            : base(_D, physicalParameters, phaseA, phaseB) {

            this.m_d = _d;
            this.movingMesh = _movingMesh;
            if (m_d >= m_D)
                throw new ArgumentOutOfRangeException();
        }

        bool movingMesh;
        int m_d;

        /// <summary>
        /// 
        /// </summary>
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

            // moving-mesh-contribution
            // ========================

            if (movingMesh) {
                double s = ComputeInterfaceNormalVelocity(ref cp);
                Console.WriteLine("interface normal velocity = {0}", s);
                double movingFlux;
                if (M < 0) { // select DOWN-wind! (sign of M is equal to relative interface velocity (-s + u * n))
                    Console.WriteLine($"Velocity OUT: {cp.Parameters_OUT[1 + m_d]}");
                    movingFlux = (-s) * m_rhoB * cp.Parameters_OUT[1 + m_d]; // uB[0];
                } else {
                    Console.WriteLine($"Velocity IN: {cp.Parameters_IN[1 + m_d]}");
                    movingFlux = (-s) * m_rhoA * cp.Parameters_IN[1 + m_d]; // uA[0];
                }
                Console.WriteLine($"{Ret}   |   {movingFlux}");
                Ret -= movingFlux * Normal[m_d] * 0.5 * (vA + vB);
            }

            return Ret;
        }

        private double ComputeInterfaceNormalVelocity(ref CommonParams cp) {

            double M = MassFlux(cp);

            double sNeg = 0.0;
            for (int d = 0; d < m_D; d++)
                sNeg += (cp.Parameters_IN[1 + d]) * cp.Normal[d];
            sNeg -= (M / m_rhoA);

            double sPos = 0.0;
            for (int d = 0; d < m_D; d++)
                sPos += (cp.Parameters_OUT[1 + d]) * cp.Normal[d];
            sPos -= (M / m_rhoB);

            double s = (m_rhoA * sNeg + m_rhoB * sPos) / (m_rhoA + m_rhoB);     // density averaged, corresponding to the mean evo velocity 

            return s;
        }

        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0Vector(m_D));
            }
        }

    }

    /// <summary>
    /// Interfaceflux Extension for jump conditions in Navierstokes with Massflux
    /// </summary>
    public class MassFluxAtLevelSet_Evaporation_StrongCoupling : MassFluxAtLevelSet_StrongCoupling {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        public MassFluxAtLevelSet_Evaporation_StrongCoupling(int _d, int _D, ThermalParameters thermalParameters, bool _movingMesh, string phaseA, string phaseB)
            : base(_D, thermalParameters, phaseA, phaseB) {

            this.m_d = _d;
            this.movingMesh = _movingMesh;
            if (m_d >= m_D)
                throw new ArgumentOutOfRangeException();
        }

        bool movingMesh;
        int m_d;

        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.Normal;

            double M = MassFlux(cp, Grad_uA, Grad_uB);           

            if (M == 0.0)
                return 0.0;

            double massFlux = 0.0;

            // moving-mesh-contribution
            // ========================
            double Ret = 0.0;
            double FlxNeg = 0.0;
            double FlxPos = 0.0;

            if (movingMesh) {
                massFlux = M * 0.5 * (uA[1 + m_d]+ uB[1 + m_d]);

                FlxNeg = massFlux;
                FlxPos = massFlux;

                Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
                Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

                Ret = FlxNeg * vA - FlxPos * vB;
            } else {
                massFlux = M.Pow2() * ((1 / m_rhoA) - (1 / m_rhoB)) * Normal[m_d];

                FlxNeg = -0.5 * massFlux;
                FlxPos = +0.5 * massFlux;

                Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
                Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

                Ret = FlxNeg * vA - FlxPos * vB;
            }

            return Ret;
        }

        public override TermActivationFlags LevelSetTerms { 
            get {
                if (this.movingMesh) return base.LevelSetTerms | TermActivationFlags.UxV;
                else return base.LevelSetTerms;
            } 
        }

        public override IList<string> ArgumentOrdering => base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    /// <summary>
    /// Interfaceflux Extension for jump conditions in Navierstokes with Massflux,
    /// This is "Tangential Recoil pressure", only used when slip on interface, i.e., a velocity jump in tangential direction exists
    /// </summary>
    public class MassFluxAtLevelSet_Evaporation_StrongCoupling_Tangential : MassFluxAtLevelSet_StrongCoupling {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        public MassFluxAtLevelSet_Evaporation_StrongCoupling_Tangential(int _d, int _D, PhysicalParameters physicalParameters, ThermalParameters thermalParameters, bool _movingMesh, string phaseA, string phaseB)
            : base(_D, thermalParameters, phaseA, phaseB) {

            this.m_d = _d;
            this.movingMesh = _movingMesh;
            this.physParams = physicalParameters;
            if (m_d >= m_D)
                throw new ArgumentOutOfRangeException();
        }

        PhysicalParameters physParams;
        bool movingMesh;
        int m_d;

        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.Normal;

            double M = MassFlux(cp, Grad_uA, Grad_uB);

            if (M == 0.0)
                return 0.0;

            // moving-mesh-contribution
            // ========================
            double Ret = 0.0;
            double beta = physParams.slipI.IsInfinity() ? 0.0 : (physParams.mu_A + physParams.mu_B) / (2 * physParams.slipI);

            if (movingMesh) {
                throw new NotImplementedException();
            } else {
                var P = SurfaceProjection(Normal);
                for (int i = 0; i < m_D; i++) {
                    for (int j = 0; j < m_D; j++) {
                        Ret -= M * P[i, j] * (uA[1 + j] - uB[1 + j]) * 0.5 * P[i, m_d] * (vA + vB); // using the velocity jump directly
                    }
                }
                Debug.Assert(!(double.IsNaN(Ret) || double.IsInfinity(Ret)));
            }

            return Ret;
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        public override TermActivationFlags LevelSetTerms {
            get {
                if (this.movingMesh) return base.LevelSetTerms | TermActivationFlags.UxV;
                else return base.LevelSetTerms | TermActivationFlags.UxV;
            }
        }

        public override IList<string> ArgumentOrdering => base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    /// <summary>
    /// Interfaceflux Extension for jump conditions in Heatequation with Massflux
    /// </summary>
    public class HeatFluxAtLevelSet_Evaporation_StrongCoupling : MassFluxAtLevelSet_StrongCoupling {


        protected double m_cA;
        protected double m_cB;
        protected double m_Tsat;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        public HeatFluxAtLevelSet_Evaporation_StrongCoupling(int _D, ThermalParameters thermalParameters, bool _movingMesh, string phaseA, string phaseB)
            : base(_D, thermalParameters, phaseA, phaseB) {

            this.m_cA = thermalParameters.c_A;
            this.m_cB = thermalParameters.c_B;
            this.m_Tsat = thermalParameters.T_sat;
            this.movingMesh = _movingMesh;
            if (m_d >= m_D)
                throw new ArgumentOutOfRangeException();
        }

        bool movingMesh;
        int m_d;

        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.Normal;

            double M = MassFlux(cp, Grad_uA, Grad_uB);

            if (M == 0.0)
                return 0.0;

            double massFlux = 0.0;

            // moving-mesh-contribution
            // ========================
            double Ret = 0.0;
            double FlxNeg = 0.0;
            double FlxPos = 0.0;

            if (movingMesh) {
                massFlux = M * m_Tsat;

                FlxNeg = massFlux;
                FlxPos = massFlux;

                Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
                Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

                Ret = FlxNeg * vA * m_cA - FlxPos * vB * m_cB;
            } else { 
                Ret = 0.0;
            }

            return Ret;
        }

        public override TermActivationFlags LevelSetTerms {
            get {                
                return base.LevelSetTerms;
            }
        }

        public override IList<string> ArgumentOrdering => base.ArgumentOrdering;
    }

    /// <summary>
    /// Mass flux correction for the LLF flux in the convection terms for non material interface (only in case of splitting)
    /// </summary>
    public class ConvectionAtLevelSet_nonMaterialLLF_withMassFlux : MassFluxAtLevelSet {

        public ConvectionAtLevelSet_nonMaterialLLF_withMassFlux(int _d, int _D, PhysicalParameters physicalParameters, string phaseA, string phaseB)
            : base(_D, physicalParameters, phaseA, phaseB) {

            this.m_d = _d;
        }

        int m_d;


        public override double InnerEdgeForm(ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {



            double M = MassFlux(cp);
            if (M == 0.0)
                return 0.0;

            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = cp.Parameters_IN[1 + d];
                VelocityMeanOut[d] = cp.Parameters_OUT[1 + d];

            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);


            double uJump = -M * ((1 / m_rhoA) - (1 / m_rhoB)) * cp.Normal[m_d];


            double flx = Lambda * uJump * 0.8;

            return flx * (m_rhoA * vA - m_rhoB * vB);
        }



        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0MeanVector(m_D));
            }
        }
    }

    /// <summary>
    /// Mass flux correction for the LLF flux in the convection terms for non material interface (only in case of splitting)
    /// </summary>
    public class ConvectionAtLevelSet_nonMaterialLLF_withMassFlux_StrongCoupling : MassFluxAtLevelSet_StrongCoupling {

        public ConvectionAtLevelSet_nonMaterialLLF_withMassFlux_StrongCoupling(int _d, int _D, ThermalParameters thermParams, string phaseA, string phaseB)
            : base(_D, thermParams, phaseA, phaseB) {

            this.m_d = _d;
        }

        int m_d;


        public override double InnerEdgeForm(ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {



            double M = MassFlux(cp, Grad_uA, Grad_uB);
            if (M == 0.0)
                return 0.0;

            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = cp.Parameters_IN[d];
                VelocityMeanOut[d] = cp.Parameters_OUT[d];

            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);


            double uJump = -M * ((1 / m_rhoA) - (1 / m_rhoB)) * cp.Normal[m_d];


            double flx = Lambda * uJump * 0.8;

            return flx * (m_rhoA * vA - m_rhoB * vB);
        }



        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0MeanVector(m_D));
            }
        }
    }

    /// <summary>
    /// Same as <see cref="ConvectionAtLevelSet_Consistency_withMassFlux_StrongCoupling"/>, but using Newton Solver compatible fluxes
    /// </summary>
    public class ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton : MassFluxAtLevelSet_StrongCoupling {

        public ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton(int _d, int _D, PhysicalParameters physParams, ThermalParameters thermParams, string phaseA, string phaseB)
            : base(_D, thermParams, phaseA, phaseB) {

            this.m_d = _d;
            this.physParams = physParams;
        }

        int m_d;
        PhysicalParameters physParams;

        public override double InnerEdgeForm(ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {



            double M = MassFlux(cp, Grad_uA, Grad_uB);
            if (M == 0.0)
                return 0.0;

            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = cp.Parameters_IN[d]; //U_Neg[1+d];
                VelocityMeanOut[d] = cp.Parameters_OUT[d]; //U_Pos[1+d];

            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            double Lambda = 0.8 * Math.Max(LambdaIn, LambdaOut);

            double uJump = -M * ((1 / m_rhoA) - (1 / m_rhoB)) * cp.Normal[m_d];

            double slip = 0.0; // if there is slip, the tangential velocity is not continuous anymore!
            if(physParams.slipI != 0.0) {
                var P = SurfaceProjection(cp.Normal);
                double UxN = 0.0;
                for (int j = 0; j < m_D; j++) {
                    UxN += 0.5 * (U_Neg[j+1] + U_Pos[j+1]) * cp.Normal[j];
                }
                for (int i = 0; i < m_D; i++) {
                    slip -= P[m_d, i] * (U_Neg[i+1] - U_Pos[i+1]) * (Lambda*(m_rhoA * vA - m_rhoB * vB) - UxN * 0.5 * (m_rhoA * vA + m_rhoB * vB));
                }
            }

            double flx = Lambda * uJump;

            return flx * (m_rhoA * vA - m_rhoB * vB) + slip;
        }
        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        public override TermActivationFlags LevelSetTerms => base.LevelSetTerms | TermActivationFlags.UxV;

        public override IList<string> ArgumentOrdering {
            get {
                return base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    /// <summary>
    /// Mass flux correction for the LLF flux (central part) in the convection terms for non material interface (only in case of splitting)
    /// </summary>
    public class ConvectionAtLevelSet_Consistency_withMassFlux : MassFluxAtLevelSet {


        public ConvectionAtLevelSet_Consistency_withMassFlux(int _d, int _D, double vorZeichen, bool RescaleConti, PhysicalParameters physParams, string phaseA, string phaseB)
            : base(_D, physParams, phaseA, phaseB) {

            this.m_d = _d;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= m_rhoA;
                scaleB /= m_rhoB;
            }

        }

        int m_d;

        double scaleA;
        double scaleB;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }



        public override double InnerEdgeForm(ref CommonParams cp,

            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = MassFlux(cp);

            if (M == 0.0)
                return 0.0;

            double Ucentral = 0.0;
            for (int d = 0; d < m_D; d++) {
                Ucentral += 0.5 * (cp.Parameters_IN[1 + d] + cp.Parameters_OUT[1 + d]) * cp.Normal[d];
            }

            double uAxN = Ucentral * (-M * (1 / m_rhoA) * cp.Normal[m_d]);
            double uBxN = Ucentral * (-M * (1 / m_rhoB) * cp.Normal[m_d]);


            uAxN += -M * (1 / m_rhoA) * 0.5 * (U_Neg[0] + U_Pos[0]);
            uBxN += -M * (1 / m_rhoB) * 0.5 * (U_Neg[0] + U_Pos[0]);

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

            FlxNeg *= m_rhoA;
            FlxPos *= m_rhoB;

            double Ret = FlxNeg * vA - FlxPos * vB;

            return Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }



        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }


        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0Vector(m_D));
            }
        }

    }

    /// <summary>
    /// Mass flux correction for the LLF flux (central part) in the convection terms for non material interface (only in case of splitting)
    /// </summary>
    public class ConvectionAtLevelSet_Consistency_withMassFlux_StrongCoupling : MassFluxAtLevelSet_StrongCoupling {


        public ConvectionAtLevelSet_Consistency_withMassFlux_StrongCoupling(int _d, int _D,
            double vorZeichen, bool RescaleConti, ThermalParameters thermParams, string phaseA, string phaseB)
            : base(_D, thermParams, phaseA, phaseB) {

            this.m_d = _d;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= m_rhoA;
                scaleB /= m_rhoB;
            }

        }

        int m_d;

        double scaleA;
        double scaleB;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.GradUxV | TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }



        public override double InnerEdgeForm(ref CommonParams cp,

            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = MassFlux(cp, Grad_uA, Grad_uB);

            if (M == 0.0)
                return 0.0;

            double Ucentral = 0.0;
            double Ujump = 0.0;
            for (int d = 0; d < m_D; d++) {
                Ucentral += 0.5 * (cp.Parameters_IN[d] + cp.Parameters_OUT[d]) * cp.Normal[d];
                Ujump += (cp.Parameters_IN[d] - cp.Parameters_OUT[d]) * cp.Normal[d];
            }

            //double uAxN = Ucentral * (-M * (1 / m_rhoA) * cp.Normal[m_d]);
            //double uBxN = Ucentral * (-M * (1 / m_rhoB) * cp.Normal[m_d]);


            //uAxN += -M * (1 / m_rhoA) * 0.5 * (U_Neg[1] + U_Pos[1]);
            //uBxN += -M * (1 / m_rhoB) * 0.5 * (U_Neg[1] + U_Pos[1]);

            //// transform from species B to A: we call this the "A-fictitious" value
            //double uAxN_fict;
            ////uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            //uAxN_fict = uBxN;

            //// transform from species A to B: we call this the "B-fictitious" value
            //double uBxN_fict;
            ////uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            //uBxN_fict = uAxN;

            //// compute the fluxes: note that for the continuity equation, we use not a real flux,
            //// but some kind of penalization, therefore the fluxes have opposite signs!
            //double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            //double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side

            //FlxNeg *= m_rhoA;
            //FlxPos *= m_rhoB;

            //double Ret = FlxNeg * vA - FlxPos * vB;

            double Ret = (Ucentral * M * (1 / m_rhoA - 1 / m_rhoB) * cp.Normal[m_d] + 0.5 * (U_Neg[1] + U_Pos[1]) * Ujump) * 0.5 * (m_rhoA * vA + m_rhoB * vB);

            return Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }



        public override IList<string> ArgumentOrdering {
            get {
                return base.ArgumentOrdering.Cat(VariableNames.Velocity_d(m_d));
            }
        }


        public override IList<string> ParameterOrdering => base.ParameterOrdering.Cat(VariableNames.Velocity0Vector(m_D));

    }

    /// <summary>
    /// Same as <see cref="ConvectionAtLevelSet_Consistency_withMassFlux_StrongCoupling"/>, but using Newton solver compatible fluxes
    /// </summary>
    public class ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton : MassFluxAtLevelSet_StrongCoupling {


        public ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton(int _d, int _D,
            double vorZeichen, bool RescaleConti, ThermalParameters thermParams, string phaseA, string phaseB)
            : base(_D, thermParams, phaseA, phaseB) {

            this.m_d = _d;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= m_rhoA;
                scaleB /= m_rhoB;
            }

        }

        int m_d;

        double scaleA;
        double scaleB;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.GradUxV | TermActivationFlags.UxV;
            }
        }



        public override double InnerEdgeForm(ref CommonParams cp,

            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = MassFlux(cp, Grad_uA, Grad_uB);

            if (M == 0.0)
                return 0.0;

            double Ucentral = 0.0;
            double Ujump = 0.0;
            for (int d = 0; d < m_D; d++) {
                Ucentral += 0.5 * (U_Neg[1 + d] + U_Pos[1 + d]) * cp.Normal[d];
                Ujump += (U_Neg[1 + d] - U_Pos[1 + d]) * cp.Normal[d];
            }

            //double uAxN = Ucentral * (-M * (1 / m_rhoA) * cp.Normal[m_d]);
            //double uBxN = Ucentral * (-M * (1 / m_rhoB) * cp.Normal[m_d]);


            //uAxN += -M * (1 / m_rhoA) * 0.5 * (U_Neg[1] + U_Pos[1]);
            //uBxN += -M * (1 / m_rhoB) * 0.5 * (U_Neg[1] + U_Pos[1]);

            //// transform from species B to A: we call this the "A-fictitious" value
            //double uAxN_fict;
            ////uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            //uAxN_fict = uBxN;

            //// transform from species A to B: we call this the "B-fictitious" value
            //double uBxN_fict;
            ////uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            //uBxN_fict = uAxN;

            //// compute the fluxes: note that for the continuity equation, we use not a real flux,
            //// but some kind of penalization, therefore the fluxes have opposite signs!
            //double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            //double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side

            //FlxNeg *= m_rhoA;
            //FlxPos *= m_rhoB;

            //double Ret = FlxNeg * vA - FlxPos * vB;

            double Ret = (Ucentral * M * (1 / m_rhoA - 1 / m_rhoB) * cp.Normal[m_d] + 0.5 * (U_Neg[1 + m_d] + U_Pos[1 + m_d]) * Ujump) * 0.5 * (m_rhoA * vA + m_rhoB * vB);

            return Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }



        public override IList<string> ArgumentOrdering {
            get {
                return base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));
            }
        }
        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }

    }

    /// <summary>
    /// Mass flux correction for the convection terms and moving mesh terms for non material interface
    /// </summary>
    public class ConvectionAtLevelSet_MovingMesh_withMassFlux : MassFluxAtLevelSet {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        public ConvectionAtLevelSet_MovingMesh_withMassFlux(int _d, int _D, PhysicalParameters physicalParameters, string phaseA, string phaseB)
            : base(_D, physicalParameters, phaseA, phaseB) {

            this.m_d = _d;
            if (m_d >= m_D)
                throw new ArgumentOutOfRangeException();
        }

        int m_d;

        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.Normal;

            double M = MassFlux(cp);
            if (M == 0.0)
                return 0.0;

            // moving-mesh-contribution
            // ========================
            double Ret = 0.0;
            double movingFlux;

            movingFlux = M * 0.5 * (cp.Parameters_OUT[1 + m_d] + cp.Parameters_IN[1 + m_d]);
            Ret = movingFlux * (vA - vB);
            

            return Ret;
        }

        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0Vector(m_D));
            }
        }
    }   

}
