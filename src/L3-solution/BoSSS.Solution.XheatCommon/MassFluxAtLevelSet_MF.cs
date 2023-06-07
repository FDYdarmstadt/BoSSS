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
    /// Temperature gradient is calculated from mixture fraction field
    /// </summary>
    public abstract class MassFluxAtLevelSet_StrongCoupling_MF : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        // for micro regions
        protected double m_sigma;

        //protected double rho;     // density of liquid phase 
        protected double m_rhoA;
        protected double m_rhoB;
        protected double m_kA;
        protected double m_kB;
        protected double m_hvap;

        protected int m_D;


        public MassFluxAtLevelSet_StrongCoupling_MF(int _D, ThermalParameters thermalParameters, string phaseA, string phaseB) {

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
                    double TF0 = 1; 
                    double TO0 = 1;
                    double Q = 1;
                    double dTdZ = TF0 - TO0 - Q;
                    M += -(m_kA * Grad_uA[0, d] - m_kB * Grad_uB[0, d]) * cp.Normal[d]* dTdZ;
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
                return new string[] { VariableNames.MixtureFraction };
            }
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return new List<string> { };
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
    /// velocity jump penalty of the low mach equations for the divergence operator (continuity equation), on the level set 
    /// </summary>
    public class DivergenceAtLevelSet_Evaporation_StrongCoupling_LowMach_MF : MassFluxAtLevelSet_StrongCoupling_MF {

        public DivergenceAtLevelSet_Evaporation_StrongCoupling_LowMach_MF(int _D,
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

            // double uAxN = -M * (1 / m_rhoA);
            // double uBxN = -M * (1 / m_rhoB);

            // // transform from species B to A: we call this the "A-fictitious" value
            // double uAxN_fict;
            // //uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            // uAxN_fict = uBxN;

            // // transform from species A to B: we call this the "B-fictitious" value
            // double uBxN_fict;
            // //uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            // uBxN_fict = uAxN;


            // // compute the fluxes: note that for the continuity equation, we use not a real flux,
            // // but some kind of penalization, therefore the fluxes have opposite signs!
            // double FlxNeg = -Flux(uAxN, uAxN_fict) * -1 * m_rhoA; // flux on A-side
            // double FlxPos = +Flux(uBxN_fict, uBxN) * -1 * m_rhoB;  // flux on B-side

            // FlxNeg *= scaleA;
            // FlxPos *= scaleB;

            //double Ret = FlxNeg * vA - FlxPos * vB;
            double uAxN = GenericBlas.InnerProd(U_Neg.GetSubVector(1, m_D), cp.Normal);
            double uBxN = GenericBlas.InnerProd(U_Pos.GetSubVector(1, m_D), cp.Normal);
            double Ret = -(M * 0.5 * (1.0 / m_rhoA + 1.0 / m_rhoB) * (m_rhoA - m_rhoB) - (m_rhoA - m_rhoB) * 0.5 * (uAxN + uBxN)) * 0.5 * (vA + vB);
            double res2 = 0.5 * (m_rhoA * uAxN + m_rhoB * uBxN) * (vA - vB); // This is correct WITH mass evaporation flux


            return Ret + res2;
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
    /// Same as <see cref="ConvectionAtLevelSet_Consistency_withMassFlux_StrongCoupling"/>, but using Newton Solver compatible fluxes
    /// </summary>
    public class ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton_MF : MassFluxAtLevelSet_StrongCoupling_MF {

        public ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton_MF(int _d, int _D, ThermalParameters thermParams, string phaseA, string phaseB)
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
                VelocityMeanIn[d] = U_Neg[1 + d];
                VelocityMeanOut[d] = U_Pos[1 + d];

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

        public override TermActivationFlags LevelSetTerms => base.LevelSetTerms | TermActivationFlags.UxV;

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
    /// Same as <see cref="ConvectionAtLevelSet_Consistency_withMassFlux_StrongCoupling"/>, but using Newton solver compatible fluxes
    /// </summary>
    public class ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton_MF : MassFluxAtLevelSet_StrongCoupling_MF {


        public ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton_MF(int _d, int _D,
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
    /// Interfaceflux Extension for jump conditions in Navierstokes with Massflux
    /// </summary>
    public class MassFluxAtLevelSet_Evaporation_StrongCoupling_MF : MassFluxAtLevelSet_StrongCoupling_MF {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        public MassFluxAtLevelSet_Evaporation_StrongCoupling_MF(int _d, int _D, ThermalParameters thermalParameters, bool _movingMesh, string phaseA, string phaseB)
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
                massFlux = M * 0.5 * (uA[1 + m_d] + uB[1 + m_d]);

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
    /// Correction terms at level set for a non material interface, for the viscous terms.
    /// </summary>
    public class ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling_MF : MassFluxAtLevelSet_StrongCoupling_MF {


        public ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling_MF(double _penalty, int _component, int D,
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

        public override TermActivationFlags LevelSetTerms {
            get {
                var terms = base.LevelSetTerms | TermActivationFlags.GradUxGradV;
                if (MEvapIsPrescribd)
                    terms |= TermActivationFlags.V;
                return TermActivationFlags.AllOn;

            }
        }
        //public override TermActivationFlags LevelSetTerms => base.LevelSetTerms | TermActivationFlags.GradUxGradV /*|TermActivationFlags.V*/;

    }

}
