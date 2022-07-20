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
    /// Moving mesh consideration of convection and moving mesh terms in case of non material interface
    /// </summary>
    public class HeatConvectionAtLevelSet_MovingMesh_withMassflux : MassFluxAtLevelSet {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="physicalParameters"></param>
        /// <param name="_movingMesh"></param>
        public HeatConvectionAtLevelSet_MovingMesh_withMassflux(int _D, double _Tsat, PhysicalParameters physicalParameters, ThermalParameters thermalParameters, string phaseA, string phaseB)
            : base(_D, physicalParameters, phaseA, phaseB) {

            this.Tsat = _Tsat;
            c_A = thermalParameters.c_A;
            c_B = thermalParameters.c_B;
        }

        double c_A, c_B, Tsat;

        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double M = MassFlux(cp);
            if (M == 0.0)
                return 0.0;

            // moving-mesh-contribution
            // ========================
            double Ret = 0.0;
            double movingFlux;

            movingFlux = M * Tsat;
            Ret = movingFlux * (c_A * vA - c_B * vB);

            return Ret;
        }

    }    

    /// <summary>
    /// 
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_material: ILevelSetForm, ILevelSetEquationComponentCoefficient {

        //LevelSetTracker m_LsTrk;

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_material(int _D, double _capA, double _capB, double _LFFA, double _LFFB,
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh ) {

            m_D = _D;

            //m_LsTrk = LsTrk;

            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, _LFFA, double.NaN);
            NegFlux.SetParameter("A");
            PosFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, double.NaN, _LFFB);
            PosFlux.SetParameter("B");

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

        }

        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        // Use Fluxes as in Bulk Convection
        HeatConvectionInBulk NegFlux;
        HeatConvectionInBulk PosFlux;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }


        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;


            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);




            //Flux for negative side
            double FlxNeg, FlxPos, LLF;

            // NEGATIVE
            // central part
            double Tavg = U_Neg[0];
            FlxNeg = Tavg * (ParamsNeg[0] * cp.Normal[0] + ParamsNeg[1] * cp.Normal[1]);
            if (m_D == 3) {
                FlxNeg += Tavg * ParamsNeg[2] * cp.Normal[2];
            }

            // dissipative part
            double[] VelocityMeanIn = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = ParamsNeg[m_D + d];
            }
            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);
            double uJump = U_Neg[0] - U_Pos[0];
            //FlxNeg += uJump * LambdaIn * LFFA;
            //FlxNeg *= capA;

            // POSITIVE
            Tavg = U_Pos[0];
            FlxPos = Tavg * (ParamsPos[0] * cp.Normal[0] + ParamsPos[1] * cp.Normal[1]);
            if (m_D == 3) {
                FlxPos += Tavg * ParamsPos[2] * cp.Normal[2];
            }

            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanOut[d] = ParamsPos[m_D + d];
            }
            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);
            double Lambda = Math.Max(LambdaIn, LambdaOut);
            uJump = U_Neg[0] - U_Pos[0];
            //FlxPos += uJump * LambdaOut * LFFB;
            //FlxPos *= capB;

            //LLF = FlxNeg * v_Neg - FlxPos * v_Pos;
            LLF = 0.5 * (FlxPos + FlxNeg);
            LLF += uJump * Lambda;

            //LLF = 0.5 * (FlxNeg + FlxPos) * (capA * v_Neg - capB * v_Neg);

            if (movingmesh)
                return 0.0;
            else
                return LLF * (capA * v_Neg - capB * v_Pos);
        }


        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            this.NegFlux.CoefficientUpdate(csA, DomainDGdeg, TestDGdeg);
            this.PosFlux.CoefficientUpdate(csB, DomainDGdeg, TestDGdeg);

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        BitArray evapMicroRegion;


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_material_Newton : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent, IEquationComponentCoefficient {

        //LevelSetTracker m_LsTrk;

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_material_Newton(int _D, double _capA, double _capB, double _LFFA, double _LFFB,
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh, string phaseA, string phaseB) {

            m_D = _D;

            //m_LsTrk = LsTrk;

            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk_Newton(_D, _bcmap, _capA, _capB, _LFFA, double.NaN);
            NegFlux.SetParameter("A");
            PosFlux = new HeatConvectionInBulk_Newton(_D, _bcmap, _capA, _capB, double.NaN, _LFFB);
            PosFlux.SetParameter("B");

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

            this.NegativeSpecies = phaseA;
            this.PositiveSpecies = phaseB;
        }

        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        // Use Fluxes as in Bulk Convection
        HeatConvectionInBulk_Newton NegFlux;
        HeatConvectionInBulk_Newton PosFlux;

        /// <summary>
        /// Scale of this component, used for homotopy
        /// </summary>
        double Scale = 1.0;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }


        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            double FlxNeg;
            {
                BoSSS.Foundation.CommonParams inp = cp;
                inp.Parameters_OUT = ParamsNegFict;

                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            double FlxPos;
            {
                BoSSS.Foundation.CommonParams inp = cp;
                inp.Parameters_IN = ParamsPosFict;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }

            if (movingmesh)
                return 0.0;
            else
                return Scale * (FlxNeg * v_Neg - FlxPos * v_Pos);
        }


        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            this.NegFlux.CoefficientUpdate(csA, DomainDGdeg, TestDGdeg);
            this.PosFlux.CoefficientUpdate(csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        BitArray evapMicroRegion;


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature }.Cat(VariableNames.VelocityVector(m_D));
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
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

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_material_Newton_Hamiltonian : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent, IEquationComponentCoefficient {

        //LevelSetTracker m_LsTrk;

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_material_Newton_Hamiltonian(int _D, double _capA, double _capB, double _LFFA, double _LFFB,
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh, string phaseA, string phaseB) {

            m_D = _D;

            //m_LsTrk = LsTrk;

            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk_Newton(_D, _bcmap, _capA, _capB, _LFFA, double.NaN);
            NegFlux.SetParameter("A");
            PosFlux = new HeatConvectionInBulk_Newton(_D, _bcmap, _capA, _capB, double.NaN, _LFFB);
            PosFlux.SetParameter("B");

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

            this.NegativeSpecies = phaseA;
            this.PositiveSpecies = phaseB;

        }

        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        // Use Fluxes as in Bulk Convection
        HeatConvectionInBulk_Newton NegFlux;
        HeatConvectionInBulk_Newton PosFlux;

        /// <summary>
        /// Scale of this component, used for homotopy
        /// </summary>
        double Scale = 1.0;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }


        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double flx = 0.0;

            //===========================================================================================================
            //===========================================================================================================
            // First variant, using central flux for temperature            
            /*
            double FlxNeg = U_Neg[0] * (U_Neg[1] * cp.Normal[0] + U_Neg[2] * cp.Normal[1]);
            double FlxPos = U_Pos[0] * (U_Pos[1] * cp.Normal[0] + U_Pos[2] * cp.Normal[1]);
            if (m_D == 3) {
                FlxNeg += U_Neg[0] * U_Neg[3] * cp.Normal[2];
                FlxPos += U_Pos[0] * U_Pos[3] * cp.Normal[2];
            }

            // Term from partial integration back to strong form
            double sflx = capA * FlxNeg * v_Neg - capB * FlxPos * v_Pos;
            //funktioniert mit starker Form
            flx = (0.5 * (U_Neg[0] + U_Pos[0])) * ((U_Neg[1] * cp.Normal[0] + U_Neg[2] * cp.Normal[1]) * capA * v_Neg - (U_Pos[1] * cp.Normal[0] + U_Pos[2] * cp.Normal[1]) * capB * v_Pos) - sflx;
            */
            //===========================================================================================================
            //===========================================================================================================



            //===========================================================================================================
            //===========================================================================================================
            // Second variant using Roe-Type Scheme
            // Normal velocities
            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            double vINxN = 0.0, vOUTxN = 0.0;
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = U_Neg[1 + d];
                vINxN += VelocityMeanIn[d] * cp.Normal[d];
                VelocityMeanOut[d] = U_Pos[1 + d];
                vOUTxN += VelocityMeanOut[d] * cp.Normal[d];
            }

            vINxN *= capA;
            vOUTxN *= capB;
            double uJump = U_Neg[0] - U_Pos[0];

            flx = 0.5 * (Math.Min(vINxN, vOUTxN) - Math.Abs(Math.Min(vINxN, vOUTxN))) * -uJump * v_Neg;
            flx += 0.5 * (Math.Max(vINxN, vOUTxN) + Math.Abs(Math.Max(vINxN, vOUTxN))) * -uJump * v_Pos;
            //===========================================================================================================
            //===========================================================================================================

            if (movingmesh) {
                return 0.0;
            } else {
                return Scale * flx;
            }
        }


        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            this.NegFlux.CoefficientUpdate(csA, DomainDGdeg, TestDGdeg);
            this.PosFlux.CoefficientUpdate(csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        BitArray evapMicroRegion;


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature }.Cat(VariableNames.VelocityVector(m_D));
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
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

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }
    }

    /// <summary>
    /// Extension for dissipative part of LLF to ensure $T=T_sat$ at the interface
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_withMassflux : MassFluxAtLevelSet {

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_withMassflux(int _D, LevelSetTracker LsTrk, double _capA, double _capB, double _LFFA, double _LFFB,
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh, double _Tsat, PhysicalParameters _physicalParameters, string phaseA, string phaseB) : base(_D, _physicalParameters, phaseA, phaseB) {

            m_D = _D;

            //m_LsTrk = LsTrk;

            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, _LFFA, double.NaN);
            NegFlux.SetParameter("A");
            PosFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, double.NaN, _LFFB);
            PosFlux.SetParameter("B");


            //DirichletCond = _DiriCond;
            Tsat = _Tsat;

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

        }

        //bool MaterialInterface;
        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        //bool DirichletCond;
        double Tsat;

        // Use Fluxes as in Bulk Convection
        HeatConvectionInBulk NegFlux;
        HeatConvectionInBulk PosFlux;



        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }


        public override double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;


            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);





            // Calculate dissipative part
            // ==========================
            double[] VelocityMeanIn = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = ParamsNeg[1 + m_D + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);

            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanOut[d] = ParamsPos[1 + m_D + d];
            }

            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);


            double Flx = 0;
            //Flx += Lambda * (U_Pos[0] - Tsat) * capA * v_Neg - Lambda * (Tsat - U_Neg[0]) * capB * v_Pos;

            //double FlxNeg, FlxPos;
            //{
            //    double r = 0.0;


            //    // Calculate dissipative part
            //    // ==========================

            //    double[] VelocityMeanIn = new double[m_D];
            //    for (int d = 0; d < m_D; d++) {
            //        VelocityMeanIn[d] = ParamsNeg[m_D + d];
            //    }

            //    double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);

            //    double uJump = U_Neg[0] - Tsat;

            //    r += LambdaIn * uJump * LFFA;

            //    FlxNeg = capA * r;

            //}

            ////Flux for positive side
            //{
            //    double r = 0.0;

            //    // Calculate dissipative part
            //    // ==========================

            //    double[] VelocityMeanOut = new double[m_D];
            //    for (int d = 0; d < m_D; d++) {
            //        VelocityMeanOut[d] = ParamsPos[m_D + d];
            //    }


            //    double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            //    double uJump = Tsat - U_Pos[0];

            //    r += LambdaOut * uJump * LFFB;

            //    FlxPos = capB * r;

            //}

            //double Flx = FlxNeg * v_Neg - FlxPos * v_Pos;

            // contribution due to massflux
            {
                double M = MassFlux(cp);
                Flx += Tsat * M * (1 / base.m_rhoA - 1 / base.m_rhoB) * 0.5 * (capA * v_Neg + capB * v_Pos);

                //double FlxNeg = Tsat * (ParamsNeg[0] * cp.Normal[0] + ParamsNeg[1] * cp.Normal[1]);
                //if (m_D == 3) {
                //    FlxNeg += Tsat * ParamsNeg[2] * cp.Normal[2];
                //}
                //double FlxPos = Tsat * (ParamsPos[0] * cp.Normal[0] + ParamsPos[1] * cp.Normal[1]);
                //if (m_D == 3) {
                //    FlxPos += Tsat * ParamsPos[2] * cp.Normal[2];
                //}
                //Flx -= 0.5 * (FlxNeg + FlxPos) * (capA * v_Neg - capB * v_Pos);

                //double FlxNeg = 0, FlxPos = 0;
                //for(int d = 0; d < m_D; d++) {
                //    FlxNeg += ParamsNeg[1 + d] * cp.Normal[d];
                //    FlxPos += ParamsPos[1 + d] * cp.Normal[d];
                //}
                //FlxNeg *= U_Neg[0];
                //FlxPos *= U_Pos[0];
                //Flx += (FlxNeg - FlxPos) * 0.5 * (capA * v_Neg + capB * v_Pos);

            }
            if (movingmesh)
                return 0.0;
            else
                return Flx;
        }


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            this.NegFlux.CoefficientUpdate(csA, DomainDGdeg, TestDGdeg);
            this.PosFlux.CoefficientUpdate(csB, DomainDGdeg, TestDGdeg);
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        BitArray evapMicroRegion;


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return base.ParameterOrdering.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
    }

    /// <summary>
    /// Extension for dissipative part of LLF to ensure $T=T_sat$ at the interface
    /// </summary>
    public class HeatConvectionAtLevelSet_Direct : ILevelSetForm {

        bool movingmesh;
        string m_spcId;
        public HeatConvectionAtLevelSet_Direct(int _D, LevelSetTracker LsTrk, double _capA, double _capB,
            double _Tsat, PhysicalParameters _physicalParameters, string spcId)  {

            m_D = _D;

            m_spcId = spcId;

            Tsat = _Tsat;

            capA = _capA;
            capB = _capB;
        }

        int m_D;

        double capA;
        double capB;

        double Tsat;     


        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double Flx = 0;
            if(m_spcId == this.NegativeSpecies) {
                Flx += cp.Parameters_IN[0] * cp.Normal[0] + cp.Parameters_IN[1] * cp.Normal[1];
                if(m_D == 3)
                    Flx += cp.Parameters_IN[2] * cp.Normal[2];
                Flx *= U_Neg[0];
                Flx *= capA * v_Neg;
            } else {
                Flx -= cp.Parameters_OUT[0] * cp.Normal[0] + cp.Parameters_OUT[1] * cp.Normal[1];
                if (m_D == 3)
                    Flx -= cp.Parameters_OUT[2] * cp.Normal[2];
                Flx *= U_Pos[0];
                Flx *= capB * v_Pos;
            }
            return Flx;
        }  

       

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public int LevelSetIndex => 0;

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }
    }

    /// <summary>
    /// Extension for dissipative part of LLF to ensure $T=T_sat$ at the interface
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_withMassflux_StrongCoupling : MassFluxAtLevelSet_StrongCoupling, IEquationComponentCoefficient {

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_withMassflux_StrongCoupling(int _D, double _capA, double _capB, double _LFFA, double _LFFB, bool _movingmesh, double _Tsat, ThermalParameters thermParams, string phaseA, string phaseB) : base(_D, thermParams, phaseA, phaseB) {

            m_D = _D;

            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;
            
            //DirichletCond = _DiriCond;
            Tsat = _Tsat;

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

        }

        //bool MaterialInterface;
        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        //bool DirichletCond;
        double Tsat;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }

        /// <summary>
        /// Scale of the component, used for homotopy
        /// </summary>
        double Scale = 1.0;

        public override double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;


            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);





            // Calculate dissipative part
            // ==========================
            double[] VelocityMeanIn = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = U_Neg[1 + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);

            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanOut[d] = U_Pos[1 + d];
            }


            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);
            double Lambda = Math.Max(LambdaIn, LambdaOut);


            //double Flx = Lambda * (U_Pos[0] - Tsat) * capA * v_Neg - Lambda * (Tsat - U_Neg[0]) * capB * v_Pos;

            //double FlxNeg, FlxPos;
            //{
            //    double r = 0.0;


            //    // Calculate dissipative part
            //    // ==========================

            //    double[] VelocityMeanIn = new double[m_D];
            //    for (int d = 0; d < m_D; d++) {
            //        VelocityMeanIn[d] = ParamsNeg[m_D + d];
            //    }

            //    double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);

            //    double uJump = U_Neg[0] - Tsat;

            //    r += LambdaIn * uJump * LFFA;

            //    FlxNeg = capA * r;

            //}

            ////Flux for positive side
            //{
            //    double r = 0.0;

            //    // Calculate dissipative part
            //    // ==========================

            //    double[] VelocityMeanOut = new double[m_D];
            //    for (int d = 0; d < m_D; d++) {
            //        VelocityMeanOut[d] = ParamsPos[m_D + d];
            //    }


            //    double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            //    double uJump = Tsat - U_Pos[0];

            //    r += LambdaOut * uJump * LFFB;

            //    FlxPos = capB * r;

            //}

            //double Flx = FlxNeg * v_Neg - FlxPos * v_Pos;

            // contribution due to massflux
            //{
            //    double M = MassFlux(cp, Grad_uA, Grad_uB);
            //    Flx += Tsat * M * (1 / base.m_rhoA - 1 / base.m_rhoB) * 0.5 * (capA * v_Neg + capB * v_Pos);
            //}

            // ++++++++++
            // ++++++++++
            double UxN_Neg = 0, UxN_Pos = 0;
            for(int d = 0; d < m_D; d++) {
                UxN_Neg += U_Neg[1 + d] * cp.Normal[d];
                UxN_Pos += U_Pos[1 + d] * cp.Normal[d];
            }

            double Flx = Tsat * (UxN_Neg - UxN_Pos) * 0.5 * (capA * v_Neg + capB * v_Pos);
            Flx += (0.5 * (UxN_Neg * U_Neg[0] + UxN_Pos * U_Pos[0]) + Lambda * (U_Neg[0] - U_Pos[0])) * (capA * v_Neg - capB * v_Pos);
            //double Flx = 0.0;
            //Flx += 1.0 * (U_Neg[0] - Tsat) * capA * v_Neg;
            //Flx -= 1.0 * (Tsat - U_Pos[0]) * capB * v_Pos;
            //Flx += -(capA * v_Neg - capB * v_Pos) * Lambda * (U_Neg[0] - U_Pos[0]) + capA * v_Neg * Lambda * (U_Neg[0] - Tsat) - capB * v_Pos * Lambda * (Tsat - U_Pos[0]);
            //double Flx = Tsat * UxN_Neg * capA * v_Neg - Tsat * UxN_Pos * capB * v_Pos;
            // ++++++++++
            // ++++++++++

            
            double FlxNeg = Tsat * UxN_Neg;
            FlxNeg += LambdaIn * (U_Neg[0] - Tsat);
            double sFlxNeg = U_Neg[0] * UxN_Neg;
            FlxNeg += -sFlxNeg;

            double FlxPos = Tsat * UxN_Pos;
            FlxPos += LambdaOut * (Tsat - U_Pos[0]);
            double sFlxPos = U_Pos[0] * UxN_Pos;
            FlxPos += -sFlxPos;
            

            if (movingmesh)
                return 0.0;
            else
                return Scale * (capA * FlxNeg * v_Neg - capB * FlxPos * v_Pos);// Scale * Flx; //
        }


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        BitArray evapMicroRegion;


        public override IList<string> ArgumentOrdering {
            get {
                return base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D)) ;
            }
        }



        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    /// <summary>
    /// Extension for <see cref="HeatConvectionInSpeciesBulk_Hamiltonian_Newton"/> on interface
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian : MassFluxAtLevelSet_StrongCoupling, IEquationComponentCoefficient {

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian(int _D, double _capA, double _capB, double _LFFA, double _LFFB, bool _movingmesh, double _Tsat, ThermalParameters thermParams, string phaseA, string phaseB) : base(_D, thermParams, phaseA, phaseB) {

            m_D = _D;

            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            //DirichletCond = _DiriCond;
            Tsat = _Tsat;

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

        }

        //bool MaterialInterface;
        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        //bool DirichletCond;
        double Tsat;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }

        /// <summary>
        /// Scale of the component, used for homotopy
        /// </summary>
        double Scale = 1.0;

        public override double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);

            // Normal velocities
            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            double vINxN = 0.0, vOUTxN = 0.0;
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = U_Neg[1 + d];
                vINxN += VelocityMeanIn[d] * cp.Normal[d];
                VelocityMeanOut[d] = U_Pos[1 + d];
                vOUTxN += VelocityMeanOut[d] * cp.Normal[d];
            }

            //===========================================================================================================
            //===========================================================================================================
            // First variant, using central flux for temperature       
            /*
            double FlxNeg = Tsat * vINxN;
            double sFlxNeg = U_Neg[0] * vINxN;
            FlxNeg += -sFlxNeg;

            double FlxPos = Tsat * vOUTxN;
            double sFlxPos = U_Pos[0] * vOUTxN;
            FlxPos += -sFlxPos;           
            */
            //===========================================================================================================
            //===========================================================================================================


            //===========================================================================================================
            //===========================================================================================================
            // Second variant using Roe-Type Scheme
            // if VxN < 0 (Inflow) enforce Dirichlet condition, if > 0 (outflow), we still want to enforce saturation temperature?
            double FlxNeg = 0.5 * (vINxN - Math.Abs(vINxN)) * (Tsat - U_Neg[0]);
            FlxNeg -= 0.5 * (vINxN + Math.Abs(vINxN)) * (Tsat - U_Neg[0]);

            double FlxPos = 0.5 * (vOUTxN - Math.Abs(vOUTxN)) * (U_Pos[0] - Tsat);
            FlxPos -= 0.5 * (vOUTxN + Math.Abs(vOUTxN)) * (U_Pos[0] - Tsat);


            //===========================================================================================================
            //===========================================================================================================

            if (movingmesh)
                return 0.0;
            else
                return Scale * (capA * FlxNeg * v_Neg - capB * FlxPos * v_Pos);
        }


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        BitArray evapMicroRegion;


        public override IList<string> ArgumentOrdering {
            get {
                return base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));
            }
        }



        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }




    /// <summary>
    /// Extension for <see cref="HeatConvectionInSpeciesBulk_Hamiltonian_Newton"/> on interface
    /// Difference with <see cref="HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian"/> are the extra terms comming from the 
    /// weak formulation of the discretization
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian_LowMach : MassFluxAtLevelSet_StrongCoupling, IEquationComponentCoefficient {

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian_LowMach(int _D, double _capA, double _capB, double _LFFA, double _LFFB, bool _movingmesh, double _Tsat, ThermalParameters thermParams, string phaseA, string phaseB) : base(_D, thermParams, phaseA, phaseB) {

            m_D = _D;

            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            //DirichletCond = _DiriCond;
            Tsat = _Tsat;

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

        }

        //bool MaterialInterface;
        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        //bool DirichletCond;
        double Tsat;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }

        /// <summary>
        /// Scale of the component, used for homotopy
        /// </summary>
        double Scale = 1.0;

        public override double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);

            // Normal velocities
            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            double vINxN = 0.0, vOUTxN = 0.0;
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = U_Neg[1 + d];
                vINxN += VelocityMeanIn[d] * cp.Normal[d];
                VelocityMeanOut[d] = U_Pos[1 + d];
                vOUTxN += VelocityMeanOut[d] * cp.Normal[d];
            }

    

            //===========================================================================================================
            //===========================================================================================================
            // Second variant using Roe-Type Scheme
            // if VxN < 0 (Inflow) enforce Dirichlet condition, if > 0 (outflow), we still want to enforce saturation temperature?
            double FlxNeg = 0.5 * (vINxN - Math.Abs(vINxN)) * (Tsat - U_Neg[0]);
            FlxNeg -= 0.5 * (vINxN + Math.Abs(vINxN)) * (Tsat - U_Neg[0]);
            FlxNeg += vINxN * U_Neg[0]; // Extra term from weak formulation

            double FlxPos = 0.5 * (vOUTxN - Math.Abs(vOUTxN)) * (U_Pos[0] - Tsat);
            FlxPos -= 0.5 * (vOUTxN + Math.Abs(vOUTxN)) * (U_Pos[0] - Tsat);
            FlxPos -= vOUTxN * U_Pos[0];// Extra term from weak formulation         


            if (movingmesh)
                return 0.0;
            else
                return Scale * (capA * FlxNeg * v_Neg - capB * FlxPos * v_Pos);
        }


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        BitArray evapMicroRegion;


        public override IList<string> ArgumentOrdering {
            get {
                return base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));
            }
        }



        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }



    /// <summary>
    /// Extension for <see cref="HeatConvectionInSpeciesBulk_Hamiltonian_Newton"/> on interface
    /// </summary>
    public class SpeciesConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian : MassFluxAtLevelSet_StrongCoupling, IEquationComponentCoefficient {

        bool movingmesh;
        int component; //chemical species index
        public SpeciesConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian(int _D, double _rhoA, double _rhoB, double _LFFA, double _LFFB, bool _movingmesh, double _interfaceValue, ThermalParameters thermParams, string phaseA, string phaseB, int _component) : base(_D, thermParams, phaseA, phaseB) {

            m_D = _D;
            interfaceValue = _interfaceValue;
            capA = _rhoA;
            capB = _rhoB;
            LFFA = _LFFA;
            LFFB = _LFFB;
            component = _component;
        }

        int m_D;
        double capA;
        double capB;
        double LFFA;
        double LFFB;

        //bool DirichletCond;
        double interfaceValue;

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }

        /// <summary>
        /// Scale of the component, used for homotopy
        /// </summary>
        double Scale = 1.0;

        public override double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);

            // Normal velocities
            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            double vINxN = 0.0, vOUTxN = 0.0;
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = U_Neg[1 + d];
                vINxN += VelocityMeanIn[d] * cp.Normal[d];
                VelocityMeanOut[d] = U_Pos[1 + d];
                vOUTxN += VelocityMeanOut[d] * cp.Normal[d];
            }

            // Second variant using Roe-Type Scheme
            // if VxN < 0 (Inflow) enforce Dirichlet condition, if > 0 (outflow), we still want to enforce saturation temperature?
            double FlxNeg = 0.5 * (vINxN - Math.Abs(vINxN)) * (interfaceValue - U_Neg[3]);
            FlxNeg -= 0.5 * (vINxN + Math.Abs(vINxN)) * (interfaceValue - U_Neg[3]);
            FlxNeg += vINxN * U_Neg[3]; // Extra term from weak formulation

            double FlxPos = 0.5 * (vOUTxN - Math.Abs(vOUTxN)) * (U_Pos[3] - interfaceValue);
            FlxPos -= 0.5 * (vOUTxN + Math.Abs(vOUTxN)) * (U_Pos[3] - interfaceValue);
            FlxPos -= vOUTxN * U_Pos[3];// Extra term from weak formulation         
            return Scale * (capA * FlxNeg * v_Neg - capB * FlxPos * v_Pos);
        }


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);
            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        BitArray evapMicroRegion;


        public override IList<string> ArgumentOrdering {
            get {
                var bla = base.ArgumentOrdering.Cat(VariableNames.VelocityVector(m_D));
                return bla.Cat(VariableNames.MassFraction_n(component));
            }
        }



        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }



    /// <summary>
    /// 
    /// </summary>
    public class ConductivityAtLevelSet_material : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;
        string phaseA, phaseB;
        public ConductivityAtLevelSet_material(int SpatialDim, double _kA, double _kB, double _penalty, double _Tsat, string phaseA, string phaseB, int iLevSet = 0) {
            this.kA = _kA;
            this.kB = _kB;
            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;
            this.iLevSet = iLevSet;
            //m_LsTrk = lstrk;

        }

        double kA;
        double kB;

        double penalty;
        int m_D;
        int iLevSet;

        //bool DirichletCond;
        double Tsat;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(inp.jCellIn == inp.jCellOut);

            double[] N = inp.Normal;
            //double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];

            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == m_D);
            Debug.Assert(Grad_uB.GetLength(1) == m_D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = (Math.Abs(kA) > Math.Abs(kB)) ? kA : kB;

            double Ret = 0.0;

            Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - uB[0]);                     // symmetry term

            Ret += pnlty * wPenalty * (uA[0] - uB[0]) * (vA - vB); // penalty term


            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        BitArray evapMicroRegion;


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
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public int LevelSetIndex {
            get { return iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public string PositiveSpecies {
            get { return phaseB; }
        }

        public string NegativeSpecies {
            get { return phaseA; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }

    /// <summary>
    /// 
    /// </summary>
    public class ConductivityAtLevelSet_withMassflux : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;

        public ConductivityAtLevelSet_withMassflux(int SpatialDim, double _kA, double _kB, double _penalty, double _Tsat, bool ZeroGradientAtContactline) {
            this.kA = _kA;
            this.kB = _kB;
            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;

            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;
            this.ZeroGradientAtContactline = ZeroGradientAtContactline;
            //m_LsTrk = lstrk;

        }

        double kA;
        double kB;

        double penalty;
        int m_D;

        //bool DirichletCond;
        double Tsat;
        bool ZeroGradientAtContactline;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(inp.jCellIn == inp.jCellOut);

            double[] N = inp.Normal;
            //double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];

            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == m_D);
            Debug.Assert(Grad_uB.GetLength(1) == m_D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = (Math.Abs(kA) > Math.Abs(kB)) ? kA : kB;

            double Ret = 0.0;

            // old, extension to material form
            /*
            // symmetry term
            Ret += (Tsat - 0.5 * (uA[0] + uB[0])) * (kA * Grad_vA_xN - kB * Grad_vB_xN);

            // consistency term
            Ret -= 0.5 * (vA + vB) * (kA * Grad_uA_xN - kB * Grad_uB_xN);

            // penalty
            Ret += pnlty * wPenalty * vA * (uB[0] - Tsat) - pnlty * wPenalty * vB * (Tsat - uA[0]);   
            //Ret += pnlty * wPenalty * vA * (uA[0] - Tsat) - pnlty * wPenalty * vB * (Tsat - uB[0]);
            */

            /*
            // new standalone non-material form
            // symmetry term
            Ret -=  kA * Grad_vA_xN * (uA[0] - Tsat) - kB * Grad_vB_xN * (Tsat - uB[0]);

            // consistency term
            Ret -= (vA * kA * Grad_uA_xN - vB * kB * Grad_uB_xN);

            // penalty, like a dirichlet condition /for temperature from both sides
            Ret += pnlty * wPenalty * vA * (uA[0] - Tsat) - pnlty * wPenalty * vB * (Tsat - uB[0]);
            */


            double g_D = Tsat;

            if (!ZeroGradientAtContactline) {
                // dirichlet condition from A-side
                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Ret += (kA * Grad_uA[0, d]) * (vA) * nd;
                    Ret += (kA * Grad_vA[d]) * (uA[0] - g_D) * nd;
                }

                Ret -= wPenalty * (uA[0] - g_D) * (vA - 0) * pnlty;

                // dirichlet condition from B-side
                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Ret += (kB * Grad_uB[0, d]) * (-vB) * nd;
                    Ret += (kB * Grad_vB[d]) * (g_D - uB[0]) * nd;
                }

                Ret -= wPenalty * (g_D - uB[0]) * (0 - vB) * pnlty;
                Ret *= -1.0;
            } else {
                double a = 1.0 * inp.X[0].Pow2();
                // Robin condition from A-side
                Ret += a * (uA[0] - g_D) * (vA - 0);
                // Robin condition from B-side
                Ret += a * (g_D - uB[0]) * (0 - vB);
            }

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        BitArray evapMicroRegion;


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
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

    }


    /// <summary>
    /// 
    /// </summary>
    public class MassDifusivityAtLevelSet_withMassflux : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {


        public MassDifusivityAtLevelSet_withMassflux(int SpatialDim, double _kA, double _kB, double _penalty, double _Tsat, int _component, bool ZeroGradientAtContactline) {
            this.kA = _kA;
            this.kB = _kB;
            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;
            this.component = _component;
            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;
            this.ZeroGradientAtContactline = ZeroGradientAtContactline;
            //m_LsTrk = lstrk;

        }

        double kA;
        double kB;
        int component;
        double penalty;
        int m_D;

        //bool DirichletCond;
        double Tsat;
        bool ZeroGradientAtContactline;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(inp.jCellIn == inp.jCellOut);

            double[] N = inp.Normal;
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == m_D);
            Debug.Assert(Grad_uB.GetLength(1) == m_D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = (Math.Abs(kA) > Math.Abs(kB)) ? kA : kB;

            double Ret = 0.0;

            double g_D = Tsat;

            if (!ZeroGradientAtContactline) {
                // dirichlet condition from A-side
                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Ret += (kA * Grad_uA[0, d]) * (vA) * nd;
                    Ret += (kA * Grad_vA[d]) * (uA[0] - g_D) * nd;
                }

                Ret -= wPenalty * (uA[0] - g_D) * (vA - 0) * pnlty;

                // dirichlet condition from B-side
                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Ret += (kB * Grad_uB[0, d]) * (-vB) * nd;
                    Ret += (kB * Grad_vB[d]) * (g_D - uB[0]) * nd;
                }

                Ret -= wPenalty * (g_D - uB[0]) * (0 - vB) * pnlty;
                Ret *= -1.0;
            } else {
                double a = 1.0 * inp.X[0].Pow2();
                // Robin condition from A-side
                Ret += a * (uA[0] - g_D) * (vA - 0);
                // Robin condition from B-side
                Ret += a * (g_D - uB[0]) * (0 - vB);
            }

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        BitArray evapMicroRegion;


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
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.MassFraction_n(component)}; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

    }
}
