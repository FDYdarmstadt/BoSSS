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
        public HeatConvectionAtLevelSet_MovingMesh_withMassflux(int _D, LevelSetTracker LsTrk, double _Tsat, PhysicalParameters physicalParameters, ThermalParameters thermalParameters)
            : base(_D, LsTrk, physicalParameters) {

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
    public class HeatConvectionAtLevelSet_LLF_material : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_material(int _D, LevelSetTracker LsTrk, double _capA, double _capB, double _LFFA, double _LFFB,
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {

            m_D = _D;

            m_LsTrk = LsTrk;

            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, _LFFA, double.NaN, LsTrk);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, double.NaN, _LFFB, LsTrk);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));

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

            // central part

            FlxNeg = U_Neg[0] * (ParamsNeg[0] * cp.Normal[0] + ParamsNeg[1] * cp.Normal[1]);
            if (m_D == 3) {
                FlxNeg += U_Neg[0] * ParamsNeg[2] * cp.Normal[2];
            }

            FlxPos = U_Pos[0] * (ParamsPos[0] * cp.Normal[0] + ParamsPos[1] * cp.Normal[1]);
            if (m_D == 3) {
                FlxPos += U_Pos[0] * ParamsPos[2] * cp.Normal[2];
            }

            LLF = 0.5 * (FlxNeg + FlxPos);

            // dissipative part
            double[] VelocityMeanIn = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = ParamsNeg[m_D + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);

            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanOut[d] = ParamsPos[m_D + d];
            }

            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = U_Neg[0] - U_Pos[0];

            LLF += uJump * Lambda;

            if (movingmesh)
                return 0.0;
            else
                return LLF * ( capA * v_Neg - capB * v_Pos);
        }


        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

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

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
    }

    /// <summary>
    /// Extension for dissipative part of LLF to ensure $T=T_sat$ at the interface
    /// </summary>
    public class HeatConvectionAtLevelSet_LLF_withMassflux : MassFluxAtLevelSet {

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF_withMassflux(int _D, LevelSetTracker LsTrk, double _capA, double _capB, double _LFFA, double _LFFB,
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh, double _Tsat, PhysicalParameters _physicalParameters) : base(_D, LsTrk, _physicalParameters){

            m_D = _D;

            m_LsTrk = LsTrk;

            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, _LFFA, double.NaN, LsTrk);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, double.NaN, _LFFB, LsTrk);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));


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

            //Flux for negativ side
            double FlxNeg, FlxPos;
            {
                double r = 0.0;


                // Calculate dissipative part
                // ==========================

                double[] VelocityMeanIn = new double[m_D];
                for (int d = 0; d < m_D; d++) {
                    VelocityMeanIn[d] = ParamsNeg[m_D + d];
                }

                double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);

                double uJump = U_Neg[0] - Tsat;

                r += LambdaIn * uJump * LFFA;

                FlxNeg = capA * r;

            }

            //Flux for positive side
            {
                double r = 0.0;

                // Calculate dissipative part
                // ==========================

                double[] VelocityMeanOut = new double[m_D];
                for (int d = 0; d < m_D; d++) {
                    VelocityMeanOut[d] = ParamsPos[m_D + d];
                }


                double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

                double uJump = Tsat - U_Pos[0];

                r += LambdaOut * uJump * LFFB;

                FlxPos = capB * r;

            }

            double Flx = FlxNeg * v_Neg - FlxPos * v_Pos;

            // contribution due to massflux
            {
                double M = MassFlux(cp);
                Flx += Tsat * M * (1/base.m_rhoA - 1/base.m_rhoB) * 0.5 * (capA * v_Neg + capB * v_Pos);
            }

            if (movingmesh)
                return 0.0;
            else
                return Flx;
        }


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
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
    /// 
    /// </summary>
    public class ConductivityAtLevelSet_material : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ConductivityAtLevelSet_material(LevelSetTracker lstrk, double _kA, double _kB, double _penalty, double _Tsat) {
            this.kA = _kA;
            this.kB = _kB;
            this.m_penalty_base = _penalty;
            this.m_D = lstrk.GridDat.SpatialDimension;

            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;

            m_LsTrk = lstrk;

        }

        double kA;
        double kB;

        double penalty;
        int m_D;

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


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
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
    public class ConductivityAtLevelSet_withMassflux : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ConductivityAtLevelSet_withMassflux(LevelSetTracker lstrk, double _kA, double _kB, double _penalty, double _Tsat) {
            this.kA = _kA;
            this.kB = _kB;
            this.m_penalty_base = _penalty;
            this.m_D = lstrk.GridDat.SpatialDimension;

            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;

            m_LsTrk = lstrk;

        }

        double kA;
        double kB;

        double penalty;
        int m_D;

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

            // enforce saturation temperature
            Ret += (Tsat - 0.5 * (uA[0] + uB[0])) * (kA * Grad_vA_xN - kB * Grad_vB_xN);

            // incorporate massflux
            Ret -= 0.5 * (vA + vB) * (kA * Grad_uA_xN - kB * Grad_uB_xN);

            // extension of penalty
            Ret += pnlty * wPenalty * vA * (uB[0] - Tsat) - pnlty * wPenalty * vB * (Tsat - uA[0]);   

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


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
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

}
