using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

//namespace BoSSS.Solution.NSECommon {
namespace BoSSS.Application.XNSEC {

    public class SIPDiffusionMixtureFraction : SIPDiffusionBase, ISupportsJacobianComponent {

        /// <summary>
        ///
        /// </summary>
        /// <param name="PenaltyBase"></param>
        /// <param name="PenaltyLengthScales"></param>
        /// <param name="BcMap"></param>
        /// <param name="EoS"></param>
        /// <param name="Reynolds"></param>
        /// <param name="Prandtl"></param>
        /// <param name="prmsOK"></param>
        public SIPDiffusionMixtureFraction(double PenaltyBase, IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS, double Reynolds, double Prandtl, bool prmsOK) : base(PenaltyBase, prmsOK) {
            this.EoS = EoS;
            this.BcMap = BcMap;
            this.m_reynolds = Reynolds;
            this.m_Prandtl = Prandtl;
            varname = VariableNames.MixtureFraction;
        }

        /// <summary>
        /// The equation of state used
        /// </summary>
        private MaterialLaw EoS;

        /// <summary>
        ///
        /// </summary>
        private IncompressibleBoundaryCondMap BcMap;

        /// <summary>
        /// Name of the variable
        /// </summary>
        private string varname;

        public override IList<string> ArgumentOrdering {
            get {
                return new List<string> { varname };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return new List<string> { }; // no parameters
            }
        }

        /// <summary>
        /// Reynolds number
        /// </summary>
        private double m_reynolds;

        /// <summary>
        /// Prandtl number
        /// </summary>
        private double m_Prandtl;

        /// <summary>
        ///
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="_uA"></param>
        /// <param name="_Grad_uA"></param>
        /// <param name="_vA"></param>
        /// <param name="_Grad_vA"></param>
        /// <returns></returns>

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;
            double u_D;

            IncompressibleBcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall:
                    Acc = 0; // Zero-Neumann
                    break;

                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet:
                    // inhom. Dirichlet b.c.
                    // =====================
                    u_D = BcMap.bndFunction[varname][inp.EdgeTag](inp.X, inp.time);
                    Acc = -1 * base.BoundaryEdgeFormDirichlet(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, u_D);
                    break;

                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.NoSlipNeumann:
                    Acc = base.BoundaryEdgeFormNeumann();
                    break;

                default:
                    throw new NotSupportedException();
            }

            return -Acc;
        }

        /// <summary>
        /// For the
        /// </summary>
        /// <param name="Parameters"></param>
        /// <returns></returns>
        protected override double Diffusivity(double[] U, double[,] GradU, Vector NodeCoordinates) {
            //double Diffusivity = ((MaterialLawLowMach)EoS).GetHeatConductivity(U[0]); // Just a Temperature dependence
            double Diffusivity = ((MaterialLaw_MultipleSpecies)EoS).GetHeatConductivity(U[0]); // Just a Temperature dependence

            Debug.Assert(Diffusivity >= 0.0);
            Diffusivity *= 1 / (m_reynolds * m_Prandtl);
            return Diffusivity;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivEdg, DerivVol };
        }

        public override void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);
            // Set the Reynolds number to a user defined value contained in the CoefficientSet cs
            // Useful in case that the Reynolds number changes during a simulation...
            if (cs.UserDefinedValues.Keys.Contains("Reynolds"))
                m_reynolds = (double)cs.UserDefinedValues["Reynolds"];
        }
    }

    /// <summary>
    ///
    /// </summary>
    public class MixtureFractionDiffusivity_AtLevelSet_material_LowMach : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;
        private string phaseA, phaseB;

        public MixtureFractionDiffusivity_AtLevelSet_material_LowMach(int SpatialDim, MaterialLaw EoS, double Reynolds, double Prandtl, double _penalty, string phaseA, string phaseB, int iLevSet = 0) {
            this.EoS = EoS;

            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;

            this.phaseA = phaseA;
            this.phaseB = phaseB;
            this.m_reynolds = Reynolds;
            this.m_Prandtl = Prandtl;

            this.iLevSet = iLevSet;
            varname = VariableNames.MixtureFraction;
        }

        private int m_D;
        private int iLevSet;
        private double T_boundary;

        /// <summary>
        ///
        /// </summary>
        private int massFractionComponent;

        /// <summary>
        /// Number of total species
        /// </summary>
        private int numOfSpecies;

        /// <summary>
        /// The equation of state used
        /// </summary>
        private MaterialLaw EoS;

        /// <summary>
        /// Reynolds number
        /// </summary>
        private double m_reynolds;

        /// <summary>
        /// Prandtl number
        /// </summary>
        private double m_Prandtl;

        /// <summary>
        /// Components Lewis numbers
        /// </summary>
        private double[] m_lewis;

        /// <summary>
        /// For the
        /// </summary>
        protected double Diffusivity(double[] U, double[,] GradU, Vector NodeCoordinates) {
            double Diffusivity = ((MaterialLaw_MultipleSpecies)EoS).GetHeatConductivity(U[0]); // Just a Temperature dependence

            Debug.Assert(Diffusivity >= 0.0);
            Diffusivity *= 1 / (m_reynolds * m_Prandtl);
            return Diffusivity;
        }

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

            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
            }
            double[] difusivityArguments_IN = uA;

            double DiffusivityA = Diffusivity(difusivityArguments_IN, Grad_uA, inp.X);

            double Ret = 0.0;
            double pnlty = 2 * this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = DiffusivityA;

            double Z_BOUND = 1.0; //Hardcoded mixture fraction at boundary

            Ret -= (DiffusivityA * Grad_uA_xN) * (vA);                           // consistency term
            Ret -= (DiffusivityA * Grad_vA_xN) * (uA[0] - Z_BOUND);                     // symmetry term

            Ret += pnlty * wPenalty * (uA[0] - Z_BOUND) * (vA); // penalty term

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        private BitArray evapMicroRegion;

        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        private double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double Penalty(int jCellIn, int jCellOut) {
            double penaltySizeFactor = 1.0 / NegLengthScaleS[jCellIn];

            Debug.Assert(!double.IsNaN(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }

        private MultidimensionalArray NegLengthScaleS;

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
            get {
                return new List<string> { varname };
            }
        }

        /// <summary>
        /// Name of the variable
        /// </summary>
        private string varname;

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
    public class Interface_MixtureFractionDiffusivity_LowMach : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;
        private string phaseA, phaseB;

        public Interface_MixtureFractionDiffusivity_LowMach(int SpatialDim, MaterialLaw EoS, double Reynolds, double Prandtl, double _penalty, string phaseA, string phaseB, int iLevSet = 0) {
            this.EoS = EoS;

            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;

            this.phaseA = phaseA;
            this.phaseB = phaseB;
            this.m_reynolds = Reynolds;
            this.m_Prandtl = Prandtl;

            this.iLevSet = iLevSet;
            varname = VariableNames.MixtureFraction;
        }

        private int m_D;
        private int iLevSet;

        /// <summary>
        ///
        /// </summary>
        private int massFractionComponent;

        /// <summary>
        /// Number of total species
        /// </summary>
        private int numOfSpecies;

        /// <summary>
        /// The equation of state used
        /// </summary>
        private MaterialLaw EoS;

        /// <summary>
        /// Reynolds number
        /// </summary>
        private double m_reynolds;

        /// <summary>
        /// Prandtl number
        /// </summary>
        private double m_Prandtl;

        /// <summary>
        /// Components Lewis numbers
        /// </summary>
        private double[] m_lewis;

        /// <summary>
        /// For the
        /// </summary>
        protected double Diffusivity(double[] U, double[,] GradU, Vector NodeCoordinates) {
            double Diffusivity = ((MaterialLaw_MultipleSpecies)EoS).GetHeatConductivity(U[0]); // Just a Temperature dependence

            Debug.Assert(Diffusivity >= 0.0);
            Diffusivity *= 1 / (m_reynolds * m_Prandtl);
            return Diffusivity;
        }

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

            double[] difusivityArguments_IN = uA;
            double[] difusivityArguments_OT = uB;

            double DiffusivityA = Diffusivity(difusivityArguments_IN, Grad_uA, inp.X);
            double DiffusivityB = Diffusivity(difusivityArguments_OT, Grad_uB, inp.X);

            double Ret = 0.0;
            double pnlty =/* 2 * */this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;

            double g_D = 1.0;//Hardcoded mixture fraction at boundary
            bool ZeroGradientAtContactline = false;
            if (!ZeroGradientAtContactline) {
                // dirichlet condition from A-side
                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Ret += (DiffusivityA * Grad_uA[0, d]) * (vA) * nd;
                    Ret += (DiffusivityA * Grad_vA[d]) * (uA[0] - g_D) * nd;
                }

                Ret -= wPenalty * (uA[0] - g_D) * (vA - 0) * pnlty;

                // dirichlet condition from B-side
                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Ret += (DiffusivityB * Grad_uB[0, d]) * (-vB) * nd;
                    Ret += (DiffusivityB * Grad_vB[d]) * (g_D - uB[0]) * nd;
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

        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        private double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double Penalty(int jCellIn, int jCellOut) {
            double penaltySizeFactor = 1.0 / NegLengthScaleS[jCellIn];

            Debug.Assert(!double.IsNaN(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }

        private MultidimensionalArray NegLengthScaleS;

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

            //if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
            //    evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public int LevelSetIndex {
            get { return iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new List<string> { varname };
            }
        }

        /// <summary>
        /// Name of the variable
        /// </summary>
        private string varname;

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



}