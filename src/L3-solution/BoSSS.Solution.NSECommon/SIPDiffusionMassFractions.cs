using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.NSECommon {
    public class SIPDiffusionMassFractions : SIPDiffusionBase, ISupportsJacobianComponent {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="PenaltyBase"></param>
        /// <param name="PenaltyLengthScales"></param>
        /// <param name="BcMap"></param>
        /// <param name="EoS"></param>
        /// <param name="Reynolds"></param>
        /// <param name="Prandtl"></param>
        /// <param name="prefactor"></param>
        /// <param name="prmsOK"></param>
        /// <param name="massFractionComponent"></param>
        /// <param name="numOfSpecies"></param>
        public SIPDiffusionMassFractions(double PenaltyBase, IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS, double Reynolds, double Prandtl, double[] Lewis, int massFractionComponent, int numOfSpecies) : base(PenaltyBase,  false, massFractionComponent + 1) {
            this.EoS = EoS;
            this.numOfSpecies = numOfSpecies;
            this.BcMap = BcMap;
            this.m_reynolds = Reynolds;
            this.m_prandtl = Prandtl;
            this.m_lewis = Lewis;
        }

        /// <summary>
        /// The equation of state used
        /// </summary>
        MaterialLaw EoS;

        /// <summary>
        /// 
        /// </summary>
        IncompressibleBoundaryCondMap BcMap;
        /// <summary>
        /// Number of total species 
        /// </summary>
        int numOfSpecies;

        /// <summary>
        /// Reynolds number
        /// </summary>
        double m_reynolds;

        /// <summary>
        /// Schmidt number
        /// </summary>
        double m_prandtl;

        /// <summary>
        /// Components Lewis numbers
        /// </summary>
        double[] m_lewis;

        ///// <summary>
        ///// Density multiplied by Diffusivity.
        ///// </summary>
        //double rhoD;


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;


            IncompressibleBcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall:
                //Neumann boundary condition
                Acc = 0.0;
                break;
                //case IncompressibleBcType.Wall:
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet:
                // inhom. Dirichlet b.c.
                // =====================
                double u_D = BcMap.bndFunction[VariableNames.MassFraction_n(i - 1)][inp.EdgeTag](inp.X, inp.time);
                Acc = -1 * base.BoundaryEdgeFormDirichlet(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, u_D);
                break;
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.NoSlipNeumann: {
                    Acc = 0.0;
                    break;
                }
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
            //double rhoDAd = ((MaterialLawLowMach)EoS).GetDiffusivity(U[0]);
            double rhoDAd = ((MaterialLaw_MultipleSpecies)EoS).GetDiffusivity(U[0]);

            int massfractionIndex = this.i - 1;
            double Lewis = m_lewis[massfractionIndex];
            double res = rhoDAd / (m_reynolds * m_prandtl * Lewis);
            return res;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivEdg, DerivVol };
        }

        public override IList<string> ArgumentOrdering {
            get {
                return (ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(numOfSpecies)));
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return new List<string> { }; // no parameters
            }
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
    public class SpeciesMassDiffusivity_AtLevelSet_material_LowMach : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;
        string phaseA, phaseB;
        public SpeciesMassDiffusivity_AtLevelSet_material_LowMach(int SpatialDim, MaterialLaw EoS, double Reynolds, double Prandtl, double[] Lewis, int massFractionComponent, int numOfSpecies, double _penalty, double _T_boundary, string phaseA, string phaseB, int iLevSet = 0) {
            this.EoS = EoS;


            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;

            this.phaseA = phaseA;
            this.phaseB = phaseB;
            this.m_reynolds = Reynolds;
            this.m_Prandtl = Prandtl;
            this.m_lewis = Lewis;
            this.numOfSpecies = numOfSpecies;

            //this.DirichletCond = _DiriCond;
            this.iLevSet = iLevSet;
            this.T_boundary = _T_boundary;

            this.massFractionComponent = massFractionComponent;
        }


        int m_D;
        int iLevSet;
        double T_boundary;

        /// <summary>
        /// 
        /// </summary>
        int massFractionComponent;

        /// <summary>
        /// Number of total species 
        /// </summary>
        int numOfSpecies;

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
        double[] m_lewis;

        /// <summary>
        /// For the
        /// </summary>
        protected double Diffusivity(double[] U, double[,] GradU, Vector NodeCoordinates) {
            double rhoDAd = ((MaterialLaw_MultipleSpecies)EoS).GetDiffusivity(U[0]);

            double Lewis = m_lewis[massFractionComponent];
            double res = rhoDAd / (m_reynolds * m_Prandtl * Lewis);
            return res;
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
                Grad_uA_xN += Grad_uA[1+ massFractionComponent, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];

            }
            double[] difusivityArguments_IN = uA;

            double DiffusivityA = Diffusivity(difusivityArguments_IN, Grad_uA, inp.X);

            double Ret = 0.0;
            double pnlty = 2 * this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = DiffusivityA;

            double Y_BOUND = massFractionComponent == 0 ? 1.0 : 0.0; //Hardcoded species at boundary


            Ret -= (DiffusivityA * Grad_uA_xN) * (vA);                           // consistency term
            Ret -= (DiffusivityA * Grad_vA_xN) * (uA[1+ massFractionComponent] - Y_BOUND);                     // symmetry term

            Ret += pnlty * wPenalty * (uA[1+ massFractionComponent] - Y_BOUND) * (vA); // penalty term


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

            double penaltySizeFactor = 1.0 / NegLengthScaleS[jCellIn];


            Debug.Assert(!double.IsNaN(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }


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
                return (ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(numOfSpecies)));
            }
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
    public class Interface_MassDiffusion_LowMach : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;
        private string phaseA, phaseB;

        public Interface_MassDiffusion_LowMach(int SpatialDim, MaterialLaw EoS_A, MaterialLaw EoS_B, double Reynolds, double Prandtl, double _penalty, double _Yinterface, string phaseA, string phaseB,  int Y_component, int iLevSet = 0) {

            this.m_penalty_base = _penalty;
            this.m_D = SpatialDim;
            this.phaseA = phaseA;
            this.phaseB = phaseB;
            this.m_reynolds = Reynolds;
            this.m_Prandtl = Prandtl;
            this.Y_interface = _Yinterface;
            this.iLevSet = iLevSet;
            this.EoS_A = EoS_A;
            this.EoS_B = EoS_B;

            varname = VariableNames.MassFraction_n(Y_component);
        }

        private int m_D;
        private int iLevSet;

        private double Y_interface;
 
        /// <summary>
        /// Reynolds number
        /// </summary>
        private double m_reynolds;

        /// <summary>
        /// Prandtl number
        /// </summary>
        private double m_Prandtl;

        /// <summary>
        /// The equation of state in phase A
        /// </summary>
        private MaterialLaw EoS_A;
        /// <summary>
        /// The equation of state in phase B
        /// </summary>
        private MaterialLaw EoS_B;

        /// <summary>
        /// For the
        /// </summary>
        protected double Diffusivity_A(double[] U, double[,] GradU, Vector NodeCoordinates) {
            double Diffusivity = ((MaterialLaw_MultipleSpecies)EoS_A).GetHeatConductivity(U[1]); // Just a Temperature dependence
            double Le = 1; //TODO

            Debug.Assert(Diffusivity >= 0.0);
            Diffusivity *= 1 / (m_reynolds * m_Prandtl * Le);
            return Diffusivity;
        }


        /// <summary>
        /// For the
        /// </summary>
        protected double Diffusivity_B(double[] U, double[,] GradU, Vector NodeCoordinates) {
            double Diffusivity = ((MaterialLaw_MultipleSpecies)EoS_B).GetHeatConductivity(U[0]); // Just a Temperature dependence

            Debug.Assert(Diffusivity >= 0.0);
            double Le = 1; //TODO
            Diffusivity *= 1 / (m_reynolds * m_Prandtl * Le);
            return Diffusivity;
        }


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(inp.jCellIn == inp.jCellOut);

            double[] N = inp.Normal; ;
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

            double[] HeatConducitityArguments_IN = uA;
            double[] HeatConducitityArguments_OT = uB;

            double DiffusivityA = Diffusivity_A(HeatConducitityArguments_IN, Grad_uA, inp.X);
            double DiffusivityB = Diffusivity_B(HeatConducitityArguments_OT, Grad_uB, inp.X);

            double Ret = 0.0;
            double pnlty =/* 2 * */this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;

            double g_D = this.Y_interface; ;//Hardcoded temperature at boundary
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
                return new List<string> { varname, VariableNames.Temperature };// transport parameter is temperature dependent
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
