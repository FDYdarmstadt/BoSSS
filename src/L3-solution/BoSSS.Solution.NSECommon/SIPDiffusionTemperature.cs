using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution.NSECommon {

    public class SIPDiffusionTemperature : SIPDiffusionBase, ISupportsJacobianComponent {

        public enum ThermalWallType {

            /// <summary>
            /// Dirichlet value for T
            /// </summary>
            fixedTemperature,

            /// <summary>
            /// homogen Neumann value for T
            /// </summary>
            Adiabatic
        }

        private ThermalWallType myWallType;

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
        public SIPDiffusionTemperature(double PenaltyBase, IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS, double Reynolds, double Prandtl, bool prmsOK, ThermalWallType _myType) : base(PenaltyBase, prmsOK) {
            this.EoS = EoS;
            this.BcMap = BcMap;
            this.m_reynolds = Reynolds;
            this.m_Prandtl = Prandtl;
            varname = VariableNames.Temperature;
            myWallType = _myType;
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
                switch (myWallType) {
                    case ThermalWallType.Adiabatic:
                    Acc = 0.0;
                    break;

                    case ThermalWallType.fixedTemperature:
                    u_D = BcMap.bndFunction[varname][inp.EdgeTag](inp.X, inp.time);
                    Acc = -1 * BoundaryEdgeFormDirichlet2(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, u_D);
                    break;
                }

                //u_D = BcMap.bndFunction[varname][inp.EdgeTag](inp.X, inp.time);
                //Acc = -1 * base.BoundaryEdgeFormDirichlet(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, u_D);
                break;

                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Velocity_Inlet:
                    // inhom. Dirichlet b.c.
                    // =====================
                    u_D = BcMap.bndFunction[varname][inp.EdgeTag](inp.X, inp.time);
                    Acc = -1 * BoundaryEdgeFormDirichlet2(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, u_D);
                    break;
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet:
                    // inhom. Dirichlet b.c.
                    // =====================
                    u_D = BcMap.bndFunction[varname][inp.EdgeTag](inp.X, inp.time);
                    Acc = -1 * BoundaryEdgeFormDirichlet2(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, u_D);
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
        ///   The BoundaryEdgeForm for Dirichlet b.c. with value u_D
        /// </summary>
        protected double BoundaryEdgeFormDirichlet2(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA, double u_D) {
            double Acc = 0.0;

            double pnlty = 2 * GetPenalty(inp.jCellIn, -1);
            double[] difusivityArguments_IN = _uA;

            double DiffusivityA = Diffusivity(difusivityArguments_IN, _Grad_uA, inp.X);
            Debug.Assert(!double.IsNaN(DiffusivityA));
            Debug.Assert(!double.IsInfinity(DiffusivityA));
            double DiffusivityB = Diffusivity(new double[] { u_D}, _Grad_uA, inp.X);


            // penalty term          
            double  DiffusivityMax = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;


            // inhom. Dirichlet b.c.
            // =====================
            for (int d = 0; d < m_D; d++) {
                Acc += (DiffusivityA * _Grad_uA[i, d]) * (_vA) * inp.Normal[d];
                Acc += (DiffusivityA * _Grad_vA[d]) * (_uA[i] - u_D) * inp.Normal[d];
            }

            Acc -= DiffusivityMax * (_uA[i] - u_D) * (_vA - 0) * pnlty;
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
}