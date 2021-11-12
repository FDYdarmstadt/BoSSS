using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
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
}