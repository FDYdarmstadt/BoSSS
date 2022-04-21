using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
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
}
