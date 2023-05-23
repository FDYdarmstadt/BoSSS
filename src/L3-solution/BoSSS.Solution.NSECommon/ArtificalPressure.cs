using BoSSS.Foundation;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// SIP discretization of Laplace operator for pressure correction,
    /// with penalty for pressure correction.
    /// </summary>
    public class ArtificalPressure : SIPLaplace {
        IncompressibleBoundaryCondMap m_BcMap;

        /// <summary>
        /// Ctor.
        /// </summary>
        public ArtificalPressure(double penalty_base, IncompressibleBoundaryCondMap _BcMap, double artviscosity)
            : base(penalty_base, VariableNames.Pressure) {
            m_artviscosity = artviscosity;
            m_BcMap = _BcMap;
        }

        private double m_artviscosity = 1.0;

        public override double Nu(double[] x, double[] p, int jCell) {
            double h = base.LengthScales[jCell];
            return m_artviscosity * h * h;
        }

        /// <summary>
        /// this is the place to set pressure Dirichlet
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        protected override double g_Diri(ref Foundation.CommonParamsBnd inp) { return 0; }

        /// <summary>
        /// 
        /// </summary>
        protected override bool IsDirichlet(ref BoSSS.Foundation.CommonParamsBnd inp) {
            IncompressibleBcType edgType = m_BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Outflow:
                return true;
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                return false;
                default:
                throw new NotImplementedException("unsupported/unknown b.c. - missing implementation;");
            }
        }

        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);
                //double g_D = _uA[0];

                for (int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                    Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                }
                Acc *= this.m_alpha;

                Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty

            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }
}