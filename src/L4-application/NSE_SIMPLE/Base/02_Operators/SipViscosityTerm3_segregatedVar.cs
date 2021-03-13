using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS;
using BoSSS.Foundation;

namespace NSE_SIMPLE {

     /// <summary>
    /// \f[ 
    ///   \frac{2}{3} \operatorname{div} \left( \mu \myMatrix{I} \operatorname{div} ( \vec{u} )  \right)
    /// \f]
    /// </summary>
    public class SipViscosity_Term3_segregatedVar : SipViscosityBase {

        //private ViscositySolverMode ViscSolverMode;

        /// <summary>
        /// ctor; parameter documentation see <see cref="SipViscosityBase.SipViscosityBase"/>.
        /// </summary>
        public SipViscosity_Term3_segregatedVar(double _penalty, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
                                   ViscosityOption _ViscosityMode/*, ViscositySolverMode ViscSolverMode = ViscositySolverMode.FullyCoupled*/,
                                   double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null, bool ignoreVectorized = false)
            : base(_penalty, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS, ignoreVectorized) {

            //this.ViscSolverMode = ViscSolverMode;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double visc = Viscosity(cpv.Parameters);
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[d, d] * GradV[base.m_iComp] * visc * base.m_alpha;
            return acc * (2.0 / 3.0);
        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            for(int i = 0; i < inp.D; i++) {
                // consistency term
                Acc += 0.5 * (muA * _Grad_uA[i, i] + muB * _Grad_uB[i, i]) * (_vA - _vB) * inp.Normal[m_iComp];
                // symmetry term
                //switch(ViscSolverMode) {
                //    case ViscositySolverMode.FullyCoupled:
                //        Acc += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normal[i];
                //        break;
                //    case ViscositySolverMode.Segregated:
                        if(i == m_iComp)
                            Acc += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normal[i];
                //        break;
                //    default:
                //        throw new NotImplementedException();
                //}
            }
            Acc *= base.m_alpha;

            // penalty term
            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;

            return Acc * (2.0 / 3.0);
        }



        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag) {
            if(base.g_Neu_Override == null) {
                //return 0.0;

                throw new NotSupportedException("Neumann BC. for the \\/U^T -- term is problematic!");

            } else {
                double Acc = 0;
                for(int i = 0; i < base.m_D; i++) {
                    Acc += N[m_iComp] * g_Neu_Override(i, X, i);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;
            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++                      

                    for(int i = 0; i < inp.D; i++) {
                        // consistency
                        Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normal[m_iComp];
                        // symmetry
                        //switch(ViscSolverMode) {
                        //    case ViscositySolverMode.FullyCoupled:
                        //        Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[i];
                        //        break;
                        //    case ViscositySolverMode.Segregated:
                                if(i == m_iComp)
                                    Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[i];
                        //        break;
                        //    default:
                        //        throw new NotImplementedException();
                        //}
                    }
                    Acc *= base.m_alpha;

                    // penalty
                    Acc -= muA * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;

                    break;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {

                    if(base.g_Neu_Override == null) {
                        // Inner values of velocity gradient are taken, i.e.
                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                        for(int i = 0; i < inp.D; i++) {
                            Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normal[m_iComp];
                        }
                    } else {
                        double g_N = g_Neu(inp.X, inp.Normal, inp.EdgeTag);
                        Acc += muA * g_N * _vA;
                    }
                    Acc *= base.m_alpha;

                    break;
                }
                default:
                    throw new NotSupportedException();
            }

            return Acc * (2.0 / 3.0);
        }
    }
}
