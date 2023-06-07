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
    ///   - \operatorname{div} \left( \mu (\partial_d \vec{u})^T \right)
    /// \f]
    /// </summary>
    public class SipViscosity_Term2_segregatedVar : SipViscosityBase {

        //private ViscositySolverMode ViscSolverMode;

        /// <summary>
        /// ctor; parameter documentation see <see cref="SipViscosityBase.SipViscosityBase"/>.
        /// </summary>
        public SipViscosity_Term2_segregatedVar(double _penalty, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
                                   ViscosityOption _ViscosityMode, /*ViscositySolverMode ViscSolverMode = ViscositySolverMode.FullyCoupled,*/
                                   double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null)
            : base(_penalty, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS) {

            //this.ViscSolverMode = ViscSolverMode;
        }

        public override double VolumeForm(ref BoSSS.Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            double visc = Viscosity(cpv.Parameters, U, GradU);
            for(int d = 0; d < cpv.D; d++)
                // we want to:
                //    sum(  \partial_{m_iComp} u_d  * \partial_{d} v, d=0..D-1)
                acc += GradU[d, base.m_iComp] * GradV[d] * visc * base.m_alpha;
            return acc;
        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            double muB = this.Viscosity(inp.Parameters_OUT, _uB, _Grad_uB);


            for(int i = 0; i < inp.D; i++) {
                // consistency term
                Acc += 0.5 * (muA * _Grad_uA[i, m_iComp] + muB * _Grad_uB[i, m_iComp]) * (_vA - _vB) * inp.Normal[i];
                // symmetry term
                //switch(ViscSolverMode) {
                //    case ViscositySolverMode.FullyCoupled:
                //        Acc += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normal[m_iComp];
                //        break;
                //    case ViscositySolverMode.Segregated:
                        if(i == m_iComp)
                            Acc += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normal[m_iComp];
                //        break;
                //    default:
                //        throw new NotImplementedException();
                //}
            }
            Acc *= base.m_alpha;

            // penalty term
            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;

            return -Acc;
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
                    Acc += N[i] * g_Neu_Override(i, X, base.m_iComp);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++
                     

                    for(int i = 0; i < inp.D; i++) {
                        // consistency
                        Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                        // symmetry
                        //switch(ViscSolverMode) {
                        //    case ViscositySolverMode.FullyCoupled:
                        //        Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[m_iComp];
                        //        break;
                        //    case ViscositySolverMode.Segregated:
                                if(i == m_iComp)
                                    Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[m_iComp];
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
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry: {
                    throw new NotImplementedException();
                    /*
                    int D = inp.D;
                    double g_D;

                    for(int dN = 0; dN < D; dN++) {
                        for(int dD = 0; dD < D; dD++) {
                            
                            // consistency
                            Acc += muA * (inp.Normal[dN] * _Grad_uA[dD, dN] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                            // symmetry
                            //switch(ViscSolverMode) {
                            //    case ViscositySolverMode.FullyCoupled:
                            //        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dD);
                            //        Acc += muA * (inp.Normal[dN] * _Grad_vA[dN] * inp.Normal[m_iComp]) * (_uA[dD] - g_D) * inp.Normal[dD] * base.m_alpha;
                            //        break;
                            //    case ViscositySolverMode.Segregated:
                            //    default:
                            //}
                        }
                        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);
                        // penalty
                        Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                    }

                    break;*/
                }
                case IncompressibleBcType.NavierSlip_Linear: {

                    double ls = Lslip[inp.jCellIn];
                    if(ls == 0.0)
                        goto case IncompressibleBcType.Velocity_Inlet;
                    else
                        goto case IncompressibleBcType.FreeSlip;


                }
                //case IncompressibleBcType.InnerValues: {

                //        for (int i = 0; i < inp.D; i++) {
                //            Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                //        }
                //        Acc *= base.m_alpha;

                //        break;
                //    }
                case IncompressibleBcType.Pressure_Dirichlet: 
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {

                    if (base.g_Neu_Override == null) {
                        // Inner values of velocity gradient are taken, i.e.
                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                        for(int i = 0; i < inp.D; i++) {
                            Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
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

            return -Acc;
        }

    }


}
