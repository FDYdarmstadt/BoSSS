using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS;

namespace NSE_SIMPLE {

    /// <summary>
    /// \f[ 
    ///   -\operatorname{div} \left( \mu \nabla \vec{u} \right)
    /// \f]
    /// </summary>
    public class swipViscosity_Term1_variante : swipViscosityBase {

        /// <summary>
        /// ctor; parameter documentation see <see cref="swipViscosityBase.swipViscosityBase"/>.
        /// </summary>
        public swipViscosity_Term1_variante(double _penalty, MultidimensionalArray PenaltyLengthScales, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
            ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null,
            Func<double, int, int, MultidimensionalArray, double> ComputePenalty = null)
            : base(_penalty, PenaltyLengthScales, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS, ComputePenalty) {

        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_iComp) };
            }
        }

        public override double VolumeForm(ref BoSSS.Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[0, d] * GradV[d] * Viscosity(cpv.Parameters) * base.m_alpha;
            //acc -= GradU[m_iComp, d] * GradV[d] * Viscosity(cpv.Parameters) * base.m_alpha;
            return -acc;
        }



        public override double InnerEdgeForm(ref BoSSS.Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            //switch (base.m_implMode) {
            //case ViscosityImplementation.H: //
            //{
            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (muA * _Grad_uA[0, d] + muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
                                                                                                            //Acc += 0.5 * (muA * _Grad_uA[m_iComp, d] + muB * _Grad_uB[m_iComp, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                                                                                                            //Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[m_iComp] - _uB[m_iComp]) * inp.Normale[d];  // symmetry term
            }
            Acc *= base.m_alpha;

            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term
                                                                    //Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax; // penalty term

            return -Acc;

            //}

            //    case ViscosityImplementation.SWIP: //
            //    {
            //            for (int d = 0; d < inp.D; d++) {
            //                Acc += (muB * muA * _Grad_uA[0, d] + muA * muB * _Grad_uB[0, d]) / (muA + muB) * (_vA - _vB) * inp.Normale[d];  // consistency term
            //                Acc += (muB * muA * _Grad_vA[d] + muA * muB * _Grad_vB[d]) / (muA + muB) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            //                //Acc += (muB * muA * _Grad_uA[m_iComp, d] + muA * muB * _Grad_uB[m_iComp, d]) / (muA + muB) * (_vA - _vB) * inp.Normale[d];  // consistency term
            //                //Acc += (muB * muA * _Grad_vA[d] + muA * muB * _Grad_vB[d]) / (muA + muB) * (_uA[m_iComp] - _uB[m_iComp]) * inp.Normale[d];  // symmetry term
            //            }
            //            Acc *= base.m_alpha;
            //            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * (2 * muA * muB) / (muA + muB); // penalty term
            //            //Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * (2 * muA * muB) / (muA + muB); // penalty term

            //        return -Acc;
            //    }
            //    default:
            //    throw new NotImplementedException();
            //}
        }


        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag) {
            if(base.g_Neu_Override == null) {
                return 0.0;
            } else {
                double Acc = 0;
                for(int i = 0; i < base.m_D; i++) {
                    Acc += N[i] * g_Neu_Override(base.m_iComp, X, i);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref BoSSS.Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++

                    double g_D = base.g_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    for(int d = 0; d < inp.D; d++) {
                        double nd = inp.Normale[d];
                        Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                        Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        //Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * nd;
                        //Acc += (muA * _Grad_vA[d]) * (_uA[m_iComp] - g_D) * nd;
                    }
                    Acc *= base.m_alpha;

                    Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                    //Acc -= muA * (_uA[m_iComp] - g_D) * (_vA - 0) * pnlty;
                    //            break;
                    //        }
                    //    default:
                    //        throw new NotImplementedException();

                    //}
                    break;
                }
                case IncompressibleBcType.FreeSlip: {
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: { 
                    // consistency
                    for(int d = 0; d < inp.D; d++) {
                        Acc += (inp.Normale[m_iComp] * muA * _Grad_uA[0, d] * inp.Normale[d]) * (_vA * inp.Normale[m_iComp]) * base.m_alpha;
                    }

                    // penalty
                    Acc -= muA * (_uA[0] - 0) * inp.Normale[m_iComp] * inp.Normale[m_iComp] * (_vA - 0) * pnlty;

                    break;
                    //    }
                    //    default:
                    //        throw new NotImplementedException();
                    //}
                    //break;
                }
                case IncompressibleBcType.NavierSlip_Linear: {

                    double[,] P = new double[inp.D, inp.D];
                    for(int d1 = 0; d1 < inp.D; d1++) {
                        for(int d2 = 0; d2 < inp.D; d2++) {
                            double nn = inp.Normale[d1] * inp.Normale[d2];
                            if(d1 == d2) {
                                P[d1, d2] = 1 - nn;
                            } else {
                                P[d1, d2] = -nn;
                            }
                        }
                    }

                    double g_D = base.g_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    for(int d = 0; d < inp.D; d++) {
                        // consistency
                        Acc += muA * (inp.Normale[m_iComp] * _Grad_uA[0, d] * inp.Normale[d]) * (_vA * inp.Normale[m_iComp]);
                        // symmetry
                        Acc += muA * (inp.Normale[m_iComp] * _Grad_vA[d] * inp.Normale[d]) * (_uA[0] - g_D) * inp.Normale[m_iComp];
                    }
                    Acc *= base.m_alpha;
                    // penalty
                    Acc -= muA * ((_uA[0] - g_D) * inp.Normale[m_iComp]) * ((_vA - 0) * inp.Normale[m_iComp]) * pnlty;

                    // tangential dissipation force term
                    Acc -= (m_beta * P[0, m_iComp] * (_uA[0] - g_D)) * (P[0, m_iComp] * _vA) * base.m_alpha;

                    //for (int d = 0; d < inp.D; d++) {
                    //    // no-penetration term
                    //    Acc += (inp.Normale[m_iComp] * muA * _Grad_uA[0, d] * inp.Normale[d]) * (_vA * inp.Normale[m_iComp]) * base.m_alpha;
                    //    // tangential dissipation force term
                    //    Acc -= (m_beta * P[0, d] * (_uA[0] - g_D)) * (P[0, d] * _vA) * base.m_alpha;
                    //}
                    //        break;
                    //    }
                    //default:
                    //    throw new NotImplementedException();
                    //}
                    break;
                }
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                    // Atmospheric outlet/pressure outflow: hom. Neumann
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    double g_N = g_Neu(inp.X, inp.Normale, inp.EdgeTag);
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    Acc += muA * g_N * _vA * base.m_alpha;
                    break;
                    //        }

                    //    default:
                    //        throw new NotImplementedException();

                    //}

                    //break;
                }
                case IncompressibleBcType.Pressure_Dirichlet: {
                    // Dirichlet boundary condition for pressure.
                    // Inner values of velocity gradient are taken, i.e.
                    // no boundary condition for the velocity (resp. velocity gradient) is imposed.                        
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    for(int d = 0; d < inp.D; d++) {
                        Acc += (muA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];
                        //Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * inp.Normale[d];
                    }
                    Acc *= base.m_alpha;
                    //          break;
                    //        }
                    //    default:
                    //        throw new NotImplementedException();
                    //}
                    break;
                }
                default:
                throw new NotImplementedException();
            }

            return -Acc;
        }
    }

}
