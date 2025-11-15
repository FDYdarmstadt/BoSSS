using BoSSS.Foundation;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {
    class ExtensionSIP : SipViscosity_GradU {

        public ExtensionSIP(double _penalty_safety, int iComp, int D, IncompressibleBoundaryCondMap bcmap, ViscosityOption _ViscosityMode, double constantViscosityValue, StokesExtentionBoundaryOption useBCMap)
            : base(_penalty_safety, iComp, D, bcmap, _ViscosityMode, constantViscosityValue) {
            this.m_useBCMap = useBCMap;
        }

        StokesExtentionBoundaryOption m_useBCMap;

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            switch(m_useBCMap) {
                case StokesExtentionBoundaryOption.useBcMap: {
                    return base.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
                }
                case StokesExtentionBoundaryOption.FreeSlipAtWall: {

                    if(base.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.Wall) {  // redefines as slip wall
                                                                                       //return base.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
                        double Acc = 0.0;
                        double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);

                        double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
                        int D = inp.D;

                        for(int dN = 0; dN < D; dN++) {

                            for(int dD = 0; dD < D; dD++) {
                                // consistency
                                Acc += muA * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                                // symmetry
                                Acc += muA * (inp.Normal[m_iComp] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - 0.0) * inp.Normal[dN] * base.m_alpha;
                            }

                            // penalty
                            Acc -= muA * ((_uA[dN] - 0.0) * inp.Normal[dN]) * ((_vA - 0.0) * inp.Normal[m_iComp]) * pnlty;
                            //Acc = 0;
                        }
                        return Acc;
                    } else
                        return 0;
                }
                case StokesExtentionBoundaryOption.ZeroFlux:
                return 0;
                default:
                throw new NotImplementedException();
            }
        }
    }

    class ExtensionPressureGradient : PressureGradientLin_d {

        StokesExtentionBoundaryOption m_useBCMap;

        public ExtensionPressureGradient(int _d, IncompressibleBoundaryCondMap bcmap, StokesExtentionBoundaryOption useBCMap) : base(_d, bcmap) {
            this.m_useBCMap = useBCMap;
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            switch(m_useBCMap) {
                case StokesExtentionBoundaryOption.useBcMap: {
                    return base.BorderEdgeFlux(ref inp, Uin);
                }
                case StokesExtentionBoundaryOption.FreeSlipAtWall: {

                    if(base.m_bcmap.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.Wall) {
                        return base.BorderEdgeFlux(ref inp, Uin);
                    } else {
                        return 0;
                    }

                }
                case StokesExtentionBoundaryOption.ZeroFlux:
                return 0;
                default:
                throw new NotImplementedException();
            }
        }
    }

    class ExtensionDivergenceFlux : Divergence_DerivativeSource_Flux {

        StokesExtentionBoundaryOption m_useBCMap;

        public ExtensionDivergenceFlux(int component, IncompressibleBoundaryCondMap bcmap, StokesExtentionBoundaryOption useBCMap) : base(component, bcmap) {
            this.m_useBCMap = useBCMap;
        }

        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            switch(m_useBCMap) {
                case StokesExtentionBoundaryOption.useBcMap: {
                    base.BorderEdgeFlux_(ref inp, Uin, out FluxInCell);
                    return;
                }
                case StokesExtentionBoundaryOption.FreeSlipAtWall: {

                    if(base.bcmap.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.Wall) {
                        base.BorderEdgeFlux_(ref inp, Uin, out FluxInCell);
                    } else {
                        FluxInCell = 0;
                    }
                    return;
                }
                case StokesExtentionBoundaryOption.ZeroFlux: {
                    FluxInCell = 0;
                    return;
                }
                default:
                throw new NotImplementedException();
            }
        }
    }

}
