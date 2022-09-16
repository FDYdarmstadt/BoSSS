using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {
    class ExtensionSIP : SipViscosity_GradU {

        public ExtensionSIP(double _penalty_safety, int iComp, int D, IncompressibleBoundaryCondMap bcmap, ViscosityOption _ViscosityMode, double constantViscosityValue) 
            : base(_penalty_safety, iComp, D, bcmap,_ViscosityMode, constantViscosityValue) {

        }

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            //// freeslip, for now only non-moving walls
            switch (EdgeTag2Type[inp.EdgeTag]) {
                case IncompressibleBcType.NavierSlip_Linear:
                    double Acc = 0.0;
                    double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);

                    double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);

                    int D = inp.D;

                    for (int dN = 0; dN < D; dN++) {

                        for (int dD = 0; dD < D; dD++) {
                            // consistency
                            Acc += muA * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                            // symmetry
                            Acc += muA * (inp.Normal[m_iComp] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - 0.0) * inp.Normal[dN] * base.m_alpha;
                        }

                        // penalty
                        Acc -= muA * ((_uA[dN] - 0.0) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                    }
                    return -Acc;                    
                default:
                    return 0.0;            
            }

        }

    }

    class ExtensionPressureGradient : PressureGradientLin_d {

        int _d;
        IncompressibleBcType[] EdgeTag2Type;

        public ExtensionPressureGradient(int _d, IncompressibleBoundaryCondMap bcmap) : base(_d, bcmap) {
            EdgeTag2Type = bcmap.EdgeTag2Type;
            this._d = _d;
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            // freeslip, no condition for pressure
            //return Uin[0] * inp.Normal[m_d];
            switch (EdgeTag2Type[inp.EdgeTag]) {
                case IncompressibleBcType.NavierSlip_Linear:
                    return Uin[0] * inp.Normal[_d]; 
                default:
                    return 0.0; // pressure outlet
            }                   
        }
    }

    class ExtensionDivergenceFlux : Divergence_DerivativeSource_Flux {

        IncompressibleBcType[] EdgeTag2Type;

        public ExtensionDivergenceFlux(int component, IncompressibleBoundaryCondMap bcmap) : base(component, bcmap) {
            EdgeTag2Type = bcmap.EdgeTag2Type;
        }

        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            switch (EdgeTag2Type[inp.EdgeTag]) {
                case IncompressibleBcType.NavierSlip_Linear:
                    double u_j_In = Uin[0];
                    double u_j_Out = 0.0;// this.bndFunction[inp.EdgeTag](inp.X, inp.time);
                    FluxInCell = -(u_j_In - u_j_Out) * inp.Normal[component];
                    break;
                default:
                    FluxInCell = 0.0;
                    break;
            }
        }
    }
}
