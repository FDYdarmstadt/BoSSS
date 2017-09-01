using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Heat conduction in temperature equation.
    /// Analog to swipViscosity_Term1.
    /// </summary>
    public class swipHeatConduction : BoSSS.Foundation.IEdgeComponent, BoSSS.Foundation.IVolumeForm {

        double Reynolds;
        double Prandtl;
        MaterialLaw EoS;

        double PenaltyBase;
        IncompressibleBoundaryCondMap BcMap;
        Func<double[],double>[] TemperatureFunction;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Reynolds"></param>
        /// <param name="Prandtl"></param>
        /// <param name="EoS"></param>
        /// <param name="PenaltyBase"></param>
        /// <param name="BcMap"></param>        
        public swipHeatConduction(double Reynolds, double Prandtl, MaterialLaw EoS, double PenaltyBase, IncompressibleBoundaryCondMap BcMap) {
            this.Reynolds = Reynolds;
            this.Prandtl = Prandtl;
            this.EoS = EoS;
            this.PenaltyBase = PenaltyBase;
            this.BcMap = BcMap;
            this.TemperatureFunction = BcMap.bndFunction[VariableNames.Temperature];
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        private double GetPenalty(int jCellIn, int jCellOut, double[] cj) {
            double cj_in = cj[jCellIn];
            double mu = PenaltyBase * cj_in;
            if (jCellOut >= 0) {
                double cj_out = cj[jCellOut];
                mu = Math.Max(mu, PenaltyBase * cj_out);
            }
            return mu;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = GetPenalty(inp.jCellIn, inp.jCellOut, inp.GridDat.Cells.cj);
            double lambdaA = EoS.GetViscosity(inp.Parameters_IN[0]);
            double lambdaB = EoS.GetViscosity(inp.Parameters_OUT[0]);

            for (int d = 0; d < inp.D; d++) {
                // consistency term
                Acc += 0.5 * (lambdaA * _Grad_uA[0, d] + lambdaB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];
                // symmetry term                
                Acc += 0.5 * (lambdaA * _Grad_vA[d] + lambdaB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];
            }

            // penalty term            
            double lambdaMax = (Math.Abs(lambdaA) > Math.Abs(lambdaB)) ? lambdaA : lambdaB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * lambdaMax;

            Acc *= 1.0 / (Reynolds * Prandtl);

            return -Acc;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * GetPenalty(inp.jCellIn, -1, inp.GridDat.Cells.cj);
            double lambdaA = EoS.GetViscosity(inp.Parameters_IN[0]);
            IncompressibleBcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall: {
                        // inhom. Dirichlet b.c.
                        // =====================

                        double T_D = TemperatureFunction[inp.EdgeTag](inp.X);
                        for (int d = 0; d < inp.D; d++) {
                            double nd = inp.Normale[d];
                            Acc += (lambdaA * _Grad_uA[0, d]) * (_vA) * nd;
                            Acc += (lambdaA * _Grad_vA[d]) * (_uA[0] - T_D) * nd;
                        }

                        Acc -= lambdaA * (_uA[0] - T_D) * (_vA - 0) * pnlty;

                        break;
                    }
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.NoSlipNeumann: {
                        Acc = 0.0;
                        break;
                    }
                default:
                    throw new NotSupportedException();
            }

            Acc *= 1.0 / (Reynolds * Prandtl);

            return -Acc;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0;
            double HeatConductivity = EoS.GetViscosity(cpv.Parameters[0]);
            for (int d = 0; d < cpv.D; d++)
                Acc -= HeatConductivity * GradU[0, d] * GradV[d];
            Acc *= 1.0 / (Reynolds * Prandtl);
            return -Acc;
        }

        /// <summary>
        /// Temperature
        /// </summary>
        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        /// <summary>
        /// Linearization point for temperature to calculate thermal conductivity.
        /// </summary>
        public IList<string> ParameterOrdering {
            get { return new string[] { VariableNames.Temperature0 }; }
        }
    }
}
