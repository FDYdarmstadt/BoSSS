/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Mode for the diffusion operator
    /// i.e. operator part of energy balance (temperature equation) or species mass balance
    /// </summary>
    public enum DiffusionMode {

        /// <summary>
        /// Diffusion of heat
        /// </summary>
        Temperature,

        /// <summary>
        /// Diffusion of mass (molecules)
        /// </summary>
        MassFraction,
    }

    /// <summary>
    /// SIP discretization of diffusion operators for scalar transport equations (i.e. species mass transport and temperature). Analog to swipViscosity_Term1.
    /// </summary>
    public class SIPDiffusion : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {

        double Reynolds;
        double Schmidt;
        MaterialLaw EoS;
        DiffusionMode Mode;

        string Argument;

        double PenaltyBase;
        IncompressibleBoundaryCondMap BcMap;
        Func<double[], double, double>[] ArgumentFunction;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Reynolds"></param>
        /// <param name="Schmidt"></param>
        /// <param name="EoS">Material law</param>
        /// <param name="PenaltyBase">C.f. Calculation of SIP penalty base, cf. Chapter 3 in 
        /// K. Hillewaert, “Development of the discontinuous Galerkin method for high-resolution, large scale CFD and acoustics in industrial geometries”,
        /// Université catholique de Louvain, 2013.</param>
        /// <param name="BcMap">Boundary condition map</param>
        /// <param name="Mode">Equation type. Can be Temperature or MassFraction</param>
        /// <param name="Argument">The argument of the flux. Must be compatible with the DiffusionMode.</param>
        /// <param name="PenaltyLengthScales"></param>
        public SIPDiffusion(double Reynolds, double Schmidt, MaterialLaw EoS, double PenaltyBase, MultidimensionalArray PenaltyLengthScales, IncompressibleBoundaryCondMap BcMap, DiffusionMode Mode, string Argument) {
            this.Reynolds = Reynolds;
            this.Schmidt = Schmidt;
            this.EoS = EoS;
            this.PenaltyBase = PenaltyBase;
            this.BcMap = BcMap;
            this.ArgumentFunction = BcMap.bndFunction[Argument];
            this.Mode = Mode;
            this.Argument = Argument;
            this.cj = PenaltyLengthScales;
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

        MultidimensionalArray cj;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        private double GetPenalty(int jCellIn, int jCellOut) {
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
            double pnlty = GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);

            double DiffusivityA;
            double DiffusivityB;
            double DiffusivityMax;
            double rhoMax;
            switch (Mode) {
                case DiffusionMode.Temperature:
                    DiffusivityA = ((MaterialLawCombustion)EoS).GetHeatConductivity(inp.Parameters_IN[0]);
                    DiffusivityB = ((MaterialLawCombustion)EoS).GetHeatConductivity(inp.Parameters_OUT[0]);

                    for (int d = 0; d < inp.D; d++) {
                        // consistency term
                        Acc += 0.5 * (DiffusivityA * _Grad_uA[0, d] + DiffusivityB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];
                        // symmetry term                
                        Acc += 0.5 * (DiffusivityA * _Grad_vA[d] + DiffusivityB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];
                    }
                    // penalty term          
                    DiffusivityMax = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;
                    Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * DiffusivityMax;
                    break;
                case DiffusionMode.MassFraction:
                    double rhoA = 0.0;
                    double rhoB = 0.0;
                    rhoA = EoS.GetDensity(inp.Parameters_IN);
                    rhoB = EoS.GetDensity(inp.Parameters_OUT);
                    DiffusivityA = ((MaterialLawCombustion)EoS).GetDiffusivity(inp.Parameters_IN[0]);
                    DiffusivityB = ((MaterialLawCombustion)EoS).GetDiffusivity(inp.Parameters_OUT[0]);

                    for (int d = 0; d < inp.D; d++) {
                        // consistency term
                        Acc += 0.5 * (DiffusivityA * rhoA * _Grad_uA[0, d] + rhoB * DiffusivityB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];
                        // symmetry term                
                        Acc += 0.5 * (DiffusivityA * rhoA * _Grad_vA[d] + DiffusivityB * rhoB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];
                    }
                    // penalty term       
                    DiffusivityMax = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;
                    rhoMax = (rhoA > rhoB) ? rhoA : rhoB;
                    Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * rhoMax * DiffusivityMax;
                    break;
                default:
                    throw new NotImplementedException();
            }

            Acc *= 1.0 / (Reynolds * Schmidt);
            return -Acc;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double DiffusivityA;
            switch (Mode) {
                case DiffusionMode.Temperature:
                    DiffusivityA = ((MaterialLawCombustion)EoS).GetHeatConductivity(inp.Parameters_IN[0]);
                    break;
                case DiffusionMode.MassFraction:
                    DiffusivityA = ((MaterialLawCombustion)EoS).GetDiffusivity(inp.Parameters_IN[0]);
                    break;
                default:
                    throw new NotImplementedException();
            }
            IncompressibleBcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall: {
                        double u_D;
                        switch (Mode) {
                            case DiffusionMode.Temperature:
                                // inhom. Dirichlet b.c.
                                // =====================
                                u_D = ArgumentFunction[inp.EdgeTag](inp.X, 0);
                                for (int d = 0; d < inp.D; d++) {
                                    Acc += (DiffusivityA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];
                                    Acc += (DiffusivityA * _Grad_vA[d]) * (_uA[0] - u_D) * inp.Normale[d];
                                }

                                Acc -= DiffusivityA * (_uA[0] - u_D) * (_vA - 0) * pnlty;
                                break;
                            case DiffusionMode.MassFraction:
                                //Neumann boundary condition
                                Acc = 0.0;
                                break;
                            default:
                                throw new NotImplementedException();
                        }
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        // inhom. Dirichlet b.c.
                        // =====================

                        double u_D;
                        switch (Mode) {
                            case DiffusionMode.Temperature:
                                u_D = ArgumentFunction[inp.EdgeTag](inp.X, inp.time);
                                for (int d = 0; d < inp.D; d++) {
                                    Acc += (DiffusivityA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];
                                    Acc += (DiffusivityA * _Grad_vA[d]) * (_uA[0] - u_D) * inp.Normale[d];
                                }

                                Acc -= DiffusivityA * (_uA[0] - u_D) * (_vA - 0) * pnlty;
                                break;
                            case DiffusionMode.MassFraction:
                                double rhoA = 0.0;
                                rhoA = EoS.GetDensity(inp.Parameters_IN);
                                u_D = ArgumentFunction[inp.EdgeTag](inp.X, inp.time);
                                for (int d = 0; d < inp.D; d++) {
                                    Acc += (DiffusivityA * rhoA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];
                                    Acc += (DiffusivityA * rhoA * _Grad_vA[d]) * (_uA[0] - u_D) * inp.Normale[d];
                                }

                                Acc -= DiffusivityA * rhoA * (_uA[0] - u_D) * (_vA - 0) * pnlty;
                                break;
                            default:
                                throw new NotImplementedException();
                        }
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

            Acc *= 1.0 / (Reynolds * Schmidt);
            return -Acc;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0;
            double Diffusivity;
            double rho = 0.0;
            rho = EoS.GetDensity(cpv.Parameters);
            switch (Mode) {
                case DiffusionMode.Temperature:
                    Diffusivity = ((MaterialLawCombustion)EoS).GetHeatConductivity(cpv.Parameters[0]);
                    for (int d = 0; d < cpv.D; d++)
                        Acc -= Diffusivity * GradU[0, d] * GradV[d];
                    break;
                case DiffusionMode.MassFraction:
                    Diffusivity = ((MaterialLawCombustion)EoS).GetDiffusivity(cpv.Parameters[0]);
                    for (int d = 0; d < cpv.D; d++)
                        Acc -= Diffusivity * rho * GradU[0, d] * GradV[d];
                    break;
                default:
                    throw new NotImplementedException();
            }
            Acc *= 1.0 / (Reynolds * Schmidt);
            return -Acc;
        }

        /// <summary>
        /// Arguments
        /// </summary>
        public IList<string> ArgumentOrdering {
            get { return new string[] { Argument }; }
        }

        /// <summary>
        /// Parameters at linearization point to calculate material properties.
        /// </summary>
        public IList<string> ParameterOrdering {
            get { return new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 }; }
        }
    }
}
