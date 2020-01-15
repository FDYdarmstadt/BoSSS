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
using System.Diagnostics;
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

        /// <summary>
        /// Viscous term in k-omega model
        /// </summary>
        k_omegaModel
    }

    /// <summary>
    /// SIP discretization of diffusion operators for scalar transport equations (i.e. species mass transport and temperature). Analog to swipViscosity_Term1.
    /// </summary>
    public class SIPDiffusion : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IEquationComponentCoefficient, ISupportsJacobianComponent    {
 
        Func<double[],double> Coefficient;
        DiffusionMode Mode;

        string Argument;

        double PenaltyBase;
        string[] m_ParameterOrdering;
        IncompressibleBoundaryCondMap BcMap;
        Func<double[], double, double>[] ArgumentFunction;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Coefficient">Coefficient Function in \nabla \dot (Coefficient \nabla u)</param>
        /// <param name="PenaltyBase">C.f. Calculation of SIP penalty base, cf. Chapter 3 in 
        /// K. Hillewaert, “Development of the discontinuous Galerkin method for high-resolution, large scale CFD and acoustics in industrial geometries”,
        /// Université catholique de Louvain, 2013.</param>
        /// <param name="BcMap">Boundary condition map</param>
        /// <param name="Argument">The argument of the flux. Must be compatible with the DiffusionMode.</param>
        /// <param name="PenaltyLengthScales"></param>
        public SIPDiffusion(Func<double[],double> Coefficient, double PenaltyBase, MultidimensionalArray PenaltyLengthScales, IncompressibleBoundaryCondMap BcMap, string Argument) {
            this.Coefficient = Coefficient;
            this.PenaltyBase = PenaltyBase;
            this.BcMap = BcMap;
            this.ArgumentFunction = BcMap.bndFunction[Argument];
            this.Argument = Argument;
            this.cj = PenaltyLengthScales;
            switch (BcMap.PhysMode) {
                case PhysicsMode.LowMach:
                    this.m_ParameterOrdering = new string[] { VariableNames.Temperature0};
                    break;
                case PhysicsMode.Combustion:
                    this.m_ParameterOrdering = new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 };
                    break;
                case PhysicsMode.RANS:
                    this.m_ParameterOrdering = new string[] { "k0", "omega0" };
                    break;
                default:
                    throw new ApplicationException("Wrong physicsMode");
            }
        }

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="ConstantCoefficient">Constant Coefficient in \nabla \dot (Coefficient \nabla u)</param>
        /// <param name="PenaltyBase">C.f. Calculation of SIP penalty base, cf. Chapter 3 in 
        /// K. Hillewaert, “Development of the discontinuous Galerkin method for high-resolution, large scale CFD and acoustics in industrial geometries”,
        /// Université catholique de Louvain, 2013.</param>
        /// <param name="BcMap">Boundary condition map</param>
        /// <param name="Argument">The argument of the flux. Must be compatible with the DiffusionMode.</param>
        /// <param name="PenaltyLengthScales"></param>
        public SIPDiffusion(double ConstantCoefficient, double PenaltyBase, MultidimensionalArray PenaltyLengthScales, IncompressibleBoundaryCondMap BcMap, string Argument) {
            this.Coefficient = (x => ConstantCoefficient); //TODO does this work when inp.params is null?
            this.PenaltyBase = PenaltyBase;
            this.BcMap = BcMap;
            this.ArgumentFunction = BcMap.bndFunction[Argument];
            this.Argument = Argument;
            this.cj = PenaltyLengthScales;
            switch (BcMap.PhysMode) {
                case PhysicsMode.LowMach:
                    this.m_ParameterOrdering = new string[] { VariableNames.Temperature0};
                    break;
                case PhysicsMode.Combustion:
                    this.m_ParameterOrdering = new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 };
                    break;
                case PhysicsMode.RANS:
                    this.m_ParameterOrdering = new string[] { "k0", "omega0" };
                    break;
                default:
                    throw new ApplicationException("Wrong physicsMode");
            }
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
        //private double GetPenalty(int jCellIn, int jCellOut) {
        //    double cj_in = cj[jCellIn];
        //    double mu = PenaltyBase * cj_in;
        //    if (jCellOut >= 0) {
        //        double cj_out = cj[jCellOut];
        //        mu = Math.Max(mu, PenaltyBase * cj_out);
        //    }
        //    return mu;
        //}

        protected double GetPenalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor = 1;
            if (cj != null) {
                double penaltySizeFactor_A = 1.0 / cj[jCellIn];
                double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;
                Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
                Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
                Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
                Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
                penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);
            }

            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * PenaltyBase;
        }

        /// <summary>
        /// spatial dimension
        /// </summary>
        protected int m_D;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        /// <param name="cs"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_D = cs.GrdDat.SpatialDimension;
            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;
            double pnlty = GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);

            double DiffusivityA;
            double DiffusivityB;
            double DiffusivityMax;
            DiffusivityA = Coefficient(inp.Parameters_IN);
            DiffusivityB = Coefficient(inp.Parameters_OUT);

            foreach (var Diffusivity in new double[]{DiffusivityA, DiffusivityB})
            {
                Debug.Assert(!double.IsNaN(Diffusivity));
                Debug.Assert(!double.IsInfinity(Diffusivity));
            }

            for (int d = 0; d < inp.D; d++) {
                // consistency term
                Acc += 0.5 * (DiffusivityA * _Grad_uA[0, d] + DiffusivityB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];
                // symmetry term                
                Acc += 0.5 * (DiffusivityA * _Grad_vA[d] + DiffusivityB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];
            }
            // penalty term          
            DiffusivityMax = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * DiffusivityMax;
            return -Acc;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double DiffusivityA = Coefficient(inp.Parameters_IN);
            Debug.Assert(!double.IsNaN(DiffusivityA));
            Debug.Assert(!double.IsInfinity(DiffusivityA));
            IncompressibleBcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            // TODO
            switch (edgType) {
                case IncompressibleBcType.Wall: {
                    double u_D;
                    switch (Mode) {
                        case DiffusionMode.Temperature:
                            // inhom. Dirichlet b.c.
                            // =====================
                            u_D = ArgumentFunction[inp.EdgeTag](inp.X, 0);
                            for (int d = 0; d < inp.D; d++) {
                                Acc += (DiffusivityA * _Grad_uA[0, d]) * (_vA) * inp.Normal[d];
                                Acc += (DiffusivityA * _Grad_vA[d]) * (_uA[0] - u_D) * inp.Normal[d];
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
                                Acc += (DiffusivityA * _Grad_uA[0, d]) * (_vA) * inp.Normal[d];
                                Acc += (DiffusivityA * _Grad_vA[d]) * (_uA[0] - u_D) * inp.Normal[d];
                            }

                            Acc -= DiffusivityA * (_uA[0] - u_D) * (_vA - 0) * pnlty;
                            break;
                        case DiffusionMode.MassFraction:
                            double rhoA = 0.0;
                            rhoA = EoS.GetDensity(inp.Parameters_IN);
                            u_D = ArgumentFunction[inp.EdgeTag](inp.X, inp.time);
                            for (int d = 0; d < inp.D; d++) {
                                Acc += (DiffusivityA * rhoA * _Grad_uA[0, d]) * (_vA) * inp.Normal[d];
                                Acc += (DiffusivityA * rhoA * _Grad_vA[d]) * (_uA[0] - u_D) * inp.Normal[d];
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

            return -Acc;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0;
            double Diffusivity = Coefficient(cpv.Parameters);
            Debug.Assert(!double.IsNaN(Diffusivity));
            Debug.Assert(!double.IsInfinity(Diffusivity));

            for(int d = 0; d < cpv.D; d++)
                Acc -= Diffusivity * GradU[0, d] * GradV[d];
            return -Acc;
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
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
            get { return m_ParameterOrdering; }
        }

    }
}
