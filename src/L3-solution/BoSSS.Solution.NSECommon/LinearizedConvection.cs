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
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Linearized convective operator (used e.g. in the SIMPLE-algorithm) according to:
    /// K. Shahbazi, P. F. Fischer, C. R. Ethier,
    /// A high-order discontinuous Galerkin method for the unsteady incompressible Navier-Stokes equations,
    /// J. Comput. Phys. 222 (2007) 391–407.
    /// (see BoSSS\doc\notes\0010-BoSSS-References\ShahbaziEtAl2007)
    /// </summary>
    public class LinearizedConvection : LinearFlux {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;
        //BoundaryCondMap<IncompressibleBcType> m_bcmap;
        IncompressibleBoundaryCondMap m_bcmap;

        /// <summary>
        /// Component index of the momentum equation.
        /// </summary>
        protected int m_component;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        PhysicsMode PhysMode;

        /// <summary>
        /// Ctor for common part of incompressible and low Mach number flows.
        /// </summary>
        /// <param name="SpatDim"></param>
        /// <param name="_bcmap"></param>
        /// <param name="_component"></param>
        private LinearizedConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);

            PhysMode = _bcmap.PhysMode;
        }

        bool m_UseBoundaryVelocityParameter = false;

        /// <summary>
        /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
        /// </summary>
        protected double LaxFriedrichsSchemeSwitch = 1.0;

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>
        /// <param name="SpatDim">
        /// Spatial dimension (either 2 or 3).
        /// </param>
        /// <param name="_bcmap"></param>
        /// <param name="_component"></param>
        /// <param name="UseBoundaryVelocityParameter">
        /// True, if (an offset to) the boundary velocity is supplied in form of parameter variables.
        /// </param>
        public LinearizedConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, bool UseBoundaryVelocityParameter = false)
            : this(SpatDim, _bcmap, _component) {

            m_UseBoundaryVelocityParameter = UseBoundaryVelocityParameter;

            m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim));
            if (m_UseBoundaryVelocityParameter)
                m_ParameterOrdering = ArrayTools.Cat(m_ParameterOrdering, VariableNames.BoundaryVelocityVector(SpatDim));
        }

        //bool m_VariableDensity = false;
        MaterialLaw EoS = null;
        Func<double[], double, double>[] scalarFunction = null;
        int NumberOfReactants;

        /// <summary>
        /// Ctor for variable density flows,
        /// i.e. low Mach number flows and
        /// multiphase flows with smooth interface.
        /// </summary>
        /// <param name="SpatDim">
        /// Spatial dimension (either 2 or 3).
        /// </param>
        /// <param name="_bcmap"></param>
        /// <param name="_component"></param>
        /// <param name="EoS">
        /// Material law for variable density flows.
        /// </param>
        public LinearizedConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, MaterialLaw EoS, int NumberOfReactants = -1)
            : this(SpatDim, _bcmap, _component) {

            //m_VariableDensity = true;
            this.EoS = EoS;
            this.NumberOfReactants = NumberOfReactants;

            switch (_bcmap.PhysMode) {
                case PhysicsMode.LowMach:
                    scalarFunction = m_bcmap.bndFunction[VariableNames.Temperature];
                    m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim),
                        VariableNames.Velocity0MeanVector(SpatDim),
                        VariableNames.Temperature0,
                        VariableNames.Temperature0Mean);
                    break;
                case PhysicsMode.Multiphase:
                    scalarFunction = m_bcmap.bndFunction[VariableNames.LevelSet];
                    m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim),
                        VariableNames.Velocity0MeanVector(SpatDim),
                        VariableNames.Phi0,
                        VariableNames.Phi0Mean);
                    break;
                case PhysicsMode.Combustion:
                    if (NumberOfReactants == -1)
                        throw new ArgumentException("NumberOfReactants needs to be specified!");
                    m_ParameterOrdering = ArrayTools.Cat(
                        VariableNames.Velocity0Vector(SpatDim),
                        VariableNames.Velocity0MeanVector(SpatDim),
                        VariableNames.Temperature0,
                        VariableNames.MassFractions0(NumberOfReactants),
                        VariableNames.Temperature0Mean,
                        VariableNames.MassFractionsMean(NumberOfReactants));
                    break;
                case PhysicsMode.Viscoelastic:
                    throw new ApplicationException("Using of wrong constructor for viscoelastic flows.");
                case PhysicsMode.Incompressible:
                    throw new ApplicationException("Using of wrong constructor for incompressible flows.");
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// flux at the boundary
        /// </summary>
        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.NavierSlip_localized:
                case IncompressibleBcType.Velocity_Inlet: {

                        // Fluss am Rand: f(u[d]) = n∙v∙u[d]
                        // wobei n der Normalenvektor, v=(v1,v2) resp. v=(v1,v2,v3) der Linearisierungspunkt.
                        //
                        // Begründung: im Gegensatz zu obigem Code scheint dies besser zu funktionieren,
                        // wenn ein Offset (m_UseBoundaryVelocityParameter == true) addiert wird.
                        // Details: siehe Note 0022;

                        double r = 0.0;
                        double v1, v2, v3 = 0.0, u_d;

                        if (m_UseBoundaryVelocityParameter) {

                            Debug.Assert(m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Combustion || m_bcmap.PhysMode == PhysicsMode.Multiphase, "A boundary velocity is not implemented for variable density!");

                            u_d = Uin[0];

                            v1 = velFunction[inp.EdgeTag, 0](inp.X, inp.time) + inp.Parameters_IN[0 + 2 * m_SpatialDimension];
                            v2 = velFunction[inp.EdgeTag, 1](inp.X, inp.time) + inp.Parameters_IN[1 + 2 * m_SpatialDimension];
                            if (m_SpatialDimension == 3)
                                v3 = velFunction[inp.EdgeTag, 2](inp.X, inp.time) + inp.Parameters_IN[2 + 2 * m_SpatialDimension];

                            r += u_d * (v1 * inp.Normale[0] + v2 * inp.Normale[1]);
                            if (m_SpatialDimension == 3) {
                                r += u_d * v3 * inp.Normale[2];
                            }
                        }
                        else {
                            // Setup params
                            // ============
                            Foundation.CommonParams inp2;
                            inp2.GridDat = inp.GridDat;
                            inp2.Normale = inp.Normale;
                            inp2.iEdge = inp.iEdge;
                            inp2.Parameters_IN = inp.Parameters_IN;
                            inp2.X = inp.X;
                            inp2.time = inp.time;

                            // Dirichlet value for velocity
                            // ============================
                            double Uout = velFunction[inp.EdgeTag, m_component](inp.X, inp.time);


                            // Specify Parameters_OUT
                            // ======================
                            inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                            // Outer values for Velocity and VelocityMean
                            for (int j = 0; j < m_SpatialDimension; j++) {

                                inp2.Parameters_OUT[j] = velFunction[inp.EdgeTag, j](inp.X, inp.time);

                                // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                                inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                            }

                            // Outer values for Scalar and ScalarMean
                            switch (m_bcmap.PhysMode) {
                                case PhysicsMode.Viscoelastic:
                                case PhysicsMode.Incompressible:
                                    break;
                                case PhysicsMode.LowMach: 
                                case PhysicsMode.Multiphase: {
                                        // opt1:
                                        switch (edgeType) {
                                            case IncompressibleBcType.Velocity_Inlet:
                                            case IncompressibleBcType.Wall:
                                                inp2.Parameters_OUT[2 * m_SpatialDimension] = scalarFunction[inp.EdgeTag](inp.X, inp.time);
                                                break;
                                            case IncompressibleBcType.NoSlipNeumann:
                                                inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                                break;
                                            default:
                                                throw new ApplicationException();
                                        }
                                        // opt2:
                                        // Inner values are used for the Scalar variable (even at Dirichlet boundaries of the Scalar variable).                                
                                        // The Dirichlet value for the Scalar variable will be used while solving the Scalar equation, but not in the momentum equation.
                                        //inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                        // Use inner value for ScalarMean, i.e. LambdaIn is used.
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + 1] = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                                        break;
                                    }
                                case PhysicsMode.Combustion: {
                                        switch (edgeType) {
                                            case IncompressibleBcType.Velocity_Inlet:
                                                // opt1: (using Dirichlet values)
                                                inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                                for (int n = 1; n < NumberOfReactants + 1; n++) {
                                                    // opt1: (using Dirichlet values)
                                                    inp2.Parameters_OUT[2 * m_SpatialDimension + n] = m_bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, inp.time);
                                                }
                                                break;
                                            case IncompressibleBcType.Wall:
                                                // opt1: (using Dirichlet values for the temperature)
                                                inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                                for (int n = 1; n < NumberOfReactants + 1; n++) {
                                                    // using inner values for the mass the mass fractions
                                                    inp2.Parameters_OUT[2 * m_SpatialDimension + n] = inp2.Parameters_IN[2 * m_SpatialDimension + n];
                                                }
                                                break;
                                            case IncompressibleBcType.NoSlipNeumann:
                                                for (int n = 0; n < NumberOfReactants + 1; n++) {
                                                    // using inner values
                                                    inp2.Parameters_OUT[2 * m_SpatialDimension + n] = inp2.Parameters_IN[2 * m_SpatialDimension + n];
                                                }
                                                break;
                                            default:
                                                throw new ApplicationException();
                                        }
                                        for (int n = 0; n < NumberOfReactants + 1; n++) {
                                            // Use inner value for mean scalar input parameters, i.e. LambdaIn is used.
                                            inp2.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                                        }
                                        break;
                                    }
                                default:
                                    throw new NotImplementedException("PhysicsMode not implemented");
                            }

                            // Calculate BorderEdgeFlux as InnerEdgeFlux
                            // =========================================
                            r = InnerEdgeFlux(ref inp2, Uin, new double[] { Uout });
                        }

                        return r;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                        double r = 0.0;
                        double u1, u2, u3 = 0, u_d;

                        if (m_UseBoundaryVelocityParameter) {

                            Debug.Assert(m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Combustion || m_bcmap.PhysMode == PhysicsMode.Multiphase, "A boundary velocity is not implemented for variable density!");

                            u_d = Uin[0];
                            u1 = inp.Parameters_IN[0] + inp.Parameters_IN[0 + 2 * m_SpatialDimension];
                            u2 = inp.Parameters_IN[1] + inp.Parameters_IN[1 + 2 * m_SpatialDimension];
                            if (m_SpatialDimension == 3)
                                u3 = inp.Parameters_IN[2] + inp.Parameters_IN[2 + 2 * m_SpatialDimension];
                        }
                        else {
                            u_d = Uin[0];
                            u1 = inp.Parameters_IN[0];
                            u2 = inp.Parameters_IN[1];
                            if (m_SpatialDimension == 3)
                                u3 = inp.Parameters_IN[2];
                        }

                        r += u_d * (u1 * inp.Normale[0] + u2 * inp.Normale[1]);
                        if (m_SpatialDimension == 3) {
                            r += u_d * u3 * inp.Normale[2];
                        }

                        if (m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Multiphase) {
                            double rho = EoS.GetDensity(inp.Parameters_IN[2 * m_SpatialDimension]);
                            r *= rho;
                        }

                        if (m_bcmap.PhysMode == PhysicsMode.Combustion) {
                            double[] args = new double[NumberOfReactants + 1];
                            for (int n = 0; n < NumberOfReactants + 1; n++) {
                                args[n] = inp.Parameters_IN[2 * m_SpatialDimension + n];
                            }
                            double rho = EoS.GetDensity(args);
                            r *= rho;
                        }

                        return r;
                    }
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }
        }

        /// <summary>
        /// bla bla bla
        /// </summary>
        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            double rhoIn = 1.0;
            double rhoOut = 1.0;

            switch (m_bcmap.PhysMode) {
                case PhysicsMode.Viscoelastic:
                case PhysicsMode.Incompressible:
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    rhoIn = EoS.GetDensity(inp.Parameters_IN[2 * m_SpatialDimension]);
                    rhoOut = EoS.GetDensity(inp.Parameters_OUT[2 * m_SpatialDimension]);
                    break;
                case PhysicsMode.Combustion:
                    double[] args_IN = new double[NumberOfReactants + 1];
                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                        args_IN[n] = inp.Parameters_IN[2 * m_SpatialDimension + n];
                    }
                    double[] args_OUT = new double[NumberOfReactants + 1];
                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                        args_OUT[n] = inp.Parameters_OUT[2 * m_SpatialDimension + n];
                    }
                    rhoIn = EoS.GetDensity(args_IN);
                    rhoOut = EoS.GetDensity(args_OUT);
                    break;
                default:
                    throw new NotImplementedException("PhysicsMode not implemented");
            }

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
            if (m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
            }

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for (int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
            }

            double LambdaIn;
            double LambdaOut;

            switch (m_bcmap.PhysMode) {
                case PhysicsMode.Viscoelastic:
                case PhysicsMode.Incompressible:
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, true);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, true);
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    double TemperatureMeanIn = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                    double TemperatureMeanOut = inp.Parameters_OUT[2 * m_SpatialDimension + 1];
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, EoS, true, TemperatureMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, EoS, true, TemperatureMeanOut);
                    break;
                case PhysicsMode.Combustion:
                    double[] ScalarMeanIn = new double[NumberOfReactants + 1];
                    double[] ScalarMeanOut = new double[NumberOfReactants + 1];
                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                        ScalarMeanIn[n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                        ScalarMeanOut[n] = inp.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                    }
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, EoS, true, ScalarMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, EoS, true, ScalarMeanOut);
                    break;
                default:
                    throw new NotImplementedException();
            }

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = Uin[0] - Uout[0];

            r += Lambda * uJump * LaxFriedrichsSchemeSwitch;

            r *= 0.5;
            return r;
        }

        /// <summary>
        /// returns
        /// \f[ 
        ///   \vec{v} \cdot u_d,
        /// \f]
        /// where \f$ \vec{v}\f$  is the linearization point.
        /// For variable density the result is multiplied by \f$ \rho\f$ .
        /// </summary>
        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            output[0] = U[0] * inp.Parameters[0];
            output[1] = U[0] * inp.Parameters[1];
            if (m_SpatialDimension == 3) {
                output[2] = U[0] * inp.Parameters[2];
            }

            if (m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Multiphase) {

                double rho = EoS.GetDensity(inp.Parameters[2 * m_SpatialDimension]);
                for (int d = 0; d < m_SpatialDimension; d++)
                    output[d] *= rho;
            }

            if (m_bcmap.PhysMode == PhysicsMode.Combustion) {
                double[] args = new double[NumberOfReactants + 1];
                for (int n = 0; n < NumberOfReactants + 1; n++) {
                    args[n] = inp.Parameters[2 * m_SpatialDimension + n];
                }
                double rho = EoS.GetDensity(args);
                for (int d = 0; d < m_SpatialDimension; d++)
                    output[d] *= rho;
            }


        }

        /// <summary>
        /// name of the <em>d</em>-th velocity component
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Velocity_d(m_component) }; }
        }

        string[] m_ParameterOrdering;

        /// <summary>
        /// Parameters are the velocity vector (see <see cref="VariableNames.VelocityVector"/>),
        /// concatenated with the mean velocity vector (see <see cref="VariableNames.Velocity0MeanVector"/>) and
        /// optionally (controlled by parameter in constructor) concatenated with
        /// scalar (see <see cref="VariableNames.Phi0"/>) and scalar mean (see <see cref="VariableNames.Phi0Mean"/>) for variable density flows and
        /// optionally (controlled by parameter in constructor) concatenated with the 
        /// boundary velocity offset vector (see <see cref="VariableNames.BoundaryVelocityVector"/>);
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }
    }
}