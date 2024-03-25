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
using ilPSP;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Linearized convective operator (used e.g. in the SIMPLE-algorithm) according to:
    /// K. Shahbazi, P. F. Fischer, C. R. Ethier,
    /// A high-order discontinuous Galerkin method for the unsteady incompressible Navier-Stokes equations,
    /// J. Comput. Phys. 222 (2007) 391–407.
    /// (see BoSSS\doc\notes\0010-BoSSS-References\ShahbaziEtAl2007)
    /// </summary>
    //[Obsolete]
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
            for(int d = 0; d < SpatDim; d++)
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
            if(m_UseBoundaryVelocityParameter)
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

            switch(_bcmap.PhysMode) {
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
                if(NumberOfReactants == -1)
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

            switch(edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.Velocity_Inlet: {

                    // Fluss am Rand: f(u[d]) = n∙v∙u[d]
                    // wobei n der Normalenvektor, v=(v1,v2) resp. v=(v1,v2,v3) der Linearisierungspunkt.
                    //
                    // Begründung: im Gegensatz zu obigem Code scheint dies besser zu funktionieren,
                    // wenn ein Offset (m_UseBoundaryVelocityParameter == true) addiert wird.
                    // Details: siehe Note 0022;

                    double r = 0.0;
                    double v1, v2, v3 = 0.0, u_d;

                    if(m_UseBoundaryVelocityParameter) {

                        Debug.Assert(m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Combustion || m_bcmap.PhysMode == PhysicsMode.Multiphase, "A boundary velocity is not implemented for variable density!");

                        u_d = Uin[0];

                        v1 = velFunction[inp.EdgeTag, 0](inp.X, inp.time) + inp.Parameters_IN[0 + 2 * m_SpatialDimension];
                        v2 = velFunction[inp.EdgeTag, 1](inp.X, inp.time) + inp.Parameters_IN[1 + 2 * m_SpatialDimension];
                        if(m_SpatialDimension == 3)
                            v3 = velFunction[inp.EdgeTag, 2](inp.X, inp.time) + inp.Parameters_IN[2 + 2 * m_SpatialDimension];

                        r += u_d * (v1 * inp.Normal[0] + v2 * inp.Normal[1]);
                        if(m_SpatialDimension == 3) {
                            r += u_d * v3 * inp.Normal[2];
                        }
                    } else {
                        // Setup params
                        // ============
                        Foundation.CommonParams inp2;
                        inp2.GridDat = inp.GridDat;
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;
                        inp2.jCellIn = inp.jCellIn;
                        inp2.jCellOut = int.MinValue;
                        inp2.EdgeTag = inp.EdgeTag;

                        // Dirichlet value for velocity
                        // ============================
                        double Uout = velFunction[inp.EdgeTag, m_component](inp.X, inp.time);

                        // Specify Parameters_OUT
                        // ======================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Outer values for Velocity and VelocityMean
                        for(int j = 0; j < m_SpatialDimension; j++) {

                            inp2.Parameters_OUT[j] = velFunction[inp.EdgeTag, j](inp.X, inp.time);

                            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                            inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                        }

                        // Outer values for Scalar and ScalarMean
                        switch(m_bcmap.PhysMode) {
                            case PhysicsMode.Viscoelastic:
                            case PhysicsMode.Incompressible:
                            case PhysicsMode.RANS:
                            break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase: {
                                // opt1:
                                switch(edgeType) {
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
                                switch(edgeType) {
                                    case IncompressibleBcType.Velocity_Inlet:
                                    // opt1: (using Dirichlet values)
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    for(int n = 1; n < NumberOfReactants + 1; n++) {
                                        // opt1: (using Dirichlet values)
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = m_bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, inp.time);
                                    }
                                    break;
                                    case IncompressibleBcType.Wall:
                                    // opt1: (using Dirichlet values for the temperature)
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    for(int n = 1; n < NumberOfReactants + 1; n++) {
                                        // using inner values for the mass fractions
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = inp2.Parameters_IN[2 * m_SpatialDimension + n];
                                    }
                                    break;
                                    case IncompressibleBcType.NoSlipNeumann:
                                    for(int n = 0; n < NumberOfReactants + 1; n++) {
                                        // using inner values for the temperature and the mass fractions
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = inp2.Parameters_IN[2 * m_SpatialDimension + n];
                                    }
                                    break;
                                    default:
                                    throw new ApplicationException();
                                }
                                for(int n = 0; n < NumberOfReactants + 1; n++) {
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
                        r = InnerEdgeFlux_impl(ref inp2, Uin, new double[] { Uout });
                    }

                    return r;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                    double r = 0.0;
                    double u1, u2, u3 = 0, u_d;

                    if(m_UseBoundaryVelocityParameter) {

                        Debug.Assert(m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Combustion || m_bcmap.PhysMode == PhysicsMode.Multiphase, "A boundary velocity is not implemented for variable density!");

                        u_d = Uin[0];
                        u1 = inp.Parameters_IN[0] + inp.Parameters_IN[0 + 2 * m_SpatialDimension];
                        u2 = inp.Parameters_IN[1] + inp.Parameters_IN[1 + 2 * m_SpatialDimension];
                        if(m_SpatialDimension == 3)
                            u3 = inp.Parameters_IN[2] + inp.Parameters_IN[2 + 2 * m_SpatialDimension];
                    } else {
                        u_d = Uin[0];
                        u1 = inp.Parameters_IN[0];
                        u2 = inp.Parameters_IN[1];
                        if(m_SpatialDimension == 3)
                            u3 = inp.Parameters_IN[2];
                    }

                    r += u_d * (u1 * inp.Normal[0] + u2 * inp.Normal[1]);
                    if(m_SpatialDimension == 3) {
                        r += u_d * u3 * inp.Normal[2];
                    }

                    if(m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Multiphase) {
                        double rho = EoS.GetDensity(inp.Parameters_IN[2 * m_SpatialDimension]);
                        r *= rho;
                    }

                    if(m_bcmap.PhysMode == PhysicsMode.Combustion) {
                        double[] args = new double[NumberOfReactants + 1];
                        for(int n = 0; n < NumberOfReactants + 1; n++) {
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
            return InnerEdgeFlux_impl(ref inp, Uin, Uout);
        }

        /// <summary>
        /// implemented in a separate function to prevent from overloading, i.e. make sure that <see cref="BorderEdgeFlux"/> calls this implementation, not an overloaded one
        /// </summary>
        protected double InnerEdgeFlux_impl(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            double rhoIn = 1.0;
            double rhoOut = 1.0;

            switch(m_bcmap.PhysMode) {
                case PhysicsMode.Viscoelastic:
                case PhysicsMode.Incompressible:
                case PhysicsMode.RANS:
                break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                rhoIn = EoS.GetDensity(inp.Parameters_IN[2 * m_SpatialDimension]);
                rhoOut = EoS.GetDensity(inp.Parameters_OUT[2 * m_SpatialDimension]);
                break;
                case PhysicsMode.Combustion:
                //double[] args_IN = new double[NumberOfReactants + 1];
                //for(int n = 0; n < NumberOfReactants + 1; n++) {
                //    args_IN[n] = inp.Parameters_IN[2 * m_SpatialDimension + n];
                //}
                //double[] args_OUT = new double[NumberOfReactants + 1];
                //for(int n = 0; n < NumberOfReactants + 1; n++) {
                //    args_OUT[n] = inp.Parameters_OUT[2 * m_SpatialDimension + n];
                //}
                //rhoIn = EoS.GetDensity(args_IN);
                //rhoOut = EoS.GetDensity(args_OUT);
                break;
                default:
                throw new NotImplementedException("PhysicsMode not implemented");
            }

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normal[0] + inp.Parameters_OUT[1] * inp.Normal[1]);
            if(m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normal[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for(int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
            }

            double LambdaIn;
            double LambdaOut;

            switch(m_bcmap.PhysMode) {
                case PhysicsMode.Combustion:
                case PhysicsMode.Viscoelastic:
                case PhysicsMode.Incompressible:
                case PhysicsMode.RANS:
                LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, true);
                LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, true);
                break;

                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                double TemperatureMeanIn = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                double TemperatureMeanOut = inp.Parameters_OUT[2 * m_SpatialDimension + 1];
                LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, true, TemperatureMeanIn);
                LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, true, TemperatureMeanOut);
                break;

                //case PhysicsMode.Combustion:
                //double[] ScalarMeanIn = new double[NumberOfReactants + 1];
                //double[] ScalarMeanOut = new double[NumberOfReactants + 1];
                //for(int n = 0; n < NumberOfReactants + 1; n++) {
                //    ScalarMeanIn[n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                //    ScalarMeanOut[n] = inp.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                //}
                //LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, true, ScalarMeanIn);
                //LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, true, ScalarMeanOut);
                //break;
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
            if(m_SpatialDimension == 3) {
                output[2] = U[0] * inp.Parameters[2];
            }

            if(m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Multiphase) {

                double rho = EoS.GetDensity(inp.Parameters[2 * m_SpatialDimension]);
                for(int d = 0; d < m_SpatialDimension; d++)
                    output[d] *= rho;
            }

            if(m_bcmap.PhysMode == PhysicsMode.Combustion) {
                double[] args = new double[NumberOfReactants + 1];
                for(int n = 0; n < NumberOfReactants + 1; n++) {
                    args[n] = inp.Parameters[2 * m_SpatialDimension + n];
                }
                double rho = EoS.GetDensity(args);
                for(int d = 0; d < m_SpatialDimension; d++)
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

    public class LinearizedConvectionJacobi : IVolumeForm, IEdgeForm, ISupportsJacobianComponent, IEquationComponentCoefficient {

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

        int argumentIndex;
        /// <summary>
        /// Ctor for common part of incompressible and low Mach number flows.
        /// </summary>
        /// <param name="SpatDim"></param>
        /// <param name="_bcmap"></param>
        /// <param name="_component"></param>
        private LinearizedConvectionJacobi(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;
            argumentIndex = _component;


            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for(int d = 0; d < SpatDim; d++)
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
        public LinearizedConvectionJacobi(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, bool UseBoundaryVelocityParameter = false)
            : this(SpatDim, _bcmap, _component) {

            m_UseBoundaryVelocityParameter = UseBoundaryVelocityParameter;
            //m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim));
            m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0MeanVector(SpatDim));
            if(m_UseBoundaryVelocityParameter)
                m_ParameterOrdering = ArrayTools.Cat(m_ParameterOrdering, VariableNames.BoundaryVelocityVector(SpatDim));
            m_ArgumentOrdering = VariableNames.VelocityVector(SpatDim); // VelocityX,VelocityY,(VelocityZ) as variables. 

        }

        //bool m_VariableDensity = false;
        MaterialLaw EoS = null;
        Func<double[], double, double>[] scalarFunction = null;
        int NumberOfReactants;
        //int idx;
      
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
        public LinearizedConvectionJacobi(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, MaterialLaw EoS, int NumberOfComponents = -1)
            : this(SpatDim, _bcmap, _component) {

            //m_VariableDensity = true;
            this.EoS = EoS;
            this.NumberOfReactants = NumberOfComponents;
            //idx = _component; // Velocity-i as argument...
            m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0MeanVector(SpatDim)); // used for computation of penalties

            switch(_bcmap.PhysMode) {
                case PhysicsMode.MixtureFraction:
                scalarFunction = m_bcmap.bndFunction[VariableNames.MixtureFraction];
                m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.MixtureFraction); // VelocityX,VelocityY,(VelocityZ), Temperature as variables. 
                m_ParameterOrdering = m_ParameterOrdering.Cat(new string[] { /*VariableNames.Rho*/ });
                break;
                case PhysicsMode.LowMach:
                scalarFunction = m_bcmap.bndFunction[VariableNames.Temperature];
                m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature); // VelocityX,VelocityY,(VelocityZ), Temperature as variables. 
                break;
                case PhysicsMode.Multiphase: //TODO
                                             //scalarFunction = m_bcmap.bndFunction[VariableNames.LevelSet];
                                             //m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim),
                                             //                                     VariableNames.Velocity0MeanVector(SpatDim),
                                             //                                     VariableNames.Phi0,
                                             //                                     VariableNames.Phi0Mean);
                break;
                case PhysicsMode.Combustion:
                    if(NumberOfComponents == -1)
                        throw new ArgumentException("NumberOfReactants needs to be specified!");
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature, VariableNames.MassFractions(NumberOfComponents)); // VelocityX,VelocityY,(VelocityZ), Temperature  and MassFractions as variables. 
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
        protected virtual double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch(edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.Velocity_Inlet: {


                    // Fluss am Rand: f(u[d]) = n∙v∙u[d]
                    // wobei n der Normalenvektor, v=(v1,v2) resp. v=(v1,v2,v3) der Linearisierungspunkt.
                    //
                    // Begründung: im Gegensatz zu obigem Code scheint dies besser zu funktionieren,
                    // wenn ein Offset (m_UseBoundaryVelocityParameter == true) addiert wird.
                    // Details: siehe Note 0022;

                    double r = 0.0;

                    // Setup params
                    // ============
                    Foundation.CommonParams inp2;
                    inp2.GridDat = inp.GridDat;
                    inp2.Normal = inp.Normal;
                    inp2.iEdge = inp.iEdge;
                    inp2.Parameters_IN = inp.Parameters_IN;
                    inp2.X = inp.X;
                    inp2.time = inp.time;
                    inp2.jCellIn = inp.jCellIn;
                    inp2.jCellOut = int.MinValue;
                    inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];
                    inp2.EdgeTag = inp.EdgeTag;

                    // Dirichlet value for velocity
                    // ============================
                    double[] Uout = new double[Uin.Length];
                    for(int i = 0; i < m_SpatialDimension; i++) {
                        Uout[i] = velFunction[inp.EdgeTag, i](inp.X, inp.time) * VelocityMultiplier;
                    }

                        // Outer values for Scalar and ScalarMean
                        switch(m_bcmap.PhysMode) {
                            case PhysicsMode.Viscoelastic:
                            case PhysicsMode.MixtureFraction: 
                            case PhysicsMode.Incompressible:
                                break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase:
                                //Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                // opt1:
                                switch(edgeType) {
                                    case IncompressibleBcType.Velocity_Inlet:
                                    case IncompressibleBcType.Wall:
                                         Uout[m_SpatialDimension] = scalarFunction[inp.EdgeTag](inp.X, inp.time);
                                        break;
                                    case IncompressibleBcType.NoSlipNeumann:
                                        Uout[m_SpatialDimension] = Uin[m_SpatialDimension];
                                        break;
                                    default:
                                        throw new ApplicationException();
                                }                             
                                break;
                            case PhysicsMode.Combustion: {
                                    switch(edgeType) {
                                        case IncompressibleBcType.Velocity_Inlet:
                                            // opt1: (using Dirichlet values)
                                            Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                            for(int n = 0; n < NumberOfReactants ; n++) {
                                                // opt1: (using Dirichlet values)
                                                Uout[m_SpatialDimension + n + 1] =  m_bcmap.bndFunction[VariableNames.MassFraction_n(n)][inp.EdgeTag](inp.X, inp.time);
                                            }
                                            break;
                                        case IncompressibleBcType.Wall: // Dirichlet for temperature and neumann = 0 for mass Fractions
                                            // opt1: (using Dirichlet values for the temperature)
                                            Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                            for(int n = 1; n < NumberOfReactants +1 ; n++) {
                                                // using inner values for the mass fractions
                                                Uout[m_SpatialDimension + n] = Uin[m_SpatialDimension + n];
                                            }
                                            break;
                                        case IncompressibleBcType.NoSlipNeumann:
                                            for(int n = 0; n < NumberOfReactants + 1 ; n++) {
                                                // using inner values for the temperature and the mass fractions
                                                Uout[m_SpatialDimension + n] = Uin[m_SpatialDimension + n];
                                            }
                                            break;
                                        default:
                                            throw new ApplicationException();
                                    }
                                    break;
                                }
                            default:
                                throw new NotImplementedException("PhysicsMode not implemented");
                        }
                    
                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================
                        r = InnerEdgeFlux_impl(ref inp2, Uin, Uout);
                        return r;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet: {
                    double r = 0.0;
                    double u1, u2, u3 = 0, u_d;

                    u_d = Uin[argumentIndex];
                    u1 = Uin[0];
                    u2 = Uin[1];

                    if(m_SpatialDimension == 3)
                        u3 = Uin[2];

                    r += u_d * (u1 * inp.Normal[0] + u2 * inp.Normal[1]);
                    if(m_SpatialDimension == 3) {
                        r += u_d * u3 * inp.Normal[2];
                    }

                        double rho = 1.0;
                        double[] densityArguments;
                    switch (PhysMode) {
                            case PhysicsMode.Incompressible:
                                break;
                            case PhysicsMode.MixtureFraction:
                                rho = EoS.getDensityFromZ(Uin[inp.D]);
                                break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase:
                                //double rho = EoS.GetDensity(inp.Parameters_IN[2 * m_SpatialDimension]);
                                densityArguments = Uin.GetSubVector(m_SpatialDimension, 1); // Only Temp
                                rho = EoS.GetDensity(densityArguments);
                                break;
                            case PhysicsMode.Combustion:
                                densityArguments = Uin.GetSubVector(m_SpatialDimension, NumberOfReactants + 1); // T,Y0,Y1,Y2 and ??Y3??
                                rho = EoS.GetDensity(densityArguments);
                                break;
                            default:
                                throw new NotImplementedException();
                        }
                        r *= rho;

                    return r;
                }
                default:
                throw new NotImplementedException("Boundary condition not implemented!");
            }
        }

        /// <summary>
        /// bla bla bla
        /// </summary>
        protected virtual double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return InnerEdgeFlux_impl(ref inp, Uin, Uout);
        }

        /// <summary>
        /// implemented in a separate function to prevent from overloading, i.e. make sure that <see cref="BorderEdgeFlux"/> calls this implementation, not an overloaded one
        /// </summary>
        protected double InnerEdgeFlux_impl(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            double rhoIn = 1.0;
            double rhoOut = 1.0;
            double[] DensityArgumentsIn;
            double[] DensityArgumentsOut;
            switch(m_bcmap.PhysMode) {
                case PhysicsMode.Viscoelastic:
                case PhysicsMode.Incompressible:
                // Constant density
                break;
                case PhysicsMode.MixtureFraction:
                // opt1:
                rhoIn = EoS.getDensityFromZ(Uin[inp.D]);
                rhoOut = EoS.getDensityFromZ(Uout[inp.D]);
                break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                DensityArgumentsIn = Uin.GetSubVector(m_SpatialDimension, 1);
                DensityArgumentsOut = Uout.GetSubVector(m_SpatialDimension, 1);
                rhoIn = EoS.GetDensity(DensityArgumentsIn);
                rhoOut = EoS.GetDensity(DensityArgumentsOut);
                break;
                case PhysicsMode.Combustion:
                    DensityArgumentsIn = Uin.GetSubVector(m_SpatialDimension, (NumberOfReactants) + 1); 
                    DensityArgumentsOut = Uout.GetSubVector(m_SpatialDimension, (NumberOfReactants) + 1);

                rhoIn = EoS.GetDensity(DensityArgumentsIn);
                rhoOut = EoS.GetDensity(DensityArgumentsOut);
                break;
                default:
                throw new NotImplementedException("PhysicsMode not implemented");
            }

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            int idx = argumentIndex;

            r += rhoIn * Uin[idx] * (Uin[0] * inp.Normal[0] + Uin[1] * inp.Normal[1]);
            r += rhoOut * Uout[idx] * (Uout[0] * inp.Normal[0] + Uout[1] * inp.Normal[1]);
            if(m_SpatialDimension == 3) {
                r += rhoIn * Uin[idx] * Uin[2] * inp.Normal[2] + rhoOut * Uout[idx] * Uout[2] * inp.Normal[2];
            }


            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for(int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[d];// Uin[d];
                VelocityMeanOut[d] = inp.Parameters_OUT[d];// Uout[d];
            }

            double LambdaIn;
            double LambdaOut;

            switch(m_bcmap.PhysMode) {
                case PhysicsMode.Viscoelastic:
                case PhysicsMode.Incompressible:
                LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, true);
                LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, true);
                break;
                case PhysicsMode.MixtureFraction:
                rhoIn = EoS.getDensityFromZ(Uin[inp.D]);
                rhoOut = EoS.getDensityFromZ(Uout[inp.D]);
                LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, true, rhoIn);
                LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, true, rhoOut);
                break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                double TemperatureMeanIn = Uin[m_SpatialDimension];  /* inp.Parameters_IN[2 * m_SpatialDimension + 1];*/
                double TemperatureMeanOut = Uout[m_SpatialDimension]; // inp.Parameters_OUT[2 * m_SpatialDimension + 1];
                LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, true, TemperatureMeanIn);
                LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, true, TemperatureMeanOut);
                break;
                case PhysicsMode.Combustion:
                    double[] ScalarMeanIn = Uin.GetSubVector(m_SpatialDimension, (NumberOfReactants ) + 1);
                    double[] ScalarMeanOut = Uout.GetSubVector(m_SpatialDimension, (NumberOfReactants ) + 1);
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, true, ScalarMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, true, ScalarMeanOut);
                    break;
                default:
                throw new NotImplementedException();
            }

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = Uin[idx] - Uout[idx];
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
        protected virtual void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            int idx = m_component;
            double rho;
            double[] DensityArguments;
            output[0] = U[idx] * U[0];
            output[1] = U[idx] * U[1];
            if(m_SpatialDimension == 3) {
                output[2] = U[idx] * U[2];
            }

            switch(PhysMode) {
                case PhysicsMode.Incompressible:
                rho = 1.0;
                break;
                case PhysicsMode.MixtureFraction:
                rho = EoS.getDensityFromZ(U[inp.D]);
                break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                DensityArguments = U.GetSubVector(m_SpatialDimension, 1); // Only Temperature
                rho = EoS.GetDensity(DensityArguments);
                break;
                case PhysicsMode.Combustion:
                    DensityArguments = U.GetSubVector(m_SpatialDimension, (NumberOfReactants ) + 1); // Temperature, Y0,Y1 and Y2. Y3 is calculated in GetDensity.
                    rho = EoS.GetDensity(DensityArguments);
                    break;
                default:
                throw new NotImplementedException("not implemented physmode");
            }

            for(int d = 0; d < m_SpatialDimension; d++)
                output[d] *= rho;

        }




        string[] m_ArgumentOrdering;

        /// <summary>
        /// name of the <em>d</em>-th velocity component
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        string[] m_ParameterOrdering;

        /// <summary>
        ///
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return this.InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return this.BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = GradV.Length;
            double acc = 0;
            double[] buf = new double[D];
            this.Flux(ref cpv, U, buf);
            for(int d = 0; d < D; d++)
                acc += buf[d] * GradV[d];
            return -acc;
        }

        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }

        public virtual void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("VelocityMultiplier"))
                VelocityMultiplier = (double)cs.UserDefinedValues["VelocityMultiplier"];
        }

        /// <summary>
        /// Multiplier used with every velocity BC
        /// </summary>
        double VelocityMultiplier = 1.0;


        /// <summary>
        /// <see cref="IEdgeForm.BoundaryEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// <see cref="IEdgeForm.InnerEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// <see cref="IVolumeForm.VolTerms"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV | TermActivationFlags.GradV;
            }
        }

    }
}
