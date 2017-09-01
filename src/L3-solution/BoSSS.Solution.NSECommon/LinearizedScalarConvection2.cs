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
using System.Diagnostics;

namespace BoSSS.Solution.NSECommon {

    //Remark: Density and Lambda computation currently not correct for the PhysicsMode.Multphase! See LinearizedConvection.cs class for the correct switches.

    /// <summary>
    /// Linearized convection for scalar quantities like e.g. temperature for low Mach number flows.
    /// In fact this operator is the analog to <see cref="LinearizedConvection"/> for velocity.
    /// </summary>
    public class LinearizedScalarConvection2 : LinearFlux {

        int m_SpatialDimension;
        IncompressibleBoundaryCondMap m_bcmap;

        string Argument;
        int NumberOfReactants;

        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;

        MaterialLaw EoS;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="SpatDim">Spatial dimension (either 2 or 3)</param>
        /// <param name="BcMap"></param>
        /// <param name="EoS">Null for multiphase. Has to be given for Low-Mach and combustion to calculate density.</param>
        /// <param name="Argument">Variable name of the argument (e.g. "Temperature" or "MassFraction0")</param>
        public LinearizedScalarConvection2(int SpatDim, int NumberOfReactants, IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS, string Argument = null) {

            this.NumberOfReactants = NumberOfReactants;
            m_SpatialDimension = SpatDim;
            m_bcmap = BcMap;

            switch (BcMap.PhysMode) {
                case PhysicsMode.Multiphase:
                    this.Argument = VariableNames.LevelSet;
                    m_ArgumentOrdering = new string[] { VariableNames.LevelSet };
                    m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim));
                    break;
                case PhysicsMode.LowMach:
                    this.Argument = VariableNames.Temperature;
                    m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim),
                        VariableNames.Temperature0, VariableNames.Temperature0Mean);
                    m_ArgumentOrdering = new string[] { VariableNames.Temperature };
                    if (EoS == null)
                        throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                    else
                        this.EoS = EoS;
                    break;
                case PhysicsMode.Combustion:
                    this.Argument = Argument;
                    m_ParameterOrdering = ArrayTools.Cat(
                        VariableNames.Velocity0Vector(SpatDim),
                        VariableNames.Velocity0MeanVector(SpatDim),
                        VariableNames.Temperature0,
                        VariableNames.MassFractions0(NumberOfReactants),
                        VariableNames.Temperature0Mean,
                        VariableNames.MassFractionsMean(NumberOfReactants));
                    m_ArgumentOrdering = new string[] { Argument };
                    if (EoS == null)
                        throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                    else
                        this.EoS = EoS;
                    break;
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// flux at inner edges
        /// </summary>        
        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            // 2 * (rho u_j * scalar) * n_j
            double rhoIn = 0.0;
            double rhoOut = 0.0;
            switch (m_bcmap.PhysMode) {
                case PhysicsMode.LowMach:
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

            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
            if (m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
            }

            // Calculate dissipative part
            // ==========================

            IList<double> _VelocityMeanIn = new List<double>();
            IList<double> _VelocityMeanOut = new List<double>();
            for (int d = 0; d < m_SpatialDimension; d++) {
                _VelocityMeanIn.Add(inp.Parameters_IN[m_SpatialDimension + d]);
                _VelocityMeanOut.Add(inp.Parameters_OUT[m_SpatialDimension + d]);
            }
            double[] VelocityMeanIn = _VelocityMeanIn.ToArray();
            double[] VelocityMeanOut = _VelocityMeanOut.ToArray();

            double LambdaIn;
            double LambdaOut;
            switch (m_bcmap.PhysMode) {
                case PhysicsMode.Multiphase:
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, false);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, false);
                    break;
                case PhysicsMode.LowMach:
                    double TemperatureMeanIn = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                    double TemperatureMeanOut = inp.Parameters_OUT[2 * m_SpatialDimension + 1];
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, EoS, false, TemperatureMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, EoS, false, TemperatureMeanOut);
                    break;
                case PhysicsMode.Combustion:
                    double[] ScalarMeanIn = new double[NumberOfReactants + 1];
                    double[] ScalarMeanOut = new double[NumberOfReactants + 1];
                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                        ScalarMeanIn[n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                        ScalarMeanOut[n] = inp.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                    }
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, EoS, false, ScalarMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, EoS, false, ScalarMeanOut);
                    break;
                default:
                    throw new NotImplementedException();
            }

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            r += Lambda * (Uin[0] - Uout[0]);

            r *= 0.5;
            return r;
        }

        /// <summary>
        /// flux at the boundary
        /// </summary>
        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall: {

                        if ((edgeType == IncompressibleBcType.Wall) && (m_bcmap.PhysMode == PhysicsMode.Multiphase))
                            throw new ApplicationException("Use NoSlipNeumann boundary condition for multiphase flows instead of Wall.");

                        double r = 0.0;

                        // Setup params
                        // ============
                        Foundation.CommonParams inp2;
                        inp2.GridDat = inp.GridDat;
                        inp2.Normale = inp.Normale;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;

                        // Boundary values for Parameters
                        // ==============================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Velocity
                        for (int j = 0; j < m_SpatialDimension; j++) {
                            // opt1:
                            inp2.Parameters_OUT[j] = m_bcmap.bndFunction[VariableNames.Velocity_d(j)][inp.EdgeTag](inp.X, inp.time);
                            // opt2: inner values
                            //inp2.Parameters_OUT[j] = inp2.Parameters_IN[j];
                            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.                            
                            inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                        }

                        // Skalar (e.g. temperature or MassFraction)
                        switch (m_bcmap.PhysMode) {
                            case PhysicsMode.LowMach: {
                                    // opt1:
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, 0);
                                    // opt2: inner values
                                    //inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                    // Use inner value for TemperatureMean, i.e. LambdaIn is used.
                                    inp2.Parameters_OUT[2 * m_SpatialDimension + 1] = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1: (using Dirichlet values)
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, 0);
                                    for (int n = 1; n < NumberOfReactants + 1; n++) {
                                        //Using inner values
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = inp.Parameters_IN[2 * m_SpatialDimension + n];
                                    }
                                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                                        // Use inner value for mean scalar input parameters, i.e. LambdaIn is used.
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                                    }
                                    break;
                                }
                            case PhysicsMode.Multiphase:
                                break;
                            default:
                                throw new NotImplementedException();
                        }

                        // Dirichlet value for scalar
                        // ==========================
                        double[] Uout = new double[] { m_bcmap.bndFunction[Argument][inp.EdgeTag](inp.X, inp.time) };

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================    
                        r = InnerEdgeFlux(ref inp2, Uin, Uout);

                        return r;
                    }
                case IncompressibleBcType.Velocity_Inlet: {

                        if ((edgeType == IncompressibleBcType.Wall) && (m_bcmap.PhysMode == PhysicsMode.Multiphase))
                            throw new ApplicationException("Use NoSlipNeumann boundary condition for multiphase flows instead of Wall.");

                        double r = 0.0;

                        // Setup params
                        // ============
                        Foundation.CommonParams inp2;
                        inp2.GridDat = inp.GridDat;
                        inp2.Normale = inp.Normale;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;

                        // Boundary values for Parameters
                        // ==============================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Velocity
                        for (int j = 0; j < m_SpatialDimension; j++) {
                            // opt1:
                            inp2.Parameters_OUT[j] = m_bcmap.bndFunction[VariableNames.Velocity_d(j)][inp.EdgeTag](inp.X, inp.time);
                            // opt2: inner values
                            //inp2.Parameters_OUT[j] = inp2.Parameters_IN[j];
                            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.                            
                            inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                        }

                        // Skalar (e.g. temperature or MassFraction)
                        switch (m_bcmap.PhysMode) {
                            case PhysicsMode.LowMach: {
                                    // opt1:
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, 0);
                                    // opt2: inner values
                                    //inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                    // Use inner value for TemperatureMean, i.e. LambdaIn is used.
                                    inp2.Parameters_OUT[2 * m_SpatialDimension + 1] = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1: (using Dirichlet values)
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, 0);
                                    for (int n = 1; n < NumberOfReactants + 1; n++) {
                                        // opt1: (using Dirichlet values)
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = m_bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, 0);
                                    }
                                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                                        // Use inner value for mean scalar input parameters, i.e. LambdaIn is used.
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                                    }
                                    break;
                                }
                            case PhysicsMode.Multiphase:
                                break;
                            default:
                                throw new NotImplementedException();
                        }

                        // Dirichlet value for scalar
                        // ==========================
                        double[] Uout = new double[] { m_bcmap.bndFunction[Argument][inp.EdgeTag](inp.X, inp.time) };

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================    
                        r = InnerEdgeFlux(ref inp2, Uin, Uout);

                        return r;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.NoSlipNeumann: {
                        double r = 0.0;
                        double u1, u2, u3 = 0;

                        u1 = inp.Parameters_IN[0];
                        u2 = inp.Parameters_IN[1];

                        r += Uin[0] * (u1 * inp.Normale[0] + u2 * inp.Normale[1]);
                        if (m_SpatialDimension == 3) {
                            u3 = inp.Parameters_IN[2];
                            r += Uin[0] * u3 * inp.Normale[2];
                        }

                        if (m_bcmap.PhysMode == PhysicsMode.LowMach) {
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
        /// returns
        /// \f[ 
        ///   \vec{v} \cdot \phi,
        /// \f]
        /// where \f$ \vec{v}\f$  is the linearization point.
        /// </summary>
        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {

            double u = inp.Parameters[0];
            double v = inp.Parameters[1];

            output[0] = u * U[0];
            output[1] = v * U[0];
            if (m_SpatialDimension == 3) {
                double w = inp.Parameters[2];
                output[2] = w * U[0];
            }
            if (m_bcmap.PhysMode == PhysicsMode.LowMach) {

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
        /// Level-Set for multiphase.
        /// Temperature for convection in the temperature equation
        /// Name of the computed MassFraction in MassFraction balances (e.g. MassFraction0)
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }
    }
}
