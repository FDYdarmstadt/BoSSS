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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution.NSECommon {
    //Remark: Density and Lambda computation currently not correct for the PhysicsMode.Multphase! See LinearizedConvection.cs class for the correct switches.

    /// <summary>
    /// Linearized convection for scalar quantities like e.g. temperature for low Mach number flows.
    /// In fact this operator is the analog to <see cref="LinearizedConvection"/> for velocity.
    /// </summary>
    public class LinearizedScalarConvection2 : LinearFlux {
        private int m_SpatialDimension;
        private IncompressibleBoundaryCondMap m_bcmap;

        private string Argument;
        private int NumberOfReactants;

        private string[] m_ArgumentOrdering;
        private string[] m_ParameterOrdering;

        private MaterialLaw EoS;

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

            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normal[0] + inp.Parameters_OUT[1] * inp.Normal[1]);
            if (m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normal[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normal[2];
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
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, false);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, false);
                    break;

                case PhysicsMode.LowMach:
                    double TemperatureMeanIn = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                    double TemperatureMeanOut = inp.Parameters_OUT[2 * m_SpatialDimension + 1];

                    Debug.Assert(TemperatureMeanOut != 0);

                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, false, TemperatureMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, false, TemperatureMeanOut);

                    if (double.IsNaN(LambdaIn) || double.IsInfinity(LambdaIn))
                        throw new NotFiniteNumberException();

                    if (double.IsNaN(LambdaOut) || double.IsInfinity(LambdaOut))
                        throw new NotFiniteNumberException();

                    break;

                case PhysicsMode.Combustion:
                    double[] ScalarMeanIn = new double[NumberOfReactants + 1];
                    double[] ScalarMeanOut = new double[NumberOfReactants + 1];
                    for (int n = 0; n < NumberOfReactants + 1; n++) {
                        ScalarMeanIn[n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                        ScalarMeanOut[n] = inp.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                    }
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, false, ScalarMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, false, ScalarMeanOut);
                    break;

                default:
                    throw new NotImplementedException();
            }

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            r += Lambda * (Uin[0] - Uout[0]);

            r *= 0.5;
            if (double.IsNaN(r))
                throw new NotFiniteNumberException();

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
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;
                        inp2.jCellIn = inp.jCellIn;
                        inp2.jCellOut = int.MinValue;
                        inp2.EdgeTag = inp.EdgeTag;

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
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    // opt2: inner values
                                    //inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                    // Use inner value for TemperatureMean, i.e. LambdaIn is used.
                                    inp2.Parameters_OUT[2 * m_SpatialDimension + 1] = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1: (using Dirichlet values)
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
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
                        if (double.IsNaN(r))
                            throw new NotFiniteNumberException();
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
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;
                        inp2.jCellIn = inp.jCellIn;
                        inp2.jCellOut = int.MinValue;
                        inp2.EdgeTag = inp.EdgeTag;

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
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    // opt2: inner values
                                    //inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                    // Use inner value for TemperatureMean, i.e. LambdaIn is used.
                                    inp2.Parameters_OUT[2 * m_SpatialDimension + 1] = inp.Parameters_IN[2 * m_SpatialDimension + 1];
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1: (using Dirichlet values)
                                    inp2.Parameters_OUT[2 * m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    for (int n = 1; n < NumberOfReactants + 1; n++) {
                                        // opt1: (using Dirichlet values)
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = m_bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, inp.time);
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
                        if (double.IsNaN(r))
                            throw new NotFiniteNumberException();
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

                        r += Uin[0] * (u1 * inp.Normal[0] + u2 * inp.Normal[1]);
                        if (m_SpatialDimension == 3) {
                            u3 = inp.Parameters_IN[2];
                            r += Uin[0] * u3 * inp.Normal[2];
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
                        if (double.IsNaN(r))
                            throw new NotFiniteNumberException();
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
            for (int d = 0; d < m_SpatialDimension; d++) {
                if (double.IsNaN(output[d]))
                    throw new NotFiniteNumberException();
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

    public class LinearizedScalarConvection2Jacobi : IVolumeForm, IEdgeForm, ISupportsJacobianComponent, IEquationComponentCoefficient {
        private int m_SpatialDimension;
        private IncompressibleBoundaryCondMap m_bcmap;

        //string Argument;
        private int NumberOfReactants;

        private int argumentIndex;

        private double LaxFriedrichsSchemeSwitch = 1.0;

        private string[] m_ArgumentOrdering;
        //string[] m_ParameterOrdering;

        private MaterialLaw EoS;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        private bool multiplyBy_cp = false;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="SpatDim">Spatial dimension (either 2 or 3)</param>
        /// <param name="NumberOfReactants"></param>
        /// <param name="BcMap"></param>
        /// <param name="EoS">Null for multiphase. Has to be given for Low-Mach and combustion to calculate density.</param>
        /// <param name="idx">Index of scalar in array, </param>
        /// <param name="Argument">Variable name of the argument (e.g. "Temperature" or "MassFraction0")</param>
        public LinearizedScalarConvection2Jacobi(int SpatDim, int NumberOfReactants, IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS, int idx) {
            this.NumberOfReactants = NumberOfReactants;
            m_SpatialDimension = SpatDim;
            m_bcmap = BcMap;
            argumentIndex = m_SpatialDimension + idx;

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];

            for (int d = 0; d < SpatDim; d++)
                velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);

            //if idx == 0, the heat capacity should multiply this term.
            if (idx == 0)
                multiplyBy_cp = true;

            switch (BcMap.PhysMode) {
                case PhysicsMode.Multiphase:
                    //this.Argument = VariableNames.LevelSet;
                    m_ArgumentOrdering = new string[] { VariableNames.LevelSet };
                    break;

                case PhysicsMode.LowMach:

                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature); // VelocityX,VelocityY,(VelocityZ), Temperature as variables.

                    if (EoS == null)
                        throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                    else
                        this.EoS = EoS;
                    break;

                case PhysicsMode.MixtureFraction:
                    
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.MixtureFraction); // VelocityX,VelocityY,(VelocityZ), MixtureFraction as variables.

                    if (EoS == null)
                        throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                    else
                        this.EoS = EoS;
                    break;

                case PhysicsMode.Combustion:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature, VariableNames.MassFractions(NumberOfReactants)); // u,v,w,T, Y0,Y1,Y2,Y3  as variables (Y4 is calculated as Y4 = 1- (Y0+Y1+Y2+Y3)
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
        private double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;
            int idx = argumentIndex;

            // Calculate central part
            // ======================
            // 2 * (rho u_j * scalar) * n_j
            double rhoIn = 0.0;
            double rhoOut = 0.0;
            double cpIn = 1.0;
            double cpOut = 1.0;

            switch (m_bcmap.PhysMode) {
                case PhysicsMode.MixtureFraction:
                    rhoIn = EoS.getDensityFromZ(Uin[inp.D]);
                    rhoOut = EoS.getDensityFromZ(Uout[inp.D]);
                    break;

                case PhysicsMode.LowMach:
                    double[] DensityArgumentsIn = Uin.GetSubVector(m_SpatialDimension, 1);
                    double[] DensityArgumentsOut = Uout.GetSubVector(m_SpatialDimension, 1);
                    rhoIn = EoS.GetDensity(DensityArgumentsIn);
                    rhoOut = EoS.GetDensity(DensityArgumentsOut);
                    break;

                case PhysicsMode.Combustion:
                    double[] DensityArgumentsIn2 = Uin.GetSubVector(m_SpatialDimension, NumberOfReactants + 1);
                    double[] DensityArgumentsOut2 = Uout.GetSubVector(m_SpatialDimension, NumberOfReactants + 1);
                    rhoIn = EoS.GetDensity(DensityArgumentsIn2);
                    rhoOut = EoS.GetDensity(DensityArgumentsOut2);
                    if (multiplyBy_cp) {
                        cpIn = EoS.GetMixtureHeatCapacity(DensityArgumentsIn2);
                        cpOut = EoS.GetMixtureHeatCapacity(DensityArgumentsOut2);
                    }
                    break;

                default:
                    throw new NotImplementedException("PhysicsMode not implemented");
            }

            r += cpIn * rhoIn * Uin[idx] * (Uin[0] * inp.Normal[0] + Uin[1] * inp.Normal[1]);
            r += cpOut * rhoOut * Uout[idx] * (Uout[0] * inp.Normal[0] + Uout[1] * inp.Normal[1]);
            if (m_SpatialDimension == 3) {
                r += cpIn * rhoIn * Uin[idx] * Uin[2] * inp.Normal[2] + cpOut * rhoOut * Uout[idx] * Uout[2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = Uin.GetSubVector(0, m_SpatialDimension);
            double[] VelocityMeanOut = Uout.GetSubVector(0, m_SpatialDimension);

            double LambdaIn;
            double LambdaOut;
            switch (m_bcmap.PhysMode) {
                case PhysicsMode.MixtureFraction:
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, false, rhoIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, false, rhoOut);
                    break;

                case PhysicsMode.Multiphase:
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, false);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, false);
                    break;

                case PhysicsMode.LowMach:
                    double ScalarMeanIn = Uin[m_SpatialDimension];
                    double ScalarMeanOut = Uout[m_SpatialDimension];
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, false, multiplyBy_cp, ScalarMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, false, multiplyBy_cp, ScalarMeanOut);
                    break;

                case PhysicsMode.Combustion:
                    double[] ScalarMeanIn_vec = Uin.GetSubVector(m_SpatialDimension, NumberOfReactants + 1);
                    double[] ScalarMeanOut_vec = Uout.GetSubVector(m_SpatialDimension, NumberOfReactants + 1);
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, false, multiplyBy_cp, ScalarMeanIn_vec);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, false, multiplyBy_cp, ScalarMeanOut_vec);
                    break;

                default:
                    throw new NotImplementedException();
            }
            if (double.IsNaN(LambdaIn) || double.IsInfinity(LambdaIn) || double.IsNaN(LambdaOut) || double.IsInfinity(LambdaOut))
                throw new NotFiniteNumberException();
            double Lambda = Math.Max(LambdaIn, LambdaOut);

            r += Lambda * (Uin[idx] - Uout[idx]) * LaxFriedrichsSchemeSwitch;
            r *= 0.5;
            if (double.IsNaN(r))
                throw new NotFiniteNumberException();

            return r;
        }

        /// <summary>
        /// flux at the boundary
        /// </summary>
        private double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
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
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;
                        inp2.jCellIn = inp.jCellIn;
                        inp2.jCellOut = int.MinValue;
                        inp2.EdgeTag = inp.EdgeTag;

                        // Boundary values for Parameters
                        // ==============================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Dirichlet values for scalars and Velocities
                        // ==========================
                        double[] Uout = new double[Uin.Length];
                        //Velocity
                        for (int i = 0; i < m_SpatialDimension; i++) {
                            Uout[i] = m_bcmap.bndFunction[VariableNames.Velocity_d(i)][inp.EdgeTag](inp.X, inp.time)* VelocityMultiplier;
                        }

                        switch (m_bcmap.PhysMode) {
                            case PhysicsMode.MixtureFraction:
                                // opt1:
                                //inp2.Parameters_OUT = inp.Parameters_IN;
                                //return 0.0;
                                Uout[m_SpatialDimension] = Uin[m_SpatialDimension]; // m_bcmap.bndFunction[VariableNames.MixtureFraction][inp.EdgeTag](inp.X, inp.time);
                                break;

                            case PhysicsMode.LowMach: {
                                    Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1: (using Dirichlet values)
                                    Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    for (int n = 0; n < NumberOfReactants; n++) {
                                        //Using inner values for species
                                        Uout[m_SpatialDimension + 1 + n] = Uin[m_SpatialDimension + 1 + n];
                                    }
                                    break;
                                }
                            case PhysicsMode.Multiphase:
                                break;

                            default:
                                throw new NotImplementedException();
                        }
                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================

                        r = InnerEdgeFlux(ref inp2, Uin, Uout);
                        if (double.IsNaN(r))
                            throw new NotFiniteNumberException();
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
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;
                        inp2.jCellIn = inp.jCellIn;
                        inp2.jCellOut = int.MinValue;
                        inp2.EdgeTag = inp.EdgeTag;

                        // Boundary values for Parameters
                        // ==============================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];
                        double[] Uout = new double[Uin.Length];

                        // Velocity
                        for (int j = 0; j < m_SpatialDimension; j++) {
                            // opt1:
                            Uout[j] = m_bcmap.bndFunction[VariableNames.Velocity_d(j)][inp.EdgeTag](inp.X, inp.time) * VelocityMultiplier;
                        }

                        // Skalar (e.g. temperature or MassFraction)
                        switch (m_bcmap.PhysMode) {
                            case PhysicsMode.MixtureFraction: {
                                    // opt1:
                                    Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.MixtureFraction][inp.EdgeTag](inp.X, inp.time);
                                    break;
                                }
                            case PhysicsMode.LowMach: {
                                    // opt1:
                                    Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    // opt2: inner values
                                    //inp2.Parameters_OUT[2 * m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension];
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1: (using Dirichlet values)
                                    Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    for (int n = 0; n < NumberOfReactants; n++) {
                                        // opt1: (using Dirichlet values)
                                        Uout[m_SpatialDimension + 1 + n] = m_bcmap.bndFunction[VariableNames.MassFraction_n(n)][inp.EdgeTag](inp.X, inp.time);
                                    }
                                    break;
                                }
                            case PhysicsMode.Multiphase:
                                break;

                            default:
                                throw new NotImplementedException();
                        }

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================
                        r = InnerEdgeFlux(ref inp2, Uin, Uout);
                        if (double.IsNaN(r))
                            throw new NotFiniteNumberException();
                        return r;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet:
                case IncompressibleBcType.NoSlipNeumann: {
                        double r = 0.0;
                        double u1, u2, u3 = 0;

                        u1 = Uin[0];
                        u2 = Uin[1];

                        r += Uin[argumentIndex] * (u1 * inp.Normal[0] + u2 * inp.Normal[1]);
                        if (m_SpatialDimension == 3) {
                            u3 = Uin[2];
                            r += Uin[argumentIndex] * u3 * inp.Normal[2];
                        }
                        double rho = 1.0;
                        double cp = 1.0;
                        switch (m_bcmap.PhysMode) {
                            case PhysicsMode.Incompressible:
                                break;

                            case PhysicsMode.MixtureFraction:
                                rho = EoS.getDensityFromZ(Uin[inp.D]);
                                break;

                            case PhysicsMode.LowMach:
                                rho = EoS.GetDensity(Uin[argumentIndex]);
                                break;

                            case PhysicsMode.Combustion:
                                double[] args = ArrayTools.GetSubVector(Uin, m_SpatialDimension, NumberOfReactants + 1);
                                rho = EoS.GetDensity(args);
                                if (multiplyBy_cp == true)
                                    cp = EoS.GetMixtureHeatCapacity(args);
                                break;

                            default:
                                throw new NotImplementedException("not implemented physmode");
                        }
                        r *= (cp * rho);

                        if (double.IsNaN(r))
                            throw new NotFiniteNumberException();
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
        private void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            double u = U[0];
            double v = U[1];
            output[0] = u * U[argumentIndex];
            output[1] = v * U[argumentIndex];

            if (m_SpatialDimension == 3) {
                double w = U[2];
                output[2] = w * U[argumentIndex];
            }

            double rho;
            double cp = 1.0;
            switch (m_bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    rho = 1.0;
                    break;

                case PhysicsMode.MixtureFraction:
                    rho = EoS.getDensityFromZ(U[inp.D]);
                    break;

                case PhysicsMode.LowMach:
                    double[] DensityArguments = U.GetSubVector(m_SpatialDimension, 1);
                    rho = EoS.GetDensity(DensityArguments);
                    break;

                case PhysicsMode.Combustion:
                    double[] arguments = U.GetSubVector(m_SpatialDimension, NumberOfReactants + 1); // T, Y0,Y1,Y2, Y3
                    rho = EoS.GetDensity(arguments);
                    if (multiplyBy_cp)
                        cp = EoS.GetMixtureHeatCapacity(arguments);
                    break;

                default:
                    throw new NotImplementedException("not implemented physics mode");
            }

            for (int d = 0; d < m_SpatialDimension; d++)
                output[d] *= (cp * rho);

            for (int d = 0; d < m_SpatialDimension; d++) {
                if (double.IsNaN(output[d]))
                    throw new NotFiniteNumberException();
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
            var buf = new double[D];
            this.Flux(ref cpv, U, buf);
            for (int d = 0; d < D; d++)
                acc += buf[d] * GradV[d];
            return -acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivEdg, DerivVol };
        }
        /// <summary>
        /// Multiplier used with every velocity BC
        /// </summary>
        double VelocityMultiplier = 1.0;
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("VelocityMultiplier"))
                VelocityMultiplier = (double)cs.UserDefinedValues["VelocityMultiplier"];
        }

        /// <summary>
        /// Level-Set for multiphase.
        /// Temperature for convection in the temperature equation
        /// Name of the computed MassFraction in MassFraction balances (e.g. MassFraction0)
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        public virtual IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

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