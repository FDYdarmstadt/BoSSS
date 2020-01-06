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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;

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
    


    public class LinearizedScalarConvection2Jacobi : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {

        int m_SpatialDimension;
        IncompressibleBoundaryCondMap m_bcmap;

        string Argument;
        int NumberOfReactants;
        int argumentIndex;

        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;

        MaterialLaw EoS;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="SpatDim">Spatial dimension (either 2 or 3)</param>
        /// <param name="BcMap"></param>
        /// <param name="EoS">Null for multiphase. Has to be given for Low-Mach and combustion to calculate density.</param>
        /// <param name="Argument">Variable name of the argument (e.g. "Temperature" or "MassFraction0")</param>
        public LinearizedScalarConvection2Jacobi(int SpatDim, int NumberOfReactants, IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS, string Argument = null) {

            this.NumberOfReactants = NumberOfReactants;
            m_SpatialDimension = SpatDim;
            m_bcmap = BcMap;
            int idx = 0; // 0 for temperature, 1 for Y0 and so on... TODO
            argumentIndex =  m_SpatialDimension + idx;

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];

            for (int d = 0; d < SpatDim; d++)
                velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);


            switch (BcMap.PhysMode) {
                case PhysicsMode.Multiphase:
                    this.Argument = VariableNames.LevelSet;
                    m_ArgumentOrdering = new string[] { VariableNames.LevelSet };
                    m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim));
                    break;
                case PhysicsMode.LowMach:

                    m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim),
                        VariableNames.Temperature0, VariableNames.Temperature0Mean);
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature); // VelocityX,VelocityY,(VelocityZ), Temperature as variables. 

                    if(EoS == null)
                        throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                    else
                        this.EoS = EoS;
                    break;
                case PhysicsMode.Combustion: //TODO
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
        double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;
            int idx = argumentIndex; // Valid for Temperature only. has to be changed 
                                     // Calculate central part
                                     // ======================

            // 2 * (rho u_j * scalar) * n_j
            double rhoIn = 0.0;
            double rhoOut = 0.0;
            switch(m_bcmap.PhysMode) {
                case PhysicsMode.LowMach:
                    double[] DensityArgumentsIn = Uin.GetSubVector(m_SpatialDimension, 1);
                    double[] DensityArgumentsOut = Uout.GetSubVector(m_SpatialDimension, 1);
                    rhoIn = EoS.GetDensity(DensityArgumentsIn);
                    rhoOut = EoS.GetDensity(DensityArgumentsOut);
                    break;
                case PhysicsMode.Combustion: {// TODOOOOOOO 
                        double[] args_IN = new double[NumberOfReactants + 1];
                        for(int n = 0; n < NumberOfReactants + 1; n++) {
                            args_IN[n] = inp.Parameters_IN[2 * m_SpatialDimension + n];
                        }
                        double[] args_OUT = new double[NumberOfReactants + 1];
                        for(int n = 0; n < NumberOfReactants + 1; n++) {
                            args_OUT[n] = inp.Parameters_OUT[2 * m_SpatialDimension + n];
                        }
                        rhoIn = EoS.GetDensity(args_IN);
                        rhoOut = EoS.GetDensity(args_OUT);
                        break;
                    }
                default:
                    throw new NotImplementedException("PhysicsMode not implemented");
            }

            r += rhoIn * Uin[idx] * (Uin[0] * inp.Normal[0] + Uin[1] * inp.Normal[1]);
            r += rhoOut * Uout[idx] * (Uout[0] * inp.Normal[0] + Uout[1] * inp.Normal[1]);
            if(m_SpatialDimension == 3) {
                //r += rhoIn * Uin[idx] * (Uin[2] * inp.Normal[2] + Uout[2] * inp.Normal[2]);
                r += rhoIn * Uin[idx] * Uin[2] * inp.Normal[2] + rhoOut * Uout[idx] * Uout[2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================


            //double[] VelocityMeanIn = new double[m_SpatialDimension];
            //double[] VelocityMeanOut = new double[m_SpatialDimension];
            //for(int d = 0; d < m_SpatialDimension; d++) {
            //    VelocityMeanIn[d] = Uin[d];
            //    VelocityMeanOut[d] = Uout[d];
            //}
            double[] VelocityMeanIn = Uin.GetSubVector(0, m_SpatialDimension); ////////////////////////////TODO CHECK!!!!!!!!!!!!!!!!!!!!
            double[] VelocityMeanOut = Uout.GetSubVector(0, m_SpatialDimension);

            double LambdaIn;
            double LambdaOut;
            switch(m_bcmap.PhysMode) {
                case PhysicsMode.Multiphase:
                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, false);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, false);
                    break;
                case PhysicsMode.LowMach:
                    double TemperatureMeanIn = Uin[m_SpatialDimension];
                    double TemperatureMeanOut = Uout[m_SpatialDimension];

                    LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, false, TemperatureMeanIn);
                    LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, false, TemperatureMeanOut);

                    if(double.IsNaN(LambdaIn) || double.IsInfinity(LambdaIn) || double.IsNaN(LambdaOut) || double.IsInfinity(LambdaOut))
                        throw new NotFiniteNumberException();


                    break;
                case PhysicsMode.Combustion: {
                        //double TemperatureMeanIn = Uin.GetSubVector(m_SpatialDimension, 1);
                        //double TemperatureMeanOut = Uout.GetSubVector(m_SpatialDimension, 1);
                        double[] ScalarMeanIn = new double[NumberOfReactants + 1];
                        double[] ScalarMeanOut = new double[NumberOfReactants + 1];
                        for(int n = 0; n < NumberOfReactants + 1; n++) {
                            ScalarMeanIn[n] = inp.Parameters_IN[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                            ScalarMeanOut[n] = inp.Parameters_OUT[2 * m_SpatialDimension + NumberOfReactants + 1 + n];
                        }
                        LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, EoS, false, ScalarMeanIn);
                        LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, EoS, false, ScalarMeanOut);
                        break;
                    }
                default:
                    throw new NotImplementedException();
            }

            double Lambda = Math.Max(LambdaIn, LambdaOut)*1;

            r += Lambda * (Uin[idx] - Uout[idx]);
            r *= 0.5;
            if(double.IsNaN(r))
                throw new NotFiniteNumberException();

            return r;
        }

        /// <summary>
        /// flux at the boundary
        /// </summary>
        double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];
            switch(edgeType) {
                case IncompressibleBcType.Wall: {

                        if((edgeType == IncompressibleBcType.Wall) && (m_bcmap.PhysMode == PhysicsMode.Multiphase))
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


                        // Boundary values for Parameters
                        // ==============================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Dirichlet value for scalar
                        // ==========================
                        double[] Uout = new double[Uin.Length];
                        for(int i = 0; i < m_SpatialDimension; i++) {
                            Uout[i] = m_bcmap.bndFunction[VariableNames.Velocity_d(i)][inp.EdgeTag](inp.X, inp.time);
                        }
                        Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================                         

                        r = InnerEdgeFlux(ref inp2, Uin, Uout);
                        if(double.IsNaN(r))
                            throw new NotFiniteNumberException();
                        return r;
                    }
                case IncompressibleBcType.Velocity_Inlet: {

                        if((edgeType == IncompressibleBcType.Wall) && (m_bcmap.PhysMode == PhysicsMode.Multiphase))
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

                        // Boundary values for Parameters
                        // ==============================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Velocity
                        for(int j = 0; j < m_SpatialDimension; j++) {
                            // opt1:
                            inp2.Parameters_OUT[j] = m_bcmap.bndFunction[VariableNames.Velocity_d(j)][inp.EdgeTag](inp.X, inp.time);
                            // opt2: inner values
                            //inp2.Parameters_OUT[j] = inp2.Parameters_IN[j];
                            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.                            
                            inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                        }

                        // Skalar (e.g. temperature or MassFraction)
                        switch(m_bcmap.PhysMode) {
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
                                    for(int n = 1; n < NumberOfReactants + 1; n++) {
                                        // opt1: (using Dirichlet values)
                                        inp2.Parameters_OUT[2 * m_SpatialDimension + n] = m_bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, 0);
                                    }
                                    for(int n = 0; n < NumberOfReactants + 1; n++) {
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
                        double[] Uout = new double[Uin.Length];
                        for(int i = 0; i < m_SpatialDimension; i++) {
                            Uout[i] = velFunction[inp.EdgeTag, i](inp.X, inp.time);
                        }
                        Uout[m_SpatialDimension] = m_bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);


                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================    
                        r = InnerEdgeFlux(ref inp2, Uin, Uout);
                        if(double.IsNaN(r))
                            throw new NotFiniteNumberException();
                        return r;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.NoSlipNeumann: {
                        double r = 0.0;
                        double u1, u2, u3 = 0;

                        u1 = Uin[0];
                        u2 = Uin[1];

                        r += Uin[argumentIndex] * (u1 * inp.Normal[0] + u2 * inp.Normal[1]);
                        if(m_SpatialDimension == 3) {
                            u3 = Uin[2];
                            r += Uin[argumentIndex] * u3 * inp.Normal[2];
                        }
                        double rho = 1.0;
                        switch(m_bcmap.PhysMode) {
                            case PhysicsMode.Incompressible:
                                break;
                            case PhysicsMode.LowMach:
                                rho = EoS.GetDensity(Uin[argumentIndex]);
                                break;

                            case PhysicsMode.Combustion: // TODO!
                                double[] args = new double[NumberOfReactants + 1];
                                for(int n = 0; n < NumberOfReactants + 1; n++) {
                                    args[n] = inp.Parameters_IN[2 * m_SpatialDimension + n];
                                }
                                rho = EoS.GetDensity(args);
                                break;

                            default:
                                throw new NotImplementedException("not implemented physmode");
                        }
                        r *= rho;

                        if(double.IsNaN(r))
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
        void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {

            double u = U[0];
            double v = U[1];
            int idx = argumentIndex;
            output[0] = u * U[idx];
            output[1] = v * U[idx];
            if (m_SpatialDimension == 3) {
                double w = U[2];
                output[2] = w * U[idx];
            }


            double rho;
            switch(m_bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    break;
                case PhysicsMode.LowMach:
                    double[] DensityArguments = U.GetSubVector(m_SpatialDimension, 1);

                     rho = EoS.GetDensity(DensityArguments);
                    for(int d = 0; d < m_SpatialDimension; d++)
                        output[d] *= rho;
                    break;
                case PhysicsMode.Combustion:
                    double[] args = new double[NumberOfReactants + 1];
                    for(int n = 0; n < NumberOfReactants + 1; n++) {
                        args[n] = inp.Parameters[2 * m_SpatialDimension + n];
                    }
                     rho = EoS.GetDensity(args);
                    for(int d = 0; d < m_SpatialDimension; d++)
                        output[d] *= rho;
                    break;

                default:
                    throw new NotImplementedException("not implemented physics mode");



            }
            for(int d = 0; d < m_SpatialDimension; d++) {
                if(double.IsNaN(output[d]))
                    throw new NotFiniteNumberException();
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return this.InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return this.BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        double[] buf = null;
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = GradV.Length;
            double acc = 0;
            if(buf == null)
                buf = new double[D];
            this.Flux(ref cpv, U, buf);
            for(int d = 0; d < D; d++)
                acc += buf[d] * GradV[d];
            return -acc;
        }



        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivEdg, DerivVol };
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
                return m_ParameterOrdering;
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
