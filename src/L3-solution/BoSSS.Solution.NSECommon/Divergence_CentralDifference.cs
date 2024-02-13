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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP.Utils;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Central difference scheme for divergence operator.
    /// </summary>
    public class Divergence_CentralDifference : LinearFlux {

        int Component;
        IncompressibleBoundaryCondMap Bcmap;
        Func<double[], double, double>[] VelFunction;


        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>
        /// <param name="Component"></param>
        /// <param name="Bcmap"></param>
        public Divergence_CentralDifference(int Component, IncompressibleBoundaryCondMap Bcmap) {
            this.Component = Component;
            this.Bcmap = Bcmap;
            this.VelFunction = Bcmap.bndFunction[VariableNames.Velocity_d(Component)];
        }

        MaterialLaw EoS = null;
        string[] m_ParamterOrdering = null;
        int NumberOfReactants = -1;

        /// <summary>
        /// Ctor for low Mach number flows.
        /// </summary>
        /// <param name="Component"></param>
        /// <param name="Bcmap"></param>
        /// <param name="EoS"></param>
        public Divergence_CentralDifference(int Component, IncompressibleBoundaryCondMap Bcmap, MaterialLaw EoS, int NumberOfReactants = -1)
            : this(Component, Bcmap) {
            this.EoS = EoS;
            switch(Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                case PhysicsMode.Multiphase:
                    throw new ApplicationException("Wrong constructor");
                case PhysicsMode.LowMach:
                    this.m_ParamterOrdering =new string[] { VariableNames.Temperature0 };
                    break;
                case PhysicsMode.Combustion:
                    this.m_ParamterOrdering =  ArrayTools.Cat(new string[] { VariableNames.Temperature0 }, VariableNames.MassFractions0(NumberOfReactants));
                    this.NumberOfReactants = NumberOfReactants;
                    if(NumberOfReactants == -1)
                        throw new ArgumentException("NumberOfReactants must be specified for combustion flows.");
                    break;
                default:
                    throw new NotImplementedException();
            }
        }

        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            double res = 0.0;

            IncompressibleBcType edgeType = Bcmap.EdgeTag2Type[inp.EdgeTag];

            switch(edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                        res = 0.0;
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        double TemperatureOut = 0.0;
                        double Uout = VelFunction[inp.EdgeTag](inp.X, inp.time);
                        switch(Bcmap.PhysMode) {
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase: {
                                    // opt1:
                                    TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    res = EoS.GetDensity(TemperatureOut) * Uout * inp.Normal[Component];
                                    // opt2:
                                    //double TemperatureIN = inp.Parameters_IN[0];
                                    //double rhoIn = EoS.GetDensity(TemperatureIN);
                                    //res = rhoIn * Uout * inp.Normale[Component];
                                    break;
                                }
                            case PhysicsMode.Combustion:
                                // opt1:
                                TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                double[] args = new double[NumberOfReactants + 1];
                                args[0] = TemperatureOut;
                                for(int n = 1; n < NumberOfReactants + 1; n++) {
                                    args[n] = Bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, inp.time);
                                }
                                res = EoS.GetDensity(args) * Uout * inp.Normal[Component];
                                break;
                            case PhysicsMode.Incompressible:
                                res = Uout * inp.Normal[Component];
                                break;
                            default:
                                throw new ApplicationException("PhysicsMode not implemented");
                        }
                    }
                    break;
                case IncompressibleBcType.Pressure_Outlet: {
                        switch(Bcmap.PhysMode) {
                            case PhysicsMode.Incompressible:
                                res = Uin[0] * inp.Normal[Component];
                                break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase:
                                res = EoS.GetDensity(inp.Parameters_IN[0]) * Uin[0] * inp.Normal[Component];
                                break;
                            case PhysicsMode.Combustion:
                                res = EoS.GetDensity(inp.Parameters_IN) * Uin[0] * inp.Normal[Component];
                                break;
                            default:
                                throw new ApplicationException("PhysicsMode not implemented");
                        }
                    }
                    break;
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }

            return res;
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double res = 0.0;

            switch(Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    res = 0.5 * (Uin[0] + Uout[0]) * inp.Normal[Component];
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    double densityIn = (EoS.GetDensity(inp.Parameters_IN[0]));
                    double densityOut = (EoS.GetDensity(inp.Parameters_OUT[0]));

                    if(double.IsNaN(densityIn) || double.IsInfinity(densityIn))
                        throw new NotFiniteNumberException();

                    if(double.IsNaN(densityOut) || double.IsInfinity(densityOut))
                        throw new NotFiniteNumberException();

                    res = 0.5 * (densityIn * Uin[0] + densityOut * Uout[0]) * inp.Normal[Component];
                    break;
                case PhysicsMode.Combustion:
                    res = 0.5 * (EoS.GetDensity(inp.Parameters_IN) * Uin[0] + EoS.GetDensity(inp.Parameters_OUT) * Uout[0]) * inp.Normal[Component];
                    break;
                default:
                    throw new ApplicationException("PhysicsMode not implemented");
            }

            return res;
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            int D = output.Length;
            Array.Clear(output, 0, D);

            switch(Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    output[Component] = U[0];
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    output[Component] = EoS.GetDensity(inp.Parameters[0]) * U[0];
                    break;
                case PhysicsMode.Combustion:
                    output[Component] = EoS.GetDensity(inp.Parameters) * U[0];
                    break;
                default:
                    throw new ApplicationException("PhysicsMode not implemented");
            }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(Component) };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return m_ParamterOrdering;
            }
        }
    }

    /// <summary>
    /// Central difference scheme for divergence operator.
    /// </summary>
    public class Divergence_CentralDifferenceJacobian : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, ISupportsJacobianComponent, IEquationComponentCoefficient {

        int Component;
        IncompressibleBoundaryCondMap Bcmap;
        Func<double[], double, double>[] VelFunction;

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>
        /// <param name="Component"></param>
        /// <param name="Bcmap"></param>
        public Divergence_CentralDifferenceJacobian(int Component, IncompressibleBoundaryCondMap Bcmap, int SpatDim) {
            this.Component = Component;
            this.Bcmap = Bcmap;
            this.VelFunction = Bcmap.bndFunction[VariableNames.Velocity_d(Component)];
            m_SpatialDimension = SpatDim;
            m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim));
        }

        MaterialLaw EoS = null;
        int NumberOfSpecies = -1;
        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;
        /// <summary>
        /// Ctor for low Mach number flows.
        /// </summary>
        /// <param name="Component"></param>
        /// <param name="Bcmap"></param>
        /// <param name="EoS"></param>
        public Divergence_CentralDifferenceJacobian(int Component, IncompressibleBoundaryCondMap Bcmap, int SpatDim, MaterialLaw EoS, int NumberOfSpecies = -1)
            : this(Component, Bcmap, SpatDim) {
            this.EoS = EoS;

            switch (Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                case PhysicsMode.Multiphase:
                    throw new ApplicationException("Wrong constructor");
                case PhysicsMode.MixtureFraction:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), new string[] { VariableNames.MixtureFraction});
                    m_ParameterOrdering = null;
                    this.NumberOfSpecies = NumberOfSpecies;
                    break;
                case PhysicsMode.LowMach:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), new string[] { VariableNames.Temperature });
                    break;
                case PhysicsMode.Combustion:
                    if(NumberOfSpecies == -1)
                        throw new ArgumentException("NumberOfSpecies must be specified for combustion flows.");
                    string[] Parameters = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(NumberOfSpecies));
                    this.NumberOfSpecies = NumberOfSpecies;
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), Parameters);
                    break;
                default:
                    throw new NotImplementedException();
            }


        }

        double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            double res = 0.0;
            
            IncompressibleBcType edgeType = Bcmap.EdgeTag2Type[inp.EdgeTag];
            double[] DensityArgumentsIn;
            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                    res = 0.0;
                    break;
                case IncompressibleBcType.Velocity_Inlet: {
                        double TemperatureOut = 0.0;
                        double Uout = VelFunction[inp.EdgeTag](inp.X, inp.time) * VelocityMultiplier;                       

                        switch (Bcmap.PhysMode) {
                            case PhysicsMode.Incompressible:
                                res = Uout * inp.Normal[Component];
                                break;
                            case PhysicsMode.MixtureFraction:
                                // opt1:
                                double Zout = Bcmap.bndFunction[VariableNames.MixtureFraction][inp.EdgeTag](inp.X, inp.time);                                
                                res = EoS.getDensityFromZ(Zout) * Uout * inp.Normal[Component];
                                break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase: {
                                    // opt1:
                                    TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);

                                   res = EoS.GetDensity(TemperatureOut) * Uout * inp.Normal[Component];

                                    // opt2:
                                    //double TemperatureIN = inp.Parameters_IN[0];
                                    //double rhoIn = EoS.GetDensity(TemperatureIN);
                                    //res = rhoIn * Uout * inp.Normale[Component];
                                    break;
                                }
                            case PhysicsMode.Combustion: {
                                    // opt1:
                                    TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    double[] argsa = new double[NumberOfSpecies + 1];
                                    argsa[0] = TemperatureOut;
                                    for (int n = 0; n < NumberOfSpecies; n++) {
                                        argsa[n + 1] = Bcmap.bndFunction[VariableNames.MassFraction_n(n)][inp.EdgeTag](inp.X, inp.time);
                                    }
                                    res = EoS.GetDensity(argsa) * Uout * inp.Normal[Component];


                                    // todo: should i also use central difference flux here?

                                    //double TemperatureIn = Uin[m_SpatialDimension];
                                    //TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    //double[] argsa = new double[NumberOfSpecies + 1];
                                    //double[] argsb= new double[NumberOfSpecies + 1];
                                    //argsa[0] = TemperatureOut;
                                    //argsb[0] = TemperatureIn;
                                    //for (int n = 0; n < NumberOfSpecies; n++) {
                                    //    argsa[n+1] = Bcmap.bndFunction[VariableNames.MassFraction_n(n)][inp.EdgeTag](inp.X, inp.time);
                                    //    argsb[n + 1] = Uin[m_SpatialDimension + n+1];
                                    //}
                                    //res = EoS.GetDensity(argsa) * Uout * inp.Normal[Component] + EoS.GetDensity(argsb)*Uin[Component]*inp.Normal[Component];
                                    //res *= 0.5;
                                    break;
                                }
                            default:
                                throw new ApplicationException("PhysicsMode not implemented");
                        }
                    }
                    break;
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet: {
                        switch (Bcmap.PhysMode) {
                            case PhysicsMode.Incompressible:
                                res = Uin[Component] * inp.Normal[Component];
                                break;
                            case PhysicsMode.MixtureFraction:
                                double Zin = Uin[inp.D];
                                double rho = EoS.getDensityFromZ(Zin);
                                res = rho * Uin[Component] * inp.Normal[Component];
                                break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase:
                                DensityArgumentsIn = Uin.GetSubVector(m_SpatialDimension, 1); // Only TemperatureIn
                                res = EoS.GetDensity(DensityArgumentsIn) * Uin[Component] * inp.Normal[Component];
                                break;
                            case PhysicsMode.Combustion:
                                DensityArgumentsIn = Uin.GetSubVector(m_SpatialDimension, NumberOfSpecies+1);
                                res = EoS.GetDensity(DensityArgumentsIn) * Uin[Component] * inp.Normal[Component];
                                break;
                            default:
                                throw new ApplicationException("PhysicsMode not implemented");
                        }
                    }
                    break;
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }

            return res;
        }

        double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double res = 0.0;
            double[] DensityArguments_In; // Arguments used to calculate the density with the EoS
            double[] DensityArguments_Out; // Arguments used to calculate the density with the EoS
            double densityIn;
            double densityOut;
            switch (Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    res = 0.5 * (Uin[Component] + Uout[Component]) * inp.Normal[Component];
                    break;
                case PhysicsMode.MixtureFraction:
                    densityIn = EoS.getDensityFromZ(Uin[inp.D]);
                    densityOut = EoS.getDensityFromZ(Uout[inp.D]);
                    res = 0.5 * (densityIn * Uin[Component] + densityOut * Uout[Component]) * inp.Normal[Component];
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    DensityArguments_In = Uin.GetSubVector(m_SpatialDimension, 1);
                    DensityArguments_Out = Uout.GetSubVector(m_SpatialDimension, 1);
                    densityIn = (EoS.GetDensity(DensityArguments_In));
                    densityOut = (EoS.GetDensity(DensityArguments_Out));
                    if (double.IsNaN(densityIn) || double.IsInfinity(densityIn) || double.IsNaN(densityOut) || double.IsInfinity(densityOut))
                        throw new NotFiniteNumberException();
                    res = 0.5 * (densityIn * Uin[Component] + densityOut * Uout[Component]) * inp.Normal[Component];
                    break;
                case PhysicsMode.Combustion:
                    DensityArguments_In = Uin.GetSubVector(m_SpatialDimension, NumberOfSpecies+1);
                    DensityArguments_Out = Uout.GetSubVector(m_SpatialDimension, NumberOfSpecies+1);
                    densityIn = (EoS.GetDensity(DensityArguments_In));
                    densityOut = (EoS.GetDensity(DensityArguments_Out));
                    if (double.IsNaN(densityIn) || double.IsInfinity(densityIn) || double.IsNaN(densityOut) || double.IsInfinity(densityOut))
                        throw new NotFiniteNumberException();
                    res = 0.5 * (densityIn * Uin[Component] + densityOut * Uout[Component]) * inp.Normal[Component];
                    break;
                default:
                    throw new ApplicationException("PhysicsMode not implemented");
            }

            return res;
        }

        void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            int D = output.Length;
            Array.Clear(output, 0, D);
            double[] DensityArguments; // Arguments used to calculate the density with the EoS
            switch(Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    output[Component] = U[Component];
                    break;
                case PhysicsMode.MixtureFraction:
                    output[Component] = EoS.getDensityFromZ(U[inp.D]) * U[Component];
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    DensityArguments = U.GetSubVector(m_SpatialDimension, 1);
                    output[Component] = EoS.GetDensity(DensityArguments) * U[Component];
                    break;
                case PhysicsMode.Combustion:
                    DensityArguments = U.GetSubVector(m_SpatialDimension, NumberOfSpecies+1 ); //MassFraction4 does not exist as a variable, because it is just calculated at the end of each iteration
                    output[Component] = EoS.GetDensity(DensityArguments) * U[Component];
                    break;
                default:
                    throw new ApplicationException("PhysicsMode not implemented");
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
            //for(int d = 0; d < D; d++)
                acc += buf[Component] * GradV[Component];
            return -acc;
        }

        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }
        /// <summary>
        /// Multiplier used with every velocity BC
        /// </summary>
        double VelocityMultiplier = 1.0;
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("VelocityMultiplier"))
                VelocityMultiplier = (double)cs.UserDefinedValues["VelocityMultiplier"];
        }

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
                //return TermActivationFlags.AllOn;

            }
        }

        /// <summary>
        /// <see cref="IEdgeForm.InnerEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
                //return TermActivationFlags.AllOn;
            }
        }

        /// <summary>
        /// <see cref="IVolumeForm.VolTerms"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV | TermActivationFlags.GradV;
                //return TermActivationFlags.AllOn;
            }
        }
    }
}
