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
            switch (Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                case PhysicsMode.Multiphase:
                    throw new ApplicationException("Wrong constructor");
                case PhysicsMode.LowMach:
                    this.m_ParamterOrdering = new string[] { VariableNames.Temperature0 };
                    break;
                case PhysicsMode.Combustion:
                    this.m_ParamterOrdering = ArrayTools.Cat(new string[] { VariableNames.Temperature0 }, VariableNames.MassFractions0(NumberOfReactants));
                    this.NumberOfReactants = NumberOfReactants;
                    if (NumberOfReactants == -1)
                        throw new ArgumentException("NumberOfReactants must be specified for combustion flows.");
                    break;
                default:
                    throw new NotImplementedException();
            }
        }

        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            double res = 0.0;

            IncompressibleBcType edgeType = Bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                        res = 0.0;
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        double TemperatureOut = 0.0;
                        double Uout = VelFunction[inp.EdgeTag](inp.X, inp.time);
                        switch (Bcmap.PhysMode) {
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase: {
                                    // opt1:
                                    TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, inp.time);
                                    res = EoS.GetDensity(TemperatureOut) * Uout * inp.Normale[Component];
                                    // opt2:
                                    //double TemperatureIN = inp.Parameters_IN[0];
                                    //double rhoIn = EoS.GetDensity(TemperatureIN);
                                    //res = rhoIn * Uout * inp.Normale[Component];
                                    break;
                                }
                            case PhysicsMode.Combustion:
                                // opt1:
                                TemperatureOut = Bcmap.bndFunction[VariableNames.Temperature][inp.EdgeTag](inp.X, 0);
                                double[] args = new double[NumberOfReactants + 1];
                                args[0] = TemperatureOut;
                                for (int n = 1; n < NumberOfReactants + 1; n++) {
                                    args[n] = Bcmap.bndFunction[VariableNames.MassFraction_n(n - 1)][inp.EdgeTag](inp.X, 0);
                                }
                                res = EoS.GetDensity(args) * Uout * inp.Normale[Component];
                                break;
                            case PhysicsMode.Incompressible:
                                res = Uout * inp.Normale[Component];
                                break;
                            default:
                                throw new ApplicationException("PhysicsMode not implemented");
                        }
                    }
                    break;
                case IncompressibleBcType.Pressure_Outlet: {
                        switch (Bcmap.PhysMode) {
                            case PhysicsMode.Incompressible:
                                res = Uin[0] * inp.Normale[Component];
                                break;
                            case PhysicsMode.LowMach:
                            case PhysicsMode.Multiphase:
                                res = EoS.GetDensity(inp.Parameters_IN[0]) * Uin[0] * inp.Normale[Component];
                                break;
                            case PhysicsMode.Combustion:
                                res = EoS.GetDensity(inp.Parameters_IN) * Uin[0] * inp.Normale[Component];
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

            switch (Bcmap.PhysMode) {
                case PhysicsMode.Incompressible:
                    res = 0.5 * (Uin[0] + Uout[0]) * inp.Normale[Component];
                    break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Multiphase:
                    res = 0.5 * (EoS.GetDensity(inp.Parameters_IN[0]) * Uin[0] + EoS.GetDensity(inp.Parameters_OUT[0]) * Uout[0]) * inp.Normale[Component];
                    break;
                case PhysicsMode.Combustion:
                    res = 0.5 * (EoS.GetDensity(inp.Parameters_IN) * Uin[0] + EoS.GetDensity(inp.Parameters_OUT) * Uout[0]) * inp.Normale[Component];
                    break;
                default:
                    throw new ApplicationException("PhysicsMode not implemented");
            }

            return res;
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            int D = output.Length;
            Array.Clear(output, 0, D);

            switch (Bcmap.PhysMode) {
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
}
