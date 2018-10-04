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
using System.Diagnostics;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Volume part of the divergence according to:
    /// K. Shahbazi, P. F. Fischer, C. R. Ethier,
    /// A high-order discontinuous Galerkin method for the unsteady incompressible Navier-Stokes equations,
    /// J. Comput. Phys. 222 (2007) 391–407.
    /// (see BoSSS\doc\notes\0010-BoSSS-References\ShahbaziEtAl2007)
    /// </summary>
    /// <remarks>
    /// Note: Here we use the positive divergence in contrast to the cited paper.
    /// The surface part of the divergence is implemented in <see cref="Divergence_DerivativeSource_Flux"/>.
    /// </remarks>
    public class Divergence_DerivativeSource : LinearDerivativeSource {

        int component;
        int D;

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>
        /// <param name="component">
        /// component of the divergence
        /// </param>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public Divergence_DerivativeSource(int component, int D) {
            this.component = component;
            this.D = D;
        }

        /// <summary>
        /// bla bla bla.
        /// </summary>
        public override double _DerivativeSource(double[] x, double[] Parameters, double[,] GradientU) {
            return GradientU[0, this.component];
        }


        /// <summary>
        /// name of the d-th velocity component
        /// </summary>
        override public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(component) };
            }
        }



    }

    /// <summary>
    /// Surface part of the divergence according to:
    /// K. Shahbazi, P. F. Fischer, C. R. Ethier,
    /// A high-order discontinuous Galerkin method for the unsteady incompressible Navier-Stokes equations,
    /// J. Comput. Phys. 222 (2007) 391–407.
    /// (see BoSSS\doc\notes\0010-BoSSS-References\ShahbaziEtAl2007)
    /// </summary>
    /// <remarks>
    /// Note: Here we use the positive divergence in contrast to the cited paper.
    /// The volume part of the divergence is implemented in <see cref="Divergence_DerivativeSource"/>.
    /// </remarks>
    public class Divergence_DerivativeSource_Flux : LinearDualValueFlux {

        /// <summary>
        /// spatial direction of derivative.
        /// </summary>
        protected int component;


        IncompressibleBoundaryCondMap bcmap;

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>
        /// <param name="component">component of the divergence</param>
        /// <param name="bcmap"></param>        
        public Divergence_DerivativeSource_Flux(int component, IncompressibleBoundaryCondMap bcmap){//, double PresPenalty2) {
            this.component = component;
            this.bcmap = bcmap;
            this.bndFunction = bcmap.bndFunction[VariableNames.Velocity_d(component)];
        }

        /// <summary>
        /// Mapping from edge tags to boundary values.
        /// </summary>
        protected Func<double[], double, double>[] bndFunction;

        /// <summary>
        /// name of the d-th velocity component
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(component)};
            }
        }

        /// <summary>
        /// returns a penaltization of the velocity jump, i.e.
        /// \f[ 
        ///   -\frac{1}{2} \llbracket u_d \rrbracket \cdot \vec{n} \cdot \vec{e}_d,
        /// \f]
        /// where d is the velocity component index provided in the constructor
        /// </summary>
        protected override void InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOutCell) {
            double u_j_In = Uin[0];
            double u_j_Out = Uout[0];

            FluxInCell = -0.5 * (u_j_In - u_j_Out) * inp.Normale[component];
            FluxOutCell = FluxInCell;
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            IncompressibleBcType edgeType = bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Outflow:
                throw new ArithmeticException("Tests on channel flow indicate that b.c. " + edgeType + " is ill-posed, fk 25may16.");
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Pressure_Outlet:
                {
                    FluxInCell = 0.0;
                    break;
                }
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.NavierSlip_Linear: {
                    double u_j_In = Uin[0];
                    double u_j_Out = this.bndFunction[inp.EdgeTag](inp.X, inp.time);

                    FluxInCell = -(u_j_In - u_j_Out) * inp.Normale[component];
                    break;
                }
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NoSlipNeumann:
                {
                    double u_j_In = Uin[0];
                    FluxInCell = -u_j_In * inp.Normale[component];
                    break;
                }
                default:
                throw new NotImplementedException("Boundary condition not implemented!");
            }
        }
    }
}