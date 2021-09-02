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
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Derivative in <em>d</em>-th coordinate direction, using central difference;
    /// </summary>
    /// <remarks>
    /// This is a linear flux, it can be used to assemble the matrix of the spatial operator.
    /// A mathematical equal, nonlinear version,
    /// which can be used for direct evaluation by quadrature,
    /// is given by <see cref="PressureGradientNonlin_d"/>.
    /// </remarks>
    public class PressureGradientLin_d : LinearFlux {

        protected int m_d = -1;

        BoundaryCondMap<IncompressibleBcType> m_bcmap;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_d">
        /// spatial direction of derivative
        /// </param>
        /// <param name="bcmap"></param>
        public PressureGradientLin_d(int _d, IncompressibleBoundaryCondMap bcmap){
            m_d = _d;
            m_bcmap = bcmap;
            pressureFunction = bcmap.bndFunction[VariableNames.Pressure];
        }

        /// <summary>
        /// a mapping from edge tags to Dirichlet boundary values.
        /// </summary>
        protected Func<double[], double, double>[] pressureFunction;

        /// <summary>
        /// A central difference at Dirichlet boundary regions (<see cref="IncompressibleBcType.Pressure_Outlet"/>),
        /// an open end everywhere else.
        /// </summary>
        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            IncompressibleBcType edgType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet:
                    // Atmospheric outlet/pressure outlet: inhom. Dirichlet
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
                    return pressureFunction[inp.EdgeTag](inp.X, inp.time) * inp.Normal[m_d];

                case IncompressibleBcType.Outflow:
                    throw new ArithmeticException("Tests on channel flow indicate that b.c. " + edgType + " is ill-posed, fk 25may16.");
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.NoSlipNeumann:
                    // hom. Neumann b.c.
                    // +++++++++++++++++
                    return Uin[0] * inp.Normal[m_d];
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// central difference Riemannian
        /// </summary>
        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return 0.5 * (Uin[0] + Uout[0]) * inp.Normal[m_d];
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            int D = output.Length;
            Array.Clear(output, 0, D);
            output[m_d] = U[0];
        }

       

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Pressure };
            }
        }
    }

    /// <summary>
    /// An optional (constant) source term for the <em>d</em>-th component of the pressure gradient <see cref="PressureGradientLin_d"/>.
    /// </summary>
    /// <remarks>
    /// Can be used to implement periodic boundary conditions.
    /// </remarks>
    public class SrcPressureGradientLin_d : LinearSource {

        double m_SrcGradPressure_d;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_SrcGradPressure_d">
        /// Value of the source term for the <em>d</em>-th component of the pressure gradient.
        /// </param>
        public SrcPressureGradientLin_d(double _SrcGradPressure_d) {
            m_SrcGradPressure_d = _SrcGradPressure_d;
        }

        /// <summary>
        /// the variable name for the pressure
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Pressure };
            }
        }

        /// <summary>
        /// returns the value of the constant pressure gradient (<em>d</em>-th component),
        /// which is provided by the constructor.
        /// </summary>
        protected override double Source(double[] x, double[] parameters, double[] U) {
            return m_SrcGradPressure_d;
        }
    }


}
