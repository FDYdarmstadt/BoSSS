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
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid;
using ilPSP;

namespace BoSSS.Foundation.XDG.Quadrature {

    /// <summary>
    /// Implementation of <see cref="LevelSetIntegrator"/> that is able to
    /// integrate an arbitrary scalar field over a level set. The idea behind
    /// this class is write the integral of a scalar integrand g over the 
    /// zero level set (given by \f$ \Phi = 0\f$ ) I as
    /// \f{multline*}
    /// \int \limits_I g ds \\
    /// = \int \limits_I g \vec{n}_I \vec{n}_I ds \\
    /// = \int \limits_I g \frac{\nabla \Phi}{|\nabla \Phi|} \vec{n}_I ds
    /// \f{multline*}
    /// and using Gauss' theorem on the modified integrand
    /// \f$ 
    /// \vec{g} = g \frac{\nabla \Phi}{|\nabla \Phi|}
    /// \f$ 
    /// </summary>
    public class ScalarFieldLevelSetIntegrator : LevelSetIntegrator {

        /// <summary>
        /// The level set to be integrated over
        /// </summary>
        private ILevelSet m_levSet;

        /// <summary>
        /// Buffer for the gradient of the level set, see
        /// <see cref="ILevelSet.EvaluateGradient"/>
        /// </summary>
        private MultidimensionalArray m_LevelSetGradientBuffer = new MultidimensionalArray(3);

        /// <summary>
        /// Buffer for the second derivatives of the levels et,
        /// <see cref="ILevelSet.EvaluateHessian"/>
        /// </summary>
        private MultidimensionalArray m_LevelSetHessianBuffer = new MultidimensionalArray(4);

        /// <summary>
        /// Buffer for the weight function, see
        /// <see cref="LevelSetIntegrator.EvaluateWeightFunction"/>
        /// </summary>
        private MultidimensionalArray m_WeighFunctionBuffer = new MultidimensionalArray(2);

        /// <summary>
        /// Buffer for the values of the field to be integrated, see
        /// <see cref="m_Field"/>
        /// </summary>
        private MultidimensionalArray m_IntegrandBuffer = new MultidimensionalArray(2);

        /// <summary>
        /// Buffer for the gradients of the field to be integrated, see
        /// <see cref="m_Field"/>
        /// </summary>
        private MultidimensionalArray m_IntegrandGradientBuffer = new MultidimensionalArray(3);

        /// <summary>
        /// The field to be integrated
        /// </summary>
        private SinglePhaseField m_Field;

        /// <summary>
        /// Constructs an integrator for the given field
        /// <paramref name="field"/>.
        /// </summary>
        /// <param name="trk">
        /// The tracker containing the level set to be integrated over
        /// </param>
        /// <param name="field">
        /// The field to be integrated
        /// </param>
        /// <param name="volumeFactory">
        /// <see cref="LevelSetIntegrator.LevelSetIntegrator"/>
        /// </param>
        /// <param name="boundaryFactory">
        /// <see cref="LevelSetIntegrator.LevelSetIntegrator"/>
        /// </param>
        /// <param name="quadOrder">
        /// <see cref="LevelSetIntegrator.LevelSetIntegrator"/>
        /// </param>
        /// <param name="sgrd"></param>
        /// <param name="LevSetIdx"></param>
        public ScalarFieldLevelSetIntegrator(
            LevelSetTracker trk,
            SinglePhaseField field,
            ICompositeQuadRule<QuadRule> volumeFactory,
            ICompositeQuadRule<CellBoundaryQuadRule> boundaryFactory,
            int quadOrder,
            SubGrid sgrd,
            int LevSetIdx)
            : base(trk, 1, volumeFactory, boundaryFactory, quadOrder, sgrd, LevSetIdx) {
            m_levSet = trk.LevelSets[0];
            m_Field = field;
        }

        /// <summary>
        /// Evaluates the modified integrand which is
        /// \f$ 
        /// \vec{g} = g \frac{\nabla \Phi}{|\nabla \Phi|}
        /// \f$ 
        /// where g represents <see cref="m_Field"/>.
        /// </summary>
        /// <param name="j0">
        /// <see cref="LevelSetIntegrator.EvaluateIntegrand"/>
        /// </param>
        /// <param name="Length">
        /// <see cref="LevelSetIntegrator.EvaluateIntegrand"/>
        /// </param>
        /// <param name="EvalResult">
        /// <see cref="LevelSetIntegrator.EvaluateIntegrand"/>
        /// </param>
        public override void EvaluateIntegrand(NodeSet N, int j0, int Length, MultidimensionalArray EvalResult) {
            using (new FuncTrace()) {
                int D = base.m_LevSetTrk.GridDat.SpatialDimension; // spatial dimension
                int noOfNodes = EvalResult.GetLength(1);  // number of nodes

                if (m_LevelSetGradientBuffer.GetLength(0) != Length || m_LevelSetGradientBuffer.GetLength(1) != noOfNodes) {
                    m_IntegrandBuffer.Allocate(Length, noOfNodes);
                    m_LevelSetGradientBuffer.Allocate(Length, noOfNodes, D);
                }

                m_Field.Evaluate(j0, Length, N, m_IntegrandBuffer, 0.0);
                m_levSet.EvaluateGradient(j0, Length, N, m_LevelSetGradientBuffer);

                for (int i = 0; i < Length; i++) {
                    for (int j = 0; j < noOfNodes; j++) {
                        double AbsGradPhi = 0;

                        for (int d = 0; d < D; d++) {
                            double GradPhi_d = m_LevelSetGradientBuffer[i, j, d];
                            AbsGradPhi += GradPhi_d * GradPhi_d;
                        }
                        AbsGradPhi = Math.Sqrt(AbsGradPhi);
                        double ooAbsGradPhi = 1.0 / AbsGradPhi;

                        for (int d = 0; d < D; d++) {
                            EvalResult[i, j, 0, d] = ooAbsGradPhi * m_LevelSetGradientBuffer[i, j, d] * m_IntegrandBuffer[i, j];
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Evaluates the divergence of the modified integrand 
        /// \f$ 
        /// \nabla \cdot \vec{g} = \nabla \cdot (g \frac{\nabla \Phi}{|\nabla \Phi|})
        /// \f$ 
        /// which can be expanded to
        /// \f$ 
        /// \nabla \cdot \vec{g} = \nabla g \frac{\nabla \Phi}{|\nabla \Phi|} + g (
        /// \frac{\Delta \Phi}{|\nabla \Phi|}
        /// - \frac{\nabla \Phi}{|\nabla \Phi|} \frac{H(\Phi)}{|\nabla \Phi|} \frac{\nabla \Phi}{|\nabla \Phi|})
        /// \f$ 
        /// using the chain rule. HEre, \f$ H(\Phi)\f$ 
        /// denotes the Hessian of the level set (i.e., the second derivatives)
        /// </summary>
        /// <param name="j0">
        /// <see cref="LevelSetIntegrator.EvaluateDivergenceOfIntegrand"/>
        /// </param>
        /// <param name="Length">
        /// <see cref="LevelSetIntegrator.EvaluateDivergenceOfIntegrand"/>
        /// </param>
        /// <param name="EvalResult">
        /// <see cref="LevelSetIntegrator.EvaluateDivergenceOfIntegrand"/>
        /// </param>
        public override void EvaluateDivergenceOfIntegrand(NodeSet nodes, int j0, int Length, MultidimensionalArray EvalResult) {
            using (new FuncTrace()) {
                int D = base.m_LevSetTrk.GridDat.SpatialDimension; // spatial dimension
                int noOfNodes = EvalResult.GetLength(1);  // number of nodes

                if (m_WeighFunctionBuffer.GetLength(0) != Length || m_WeighFunctionBuffer.GetLength(1) != noOfNodes) {
                    m_IntegrandBuffer.Allocate(Length, noOfNodes);
                    m_IntegrandGradientBuffer.Allocate(Length, noOfNodes, D);
                    m_WeighFunctionBuffer.Allocate(Length, noOfNodes);
                    m_LevelSetGradientBuffer.Allocate(Length, noOfNodes, D);
                    m_LevelSetHessianBuffer.Allocate(Length, noOfNodes, D, D);
                }

                m_Field.Evaluate(j0, Length, nodes, m_IntegrandBuffer);
                m_Field.EvaluateGradient(j0, Length, nodes, m_IntegrandGradientBuffer, 0, 0.0);
                m_levSet.EvaluateGradient(j0, Length, nodes, m_LevelSetGradientBuffer);
                m_levSet.EvaluateHessian(j0, Length, nodes, m_LevelSetHessianBuffer);

                EvaluateWeightFunction(nodes, j0, Length, m_WeighFunctionBuffer);

                double[] Buf = new double[D];

                for (int i = 0; i < Length; i++) {
                    for (int j = 0; j < noOfNodes; j++) {
                        if (m_WeighFunctionBuffer[i, j] == 0.0) {
                            continue;
                        }

                        double AbsGradPhi = 0.0;

                        for (int d = 0; d < D; d++) {
                            double GradPhi_d = m_LevelSetGradientBuffer[i, j, d];
                            AbsGradPhi += GradPhi_d * GradPhi_d;
                        }
                        AbsGradPhi = Math.Sqrt(AbsGradPhi);
                        if (AbsGradPhi < 1e-11) {
                            // Assume zero since gradient is nearly zero
                            continue;
                        }

                        double ooAbsGradPhi = 1.0 / AbsGradPhi;

                        double term1 = 0.0;
                        for (int d = 0; d < D; d++) {
                            term1 += m_IntegrandGradientBuffer[i, j, d] * m_LevelSetGradientBuffer[i, j, d];
                        }
                        term1 *= ooAbsGradPhi;

                        double term2 = 0;
                        for (int d = 0; d < D; d++) {
                            term2 += m_LevelSetHessianBuffer[i, j, d, d];
                        }
                        term2 *= ooAbsGradPhi;

                        double term3 = 0;
                        {
                            for (int d1 = 0; d1 < D; d1++) {
                                double a = 0;
                                for (int d2 = 0; d2 < D; d2++) {
                                    a += m_LevelSetHessianBuffer[i, j, d1, d2] * m_LevelSetGradientBuffer[i, j, d2];
                                }
                                Buf[d1] = a;
                            }

                            for (int d = 0; d < D; d++)
                                term3 += Buf[d] * m_LevelSetGradientBuffer[i, j, d];

                            term3 *= (ooAbsGradPhi * ooAbsGradPhi * ooAbsGradPhi);
                        }

                        EvalResult[i, j, 0] = term1 + m_IntegrandBuffer[i, j] * (term2 - term3);
                    }
                }
            }
        }
    }
}