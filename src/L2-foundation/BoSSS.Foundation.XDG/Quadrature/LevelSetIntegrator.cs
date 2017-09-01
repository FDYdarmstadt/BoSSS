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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Integrates some vector field over the manifold defined by a zero
    /// level set. To do so, Gauss' theorem is applied to the surface integral
    /// over the level set. This results into
    /// <list type="bullet">
    ///     <item>
    ///     a surface integral over the parts of the edges of the integration
    ///     cell associated with species A and
    ///     </item>
    ///     <item>
    ///     a volume integral over the sub-volume of the integration cell
    ///     associated with species A
    ///     </item>
    /// </list>
    /// which will be evaluated by this class.
    /// </summary>
    public abstract class LevelSetIntegrator {

        /// <summary>
        /// The level set tracker tracking the level set to be integrated over
        /// </summary>
        protected LevelSetTracker m_LevSetTrk;

        /// <summary>
        /// The number of integrals to perform simultaneously
        /// </summary>
        private int m_NoOfIntegrands;

        /// <summary>
        /// Integrator for the surface integrals
        /// </summary>
        private SurfaceIntegrator m_SurfaceIntegrator;

        /// <summary>
        /// Integrator for the volume integrals
        /// </summary>
        private VolumeIntegrator m_VoumeIntegrator;

        /// <summary>
        /// Index of the level set to be integrated over
        /// </summary>
        private int m_LevSetIdx = -1;

        /// <summary>
        /// Stores the integration result
        /// </summary>
        private MultidimensionalArray m_IntegrationResult;

        /*
        /// <summary>
        /// Cache for the volume fractions occupied by a given species in each
        /// cell. Required by <see cref="EvaluateWeightFunction"/>.
        /// </summary>
        private double[,] m_VolumeFractions;
        */

        /// <summary>
        /// Restriction of the computational domain.
        /// </summary>
        private SubGrid m_SubGrid;

        /// <summary>
        /// Constructs an integrator.
        /// </summary>
        /// <param name="lsTrk">
        /// The level set tracker tracking the level set to be integrated over
        /// </param>
        /// <param name="NoOfIntegrands">
        /// The number of integrals to perform simultaneously
        /// </param>
        /// <param name="volumeFactory">
        /// The quadrature rule factory for the volume integrals. Note that the
        /// integrand is discontinuous and that special quadrature rules should
        /// be used for this case.
        /// </param>
        /// <param name="boundaryFactory">
        /// The quadrature rule factory for the surface integrals. Note that the
        /// integrand is discontinuous and that special quadrature rules should
        /// be used for this case.
        /// </param>
        /// <param name="quadOrder">
        /// The desired quadrature rule (for surface and volume integrals)
        /// </param>
        /// <param name="LevSetIdx">
        /// The id of the level set to be integrated over.
        /// </param>
        /// <param name="sgrd">
        /// the subgrid which is used for the mapping the results
        /// (see <see cref="Execute"/>).
        /// </param>
        public LevelSetIntegrator(
            LevelSetTracker lsTrk, int NoOfIntegrands,
            ICompositeQuadRule<QuadRule> volumeFactory,
            ICompositeQuadRule<CellBoundaryQuadRule> boundaryFactory,
            int quadOrder,
            SubGrid sgrd,
            int LevSetIdx) {
            m_LevSetTrk = lsTrk;
            m_NoOfIntegrands = NoOfIntegrands;
            m_SubGrid = sgrd;
            m_LevSetIdx = LevSetIdx;

            m_SurfaceIntegrator = new SurfaceIntegrator(this, boundaryFactory);
            m_VoumeIntegrator = new VolumeIntegrator(this, volumeFactory);
        }

        /// <summary>
        /// Override this method by evaluating the integrand of the surface
        /// integral on all edges of all <b>cells</b> between
        /// <paramref name="j0"/> and <paramref name="j0"/> +
        /// <paramref name="Length"/>.
        /// </summary>
        /// <param name="j0">
        /// The index of the first cell to be integrated over.
        /// </param>
        /// <param name="Length">
        /// The number of cells to be integrated over.
        /// </param>
        /// <param name="EvalResult">
        /// On exit: Contains the result of the evaluation
        /// <list type="bullet">
        ///     <item>1st index: Cell index - <paramref name="j0"/></item>
        ///     <item>2nd index: Node index</item>
        ///     <item>3rd index: Integrand index</item>
        ///     <item>4th index: Spatial dimension (0, 1 or 2)</item>
        /// </list>
        /// </param>
        /// <param name="QuadNodes">quadrature nodes</param>
        public abstract void EvaluateIntegrand(NodeSet QuadNodes, int j0, int Length, MultidimensionalArray EvalResult);

        /// <summary>
        /// Override this method by evaluating the divergence of the integrand
        /// of the surface integral (i.e., by evaluating the integrand of the
        /// volume integral) in all cells between <paramref name="j0"/> and
        /// <paramref name="j0"/> + <paramref name="Length"/>.
        /// </summary>
        /// <param name="j0">
        /// The index of the first cell to be integrated over.
        /// </param>
        /// <param name="Length">
        /// The number of cells to be integrated over.
        /// </param>
        /// <param name="EvalResult">
        /// On exit: Contains the result of the evaluation
        /// <list type="bullet">
        ///     <item>1st index: Cell index - <paramref name="j0"/></item>
        ///     <item>2nd index: Node index</item>
        ///     <item>3rd index: Integrand index</item>
        /// </list>
        /// </param>
        /// <param name="QuadNodes">quadrature nodes</param>
        public abstract void EvaluateDivergenceOfIntegrand(NodeSet QuadNodes, int j0, int Length, MultidimensionalArray EvalResult);



        /// <summary>
        /// Calculates the integral(s) over the level set. Version with memory
        /// allocation.
        /// </summary>
        /// <returns>
        /// The integral of the integrand <see cref="EvaluateIntegrand"/> over the
        /// level set.
        /// </returns>
        public MultidimensionalArray ExecuteA() {
            MultidimensionalArray ret = MultidimensionalArray.Create(
                m_SubGrid.LocalNoOfCells, m_NoOfIntegrands);
            Execute(ret);
            return ret;
        }

        /// <summary>
        /// Calculates the integral(s) over the level set. Version without
        /// memory allocation.
        /// </summary>
        /// <param name="IntegrationResult">
        /// On Exit: Contains the integral of the integrand
        /// <see cref="EvaluateIntegrand"/> over the level set.
        /// </param>
        public void Execute(MultidimensionalArray IntegrationResult) {
            using (new FuncTrace()) {
                if (IntegrationResult.Dimension != 2)
                    throw new ArgumentException("dimension must be 2.", "IntegrationResult");
                if (IntegrationResult.GetLength(1) != m_NoOfIntegrands)
                    throw new ArgumentException("mismatch in length of second dimension.", "IntegrationResult");
                if (IntegrationResult.GetLength(0) < m_SubGrid.VolumeMask.NoOfItemsLocally)
                    throw new ArgumentException("first dimension too short.", "IntegrationResult");

                m_IntegrationResult = IntegrationResult;

                //m_VolumeFractions = m_LevSetTrk.GetLevSetVolume(m_LevSetIdx);

                m_SurfaceIntegrator.Execute();
                m_VoumeIntegrator.Execute();
                m_LevSetIdx = int.MinValue;
            }
        }

        /// <summary>
        /// This weight function is used to combine the result of the
        /// integration over the volume with negative level set values with
        /// the result of the integration over the volume with positive level
        /// set values. Used by <see cref="SurfaceIntegrator"/> and
        /// <see cref="VolumeIntegrator"/>.
        /// </summary>
        /// <param name="i0">
        /// The first cell to be integrated over.
        /// </param>
        /// <param name="Length">
        /// The number of cells to be integrated over.
        /// </param>
        /// <param name="outp">
        /// On exit: The weight factor for each node.
        /// <list type="bullet">
        ///     <item>1st index: Cell index - <paramref name="i0"/></item>
        ///     <item>2nd index: Node index</item>
        /// </list>
        /// </param>
        /// <remarks>
        /// The current implementation is built such that only the integral
        /// over the larger sub-volume contributes to the integration result.
        /// This avoids bad results in cells with very small sub-volumes (which
        /// only contain very few integration points)
        /// </remarks>
        /// <param name="nodes">
        /// Quadrature nodes.
        /// </param>
        protected void EvaluateWeightFunction(NodeSet nodes, int i0, int Length, MultidimensionalArray outp) {

            throw new NotImplementedException("Killed due to refactoring/cleanup of LevelSetTracker -- is anyone actually using this stuff???");
            
            /*
            MultidimensionalArray LevSetVal = m_LevSetTrk.GetLevSetValues(m_LevSetIdx, nodes, i0, Length);
            int noOfNodes = LevSetVal.GetLength(1);

            // loop over cells
            for (int i = 0; i < Length; i++) {
                int cellIndex = i + i0;
                int subgridIndex = m_SubGrid.LocalCellIndex2SubgridIndex[cellIndex];

                double NegVol = m_VolumeFractions[subgridIndex, 0];
                double PosVol = m_VolumeFractions[subgridIndex, 1];

                // Only consider the bigger sub-volume
                if (PosVol >= NegVol) {
                    for (int j = 0; j < noOfNodes; j++) {
                        if (LevSetVal[i, j] >= 0.0) {
                            outp[i, j] = -1.0;
                        }
                    }
                } else {
                    for (int j = 0; j < noOfNodes; j++) {
                        if (LevSetVal[i, j] < 0.0) {
                            outp[i, j] = 1.0;
                        }
                    }
                }
            }
            */
        }

        /// <summary>
        /// Integrator for the surface integral.
        /// </summary>
        private class SurfaceIntegrator : CellBoundaryQuadrature<CellBoundaryQuadRule> {

            /// <summary>
            /// Creator of this object.
            /// </summary>
            private LevelSetIntegrator m_Owner;

            /// <summary>
            /// Buffer for <see cref="EvaluateWeightFunction"/>
            /// </summary>
            private MultidimensionalArray m_WeightBuffer;

            /// <summary>
            /// Buffer for the integration results.
            /// </summary>
            private MultidimensionalArray m_ResultBuffer;

            /// <summary>
            /// Constructs a surface integrator.
            /// </summary>
            /// <param name="owner">
            /// Creator of this object.
            /// </param>
            /// <param name="rule">
            /// The quadrature rule to be used
            /// </param>
            public SurfaceIntegrator(LevelSetIntegrator owner, ICompositeQuadRule<CellBoundaryQuadRule> rule)
                : base(new int[] { owner.m_NoOfIntegrands }, owner.m_LevSetTrk.GridDat, rule) {
                m_Owner = owner;
                m_WeightBuffer = new MultidimensionalArray(2);
                m_ResultBuffer = new MultidimensionalArray(4);
            }

            /// <summary>
            /// Allocates memory for <see cref="m_WeightBuffer"/> and
            /// <see cref="m_ResultBuffer"/>
            /// </summary>
            /// <param name="NoOfItems">
            /// <see cref="CellBoundaryQuadrature{T}.AllocateBuffers"/>
            /// </param>
            /// <param name="ruleNodes">
            /// <see cref="CellBoundaryQuadrature{T}.AllocateBuffers"/>
            /// </param>
            protected override void AllocateBuffers(int NoOfItems, NodeSet ruleNodes) {
                base.AllocateBuffers(NoOfItems, ruleNodes);
                int NoOfNodes = ruleNodes.GetLength(0);
                m_WeightBuffer.Allocate(NoOfItems, NoOfNodes);
                m_ResultBuffer.Allocate(NoOfItems, NoOfNodes, m_Owner.m_NoOfIntegrands, gridData.SpatialDimension);
            }

            
            /// <summary>
            /// Evaluates the integrand of the surface integral by making use
            /// of <see cref="EvaluateIntegrand"/> and
            /// <see cref="EvaluateWeightFunction"/>.
            /// </summary>
            protected override void Evaluate(int i0, int Length, CellBoundaryQuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;

                m_Owner.EvaluateIntegrand(QuadNodes, i0, Length, m_ResultBuffer);
                m_Owner.EvaluateWeightFunction(QuadNodes, i0, Length, m_WeightBuffer);

                int D = m_Owner.m_LevSetTrk.GridDat.SpatialDimension;
                double[] Normal = new double[D];
                GridData grdDat = m_Owner.m_LevSetTrk.GridDat;

                // loop over cells
                for(int i = 0; i < Length; i++) {
                    int cellIndex = i + i0;

                    int[] cells2edges = grdDat.Cells.Cells2Edges[cellIndex];

                    // loop over the edges of a cell
                    int nodeIndex = 0;
                    for(int face = 0; face < CurrentRule.NumbersOfNodesPerFace.Length; face++) {
                        // Try to determine face
                        int edge = -1;
                        double sign = 0.0;
                        for(int e = 0; e < cells2edges.Length; e++) {
                            int trialEdge = cells2edges[e];
                            sign = Math.Sign(trialEdge);
                            trialEdge = Math.Abs(trialEdge) - 1;

                            if(grdDat.Edges.FaceIndices[trialEdge, sign > 0 ? 0 : 1] == face) {
                                edge = trialEdge;
                                break;
                            }
                        }

                        if(edge < 0) {
                            throw new Exception();
                        }

                        // load normal;
                        for(int d = 0; d < D; d++) {
                            Normal[d] = grdDat.Edges.NormalsForAffine[edge, d] * sign;
                        }

                        // loop over the nodes of an edge
                        for(int j = 0; j < CurrentRule.NumbersOfNodesPerFace[face]; j++) {
                            // loop over integrands...
                            for(int k = 0; k < m_Owner.m_NoOfIntegrands; k++) {
                                double Acc = 0;
                                for(int d = 0; d < D; d++) {
                                    Acc += Normal[d] * m_ResultBuffer[i, nodeIndex, k, d];
                                }

                                EvalResult[i, nodeIndex, k] = Acc * m_WeightBuffer[i, nodeIndex];
                            }

                            nodeIndex++;
                        }
                    }
                }
            }

            /// <summary>
            /// Stores the integration results in
            /// <see cref="m_IntegrationResult"/> using the order defined by
            /// <see cref="m_SubGrid"/>. Note that the result of the surface
            /// integration is <b>subtracted</b>.
            /// </summary>
            /// <param name="i0">
            /// </param>
            /// <param name="Length">
            /// </param>
            /// <param name="ResultsOfIntegration">
            /// </param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int j = 0; j < Length; j++) {
                    int j_subgrd = m_Owner.m_SubGrid.LocalCellIndex2SubgridIndex[j + i0];
                    for (int e = 0; e < GridDat.iGeomCells.RefElements[CurrentRuleRefElementIndex].NoOfFaces; e++) {
                        for (int k = 0; k < m_Owner.m_NoOfIntegrands; k++) {
                            m_Owner.m_IntegrationResult[j_subgrd, k] -= ResultsOfIntegration[j, e, k];
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Integrator for the volume integral.
        /// </summary>
        private class VolumeIntegrator : Foundation.Quadrature.CellQuadrature {

            /// <summary>
            /// Creator of this object.
            /// </summary>
            private LevelSetIntegrator m_Owner;

            /// <summary>
            /// Buffer for <see cref="EvaluateWeightFunction"/>
            /// </summary>
            private MultidimensionalArray m_WeightBuffer;

            /// <summary>
            /// Constructs a volume integrator.
            /// </summary>
            /// <param name="owner">
            /// Creator of this object.
            /// </param>
            /// <param name="volumeRule">
            /// quadrature rules and domain
            /// </param>
            public VolumeIntegrator(LevelSetIntegrator owner, ICompositeQuadRule<QuadRule> volumeRule)
                : base(new int[] { owner.m_NoOfIntegrands }, owner.m_LevSetTrk.GridDat, volumeRule) {
                m_Owner = owner;
                m_WeightBuffer = new MultidimensionalArray(2);
            }

            

            /// <summary>
            /// Allocates memory for <see cref="m_WeightBuffer"/>.
            /// </summary>
            protected override void AllocateBuffers(int NoOfItems, NodeSet ruleNodes) {
                base.AllocateBuffers(NoOfItems, ruleNodes);
                m_WeightBuffer.Allocate(NoOfItems, ruleNodes.GetLength(0));
            }

            /// <summary>
            /// Modulates the result of
            /// <see cref="EvaluateDivergenceOfIntegrand"/> by the weights
            /// given by <see cref="EvaluateWeightFunction"/>.
            /// </summary>
            protected override void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;
                int NoOfNodes = QuadNodes.NoOfNodes;
                m_Owner.EvaluateDivergenceOfIntegrand(QuadNodes, i0, Length, EvalResult);
                m_Owner.EvaluateWeightFunction(QuadNodes, i0, Length, m_WeightBuffer);

                for(int i = 0; i < Length; i++) {
                    for(int j = 0; j < NoOfNodes; j++) {
                        double gamma = m_WeightBuffer[i, j];
                        for(int k = 0; k < m_Owner.m_NoOfIntegrands; k++) {
                            EvalResult[i, j, k] *= gamma;
                        }
                    }
                }
            }

            /// <summary>
            /// Stores the integration results in
            /// <see cref="m_IntegrationResult"/> using the order defined by
            /// <see cref="m_SubGrid"/>.
            /// </summary>
            /// <param name="i0">
            /// </param>
            /// <param name="Length">
            /// </param>
            /// <param name="ResultsOfIntegration">
            /// </param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    int j_subgrd = m_Owner.m_SubGrid.LocalCellIndex2SubgridIndex[i + i0];
                    for (int k = 0; k < m_Owner.m_NoOfIntegrands; k++) {
                        m_Owner.m_IntegrationResult[j_subgrd, k] += ResultsOfIntegration[i, k];
                    }
                }
            }
        }
    }
}