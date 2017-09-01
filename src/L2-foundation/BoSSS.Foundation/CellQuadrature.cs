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
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// quadrature over cells.
    /// </summary>
    public abstract class CellQuadrature : Quadrature<QuadRule, CellMask> {

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="noOfIntegralsPerCell"></param>
        /// <param name="context"></param>
        /// <param name="rule"></param>
        /// <param name="cs">integrate in physical or reference coordinates?</param>
        public CellQuadrature(int[] noOfIntegralsPerCell, IGridData context, ICompositeQuadRule<QuadRule> rule, CoordinateSystem cs = Quadrature.CoordinateSystem.Physical)
            : base(noOfIntegralsPerCell, context, rule, cs) //
        {
            foreach(IChunkRulePair<QuadRule> crp in rule) {
                NodeCoordinateSystem ncs = crp.Rule.Nodes.GetNodeCoordinateSystem(context);
                if(ncs!= NodeCoordinateSystem.CellCoord) {
                    throw new ArgumentException("Illegal node set for cell volume quadrature. Found some node set defined for: " + ncs.ToString() + ".");
                }
            }
        }

        
        /// <summary>
        /// the absolute value of the Jacobian determinate, for the 
        /// transformation that translates cell reference coordinates into physical coordinates.
        /// </summary>
        protected override MultidimensionalArray GetScalingsForLinearElements(int i0, int L) {
            return gridData.iGeomCells.JacobiDet.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + L - 1 });
        }

        /// <summary>
        /// the absolute value of the Jacobian determinate, for the 
        /// transformation that translates cell reference coordinates into physical coordinates.
        /// </summary>
        protected override MultidimensionalArray GetScalingsForNonlinElements(int jCell0, int L) {
            return base.GridDat.JacobianDeterminat.GetValue_Cell(this.CurrentRule.Nodes, jCell0, L);
        }

        /// <summary>
        /// Sweeps whether cell <paramref name="i0"/> is linear/nonlinear and how many cells of the same type are going to come after it.
        /// </summary>
        protected override void NextPart(out bool Linear, out int NoOfElm, int i0, int Len) {
            var cells = base.gridData.iGeomCells;

            bool bLinear = cells.IsCellAffineLinear(i0);
            Linear = bLinear;
            int L;
            for (L = 1; L < Len; L++) {
                if (cells.IsCellAffineLinear(i0 + L) != bLinear) {
                    NoOfElm = L;
                    return;
                }
            }

            NoOfElm = Len;
        }


        /// <summary>
        /// creates a cell quadrature, where integrand evaluation (<paramref name="_Evaluate"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="CellQuadrature"/>.
        /// </summary>
        static public CellQuadrature GetQuadrature(int[] noOfIntegralsPerCell,
                              IGridData context,
                              ICompositeQuadRule<QuadRule> domNrule,
                              Del_Evaluate _Evaluate,
                              Del_SaveIntegrationResults _SaveIntegrationResults,
                              Del_AllocateBuffers _AllocateBuffers = null,
                              Del_QuadNodesChanged _PostLockNodes = null,
                              CoordinateSystem cs = CoordinateSystem.Physical) {
            var ret = new CellQuadratureImpl(noOfIntegralsPerCell, context, domNrule, cs) {
                m_Evaluate = _Evaluate,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                m_AllocateBuffers = _AllocateBuffers,
                m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }

        /// <summary>
        /// Creates a cell quadrature, where integrand evaluation (<paramref name="_Evaluate"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="CellQuadrature"/>.
        /// </summary>
        static public CellQuadrature GetQuadrature2(int[] noOfIntegralsPerCell,
                              IGridData context,
                              ICompositeQuadRule<QuadRule> domNrule,
                              Del_Evaluate _Evaluate,
                              Del_SaveIntegrationResults _SaveIntegrationResults,
                              Del_AllocateBuffers _AllocateBuffers = null,
                              Del_QuadNodesChanged _PostLockNodes = null,
                              CoordinateSystem cs = CoordinateSystem.Physical) {
            var ret = new CellQuadratureImpl(noOfIntegralsPerCell, context, domNrule, cs) {
                m_ExEvaluate = _Evaluate,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                m_AllocateBuffers = _AllocateBuffers,
                m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }


        private class CellQuadratureImpl : CellQuadrature {

            public CellQuadratureImpl(
                int[] noOfIntegralsPerCell, IGridData context, ICompositeQuadRule<QuadRule> domNrule, CoordinateSystem cs)
                : base(noOfIntegralsPerCell, context, domNrule, cs) {
            }

        }

        /// <summary>
        /// Some DEBUG checks on the current quadrature chunk.
        /// </summary>
        protected override void CheckQuadratureChunk(int j0, int Len, int iKRef) {
            // check1: Ref element index:
            int JE = Len + j0;
            for (int j = 0; j < JE; j++) {
                int _iKref = gridData.iGeomCells.GetRefElementIndex(j);
                if (iKRef != _iKref)
                    throw new Exception("Internal error: mismatch between reference element index for given quadrature rule and specific cell that should be integrated.");

            }

        }

        /// <summary>
        /// Index of the current rule's reference element into
        /// <see cref="BoSSS.Foundation.Grid.GridCommons.RefElements"/>.
        /// </summary>
        public override int CurrentRuleRefElementIndex {
            get {
                var Kref = base.CurrentRule.RefElement;
                int ret = Array.IndexOf(gridData.iGeomCells.RefElements, Kref);
                Debug.Assert(ret >= 0, "unable to identify reference element");
                return ret;
            }
        }
    }

}
