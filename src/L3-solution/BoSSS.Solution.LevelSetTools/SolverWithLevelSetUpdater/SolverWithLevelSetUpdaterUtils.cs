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
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using BoSSS.Platform;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.LevelSetTools;
using System.Collections;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    public static class LevelSetUtils {


        #region interface related properties (area, length, material points, velocities)

        public static double GetSpeciesArea(LevelSetTracker LsTrk, SpeciesId spcId, int quadRuleOrder = -1) {

            double spcArea = 0.0;

            int order = (quadRuleOrder < 0) ? 1 : quadRuleOrder;
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            CellQuadratureScheme vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        spcArea += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return spcArea;
        }


        public static double GetInterfaceLength(int iLevSet, LevelSetTracker LsTrk, int quadRuleOrder = -1) {

            double interLength = 0.0;

            int order  = (quadRuleOrder < 0) ? 1 : quadRuleOrder;
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        interLength += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return interLength;
        }

        public static MultidimensionalArray GetInterfacePoints(LevelSetTracker LsTrk, SinglePhaseField LevSet, SubGrid sgrd = null, int quadRuleOrderForNodeSet = -1) {

            int D = LsTrk.GridDat.SpatialDimension;
            int p = LevSet.Basis.Degree;
            if (sgrd == null)
                sgrd = LsTrk.Regions.GetCutCellSubgrid4LevSet(0);

            int quadRule = (quadRuleOrderForNodeSet < 0) ? p * 2 : quadRuleOrderForNodeSet;
            NodeSet[] Nodes = LsTrk.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(quadRule).Nodes).ToArray();
            int Jsub = sgrd.LocalNoOfCells;
            int K = Nodes.Max(nds => nds.NoOfNodes);
            int numP = Jsub * K;

            var cp = new BoSSS.Solution.LevelSetTools.ClosestPointFinder(LsTrk, 0, sgrd, Nodes);

            MultidimensionalArray ClosestPoints = cp.X0_global_Resorted;

            MultidimensionalArray interfaceP = new MultidimensionalArray(2);
            interfaceP.Allocate(numP, D);

            for (int d = 0; d < D; d++) {
                MultidimensionalArray cp_d = ClosestPoints.ExtractSubArrayShallow(-1, -1, d).CloneAs().ResizeShallow(numP);
                interfaceP.ExtractSubArrayShallow(-1, d).Acc(1.0, cp_d);
            }

            // sort the interface points (non-equidistant)
            int[] permutation = new int[numP];
            for (int ind = 0; ind < numP; ind++)
                permutation[ind] = ind;

            Array.Sort(permutation, (int i, int j) => Math.Sign(interfaceP[i, 0] - interfaceP[j, 0]));

            MultidimensionalArray interfaceP_temp = interfaceP.CloneAs();
            for (int ip = 0; ip < numP; ip++) {
                for (int d = 0; d < D; d++) {
                    interfaceP[ip, d] = interfaceP_temp[permutation[ip], d];
                }
            }

            return interfaceP;
        }

        #endregion      

    }
}
