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
using System.Text;
using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LevelSetTools;
using System.Linq;
using ilPSP.Tracing;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// checks the Gauss theorem on cut cells and refines if a given error threshold is exceeded
    /// </summary>
    [Serializable]
    public class AMRbasedOnLocalCurvature : LevelSetTools.SolverWithLevelSetUpdater.AMRLevelIndicatorWithLevelset {

        /// <summary>
        /// Threshold for triggering local mesh refinement:
        /// If the local radius of the interface is smaller than the local mesh width by this fraction, mesh refinement is triggered.
        /// </summary>
        public double LocalInterfaceRadius_To_MeshWidth_Ratio = 0.5;


        public int LevelSetIndex = 0;


        public override int[] DesiredCellChanges() {
            using(var tr = new FuncTrace("AMRcheckGauss")) {
                tr.InfoToConsole = true;

                int J = GridData.CellPartitioning.LocalLength;
                double ooD = 1.0 / GridData.SpatialDimension;
                
                double GetLogicalCellLenthSchale(int j) {
                    double vol = base.GridData.iLogicalCells.GetCellVolume(j);
                    return Math.Pow(vol, ooD);
                } 
                var cutCells = base.LsTrk.Regions.GetCutCellMask4LevSet(this.LevelSetIndex);


                // obtain local curvature minimum for each **logical** cut cell
                // ============================================================
                MultidimensionalArray localCurv = MultidimensionalArray.Create(J);
                {
                    int order = (LsTrk.LevelSets[this.LevelSetIndex] as LevelSet)?.Basis?.Degree ?? 4;
                    order *= 2;

                    var testNodes = GridData.iGeomCells.RefElements.Select(
                        Kref => Kref.GetQuadratureRule(order).Nodes
                    ).ToArray();


                    var geom2Log = GridData.iGeomCells.GeomCell2LogicalCell;
                    foreach(int jGeom in cutCells.ToGeometicalMask().ItemEnum) {
                        int jLog = geom2Log != null ? geom2Log[jGeom] : jGeom;
                        int iKref = GridData.iGeomCells.GetRefElementIndex(jGeom);


                        var curvValues = MultidimensionalArray.Create(1, testNodes[iKref].NoOfNodes);
                        base.LevSet.EvaluateTotalCurvature(J, 1, testNodes[iKref], curvValues);


                        localCurv[jLog] = Math.Min(localCurv[jLog], curvValues.Min());

                    }
                }

                int[] levels = new int[J];
                foreach(int jLog in cutCells.ItemEnum) {
                    throw new Exception("todo");
                }

                return levels;
            }
        }
    }



}