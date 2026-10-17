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
using BoSSS.Platform.Utils.Geom;
using ilPSP.Utils;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// checks the average curvature radius in cut cells and refines if a given error threshold is exceeded
    /// </summary>
    [Serializable]
    public class AMRbasedOnLocalCurvature : LevelSetTools.SolverWithLevelSetUpdater.AMRLevelIndicatorWithLevelset {

        /// <summary>
        /// Threshold for triggering local mesh refinement:
        /// If the local radius of the interface is smaller than the local mesh width by this fraction, mesh refinement is triggered.
        /// - a higher value (e.g., 2) in general triggers more aggressive refinement,
        /// - a lower value (e.g., 0.1) is less aggressive and allows higher curvature in cells without triggering refinement.
        /// </summary>
        public double LocalInterfaceRadius_To_MeshWidth_Ratio = 0.5;


        public int LevelSetIndex = 0;


        public override int[] DesiredCellChanges() {
            using(var tr = new FuncTrace("AMRcheckGauss")) {
                tr.InfoToConsole = true;

                int J = GridData.CellPartitioning.LocalLength;
                double ooD = 1.0 / GridData.SpatialDimension;
                
                /*
                double GetLogicalCellLenthSchale(int j) {
                    double vol = base.GridData.iLogicalCells.GetCellVolume(j);
                    return Math.Pow(vol, ooD);
                }
                */
                var cutCells = base.LsTrk.Regions.GetCutCellMask4LevSet(this.LevelSetIndex);
                

                // obtain local curvature maximum for each **logical** cut cell
                // ============================================================
                MultidimensionalArray MaxCurv_LogicalCell = MultidimensionalArray.Create(J);
                {
                    int order = (LsTrk.LevelSets[this.LevelSetIndex] as LevelSet)?.Basis?.Degree ?? 4;
                    
                    var testNodes = GridData.iGeomCells.RefElements.Select(
                        Kref => Kref.GetQuadratureRule(order).Nodes
                    ).ToArray();


                    var geom2Log = GridData.iGeomCells.GeomCell2LogicalCell;
                    foreach(int jGeom in cutCells.ToGeometicalMask().ItemEnum) {
                        int jLog = geom2Log != null ? geom2Log[jGeom] : jGeom;
                        int iKref = GridData.iGeomCells.GetRefElementIndex(jGeom);


                        var curvValues = MultidimensionalArray.Create(1, testNodes[iKref].NoOfNodes);
                        base.LevSet.EvaluateTotalCurvature(jGeom, 1, testNodes[iKref], curvValues);

                        double maxCurv = Math.Max(curvValues.Min().Abs(), curvValues.Max().Abs()); // maximum curvature in geometrical cell
                        MaxCurv_LogicalCell[jLog] = Math.Max(MaxCurv_LogicalCell[jLog], maxCurv); // maximum curvature in logical call
                    }
                }

                int[] levels = new int[J];
                var bb = new BoundingBox(GridData.SpatialDimension);
                int cellsToCoarse = 0, cellsToRefine = 0;
                levels.SetAll(-1);
                foreach(int jLog in cutCells.ItemEnum) {
                    double maxCurv = MaxCurv_LogicalCell[jLog];
                    double minRadius = Math.Min(1.0 / maxCurv, double.MaxValue);
                    if(minRadius.IsNaNorInf())
                        throw new ArithmeticException($"NAN/INF in local radius evaluation: minRadius = {minRadius}, (jLog = {jLog}, maxCurv = {maxCurv})");

                    bb.Clear();
                    GridData.iLogicalCells.GetCellBoundingBox(jLog, bb);
                    double h_local = bb.h_min;

                    double ratio = minRadius / h_local;


                    if(ratio < LocalInterfaceRadius_To_MeshWidth_Ratio) {
                        levels[jLog] = 1;
                        cellsToRefine++;
                    } else {
                        levels[jLog] = -1;
                        cellsToCoarse++;
                    }
                }

                tr.Info($"curvature AMR criterion: cells to refine: {cellsToRefine}, cells to coarsen: {cellsToCoarse}");

                return levels;
            }
        }
    }



}