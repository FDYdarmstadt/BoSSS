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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    ///
    /// </summary>
    [Serializable]
    public class AMR_onProblematicPoints : AMRLevelIndicatorWithLevelset {

        [DataMember]
        public IList<double[]> problematicPoints;

        private AMR_onProblematicPoints() {
        }

        public AMR_onProblematicPoints(IList<double[]> points, int maxRefinementLevelval) {
            if(points[0].Length != 2) {
                throw new Exception("only supported for 2d");
            }

            problematicPoints = points;
            maxRefinementLevel = maxRefinementLevelval;
        }

        public override int[] DesiredCellChanges() {
            double eps = 1e-8;
            int J = GridData.CellPartitioning.LocalLength;
            int D = GridData.SpatialDimension;
            int[] levels = new int[J];
            Cell[] cells = GridData.Grid.Cells;

            var allPoints = new List<double[]>() { };
            foreach(var point in problematicPoints) {
                allPoints.Add(new double[] { point[0] + eps, point[1] + eps });
                allPoints.Add(new double[] { point[0] - eps, point[1] + eps });
                allPoints.Add(new double[] { point[0] + eps, point[1] - eps });
                allPoints.Add(new double[] { point[0] - eps, point[1] - eps });
            }

            for(int p = 0; p < allPoints.Count; p++) {
                long GlobalID, GlobalIndex; // Index of the cell to be refined
                bool IsInside, OnThisProcess;
                double[] point = allPoints[p];
                GridData.LocatePoint(point, out GlobalID, out GlobalIndex, out IsInside, out OnThisProcess);

                // Get local cell index of current point
                long j0Grd = GridData.CellPartitioning.i0;
                int jLocal = (int)(GlobalIndex - j0Grd);

                //if(!IsInside) {
                //    Console.WriteLine("point" + point + "is problematic");
                //}
                int currentRefinement = cells[(int)jLocal].RefinementLevel;
                if(OnThisProcess && currentRefinement < maxRefinementLevel) {                    
                    levels[(int)jLocal] = 1;
                }
            }

            return levels;
        }
    }

    /// <summary>
    /// Refinement around the centerline of the counter-flow configuration
    /// </summary>
    [Serializable]
    public class AMR_AroundCenterline : AMRLevelIndicatorWithLevelset {

    
        private AMR_AroundCenterline() {
        }

        public AMR_AroundCenterline(int maxRefinementLevelval) {
            maxRefinementLevel = maxRefinementLevelval;
        }

        public override int[] DesiredCellChanges() {
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            int noOfPoints = 400;
            var xnodes = GenericBlas.Linspace(0, 1.0, noOfPoints);
            double eps = 1e-3;
            var allPoints = new List<double[]>() { };
  
            for(int i= 0; i < noOfPoints; i++) {
                allPoints.Add(new double[] { xnodes[i], +eps });
                allPoints.Add(new double[] { xnodes[i], -eps });
            }
            Cell[] cells = GridData.Grid.Cells;

            for(int p = 0; p < allPoints.Count; p++) {
                long GlobalID, GlobalIndex; // Index of the cell to be refined
                bool IsInside, OnThisProcess;
                double[] point = allPoints[p];
                GridData.LocatePoint(point, out GlobalID, out GlobalIndex, out IsInside, out OnThisProcess);

                // Get local cell index of current point
                long j0Grd = GridData.CellPartitioning.i0;
                int jLocal = (int)(GlobalIndex - j0Grd);

       
                //int currentLevel = cells[(int)jLocal].RefinementLevel;
                if(OnThisProcess/* && currentLevel<maxRefinementLevel*/ ) {
                    levels[(int)jLocal] = 1;
                }
            }

            return levels;
        }
    }


    /// <summary>
    ///
    ///
    /// </summary>
    [Serializable]
    public class AMR_RefineAroundProblematicPoints : AMRLevelIndicatorWithLevelset {

        [DataMember]
        public IList<double[]> problematicPoints;
        [DataMember]
        public double radius;

        private AMR_RefineAroundProblematicPoints() {
        }

        public AMR_RefineAroundProblematicPoints(IList<double[]> points, int maxRefinementLevelval, double _radius) {
            if(points[0].Length != 2) {
                throw new Exception("only supported for 2d");
            }

            problematicPoints = points;
            maxRefinementLevel = maxRefinementLevelval;
            radius = _radius;
        }

        public override int[] DesiredCellChanges() {
            double eps = 1e-8;
            int J = GridData.CellPartitioning.LocalLength;
            int D = GridData.SpatialDimension;
            int[] levels = new int[J];
            Cell[] cells = GridData.Grid.Cells;
            //System.Diagnostics. dbg_launch();

            for (int p = 0; p < problematicPoints.Count; p++) {
                long GlobalID, GlobalIndex; // Index of the cell to be refined
                bool IsInside, OnThisProcess;
                GridData.LocatePoint(problematicPoints[p], out GlobalID, out GlobalIndex, out IsInside, out OnThisProcess);

                // Get local cell index of current point
                long j0Grd = GridData.CellPartitioning.i0;
                int jLocal = (int)(GlobalIndex - j0Grd);
                if(OnThisProcess) {
                    try {
                        var CellsInRadius = CellsInRadiusAroundCell(jLocal, radius);

                        for (int j = 0; j < CellsInRadius.Count; j++) {

                            int jjLocal = (CellsInRadius[j] - (int)j0Grd * 0);
                            //if(jjLocal >= 0)
                            levels[jjLocal] = 1;
                        }
                    } catch (Exception) {
                        Console.WriteLine("Problem with point (" + problematicPoints[p][0] + "," + problematicPoints[p][1] + "). jLocal is" + jLocal);
                        
                        continue;
                        throw new Exception("Problem with point (" + problematicPoints[p][0] +","+ problematicPoints[p][1]+"). jLocal is" + jLocal);
                    }                  

                }


            }

            return levels;
        }
        public List<int> CellsInRadiusAroundCell(int BasicCellIndex, double Radius) {
            int SpatialDimension = 2;
            BoundingBox Box = new BoundingBox(SpatialDimension);
            BoundingBox Box_Cell = new BoundingBox(SpatialDimension);
            GridData.iGeomCells.GetCellBoundingBox(BasicCellIndex, Box);
            double[] Max_0 = Box.Max;
            double[] Min_0 = Box.Min;
            double[] Max_Cell, Min_Cell;
            List<int> CellIndexDoneList = new List<int>();
            List<int> CellIndexFrontList = new List<int> { BasicCellIndex };

            bool RadiusReached = false;
            while(!RadiusReached) {
                bool AllOutsideRadius = true;
                foreach(int FrontCell in CellIndexFrontList) {
                    GridData.iGeomCells.GetCellBoundingBox(FrontCell, Box_Cell);
                    Max_Cell = Box_Cell.Max;
                    Min_Cell = Box_Cell.Min;
                    bool OneDimOutsideRadius = false;
                    for(int dim = 0; dim < SpatialDimension; dim++) {
                        if(Math.Abs(Min_Cell[dim] - Max_0[dim]) > Radius && Math.Abs(Min_0[dim] - Max_Cell[dim]) > Radius) {
                            OneDimOutsideRadius = true;
                        }
                    }
                    if(FrontCell == BasicCellIndex)
                        OneDimOutsideRadius = false;
                    if(!OneDimOutsideRadius)
                        //Console.WriteLine(FrontCell + " is inside Radius");
                        AllOutsideRadius = (OneDimOutsideRadius) ? AllOutsideRadius : false;
                }

                RadiusReached = (AllOutsideRadius) ? true : false;

                if(!RadiusReached) {
                    foreach(int FrontCell in CellIndexFrontList.ToList()) {
                        GridData.GetCellNeighbours(FrontCell, GetCellNeighbours_Mode.ViaVertices, out int[] Neighbours, out int[] ConnectingEntities);
                        CellIndexDoneList.Add(FrontCell);
                        CellIndexFrontList.Remove(FrontCell);
                        foreach(int neighbourcell in Neighbours) {
                            if(!CellIndexDoneList.Contains(neighbourcell) && !CellIndexFrontList.Contains(neighbourcell))
                                CellIndexFrontList.Add(neighbourcell);
                        }
                    }
                }
            }
            return CellIndexDoneList;
        }
    }



    /// <summary>
    /// refinement on cells which are on the specified boundary
    /// </summary>
    [Serializable]
    public class AMRInBoundingBox : AMRLevelIndicatorWithLevelset {
        //[DataMember]
        //public BoundingBox bb;
        [DataMember]
        public double[] point1;
        [DataMember]
        public double[] point2;
        public AMRInBoundingBox(double[] p1, double[] p2) {
            point1 = p1;
            point2 = p2;
        }

        public override int[] DesiredCellChanges() {



            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            var bb = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
            bb.AddPoint(point1);
            bb.AddPoint(point2);
            for (int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;
                

                bool refine = false;
                foreach (var cell in cells) {
                    if (bb.Contains(GridData.Cells.CellCenter.ExtractSubArrayShallow(j, -1).To1DArray()))
                        refine = true;
                    //try {
                    //    if (bb.Contains(GridData.Cells.CellCenter.ExtractSubArrayShallow(j, -1).To1DArray()))
                    //        refine = true;
                    //} catch (Exception) {
                    //    continue;
                    //}
                }

                if (refine && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } else if (!refine && currentLevel > 0) {
                    levels[j] = -1;
                    cellsToCoarse++;
                }
            }

            return levels;
        }

    
    }





}