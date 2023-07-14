﻿/* =======================================================================
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
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// three-dimensional grids
    /// </summary>
    [Serializable]
    public class Grid3D : GridCommons {

        /// <summary>
        /// Constructs an empty 3D grid.
        /// </summary>
        /// <param name="_RefElement"></param>
        public Grid3D(RefElement _RefElement)
            : base(new RefElement[] { _RefElement }, new RefElement[] { _RefElement.FaceRefElement }) {
        }


        /// <summary>
        /// Private constructor to support de-serialization
        /// </summary>
        private Grid3D()
            : base(new RefElement[] { Cube.Instance }, new RefElement[] { Square.Instance }) {
        }


        /// <summary>
        /// Constructs a Cartesian 3D Grid
        /// </summary>
        /// <param name="xNodes">
        /// Array containing the nodes in x-Direction
        /// </param>
        /// <param name="yNodes">
        /// Array containing the nodes in y-Direction
        /// </param>
        /// <param name="zNodes">
        /// Array containing the nodes in z-Direction
        /// </param>
        /// <param name="periodicX">
        /// Toggle for periodic boundary condition in x-direction
        /// </param>
        /// <param name="periodicY">
        /// Toggle for periodic boundary condition in y-direction
        /// </param>
        /// <param name="periodicZ">
        /// Toggle for periodic boundary condition in z-direction
        /// </param>
        /// <param name="_CellType">
        /// The specific type of hexahedral elements to be used.
        /// </param>
        /// <param name="CutOuts">
        /// Optional regions that are not meshed
        /// </param>
        /// <returns>
        /// A Cartesian 3D grid with the given nodes.
        /// </returns>
        public static Grid3D Cartesian3DGrid(double[] xNodes, double[] yNodes, double[] zNodes, 
            CellType _CellType = CellType.Cube_Linear,
            bool periodicX = false, bool periodicY = false, bool periodicZ = false, 
            params BoundingBox[] CutOuts) {
            using (var tr = new FuncTrace()) {
                MPICollectiveWatchDog.Watch();

                // Some Checks
                // ===========
                CheckMonotonicity(xNodes);
                CheckMonotonicity(yNodes);
                CheckMonotonicity(zNodes);

                // Does not work, so far
                //if (String.Compare("Cube",0,_CellType.CompareTo(),0,4,false)) {
                //    throw new ApplicationException("Grid must consist of Cubes");
                //}

                // split along x-Axis (cells can be redistributed by ParMETIS anyway)
                // ==================================================================
                int nX = xNodes.Length - 1;
                int nY = yNodes.Length - 1;
                int nZ = zNodes.Length - 1;

                if (nX < 1)
                    throw new ArithmeticException("At least 2 node are required.");
                if (nY < 1)
                    throw new ArithmeticException("At least 2 node are required.");
                if (nZ < 1)
                    throw new ArithmeticException("At least 2 node are required.");

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                int i0 = nX * myrank / size;        // 1st x-Index on this proc.
                int iE = nX * (myrank + 1) / size;  // 1st x-Index on next proc. rank

                // Return object
                // =============

                Grid3D grid;
                using (new BlockTrace("GridInstantiation", tr)) {
                    grid = new Grid3D(Cube.Instance);
                }

                // define periodic transformations, if necessary
                // =============================================
                byte perxTag = 0;
                byte peryTag = 0;
                byte perzTag = 0;

                if (periodicX) {
                    Vector[] Inlet = { new Vector { xNodes[0], yNodes[0], zNodes[0] },
                                       new Vector { xNodes[0], yNodes[nY], zNodes[0] },
                                       new Vector { xNodes[0], yNodes[0], zNodes[nZ] }
                                       };
                    Vector[] Outlet = { new Vector { xNodes[nX], yNodes[0], zNodes[0] },
                                        new Vector { xNodes[nX], yNodes[nY], zNodes[0] },
                                        new Vector { xNodes[nX], yNodes[0], zNodes[nZ] } 
                                        };
                    grid.ConstructPeriodicEdgeTrafo(Outlet, new double[] { 1.0, 0, 0 }, Inlet, new double[] { 1.0, 0, 0 }, out perxTag);
                    grid.EdgeTagNames.Add(perxTag, "Periodic-X");
                }

                if (periodicY) {
                    Vector[] Inlet = { new Vector { xNodes[0], yNodes[0], zNodes[0] },
                                       new Vector { xNodes[nX], yNodes[0], zNodes[0] },
                                       new Vector { xNodes[0], yNodes[0], zNodes[nZ] }
                                       };
                    Vector[] Outlet = { new Vector { xNodes[0], yNodes[nY], zNodes[0] },
                                        new Vector { xNodes[nX], yNodes[nY], zNodes[0] },
                                        new Vector { xNodes[0], yNodes[nY], zNodes[nZ] }
                                        };
                    grid.ConstructPeriodicEdgeTrafo(Outlet, new double[] { 0, 1.0, 0 }, Inlet, new double[] { 0, 1.0, 0 }, out peryTag);
                    grid.EdgeTagNames.Add(peryTag, "Periodic-Y");
                }

                if (periodicZ) {
                    Vector[] Inlet = { new Vector { xNodes[0], yNodes[0], zNodes[0] },
                                       new Vector { xNodes[0], yNodes[nY], zNodes[0] },
                                       new Vector { xNodes[nX], yNodes[0], zNodes[0] }
                                       };
                    Vector[] Outlet = { new Vector { xNodes[0], yNodes[0], zNodes[nZ] },
                                        new Vector { xNodes[0], yNodes[nY], zNodes[nZ] },
                                        new Vector { xNodes[nX], yNodes[0], zNodes[nZ] } };
                    grid.ConstructPeriodicEdgeTrafo(Outlet, new Vector { 0, 0, 1.0 }, Inlet, new Vector { 0, 0, 1.0 }, out perzTag);
                    grid.EdgeTagNames.Add(perzTag, "Periodic-Z");
                }

                // set cells
                // =========
                int LocalNoOfCells = (iE - i0) * (yNodes.Length - 1) * (zNodes.Length - 1);
                List<Cell> Cells = new List<Cell>(LocalNoOfCells);

                bool IsInCutOut(int i, int j, int k) {
                    if (i < -1)
                        throw new IndexOutOfRangeException();
                    if (i > nX + 1)
                        throw new IndexOutOfRangeException();
                    if (i < 0)
                        i = nX - 1;
                    if (i >= nX)
                        i = 0;

                    if (j < -1)
                        throw new IndexOutOfRangeException();
                    if (j > nY + 1)
                        throw new IndexOutOfRangeException();
                    if (j < 0)
                        j = nY - 1;
                    if (j >= nY)
                        j = 0;

                    if (k < -1)
                        throw new IndexOutOfRangeException();
                    if (k > nZ + 1)
                        throw new IndexOutOfRangeException();
                    if (k < 0)
                        k = nZ - 1;
                    if (k >= nZ)
                        k = 0;

                    if (CutOuts != null) {
                        double xC = 0.5 * (xNodes[i] + xNodes[(i + 1)]);
                        double yC = 0.5 * (yNodes[j] + yNodes[(j + 1)]);
                        double zC = 0.5 * (zNodes[k] + zNodes[(k + 1)]);

                        return (CutOuts.Any(BB => BB.Contains(xC, yC, zC)));

                    } else {
                        return false;
                    }
                }


                int cnt = 0;
                if (Math.Abs(i0 - iE) <= 0) {
                    //throw new ApplicationException("unable to do partitioning; X-Slice on processor " + myrank + " is empty; Try to load grid from database;");
                    cnt = 0;
                } else {

                    var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Cube));
                    MultidimensionalArray InterpolationNodes = Kref.GetInterpolationNodes(_CellType);
                    int NoOfNodes = Kref.GetInterpolationNodes(_CellType).GetLength(0);

                    for (int i = i0; i < iE; i++) {
                        for (int j = 0; j < nY; j++) {
                            for (int k = 0; k < nZ; k++) {

                                // cut-out regions test
                                // ====================
                                if (IsInCutOut(i, j, k)) {
                                    continue;
                                }
                                cnt++;

                                
                                Cell C_cnt = new Cell();
                                Cells.Add(C_cnt);

                                // define cell
                                // ===========
                                C_cnt.GlobalID = i + j * nX + k * nX * nY;
                                C_cnt.Type = _CellType;

                                C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 3);
                                Vector xyzPoint = new Vector(3);
                                //  var Bild0 = Cj0.TransformationParams;




                                for (int PointNumber = 0; PointNumber < NoOfNodes; PointNumber++) {
                                    xyzPoint[0] = xNodes[i] + (xNodes[i + 1] - xNodes[i]) * 0.5 * (InterpolationNodes[PointNumber, 0] + 1);
                                    xyzPoint[1] = yNodes[j] + (yNodes[j + 1] - yNodes[j]) * 0.5 * (InterpolationNodes[PointNumber, 1] + 1);
                                    xyzPoint[2] = zNodes[k] + (zNodes[k + 1] - zNodes[k]) * 0.5 * (InterpolationNodes[PointNumber, 2] + 1);
                                    // Write Physical Coordinates to TransformParams
                                    for (int dim = 0; dim < 3; dim++) {
                                        C_cnt.TransformationParams[PointNumber, dim] = xyzPoint[dim];
                                    }
                                }

                                // cell neighbourship
                                // ==================
                                
                                int[,] Vertices = new int[3, 8]; 
                                Vertices[0, 0] = 0;
                                Vertices[0, 2] = 0;
                                Vertices[0, 3] = 0;
                                Vertices[0, 7] = 0;
                                Vertices[0, 6] = 1;
                                Vertices[0, 5] = 1;
                                Vertices[0, 4] = 1;
                                Vertices[0, 1] = 1;

                                Vertices[1, 0] = 0;
                                Vertices[1, 1] = 0;
                                Vertices[1, 3] = 0;
                                Vertices[1, 6] = 0;
                                Vertices[1, 7] = 1;
                                Vertices[1, 5] = 1;
                                Vertices[1, 4] = 1;
                                Vertices[1, 2] = 1;

                                Vertices[2, 0] = 0;
                                Vertices[2, 1] = 0;
                                Vertices[2, 2] = 0;
                                Vertices[2, 4] = 0;
                                Vertices[2, 6] = 1;
                                Vertices[2, 7] = 1;
                                Vertices[2, 5] = 1;
                                Vertices[2, 3] = 1;

                                C_cnt.NodeIndices = new long[8];
                                for (int pointnumber = 0; pointnumber < 8; pointnumber++) {
                                    C_cnt.NodeIndices[pointnumber] = i + Vertices[0, pointnumber] + (j + Vertices[1, pointnumber]) * xNodes.Length + (k + Vertices[2, pointnumber]) * xNodes.Length * yNodes.Length;
                                }



                                // Edge tags (only periodic)
                                // =========================
                                if (periodicY || periodicX || periodicZ) {

                                    int iNeigh;
                                    int jNeigh;
                                    int kNeigh;


                                    if (periodicX) {
                                        iNeigh = i + 1;
                                        jNeigh = j;
                                        kNeigh = k;
                                        if (iNeigh >= nX && !IsInCutOut(iNeigh, jNeigh, kNeigh)) {
                                            (new CellFaceTag() {
                                                EdgeTag = perxTag,
                                                PeriodicInverse = false,
                                                FaceIndex = (int)Cube.Faces.Right,
                                                NeighCell_GlobalID = 0 + (jNeigh) * (nX) + kNeigh * (nX) * (nY)
                                            }).AddToArray(ref C_cnt.CellFaceTags);
                                            //System.Diagnostics.Debugger.Break();

                                        }
                                        iNeigh = i - 1;
                                        jNeigh = j;
                                        kNeigh = k;
                                        if (iNeigh < 0 && !IsInCutOut(iNeigh, jNeigh, kNeigh)) {
                                            (new CellFaceTag() {
                                                EdgeTag = perxTag,
                                                PeriodicInverse = true,
                                                FaceIndex = (int)Cube.Faces.Left,
                                                NeighCell_GlobalID = (nX - 1) + (jNeigh) * (nX) + kNeigh * nX * nY
                                            }).AddToArray(ref C_cnt.CellFaceTags);
                                            //System.Diagnostics.Debugger.Break();
                                        }
                                    }

                                    if (periodicY) {
                                        iNeigh = i;
                                        jNeigh = j + 1;
                                        kNeigh = k;
                                        if (jNeigh >= nY && !IsInCutOut(iNeigh, jNeigh, kNeigh)) {
                                            (new CellFaceTag() {
                                                EdgeTag = peryTag,
                                                PeriodicInverse = false,
                                                FaceIndex = (int)Cube.Faces.Top,
                                                NeighCell_GlobalID = iNeigh + (0) * (nX) + kNeigh * (nX) * (nY)
                                            }).AddToArray(ref C_cnt.CellFaceTags);
                                        }

                                        iNeigh = i;
                                        jNeigh = j - 1;
                                        kNeigh = k;
                                        if (jNeigh < 0 && !IsInCutOut(iNeigh, jNeigh, kNeigh)) {
                                            (new CellFaceTag() {
                                                EdgeTag = peryTag,
                                                PeriodicInverse = true,
                                                FaceIndex = (int)Cube.Faces.Bottom,
                                                NeighCell_GlobalID = iNeigh + (nY - 1) * (nX) + kNeigh * nX * nY
                                            }).AddToArray(ref C_cnt.CellFaceTags);
                                        }
                                    }

                                    if (periodicZ) {
                                        iNeigh = i;
                                        jNeigh = j;
                                        kNeigh = k + 1;
                                        if (kNeigh >= nZ && !IsInCutOut(iNeigh, jNeigh, kNeigh)) {
                                            (new CellFaceTag() {
                                                EdgeTag = perzTag,
                                                PeriodicInverse = false,
                                                FaceIndex = (int)Cube.Faces.Front,
                                                NeighCell_GlobalID = iNeigh + jNeigh * (nX) + 0 * (nX) * (nY)
                                            }).AddToArray(ref C_cnt.CellFaceTags);
                                        }

                                        iNeigh = i;
                                        jNeigh = j;
                                        kNeigh = k - 1;
                                        if (kNeigh < 0 && !IsInCutOut(iNeigh, jNeigh, kNeigh)) {
                                            (new CellFaceTag() {
                                                EdgeTag = perzTag,
                                                PeriodicInverse = true,
                                                FaceIndex = (int)Cube.Faces.Back,
                                                NeighCell_GlobalID = iNeigh + jNeigh * (nX) + (nZ - 1) * (nX) * (nY)
                                            }).AddToArray(ref C_cnt.CellFaceTags);
                                        }
                                    }

                                }

                                
                            }
                        }
                    }
                }
                grid.Cells = Cells.ToArray();

                cnt = cnt.MPISum();
                if(cnt <= 0) {
                    throw new ArgumentException("Grid is empty - check arguments (cut-outs).");
                }

               

                // return
                // ======
                //Console.WriteLine("Returning Cartesian 3D-Grid with {0} Cells", nX * nY * nZ);

                if (CutOuts != null && CutOuts.Length > 0) {
                    grid.CompressGlobalID();
                    grid.CompressNodeIndices();
                }

                foreach(var cl in grid.Cells) {
                    if (cl.GlobalID < 0)
                        throw new ApplicationException("Internal error - illegal GlobalID.");
                    if (cl.GlobalID >= cnt)
                        throw new ApplicationException("Internal error - illegal GlobalID.");
                }

                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                return grid;
            }
        }


        internal static double[] Param2XY(double rPoint, double sPoint) {

            // Calculate Coordinates
            double[] xyPoint = new double[] { 0, 0 };
            // x-Coordinate
            xyPoint[0] = rPoint * Math.Cos(2 * Math.PI * sPoint);
            // y-Coordinate
            xyPoint[1] = rPoint * Math.Sin(2 * Math.PI * sPoint);
            return xyPoint;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="rNodes">
        /// nodes in radial direction: first node must be greater than 0.0.
        /// </param>
        /// <param name="sNodes">
        /// 
        /// </param>
        /// <param name="zNodes"></param>
        /// <param name="type"></param>
        /// <param name="PeriodicZ"></param>
        /// <param name="PeriodicS"></param>
        /// <returns></returns>
        static public Grid3D CylinderGrid(double[] rNodes, double[] sNodes, double[] zNodes, CellType type, bool PeriodicZ = false, bool PeriodicS = true) {
            using (new FuncTrace()) {
                if (!(Cube.Instance).SupportedCellTypes.Contains(type))
                    throw new ArgumentOutOfRangeException("illegal cell type.");

                int nR = rNodes.Length - 1;
                int nS = sNodes.Length - 1;
                int nZ = zNodes.Length - 1;

                MPICollectiveWatchDog.Watch();
                Grid3D grid = new Grid3D(Cube.Instance);

                CheckMonotonicity(rNodes);
                CheckMonotonicity(sNodes);
                CheckMonotonicity(zNodes);

                if (nS < 3 && PeriodicS)
                    throw new ArithmeticException("At least 3 Elements in S-Direction are required for Periodic Boundary Condition to work");
                if (nZ < 3 && PeriodicZ)
                    throw new ArithmeticException("At least 3 Elements in Z-Direction are required for Periodic Boundary Condition to work");


                // Exceptions for bounds of r and s for a segment of a circle
                if ((Math.Abs(sNodes.First() - sNodes.Last()) > 1)) {
                    throw new ArgumentException("sNodes must not exceed an interval width of 1");
                };
                if (rNodes.First() <= 0) {
                    throw new ArgumentException("rNodes must be r>0");
                };
                //bool fullcircle = false;
                //if ((Math.Abs(sNodes.First() - sNodes.Last()) == 1)) {
                //    fullcircle = true;
                //}

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);



                // define periodic transformations
                // =============================================
                byte persTag = 0;
                if (PeriodicS) {
                    Vector[] Inlet = {
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.First())[0],
                            Param2XY(rNodes.First(), sNodes.First())[1],
                            zNodes.First()
                        },
                        new Vector {
                            Param2XY(rNodes.Last(), sNodes.First())[0],
                            Param2XY(rNodes.Last(), sNodes.First())[1],
                            zNodes.First()
                        },
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.First())[0],
                            Param2XY(rNodes.First(), sNodes.First())[1],
                            zNodes.Last()
                        }
                     };
                    Vector InletNormal = new Vector { 
                            Param2XY(rNodes.Last(), sNodes.First())[1] - Param2XY(rNodes.First(), sNodes.First())[1],
                            -(Param2XY(rNodes.Last(), sNodes.First())[0] - Param2XY(rNodes.First(), sNodes.First())[0]),
                            0.0
                        };
                    Vector[] Outlet = {
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.Last())[0],
                            Param2XY(rNodes.First(), sNodes.Last())[1],
                            zNodes.First()
                        },
                        new Vector {
                            Param2XY(rNodes.Last(), sNodes.Last())[0],
                            Param2XY(rNodes.Last(), sNodes.Last())[1],
                            zNodes.First()
                        },
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.Last())[0],
                            Param2XY(rNodes.First(), sNodes.Last())[1],
                            zNodes.Last()
                        }
                    };
                    Vector OutletNormal = new Vector {
                            Param2XY(rNodes.Last(), sNodes.Last())[1] - Param2XY(rNodes.First(), sNodes.Last())[1],
                            -(Param2XY(rNodes.Last(), sNodes.Last())[0] - Param2XY(rNodes.First(), sNodes.Last())[0]),
                            0.0
                        };

                    grid.ConstructPeriodicEdgeTrafo(Outlet, OutletNormal, Inlet, InletNormal, out persTag);
                    grid.EdgeTagNames.Add(persTag, "Periodic-S");
                }
                byte perzTag = 0;
                if (PeriodicZ) {
                    Vector[] Inlet = {
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.First())[0],
                            Param2XY(rNodes.First(), sNodes.First())[1],
                            zNodes[0]
                        },
                        new Vector {
                            Param2XY(rNodes.Last(), sNodes.First())[0],
                            Param2XY(rNodes.Last(), sNodes.First())[1],
                            zNodes[0] },
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.Last())[0],
                            Param2XY(rNodes.First(), sNodes.Last())[1],
                            zNodes[0] }
                        };
                    Vector[] Outlet = {
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.First())[0],
                            Param2XY(rNodes.First(), sNodes.First())[1],
                            zNodes.Last() },
                        new Vector {
                            Param2XY(rNodes.Last(), sNodes.First())[0],
                            Param2XY(rNodes.Last(), sNodes.First())[1],
                            zNodes.Last()
                        },
                        new Vector {
                            Param2XY(rNodes.First(), sNodes.Last())[0],
                            Param2XY(rNodes.First(), sNodes.Last())[1],
                            zNodes.Last()
                        }
                    };
                    grid.ConstructPeriodicEdgeTrafo(Outlet, new double[] { 0, 0, 1.0 }, Inlet, new double[] { 0, 0, 1.0 }, out perzTag);
                    grid.EdgeTagNames.Add(perzTag, "Periodic-Z");
                }

                var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Cube));
                MultidimensionalArray InterpolationNodes = Kref.GetInterpolationNodes(type);
                // Point in Reference Element: InterpolationNodes[NodeIndex,SpatialDimension]


                if (myrank == 0) {
                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++

                    // create cell
                    int J = nR * nS * nZ; // total number of cells
                    grid.Cells = new Cell[J];

                    int j = 0;//counter

                    int NoOfNodes = Kref.GetInterpolationNodes(type).GetLength(0);

                    for (int i = 0; i < (nR); i++) {
                        for (int k = 0; k < (nS); k++) {
                            for (int l = 0; l < nZ; l++) {
                                Cell Cj0 = new Cell();
                                Cj0.GlobalID = i + k * nR + l * nR * nS;
                                ;
                                Cj0.Type = type;

                                Cj0.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 3);
                                //  var Bild0 = Cj0.TransformationParams;




                                for (int PointNumber = 0; PointNumber < NoOfNodes; PointNumber++) {

                                    // Interpolate in rs-Domain according to these Coordinates

                                    double rL = rNodes[i];
                                    double rR = rNodes[i + 1];
                                    double sL = sNodes[k];
                                    double sR = sNodes[k + 1];
                                    double zL = zNodes[l];
                                    double zR = zNodes[l + 1];
                                    Vector rsPoint = new Vector(3);

                                    rsPoint[0] = rL + (rR - rL) * 0.5 * (InterpolationNodes[PointNumber, 0] + 1);
                                    rsPoint[1] = sL + (sR - sL) * 0.5 * (InterpolationNodes[PointNumber, 1] + 1);
                                    rsPoint[2] = zL + (zR - zL) * 0.5 * (InterpolationNodes[PointNumber, 2] + 1);

                                    // Convert Point into physical Coordinates

                                    Cj0.TransformationParams[PointNumber, 0] = Param2XY(rsPoint[0], rsPoint[1])[0];
                                    Cj0.TransformationParams[PointNumber, 1] = Param2XY(rsPoint[0], rsPoint[1])[1];
                                    Cj0.TransformationParams[PointNumber, 2] = rsPoint[2];
                                }


                                //if (!(Bild0[0, 0] < Bild0[1, 0]))
                                //    throw new ArgumentException("Elements Degenerated - Aborting");
                                //if (!(Bild0[2, 0] < Bild0[3, 0]))
                                //    throw new ArgumentException("Elements Degenerated- Aborting");


                                // ------------------ //
                                // Cell Neigbourships //
                                // ------------------ //
                                int[,] Vertices = new int[3, 8];
                                Vertices[0, 0] = 0;
                                Vertices[0, 2] = 0;
                                Vertices[0, 3] = 0;
                                Vertices[0, 7] = 0;
                                Vertices[0, 6] = 1;
                                Vertices[0, 5] = 1;
                                Vertices[0, 4] = 1;
                                Vertices[0, 1] = 1;

                                Vertices[1, 0] = 0;
                                Vertices[1, 1] = 0;
                                Vertices[1, 3] = 0;
                                Vertices[1, 6] = 0;
                                Vertices[1, 7] = 1;
                                Vertices[1, 5] = 1;
                                Vertices[1, 4] = 1;
                                Vertices[1, 2] = 1;

                                Vertices[2, 0] = 0;
                                Vertices[2, 1] = 0;
                                Vertices[2, 2] = 0;
                                Vertices[2, 4] = 0;
                                Vertices[2, 6] = 1;
                                Vertices[2, 7] = 1;
                                Vertices[2, 5] = 1;
                                Vertices[2, 3] = 1;


                                Cj0.NodeIndices = new long[8];

                                for (int pointnumber = 0; pointnumber < 8; pointnumber++) {
                                    Cj0.NodeIndices[pointnumber] = i + Vertices[0, pointnumber] + (k + Vertices[1, pointnumber]) * rNodes.Length + (l + Vertices[2, pointnumber]) * rNodes.Length * sNodes.Length;
                                }



                                // Edge tags (only periodic)
                                // =========================
                                if (PeriodicS || PeriodicZ) {

                                    int iNeigh;
                                    int kNeigh;
                                    int lNeigh;

                                    if (PeriodicS) {
                                        iNeigh = i;
                                        kNeigh = k + 1;
                                        lNeigh = l;
                                        if (kNeigh >= nS) {
                                            (new CellFaceTag() {
                                                EdgeTag = persTag,
                                                PeriodicInverse = false,
                                                FaceIndex = (int)Cube.Faces.Top,
                                                NeighCell_GlobalID = iNeigh + (0) * (nR) + lNeigh * (nR) * (nS)
                                            }).AddToArray(ref Cj0.CellFaceTags);
                                        }

                                        iNeigh = i;
                                        kNeigh = k - 1;
                                        lNeigh = l;
                                        if (kNeigh < 0) {
                                            (new CellFaceTag() {
                                                EdgeTag = persTag,
                                                PeriodicInverse = true,
                                                FaceIndex = (int)Cube.Faces.Bottom,
                                                NeighCell_GlobalID = iNeigh + (nS - 1) * (nR) + lNeigh * nR * nS
                                            }).AddToArray(ref Cj0.CellFaceTags);
                                        }
                                    }

                                    if (PeriodicZ) {
                                        iNeigh = i;
                                        kNeigh = k;
                                        lNeigh = l + 1;
                                        if (PeriodicZ && lNeigh >= nZ) {
                                            (new CellFaceTag() {
                                                EdgeTag = perzTag,
                                                PeriodicInverse = false,
                                                FaceIndex = (int)Cube.Faces.Front,
                                                NeighCell_GlobalID = iNeigh + kNeigh * (nR) + 0 * (nR) * (nS)
                                            }).AddToArray(ref Cj0.CellFaceTags);
                                        }

                                        iNeigh = i;
                                        kNeigh = k;
                                        lNeigh = l - 1;
                                        if (PeriodicZ && lNeigh < 0) {
                                            (new CellFaceTag() {
                                                EdgeTag = perzTag,
                                                PeriodicInverse = true,
                                                FaceIndex = (int)Cube.Faces.Back,
                                                NeighCell_GlobalID = iNeigh + kNeigh * (nR) + (nZ - 1) * (nR) * (nS)
                                            }).AddToArray(ref Cj0.CellFaceTags);
                                        }
                                    }

                                }

                                grid.Cells[Cj0.GlobalID] = Cj0;
                                //Console.WriteLine("Cell {0} of {1} created.", j, nR * nS * nZ);
                                j++;
                            }
                        }
                    }

                } else {
                    grid.Cells = new Cell[0];
                }

                return grid;

            }
        }


        /// <summary>
        /// Creates a grid in 3D with hanging nodes. 
        /// Idea: gridBoxes defines the individual region and its resolution. All boxes are put together
        /// and if one gridBox overlaps a previous one, it is cut out of the previous box. Thus, the boxes 
        /// should be ordered from large to small in terms of physical dimension  
        /// </summary>
        /// <param name="periodicX"></param>
        /// <param name="periodicY"></param>
        /// <param name="periodicZ"></param>
        /// <param name="gridBoxes"></param>
        /// <returns></returns>
        static public GridCommons HangingNodes3D(bool periodicX, bool periodicY, bool periodicZ, params GridBox[] gridBoxes) {
            if (gridBoxes.Length < 2) {
                throw new ArgumentException("At least 2 GridBoxes are needed for a HangingNodes grid, but only " + gridBoxes.Length + " are specified");
            }
            List<Grid3D> gridList = new List<Grid3D>(gridBoxes.Length);
            GridBox box = gridBoxes[0];
            int dimension = box.boundingBox.D;
            double[][] nodes = new double[3][];

            if (CheckConnectivity3D(box, gridBoxes[1]))
                Console.WriteLine("Adjusted GridBox[1] to ensure Connectivity");

            for (int d = 0; d < dimension; d++) {
                nodes[d] = GenericBlas.Linspace(box.boundingBox.Min[d], box.boundingBox.Max[d], box.numOfCells[d] + 1);
            }
            var grd = Grid3D.Cartesian3DGrid(nodes[0], nodes[1], nodes[2], periodicX: periodicX, periodicY: periodicY, CutOuts: gridBoxes[1].boundingBox);
            gridList.Add(grd);


            for (int i = 1; i < gridBoxes.Length; i++) {
                box = gridBoxes[i];
                if (i < gridBoxes.Length - 1) {
                    if (CheckConnectivity3D(box, gridBoxes[i + 1]))
                        Console.WriteLine("Adjusted GridBox[{0}] to ensure Connectivity", i + 1);
                }

                // Checks, if smaller Box has also periodic boundaries 
                bool localPeriodicX = false;
                bool localPeriodicY = false;
                bool localPeriodicZ = false;
                if (periodicX) {
                    if (gridBoxes[i - 1].boundingBox.Min[0] == box.boundingBox.Min[0] && gridBoxes[i - 1].boundingBox.Max[0] == box.boundingBox.Max[0]) {
                        localPeriodicX = true;
                    } else if ((gridBoxes[i - 1].boundingBox.Min[0] == box.boundingBox.Min[0] && gridBoxes[i - 1].boundingBox.Max[0] != box.boundingBox.Max[0])
                        || (gridBoxes[i - 1].boundingBox.Min[0] != box.boundingBox.Min[0] && gridBoxes[i - 1].boundingBox.Max[0] == box.boundingBox.Max[0])) {
                        throw new ArgumentException("The defined HangingNodes grid with periodic boundary in x direction has two boxes at one side, but not at corresponding other side, i.e the boxes have not the same length in x direction. Such a grid is not possible!");
                    }
                }
                if (periodicY) {
                    if (gridBoxes[i - 1].boundingBox.Min[1] == box.boundingBox.Min[1] && gridBoxes[i - 1].boundingBox.Max[1] == box.boundingBox.Max[1]) {
                        localPeriodicY = true;
                    } else if ((gridBoxes[i - 1].boundingBox.Min[1] == box.boundingBox.Min[1] && gridBoxes[i - 1].boundingBox.Max[1] != box.boundingBox.Max[1])
                        || (gridBoxes[i - 1].boundingBox.Min[1] != box.boundingBox.Min[1] && gridBoxes[i - 1].boundingBox.Max[1] == box.boundingBox.Max[1])) {
                        throw new ArgumentException("The defined HangingNodes grid with periodic boundary in y direction has two boxes at one side, but not at corresponding other side, i.e the boxes have not the same length in y direction. Such a grid is not possible!");
                    }
                }

                if (periodicZ) {
                    if (gridBoxes[i - 1].boundingBox.Min[2] == box.boundingBox.Min[2] && gridBoxes[i - 1].boundingBox.Max[2] == box.boundingBox.Max[2]) {
                        localPeriodicY = true;
                    } else if ((gridBoxes[i - 1].boundingBox.Min[2] == box.boundingBox.Min[2] && gridBoxes[i - 1].boundingBox.Max[2] != box.boundingBox.Max[2])
                        || (gridBoxes[i - 1].boundingBox.Min[2] != box.boundingBox.Min[2] && gridBoxes[i - 1].boundingBox.Max[2] == box.boundingBox.Max[2])) {
                        throw new ArgumentException("The defined HangingNodes grid with periodic boundary in z direction has two boxes at one side, but not at corresponding other side, i.e the boxes have not the same length in z direction. Such a grid is not possible!");
                    }
                }


                // Creating Grid for box with possible cut out region
                for (int d = 0; d < dimension; d++) {
                    nodes[d] = GenericBlas.Linspace(box.boundingBox.Min[d], box.boundingBox.Max[d], box.numOfCells[d] + 1);
                }

                if (i < gridBoxes.Length - 1) {
                    grd = Grid3D.Cartesian3DGrid(nodes[0], nodes[1], nodes[2], periodicX: localPeriodicX, periodicY: localPeriodicY, periodicZ: localPeriodicZ, CutOuts: gridBoxes[i + 1].boundingBox);
                } else {
                    grd = Grid3D.Cartesian3DGrid(nodes[0], nodes[1], nodes[2], periodicX: localPeriodicX, periodicY: localPeriodicY);
                }

                gridList.Add(grd);
            }
            var gridMerged = GridCommons.MergeLogically(gridList.ToArray());
            var grid = GridCommons.Seal(gridMerged, 4);
            return grid;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="coarse"></param>
        /// <param name="fine"></param>
        /// <returns></returns>
        private static bool CheckConnectivity3D(GridBox coarse, GridBox fine) {
            bool madeAdjustments = false;
            int dimension = coarse.boundingBox.D;
            double[][] nodes = new double[coarse.boundingBox.D][];
            for (int d = 0; d < dimension; d++) {
                nodes[d] = GenericBlas.Linspace(coarse.boundingBox.Min[d], coarse.boundingBox.Max[d], coarse.numOfCells[d] + 1);
            }

            double[,] cornerPoints = new double[(int)Math.Pow(2, dimension), dimension];
            // FrontBottomLeft
            cornerPoints[0, 0] = fine.boundingBox.Min[0];
            cornerPoints[0, 1] = fine.boundingBox.Min[1];
            cornerPoints[0, 2] = fine.boundingBox.Min[2];
            // FrontBottomRight
            cornerPoints[1, 0] = fine.boundingBox.Max[0];
            cornerPoints[1, 1] = fine.boundingBox.Min[1];
            cornerPoints[1, 2] = fine.boundingBox.Min[2];
            // FrontTopLeft
            cornerPoints[2, 0] = fine.boundingBox.Min[0];
            cornerPoints[2, 1] = fine.boundingBox.Max[1];
            cornerPoints[2, 2] = fine.boundingBox.Min[2];
            // FrontTopRight
            cornerPoints[3, 0] = fine.boundingBox.Max[0];
            cornerPoints[3, 1] = fine.boundingBox.Max[1];
            cornerPoints[3, 2] = fine.boundingBox.Min[2];
            // BackBottomLeft
            cornerPoints[4, 0] = fine.boundingBox.Min[0];
            cornerPoints[4, 1] = fine.boundingBox.Min[1];
            cornerPoints[4, 2] = fine.boundingBox.Max[2];
            // BackBottomRight
            cornerPoints[5, 0] = fine.boundingBox.Max[0];
            cornerPoints[5, 1] = fine.boundingBox.Min[1];
            cornerPoints[5, 2] = fine.boundingBox.Max[2];
            // BackTopLeft
            cornerPoints[6, 0] = fine.boundingBox.Min[0];
            cornerPoints[6, 1] = fine.boundingBox.Max[1];
            cornerPoints[6, 2] = fine.boundingBox.Max[2];
            // BackTopRight
            cornerPoints[7, 0] = fine.boundingBox.Max[0];
            cornerPoints[7, 1] = fine.boundingBox.Max[1];
            cornerPoints[7, 2] = fine.boundingBox.Max[2];

            //check every cornerPoint
            for (int i = 0; i < (int)Math.Pow(2, dimension); i++) {
                double[] pt = new double[] { cornerPoints[i, 0], cornerPoints[i, 1], cornerPoints[i, 2] };
                // Only check and adjust if the pt is in the coarse bounding box
                if (coarse.boundingBox.Contains(pt)) {

                    for (int d = 0; d < dimension; d++) {
                        for (int j = 0; j < nodes[d].Length - 1; j++) {
                            double x = nodes[d][j];
                            // Coarse and fine are connected in that value, everything is fine
                            if (Math.Abs(x - pt[d]) <= 1.0e-10) {
                                break;
                            } else {
                                //Found right interval
                                if (x - pt[d] < 0.0 && nodes[d][j + 1] - pt[d] > 0.0) {
                                    // Adjust fine boundingBox to nearest neighbor
                                    if (Math.Abs(x - pt[d]) < Math.Abs(nodes[d][j + 1] - pt[d])) {
                                        cornerPoints[i, d] = x;
                                    } else {
                                        cornerPoints[i, d] = nodes[d][j + 1];
                                    }

                                    if (Math.Abs(nodes[d][j + 1] - pt[d]) > 1.0e-10) {
                                        madeAdjustments = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (madeAdjustments) {
                fine.boundingBox.Clear();
                fine.boundingBox.AddPoints(cornerPoints);
            }

            return madeAdjustments;
        }

        /// <summary>
        /// Cutout of a sphere
        /// </summary>
        /// <param name="rNodes">Radial nodes, currently r=0 is not supported</param>
        /// <param name="phiNodes">Azimutal nodes, decoding multiples of 2 Pi </param>
        /// <param name="thetaNodes">Polar nodes, decoding multiples of Pi,
        /// supported range is (-0.5, 0.5) excluding the polar points/// </param>
        /// <param name="type"></param>
        /// <returns></returns>
        public static Grid3D SphereCutout(double[] rNodes, double[] phiNodes, double[] thetaNodes, CellType type = CellType.Cube_27) {
            using (new FuncTrace()) {
                if (!(Cube.Instance).SupportedCellTypes.Contains(type))
                    throw new ArgumentException("illegal cell type.");
                if (CellTypeExtensions.IsLinear(type)) {
                    throw new ArgumentException("illegal cell type.");
                }
                CheckMonotonicity(rNodes);
                CheckMonotonicity(phiNodes);
                CheckMonotonicity(thetaNodes);

                if ((Math.Abs(phiNodes.First() - phiNodes.Last()) > 1)) {
                    throw new ArgumentException("phiNodes must not exceed an interval width of 1");
                };
                if ((Math.Abs(thetaNodes.First() - thetaNodes.Last()) > 1)) {
                    throw new ArgumentException("thetaNodes must not exceed an interval width of 1");
                };
                if (thetaNodes.First() <-0.5 | thetaNodes.Last() > 0.5) {
                    throw new ArgumentException("thetaNodes must be in range [-0.5, 0.5]");
                };
                if (rNodes.First() < 0) {
                    throw new ArgumentException("rNodes must be r>0"); // last segment should be filled with a pyramid, currently multiple reference elements are not supported...
                };

                bool zeroR = false;
                if ((Math.Abs(rNodes.First()) == 0.0)) {
                    zeroR = true;
                    throw new NotImplementedException();
                }
                bool fullcirclePhi = false;
                if ((Math.Abs(phiNodes.First() - phiNodes.Last()) == 1)) {
                    fullcirclePhi = true;
                }
                bool quartercircleThetaMinus = false;
                if ((Math.Abs(thetaNodes.First()) == 0.5)) {
                    throw new NotImplementedException();
                    //quartercircleThetaMinus = true;
                }
                bool quartercircleThetaPlus = false;
                if ((Math.Abs(thetaNodes.Last()) == 0.5)) {
                    throw new NotImplementedException();
                    //quartercircleThetaPlus = true;
                }

                int NoOfrNodes = rNodes.Length;
                int NoOfphiNodes = phiNodes.Length;
                int NoOfthetaNodes = thetaNodes.Length;

                MPICollectiveWatchDog.Watch();
                Grid3D grid = new Grid3D(Cube.Instance);

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (myrank == 0) {

                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++
                    List<Cell> AllCells = new List<Cell>();

                    // Reference element nodes
                    // =======================
                    var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Cube));
                    NodeSet InterpolationNodes = Kref.GetInterpolationNodes(type);
                    int NoOfNodes = InterpolationNodes.NoOfNodes;

                    // Allocate GlobalID's
                    long[,,] Gid = new long[NoOfrNodes - 1, NoOfphiNodes - 1, NoOfthetaNodes - 1];
                    long cnt = 0;
                    for (int iX = 0; iX < (NoOfrNodes - 1); iX++) {
                        for (int iY = 0; iY < (NoOfphiNodes - 1); iY++) {
                            for (int iZ = 0; iZ < (NoOfthetaNodes - 1); iZ++) {
                                Gid[iX, iY, iZ] = cnt;
                                cnt++;
                            }
                        }
                    }

                    // meshing of center section
                    // =========================

                    for (int iX = 0; iX < (NoOfrNodes - 1); iX++) {
                        for (int iY = 0; iY < (NoOfphiNodes - 1); iY++) {
                            for (int iZ = 0; iZ < (NoOfthetaNodes - 1); iZ++) {
                                // create cell
                                Cell Cj = new Cell();
                                Cj.GlobalID = Gid[iX, iY, iZ];
                                Cj.Type = type;
                                AllCells.Add(Cj);

                                // physical coordinates
                                Cj.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 3);

                                double r0 = rNodes[iX];
                                double r1 = rNodes[iX + 1];
                                double phi0 = phiNodes[iY];
                                double phi1 = phiNodes[iY + 1];
                                double theta0 = thetaNodes[iZ];
                                double theta1 = thetaNodes[iZ + 1];

                                for (int k = 0; k < NoOfNodes; k++) {
                                    double r = 0.5 * ((1-InterpolationNodes[k,0]) * r0 + (1 + InterpolationNodes[k, 0]) * r1);
                                    double phi = 0.5 * ((1 - InterpolationNodes[k, 1]) * phi0 + (1 + InterpolationNodes[k, 1]) * phi1);
                                    double theta = 0.5 * ((1 - InterpolationNodes[k, 2]) * theta0 + (1 + InterpolationNodes[k, 2]) * theta1);
                                    double x = r * Math.Cos(2 * Math.PI * phi) * Math.Cos(Math.PI * theta);
                                    double y = r * Math.Sin(2 * Math.PI * phi) * Math.Cos(Math.PI * theta);
                                    double z = r * Math.Sin(Math.PI * theta);
                                    Cj.TransformationParams[k, 0] = x;
                                    Cj.TransformationParams[k, 1] = y;
                                    Cj.TransformationParams[k, 2] = z;
                                }

                                // node indices (neighborship via cell face tags 
                                int NoRNodes = NoOfrNodes;
                                int NoPhiNodes = NoOfphiNodes;
                                int NoThetaNodes = NoOfphiNodes;
                                int offset = 0;

                                if (zeroR) { NoRNodes--; }
                                if (fullcirclePhi) { NoPhiNodes--; }
                                if (quartercircleThetaMinus) { offset = -(NoPhiNodes - 1) * NoRNodes; }

                                long[] NodeIndices = new long[] {
                                    (iZ*NoPhiNodes + iY) * NoRNodes + iX + offset,               // (-1,-1,-1)
                                    (iZ*NoPhiNodes + iY) * NoRNodes + iX + 1 + offset,           // (1,-1,-1)
                                    (iZ*NoPhiNodes + iY + 1) * NoRNodes + iX + offset,           // (-1,1,-1)
                                    ((iZ + 1)*NoPhiNodes + iY) * NoRNodes + iX + offset,         // (-1,-1,1)
                                    (iZ*NoPhiNodes + iY + 1) * NoRNodes + iX + 1 + offset,       // (1,1,-1)
                                    ((iZ + 1)*NoPhiNodes + iY + 1) * NoRNodes + iX + 1 + offset, // (1,1,1)
                                    ((iZ + 1)*NoPhiNodes + iY) * NoRNodes + iX + 1 + offset,     // (1,-1,1)
                                    ((iZ + 1)*NoPhiNodes + iY + 1) * NoRNodes + iX + offset      // (-1,1,1)
                                };
                                
                                if (fullcirclePhi && iY == NoOfphiNodes - 2) {
                                    NodeIndices[2] = (iZ * NoPhiNodes) * NoRNodes + iX + offset;
                                    NodeIndices[4] = (iZ * NoPhiNodes) * NoRNodes + iX + 1 + offset;
                                    NodeIndices[5] = ((iZ + 1) * NoPhiNodes) * NoRNodes + iX + 1 + offset;
                                    NodeIndices[7] = ((iZ + 1) * NoPhiNodes) * NoRNodes + iX + offset;
                                }
                                // theoretical node ids if center point/polar points are part of the domain
                                // however, with only cubes this leads to degenerate cells...
                                if (quartercircleThetaMinus && iZ == 0) {
                                    NodeIndices[0] = iX;
                                    NodeIndices[1] = iX + 1;
                                    NodeIndices[2] = iX;
                                    NodeIndices[4] = iX + 1;
                                }
                                if (quartercircleThetaPlus && iZ == NoOfthetaNodes - 2) {
                                    NodeIndices[3] = ((iZ + 1) * NoPhiNodes) * NoRNodes + iX + offset;
                                    NodeIndices[5] = ((iZ + 1) * NoPhiNodes) * NoRNodes + iX + 1 + offset;
                                    NodeIndices[6] = ((iZ + 1) * NoPhiNodes) * NoRNodes + iX + 1 + offset;
                                    NodeIndices[7] = ((iZ + 1) * NoPhiNodes) * NoRNodes + iX + offset;
                                }
                                if (zeroR && iX == 0) {
                                    NodeIndices[0] = 0;
                                    NodeIndices[2] = 0;
                                    NodeIndices[3] = 0;
                                    NodeIndices[7] = 0;
                                }                                

                                Cj.NodeIndices = NodeIndices;
                            }
                        }
                    }

                    grid.Cells = AllCells.ToArray();
                } else {
                    grid.Cells = new Cell[0];
                }

                return grid;
            }
        }


        /// <summary>
        /// Generates a so-called O-grid of a cylindrical domain.
        /// </summary>
        /// <param name="Radius">Radius of the circular domain.</param>
        /// <param name="CenterSectionWidth">With of center section which is meshed Cartesian.</param>
        /// <param name="NoOfCenterNodes">
        /// Number of nodes in the center section, in each direction. this also determines the 
        /// Number of nodes in rotational direction for each of the four ring segments.
        /// </param>
        /// <param name="NoOfRadialNodes">Number of nodes in the radial section.</param>
        /// <param name="type">Cell type.</param>
        /// <param name="zNodes">Nodes along the cylinder heigh.</param>
        /// <returns>
        /// A block-structured O-grid.
        /// </returns>
        public static Grid3D Ogrid(double CenterSectionWidth, double Radius, int NoOfCenterNodes, int NoOfRadialNodes, double[] zNodes, CellType type = CellType.Cube_8) {
            using(new FuncTrace()) {
                if(!(Cube.Instance).SupportedCellTypes.Contains(type))
                    throw new ArgumentException("illegal cell type.");
                if(CellTypeExtensions.IsLinear(type)) {
                    throw new ArgumentException("illegal cell type.");
                }
                if(Radius <= 0)
                    throw new ArgumentOutOfRangeException();
                if(CenterSectionWidth <= 0)
                    throw new ArgumentOutOfRangeException();
                if(Radius * Radius <= 0.5 * CenterSectionWidth * CenterSectionWidth)
                    throw new ArgumentOutOfRangeException();
                if(NoOfCenterNodes < 2)
                    throw new ArgumentOutOfRangeException();
                if(NoOfRadialNodes < 2)
                    throw new ArgumentOutOfRangeException();

                CheckMonotonicity(zNodes);
                int NoOfZNodes = zNodes.Length;
               
                MPICollectiveWatchDog.Watch();
                Grid3D grid = new Grid3D(Cube.Instance);
                
                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                

                

                if(myrank == 0) {


                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++
                    List<Cell> AllCells = new List<Cell>();
                    int PtCount = 0;

                    // Reference element nodes
                    // =======================
                    var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Cube));
                    NodeSet InterpolationNodes = Kref.GetInterpolationNodes(type);
                    int NoOfNodes = InterpolationNodes.NoOfNodes;

                    // Allocate GlobalID's
                    long[,,] CenterSection_Gid = new long[NoOfCenterNodes - 1, NoOfCenterNodes - 1, NoOfZNodes - 1];
                    long cnt = 0;
                    for (int iX = 0; iX < (NoOfCenterNodes - 1); iX++) {
                        for (int iY = 0; iY < (NoOfCenterNodes - 1); iY++) {
                            for (int iZ = 0; iZ < (NoOfZNodes - 1); iZ++) {
                                CenterSection_Gid[iX, iY, iZ] = cnt;
                                cnt++;
                            }
                        }
                    }

                    long[][,,] RingSections_Gid = new long[4][,,];
                    for(int _iRing = 0; _iRing < 4; _iRing++) {
                        RingSections_Gid[_iRing] = new long[NoOfRadialNodes - 1, NoOfCenterNodes - 1, NoOfZNodes - 1];
                        for(int iR = 0; iR < (NoOfRadialNodes - 1); iR++) {
                            for(int iPhi = 0; iPhi < (NoOfCenterNodes - 1); iPhi++) {
                                for (int iZ = 0; iZ < (NoOfZNodes - 1); iZ++) {
                                    RingSections_Gid[_iRing][iR, iPhi, iZ] = cnt;
                                    cnt++;
                                }
                            }
                        }
                    }


                     // meshing of center section
                     // =========================

                     double[] CenterNodes = GenericBlas.Linspace(-CenterSectionWidth * 0.5, CenterSectionWidth * 0.5, NoOfCenterNodes);


                    for (int iX = 0; iX < (NoOfCenterNodes - 1); iX++) {
                        for (int iY = 0; iY < (NoOfCenterNodes - 1); iY++) {
                            for (int iZ = 0; iZ < (NoOfZNodes - 1); iZ++) {
                                // create cell
                                Cell Cj = new Cell();
                                Cj.GlobalID = CenterSection_Gid[iX, iY, iZ];
                                Cj.Type = type;
                                AllCells.Add(Cj);

                                // physical coordinates
                                Cj.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 3);

                                double x0 = CenterNodes[iX];
                                double x1 = CenterNodes[iX + 1];
                                double y0 = CenterNodes[iY];
                                double y1 = CenterNodes[iY + 1];
                                double z0 = zNodes[iZ];
                                double z1 = zNodes[iZ + 1];

                                for (int k = 0; k < NoOfNodes; k++) {
                                    double x = 0.5 * (x1 - x0) * InterpolationNodes[k, 0] + 0.5 * (x1 + x0);
                                    double y = 0.5 * (y1 - y0) * InterpolationNodes[k, 1] + 0.5 * (y1 + y0);
                                    double z = 0.5 * (z1 - z0) * InterpolationNodes[k, 2] + 0.5 * (z1 + z0);
                                    Cj.TransformationParams[k, 0] = x;
                                    Cj.TransformationParams[k, 1] = y;
                                    Cj.TransformationParams[k, 2] = z;
                                }

                                // node indices (neighborship via cell face tags
                                Cj.NodeIndices = new long[] { PtCount, PtCount + 1, PtCount + 2, PtCount + 3, PtCount + 4, PtCount + 5, PtCount + 6, PtCount + 7 };
                                PtCount += 8;

                                // neigborship
                                if (iX > 0) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Left,
                                        NeighCell_GlobalID = CenterSection_Gid[iX - 1, iY, iZ],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                } else {
                                    // connection to left ring segment
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Left,
                                        NeighCell_GlobalID = RingSections_Gid[0][0, iY, iZ],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }

                                if (iX < (NoOfCenterNodes - 2)) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = CenterSection_Gid[iX + 1, iY, iZ],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                } else {
                                    // connection to right ring segment
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = RingSections_Gid[1][0, NoOfCenterNodes - 2 - iY, iZ],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }

                                if (iY > 0) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Bottom,
                                        NeighCell_GlobalID = CenterSection_Gid[iX, iY - 1, iZ]
                                    }, ref Cj.CellFaceTags);
                                } else {
                                    // connection to bottom ring segment
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Bottom,
                                        NeighCell_GlobalID = RingSections_Gid[3][0, NoOfCenterNodes - 2 - iX, iZ]
                                    }, ref Cj.CellFaceTags);
                                }

                                if (iY < (NoOfCenterNodes - 2)) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Top,
                                        NeighCell_GlobalID = CenterSection_Gid[iX, iY + 1, iZ]
                                    }, ref Cj.CellFaceTags);
                                } else {
                                    // connection to top ring segment
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Top,
                                        NeighCell_GlobalID = RingSections_Gid[2][0, iX, iZ]
                                    }, ref Cj.CellFaceTags);
                                }

                            }
                        }
                    }

                    // meshing of ring section
                    // ========================
                    double[] alphaS = GenericBlas.Linspace(0, 1, NoOfCenterNodes);
                    double[] betaS = GenericBlas.Linspace(0, 1, NoOfRadialNodes);

                    for (int iRing = 0; iRing < 4; iRing++) {
                        int iPrevRing, iNextRing;//, iPrevCon, iNextCon;
                        double phi_st, phi_en;
                        double x_st, x_en, y_st, y_en;

                        int inc_iX, inc_iY, OiX, OiY;

                        switch (iRing) {
                            case 0: // left
                            iNextRing = 2;
                            iPrevRing = 3;
                            //iPrevCon = NoOfCenterNodes - 2;
                            //iNextCon = 0;
                            inc_iX = 0;
                            inc_iY = 1;
                            OiX = 0;
                            OiY = 0;

                            phi_st = Math.PI * 0.25 + 2 * Math.PI * 0.5;
                            phi_en = Math.PI * 0.25 + 1 * Math.PI * 0.5;

                            x_st = -CenterSectionWidth * 0.5;
                            x_en = -CenterSectionWidth * 0.5;
                            y_st = -CenterSectionWidth * 0.5;
                            y_en = +CenterSectionWidth * 0.5;
                            break;

                            case 1: // right
                            iNextRing = 3; // bottom
                            iPrevRing = 2; // top
                            //iPrevCon = 0;
                            //iNextCon = NoOfCenterNodes - 2;
                            inc_iX = 0;
                            inc_iY = -1;
                            OiX = NoOfCenterNodes - 2;
                            OiY = NoOfCenterNodes - 2;

                            phi_st = +Math.PI * 0.25;
                            phi_en = -Math.PI * 0.25;

                            x_st = +CenterSectionWidth * 0.5;
                            x_en = +CenterSectionWidth * 0.5;
                            y_st = +CenterSectionWidth * 0.5;
                            y_en = -CenterSectionWidth * 0.5;
                            break;

                            case 2: // top
                            iNextRing = 1; // right
                            iPrevRing = 0; // left
                            inc_iX = 1;
                            inc_iY = 0;
                            OiX = 0;
                            OiY = NoOfCenterNodes - 2;

                            phi_st = +Math.PI * 0.25 + 1 * Math.PI * 0.5;
                            phi_en = +Math.PI * 0.25;

                            x_st = -CenterSectionWidth * 0.5;
                            x_en = +CenterSectionWidth * 0.5;
                            y_st = +CenterSectionWidth * 0.5;
                            y_en = +CenterSectionWidth * 0.5;
                            break;

                            case 3: // bottom
                            iNextRing = 0; // left
                            iPrevRing = 1; // right
                            inc_iX = -1;
                            inc_iY = 0;
                            OiX = NoOfCenterNodes - 2;
                            OiY = 0;

                            phi_st = +Math.PI * 0.25 + 3 * Math.PI * 0.5;
                            phi_en = +Math.PI * 0.25 + 2 * Math.PI * 0.5;

                            x_st = +CenterSectionWidth * 0.5;
                            x_en = -CenterSectionWidth * 0.5;
                            y_st = -CenterSectionWidth * 0.5;
                            y_en = -CenterSectionWidth * 0.5;
                            break;

                            default: throw new Exception();
                        }



                        for (int iPhi = 0; iPhi < (NoOfCenterNodes - 1); iPhi++) {
                            double phi0 = phi_st * (1.0 - alphaS[iPhi]) + phi_en * alphaS[iPhi];
                            double phi1 = phi_st * (1.0 - alphaS[iPhi + 1]) + phi_en * alphaS[iPhi + 1];

                            double xC0 = x_st * (1.0 - alphaS[iPhi]) + x_en * alphaS[iPhi];
                            double xC1 = x_st * (1.0 - alphaS[iPhi + 1]) + x_en * alphaS[iPhi + 1];
                            double yC0 = y_st * (1.0 - alphaS[iPhi]) + y_en * alphaS[iPhi];
                            double yC1 = y_st * (1.0 - alphaS[iPhi + 1]) + y_en * alphaS[iPhi + 1];

                            for (int iR = 0; iR < (NoOfRadialNodes - 1); iR++) {
                                for (int iZ = 0; iZ < NoOfZNodes - 1; iZ++) {
                                    //double beta0 = betaS[NoOfRadialNodes - iR - 1];
                                    //double beta1 = betaS[NoOfRadialNodes - iR - 2];
                                    double beta0 = betaS[iR];
                                    double beta1 = betaS[iR + 1];
                                    double z0 = zNodes[iZ];
                                    double z1 = zNodes[iZ + 1];

                                    // create cell
                                    Cell Cj = new Cell();
                                    Cj.GlobalID = RingSections_Gid[iRing][iR, iPhi, iZ];
                                    Cj.Type = type;
                                    AllCells.Add(Cj);

                                    // set nodes
                                    Cj.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 3);
                                    for (int k = 0; k < NoOfNodes; k++) {
                                        double phi = 0.5 * (phi1 - phi0) * InterpolationNodes[k, 0] + 0.5 * (phi1 + phi0);
                                        double xA = Math.Cos(phi) * Radius;
                                        double yA = Math.Sin(phi) * Radius;

                                        double xC = 0.5 * (xC1 - xC0) * InterpolationNodes[k, 0] + 0.5 * (xC1 + xC0);
                                        double yC = 0.5 * (yC1 - yC0) * InterpolationNodes[k, 0] + 0.5 * (yC1 + yC0);

                                        double beta = 0.5 * (beta1 - beta0) * InterpolationNodes[k, 1] + 0.5 * (beta1 + beta0);

                                        double x = xC * (1 - beta) + xA * beta;
                                        double y = yC * (1 - beta) + yA * beta;

                                        double z = 0.5 * (z1 - z0) * InterpolationNodes[k, 2] + 0.5 * (z1 + z0);

                                        Cj.TransformationParams[k, 0] = x;
                                        Cj.TransformationParams[k, 1] = y;
                                        Cj.TransformationParams[k, 2] = z;
                                    }

                                    // node indices (neighborship via cell face tags
                                    Cj.NodeIndices = new long[] { PtCount, PtCount + 1, PtCount + 2, PtCount + 3, PtCount + 4, PtCount + 5, PtCount + 6, PtCount + 7 };
                                    PtCount += 8;

                                    // neigborship
                                    if (iR > 0) {
                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Bottom,
                                            NeighCell_GlobalID = RingSections_Gid[iRing][iR - 1, iPhi, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    } else {
                                        // connection to center section

                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Bottom,
                                            NeighCell_GlobalID = CenterSection_Gid[OiX + iPhi * inc_iX, OiY + iPhi * inc_iY, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    }

                                    if (iR < (NoOfRadialNodes - 2)) {
                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Top,
                                            NeighCell_GlobalID = RingSections_Gid[iRing][iR + 1, iPhi, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    } else {
                                        //ArrayTools.AddToArray(new CellFaceTag() {
                                        //        FaceIndex = (int) Square.Edge.Right,
                                        //        NeighCell_GlobalID = RingSections_Gid[1][0, iY]
                                        //    }, ref Cj.CellFaceTags);
                                    }

                                    if (iPhi > 0) {
                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Left,
                                            NeighCell_GlobalID = RingSections_Gid[iRing][iR, iPhi - 1, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    } else {
                                        // connection to previous ring segment
                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Left,
                                            NeighCell_GlobalID = RingSections_Gid[iPrevRing][iR, NoOfCenterNodes - 2, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    }

                                    if (iPhi < (NoOfCenterNodes - 2)) {
                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Right,
                                            NeighCell_GlobalID = RingSections_Gid[iRing][iR, iPhi + 1, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    } else {
                                        // connection to next ring segment
                                        ArrayTools.AddToArray(new CellFaceTag() {
                                            FaceIndex = (int)Square.Faces.Right,
                                            NeighCell_GlobalID = RingSections_Gid[iNextRing][iR, 0, iZ],
                                            ConformalNeighborship = true
                                        }, ref Cj.CellFaceTags);
                                    }
                                }
                            }
                        }
                    }
                    // finalize
                    // ========

                    grid.Cells = AllCells.ToArray();
                } else {
                    grid.Cells = new Cell[0];
                }
                return grid;

            }
        }


    }
}
