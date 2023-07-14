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
using System.Diagnostics;
using System.Linq;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// a two-dimensional grid.
    /// </summary>
    [Serializable]
    public class Grid2D : GridCommons {

        /// <summary>
        /// Constructs an empty 2D grid.
        /// </summary>
        /// <param name="_RefElement">
        /// The <see cref="RefElement"/> for this grid. Currently only grids with _one_ RefElement are supported!!!
        /// </param>
        public Grid2D(RefElement _RefElement)
            : base(new RefElement[] { _RefElement }, new RefElement[] { _RefElement.FaceRefElement }) //
        {
            if (_RefElement.SpatialDimension != 2)
                throw new ArgumentException();
        }


        /// <summary>
        /// empty ctor. to support cloning.
        /// </summary>
        private Grid2D()
            : base(new RefElement[] { Square.Instance }, new RefElement[] { Line.Instance }) {
        }

        /*
        static private bool IsInCutoutRegion(int indX, int indY, double[] xNodes, double[] yNodes, Vector[] cutoutMin, Vector[] cutoutMax) {
            double Xcenter = 0.5 * (xNodes[indX] + xNodes[indX + 1]);
            double Ycenter = 0.5 * (yNodes[indY] + yNodes[indY + 1]);

            for (int l = cutoutMin.Length - 1; l >= 0; l--) {
                if(cutoutMin[l].Dim != 2)
                    throw new ArgumentException("expecting a 2D vector");
                if(cutoutMax[l].Dim != 2)
                    throw new ArgumentException("expecting a 2D vector");

                if (Xcenter >= cutoutMin[l].x && Xcenter <= cutoutMax[l].x
                    && Ycenter >= cutoutMin[l].y && Ycenter <= cutoutMax[l].y)
                    return true;
            }

            return false;
        }
        */

        /*
        /// <summary>
        /// constructs a new 2D Grid
        /// </summary>
        /// <param name="xNodes"></param>
        /// <param name="yNodes"></param>
        /// <param name="periodicX">turn on periodic boundary conditions in x direction</param>
        /// <param name="periodicY">turn on periodic boundary conditions in y direction</param>
        /// <param name="CutOuts">Optional regions that are not meshed</param>
        /// <param name="type">The type of the cells to be used</param>
        public static Grid2D Tutorium2DGrid(double[] xNodes, double[] yNodes, CellType type = CellType.Square_Linear, bool periodicX = true, bool periodicY = true, params BoundingBox[] CutOuts)
        {
            using (var tr = new FuncTrace())
            {
                MPICollectiveWatchDog.Watch();
                CheckMonotonicity(xNodes);
                CheckMonotonicity(yNodes);
                Grid2D grid;
                using (new BlockTrace("grid_object_instantiation", tr))
                {
                    grid = new Grid2D(Square.Instance);
                }

                // split along x-Axis (cells can be redistributed by ParMETIS anyway)
                // ==================================================================
                int nX = xNodes.Length - 1;
                int nY = yNodes.Length - 1;

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                int i0 = nX * myrank / size;        // 1st x-Index on this proc.
                int iE = nX * (myrank + 1) / size;  // 1st x-Index on next proc. rank

                if (Math.Abs(i0 - iE) <= 0)
                {
                    //throw new ApplicationException("unable to do partitioning; X-Slice on processor " + myrank + " is empty; Try to load grid from database;");
                    grid.Cells = new Cell[0];
                }
                else
                {

                    int D = grid.SpatialDimension;
                    var Kref = grid.RefElements.Single(re => re.GetType() == typeof(Square));
                    if (!Kref.SupportedCellTypes.Contains(type))
                        throw new ArgumentException("unsupported cell type");

                    // define periodic transformations, if necessary
                    // =============================================
                    byte perxTag = 0;
                    byte peryTag = 0;

                    if (periodicX)
                    {
                        double[][] x = { new double[] { xNodes[0], yNodes[0] }, new double[] { xNodes[0], yNodes[nY] } };
                        double[][] y = { new double[] { xNodes[nX], yNodes[0] - 1 }, new double[] { xNodes[nX], yNodes[nY] } };
                        grid.ConstructPeriodicEdgeTrafo(y, new double[] { 1.0, 0 }, x, new double[] { 1.0, 0 }, out perxTag);
                        grid.EdgeTagNames.Add(perxTag, "Periodic-X");
                    }

                    if (periodicY)
                    {
                        double[][] x = { new double[] { xNodes[0], yNodes[0] }, new double[] { xNodes[nX], yNodes[0] - 1 } };
                        double[][] y = { new double[] { xNodes[0], yNodes[nY] }, new double[] { xNodes[nX], yNodes[nY] } };
                        grid.ConstructPeriodicEdgeTrafo(y, new double[] { 0, 1.0 }, x, new double[] { 0, 1.0 }, out peryTag);
                        grid.EdgeTagNames.Add(peryTag, "Periodic-Y");
                    }

                    // set cells
                    // =========
                    int LocalNoOfCells = (iE - i0) * (yNodes.Length - 1);
                    var RefNodes = Kref.GetInterpolationNodes(type);
                    int NoOfNodes = Kref.GetInterpolationNodes(type).GetLength(0);

                    List<Cell> Cells = new List<Cell>(LocalNoOfCells);

                    int cnt = -1;
                    //int Cellcounter = 0;
                    for (int i = i0; i < iE; i++)
                    {
                        for (int j = 0; j < nY; j++)
                        {
                            cnt++;

                            // cut-out regions test
                            // ====================

                            if (CutOuts != null)
                            {
                                double xC = 0.5 * (xNodes[i] + xNodes[i + 1]);
                                double yC = 0.5 * (yNodes[j] + yNodes[j + 1]);

                                if (CutOuts.Any(BB => BB.Contains(xC, yC)))
                                    continue;
                            }

                            Cell C_cnt = new Cell();
                            Cells.Add(C_cnt);

                            // define cell
                            // ===========
                            C_cnt.GlobalID = i * (yNodes.Length - 1) + j;
                            C_cnt.Type = type;

                            // transformation
                            // ==============
                            C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, D);
                            {
                                NoOfNodes = RefNodes.GetLength(0);
                                Debug.Assert(RefNodes.GetLength(1) == D);

                                //C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, D);

                                double xL = xNodes[i];
                                double xR = xNodes[i + 1];
                                double yL = yNodes[j];
                                double yR = yNodes[j + 1];

                                for (int iNode = 0; iNode < NoOfNodes; iNode++)
                                {
                                    double xi = 0.5 * (RefNodes[iNode, 0] + 1.0);
                                    double eta = 0.5 * (RefNodes[iNode, 1] + 1.0);

                                    C_cnt.TransformationParams[iNode, 0] = xL * (1.0 - xi) + xR * xi;
                                    C_cnt.TransformationParams[iNode, 1] = yL * (1.0 - eta) + yR * eta;
                                }
                            }

                            // cell neighbourship
                            // ==================
                            C_cnt.NodeIndices = new int[4];
                            C_cnt.NodeIndices[0] = i + j * xNodes.Length;
                            C_cnt.NodeIndices[1] = (i + 1) + j * xNodes.Length;
                            C_cnt.NodeIndices[2] = i + (j + 1) * xNodes.Length;
                            C_cnt.NodeIndices[3] = (i + 1) + (j + 1) * xNodes.Length;

                            // Edge tags (only periodic)
                            // =========================
                            if (periodicY || periodicX)
                            {

                                int iNeigh = i;
                                int jNeigh = j + 1;
                                if (periodicY && jNeigh >= nY)
                                {
                                    (new CellFaceTag()
                                    {
                                        EdgeTag = peryTag,
                                        PeriodicInverse = false,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Edge.Top,

                                        NeighCell_GlobalID = iNeigh * (yNodes.Length - 1) + (0)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }

                                iNeigh = i + 1;
                                jNeigh = j;
                                if (periodicX && iNeigh >= nX)
                                {
                                    (new CellFaceTag()
                                    {
                                        EdgeTag = perxTag,
                                        PeriodicInverse = false,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Edge.Right,
                                        NeighCell_GlobalID = 0 * (yNodes.Length - 1) + (jNeigh)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }

                                iNeigh = i - 1;
                                jNeigh = j;
                                if (periodicX && iNeigh < 0)
                                {
                                    (new CellFaceTag()
                                    {
                                        EdgeTag = perxTag,
                                        PeriodicInverse = true,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Edge.Left,
                                        NeighCell_GlobalID = (nX - 1) * (yNodes.Length - 1) + (jNeigh)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }

                                iNeigh = i;
                                jNeigh = j - 1;
                                if (periodicY && jNeigh < 0)
                                {
                                    (new CellFaceTag()
                                    {
                                        EdgeTag = peryTag,
                                        PeriodicInverse = true,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Edge.Bottom,
                                        NeighCell_GlobalID = iNeigh * (yNodes.Length - 1) + (nY - 1)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }
                            }
                        }
                    }

                    grid.Cells = Cells.ToArray();
                }

                /*
                // boundary cells
                // ==============

                if (myrank == 0) {
                    List<BcCell> BoundaryCells = new List<BcCell>(2*iE + 2*nY);

                    cnt = nX*nY - 1;

                    grid.EdgeTagsNames.Add(1, "Bottom");
                    for (int i = 0; i < nX; i++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int j = 0;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = (i + 1) + j * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 1;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.EdgeTagsNames.Add(2, "Top");
                    for (int i = 0; i < nX; i++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int j = nY;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = (i + 1) + j * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 2;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.EdgeTagsNames.Add(3, "Left");
                    for (int j = 0; j < nY; j++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int i = 0;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = i + (j + 1) * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 3;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.EdgeTagsNames.Add(4, "Right");
                    for (int j = 0; j < nY; j++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int i = nX;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = i + (j + 1) * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 4;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.BcCells = BoundaryCells.ToArray();
                } else {
                    grid.BcCells = null;
                }
                

                // return
                // ======
                /*

                if (CutOuts != null && CutOuts.Length > 0)
                {
                    grid.CompressGlobalID();
                    grid.CompressNodeIndices();
                }

                return grid;
            }
        }

        */

        /// <summary>
        /// constructs a new 2D Grid
        /// </summary>
        /// <param name="xNodes"></param>
        /// <param name="yNodes"></param>
        /// <param name="periodicX">turn on periodic boundary conditions in x direction</param>
        /// <param name="periodicY">turn on periodic boundary conditions in y direction</param>
        /// <param name="CutOuts">Optional regions that are not meshed</param>
        /// <param name="type">The type of the cells to be used</param>
        /// <param name="NonlinearGridTrafo">
        /// Arbitrary transformation (i.e. a diffeomorphism, i.e. bijective and differentiable) applied to <paramref name="xNodes"/> and <paramref name="yNodes"/>, optioanl
        /// </param>
        public static Grid2D Cartesian2DGrid(double[] xNodes, double[] yNodes, 
            CellType type = CellType.Square_Linear, 
            bool periodicX = false, bool periodicY = false, 
            Func<Vector,Vector> NonlinearGridTrafo = null,
            params BoundingBox[] CutOuts) {
            using (var tr = new FuncTrace()) {
                MPICollectiveWatchDog.Watch();

                if(NonlinearGridTrafo != null) {
                    if(type == CellType.Square_Linear)
                        throw new NotSupportedException($"Not recommended to use a nonlinear transformation of the mesh together with linear square cells - use at least {CellType.Square_4}!");
                }

                
                // Some Checks
                // ===========
                CheckMonotonicity(xNodes);
                CheckMonotonicity(yNodes);

                // split along x-Axis (cells can be redistributed by ParMETIS anyway)
                // ==================================================================
                int nX = xNodes.Length - 1;
                int nY = yNodes.Length - 1;

                //if (nX < 3 && periodicX)
                //    throw new ArithmeticException("At least 3 Elements are required for Periodic Boundary Condition to work");
                //if (nY < 3 && periodicY)
                //    throw new ArithmeticException("At least 3 Elements are required for Periodic Boundary Condition to work");
                if (nX < 1)
                    throw new ArithmeticException("At least 2 node are required.");
                if (nY < 1)
                    throw new ArithmeticException("At least 2 node are required.");

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                int i0 = nX * myrank / size;        // 1st x-Index on this proc.
                int iE = nX * (myrank + 1) / size;  // 1st x-Index on next proc. rank

                // Return object 
                // =============

                Grid2D grid;
                using (new BlockTrace("grid_object_instantiation", tr)) {
                    grid = new Grid2D(Square.Instance);
                }

                // define periodic transformations, if necessary
                // =============================================
                byte perxTag = 0;
                byte peryTag = 0;

                if (periodicX) {
                    //if (NonlinearGridTrafo != null)
                    //    throw new NotSupportedException("grid transformation is not supported for periodic domains");
                    Vector[] x = { new Vector { xNodes[0], yNodes[0] }, new Vector { xNodes[0], yNodes[nY] } };
                    Vector[] y = { new Vector { xNodes[nX], yNodes[0] }, new Vector { xNodes[nX], yNodes[nY] } };
                    Vector N = new Vector() { 1.0, 0 };
                    Vector Nx, Ny;

                    if(NonlinearGridTrafo != null) {
                        //Vector xm = 0.5 * (x[0] + x[1]);
                        //Nx = NonlinearGridTrafo(xm + N) - NonlinearGridTrafo(xm);
                        //Nx.Normalize();

                        //Vector ym = 0.5 * (y[0] + y[1]);
                        //Ny = NonlinearGridTrafo(ym + N) - NonlinearGridTrafo(ym);
                        //Ny.Normalize();

                        x[0] = NonlinearGridTrafo(x[0]);
                        x[1] = NonlinearGridTrafo(x[1]);
                        y[0] = NonlinearGridTrafo(y[0]);
                        y[1] = NonlinearGridTrafo(y[1]);

                        var tx = x[1] - x[0];
                        Nx = tx.Rotate2D(-Math.PI / 2);
                        Nx.NormalizeInPlace();
                        
                        var ty = y[1] - y[0];
                        Ny = ty.Rotate2D(-Math.PI / 2);
                        Ny.NormalizeInPlace();

                    } else {
                        Nx = N;
                        Ny = N;
                    }


                    grid.ConstructPeriodicEdgeTrafo(y, Ny, x, Nx, out perxTag);
                    grid.EdgeTagNames.Add(perxTag, "Periodic-X");
                }

                if (periodicY) {
                    if (NonlinearGridTrafo != null)
                        throw new NotSupportedException("grid transformation is not supported for periodic domains");
                    Vector[] x = { new Vector { xNodes[0], yNodes[0] }, new Vector { xNodes[nX], yNodes[0] } };
                    Vector[] y = { new Vector { xNodes[0], yNodes[nY] }, new Vector { xNodes[nX], yNodes[nY] } };
                    Vector N = new Vector { 0, 1.0 };
                    Vector Nx, Ny;

                    if(NonlinearGridTrafo != null) {
                        Vector xm = 0.5 * (x[0] + x[1]);
                        Nx = NonlinearGridTrafo(xm + N) - NonlinearGridTrafo(xm);
                        Nx.NormalizeInPlace();

                        Vector ym = 0.5 * (y[0] + y[1]);
                        Ny = NonlinearGridTrafo(ym + N) - NonlinearGridTrafo(ym);
                        Ny.NormalizeInPlace();

                        x[0] = NonlinearGridTrafo(x[0]);
                        x[1] = NonlinearGridTrafo(x[1]);
                        y[0] = NonlinearGridTrafo(y[0]);
                        y[1] = NonlinearGridTrafo(y[1]);

                    } else {
                        Nx = N;
                        Ny = N;
                    }

                    grid.ConstructPeriodicEdgeTrafo(y, Nx, x, Ny, out peryTag);
                    grid.EdgeTagNames.Add(peryTag, "Periodic-Y");
                }

                // set cells
                // =========
                int LocalNoOfCells = (iE - i0) * (yNodes.Length - 1);
                List<Cell> Cells = new List<Cell>(LocalNoOfCells);

                bool IsInCutOut(int i, int j) {
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

                    if (CutOuts != null) {
                        double xC = 0.5 * (xNodes[i] + xNodes[(i + 1)]);
                        double yC = 0.5 * (yNodes[j] + yNodes[(j + 1)]);

                        return (CutOuts.Any(BB => BB.Contains(xC, yC)));

                    } else {
                        return false;
                    }
                }


                int cnt = 0;
                if (Math.Abs(i0 - iE) <= 0) {
                    //throw new ApplicationException("unable to do partitioning; X-Slice on processor " + myrank + " is empty; Try to load grid from database;");
                    cnt = 0;
                } else {

                    var Kref = grid.RefElements.Single(re => re.GetType() == typeof(Square));
                    if (!Kref.SupportedCellTypes.Contains(type))
                        throw new ArgumentException("unsupported cell type");

                    var RefNodes = Kref.GetInterpolationNodes(type);
                    int NoOfNodes = Kref.GetInterpolationNodes(type).GetLength(0);

                    for (int i = i0; i < iE; i++) {
                        for (int j = 0; j < nY; j++) {

                            // cut-out regions test
                            // ====================

                            if (IsInCutOut(i,j)) {
                                continue;
                            }
                            cnt++;

                            Cell C_cnt = new Cell();
                            Cells.Add(C_cnt);

                            // define cell
                            // ===========
                            C_cnt.GlobalID = i * (yNodes.Length - 1) + j;
                            C_cnt.Type = type;

                            // transformation
                            // ==============
                            C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 2);
                            {
                                NoOfNodes = RefNodes.GetLength(0);
                                Debug.Assert(RefNodes.GetLength(1) == 2);

                                //C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, D);

                                double xL = xNodes[i];
                                double xR = xNodes[i + 1];
                                double yL = yNodes[j];
                                double yR = yNodes[j + 1];

                                for(int iNode = 0; iNode < NoOfNodes; iNode++) {
                                    double xi = 0.5 * (RefNodes[iNode, 0] + 1.0);
                                    double eta = 0.5 * (RefNodes[iNode, 1] + 1.0);

                                    Vector A = new Vector(2);
                                    A.x = xL * (1.0 - xi) + xR * xi;
                                    A.y = yL * (1.0 - eta) + yR * eta;

                                    Vector B = new Vector(2);
                                    if(NonlinearGridTrafo != null) {
                                        B = NonlinearGridTrafo(A);
                                    } else {
                                        B = A;
                                    }

                                    C_cnt.TransformationParams[iNode, 0] = B.x;
                                    C_cnt.TransformationParams[iNode, 1] = B.y;
                                }
                            }

                            // cell neighbourship
                            // ==================
                            C_cnt.NodeIndices = new long[4];
                            C_cnt.NodeIndices[0] = i + j * xNodes.Length;
                            C_cnt.NodeIndices[1] = (i + 1) + j * xNodes.Length;
                            C_cnt.NodeIndices[2] = i + (j + 1) * xNodes.Length;
                            C_cnt.NodeIndices[3] = (i + 1) + (j + 1) * xNodes.Length;

                            // Edge tags (only periodic)
                            // =========================
                            if (periodicY || periodicX) {

                                int iNeigh = i;
                                int jNeigh = j + 1;
                                if (periodicY && jNeigh >= nY && !IsInCutOut(iNeigh, jNeigh)) {
                                    (new CellFaceTag() {
                                        EdgeTag = peryTag,
                                        PeriodicInverse = false,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Faces.Top,

                                        NeighCell_GlobalID = iNeigh * (yNodes.Length - 1) + (0)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }

                                iNeigh = i + 1;
                                jNeigh = j;
                                if (periodicX && iNeigh >= nX && !IsInCutOut(iNeigh, jNeigh)) {
                                    (new CellFaceTag() {
                                        EdgeTag = perxTag,
                                        PeriodicInverse = false,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = 0 * (yNodes.Length - 1) + (jNeigh)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }

                                iNeigh = i - 1;
                                jNeigh = j;
                                if (periodicX && iNeigh < 0 && !IsInCutOut(iNeigh, jNeigh)) {
                                    (new CellFaceTag() {
                                        EdgeTag = perxTag,
                                        PeriodicInverse = true,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Faces.Left,
                                        NeighCell_GlobalID = (nX - 1) * (yNodes.Length - 1) + (jNeigh)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }

                                iNeigh = i;
                                jNeigh = j - 1;
                                if (periodicY && jNeigh < 0 && !IsInCutOut(iNeigh, jNeigh)) {
                                    (new CellFaceTag() {
                                        EdgeTag = peryTag,
                                        PeriodicInverse = true,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Faces.Bottom,
                                        NeighCell_GlobalID = iNeigh * (yNodes.Length - 1) + (nY - 1)
                                    }).AddToArray(ref C_cnt.CellFaceTags);
                                }
                            }
                        }
                    }

                }

                grid.Cells = Cells.ToArray();

                /*
                // boundary cells
                // ==============

                if (myrank == 0) {
                    List<BcCell> BoundaryCells = new List<BcCell>(2*iE + 2*nY);

                    cnt = nX*nY - 1;

                    grid.EdgeTagsNames.Add(1, "Bottom");
                    for (int i = 0; i < nX; i++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int j = 0;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = (i + 1) + j * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 1;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.EdgeTagsNames.Add(2, "Top");
                    for (int i = 0; i < nX; i++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int j = nY;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = (i + 1) + j * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 2;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.EdgeTagsNames.Add(3, "Left");
                    for (int j = 0; j < nY; j++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int i = 0;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = i + (j + 1) * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 3;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.EdgeTagsNames.Add(4, "Right");
                    for (int j = 0; j < nY; j++) {
                        cnt++;
                        BcCell B_cnt = new BcCell();

                        int i = nX;
                        B_cnt.NodeIndices = new int[2];
                        B_cnt.NodeIndices[0] = i + j * xNodes.Length;
                        B_cnt.NodeIndices[1] = i + (j + 1) * xNodes.Length;
                        B_cnt.Type = CellType.Line_2;
                        B_cnt.GlobalID = cnt;

                        B_cnt.EdgeTag = 4;

                        BoundaryCells.Add(B_cnt);
                    }

                    grid.BcCells = BoundaryCells.ToArray();
                } else {
                    grid.BcCells = null;
                }
                 */

                cnt = cnt.MPISum();
                if(cnt <= 0) {
                    throw new ArgumentException("Grid is empty - check arguments (cut-outs).");
                }

                // return
                // ======

                if (CutOuts != null && CutOuts.Length > 0) {
                    grid.CompressGlobalID();
                    grid.CompressNodeIndices();
                }
  
                foreach (var cl in grid.Cells) {
                    if (cl.GlobalID < 0)
                        throw new ApplicationException("Internal error - illegal GlobalID.");
                    if (cl.GlobalID >= cnt)
                        throw new ApplicationException("Internal error - illegal GlobalID.");
                }

                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                return grid;
            }
        }

        /// <summary>
        /// Constructs a new 2D Grid.
        /// </summary>
        public static Grid2D Trapezoidal2dGrid(double widthBottom, double widthTop, int NoOfNodesX, double[] yNodes, CellType type = CellType.Square_4) {
            using (var tr = new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                if (widthBottom <= 0)
                    throw new ArgumentOutOfRangeException();
                if (widthTop <= 0)
                    throw new ArgumentOutOfRangeException();
                if (NoOfNodesX < 2)
                    throw new ArgumentOutOfRangeException();
                CheckMonotonicity(yNodes);
                Grid2D grid;
                using (new BlockTrace("grid_object_instantiation", tr)) {
                    grid = new Grid2D(Square.Instance);
                }

                // split along x-Axis (cells can be redistributed by ParMETIS anyway)
                // ==================================================================
                int nX = NoOfNodesX - 1;
                int nY = yNodes.Length - 1;

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                int i0 = nX * myrank / size;        // 1st x-Index on this proc.
                int iE = nX * (myrank + 1) / size;  // 1st x-Index on next proc. rank

                if (Math.Abs(i0 - iE) <= 0) {
                    //throw new ApplicationException("unable to do partitioning; X-Slice on processor " + myrank + " is empty; Try to load grid from database;");
                    grid.Cells = new Cell[0];
                }
                else {

                    int D = grid.SpatialDimension;
                    var Kref = grid.RefElements.Single(re => re.GetType() == typeof(Square));
                    if (!Kref.SupportedCellTypes.Contains(type))
                        throw new ArgumentException("unsupported cell type");


                    // set cells
                    // =========
                    int LocalNoOfCells = (iE - i0) * (yNodes.Length - 1);
                    var RefNodes = Kref.GetInterpolationNodes(type);
                    int NoOfNodes = Kref.GetInterpolationNodes(type).GetLength(0);

                    List<Cell> Cells = new List<Cell>(LocalNoOfCells);

                    int cnt = -1;
                    //int Cellcounter = 0;
                    for (int i = i0; i < iE; i++) {
                        for (int j = 0; j < nY; j++) {
                            cnt++;


                            Cell C_cnt = new Cell();
                            Cells.Add(C_cnt);

                            // define cell
                            // ===========
                            C_cnt.GlobalID = i * (yNodes.Length - 1) + j;
                            C_cnt.Type = type;

                            // transformation
                            // ==============
                            C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, D);
                            {
                                NoOfNodes = RefNodes.GetLength(0);
                                Debug.Assert(RefNodes.GetLength(1) == D);

                                //C_cnt.TransformationParams = MultidimensionalArray.Create(NoOfNodes, D);

                                //double xL_B = xNodes[i];
                                //double xR_B = xNodes[i + 1];
                                //double yL = yNodes[j];
                                //double yR = yNodes[j + 1];

                                double etaBot = (yNodes[j] - yNodes[0]) / (yNodes[nY] - yNodes[0]);
                                double etaTop = (yNodes[j + 1] - yNodes[0]) / (yNodes[nY] - yNodes[0]);
                                double wBot = widthBottom + etaBot * (widthTop - widthBottom);
                                double wTop = widthBottom + etaTop * (widthTop - widthBottom);

                                double xiL = (double)i / (double)nX;
                                double xiR = (double)(i + 1) / (double)nX;
                                double x0 = (-wBot * 0.5) * (1.0 - xiL) + (wBot * 0.5) * xiL;
                                double x1 = (-wBot * 0.5) * (1.0 - xiR) + (wBot * 0.5) * xiR;
                                double x3 = (-wTop * 0.5) * (1.0 - xiL) + (wTop * 0.5) * xiL;
                                double x4 = (-wTop * 0.5) * (1.0 - xiR) + (wTop * 0.5) * xiR;

                                for (int iNode = 0; iNode < NoOfNodes; iNode++) {
                                    double xi = 0.5 * (RefNodes[iNode, 0] + 1.0);
                                    double eta = 0.5 * (RefNodes[iNode, 1] + 1.0);

                                    double yN = yNodes[j] + eta * (yNodes[j + 1] - yNodes[j]);
                                    double xL = x0 + eta * (x3 - x0);
                                    double xR = x1 + eta * (x4 - x1);

                                    C_cnt.TransformationParams[iNode, 0] = xL * (1.0 - xi) + xR * xi;
                                    C_cnt.TransformationParams[iNode, 1] = yN;
                                }
                            }

                            // cell neighbourship
                            // ==================
                            C_cnt.NodeIndices = new long[4];
                            C_cnt.NodeIndices[0] = i + j * NoOfNodesX;
                            C_cnt.NodeIndices[1] = (i + 1) + j * NoOfNodesX;
                            C_cnt.NodeIndices[2] = i + (j + 1) * NoOfNodesX;
                            C_cnt.NodeIndices[3] = (i + 1) + (j + 1) * NoOfNodesX;

                        }
                    }

                    grid.Cells = Cells.ToArray();
                }

                // return
                // ======
                return grid;
            }
        }



        private static void Rotate(MultidimensionalArray Points, double Rotation) {
            MultidimensionalArray RotMatrix = MultidimensionalArray.Create(2, 2);
            RotMatrix[0, 0] = Math.Cos(Rotation);
            RotMatrix[0, 1] = -Math.Sin(Rotation);
            RotMatrix[1, 0] = Math.Sin(Rotation);
            RotMatrix[1, 1] = Math.Cos(Rotation);

            //int NPoints = Points.GetLength(0);

            var res = Points.CloneAs();

            res.Multiply(1.0, RotMatrix, Points, 0.0, "nd", "dl", "nl");

            Points.Clear();
            Points.Acc(1.0, res);
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="Nodes"></param>
        public delegate void PointsFunc(MultidimensionalArray Nodes);

        /// <summary>
        /// standard constructor;
        /// performs basic checks on array dimensions;
        /// </summary>
        /// <remarks>
        /// Grid creation is MPI-serial, i.e. on exit, all cells are located at MPI process rank 0.
        /// During setup, the grid must be redistributed, e.g. by ParMETIS, see <see cref="GridCommons.Redistribute"/>.
        /// This will be usually done automatically.
        /// </remarks>
        static public Grid2D BilinearSquareGrid(double[] xNodes, double[] yNodes, double factor = 0.8, double threshold = 0.0, double Rotation = 0.0, PointsFunc pf = null) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Square.Instance);

                CheckMonotonicity(xNodes);
                CheckMonotonicity(yNodes);

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);




                if (myrank == 0) {
                    //// create cell
                    int J = (xNodes.Length - 1) * (yNodes.Length - 1); // total number of cells
                    grid.Cells = new Cell[J];

                    int j = 0;//cell index
                    if (J == 1) {
                        Cell Cj0 = new Cell();
                        Cj0.GlobalID = j;
                        Cj0.Type = CellType.Square_4;
                        Cj0.TransformationParams = MultidimensionalArray.Create(2, 4);
                        var Bild0 = Cj0.TransformationParams;

                        Bild0[0, 0] = -1;
                        Bild0[1, 0] = 4;
                        Bild0[2, 0] = -1;
                        Bild0[3, 0] = 4;
                        Bild0[0, 1] = 2;
                        Bild0[1, 1] = 2;
                        Bild0[2, 1] = 4;
                        Bild0[3, 1] = 4;

                        grid.Cells[0] = Cj0;
                        grid.Cells[0].NodeIndices = new long[] { 0, 1, 2, 3 };

                    } else {


                        for (int i = 0; i < (xNodes.Length - 1); i++) {
                            for (int k = 0; k < (yNodes.Length - 1); k++) {
                                Cell Cj0 = new Cell();
                                Cj0.GlobalID = j;
                                Cj0.Type = CellType.Square_4;
                                Cj0.TransformationParams = MultidimensionalArray.Create(4, 2);
                                var Bild0 = Cj0.TransformationParams;

                                double dx = xNodes[i + 1] - xNodes[i];
                                double KK = dx / (yNodes.Length);


                                if (i == 0) {
                                    Bild0[0, 0] = xNodes[i];
                                    Bild0[1, 0] = xNodes[i + 1] - k * KK * factor;
                                    Bild0[2, 0] = xNodes[i];
                                    Bild0[3, 0] = xNodes[i + 1] - (k + 1) * KK * factor;
                                    Bild0[0, 1] = yNodes[k];
                                    Bild0[1, 1] = yNodes[k];
                                    Bild0[2, 1] = yNodes[k + 1];
                                    Bild0[3, 1] = yNodes[k + 1];
                                } else if (i == xNodes.Length - 2) {
                                    Bild0[0, 0] = (xNodes[i] > threshold) ? (xNodes[i]) : (xNodes[i] - k * KK * factor);
                                    Bild0[1, 0] = xNodes[i + 1];
                                    Bild0[2, 0] = (xNodes[i] > threshold) ? (xNodes[i]) : (xNodes[i] - (k + 1) * KK * factor);
                                    Bild0[3, 0] = xNodes[i + 1];
                                    Bild0[0, 1] = yNodes[k];
                                    Bild0[1, 1] = yNodes[k];
                                    Bild0[2, 1] = yNodes[k + 1];
                                    Bild0[3, 1] = yNodes[k + 1];
                                } else {
                                    Bild0[0, 0] = (xNodes[i] > threshold) ? (xNodes[i]) : (xNodes[i] - k * KK * factor);
                                    Bild0[1, 0] = (xNodes[i + 1] > threshold) ? (xNodes[i + 1]) : (xNodes[i + 1] - k * KK * factor);
                                    Bild0[2, 0] = (xNodes[i] > threshold) ? (xNodes[i]) : (xNodes[i] - (k + 1) * KK * factor);
                                    Bild0[3, 0] = (xNodes[i + 1] > threshold) ? (xNodes[i + 1]) : (xNodes[i + 1] - (k + 1) * KK * factor);
                                    Bild0[0, 1] = yNodes[k];
                                    Bild0[1, 1] = yNodes[k];
                                    Bild0[2, 1] = yNodes[k + 1];
                                    Bild0[3, 1] = yNodes[k + 1];
                                }

                                if (!(Bild0[0, 0] < Bild0[1, 0]))
                                    throw new ArgumentException("stretching factor chosen to large");
                                if (!(Bild0[2, 0] < Bild0[3, 0]))
                                    throw new ArgumentException("stretching factor chosen to large");
                                Rotate(Bild0, Rotation);

                                if (pf != null)
                                    pf(Bild0);

                                Cj0.NodeIndices = new long[] {
                                    k*xNodes.Length + i,          // unten, links
                                    k*xNodes.Length + (i + 1),    // unten, rechts
                                    (k + 1)*xNodes.Length + i,    // oben, links
                                    (k + 1)*xNodes.Length + i + 1 // oben, rechts
                                };



                                grid.Cells[j] = Cj0;
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
        /// This Method convert Points in the Parametric r-s-Space into Points in the x-y-Space
        /// </summary>
        /// <param name="rPoint"></param>
        /// <param name="sPoint"></param>
        /// <returns></returns>
        internal static Vector Param2XY(double rPoint, double sPoint) {

            // Calculate Coordinates
            Vector xyPoint = new Vector(2);
            // x-Coordinate
            xyPoint[0] = rPoint * Math.Cos(2 * Math.PI * sPoint);
            // y-Coordinate
            xyPoint[1] = rPoint * Math.Sin(2 * Math.PI * sPoint);
            return xyPoint;
        }

        /// <summary>
        /// Generates a Grid Consisting of a Ring. from Coordinates in
        /// Parametric Description<see cref="Grid2D.Param2XY"/>. Later
        /// Arbitrary Transformation Functions might be introduced.
        /// </summary>
        /// <param name="rNodes">
        /// Vector consisting of Points in radial direction, must be greater than 0.
        /// </param>
        /// <param name="sNodes">
        /// Vector consisting of Points in rotational direction (1.0 represents 2Pi)
        /// </param>
        /// <param name="CellType">
        /// Type of Elements. Currently only Quads supported.
        /// </param>
        /// <param name="PeriodicS">
        /// Toggle for periodicity in s-direction
        /// </param>
        /// <returns>
        /// A grid consisting of curved elements forming a (section of) a ring
        /// </returns>
        public static Grid2D CurvedSquareGrid(double[] rNodes, double[] sNodes, CellType CellType, bool PeriodicS = true) {
            using (new FuncTrace()) {
                if (!(Square.Instance).SupportedCellTypes.Contains(CellType))
                    throw new ArgumentOutOfRangeException("illegal cell type.");
                if (CellType == CellType.Square_Linear)
                    throw new ArgumentOutOfRangeException("not supported for linear squares.");

                int nR = rNodes.Length - 1;
                int nS = sNodes.Length - 1;

                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Square.Instance);

                CheckMonotonicity(rNodes);
                CheckMonotonicity(sNodes);

                //Exceptions for bounds of r and s for a segment of a circle
                if ((Math.Abs(sNodes.First() - sNodes.Last()) > 1)) {
                    throw new ArgumentException("sNodes must not exceed an interval width of 1");
                };
                if (rNodes.First() <= 0) {
                    throw new ArgumentException("rNodes must be r>0");
                };
                bool fullcircle = false;
                if ((Math.Abs(sNodes.First() - sNodes.Last()) == 1)) {
                    fullcircle = true;
                }

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);


                // define periodic transformations
                // =============================================
                byte persTag = 0;
                if (PeriodicS == true) {
                    //Periodic Boundary Inlet
                    Vector[] PerBoundIn = { Param2XY(rNodes.First(), sNodes.First()), Param2XY(rNodes.Last(), sNodes.First()) };
                    // Vector Connecting the two Points of the inlet
                    Vector PerBoundInCon = new Vector(2);
                    PerBoundInCon.x = PerBoundIn[1][0] - PerBoundIn[0][0];//(Param2XY(rNodes.First(), sNodes.First())[0] - Param2XY(rNodes.Last(), sNodes.First())[0]);
                    PerBoundInCon.y = PerBoundIn[1][1] - PerBoundIn[0][1];//(Param2XY(rNodes.First(), sNodes.First())[1] - Param2XY(rNodes.Last(), sNodes.First())[1]);
                    // Normal onto Inlet Pointing outwards
                    PerBoundInCon.NormalizeInPlace();
                    Vector PerBoundInNormal = new Vector(-PerBoundInCon.y, +PerBoundInCon.x );
                    //Periodic Boundary Inlet
                    Vector[] PerBoundOut = { Param2XY(rNodes.First(), sNodes.Last()), Param2XY(rNodes.Last(), sNodes.Last()) };
                    Vector PerBoundOutCon = new Vector(2);
                    PerBoundOutCon.x = PerBoundOut[1][0] - PerBoundOut[0][0];//(Param2XY(rNodes.First(), sNodes.Last())[0] - Param2XY(rNodes.Last(), sNodes.Last())[0]);
                    PerBoundOutCon.y = PerBoundOut[1][1] - PerBoundOut[0][1];//(Param2XY(rNodes.First(), sNodes.Last())[1] - Param2XY(rNodes.Last(), sNodes.Last())[1]);
                    PerBoundOutCon.NormalizeInPlace();
                    double[] PerBoundOutNormal = { -PerBoundOutCon.y, +PerBoundOutCon.x };
                    grid.ConstructPeriodicEdgeTrafo(PerBoundIn, PerBoundInNormal, PerBoundOut, PerBoundOutNormal, out persTag);
                    grid.EdgeTagNames.Add(persTag, "Periodic-S");
                }

                var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Square));
                MultidimensionalArray InterpolationNodes = Kref.GetInterpolationNodes(CellType);
                // Point in Reference Element: InterpolationNodes[NodeIndex,SpatialDimension]


                if (myrank == 0) {
                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++

                    // create cell
                    int J = nR * nS; // total number of cells
                    grid.Cells = new Cell[J];

                    int j = 0;//cell index

                    int NoOfNodes = Kref.GetInterpolationNodes(CellType).GetLength(0);

                    for (int i = 0; i < (nR); i++) {
                        for (int k = 0; k < (nS); k++) {
                            Cell Cj0 = new Cell();
                            Cj0.GlobalID = j;
                            Cj0.Type = CellType;

                            Cj0.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 2);
                            //  var Bild0 = Cj0.TransformationParams;




                            for (int PointNumber = 0; PointNumber < NoOfNodes; PointNumber++) {

                                // Interpolate in rs-Domain according to these Coordinates

                                Vector rsPoint = new Vector(2);
                                rsPoint[0] = rNodes[i] + (rNodes[i + 1] - rNodes[i]) * 0.5 * (InterpolationNodes[PointNumber, 0] + 1);
                                rsPoint[1] = sNodes[k] + (sNodes[k + 1] - sNodes[k]) * 0.5 * (InterpolationNodes[PointNumber, 1] + 1);


                                // Convert Point into physical Coordinates

                                double[] xyPoint = Param2XY(rsPoint[0], rsPoint[1]);
                                // Write Physical Coordinates to TransformParams
                                for (int dim = 0; dim < 2; dim++) {
                                    Cj0.TransformationParams[PointNumber, dim] = xyPoint[dim];
                                }
                            }


                            //if (!(Bild0[0, 0] < Bild0[1, 0]))
                            //    throw new ArgumentException("Elements Degenerated - Aborting");
                            //if (!(Bild0[2, 0] < Bild0[3, 0]))
                            //    throw new ArgumentException("Elements Degenerated- Aborting");


                            // ------------------ //
                            // Cell Neigbourships //
                            // ------------------ //

                            if (fullcircle == false || k != nS - 1) {
                                Cj0.NodeIndices = new long[] {
                                    k*rNodes.Length + i,          // unten, links
                                    k*rNodes.Length + (i + 1),    // unten, rechts
                                    (k + 1)*rNodes.Length + i,    // oben, links
                                    (k + 1)*rNodes.Length + i + 1 // oben, rechts
                               };
                            }
                            else {
                                Cj0.NodeIndices = new long[] {
                                        k*rNodes.Length + i,          // unten, links
                                        k*rNodes.Length + (i + 1),    // unten, rechts
                                        i,    // oben, links //  Periodic BC
                                        i + 1 // oben, rechts //  Periodic BC
                                    };
                            }


                            // ----------------------------- //
                            //  Periodic Boundary Conditions
                            // ----------------------------- //

                            if (PeriodicS == true && fullcircle == false) {
                                if (k == nS - 1) {
                                    (new CellFaceTag() {
                                        EdgeTag = persTag,
                                        PeriodicInverse = true,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Faces.Top,
                                        NeighCell_GlobalID = j - (nS - 1)//(i) * (nS)
                                    }).AddToArray(ref Cj0.CellFaceTags);
                                }

                                //iNeigh = i; kNeigh = k - 1;
                                if (k == 0) {
                                    (new CellFaceTag() {
                                        EdgeTag = persTag,
                                        PeriodicInverse = false,
                                        ConformalNeighborship = true,
                                        FaceIndex = (int)Square.Faces.Bottom,
                                        //NeighCell_GlobalID = (iNeigh) * (nS) + nR-1
                                        NeighCell_GlobalID = j + (nS - 1) //i * nS + nS - 1 //(iNeigh) * (nS) + nR - 1
                                    }).AddToArray(ref Cj0.CellFaceTags);
                                }
                            }





                            grid.Cells[j] = Cj0;
                            j++;
                        }
                    }

                }
                else {
                    grid.Cells = new Cell[0];
                }
                return grid;

            }
        }

        /// <summary>
        /// Generates a Grid Consisting of a curved Surface from Coordinates in
        /// Parametric Description<see cref="Grid2D.Param2XY"/>. Later
        /// Arbitrary Transformation Functions might be introduced.
        /// </summary>
        /// <param name="rNodes">
        /// Vector consisting of Points in r-Direction
        /// </param>
        /// <param name="sNodes">
        /// Vector consisting of Points in s-Direction
        /// </param>
        /// <param name="CellType">
        /// Type of Elements. Currently only Quads supported.
        /// </param>
        /// <param name="PeriodicR">
        /// Toggle for periodicity in r-direction
        /// </param>
        /// <param name="Topology">
        /// Convert Points in the Parametric r-s-Space (arguments) into Points in the x-y-Space (return value)
        /// </param>
        /// <returns>
        /// A grid consisting of curved elements forming a curved Surface
        /// </returns>
        public static Grid2D CurvedSquareGridChannel(double[] rNodes, double[] sNodes, CellType CellType, bool PeriodicR = true, Func<double,double,Vector> Topology = null) {
            using (new FuncTrace()) {
                if (!(Square.Instance).SupportedCellTypes.Contains(CellType))
                    throw new ArgumentOutOfRangeException("illegal cell type.");
                if (CellType == CellType.Square_Linear)
                    throw new ArgumentOutOfRangeException("not supported for linear squares.");

                int nR = rNodes.Length - 1;
                int nS = sNodes.Length - 1;

                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Square.Instance);

                CheckMonotonicity(rNodes);
                CheckMonotonicity(sNodes);

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                // define periodic transformations
                // =============================================
                byte persTag = 0;
                if (PeriodicR == true) {
                    //Periodic Boundary Inlet
                    Vector[] PerBoundIn = { Topology(rNodes.First(), sNodes.First()), Topology(rNodes.First(), sNodes.Last()) };
                    // Vector Connecting the two Points of the inlet
                    Vector PerBoundInCon = new Vector(2);
                    PerBoundInCon.x = PerBoundIn[1][0] - PerBoundIn[0][0]; //(Topology(rNodes.First(), sNodes.First())[0] - Topology(rNodes.Last(), sNodes.First())[0]);
                    PerBoundInCon.y = PerBoundIn[1][1] - PerBoundIn[0][1]; //(Topology(rNodes.First(), sNodes.First())[1] - Topology(rNodes.Last(), sNodes.First())[1]);
                    // Normal onto Inlet Pointing outwards
                    PerBoundInCon.NormalizeInPlace();
                    Vector PerBoundInNormal = new Vector(-PerBoundInCon.y, +PerBoundInCon.x);
                    //Periodic Boundary Inlet
                    Vector[] PerBoundOut = { Topology(rNodes.Last(), sNodes.First()), Topology(rNodes.Last(), sNodes.Last()) };
                    Vector PerBoundOutCon = new Vector(2);
                    PerBoundOutCon.x = PerBoundOut[1][0] - PerBoundOut[0][0];//(Topology(rNodes.First(), sNodes.Last())[0] - Topology(rNodes.Last(), sNodes.Last())[0]);
                    PerBoundOutCon.y = PerBoundOut[1][1] - PerBoundOut[0][1];//(Topology(rNodes.First(), sNodes.Last())[1] - Topology(rNodes.Last(), sNodes.Last())[1]);
                    PerBoundOutCon.NormalizeInPlace();
                    Vector PerBoundOutNormal = new Vector( -PerBoundOutCon.y, +PerBoundOutCon.x );
                    grid.ConstructPeriodicEdgeTrafo(PerBoundIn, PerBoundInNormal, PerBoundOut, PerBoundOutNormal, out persTag);
                    grid.EdgeTagNames.Add(persTag, "Periodic-R");
                }

                var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Square));
                MultidimensionalArray InterpolationNodes = Kref.GetInterpolationNodes(CellType);
                // Point in Reference Element: InterpolationNodes[NodeIndex,SpatialDimension]


                if (myrank == 0) {
                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++

                    // create cell
                    int J = nR * nS; // total number of cells
                    grid.Cells = new Cell[J];

                    int j = 0;//cell index

                    int NoOfNodes = Kref.GetInterpolationNodes(CellType).GetLength(0);

                    for (int i = 0; i < (nR); i++) {
                        for (int k = 0; k < (nS); k++) {
                            Cell Cj0 = new Cell();
                            Cj0.GlobalID = j;
                            Cj0.Type = CellType;

                            Cj0.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 2);
                            //  var Bild0 = Cj0.TransformationParams;




                            for (int PointNumber = 0; PointNumber < NoOfNodes; PointNumber++) {

                                // Interpolate in rs-Domain according to these Coordinates

                                Vector rsPoint = new Vector(2);
                                rsPoint[0] = rNodes[i] + (rNodes[i + 1] - rNodes[i]) * 0.5 * (InterpolationNodes[PointNumber, 0] + 1);
                                rsPoint[1] = sNodes[k] + (sNodes[k + 1] - sNodes[k]) * 0.5 * (InterpolationNodes[PointNumber, 1] + 1);

                                // Convert Point into physical Coordinates

                                double[] xyPoint = Topology(rsPoint[0], rsPoint[1]);
                                // Write Physical Coordinates to TransformParams
                                for (int dim = 0; dim < 2; dim++) {
                                    Cj0.TransformationParams[PointNumber, dim] = xyPoint[dim];
                                }
                            }

                            // ------------------ //
                            // Cell Neigbourships //
                            // ------------------ //

                            if (i != nR) {
                                Cj0.NodeIndices = new long[] {
                                    k*rNodes.Length + i,          // unten, links
                                    k*rNodes.Length + (i + 1),    // unten, rechts
                                    (k + 1)*rNodes.Length + i,    // oben, links
                                    (k + 1)*rNodes.Length + i + 1 // oben, rechts
                              };
                            }


                            // ----------------------------- //
                            //  Periodic Boundary Conditions
                            // ----------------------------- //

                            if (PeriodicR == true) {
                                if (i == nR - 1) {
                                    (new CellFaceTag() {
                                        EdgeTag = persTag,
                                        PeriodicInverse = true,
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = j - (nR - 1) * nS
                                    }).AddToArray(ref Cj0.CellFaceTags);
                                }

                                //iNeigh = i; kNeigh = k - 1;
                                if (i == 0) {
                                    (new CellFaceTag() {
                                        EdgeTag = persTag,
                                        PeriodicInverse = false,
                                        FaceIndex = (int)Square.Faces.Left,
                                        NeighCell_GlobalID = j + (nR - 1) * nS
                                    }).AddToArray(ref Cj0.CellFaceTags);
                                }
                            }

                            grid.Cells[j] = Cj0;
                            j++;
                        }
                    }

                }
                else {
                    grid.Cells = new Cell[0];
                }
                return grid;

            }
        }

        /// <summary>
        /// creates an unstructured triangle grid on a Cartesian mesh.
        /// </summary>
        static public Grid2D UnstructuredTriangleGrid(double[] xNodes, double[] yNodes, double JitterScale = 0.0) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Triangle.Instance);

                CheckMonotonicity(xNodes);
                CheckMonotonicity(yNodes);

                if (xNodes.Length < 2)
                    throw new ArgumentException("xNodes must have at least 2 entries.");
                if (yNodes.Length < 2)
                    throw new ArgumentException("yNodes must have at least 2 entries.");

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (myrank == 0) {
                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++

                    //var urBild = this.GridSimplex(0).Vertices;
                    // create cell
                    int J = (xNodes.Length - 1) * (yNodes.Length - 1) * 2; // total number of cells


                    grid.Cells = new Cell[J];


                    // create jitter for the unstructured grid
                    // ---------------------------------------
                    double[,] xDelta, yDelta;
                    {

                        Random rnd = new Random(5678);
                        xDelta = new double[xNodes.Length, yNodes.Length];
                        yDelta = new double[xNodes.Length, yNodes.Length];
                        for (int i = 0; i < xNodes.Length; i++) {
                            double dxPlus = double.MaxValue, dxMinus = double.MaxValue;
                            if (i > 0)
                                dxMinus = xNodes[i] - xNodes[i - 1];
                            if (i < xNodes.Length - 1)
                                dxPlus = xNodes[i + 1] - xNodes[i];
                            double dxMin = Math.Min(dxPlus, dxMinus);


                            for (int k = 0; k < yNodes.Length; k++) {
                                double dyPlus = double.MaxValue, dyMinus = double.MaxValue;
                                if (k > 0)
                                    dyMinus = yNodes[k] - yNodes[k - 1];
                                if (k < yNodes.Length - 1)
                                    dyPlus = yNodes[k + 1] - yNodes[k];
                                double dyMin = Math.Min(dyPlus, dyMinus);

                                xDelta[i, k] = (rnd.NextDouble() - 0.5) * dxMin * JitterScale;
                                yDelta[i, k] = (rnd.NextDouble() - 0.5) * dyMin * JitterScale;
                            }
                        }
                    }


                    int j = 0;
                    for (int i = 0; i < (xNodes.Length - 1); i++) {
                        for (int k = 0; k < (yNodes.Length - 1); k++) {
                            // lower Cell
                            {
                                Cell Cj0 = new Cell();

                                Cj0.GlobalID = j;
                                Cj0.Type = CellType.Triangle_3;
                                Cj0.TransformationParams = MultidimensionalArray.Create(3, 2);

                                Cj0.TransformationParams[0, 0] = xNodes[i] + xDelta[i, k];
                                Cj0.TransformationParams[0, 1] = yNodes[k] + yDelta[i, k];
                                Cj0.TransformationParams[1, 0] = xNodes[i + 1] + xDelta[i + 1, k];
                                Cj0.TransformationParams[1, 1] = yNodes[k] + yDelta[i + 1, k];
                                Cj0.TransformationParams[2, 0] = xNodes[i + 1] + xDelta[i + 1, k + 1];
                                Cj0.TransformationParams[2, 1] = yNodes[k + 1] + yDelta[i + 1, k + 1];

                                Cj0.NodeIndices = new long[] {
                                    k*xNodes.Length + i,
                                    k*xNodes.Length + i + 1,
                                    (k + 1)*xNodes.Length + i + 1
                                };

                                grid.Cells[j] = Cj0;
                                j++;
                            }

                            // upper cell
                            {
                                Cell Cj1 = new Cell();
                                Cj1.GlobalID = j;
                                Cj1.Type = CellType.Triangle_3;
                                Cj1.TransformationParams = MultidimensionalArray.Create(3, 2);

                                Cj1.TransformationParams[0, 0] = xNodes[i] + xDelta[i, k];
                                Cj1.TransformationParams[0, 1] = yNodes[k] + yDelta[i, k];
                                Cj1.TransformationParams[1, 0] = xNodes[i + 1] + xDelta[i + 1, k + 1];
                                Cj1.TransformationParams[1, 1] = yNodes[k + 1] + yDelta[i + 1, k + 1];
                                Cj1.TransformationParams[2, 0] = xNodes[i] + xDelta[i, k + 1];
                                Cj1.TransformationParams[2, 1] = yNodes[k + 1] + yDelta[i, k + 1];

                                Cj1.NodeIndices = new long[] {
                                    k*xNodes.Length + i,
                                    (k + 1)*xNodes.Length + i + 1,
                                    (k + 1)*xNodes.Length + i
                                };

                                grid.Cells[j] = Cj1;
                                j++;
                            }
                        }
                    }

                }
                else {
                    // all MPI processes with rank > 0
                    // return empty grid (use Redistribution)
                    // ++++++++++++++++++++++++++++++++++++++++

                    grid.Cells = new Cell[0];
                }


                return grid;
            }
        }


        static public Grid2D AcuteCornerTriangleGrid(double baseL, double theta, int lvl, bool inlet = false) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Triangle.Instance);

                if ((theta <= 0) || (theta > Math.PI / 2.0))
                    throw new ArgumentOutOfRangeException();

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (myrank == 0) {
                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++

                    int J = (int)Math.Pow(Math.Pow(2, lvl), 2);
                    grid.Cells = new Cell[J];

                    int stages = (int)Math.Pow(2, lvl);
                    double stageL = baseL / stages;
                    double height = baseL * Math.Sin(theta);
                    double[] yNodes = GenericBlas.Linspace(-height, 0, stages + 1);
                    yNodes.ScaleV(-1.0);
                    double heightX = -baseL * Math.Cos(theta);

                    int j = 0;
                    int nIdx1 = 0; int nIdx2 = 0;
                    for (int s = 0; s < stages; s++) {

                        double x1 = heightX + (s * stageL * Math.Cos(theta));
                        double[] xNodes1 = (s == 0) ? new double[] { heightX } : GenericBlas.Linspace(x1-(s*stageL), x1, s+1);
                        nIdx1 += xNodes1.Length - 1;
                        double x2 = heightX + ((s + 1) * stageL * Math.Cos(theta));
                        double[] xNodes2 = GenericBlas.Linspace(x2 - ((s+1) * stageL), x2, s + 2);
                        nIdx2 += xNodes2.Length - 1;

                        int stageJ = xNodes1.Length + xNodes2.Length - 2;
                        for (int sj = 0; sj < stageJ; sj++) {

                            Cell Cj0 = new Cell();

                            Cj0.GlobalID = j;
                            Cj0.Type = CellType.Triangle_3;
                            Cj0.TransformationParams = MultidimensionalArray.Create(3, 2);

                            int sjIdx = sj / 2;
                            if ((sj % 2) == 0) {

                                Cj0.TransformationParams[0, 0] = xNodes1[sjIdx];
                                Cj0.TransformationParams[0, 1] = yNodes[s];
                                Cj0.TransformationParams[1, 0] = xNodes2[sjIdx];
                                Cj0.TransformationParams[1, 1] = yNodes[s + 1];
                                Cj0.TransformationParams[2, 0] = xNodes2[sjIdx + 1];
                                Cj0.TransformationParams[2, 1] = yNodes[s + 1];

                                Cj0.NodeIndices = new long[] {
                                        nIdx1 + sjIdx,
                                        nIdx2 + sjIdx,
                                        nIdx2 + sjIdx + 1 };

                            } else {

                                Cj0.TransformationParams[0, 0] = xNodes1[sjIdx];
                                Cj0.TransformationParams[0, 1] = yNodes[s];
                                Cj0.TransformationParams[1, 0] = xNodes2[sjIdx + 1];
                                Cj0.TransformationParams[1, 1] = yNodes[s + 1];
                                Cj0.TransformationParams[2, 0] = xNodes1[sjIdx + 1];
                                Cj0.TransformationParams[2, 1] = yNodes[s];

                                Cj0.NodeIndices = new long[] {
                                        nIdx1 + sjIdx,
                                        nIdx2 + sjIdx + 1,
                                        nIdx1 + sjIdx + 1 };

                            }

                            grid.Cells[j] = Cj0;
                            j++;
                        }

                    }

                } else {
                    // all MPI processes with rank > 0
                    // return empty grid (use Redistribution)
                    // ++++++++++++++++++++++++++++++++++++++++

                    grid.Cells = new Cell[0];
                }

                return grid;
            }
        }


        static public Grid2D ObtuseCornerTriangleGrid(double baseL, double theta, double lvl) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Triangle.Instance);

                if ((theta < Math.PI / 2.0) || (theta > Math.PI))
                    throw new ArgumentOutOfRangeException();

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (myrank == 0) {
                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++

                    int J = (int)Math.Pow(Math.Pow(2, lvl), 2) * 2;
                    grid.Cells = new Cell[J];

                    int stages = (int)Math.Pow(2, lvl);
                    double stageL = baseL / stages;
                    double height = baseL * Math.Cos(theta - (Math.PI / 2.0));
                    double[] yNodes = GenericBlas.Linspace(-height, 0, stages + 1);
                    yNodes.ScaleV(-1.0);
                    double heightX = baseL * Math.Sin(theta - (Math.PI / 2.0));

                    int j = 0;
                    int nIdx1 = -(stages + 1); int nIdx2 = 0;
                    for (int s = 0; s < stages; s++) {

                        double x1 = heightX - (s * stageL * Math.Sin(theta - (Math.PI / 2.0)));
                        double[] xNodes1 = GenericBlas.Linspace(x1 - baseL, x1, stages + 1);
                        nIdx1 += xNodes1.Length;
                        double x2 = heightX - ((s + 1) * stageL * Math.Sin(theta - (Math.PI / 2.0)));
                        double[] xNodes2 = GenericBlas.Linspace(x2 - baseL, x2, stages + 1);
                        nIdx2 += xNodes2.Length;

                        int stageJ = xNodes1.Length + xNodes2.Length - 2;
                        for (int sj = 0; sj < stageJ; sj++) {

                            Cell Cj0 = new Cell();

                            Cj0.GlobalID = j;
                            Cj0.Type = CellType.Triangle_3;
                            Cj0.TransformationParams = MultidimensionalArray.Create(3, 2);

                            int sjIdx = sj / 2;
                            if ((sj % 2) == 0) {

                                Cj0.TransformationParams[0, 0] = xNodes1[sjIdx];
                                Cj0.TransformationParams[0, 1] = yNodes[s];
                                Cj0.TransformationParams[1, 0] = xNodes2[sjIdx];
                                Cj0.TransformationParams[1, 1] = yNodes[s + 1];
                                Cj0.TransformationParams[2, 0] = xNodes2[sjIdx + 1];
                                Cj0.TransformationParams[2, 1] = yNodes[s + 1];

                                Cj0.NodeIndices = new long[] {
                                        nIdx1 + sjIdx,
                                        nIdx2 + sjIdx,
                                        nIdx2 + sjIdx + 1 };

                            } else {

                                Cj0.TransformationParams[0, 0] = xNodes1[sjIdx];
                                Cj0.TransformationParams[0, 1] = yNodes[s];
                                Cj0.TransformationParams[1, 0] = xNodes2[sjIdx + 1];
                                Cj0.TransformationParams[1, 1] = yNodes[s + 1];
                                Cj0.TransformationParams[2, 0] = xNodes1[sjIdx + 1];
                                Cj0.TransformationParams[2, 1] = yNodes[s];

                                Cj0.NodeIndices = new long[] {
                                        nIdx1 + sjIdx,
                                        nIdx2 + sjIdx + 1,
                                        nIdx1 + sjIdx + 1 };

                            }

                            grid.Cells[j] = Cj0;
                            j++;
                        }

                    }

                } else {
                    // all MPI processes with rank > 0
                    // return empty grid (use Redistribution)
                    // ++++++++++++++++++++++++++++++++++++++++

                    grid.Cells = new Cell[0];
                }

                return grid;
            }
        }


        /// <summary>
        /// Creates a grid in 2D with hanging nodes. 
        /// Idea: gridBoxes defines the individual region and its resolution. All boxes are put together
        /// and if one gridBox overlaps a previous one, it is cut out of the previous box. Thus, the boxes 
        /// should be ordered from large to small in terms of physical dimension  
        /// </summary>
        /// <param name="periodicX"></param>
        /// <param name="periodicY"></param>
        /// <param name="gridBoxes"></param>
        /// <returns></returns>
        static public GridCommons HangingNodes2D(bool periodicX, bool periodicY, params GridBox[] gridBoxes) {
            if (gridBoxes.Length < 2) {
                throw new ArgumentException("At least 2 GridBoxes are needed for a HangingNodes grid, but only " + gridBoxes.Length + " are specified");
            }
            List<Grid2D> gridList = new List<Grid2D>(gridBoxes.Length);
            GridBox box = gridBoxes[0];
            int dimension = box.boundingBox.D;
            double[][] nodes = new double[2][];

            if (CheckConnectivity2D(box, gridBoxes[1]))
                Console.WriteLine("Adjusted GridBox[1] to ensure Connectivity");

            for (int d = 0; d < dimension; d++) {
                nodes[d] = GenericBlas.Linspace(box.boundingBox.Min[d], box.boundingBox.Max[d], box.numOfCells[d] + 1);
            }
            var grd = Grid2D.Cartesian2DGrid(nodes[0], nodes[1], periodicX: periodicX, periodicY: periodicY, CutOuts: new[] { gridBoxes[1].boundingBox });
            gridList.Add(grd);


            for (int i = 1; i < gridBoxes.Length; i++) {
                box = gridBoxes[i];
                if (i < gridBoxes.Length - 1) {
                    if (CheckConnectivity2D(box, gridBoxes[i + 1]))
                        Console.WriteLine("Adjusted GridBox[{0}] to ensure Connectivity", i + 1);
                }

                // Checks, if smaller Box has also periodic boundaries 
                bool localPeriodicX = false;
                bool localPeriodicY = false;
                if (periodicX) {
                    if (gridBoxes[i - 1].boundingBox.Min[0] == box.boundingBox.Min[0] && gridBoxes[i - 1].boundingBox.Max[0] == box.boundingBox.Max[0]) {
                        localPeriodicX = true;
                    }
                    else if ((gridBoxes[i - 1].boundingBox.Min[0] == box.boundingBox.Min[0] && gridBoxes[i - 1].boundingBox.Max[0] != box.boundingBox.Max[0])
                      || (gridBoxes[i - 1].boundingBox.Min[0] != box.boundingBox.Min[0] && gridBoxes[i - 1].boundingBox.Max[0] == box.boundingBox.Max[0])) {
                        throw new ArgumentException("The defined HangingNodes grid with periodic boundary in x direction has two boxes at one side, but not at corresponding other side, i.e the boxes have not the same length in x direction. Such a grid is not possible!");
                    }
                }
                if (periodicY) {
                    if (gridBoxes[i - 1].boundingBox.Min[1] == box.boundingBox.Min[1] && gridBoxes[i - 1].boundingBox.Max[1] == box.boundingBox.Max[1]) {
                        localPeriodicY = true;
                    }
                    else if ((gridBoxes[i - 1].boundingBox.Min[1] == box.boundingBox.Min[1] && gridBoxes[i - 1].boundingBox.Max[1] != box.boundingBox.Max[1])
                      || (gridBoxes[i - 1].boundingBox.Min[1] != box.boundingBox.Min[1] && gridBoxes[i - 1].boundingBox.Max[1] == box.boundingBox.Max[1])) {
                        throw new ArgumentException("The defined HangingNodes grid with periodic boundary in y direction has two boxes at one side, but not at corresponding other side, i.e the boxes have not the same length in y direction. Such a grid is not possible!");
                    }
                }


                // Creating Grid for box with possible cut out region
                for (int d = 0; d < dimension; d++) {
                    nodes[d] = GenericBlas.Linspace(box.boundingBox.Min[d], box.boundingBox.Max[d], box.numOfCells[d] + 1);
                }

                if (i < gridBoxes.Length - 1) {
                    grd = Grid2D.Cartesian2DGrid(nodes[0], nodes[1], periodicX: localPeriodicX, periodicY: localPeriodicY, CutOuts: gridBoxes[i + 1].boundingBox);
                }
                else {
                    grd = Grid2D.Cartesian2DGrid(nodes[0], nodes[1], periodicX: localPeriodicX, periodicY: localPeriodicY);
                }

                gridList.Add(grd);
            }
            var gridMerged = GridCommons.MergeLogically(gridList.ToArray());
            var grid = GridCommons.Seal(gridMerged, 4);
            return grid;
        }

        /// <summary>
        /// Overload, if no periodic boundaries is needed
        /// </summary>
        /// <param name="gridBoxes"></param>
        /// <returns></returns>
        static public GridCommons HangingNodes2D(params GridBox[] gridBoxes) {
            return HangingNodes2D(false, false, gridBoxes);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="coarse"></param>
        /// <param name="fine"></param>
        /// <returns></returns>
        private static bool CheckConnectivity2D(GridBox coarse, GridBox fine) {
            bool madeAdjustments = false;
            int dimension = coarse.boundingBox.D;
            double[][] nodes = new double[2][];
            for (int d = 0; d < dimension; d++) {
                nodes[d] = GenericBlas.Linspace(coarse.boundingBox.Min[d], coarse.boundingBox.Max[d], coarse.numOfCells[d] + 1);
            }

            double[,] cornerPoints = new double[(int)Math.Pow(2, dimension), dimension];
            // BottomLeft
            cornerPoints[0, 0] = fine.boundingBox.Min[0];
            cornerPoints[0, 1] = fine.boundingBox.Min[1];
            // BottomRight
            cornerPoints[1, 0] = fine.boundingBox.Max[0];
            cornerPoints[1, 1] = fine.boundingBox.Min[1];
            // TopLeft
            cornerPoints[2, 0] = fine.boundingBox.Min[0];
            cornerPoints[2, 1] = fine.boundingBox.Max[1];
            // TopRight
            cornerPoints[3, 0] = fine.boundingBox.Max[0];
            cornerPoints[3, 1] = fine.boundingBox.Max[1];

            //check every cornerPoint
            for (int i = 0; i < (int)Math.Pow(2, dimension); i++) {
                double[] pt = new double[] { cornerPoints[i, 0], cornerPoints[i, 1] };
                // Only check and adjust if the pt is in the coarse bounding box
                if (coarse.boundingBox.Contains(pt)) {

                    for (int d = 0; d < dimension; d++) {
                        for (int j = 0; j < nodes[d].Length - 1; j++) {
                            double x = nodes[d][j];
                            // Coarse and fine are connected in that value, everything is fine
                            if (Math.Abs(x - pt[d]) <= 1.0e-10) {
                                break;
                            }
                            else {
                                //Found right interval
                                if (x - pt[d] < 0.0 && nodes[d][j + 1] - pt[d] > 0.0) {
                                    // Adjust fine boundingBox to nearest neighbor
                                    if (Math.Abs(x - pt[d]) < Math.Abs(nodes[d][j + 1] - pt[d])) {
                                        cornerPoints[i, d] = x;
                                    }
                                    else {
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
        /// Generates a so-called O-grid of a circular domain.
        /// </summary>
        /// <param name="Radius">Radius of the circular domain.</param>
        /// <param name="CenterSectionWidth">With of center section which is meshed Cartesian.</param>
        /// <param name="NoOfCenterNodes">
        /// Number of nodes in the center section, in each direction. this also determines the 
        /// Number of nodes in rotational direction for each of the four ring segments.
        /// </param>
        /// <param name="NoOfRadialNodes">Number of nodes in the radial section.</param>
        /// <param name="type">Cell type.</param>
        /// <returns>
        /// A block-structured O-grid.
        /// </returns>
        public static Grid2D Ogrid(double CenterSectionWidth, double Radius, int NoOfCenterNodes, int NoOfRadialNodes, CellType type = CellType.Square_4) {
            using (new FuncTrace()) {
                if (!(Square.Instance).SupportedCellTypes.Contains(type))
                    throw new ArgumentException("illegal cell type.");
                if (CellTypeExtensions.IsLinear(type)) {
                    throw new ArgumentException("illegal cell type.");
                }
                if (Radius <= 0)
                    throw new ArgumentOutOfRangeException();
                if (CenterSectionWidth <= 0)
                    throw new ArgumentOutOfRangeException();
                if (Radius * Radius <= 0.5 * CenterSectionWidth * CenterSectionWidth)
                    throw new ArgumentOutOfRangeException();
                if (NoOfCenterNodes < 2)
                    throw new ArgumentOutOfRangeException();
                if (NoOfRadialNodes < 2)
                    throw new ArgumentOutOfRangeException();


                MPICollectiveWatchDog.Watch();
                Grid2D grid = new Grid2D(Square.Instance);

                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);




                if (myrank == 0) {


                    // create all cells on process 0
                    // (can be redistributed later on)
                    // ++++++++++++++++++++++++++++++++++
                    List<Cell> AllCells = new List<Cell>();
                    int PtCount = 0;

                    // Reference element nodes
                    // =======================
                    var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Square));
                    NodeSet InterpolationNodes = Kref.GetInterpolationNodes(type);
                    int NoOfNodes = InterpolationNodes.NoOfNodes;

                    // Allocate GlobalID's
                    long[,] CenterSection_Gid = new long[NoOfCenterNodes - 1, NoOfCenterNodes - 1];
                    long cnt = 0;
                    for (int iX = 0; iX < (NoOfCenterNodes - 1); iX++) {
                        for (int iY = 0; iY < (NoOfCenterNodes - 1); iY++) {
                            CenterSection_Gid[iX, iY] = cnt;
                            cnt++;
                        }
                    }

                    long[][,] RingSections_Gid = new long[4][,];
                    for (int _iRing = 0; _iRing < 4; _iRing++) {
                        RingSections_Gid[_iRing] = new long[NoOfRadialNodes - 1, NoOfCenterNodes - 1];
                        for (int iR = 0; iR < (NoOfRadialNodes - 1); iR++) {
                            for (int iPhi = 0; iPhi < (NoOfCenterNodes - 1); iPhi++) {
                                RingSections_Gid[_iRing][iR, iPhi] = cnt;
                                cnt++;
                            }
                        }
                    }


                    // meshing of center section
                    // =========================

                    double[] CenterNodes = GenericBlas.Linspace(-CenterSectionWidth * 0.5, CenterSectionWidth * 0.5, NoOfCenterNodes);


                    for (int iX = 0; iX < (NoOfCenterNodes - 1); iX++) {
                        for (int iY = 0; iY < (NoOfCenterNodes - 1); iY++) {

                            // create cell
                            Cell Cj = new Cell();
                            Cj.GlobalID = CenterSection_Gid[iX, iY];
                            Cj.Type = type;
                            AllCells.Add(Cj);

                            // physical coordinates
                            Cj.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 2);

                            double x0 = CenterNodes[iX];
                            double x1 = CenterNodes[iX + 1];
                            double y0 = CenterNodes[iY];
                            double y1 = CenterNodes[iY + 1];

                            for (int k = 0; k < NoOfNodes; k++) {
                                double x = 0.5 * (x1 - x0) * InterpolationNodes[k, 0] + 0.5 * (x1 + x0);
                                double y = 0.5 * (y1 - y0) * InterpolationNodes[k, 1] + 0.5 * (y1 + y0);
                                Cj.TransformationParams[k, 0] = x;
                                Cj.TransformationParams[k, 1] = y;
                            }

                            // node indices (neighborship via cell face tags
                            Cj.NodeIndices = new long[] { PtCount, PtCount + 1, PtCount + 2, PtCount + 3 };
                            PtCount += 4;

                            // neigborship
                            if (iX > 0) {
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Left,
                                    NeighCell_GlobalID = CenterSection_Gid[iX - 1, iY],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }
                            else {
                                // connection to left ring segment
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Left,
                                    NeighCell_GlobalID = RingSections_Gid[0][0, iY],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }

                            if (iX < (NoOfCenterNodes - 2)) {
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Right,
                                    NeighCell_GlobalID = CenterSection_Gid[iX + 1, iY],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }
                            else {
                                // connection to right ring segment
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Right,
                                    NeighCell_GlobalID = RingSections_Gid[1][0, NoOfCenterNodes - 2 - iY],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }

                            if (iY > 0) {
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Bottom,
                                    NeighCell_GlobalID = CenterSection_Gid[iX, iY - 1],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }
                            else {
                                // connection to bottom ring segment
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Bottom,
                                    NeighCell_GlobalID = RingSections_Gid[3][0, NoOfCenterNodes - 2 - iX],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }

                            if (iY < (NoOfCenterNodes - 2)) {
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Top,
                                    NeighCell_GlobalID = CenterSection_Gid[iX, iY + 1]
                                }, ref Cj.CellFaceTags);
                            }
                            else {
                                // connection to top ring segment
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Top,
                                    NeighCell_GlobalID = RingSections_Gid[2][0, iX],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
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
                                //double beta0 = betaS[NoOfRadialNodes - iR - 1];
                                //double beta1 = betaS[NoOfRadialNodes - iR - 2];
                                double beta0 = betaS[iR];
                                double beta1 = betaS[iR + 1];

                                // create cell
                                Cell Cj = new Cell();
                                Cj.GlobalID = RingSections_Gid[iRing][iR, iPhi];
                                Cj.Type = type;
                                AllCells.Add(Cj);

                                // set nodes
                                Cj.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 2);
                                for (int k = 0; k < NoOfNodes; k++) {
                                    double phi = 0.5 * (phi1 - phi0) * InterpolationNodes[k, 0] + 0.5 * (phi1 + phi0);
                                    double xA = Math.Cos(phi) * Radius;
                                    double yA = Math.Sin(phi) * Radius;

                                    double xC = 0.5 * (xC1 - xC0) * InterpolationNodes[k, 0] + 0.5 * (xC1 + xC0);
                                    double yC = 0.5 * (yC1 - yC0) * InterpolationNodes[k, 0] + 0.5 * (yC1 + yC0);

                                    double beta = 0.5 * (beta1 - beta0) * InterpolationNodes[k, 1] + 0.5 * (beta1 + beta0);

                                    double x = xC * (1 - beta) + xA * beta;
                                    double y = yC * (1 - beta) + yA * beta;

                                    Cj.TransformationParams[k, 0] = x;
                                    Cj.TransformationParams[k, 1] = y;
                                }

                                // node indices (neighborship via cell face tags
                                Cj.NodeIndices = new long[] { PtCount, PtCount + 1, PtCount + 2, PtCount + 3 };
                                PtCount += 4;

                                // neigborship
                                if (iR > 0) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Bottom,
                                        NeighCell_GlobalID = RingSections_Gid[iRing][iR - 1, iPhi],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }
                                else {
                                    // connection to center section

                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Bottom,
                                        NeighCell_GlobalID = CenterSection_Gid[OiX + iPhi * inc_iX, OiY + iPhi * inc_iY],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }

                                if (iR < (NoOfRadialNodes - 2)) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Top,
                                        NeighCell_GlobalID = RingSections_Gid[iRing][iR + 1, iPhi],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }
                                else {
                                    //ArrayTools.AddToArray(new CellFaceTag() {
                                    //        FaceIndex = (int) Square.Edge.Right,
                                    //        NeighCell_GlobalID = RingSections_Gid[1][0, iY]
                                    //    }, ref Cj.CellFaceTags);
                                }

                                if (iPhi > 0) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Left,
                                        NeighCell_GlobalID = RingSections_Gid[iRing][iR, iPhi - 1],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }
                                else {
                                    // connection to previous ring segment
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Left,
                                        NeighCell_GlobalID = RingSections_Gid[iPrevRing][iR, NoOfCenterNodes - 2],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }

                                if (iPhi < (NoOfCenterNodes - 2)) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = RingSections_Gid[iRing][iR, iPhi + 1],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }
                                else {
                                    // connection to next ring segment
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = RingSections_Gid[iNextRing][iR, 0],
                                        ConformalNeighborship = true
                                    }, ref Cj.CellFaceTags);
                                }
                            }
                        }
                    }
                    // finalize
                    // ========

                    grid.Cells = AllCells.ToArray();
                }
                else {
                    grid.Cells = new Cell[0];
                }
                return grid;

            }
        }


        /// <summary>
        /// Special-Purpose grid generator for helical solver project
        /// </summary>
        public static Grid2D HelicalHangingNodes(double rmin, double rmax, double ximin, double ximax, int NoOfRnodes, int NoOfXinodes0,
            int NoOfXiRefinements, int XiRefinementGrades) {

            Grid2D grid = new Grid2D(Square.Instance);


            int myrank;
            int size;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);


            if (myrank == 0) {
                // create all cells on process 0
                // (can be redistributed later on)
                // ++++++++++++++++++++++++++++++++++
                List<Cell> AllCells = new List<Cell>();
                int PtCount = 0;

                // Reference element nodes
                // =======================
                var Kref = grid.RefElements.Single(KK => KK.GetType() == typeof(Square));
                NodeSet InterpolationNodes = Kref.GetInterpolationNodes(CellType.Square_Linear);
                int NoOfNodes = InterpolationNodes.NoOfNodes;

                
                //Case A
                //======
                //double[] rNodes = GenericBlas.Linspace(rmin, rmax, NoOfRnodes);

                //Case B
                //======
                //double[] rNodes = new double[NoOfRnodes];
                //double[] rNodesFunc = GenericBlas.Linspace(0, Math.PI / 2.0, NoOfRnodes);
                //for (int i=0; i<NoOfRnodes; i++) {
                //    rNodes[i] = Math.Sin(rNodesFunc[i]);
                //}

                //Case C
                //======
                // Compute cell length in r direction hr
                double[] hrArray = new double[NoOfRnodes-1];

                List<double> hrList = new List<double>();
                
                for (int i=0; i < NoOfXiRefinements; i++) {
                    double hr = ximax / (NoOfXinodes0 * Math.Pow(XiRefinementGrades, i));
                    hrList.Add(hr);
                }


                //for (int k=0; k < NoOfRnodes-1; k++) {
                //    if (k < NoOfXiRefinements) {
                //        hrArray[k] = ximax / (NoOfXinodes0 * Math.Pow(XiRefinementGrades, k));
                //    }
                //    else {
                //        //hrArray[k] = ximax / (NoOfXinodes0 * Math.Pow(XiRefinementGrades, NoOfXiRefinements-1));
                //    }
                //}
                double hrFinest = ximax / (NoOfXinodes0 * Math.Pow(XiRefinementGrades, NoOfXiRefinements - 1));
                double sum = hrList.Sum();
                double numberOfMissingCells = Math.Round((rmax - sum) / hrFinest)+1;

                //for (int j = 0; j < (int)numberOfMissingCells; j++) {
                //    hrList.Add(hrFinest);
                //}
                while (sum <= rmax){
                    hrList.Add(hrFinest);
                    sum += hrFinest;
                }

                //Debug.Assert(hrList.Sum() < 1.1);


                List<double> rNodesList = new List<double>();
                rNodesList.Add(rmin);

                double[] rNodesTemp = new double[hrList.Count()+1];
                rNodesTemp[0] = rmin;
                for (int i = 1; i < hrList.Count() + 1; i++) {
                    rNodesTemp[i] = rNodesTemp[i-1] + hrList[i-1];
                    if (rNodesTemp[i] <= rmax) {

                        rNodesList.Add(rNodesTemp[i]);
                    }
                    else {
                        //
                    }
                }
                


                rNodesList.Add(rmax);
                NoOfRnodes = rNodesList.Count();

                // if there is a thin cell on the right boundary edge
                if (rmax - rNodesList.Last() < 0.5 * hrList.Last()) {
                    rNodesList.RemoveAt(rNodesList.Count() - 2);
                    NoOfRnodes--;
                }

                double[] rNodes = rNodesList.ToArray();



                // Allocate GlobalID's (umschreiben)
                long[][] CellGid = new long[NoOfRnodes - 1][];
                int[] NoOfCellsPerRcol = new int[CellGid.Length];
                int[] cOffset = new int[CellGid.Length];
                long cnt = 0;
                for(int iR = 0; iR < CellGid.Length; iR++) {
                    int NoOfCells; // no of cells in column
                    if (iR == 0) {
                        NoOfCells = NoOfXinodes0 - 1;
                    }
                    else 
                    {
                        if(iR < NoOfXiRefinements)
                            NoOfCells = NoOfCellsPerRcol[iR - 1] * XiRefinementGrades;
                        else
                            NoOfCells = NoOfCellsPerRcol[iR - 1];
                    }

                    NoOfCellsPerRcol[iR] = NoOfCells;
                    if (iR == 0)
                        cOffset[iR] = 0;
                    else
                        cOffset[iR] = cOffset[iR - 1] + NoOfCellsPerRcol[iR - 1];


                    CellGid[iR] = new long[NoOfCells];
                    for(int iXi =  0; iXi < NoOfCells; iXi++) {
                        CellGid[iR][iXi] = cnt;
                        cnt++;
                    }
                }




                for (int iR = 0; iR < CellGid.Length; iR++) {
                    int NoOfXiNodes = CellGid[iR].Length + 1;
                    double[] xiNodes = GenericBlas.Linspace(ximin, ximax, NoOfXiNodes);
                    
                    for (int iXi = 0; iXi < (NoOfXiNodes - 1); iXi++) {

                        // create cell
                        Cell Cj = new Cell();
                        Cj.GlobalID = CellGid[iR][iXi];
                        Cj.Type = CellType.Square_Linear;
                        AllCells.Add(Cj);

                        // physical coordinates
                        Cj.TransformationParams = MultidimensionalArray.Create(NoOfNodes, 2);

                        double x0 = rNodes[iR];
                        double x1 = rNodes[iR + 1];
                        double y0 = xiNodes[iXi];
                        double y1 = xiNodes[iXi + 1];

                        for (int k = 0; k < NoOfNodes; k++) {
                            double x = 0.5 * (x1 - x0) * InterpolationNodes[k, 0] + 0.5 * (x1 + x0);
                            double y = 0.5 * (y1 - y0) * InterpolationNodes[k, 1] + 0.5 * (y1 + y0);
                            Cj.TransformationParams[k, 0] = x;
                            Cj.TransformationParams[k, 1] = y;
                        }

                        // node indices (neighborship via cell face tags
                        Cj.NodeIndices = new long[] { PtCount, PtCount + 1, PtCount + 2, PtCount + 3 };
                        PtCount += 4;

                        // cell neighborship
                        if (iR > 0) {
                            if (iR >= NoOfXiRefinements) {
                                // nur ein linker Nachbar!

                                Debug.Assert(CellGid[iR].Length == CellGid[iR - 1].Length);

                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Left,
                                    NeighCell_GlobalID = CellGid[iR - 1][iXi],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            }
                            else {
                                // Anzahl an Nachbarn = 'XiRefinementGrades'

                                Debug.Assert(CellGid[iR].Length == CellGid[iR - 1].Length*XiRefinementGrades);
                                
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Left,
                                    NeighCell_GlobalID = CellGid[iR - 1][iXi/XiRefinementGrades],
                                    ConformalNeighborship = false
                                }, ref Cj.CellFaceTags);
                            }
                        } 
                        else 
                        {
                            // no left neighbour
                        }

                        if (iR < (NoOfRnodes - 2)) {
                            if (iR >= NoOfXiRefinements - 1) {
                                // nur ein rechter Nachbar!

                                Debug.Assert(CellGid[iR].Length == CellGid[iR + 1].Length);

                                ArrayTools.AddToArray(new CellFaceTag() {
                                    FaceIndex = (int)Square.Faces.Right,
                                    NeighCell_GlobalID = CellGid[iR + 1][iXi],
                                    ConformalNeighborship = true
                                }, ref Cj.CellFaceTags);
                            } else {
                                // Anzahl an Nachbarn = 'XiRefinementGrades'

                                Debug.Assert(CellGid[iR].Length * XiRefinementGrades == CellGid[iR + 1].Length);

                                for (int iNeigh = 0; iNeigh < XiRefinementGrades; iNeigh++) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        FaceIndex = (int)Square.Faces.Right,
                                        NeighCell_GlobalID = CellGid[iR + 1][iXi * XiRefinementGrades + iNeigh],
                                        ConformalNeighborship = false
                                    }, ref Cj.CellFaceTags);
                                }
                            }
                        }
                        else {
                            // no right neighbour
                        }

                        if (iXi > 0) {
                            ArrayTools.AddToArray(new CellFaceTag() {
                                FaceIndex = (int)Square.Faces.Bottom,
                                NeighCell_GlobalID = CellGid[iR][iXi - 1],
                                ConformalNeighborship = true
                            }, ref Cj.CellFaceTags);
                        }
                        else {
                            // no bottom neighbour
                        }

                        if (iXi < (NoOfXiNodes - 2)) {
                            ArrayTools.AddToArray(new CellFaceTag() {
                                FaceIndex = (int)Square.Faces.Top,
                                NeighCell_GlobalID = CellGid[iR][iXi + 1]
                            }, ref Cj.CellFaceTags);
                        }
                        else {
                            // no top neighbour
                        }

                    }
                }

                grid.Cells = AllCells.ToArray();
            }
            else {
                grid.Cells = new Cell[0];
            }
            

            return grid;
        }
    }



}
