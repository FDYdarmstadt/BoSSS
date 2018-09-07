using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding {


    /// <summary>
    /// Instantiation of BoSSS grids
    /// </summary>
    public static class Grid_ {

        /// <summary>
        /// Create a BoSSS grid from an OpenFOAM mesh
        /// </summary>
        /// <param name="_ref"></param>
        /// <param name="ierr">output, error code; on success, set to 0</param>
        /// <param name="faces">
        /// point labels (indices) of faces, according to OpenFOAM polymesh description (array of arrays)
        /// </param>
        /// <param name="vertices_per_face">
        /// number of vertices for each face
        /// </param>
        /// <param name="nCells">
        /// number of cells in OpenFOAM polymesh
        /// </param>
        /// <param name="neighbour">
        /// Neighbor cell labels for internal faces
        /// </param>
        /// <param name="owner">
        /// for each face, the owner cell label
        /// </param>
        /// <param name="nPoints">number of points/vertices</param>
        /// <param name="nFaces">
        /// number of faces
        /// </param>
        /// <param name="nInternalFaces">
        /// number of internal faces
        /// </param>
        /// <param name="points">
        /// point coordinates, array of <paramref name="nPoints"/>*3
        /// </param>
        unsafe public static void CreateGrid(out int _ref, 
            ref int nPoints, ref int nCells, ref int nFaces, ref int nInternalFaces,
            int** faces,
            int* vertices_per_face,
            int* neighbour,
            int* owner,
            double* points,
            out int ierr) {

            try {
                // copy data (unmanaged to managed)
                int[][] _faces = new int[nFaces][];
                int[] _neighbour = new int[nInternalFaces];
                int[] _owner = new int[nFaces];
                double[,] _points = new double[nPoints, 3];

                for(int i = 0; i < nFaces; i++) {
                    int N = vertices_per_face[i];
                    _faces[i] = new int[N];
                    for (int n = 0; n < N; n++)
                        _faces[i][n] = faces[i][n];
                }

                for (int i = 0; i < nInternalFaces; i++) {
                    _neighbour[i] = neighbour[i];
                }

                for (int i = 0; i < nFaces; i++) {
                    _owner[i] = owner[i];
                }

                for (int i = 0; i < nPoints; i++) {
                    _points[i, 0] = points[i * 3 + 0];
                    _points[i, 0] = points[i * 3 + 1];
                    _points[i, 0] = points[i * 3 + 2];
                }

                
                // create BoSSS grid
                GridData g = FOAMmesh_to_BoSSS(nCells, _faces, _neighbour, _owner, _points);

                // register object
                _ref = Infrastructure.RegisterObject(g);
                ierr = 0;
                return;

            } catch(Exception e) {
                ierr = Infrastructure.ErrorHandler(e);
                _ref = -1;
                return;
            }
        }


        internal static GridData FOAMmesh_to_BoSSS(int nCells, int[][] faces, int[] neighbour, int[] owner, double[,] points) { 

            // Checks
            // ======

            Common_.BoSSSInitialize();
            int nFaces = faces.Length;
            int nInternalFaces = neighbour.Length;
            Debug.Assert(nFaces == owner.Length);
            int nPoints = points.GetLength(0);

            if (owner.Length != faces.Length)
                throw new ArgumentException("mismatch between length of faces and owner array");
#if DEBUG
            // extended sanity checks
            for(int i = 0; i < faces.Length; i++) {
                int[] pts = faces[i];

                foreach(int iVtx in pts) {
                    if (i < 0)
                        throw new ArgumentException("negative vertex index for face " + i);
                    if (i >= points.GetLength(0))
                        throw new ArgumentException("vertex index out-of-range for face " + i);

                }

                if (owner[i] < 0 || owner[i] >= nCells)
                    throw new ArgumentException("face owner out-of-range for face " + i);
            }

            for(int i = 0; i < neighbour.Length; i++) {
                if (owner[i] < 0 || owner[i] >= nCells)
                    throw new ArgumentException("face neighbor out-of-range for face " + i);
            }
#endif

            // Build Cells-to-Faces correlation
            // ================================

            List<int>[] Cells2Faces = new List<int>[nCells];
            for(int j = 0; j < nCells; j++) {
                Cells2Faces[j] = new List<int>();
            }

            for(int iFace = 0; iFace < nFaces; iFace++) {
                Debug.Assert(Cells2Faces[owner[iFace]].Contains(iFace) == false);
                Cells2Faces[owner[iFace]].Add(iFace);

            }

            for(int iInternalFace = 0; iInternalFace < nInternalFaces; iInternalFace++) {
                int NeighCell = neighbour[iInternalFace];

                Debug.Assert(Cells2Faces[NeighCell].Contains(iInternalFace) == false);
                Cells2Faces[NeighCell].Add(iInternalFace);
            }

            // convert to BoSSS grid
            // =====================

            Cell[] bosss_cells = new Cell[nCells];

            for(int iCell = 0; iCell < nCells; iCell++) {
                var Cell2Faces = Cells2Faces[iCell];

                if(    Cell2Faces.Count == 6
                   && !Cell2Faces.Select(iFace => faces[iFace].Length == 4).Contains(false)) {
                    // +++++++++++++++
                    // found some Cube
                    // +++++++++++++++

                    int[][] cube_faces = Cell2Faces.Select(iFace => faces[iFace]).ToArray();
                    Debug.Assert(cube_faces.Length == 6);

                    HashSet<int> NodesUnsorted = new HashSet<int>();
                    foreach(int iFace in Cell2Faces) {
                        foreach(int iPoint in faces[iFace]) {
                            NodesUnsorted.Add(iPoint);
                        }
                    }
                    if(NodesUnsorted.Count != 8) {
                        throw new NotImplementedException("Degenerate cubes are not supported.");
                    }

                    // create cell
                    Cell cl = new Cell();
                    bosss_cells[iCell] = cl;
                    cl.GlobalID = iCell;
                    cl.Type = CellType.Cube_8;
                    cl.NodeIndices = new int[8];
                    cl.TransformationParams = MultidimensionalArray.Create(8, 3);

                    // Find point indices for bottom, i.e. 'y == -1' (in ref. coordinates):
                    cl.NodeIndices[0] = cube_faces[0][0];
                    cl.NodeIndices[1] = cube_faces[0][1];
                    cl.NodeIndices[6] = cube_faces[0][2];
                    cl.NodeIndices[3] = cube_faces[0][3];


                    // Find point indices for 'x == -1'-wall (in ref. coordinates):
                    bool found = false;
                    for(int iCF = 1; iCF < 6; iCF++) {
                        found = ContainsEdge(cube_faces[iCF], cl.NodeIndices[0], cl.NodeIndices[3], out cl.NodeIndices[2], out cl.NodeIndices[7]);
                        if (found)
                            break;
                    }
                    if (found == false)
                        throw new NotSupportedException("weired cell topology in cell " + iCell + "; assuming a cube/hex;");

                    // Find point indices for 'x == +1'-wall (in ref. coordinates):
                    found = false;
                    for(int iCF = 1; iCF < 6; iCF++) {
                        found = ContainsEdge(cube_faces[iCF], cl.NodeIndices[1], cl.NodeIndices[6], out cl.NodeIndices[4], out cl.NodeIndices[5]);
                        if (found)
                            break;
                    }
                    if (found == false)
                        throw new NotSupportedException("weired cell topology in cell " + iCell + "; assuming a cube/hex;");

                    
                    // set point coordinates
                    for(int i = 0; i < 8; i++) {
                        cl.TransformationParams.SetRow(i, points.GetRow(cl.NodeIndices[i]));
                    }

                    // Jacobian Test
                    JacobianTest(cl, out bool PositiveJacobianFlag, out bool NegativeJacobianFlag, out bool Linear);
                    if(NegativeJacobianFlag) {
                        // flip top/bottom
                        int _0 = cl.NodeIndices[0];
                        int _1 = cl.NodeIndices[1];
                        int _6 = cl.NodeIndices[6];
                        int _3 = cl.NodeIndices[3];

                        cl.NodeIndices[0] = cl.NodeIndices[2];
                        cl.NodeIndices[1] = cl.NodeIndices[4];
                        cl.NodeIndices[6] = cl.NodeIndices[5];
                        cl.NodeIndices[3] = cl.NodeIndices[7];

                        cl.NodeIndices[2] = _0;
                        cl.NodeIndices[4] = _1;
                        cl.NodeIndices[5] = _6;
                        cl.NodeIndices[7] = _3;
                        
                        for (int i = 0; i < 8; i++) {
                            cl.TransformationParams.SetRow(i, points.GetRow(cl.NodeIndices[i]));
                        }

                        JacobianTest(cl, out PositiveJacobianFlag, out NegativeJacobianFlag, out Linear);
                    }

                    if(NegativeJacobianFlag || !PositiveJacobianFlag) {
                        throw new NotSupportedException("Found degenerate Jacobian in cell " + iCell);
                    }
                    Debug.Assert(PositiveJacobianFlag == true);
                    Debug.Assert(NegativeJacobianFlag == false);

                    if (Linear) {
                        cl.Type = CellType.Cube_Linear;
                        cl.TransformationParams = cl.TransformationParams.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { 3, 2 }).CloneAs();
                    }
                    

                } else {
                    throw new NotImplementedException("At the moment, only support for hex cells.");
                }
            }

            // create BoSSS grid
            // ==================

            GridCommons grd = new GridCommons(new RefElement[] { Cube.Instance }, new RefElement[] { Square.Instance }) {
                Cells = bosss_cells,
                Description = "imported form OpenFOAM"
            };

            var gDat = new GridData(grd);


            return gDat;
            
        }

        static void JacobianTest(Cell cl, out bool PositiveJacobianFlag, out bool NegativeJacobianFlag, out bool Linear) {
            RefElement Kref = cl.Type.GetRefElement();

            int D = 3; // for OpenFOAM, we assume 3D;

            // get nodes within ref element, to test the Jacobian
            int deg = Kref.GetInterpolationDegree(cl.Type);
            if (deg > 1)
                deg--;
            deg *= 3;

            NodeSet TestNodes = Kref.GetQuadratureRule(2 * deg).Nodes;

            // evaluate derivatives of nodal polynomials for transformation
            PolynomialList[] Deriv = Kref.GetInterpolationPolynomials1stDeriv(cl.Type);
            MultidimensionalArray[] DerivEval = new MultidimensionalArray[D];
            for (int d = 0; d < D; d++) {
                DerivEval[d] = Deriv[d].Values.GetValues(TestNodes);
            }

            // evaluate Jacobian matrix
            MultidimensionalArray Jacobi = MultidimensionalArray.Create(TestNodes.NoOfNodes, D, D); // temporary storage for Jacobian matrix

            Debug.Assert(cl.TransformationParams.Dimension == 2);
            Debug.Assert(cl.TransformationParams.GetLength(1) == 3);

            for (int d1 = 0; d1 < D; d1++) {
                Debug.Assert(cl.TransformationParams.GetLength(0) == Deriv[d1].Count);
                MultidimensionalArray JacobiCol = Jacobi.ExtractSubArrayShallow(-1, -1, d1);
                JacobiCol.Multiply(1.0, DerivEval[d1], cl.TransformationParams, 0.0, "kd", "kn", "nd");
            }

            // do the tests
            NegativeJacobianFlag = false;
            PositiveJacobianFlag = false;
            double MinAbsJacobiDet = double.MaxValue;
            double MaxAbsJacobiDet = 0.0;
            for (int n = 0; n < TestNodes.NoOfNodes; n++) {
                double detJac = Jacobi.ExtractSubArrayShallow(n, -1, -1).Determinant();
                if (detJac <= 0) {
                    NegativeJacobianFlag = true;
                }
                if (detJac > 0) {
                    PositiveJacobianFlag = true;
                }

                MinAbsJacobiDet = Math.Min(MinAbsJacobiDet, detJac.Abs());
                MaxAbsJacobiDet = Math.Max(MaxAbsJacobiDet, detJac.Abs());
            }

            if ((MaxAbsJacobiDet - MinAbsJacobiDet) / MinAbsJacobiDet <= 1.0e-8)
                Linear = true;
            else
                Linear = false;

        }

        static bool ContainsEdge(int[] face, int iPt1, int iPt2, out int iPrev, out int iNext) {
            int L = face.Length;
            for (int i = 0; i < L; i++) {
                int v1 = face[i];
                int v2 = face[(i + 1) % L];
                int _iPrev = face[i > 0 ? (i - 1) : (L - 1)];
                int _iNext = face[(i + 2) % L];

                if (v1 == iPt1 && v2 == iPt2) {
                    iPrev = _iPrev;
                    iNext = _iNext;
                    return true;
                }
                if (v2 == iPt1 && v1 == iPt2) {
                    iNext = _iPrev;
                    iPrev= _iNext;
                    return true;
                }
            }
            iNext = int.MinValue;
            iPrev = int.MinValue;
            return false;
        }

    }   
}
