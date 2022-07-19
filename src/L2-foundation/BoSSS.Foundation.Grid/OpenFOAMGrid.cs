using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using ilPSP.Connectors;
using ilPSP.Utils;
using System;
using System.IO;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {
    

    /// <summary>
    /// Instantiation of BoSSS grids from an OpenFOAM polymesh structure
    /// </summary>
    public class OpenFOAMGrid : GridCommons, IForeignLanguageProxy {

        

        IntPtr m_ForeignPtr;

        public string[] BoundaryFacePatchTypes;

        /// <summary>
        /// %
        /// </summary>
        public void _SetForeignPointer(IntPtr ptr) {
            if (ptr == IntPtr.Zero) {
                m_ForeignPtr = IntPtr.Zero;
            } else {

                if (m_ForeignPtr != IntPtr.Zero) {
                    throw new ApplicationException("already registered");
                }
                m_ForeignPtr = ptr;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IntPtr _GetForeignPointer() {
            return m_ForeignPtr;
        }



        /// <summary>
        /// Create a BoSSS grid from an OpenFOAM mesh (polymesh)
        /// </summary>
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
        /// Neighbor cell labels for internal faces (length is <paramref name="nInternalFaces"/>)
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
        [CodeGenExport]
        unsafe public OpenFOAMGrid(
            int nPoints, int nCells, int nFaces, int nInternalFaces, int nNames, int* nameLenghts, int emptyTag,
            int** faces,
            int* vertices_per_face,
            int* neighbour,
            int* owner,
            double* points,
            // byte** names,
            int** names,
            int* patchIDs
        ) :
            base(new[] { Cube.Instance}, new[] { Square.Instance}) //
        {
            try {
                // copy data (unmanaged to managed)
                int[][] _faces = new int[nFaces][];
                int[] _neighbour = new int[nInternalFaces];
                int[] _owner = new int[nFaces];
                double[,] _points = new double[nPoints, 3];
                string[] _names = new string[nNames];
                int[] _patchIDs = new int[nFaces - nInternalFaces];

                //Debug.Assert(nFaces == GridImportTest.faces.Length, "mism nFaces len");
                for(int i = 0; i < nFaces; i++) {
                    int N = vertices_per_face[i];
                    //Debug.Assert(N == GridImportTest.faces[i].Length, "mism nVerts face " + i);
                    _faces[i] = new int[N];
                    for(int n = 0; n < N; n++) {
                        _faces[i][n] = faces[i][n];
                        //Debug.Assert(_faces[i][n] == GridImportTest.faces[i][n], "mism face " + i + " vert " + n);
                    }
                }

                //Debug.Assert(nInternalFaces == GridImportTest.neighbour.Length, "mismatch neighbour length");
                for(int i = 0; i < nInternalFaces; i++) {
                    _neighbour[i] = neighbour[i];
                    //Debug.Assert(_neighbour[i] == GridImportTest.neighbour[i], "neighbour " + i + " mismatch ");
                }

                for(int i = 0; i < nFaces; i++) {
                    _owner[i] = owner[i];
                    //Debug.Assert(_owner[i] == GridImportTest.owner[i], "owner " + i + " mismatch ");
                }

                for(int i = 0; i < nFaces - nInternalFaces; i++) {
                    _patchIDs[i] = patchIDs[i];
                }

                for(int i = 0; i < nPoints; i++) {
                    _points[i, 0] = points[i * 3 + 0];
                    _points[i, 1] = points[i * 3 + 1];
                    _points[i, 2] = points[i * 3 + 2];
                }

                for(int i = 0; i < nNames; i++) {
                    int nameLenght = nameLenghts[i];
                    char[] _name = new char[nameLenght];
                    for (int j = 0; j < nameLenght; j++){
                        _name[j] = (char)(names[i][j]);
                    }
                    _names[i] = String.Join("", _name);
                }

                // save boundary face information
                // =====================

                int nBoundaryFaces = nFaces - nInternalFaces;
                this.BoundaryFacePatchTypes = new string[nBoundaryFaces];
                for (int i = 0; i < nBoundaryFaces; i++)
                {
                    this.BoundaryFacePatchTypes[i] = _names[_patchIDs[i]];
                    // Console.WriteLine("BoundaryFacepatchtypes");
                    // Console.WriteLine(this.BoundaryFacePatchTypes[i]);
                }


                // create BoSSS grid
                FOAMmesh_to_BoSSS(this, nCells, _faces, _neighbour, _owner, _points, _names, _patchIDs, emptyTag);

                foreach (var name in _names){
                    this.AddEdgeTag(name);
                }
                // create grid data object
                this.GridDataObject = new GridData(this);
            } catch(Exception e) {
                Console.Error.WriteLine(e.GetType() + ": " + e.Message);
                Console.Error.WriteLine(e.StackTrace);
            }
        }

        /// <summary>
        /// Create a BoSSS grid from an OpenFOAM mesh
        /// </summary>
        /// <param name="faces">
        /// point labels (indices) of faces, according to OpenFOAM polymesh description (array of arrays)
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
        /// <param name="points">
        /// point coordinates, array of <paramref name="nPoints"/>*3
        /// </param>
        public OpenFOAMGrid(
            int nCells, 
            int[][] faces,
            int[] neighbour,
            int[] owner,
            double[,] points,
            string[] names,
            int[] patchIDs,
            int emptyTag
        ) :
            base(new[] { Cube.Instance }, new[] { Square.Instance }) //
        {


            // create BoSSS grid
            FOAMmesh_to_BoSSS(this, nCells, faces, neighbour, owner, points, names, patchIDs, emptyTag);

                Console.Write("test4");
            // create grid data object
            this.GridDataObject = new GridData(this);
                Console.Write("test5");
        }


        /// <summary>
        /// Only for testing purposes
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        [CodeGenExport]
        public int TestMethod(int a) {
            Console.WriteLine("This is a BoSSS mesh with " + CellPartitioning.TotalLength + " cells on " + CellPartitioning.MpiSize + " MPI processes.");

            //Console.WriteLine("C#: Got number " + a + " from external code; returning " + (a * 2));
            return (int)CellPartitioning.TotalLength;
        }


        /// <summary>
        /// The grid metrics
        /// </summary>
        public GridData GridDataObject {
            get;
            private set;
        }

        /// <summary>
        /// <see cref="GridDataObject"/>
        /// </summary>
        [CodeGenExport]
        public GridData GetGridData() {
            return GridDataObject;
        }


        internal static void FOAMmesh_to_BoSSS(GridCommons grid, int nCells, int[][] faces, int[] neighbour, int[] owner, double[,] points, string[] names, int[] patchIDs, int emptyTag) {

            // write everything to a file for debugging
            File.WriteAllLines("neighbour", neighbour.Select(i=>i.ToString()).ToArray());
            File.WriteAllLines("owner", owner.Select(i=>i.ToString()).ToArray());
            // System.IO.File.WriteAllLines("faces", faces.Select(i=>i.ToString()).ToArray());
            File.WriteAllLines("faces", faces.Select(line => String.Join(" ", line)));

            using (var sw = new StreamWriter("points"))
            {
                for (int i = 0; i < points.Length/3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(points[i, j] + " ");
                    }
                    sw.Write("\n");
                }

                sw.Flush();
                sw.Close();
            }

            // Checks
            // ======


            int nFaces = faces.Length;
            int nInternalFaces = neighbour.Length;
            Debug.Assert(nFaces == owner.Length);
            int nPoints = points.GetLength(0);

            if (owner.Length != faces.Length)
                throw new ArgumentException("mismatch between length of faces and owner array");
#if DEBUG
            // extended sanity checks
            for (int i = 0; i < faces.Length; i++) {
                int[] pts = faces[i];

                foreach (int iVtx in pts) {
                    if (iVtx < 0)
                        throw new ArgumentException("negative vertex index for face " + i);
                    if (iVtx >= points.GetLength(0))
                        throw new ArgumentException("vertex index out-of-range for face " + i);

                }

                if (owner[i] < 0 || owner[i] >= nCells)
                    throw new ArgumentException("face owner out-of-range for face " + i);
            }

            for (int i = 0; i < neighbour.Length; i++) {
                if (owner[i] < 0 || owner[i] >= nCells)
                    throw new ArgumentException("face neighbor out-of-range for face " + i);
            }
#endif

            // Build Cells-to-Faces correlation
            // ================================
            Console.WriteLine("faces");
            int k = 0;
            foreach (var face in faces){
                Console.WriteLine("face " + k);
                k++;
                foreach (var i in face)
                    Console.WriteLine(i);
            }

            List<int>[] Cells2Faces = new List<int>[nCells];
            for (int j = 0; j < nCells; j++) {
                Cells2Faces[j] = new List<int>();
            }

            for (int iInternalFace = 0; iInternalFace < nInternalFaces; iInternalFace++) {
                int NeighCell = neighbour[iInternalFace];

                Debug.Assert(Cells2Faces[NeighCell].Contains(iInternalFace) == false);
                Cells2Faces[NeighCell].Add(iInternalFace);
            }

            for (int iFace = 0; iFace < nFaces; iFace++) {
                Debug.Assert(Cells2Faces[owner[iFace]].Contains(iFace) == false);
                Cells2Faces[owner[iFace]].Add(iFace);

            }

            int index = 0;
            foreach (var cfs in Cells2Faces){
                Console.WriteLine("Cell " + index);
                foreach(var face in cfs){
                    Console.WriteLine(face);
                }
                index++;
            }

            // convert to BoSSS grid
            // =====================

            Cell[] bosss_cells = new Cell[nCells];

            for (int iCell = 0; iCell < nCells; iCell++) {
                var Cell2Faces = Cells2Faces[iCell];

                if (Cell2Faces.Count == 6
                   && !Cell2Faces.Select(iFace => faces[iFace].Length == 4).Contains(false)) {
                    // +++++++++++++++
                    // found some Cube
                    // +++++++++++++++

                    int[][] cube_faces = Cell2Faces.Select(iFace => faces[iFace]).ToArray();
                    Debug.Assert(cube_faces.Length == 6);

                    HashSet<int> NodesUnsorted = new HashSet<int>();
                    foreach (int iFace in Cell2Faces) {
                        foreach (int iPoint in faces[iFace]) {
                            NodesUnsorted.Add(iPoint);
                        }
                    }
                    if (NodesUnsorted.Count != 8) {
                        throw new NotImplementedException("Degenerate cubes are not supported.");
                    }
                    int iii = 0;
                    Console.WriteLine("cubefaces of cell " + iCell);
                    foreach (var cf in cube_faces){
                        Console.WriteLine("face " + iii);
                        foreach (var f in cf){
                            Console.WriteLine(f);
                        }
                        iii++;
                    }
                    // Console.WriteLine("Hello from loop for Cell " + iCell);
                    // // int index = 0;
                    // foreach (var node in NodesUnsorted){
                    //     // Console.WriteLine("Node " + index);
                    //     Console.WriteLine("Node ");
                    //     Console.WriteLine(node);
                    //     // foreach (var coord in node){
                    //     //     Console.WriteLine(coord);
                    //     // }
                    //     // index++;
                    // }

                    // create cell
                    Cell cl = new Cell();
                    bosss_cells[iCell] = cl;
                    cl.GlobalID = iCell;
                    cl.Type = CellType.Cube_8;
                    cl.NodeIndices = new long[8];
                    cl.TransformationParams = MultidimensionalArray.Create(8, 3);

                    // Find point indices for bottom, i.e. 'y == -1' (in ref. coordinates):
                    // TODO face 0 is not always the bottom one -> might be fixed by transformationparams
                    cl.NodeIndices[0] = cube_faces[0][0];
                    cl.NodeIndices[1] = cube_faces[0][1];
                    cl.NodeIndices[6] = cube_faces[0][2];
                    cl.NodeIndices[3] = cube_faces[0][3];

                    // Find point indices for 'x == -1'-wall (in ref. coordinates):
                    bool found = false;
                    for (int iCF = 1; iCF < 6; iCF++) {
                        found = ContainsEdge(cube_faces[iCF], cl.NodeIndices[0], cl.NodeIndices[3], out cl.NodeIndices[2], out cl.NodeIndices[7]);
                        if (found)
                            break;
                    }
                    if (found == false)
                        throw new NotSupportedException("weird cell topology in cell " + iCell + "; assuming a cube/hex;");

                    // Find point indices for 'x == +1'-wall (in ref. coordinates):
                    found = false;
                    for (int iCF = 1; iCF < 6; iCF++) {
                        found = ContainsEdge(cube_faces[iCF], cl.NodeIndices[1], cl.NodeIndices[6], out cl.NodeIndices[4], out cl.NodeIndices[5]);
                        if (found)
                            break;
                    }
                    if (found == false)
                        throw new NotSupportedException("weird cell topology in cell " + iCell + "; assuming a cube/hex;");


                    // set point coordinates
                    for (int i = 0; i < 8; i++) {
                        cl.TransformationParams.SetRow(i, points.GetRow(checked((int)(cl.NodeIndices[i]))));
                    }

                    cl.CellFaceTags = new CellFaceTag[6];

                    // in order to establish the neighbor relationships of empty faces, we first have to collect them
                    var emptyFacesOfCell = new List<int>();
                    for (int i = 0; i < 6; i++){
                        int faceIndex = Cells2Faces[iCell][i];
                        if (faceIndex >= nInternalFaces && patchIDs[faceIndex - nInternalFaces] == emptyTag){
                            throw new NotSupportedException("empty patches are currently not supported");
                            emptyFacesOfCell.Add(faceIndex);
                            // Console.WriteLine("Found empty face " + faceIndex);
                        }
                    }
                    int iEmpty = 0;

                    // loop over all faces in order to set the CellFaceTags appropriately.
                    // This is important for the specification of boundary conditions.
                    for (int i = 0; i < 6; i++) {
                        int faceIndex = Cells2Faces[iCell][i];
                        cl.CellFaceTags[i] = new CellFaceTag();
                        if (faceIndex < nInternalFaces) {
                            // this is an internal face
                            Console.WriteLine("Face Nr " + faceIndex + " is an internal face, local index " + i);
                            // cl.CellFaceTags[i].FaceIndex = faceIndex;
                            cl.CellFaceTags[i].FaceIndex = i;
                            cl.CellFaceTags[i].EdgeTag = 0;
                            // cl.CellFaceTags[i].NeighCell_GlobalID = neighbour[faceIndex];
                            cl.CellFaceTags[i].NeighCell_GlobalID = -1;
                            cl.CellFaceTags[i].ConformalNeighborship = true;
                        } else if (patchIDs[faceIndex - nInternalFaces] == emptyTag) {
                            throw new NotSupportedException("empty patches are currently not supported");
                            // this is an empty face
                            cl.CellFaceTags[i].FaceIndex = faceIndex;
                            // in degenerate dimensions, there can only be one cell layer. Therefore, each cell is its own neighbor
                            cl.CellFaceTags[i].NeighCell_GlobalID = iCell;
                            cl.CellFaceTags[i].ConformalNeighborship = true;
                            int D = 3; // since OpenFOAM only knows 3D grids, it seems acceptable to hardcode this

                            int faceIndexNeighbor;
                            if (iEmpty % 2 == 0){
                                faceIndexNeighbor = emptyFacesOfCell[iEmpty + 1];
                                cl.CellFaceTags[i].PeriodicInverse = true;
                            } else {
                                faceIndexNeighbor = emptyFacesOfCell[iEmpty - 1];
                                cl.CellFaceTags[i].PeriodicInverse = false;
                            }
                            iEmpty++;

                            var pointsOfFace = new List<List<double>>();
                            for (int iPoint = 0; iPoint < 4; iPoint++){
                                pointsOfFace.Add(new List<double>());
                                // Console.WriteLine("point: " + iPoint);
                                for (int iCoord = 0; iCoord < D ; iCoord++){
                                    pointsOfFace[iPoint].Add(points[faces[faceIndex][iPoint], iCoord]);
                                    // Console.WriteLine(pointsOfFace[iPoint][iCoord]);
                                }
                            }
                            var pointsOfNeighborFace = new List<List<double>>();
                            for (int iPoint = 0; iPoint < 4; iPoint++){
                                pointsOfNeighborFace.Add(new List<double>());
                                // Console.WriteLine("point of neighbor: " + iPoint);
                                for (int iCoord = 0; iCoord < D ; iCoord++){
                                    pointsOfNeighborFace[iPoint].Add(points[faces[faceIndexNeighbor][iPoint], iCoord]);
                                    // Console.WriteLine(pointsOfNeighborFace[iPoint][iCoord]);
                                }
                            }
                            byte periodicTag = 0;
                            Vector[] inlet = new Vector[D];
                            for (int j = 0; j < D; j++){ // TODO do we have to do this for the entire patch, not individual cells?
                                inlet[j] = new Vector();
                                for (int jj = 0; jj < D; jj++){
                                    inlet[j].Add(pointsOfFace[j][jj]);
                                }
                            }
                            Vector[] outlet = new Vector[D];
                            for (int j = 0; j < D; j++){ // TODO do we have to do this for the entire patch, not individual cells?
                                outlet[j] = new Vector();
                                for (int jj = 0; jj < D; jj++){
                                    outlet[j].Add(pointsOfNeighborFace[j][jj]);
                                }
                            }

                            var inletVec1 = inlet[0] - inlet[1];
                            var inletVec2 = inlet[0] - inlet[2];
                            var inletNormal = inletVec1.CrossProduct(inletVec2);

                            var outletVec1 = outlet[0] - outlet[1];
                            var outletVec2 = outlet[0] - outlet[2];
                            var outletNormal = outletVec1.CrossProduct(outletVec2);

                            // Console.WriteLine("inletNormalVec" + inletNormal);
                            // Console.WriteLine("outletNormalVec" + outletNormal);
                            grid.ConstructPeriodicEdgeTrafo(outlet, outletNormal, inlet, inletNormal, out periodicTag);
                            cl.CellFaceTags[i].EdgeTag = periodicTag;
                            // cl.CellFaceTags[i].EdgeTag = 182;
                            Console.WriteLine("Face Nr " + faceIndex + " is an empty face, GlobalID: " + iCell + " , EdgeTag: " + cl.CellFaceTags[i].EdgeTag + ", neighborFaceIndex: " + faceIndexNeighbor + ", localFaceIndex: " + i + " , periodicInverse: " + cl.CellFaceTags[i].PeriodicInverse);
                        } else {
                            // this is an actual boundary face
                            cl.CellFaceTags[i].EdgeTag = (byte)(patchIDs[faceIndex - nInternalFaces] + 1);
                            cl.CellFaceTags[i].NeighCell_GlobalID = -1;
                            // cl.CellFaceTags[i].FaceIndex = faceIndex;
                            cl.CellFaceTags[i].FaceIndex = i;
                            cl.CellFaceTags[i].ConformalNeighborship = true;
                            Console.WriteLine("Face Nr " + faceIndex + " is a boundary face, EdgeTag: " + cl.CellFaceTags[i].EdgeTag + ", local index: " + i);
                            if (faceIndex == 9){
                                Console.WriteLine("nodes: ");
                                foreach (var node in faces[faceIndex]){
                                    Console.WriteLine(node);
                                }
                                foreach (var node in faces[faceIndex]){
                                    for (int kk = 0; kk < 3; kk++)
                                    {Console.WriteLine(points[node,kk]);}
                                }
                            }
                        }
                        // Console.WriteLine("CellFaceTag: " + cl.CellFaceTag);
                    }

                    // Jacobian Test
                    JacobianTest(cl, out bool PositiveJacobianFlag, out bool NegativeJacobianFlag, out bool Linear);
                    if (NegativeJacobianFlag) {
                        // flip top/bottom
                        long _0 = cl.NodeIndices[0];
                        long _1 = cl.NodeIndices[1];
                        long _6 = cl.NodeIndices[6];
                        long _3 = cl.NodeIndices[3];

                        cl.NodeIndices[0] = cl.NodeIndices[2];
                        cl.NodeIndices[1] = cl.NodeIndices[4];
                        cl.NodeIndices[6] = cl.NodeIndices[5];
                        cl.NodeIndices[3] = cl.NodeIndices[7];

                        cl.NodeIndices[2] = _0;
                        cl.NodeIndices[4] = _1;
                        cl.NodeIndices[5] = _6;
                        cl.NodeIndices[7] = _3;

                        for (int i = 0; i < 8; i++) {
                            cl.TransformationParams.SetRow(i, points.GetRow(checked((int)(cl.NodeIndices[i]))));
                        }

                        JacobianTest(cl, out PositiveJacobianFlag, out NegativeJacobianFlag, out Linear);
                    }

                    if (NegativeJacobianFlag || !PositiveJacobianFlag) {
                        throw new NotSupportedException("Found degenerate Jacobian in cell " + iCell);
                    }
                    Debug.Assert(PositiveJacobianFlag == true);
                    Debug.Assert(NegativeJacobianFlag == false);

                    Console.WriteLine("transformationparams of cell " + iCell);
                    for (int ii = 0; ii < 8; ii++){
                        Console.WriteLine("row " + ii);
                        Console.WriteLine("nodeindex " + cl.NodeIndices[ii]);
                        for (int jj = 0; jj < 3; jj++){
                            Console.Write(cl.TransformationParams[ii,jj]);
                            Console.Write(" ");
                        }
                        Console.WriteLine();
                    }

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

            grid.Cells = bosss_cells;
            // grid.DefineEdgeTags(delegate (double[] X){ // TODO generalize
            //         if (Math.Abs(X[0] - 0) < 1e-10){
            //             return 1;
            //         }
            //         if (Math.Abs(X[0] - 5) < 1e-10){
            //             return 2;
            //         }
            //         if (Math.Abs(X[1] - 0) < 1e-10){
            //             return 3;
            //         }
            //         if (Math.Abs(X[1] - 1) < 1e-10){
            //             return 3;
            //         }
            //         if (Math.Abs(X[2] - 0) < 1e-10){
            //             return 3;
            //         }
            //         if (Math.Abs(X[2] - 1) < 1e-10){
            //             return 3;
            //         }
            //         // return 3;
            //         Console.WriteLine("Argument out of range!");
            //         throw new ArgumentOutOfRangeException();
            //     }
            // );
            // grid.BcCells = new BCElement[nFaces - nInternalFaces];
            // for (int i = 0; i < nFaces - nInternalFaces; i++){
            //     grid.BcCells[i] = new BCElement();
            //     grid.BcCells[i].EdgeTag = (byte)(patchIDs[i - nInternalFaces] + 1);
            // }
            grid.Description = "imported form OpenFOAM";
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

        static bool ContainsEdge(int[] face, long iPt1, long iPt2, out long iPrev, out long iNext) {
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
                    iPrev = _iNext;
                    return true;
                }
            }
            iNext = int.MinValue;
            iPrev = int.MinValue;
            return false;
        }

    }

}
