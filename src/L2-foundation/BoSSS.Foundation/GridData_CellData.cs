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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Comm;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using System.Collections;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {

    partial class GridData {

        /// <summary>
        /// <see cref="Cells"/>
        /// </summary>
        private CellData m_Cells;

        /// <summary>
        /// metrics and operations which are associated to one cell
        /// </summary>
        public CellData Cells {
            get {
                return m_Cells;
            }
        }

        public IGeometricalCellsData iGeomCells {
            get {
                return m_Cells;
            }
        }

        public ILogicalCellData iLogicalCells {
            get {
                return m_Cells;
            }
        }



        /// <summary>
        /// all metrics which are associated to one cell
        /// </summary>
        public class CellData : IGeometricalCellsData, ILogicalCellData {
            
            /// <summary>
            /// All reference elements for cells, see <see cref="GetRefElementIndex(int)"/> resp. <see cref="GetRefElement(int)"/>.
            /// </summary>
            public RefElement[] RefElements {
                get {
                    return m_owner.Grid.RefElements;
                }
            }

            CellMask[] m_RefElementMask;

            /// <summary>
            /// Cell-Mask of all cells which share the same reference element.
            /// </summary>
            public CellMask GetCells4Refelement(int iKref) {
                var KRefs = this.RefElements;
                if (m_RefElementMask == null) {
                    m_RefElementMask = new CellMask[KRefs.Length];
                }

                if (iKref < 0 || iKref >= KRefs.Length) {
                    throw new ArgumentOutOfRangeException();
                }

                if (m_RefElementMask[iKref] == null) {
                    int J = this.NoOfLocalUpdatedCells;
                    BitArray ba = new BitArray(J);

                    for (int j = 0; j < J; j++) {
                        ba[j] = (this.GetRefElementIndex(j) == iKref);
                    }

                    m_RefElementMask[iKref] = new CellMask(this.m_owner, ba);
                }

                return m_RefElementMask[iKref];
            }

            /// <summary>
            /// Cell-Mask of all cells which share the same reference element.
            /// </summary>
            public CellMask GetCells4Refelement(RefElement Kref) {
                var KRefs = this.RefElements;
                int iKref = KRefs.IndexOf(Kref, (A, B) => object.ReferenceEquals(A, B));
                if (iKref < 0 || iKref >= KRefs.Length)
                    throw new ArgumentException();

                if (m_RefElementMask == null) {
                    m_RefElementMask = new CellMask[KRefs.Length];
                }

                if (iKref < 0 || iKref >= KRefs.Length) {
                    throw new ArgumentOutOfRangeException();
                }

                if (m_RefElementMask[iKref] == null) {
                    int J = this.NoOfLocalUpdatedCells;
                    BitArray ba = new BitArray(J);

                    for (int j = 0; j < J; j++) {
                        ba[j] = (this.GetRefElementIndex(j) == iKref);
                    }

                    m_RefElementMask[iKref] = new CellMask(this.m_owner, ba);
                }

                return m_RefElementMask[iKref];
            }

            /// <summary>
            /// ctor
            /// </summary>
            internal CellData(GridData _owner) {
                m_owner = _owner;
            }

            bool m_ContainsNonlinearCells;

            /// <summary>
            /// true, if there is any nonlinear cell on the current MPI process
            /// </summary>
            /// <returns></returns>
            public bool ContainsNonlinearCell() {
                return m_ContainsNonlinearCells;
            }

            /// <summary>
            /// pointer to owner object
            /// </summary>
            private GridData m_owner;

            /// <summary>
            /// local indices of neighbor cells;
            ///  - 1st index: local cell index;
            ///  - 2nd index: enumeration
            /// </summary>
            public int[][] CellNeighbours {
                get {
                    return m_CellNeighbours;
                }
            }

            internal int[][] m_CellNeighbours;

            /// <summary>
            /// global indices of cell neighbors;
            /// </summary>
            internal IEnumerable<GridCommons.Neighbour>[] CellNeighbours_global_tmp;

            
            /// <summary>
            /// true if point <paramref name="pt"/> is inside cell <paramref name="j"/>
            /// </summary>
            public bool IsInCell(double[] pt, int j, double[] pt_Loc = null) {
                int D = m_owner.Grid.SpatialDimension;
                if (pt.Length != D)
                    throw new ArgumentException("length must be equal to spatial dimension", "pt");

                MultidimensionalArray _pt = MultidimensionalArray.CreateWrapper(pt, 1, D);          // point to search for
                MultidimensionalArray _pt_local = AllocHelper(ref pt_Loc, D);  // .. in cell-local coordinate

                bool[] bb = new bool[1];
                this.m_owner.TransformGlobal2Local(_pt, _pt_local, j, bb);

                if (bb[0] == false)
                    return false;

                for (int d = 0; d < D; d++) {
                    if (double.IsInfinity(pt_Loc[d]) || double.IsNaN(pt_Loc[d]))
                        return false;
                }

                return m_owner.Cells.GetRefElement(j).IsWithin(pt_Loc);
            }

            /// <summary>
            /// Helper for safely allocating a wrapped double array.
            /// </summary>
            /// <param name="pt_Loc"></param>
            /// <param name="D"></param>
            /// <returns></returns>
            private static MultidimensionalArray AllocHelper(ref double[] pt_Loc, int D) {
                if (pt_Loc != null) {
                    if (pt_Loc.Length != D)
                        throw new ArgumentException();
                } else {
                    pt_Loc = new double[D];
                }
                MultidimensionalArray _pt_local = MultidimensionalArray.CreateWrapper(
                    pt_Loc, 1, D); // .. in cell-local coordinate
                return _pt_local;
            }

            /// <summary>
            /// Computes, for point <paramref name="pt"/>, the closest point on
            /// the boundary of cell <paramref name="j"/>.
            /// </summary>
            /// <param name="j">local cell index</param>
            /// <param name="pt">
            /// input; some point in global coordinates.
            /// </param>
            /// <param name="pt_Loc">
            /// output, if unequal null: <paramref name="pt"/> in local
            /// coordinates of cell <paramref name="j"/>.
            /// </param>
            /// <param name="closestPoint_global">
            /// output, if unequal null: within the boundary of cell
            /// <paramref name="j"/>, the closest point to
            /// <paramref name="pt"/>.
            /// </param>
            /// <param name="closestPoint_local">
            /// output, if unequal null: <paramref name="closestPoint_global"/>
            /// in global coordinates.
            /// </param>
            /// <returns>
            /// the distance.
            /// </returns>
            public double ClosestPointInCell(double[] pt, int j, double[] pt_Loc = null, double[] closestPoint_global = null, double[] closestPoint_local = null) {
                int D = m_owner.SpatialDimension;
                if (pt.Length != D)
                    throw new ArgumentException("length must be equal to spatial dimension", "pt");

                // point to search for
                MultidimensionalArray _pt = MultidimensionalArray.CreateWrapper(pt, 1, D);

                MultidimensionalArray _pt_local = AllocHelper(ref pt_Loc, D);
                //MultidimensionalArray _closestPoint_local = AllocHelper(ref closestPoint_local, D);

                m_owner.TransformGlobal2Local(_pt, _pt_local, j, null);

                this.GetRefElement(j).ClosestPoint(pt_Loc, closestPoint_local);

                MultidimensionalArray _closestPoint_global = AllocHelper(ref closestPoint_global, D);
                m_owner.TransformLocal2Global(new NodeSet(this.GetRefElement(j), closestPoint_local), _closestPoint_global, j);

                return GenericBlas.L2Dist(closestPoint_global, pt);
            }

            /// <summary>
            /// returns the cell type index (index into <see cref="GridCommons.RefElements"/>).
            /// </summary>
            /// <param name="j"></param>
            /// <returns></returns>
            public int GetRefElementIndex(int j) {
                int info = (int)(InfoFlags[j]);
                Debug.Assert((((int)(CellInfo.RefElementIndex_Mask)) & 1) != 0); // ensure that the mask starts at the least significant bit
                int ret = ((int)(CellInfo.RefElementIndex_Mask)) & info;
                return ret;
            }

            /// <summary>
            /// returns the reference element for cell <paramref name="j"/>
            /// </summary>
            public RefElement GetRefElement(int j) {
                int iKref = GetRefElementIndex(j);
                return this.RefElements[iKref];
            }

            /// <summary>
            /// polynomial interpolation degree of the Reference-to-Global coordinate transformation.
            /// </summary>
            /// <param name="jCell"></param>
            /// <returns></returns>
            public int GetInterpolationDegree(int jCell) {
                return this.GetRefElement(jCell).GetInterpolationDegree(m_owner.Cells.GetCell(jCell).Type);
            }

            /// <summary>
            /// another init ....
            /// </summary>
            internal void Init() {
                using (new FuncTrace()) {
                    int JE = NoOfCells;
                    int J = NoOfLocalUpdatedCells;
                    int D = m_owner.SpatialDimension;

                    this.m_ContainsNonlinearCells = false;

                    MultidimensionalArray Trf = MultidimensionalArray.Create(1, 1, D, D);
                    MultidimensionalArray _Trf = Trf.ExtractSubArrayShallow(0, 0, -1, -1);

                    InfoFlags = new CellInfo[JE];

                    MultidimensionalArray _JacobiDet = MultidimensionalArray.Create(JE);
                    MultidimensionalArray _Transformation = MultidimensionalArray.Create(JE, D, D);
                    MultidimensionalArray _CellCenters = MultidimensionalArray.Create(JE, 1, D);
                    MultidimensionalArray _InverseTransformation = MultidimensionalArray.Create(JE, D, D);
                    MultidimensionalArray Tr = MultidimensionalArray.Create(D, D);
                    MultidimensionalArray invTr = MultidimensionalArray.Create(D, D);

                    IEnumerable<CellType>[] Types = m_owner.Grid.RefElements.Select(x => x.SupportedCellTypes).ToArray();

                    int NegativeJacobianFlag = 0;

                    for (int j = 0; j < JE; j++) {
                        var Cj = this.GetCell(j);

                        int iKref = Types.IndexOfMax(suppTypes => suppTypes.Contains(Cj.Type));
                        if (iKref < 0)
                            throw new NotSupportedException("unknown cell type;");
                        Debug.Assert(iKref == Types.IndexOfMin(suppTypes => suppTypes.Contains(Cj.Type)));
                        
                        var Kref = m_owner.Grid.GetRefElement(iKref);

                        InfoFlags[j] |= (CellInfo.RefElementIndex_Mask & ((CellInfo)iKref));
                        Debug.Assert(((((int)Cj.Type) << 8) & ((int)CellInfo.CellType_Mask)) == (((int)Cj.Type) << 8));
                        InfoFlags[j] |= (CellInfo)(((int)Cj.Type) << 8);
                        
                        // affine-linear cell ?
                        if (Cj.Type.IsLinear()) {
                            InfoFlags[j] |= CellInfo.CellIsAffineLinear;


                            //Kref.JacobianOfTransformation(vtx, Trf, 0, Cj.Type, Cj.TransformationParams);
                            m_owner.EvaluateJacobian(Kref.Center, j, 1, Trf);
                            m_owner.TransformLocal2Global(Kref.Center, j, 1, _CellCenters, j);
                            _Transformation.ExtractSubArrayShallow(j, -1, -1).Set(_Trf);

                            Tr.Set(_Trf);
                            _JacobiDet[j] = Tr.Determinant();
                            if (_JacobiDet[j] <= 0) {
                                NegativeJacobianFlag++;
                            }

                            Tr.InvertTo(invTr);
                            _InverseTransformation.ExtractSubArrayShallow(j, -1, -1).SetMatrix(invTr);

                        } else {
                            // metrics are non-constant
                            // ++++++++++++++++++++++++

                            // mark to prevent miss-use
                            _JacobiDet[j] = double.NaN;
                            _InverseTransformation.SetAll(double.NaN, j, -1, -1);
                            _Transformation.SetAll(double.NaN, j, -1, -1);
                            _CellCenters.SetAll(double.NaN, j, -1, -1);

                            m_ContainsNonlinearCells = true;

                            {
                                // test that the Jacobian determinant is positive definite

                                int deg = Kref.GetInterpolationDegree(Cj.Type);
                                if (deg > 1)
                                    deg--;
                                deg *= D;

                                NodeSet TestNodes = Kref.GetQuadratureRule(2 * deg).Nodes;
                                MultidimensionalArray tmpTrf = MultidimensionalArray.Create(1, TestNodes.NoOfNodes, D, D);
                                m_owner.EvaluateJacobian(TestNodes, j, 1, tmpTrf);

                                for (int n = 0; n < TestNodes.NoOfNodes; n++) {
                                    double detJac = tmpTrf.ExtractSubArrayShallow(0, n, -1, -1).Determinant();
                                    if (detJac <= 0) {
                                        NegativeJacobianFlag++;
                                        _JacobiDet[j] = detJac;
                                        break;
                                    }
                                }

                            }
                        }
                    }

                    if (NegativeJacobianFlag > 0) {
                        for (int j = 0; j < JE; j++) {
                            if (_JacobiDet[j] < 0)
                                throw new ArithmeticException(string.Format("Found {3} cell(s) (of {4} tested) with non-positive Jacobian determinant. First problem occurred at cell #{0}, global ID {1}, Jacobian determinant is {2}.", j, this.GetCell(j).GlobalID, _JacobiDet[j], NegativeJacobianFlag, JE));
                        }
                    }

                    this.JacobiDet = _JacobiDet;
                    this.Transformation = _Transformation;
                    this.CellCenter = _CellCenters.ResizeShallow(JE, D);
                    this.InverseTransformation = _InverseTransformation;

                }
            }

          

            /// <summary>
            /// see <see cref="CellInfo"/>
            /// </summary>
            public CellInfo[] InfoFlags {
                get;
                private set;
            }

            /// <summary>
            /// Aids the vectorization of various code parts.
            /// </summary>
            /// <param name="mask">
            /// masks which properties of the cell information (see
            /// <see cref="InfoFlags"/>) should be considered.
            /// </param>
            /// <param name="j0">start index.</param>
            /// <returns>
            /// the number of consecutive cells after cell
            /// <paramref name="j0"/>, which share the same information flags,
            /// or <paramref name="Lmax"/>, whichever is lower.
            /// </returns>
            /// <param name="Lmax">
            /// upper limit for the return value of this function
            /// </param>
            public int GetNoOfSimilarConsecutiveCells(CellInfo mask, int j0, int Lmax) {
                // this will be optimized some time, if necessary...
                uint Info_at_j0 = ((uint)InfoFlags[j0]) & (uint)mask;
                int R = 1;
                int J = this.NoOfLocalUpdatedCells;
                for (int j = j0 + 1; j < J; j++) {
                    uint Info_at_j = ((uint)InfoFlags[j0]) & (uint)mask;
                    if (Info_at_j != Info_at_j0 || R >= Lmax)
                        return R;
                    R++;
                }
                return R;
            }

            /// <summary>
            /// true if cell <paramref name="j"/> is affine-linear.
            /// </summary>
            public bool IsCellAffineLinear(int j) {
                return ((InfoFlags[j] & CellInfo.CellIsAffineLinear) != 0);
            }

            /// <summary>
            /// For affine-linear cells,
            /// the absolute value of the (Jacobi) determinant of the 
            /// transformation from local cell coordinate system to global
            /// coordinate system.
            /// (see <see cref="EvaluateJacobian(NodeSet, int, int, MultidimensionalArray)"/>
            /// 1st index: local cell index;
            /// </summary>
            public MultidimensionalArray JacobiDet {
                get;
                private set;
            }

            /// <summary>
            /// The minimal Euclidean distance between two vertices for each cell;
            /// (Can be used to compute the CFL number);
            /// 1st index: local cell index;
            /// </summary>
            public MultidimensionalArray h_min {
                get;
                private set;
            }

            /// <summary>
            /// <see cref="h_minGlobal"/>
            /// </summary>
            double m_h_minGlobal = double.NaN;

            /// <summary>
            /// Minimum over all entries in <see cref="h_min"/>, over all MPI processes.
            /// processes
            /// </summary>
            public double h_minGlobal {
                get {
                    return m_h_minGlobal;
                }
            }

            /// <summary>
            /// <see cref="h_minGlobal"/>
            /// </summary>
            double m_h_maxGlobal = double.NaN;

            /// <summary>
            /// Maximum over all entries in <see cref="h_max"/>, over all MPI processes.
            /// </summary>
            public double h_maxGlobal {
                get {
                    return m_h_maxGlobal;
                }
            }

            /// <summary>
            /// The maximal Euclidean distance between two vertices for each cell;
            /// (Can be used to compute the CFL number);
            /// 1st index: local cell index;
            /// </summary>
            public MultidimensionalArray h_max {
                get;
                private set;
            }

            /// <summary>
            /// the (<see cref="SpatialDimension"/>-1) - dimensional measure of
            /// the cell boundary;
            /// - index: local cell index;
            /// </summary>
            public double[] CellSurfaceArea;

            /// <summary>
            /// for each cell <em>j</em> this is
            /// \f[ 
            /// c_j = frac{Area(\partial K_j \\ \partial \Omega) \cdot 0.5 + Area(\partial K_j \cap \partial \Omega)}{Volume(K_j)},
            /// \f]
            /// where \f$ K_j \f$ is the cell an
            /// \f$ \Omega \f$ is the whole computational
            /// domain.
            /// </summary>
            /// <remarks>
            /// Needed for Interior Penalty method, for the Poisson equation  
            /// penalty parameter according to:
            /// An explicit expression for the penalty parameter of the
            /// interior penalty method,
            /// K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
            /// look at formula (7).
            /// </remarks>
            public MultidimensionalArray cj;

            /// <summary>
            /// Alias for cj
            /// </summary>
            public MultidimensionalArray PenaltyLengthScales
            {
                get { return cj ; }
            }



            /// <summary>
            /// Returns the volume (to be more exact: the
            /// <see cref="SpatialDimension"/> - dimensional measure) of the
            /// cell <paramref name="j"/>;
            /// </summary>
            /// <param name="j">local cell index</param>
            /// <returns></returns>
            public double GetCellVolume(int j) {
                if (IsCellAffineLinear(j)) {
                    double vol = this.GetRefElement(j).Volume;
                    return vol * JacobiDet[j];
                } else {
                    double D = m_owner.SpatialDimension;
                    double deg = this.GetInterpolationDegree(j);
                    var Cj = this.GetCell(j);

                    if (deg > 1)
                        deg -= 1;
                    var qr = m_owner.Grid.GetRefElement(this.GetRefElementIndex(j)).GetQuadratureRule((int)(deg * D));
                    double vol = 0;
                    NodeSet RefVertices = qr.Nodes;
                    MultidimensionalArray JacobianDetOut = this.m_owner.JacobianDeterminat.GetValue_Cell(RefVertices, j, 1);
                    for (int n = qr.Weights.GetLength(0) - 1; n >= 0; n--) {
                        vol += qr.Weights[n] * (double)JacobianDetOut[0, n];
                    }
                    return vol; //check for absolute value
                }
            }

            /// <summary>
            /// Mapping from cells to vertices/nodes of the grid (stored in
            /// <see cref="IVertexData.Coordinates"/>) 
            /// - content: indices into <see cref="IVertexData.Coordinates"/> 
            /// - 1st index: local cell index (externals included) 
            /// - 2nd index: cell vertex index
            /// </summary>
            public int[][] CellVertices {
                get;
                internal set;
            }

            /// <summary>
            /// Computes the bounding box of cell <paramref name="j"/>.
            /// </summary>
            /// <param name="j">local cell index.</param>
            /// <param name="bb">
            /// on exit, the bounding box of cell j.
            /// </param>
            public void GetCellBoundingBox(int j, BoundingBox bb) {
                int D = bb.D;
                if (bb.D != m_owner.SpatialDimension)
                    throw new ArgumentException("wrong dimension of bounding box.");
                bb.Clear();

                Cell Cj = this.GetCell(j);
                RefElement Kref = this.GetRefElement(j);
                NodeSet verticesLoc = Kref.GetInterpolationNodes(Cj.Type);

                MultidimensionalArray verticesGlob = MultidimensionalArray.Create(1, verticesLoc.GetLength(0), verticesLoc.GetLength(1));
                m_owner.TransformLocal2Global(verticesLoc, j, 1, verticesGlob, 0);
                bb.AddPoints(verticesGlob
                    .ExtractSubArrayShallow(0, -1, -1));
            }

            /// <summary>
            /// For all affine-linear cells, the 
            /// linear part of the affine-linear transformation
            /// from the local coordinate 
            /// system of some cell to
            /// the global coordinate system, or Jacobi-matrix.
            /// </summary>
            /// <remarks>
            /// Indices are defined as follows:
            /// <list type="bullet">
            ///   <item>1st index: local cell index (locally updated and external cells);</item>
            ///   <item>2nd index: matrix row index;</item>
            ///   <item>3rd index: matrix column index;</item>
            /// </list>
            /// </remarks>
            public MultidimensionalArray Transformation {
                get;
                private set;
            }

            /// <summary>
            /// For all affine-linear cells, the 
            /// affine part of the affine-linear transformation
            /// from the local coordinate 
            /// system of some cell to
            /// the global coordinate system, or the cell-center.
            /// </summary>
            /// <remarks>
            /// Indices are defined as follows:
            /// <list type="bullet">
            ///   <item>1st index: local cell index (locally updated and external cells);</item>
            ///   <item>2nd index: spatial dimension</item>
            /// </list>
            /// </remarks>
            public MultidimensionalArray CellCenter;


            /// <summary>
            /// inverse matrices to <see cref="Transformation"/>
            /// </summary>
            /// <remarks>
            /// Indices are defined as follows:
            /// <list type="bullet">
            ///   <item>1st index: local cell index (locally updated and external cells);</item>
            ///   <item>2nd index: matrix row index;</item>
            ///   <item>3rd index: matrix column index;</item>
            /// </list>
            /// </remarks>
            public MultidimensionalArray InverseTransformation {
                get;
                private set;
            }

            /// <summary>
            /// initializes <see cref="cj"/>
            /// </summary>
            internal void InitializeCk() {
                using (new FuncTrace()) {
                    var edgeDat = m_owner.Edges;

                    int Jtot = this.NoOfCells;
                    int Jloc = this.NoOfLocalUpdatedCells;
                    this.cj = MultidimensionalArray.Create(Jtot);
                    double[] bndArea = new double[Jtot];

                    int E = edgeDat.Count;
                    for (int e = 0; e < E; e++) {
                        int Cell1 = edgeDat.CellIndices[e, 0];
                        int Cell2 = edgeDat.CellIndices[e, 1];

                        if (Cell2 < 0)
                            bndArea[Cell1] += edgeDat.GetEdgeArea(e);
                    }
                    for (int j = 0; j < Jloc; j++) {
                        cj[j] = (CellSurfaceArea[j] + bndArea[j]) * 0.5 / GetCellVolume(j);
                    }



                    // for external cells, the result will be wrong, because edges of external cells
                    // are only tracked on this processor if they bound to an internal cell too
                    // -> send volume of external cells via MPI

                    for (int j = Jloc; j < Jtot; j++)
                        this.cj[j] = double.NaN;

                    this.cj.MPIExchange(this.m_owner);

                }
            }

            /// <summary>
            /// initializes <see cref="CellSurfaceArea"/>
            /// </summary>
            internal void InitializeCellSurfaceArea() {
                using (new FuncTrace()) {

                    var edgDat = this.m_owner.Edges;
                    int Jtot = this.NoOfCells;
                    int Jloc = this.NoOfLocalUpdatedCells;
                    this.CellSurfaceArea = new double[Jtot];

                    int E = edgDat.Count;
                    for (int e = 0; e < E; e++) {
                        this.CellSurfaceArea[edgDat.CellIndices[e, 0]] += edgDat.GetEdgeArea(e);

                        int Cell2 = edgDat.CellIndices[e, 1];
                        if (Cell2 >= 0 && Cell2 < Jloc)
                            this.CellSurfaceArea[Cell2] += edgDat.GetEdgeArea(e);

                    }

                    // for external cells, the result will be wrong, because edges of external cells
                    // are only tracked on this processor if they bound to an internal cell too
                    // -> send volume of external cells via MPI

                    for (int j = Jloc; j < Jtot; j++)
                        this.CellSurfaceArea[j] = double.NaN;

                    this.CellSurfaceArea.MPIExchange(this.m_owner);
                }
            }

            /// <summary>
            /// Number of locally updated cells - the cells which are computed on
            /// this processor (in contrast, see <see cref="NoOfExternalCells"/>);
            /// </summary>
            public int NoOfLocalUpdatedCells {
                get {
                    return this.m_owner.Grid.NoOfUpdateCells;
                }
            }

            /// <summary>
            /// Number of locally stored external cells - no computations are carried out for
            /// that cells, but their values are needed.
            /// </summary>
            public int NoOfExternalCells {
                get {
                    return m_owner.m_Parallel.GlobalIndicesExternalCells.Length;
                }
            }

            /// <summary>
            /// <see cref="NoOfExternalCells"/> plus <see cref="NoOfLocalUpdatedCells"/>;
            /// </summary>
            public int NoOfCells {
                get {
                    return (NoOfLocalUpdatedCells + NoOfExternalCells);
                }
            }

            /// <summary>
            /// returns the <see cref="Cell"/>-object for a cell index <paramref name="j"/>,
            /// </summary>
            /// <param name="j">
            /// local cell index, (can also be in the range of external/ghost cells)
            /// </param>
            public Cell GetCell(int j) {
                Debug.Assert(j >= 0);
                Debug.Assert(j < NoOfCells);
                int J = NoOfLocalUpdatedCells;
                return ((j < J) ? m_owner.m_Grid.Cells[j] : m_owner.m_Parallel.ExternalCells[j - J]);
            }
            
            /// <summary>
            /// Cell type for cell <paramref name="jCell"/>.
            /// </summary>
            public CellType GetCellType(int j) {
                int iType = (((int)(InfoFlags[j])) & ((int)(CellInfo.CellType_Mask))) >> 8;
                Debug.Assert(iType == ((int)(GetCell(j).Type)));
                return ((CellType)iType);
            }


            /// <summary>
            /// fills <see cref="h_min"/>, <see cref="h_max"/>;
            /// </summary>
            internal void Initialize_h() {
                using (new FuncTrace()) {
                    var GridSimplices = m_owner.Grid.RefElements;
                    int L = GridSimplices.Length;
                    NodeSet[] _vertices = new NodeSet[L];
                    MultidimensionalArray[] _vertGlob = new MultidimensionalArray[L];
                    for (int l = 0; l < L; l++) {
                        _vertices[l] = GridSimplices[l].Vertices;
                        _vertGlob[l] = MultidimensionalArray.Create(1, _vertices[l].GetLength(0), _vertices[l].GetLength(1));
                    }

                    int JE = NoOfCells;
                    int J = NoOfLocalUpdatedCells;
                    double __m_h_minGlobal = double.MaxValue;
                    double __m_h_max_Global = 0;
                    h_min = MultidimensionalArray.Create(JE);
                    h_max = MultidimensionalArray.Create(JE);
                    int D = m_owner.SpatialDimension;

                    double[] vk = new double[D];
                    double[] vi = new double[D];

                    for (int j = 0; j < J; j++) {

                        int l = this.GetRefElementIndex(j);
                        var vertices = _vertices[l];
                        var vertGlob = _vertGlob[l];

                        int N = vertices.GetLength(0);

                        m_owner.TransformLocal2Global(vertices, j, 1, vertGlob, 0);

                        double dist_min = double.MaxValue;
                        double dist_max = 0;
                        for (int i = 0; i < N; i++) {
                            for (int d = 0; d < D; d++)
                                vi[d] = vertGlob[0, i, d];

                            for (int k = i + 1; k < N; k++) {
                                for (int d = 0; d < D; d++)
                                    vk[d] = vertGlob[0, k, d];

                                double dist_ik = 0;
                                for (int d = 0; d < D; d++) {
                                    double hh = vi[d] - vk[d];
                                    dist_ik += hh * hh;
                                }
                                if (dist_ik < dist_min)
                                    dist_min = dist_ik;
                                if (dist_ik > dist_max)
                                    dist_max = dist_ik;
                            }
                        }

                        h_min[j] = Math.Sqrt(dist_min);
                        h_max[j] = Math.Sqrt(dist_max);

                        if (h_min[j] <= 0) {
                            throw new ArithmeticException("Error in grid. Found a degenerate cell (i.e. distance between tow vertices is 0.0)"
                                + "(GlobalId = " + GetCell(j).GlobalID + " );");
                        }

                        __m_h_minGlobal = Math.Min(h_min[j], __m_h_minGlobal);
                        __m_h_max_Global = Math.Max(h_max[j], __m_h_max_Global);
                    }

                    h_min.MPIExchange(m_owner);
                    h_max.MPIExchange(m_owner);
                    m_h_minGlobal = __m_h_minGlobal.MPIMin();
                    m_h_maxGlobal = __m_h_max_Global.MPIMax();
                }
            }

            /// <summary>
            /// As defined by <see cref="iLogicalCells"/>, the Global Id for the <paramref name="j"/>-th cell.
            /// </summary>
            public long GetGlobalID(int j) {
                return this.GetCell(j).GlobalID;
            }

            /// <summary>
            /// Which edges (see <see cref="EdgeData.CellIndices"/>) bound to
            /// which cells? <br/>
            /// - 1st index: local cell index <em>j</em>, only local updated<br/>
            /// - 2nd index: collection, order is arbitrary; <br/>
            /// - content <em>e</em>: 
            ///   If <em>e</em> is positive, then cell <em>j</em> is the first
            ///   (IN) cell of edge <em>e - 1</em>. Otherwise, if <em>e</em> is
            ///   negative, then cell <em>j</em> is the second (OUT) cell of edge
            ///   <em>-e - 1</em>.
            /// </summary>
            /// <remarks>
            /// Note: the second index does NOT correlate with the face index
            /// of the cell. This is because, in the case of hanging nodes, the
            /// number of edges that bound to one cell is not equal to the
            /// number of faces, i.e., more than one edge is associated with
            /// one face.
            /// </remarks>
            public int[][] Cells2Edges {
                get;
                internal set;
            }

            /// <summary>
            /// Since each cell is already elementary (i.e. it maps to one reference element)
            /// this is not required and therefore, equal to null.
            /// </summary>
            public int[][] AggregateCellToParts {
                get {
                    return null;
                }
            }
        }
    }
}
