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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// Contains extended information about the computational grid 
    /// (stored in <see cref="GridCommons"/>-objects), which are not stored on
    /// disk but computed at run-time, i.e. load-time. This information is
    /// essential for the Discontinuous Galerkin algorithms.
    /// </summary>
    sealed public partial class GridData : IGridData {

        private static ILog Logger = LogManager.GetLogger(typeof(GridData));

        private GridCommons m_Grid;

        /// <summary>
        /// The grid for which information is provided
        /// </summary>
        public GridCommons Grid {
            get {
                return m_Grid;
            }
        }

        /// <summary>
        /// The grid for which information is provided
        /// </summary>
        IGrid IGridData.Grid {
            get {
                return m_Grid;
            }
        }

        /// <summary>
        /// Identification of the grid in the BoSSS database, 
        /// equal to the <see cref="BoSSS.Foundation.IO.IDatabaseEntityInfo{T}.ID"/>.
        /// </summary>
        public Guid GridID {
            get {
                return Grid.ID;
            }
        }

        /// <summary>
        /// Equal to <see cref="GridCommons.EdgeTagNames"/>.
        /// </summary>
        public IDictionary<byte, string> EdgeTagNames {
            get {
                return Grid.EdgeTagNames;
            }
        }

        /// <summary>
        /// MPI process rank (within world communicator)
        /// </summary>
        public int MpiRank {
            get {
                return this.CellPartitioning.MpiRank;
            }
        }

        /// <summary>
        /// MPI world communicator size 
        /// </summary>
        public int MpiSize {
            get {
                return this.CellPartitioning.MpiSize;
            }
        }

        /// <summary>
        /// constructor 
        /// </summary>
        public GridData(GridCommons grd) {

            if(grd.RefElements.Length != 1)
                throw new ApplicationException("Currently only grids with _one_ RefElement are supported!!!");

            {
                // test that each face reference element, 
                //     for each cell/volume reference element, 
                //         is contained in the list of edge reference elements, 
                // by reference-equality!

                RefElement[] Edge_KrefS = grd.EdgeRefElements;
                foreach(var Cell_Kref in grd.RefElements) {
                    RefElement CellFaceKref = Cell_Kref.FaceRefElement;

                    int FoundCount = 0;
                    for(int i = 0; i < Edge_KrefS.Count(); i++) {
                        if(object.ReferenceEquals(CellFaceKref, Edge_KrefS.ElementAt(i)))
                            FoundCount++;
                    }
                    if(FoundCount != 1) {
                        throw new ArgumentException("Insane reference element structure.");
                    }

                    //Debug.Assert(Edge_KrefS.Where(A => object.ReferenceEquals(Cell_Kref.FaceRefElement, A)).Count() == 1);
                }
            }
            m_Grid = grd;


            using(var ft = new ilPSP.Tracing.FuncTrace()) {
                // Test Grid
                // ---------
                //grd.TestJacobianForAllCells();
                CheckEdgeTagNames();
                Grid.CheckGridForNANorINF();
                Grid.CheckCellTypes();
                int myRank = this.MpiRank;
                if (m_Grid.NoOfUpdateCells <= 0)
                    throw new ApplicationException("grid contains no cells on processor " + myRank + ";");

                // start init
                // ----------
                m_Edges = new EdgeData(this);
                m_Cells = new CellData(this);
                m_VerticeData = new VertexData(this);
                m_Parallel = new Parallelization(this);

                // caches
                // ------
                int D = this.m_Grid.SpatialDimension;
                this.Jacobian = new Caching.CacheLogicImplBy_CNs(this, this.EvaluateJacobian,
                    (Len, NoOfNodes) => new int[] { Len, NoOfNodes, D, D });
                this.AdjungateJacobian = new Caching.CacheLogicImplBy_CNs(this, this.EvaluateAdjungateJacobian,
                    (Len, NoOfNodes) => new int[] { Len, NoOfNodes, D, D });
                this.InverseJacobian = new Caching.CacheLogicImplBy_CNs(this, this.EvaluateInverseJacobian,
                    (Len, NoOfNodes) => new int[] { Len, NoOfNodes, D, D });
                this.JacobianDeterminat = new Caching.CacheLogicImplBy_CNs(this, this.EvaluateJacobianDeterminat,
                    (Len, NoOfNodes) => new int[] { Len, NoOfNodes});

                // cell init
                // ---------
                ParallelSetup();
                m_Cells.Init();

                // collect edges
                // -------------

                
                m_Edges.CollectEdges();
                m_Edges.DetermineEdgeTrafo();
                m_Edges.CollectBoundaryEdges();
                m_Edges.SetEdgeTags();
                m_Edges.NegogiateNeighbourship();

                m_Edges.InitCells2Edges();
                m_Edges.FinalizeAssembly();
                m_Cells.CellNeighbours_global_tmp = null;
                this.m_BcCells_tmp = null;

                // some edges metrics
                // ------------------
                //Logger.Info("edge metrics...");

                Edges.ComputeSqrtGramian();
                Edges.Initialize_h_Edge();
                Edges.InitNormals();

                // some other cell metrics
                // -----------------------

                Cells.InitializeCellSurfaceArea();
                Cells.InitializeCk();
                Cells.Initialize_h();

                InitNoOfCellsPerRefElement();


                // vertices
                // --------
                MultidimensionalArray Vertices_Coordinates;
                int[][] Cells_CellVertices;
                VertexData.CollectVertices(this, out Cells_CellVertices, out Vertices_Coordinates, Cells.h_min.To1DArray(), true);
                Vertices.Coordinates = Vertices_Coordinates;
                Cells.CellVertices = Cells_CellVertices;

                this.Vertices.Init_VerticeToCell();
                this.Vertices.Init_PeriodicPairings();
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                this.Vertices.NegogiateSharing();
                
                // bla
                // ===

                DefineBoundaryMasks();

                // bla bla
                // =======

                this.ChefBasis = new GridData._BasisData(this);

                // a good point for garbage collection
                // -----------------------------------
                ft.LogMemoryStat();
                Logger.Info("calling garbage collector...");
                GC.Collect();
                ft.LogMemoryStat();
            }
        }


        /// <summary>
        /// Clears all internal references for this object, to make sure that any attempt to use it leads to an exception.
        /// </summary>
        public void Invalidate() {
            this.m_Cells = null;
            this.m_Edges = null;
            this.m_GlobalNodes = null;
            this.m_Grid = null;
            this.m_Parallel = null;
            this.m_VerticeData = null;
            this.m_Edges = null;
        }



        private void InitNoOfCellsPerRefElement() {
            int JE = Cells.Count;
            int J = Cells.NoOfLocalUpdatedCells;

            var KrefS = Grid.RefElements;

            m_NoOfCellsPerRefElement_Local = new int[KrefS.Length];
            m_NoOCellsPerRefElement_External = new int[KrefS.Length];

            // Beware of correct incrementation
            int j;
            for (j = 0; j < J; j++)
                m_NoOfCellsPerRefElement_Local[Cells.GetRefElementIndex(j)]++;
            for (; j < JE; j++)
                m_NoOCellsPerRefElement_External[Cells.GetRefElementIndex(j)]++;
        }

        int[] m_NoOfCellsPerRefElement_Local;

        int[] m_NoOCellsPerRefElement_External;

        /// <summary>
        /// number of cells
        /// of type <paramref name="r"/>
        /// </summary>
        /// <param name="Ext">number of <paramref name="r"/>-type cells in external/ghost cells</param>
        /// <param name="Loc">number of <paramref name="r"/>-type cells in locally updated cells</param>
        /// <param name="r">cell type</param>
        public void GetLocalNoOfCellsRefElement(RefElement r, out int Loc, out int Ext) {
            int iKref = Array.IndexOf(Grid.RefElements, r);
            GetLocalNoOfCellsRefElement(iKref, out Loc, out Ext);
        }

        /// <summary>
        /// number of cells (including ghost/external cells)
        /// of type <see cref="GridCommons.RefElements"/>[<paramref name="iKref"/>].
        /// </summary>
        /// <param name="Ext">
        /// on exit, number of external/ghost cells of specified type
        /// </param>
        /// <param name="Loc">
        /// on exit, number of locally updated cells of specified type
        /// </param>
        /// <param name="iKref">
        /// reference element index
        /// </param>
        public void GetLocalNoOfCellsRefElement(int iKref, out int Loc, out int Ext) {
            Loc = m_NoOfCellsPerRefElement_Local[iKref];
            Ext = m_NoOCellsPerRefElement_External[iKref];
        }

        /// <summary>
        /// transforms vertices from the local coordinate system of cells <paramref name="jCell"/>
        /// to global coordinates;
        /// </summary>
        /// <param name="LocalVerticesIn">
        /// Input; vertices in the local coordinate system of a cell;
        /// <list type="bullet">
        ///   <item>1st index: vertex index;</item>
        ///   <item>2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and 0,1,2 for 3D;</item>
        /// </list>
        /// </param>
        /// <param name="GlobalVerticesOut">
        /// Output; the vertices form <paramref name="LocalVerticesIn"/>, transformed to global
        /// coordinates;
        /// <list type="bullet">
        ///   <item>
        ///     1st index: vertex index, corresponds with the 1st index of <paramref name="LocalVerticesIn"/>;
        ///   </item>
        ///   <item>
        ///     2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and 0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="jCell">local cell index of the cell to transform</param>
        public void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, MultidimensionalArray GlobalVerticesOut, int jCell) {
            int N = LocalVerticesIn.GetLength(0);
            int D = SpatialDimension;
            if (GlobalVerticesOut.GetLength(1) != D)
                throw new ArgumentException("wrong spatial dimension of GlobalVerticesOut");
            if (LocalVerticesIn.GetLength(1) != D)
                throw new ArgumentException("wrong spatial dimension of LocalVerticesIn");
            if (LocalVerticesIn.GetLength(0) != GlobalVerticesOut.GetLength(0))
                throw new ArgumentException("mismatch in number of vertices per cell.");


            var Cl = m_Cells.GetCell(jCell);

            var _GlobalVerticesOut = MultidimensionalArray.Create(1, N, D);
            this.TransformLocal2Global(LocalVerticesIn, jCell, 1, _GlobalVerticesOut, 0);
            GlobalVerticesOut.Set(_GlobalVerticesOut.ExtractSubArrayShallow(0, -1, -1));
        }

        /// <summary>
        /// transforms vertices from the local coordinate system of cells
        /// <paramref name="j0"/> to <paramref name="j0"/>+<paramref name="Len"/>-1
        /// to global coordinates;
        /// </summary>
        /// <param name="LocalVerticesIn">
        /// Input; vertices in the local coordinate system of a cell;
        /// <list type="bullet">
        ///   <item>1st index: vertex index;</item>
        ///   <item>
        ///   2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///   0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="GlobalVerticesOut">
        /// Output; the vertices form <paramref name="LocalVerticesIn"/>,
        /// transformed to global coordinates;
        /// <list type="bullet">
        ///   <item>
        ///     1st index: local cell index minus <paramref name="j0"/>, in the
        ///     range of 0 to <paramref name="Len"/>-1;
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="LocalVerticesIn"/>;
        ///   </item>
        ///   <item>
        ///     3rd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///     0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="j0">local cell index of the first cell to transform</param>
        /// <param name="Len">Number of cells to transform</param>
        /// <param name="OutArrayOffset">
        /// an offset into the first index of <paramref name="GlobalVerticesOut"/>;
        /// </param>
        public void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, int j0, int Len, MultidimensionalArray GlobalVerticesOut, int OutArrayOffset) {
            int NoOfNodes = LocalVerticesIn.GetLength(0);
            int D = SpatialDimension;
            if(GlobalVerticesOut.GetLength(2) != D)
                throw new ArgumentException("wrong spatial dimension of GlobalVerticesOut");
            if(LocalVerticesIn.GetLength(1) != D)
                throw new ArgumentException("wrong spatial dimension of LocalVerticesIn");
            if(LocalVerticesIn.GetLength(0) != GlobalVerticesOut.GetLength(1))
                throw new ArgumentException("mismatch in number of vertices per cell.");
            if(Len > GlobalVerticesOut.GetLength(0) - OutArrayOffset) {
                throw new ArgumentException("first dimension of GlobalVerticesOut to short.");
            }

            bool AffineLinear = this.Cells.IsCellAffineLinear(j0);
#if DEBUG
            for(int j = 1; j < Len; j++) {
                Debug.Assert(this.Cells.IsCellAffineLinear(j0 + j) == AffineLinear);
            }
#endif
            if(GlobalVerticesOut.GetLength(0) > Len || OutArrayOffset != 0)
                GlobalVerticesOut = GlobalVerticesOut.ExtractSubArrayShallow(new int[] { OutArrayOffset, 0, 0 }, new int[] { OutArrayOffset + Len - 1, NoOfNodes - 1, D - 1 });

            if (AffineLinear && this.Cells.Transformation != null) {
                var Trf = this.Cells.Transformation.ExtractSubArrayShallow(new int[] { j0, 0, 0 }, new int[] { j0 + Len - 1, D - 1, D - 1 });
                var Aff = this.Cells.CellCenter;//.ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, D - 1 });

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < NoOfNodes; k++) {
                        for (int d = 0; d < D; d++) {
                            GlobalVerticesOut[j, k, d] = Aff[j + j0, d];
                        }
                    }
                }

                GlobalVerticesOut.Multiply(1.0, Trf, LocalVerticesIn, 1.0, ref mp_jkd_jde_ke);
            } else {

                // work with polynomial caches

                var Krefs = this.Cells.RefElements;
                NodeSet[] _LocalVerticesIn = new NodeSet[Krefs.Length];        

                // loop over cells ...
                for (int j = 0; j < Len; j++) {
                    int jCell = j0 + j;
                    var Cl = m_Cells.GetCell(jCell);
                    int iKref = m_Cells.GetRefElementIndex(jCell);
                    var Kref = Krefs[iKref];

                    if(_LocalVerticesIn[iKref] == null) {
                        if (LocalVerticesIn is NodeSet) {
                            _LocalVerticesIn[iKref] = (NodeSet)LocalVerticesIn; 
                        } else {
                            _LocalVerticesIn[iKref] = new NodeSet(Kref, LocalVerticesIn); // convert to node set
                        }
                    }

                    Debug.Assert(object.ReferenceEquals(Kref, _LocalVerticesIn[iKref].RefElement));

                    PolynomialList polys = Kref.GetInterpolationPolynomials(Cl.Type);
                    MultidimensionalArray polyVals = polys.Values.GetValues(_LocalVerticesIn[iKref]);

                    for (int d = 0; d < D; d++) {
                        GlobalVerticesOut.ExtractSubArrayShallow(j, -1, d)
                            .Multiply(1.0, polyVals, Cl.TransformationParams.ExtractSubArrayShallow(-1, d), 0.0, ref mp_k_kn_n);
                    }
                }
            }
        }

        static MultidimensionalArray.MultiplyProgram mp_k_kn_n = MultidimensionalArray.MultiplyProgram.Compile("k", "kn", "n");
        static MultidimensionalArray.MultiplyProgram mp_jkd_jde_ke = MultidimensionalArray.MultiplyProgram.Compile("jkd", "jde", "ke");
        
        Caching.CacheLogic_CNs m_GlobalNodes;

        /// <summary>
        /// Cached transformation of node sets to global coordinates.
        /// </summary>
        public Caching.CacheLogic_CNs GlobalNodes {
            get {
                if(m_GlobalNodes == null) {
                    m_GlobalNodes = new Caching.CacheLogicImplBy_CNs(this,
                        this.TransformLocal2Global,
                        (Len, NoNodes) => new int[] { Len, NoNodes, this.SpatialDimension });
                }
                return m_GlobalNodes;
            }
        }

        /// <summary>
        /// Jacobian of transformation from reference to physical space, \f$ (\nabla T_j) \f$.
        /// </summary>
        public Caching.CacheLogicImplBy_CNs Jacobian {
            get;
            private set;
        }

        
        void EvaluateJacobian(NodeSet NS, int j0, int Len, MultidimensionalArray output) {

            // init + argcheck
            // ===============

            RefElement Kref = this.Cells.GetRefElement(j0);
            CellType ct = this.Cells.GetCell(j0).Type;
            int D = this.SpatialDimension;

            if(!object.ReferenceEquals(NS.RefElement, Kref))
                throw new ArgumentException("Mismatch between ref element for node set and cell.");
#if DEBUG
            for(int i = 0; i < Len; i++) {
                int jCell = i + j0;
                if(this.Cells.GetCell(jCell).Type != ct)
                    throw new NotSupportedException("One evaluation junk may contain only type of cells.");
            }
#endif

            PolynomialList[] Deriv = Kref.GetInterpolationPolynomials1stDeriv(ct);
            Debug.Assert(Deriv.Length == D);

            Debug.Assert(output.Dimension == 4);
            Debug.Assert(output.GetLength(0) == Len);
            Debug.Assert(output.GetLength(1) == NS.NoOfNodes);
            Debug.Assert(output.GetLength(2) == D);
            Debug.Assert(output.GetLength(3) == D);

            int K = NS.NoOfNodes;

            if(ct.IsLinear() && this.Cells.Transformation != null) {
                // evaluate linear
                // ===============

                MultidimensionalArray Trf = this.Cells.Transformation;


                for(int i = 0; i < Len; i++) { // loop over cells...
                    int jCell = i + j0;
                    MultidimensionalArray Trf_jCell = Trf.ExtractSubArrayShallow(jCell, -1, -1);

                    for(int k = 0; k < K; k++) { // loop over nodes...
                        output.SetSubArray(Trf_jCell, i, k, -1, -1); // Jacobian is equal for all nodes.
                    }
                }


            } else {

                MultidimensionalArray[] DerivEval = new MultidimensionalArray[D];
                for(int d = 0; d < D; d++) {
                    DerivEval[d] = Deriv[d].Values.GetValues(NS);
                }

                // evaluate nonlinear
                // ==================

                for(int i = 0; i < Len; i++) { // loop over cells...
                    int jCell = i + j0;

                    Cell Cj = this.Cells.GetCell(jCell);
                    MultidimensionalArray TrafoParams = Cj.TransformationParams;
                    Debug.Assert(TrafoParams.Dimension == 2);
                    Debug.Assert(TrafoParams.GetLength(1) == D);

                    for(int d1 = 0; d1 < D; d1++) {
                        Debug.Assert(TrafoParams.GetLength(0) == Deriv[d1].Count);
                        MultidimensionalArray JacobiCol = output.ExtractSubArrayShallow(i, -1, -1, d1);
                        JacobiCol.Multiply(1.0, DerivEval[d1], TrafoParams, 0.0, "kd", "kn", "nd");
                    }

                }
            }
        }

        /// <summary>
        /// Adjungate of the Jacobian of the reference-to-physical coordinate transformation, 
        /// \f$ \mathrm{Adj}( \nabla T_j ) =  \determinant{ \nabla T_j } ( \nabla T_j )^{-1} \f$.
        /// </summary>
        public Caching.CacheLogicImplBy_CNs AdjungateJacobian {
            get;
            private set;
        }

        void EvaluateAdjungateJacobian(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
            int K = NS.NoOfNodes;
            int D = this.SpatialDimension;

            MultidimensionalArray Jacobian = this.Jacobian.GetValue_Cell(NS, j0, Len);

            Debug.Assert(output.Dimension == 4);
            Debug.Assert(output.GetLength(0) == Len);
            Debug.Assert(output.GetLength(1) == K);
            Debug.Assert(output.GetLength(2) == D);
            Debug.Assert(output.GetLength(3) == D);
            Debug.Assert(Jacobian.Dimension == 4);
            Debug.Assert(Jacobian.GetLength(0) == Len);
            Debug.Assert(Jacobian.GetLength(1) == K);
            Debug.Assert(Jacobian.GetLength(2) == D);
            Debug.Assert(Jacobian.GetLength(3) == D);

            output.Set(Jacobian);

            for(int i = 0; i < Len; i++) {
                Adjugate(output, i, D, K);
            }
        }

        private static void Adjugate(MultidimensionalArray AdjugateJacobianOut, int jOffset, int D, int L) {
            if(D == 1) {
                // 1D - adjungate
                for(int l = 0; l < L; l++)
                    AdjugateJacobianOut[jOffset, l, 0, 0] = 1.0;
            } else if(D == 2) {
                // 2D - adjungate

                for(int l = 0; l < L; l++) {
                    double m00 = AdjugateJacobianOut[jOffset, l, 0, 0];
                    double m01 = AdjugateJacobianOut[jOffset, l, 0, 1];
                    double m10 = AdjugateJacobianOut[jOffset, l, 1, 0];
                    double m11 = AdjugateJacobianOut[jOffset, l, 1, 1];


                    AdjugateJacobianOut[jOffset, l, 0, 0] = m11;
                    AdjugateJacobianOut[jOffset, l, 0, 1] = -m01;
                    AdjugateJacobianOut[jOffset, l, 1, 0] = -m10;
                    AdjugateJacobianOut[jOffset, l, 1, 1] = m00;
                }
            } else if(D == 3) {
                // 3D - adjungate

                for(int l = 0; l < L; l++) {
                    double m00 = AdjugateJacobianOut[jOffset, l, 0, 0];
                    double m01 = AdjugateJacobianOut[jOffset, l, 0, 1];
                    double m02 = AdjugateJacobianOut[jOffset, l, 0, 2];
                    double m10 = AdjugateJacobianOut[jOffset, l, 1, 0];
                    double m11 = AdjugateJacobianOut[jOffset, l, 1, 1];
                    double m12 = AdjugateJacobianOut[jOffset, l, 1, 2];
                    double m20 = AdjugateJacobianOut[jOffset, l, 2, 0];
                    double m21 = AdjugateJacobianOut[jOffset, l, 2, 1];
                    double m22 = AdjugateJacobianOut[jOffset, l, 2, 2];


                    AdjugateJacobianOut[jOffset, l, 0, 0] = m11 * m22 - m12 * m21;
                    AdjugateJacobianOut[jOffset, l, 0, 1] = -m01 * m22 + m02 * m21;
                    AdjugateJacobianOut[jOffset, l, 0, 2] = m01 * m12 - m02 * m11;
                    AdjugateJacobianOut[jOffset, l, 1, 0] = -m10 * m22 + m12 * m20;
                    AdjugateJacobianOut[jOffset, l, 1, 1] = m00 * m22 - m02 * m20;
                    AdjugateJacobianOut[jOffset, l, 1, 2] = -m00 * m12 + m02 * m10;
                    AdjugateJacobianOut[jOffset, l, 2, 0] = m10 * m21 - m11 * m20;
                    AdjugateJacobianOut[jOffset, l, 2, 1] = -m00 * m21 + m01 * m20;
                    AdjugateJacobianOut[jOffset, l, 2, 2] = m00 * m11 - m01 * m10;
                }
            } else {
                throw new NotSupportedException("spatial dimension of " + D + " is not supported yet.");
            }
        }

        /// <summary>
        /// Inverse of the Jacobian of the reference-to-physical coordinate transformation, \f$ ( \nabla T_j )^{-1} \f$.
        /// </summary>
        public Caching.CacheLogicImplBy_CNs InverseJacobian {
            get;
            private set;
        }

        void EvaluateInverseJacobian(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
            MultidimensionalArray AdjJac = this.AdjungateJacobian.GetValue_Cell(NS, j0, Len);
            MultidimensionalArray DetJac = this.JacobianDeterminat.GetValue_Cell(NS, j0, Len);
            int K = NS.NoOfNodes;
            int D = this.SpatialDimension;

            Debug.Assert(output.Dimension == 4);
            Debug.Assert(output.GetLength(0) == Len);
            Debug.Assert(output.GetLength(1) == K);
            Debug.Assert(output.GetLength(2) == D);
            Debug.Assert(output.GetLength(3) == D);
            Debug.Assert(AdjJac.Dimension == 4);
            Debug.Assert(AdjJac.GetLength(0) == Len);
            Debug.Assert(AdjJac.GetLength(1) == K);
            Debug.Assert(AdjJac.GetLength(2) == D);
            Debug.Assert(AdjJac.GetLength(3) == D);
            Debug.Assert(DetJac.Dimension == 2);
            Debug.Assert(DetJac.GetLength(0) == Len);
            Debug.Assert(DetJac.GetLength(1) == K);

            output.Set(AdjJac);

            for(int i = 0; i < Len; i++) {
                for(int k = 0; k < K; k++) {
                    output.ExtractSubArrayShallow(i, k, -1, -1).Scale(1.0 / DetJac[i, k]);
                }
            }
        }

        /// <summary>
        /// Determinant of the Jacobian of the reference-to-physical coordinate transformation, \f$ \determinant{ \nabla T_j } \f$.
        /// </summary>
        public Caching.CacheLogicImplBy_CNs JacobianDeterminat {
            get;
            private set;
        }

        void EvaluateJacobianDeterminat(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
            int K = NS.NoOfNodes;
            int D = this.SpatialDimension;

            MultidimensionalArray Jacobian = this.Jacobian.GetValue_Cell(NS, j0, Len);

            Debug.Assert(output.Dimension == 2);
            Debug.Assert(output.GetLength(0) == Len);
            Debug.Assert(output.GetLength(1) == K);
            Debug.Assert(Jacobian.Dimension == 4);
            Debug.Assert(Jacobian.GetLength(0) == Len);
            Debug.Assert(Jacobian.GetLength(1) == K);
            Debug.Assert(Jacobian.GetLength(2) == D);
            Debug.Assert(Jacobian.GetLength(3) == D);

            unsafe {
                if(!Jacobian.IsContinious)
                    throw new NotSupportedException();

                fixed(double *pJacobian = Jacobian.Storage) {

                    for(int i = 0; i < Len; i++) {
                        for(int k = 0; k < K; k++) {
                            int offset = Jacobian.Index(i, k, 0, 0);
                            output[i, k] = Determinant(D, pJacobian + offset);
                        }
                    }
                }
            }
        }

        unsafe private static double Determinant(int D, double* Jacobian) {
            double det;
            switch(D) {
                case 1:
                det = *Jacobian;
                break;

                case 2:
                det = (Jacobian[0 * 2 + 0] * Jacobian[1 * 2 + 1] - Jacobian[0 * 2 + 1] * Jacobian[1 * 2 + 0]);
                break;

                case 3:
                double m00 = Jacobian[0 * 3 + 0];
                double m01 = Jacobian[0 * 3 + 1];
                double m02 = Jacobian[0 * 3 + 2];
                double m10 = Jacobian[1 * 3 + 0];
                double m11 = Jacobian[1 * 3 + 1];
                double m12 = Jacobian[1 * 3 + 2];
                double m20 = Jacobian[2 * 3 + 0];
                double m21 = Jacobian[2 * 3 + 1];
                double m22 = Jacobian[2 * 3 + 2];
                det = m00 * m11 * m22 - m00 * m12 * m21 - m01 * m10 * m22 + m01 * m12 * m20 + m02 * m10 * m21 - m02 * m11 * m20;
                break;

                default:
                throw new NotSupportedException();
            }
            return det;
        }


        /// <summary>
        /// transforms the <paramref name="NS"/>-th node set to global coordinates of
        /// cell <paramref name="j0"/> to <paramref name="j0"/>+<paramref name="Len"/>-1
        /// </summary>
        /// <param name="j0">first cell to transform</param>
        /// <param name="Len">number of cells to transform</param>
        /// <param name="NS"></param>
        /// <param name="Nodesglob">
        /// output, global coordinates;
        /// <list type="bullet">
        ///   <item>
        ///     1st index: local cell index minus <paramref name="j0"/>, in the
        ///     range of 0 to <paramref name="Len"/>-1;
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of the
        ///     local nodes references by <paramref name="NS"/>
        ///   </item>
        ///   <item>
        ///     3rd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///     0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        public void TransformLocal2Global(MultidimensionalArray NS, int j0, int Len, MultidimensionalArray Nodesglob) {
            this.TransformLocal2Global(NS, j0, Len, Nodesglob, 0);
        }

        
        /// <summary>
        /// transforms vertices from the global coordinate system
        /// the local coordinate systems 
        /// of cells <paramref name="j0"/> to <paramref name="j0"/>+<paramref name="Len"/>-1.
        /// </summary>
        /// <param name="GlobalVerticesIn">
        /// Input; vertices in the global coordinate system;
        /// <list type="bullet">
        ///   <item>1st index: vertex index;</item>
        ///   <item>
        ///   2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///   0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="LocalVerticesOut">
        /// Output; the vertices form <paramref name="GlobalVerticesIn"/>,
        /// transformed to local coordinates;
        /// <list type="bullet">
        ///   <item> 
        ///     1st index: local cell index minus <paramref name="j0"/>, in the
        ///     range of 0 to <paramref name="Len"/>-1;
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="GlobalVerticesIn"/>;
        ///   </item>
        ///   <item>
        ///     3rd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///     0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="j0">local cell index of the first cell to transform</param>
        /// <param name="Len">Number of cells to transform</param>
        /// <param name="OutArrayOffset">
        /// an offset into the first index of <paramref name="LocalVerticesOut"/>;
        /// </param>
        public void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int j0, int Len, int OutArrayOffset) {
            int N = GlobalVerticesIn.GetLength(0);
            int D = SpatialDimension;
            if (GlobalVerticesIn.Dimension != 2)
                throw new ArgumentException("wrong dimension", "GlobalVerticesIn");
            if (LocalVerticesOut.Dimension != 3)
                throw new ArgumentException("wrong dimension", "LocalVerticesOut");
            if (GlobalVerticesIn.GetLength(1) != D)
                throw new ArgumentException("wrong spatial dimension of GlobalVerticesIn", "GlobalVerticesIn");
            if (LocalVerticesOut.GetLength(2) != D)
                throw new ArgumentException("wrong spatial dimension of LocalVerticesOut", "LocalVerticesOut");
            if (GlobalVerticesIn.GetLength(0) != LocalVerticesOut.GetLength(1))
                throw new ArgumentException("mismatch in number of vertices per cell.", "GlobalVerticesIn,LocalVerticesOut");
            if (LocalVerticesOut.GetLength(0) < Len + OutArrayOffset)
                throw new ArgumentException("insufficient number of cells in output array.", "OutArrayOffset");

            for (int j = 0; j < Len; j++) {
                var _LocalVerticesOut = LocalVerticesOut.ExtractSubArrayShallow(j + OutArrayOffset, -1, -1);
                TransformGlobal2Local(GlobalVerticesIn, _LocalVerticesOut, j0 + j, null);
            }
        }

#if NEWTON_DIAGNOSIS
        bool printed_NEWTON_DIAGNOSIS_reminder = false;
#endif

        /// <summary>
        /// transforms vertices from the global coordinate system
        /// the local coordinate systems 
        /// of cells
        /// <paramref name="jCell"/>.
        /// </summary>
        /// <param name="GlobalVerticesIn">
        /// Input; vertices in the global coordinate system;
        /// <list type="bullet">
        ///   <item>1st index: vertex index;</item>
        ///   <item>
        ///   2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///   0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="LocalVerticesOut">
        /// Output; the vertices form <paramref name="GlobalVerticesIn"/>,
        /// transformed to local coordinates;
        /// <list type="bullet">
        ///   <item>
        ///     1st index: vertex index, corresponds with the 1st index of
        ///     <paramref name="GlobalVerticesIn"/>;
        ///   </item>
        ///   <item>
        ///     2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///     0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="jCell">local cell index of the cell to transform</param>
        public void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int jCell, bool[] NewtonConvergence) {
            int N = GlobalVerticesIn.GetLength(0); // number of nodes/points
            int D = SpatialDimension;
            if (GlobalVerticesIn.Dimension != 2)
                throw new ArgumentException("wrong dimension", "GlobalVerticesIn");
            if (LocalVerticesOut.Dimension != 2)
                throw new ArgumentException("wrong dimension", "LocalVerticesOut");
            if (GlobalVerticesIn.GetLength(1) != D)
                throw new ArgumentException("wrong spatial dimension of GlobalVerticesIn", "GlobalVerticesIn");
            if (LocalVerticesOut.GetLength(1) != D)
                throw new ArgumentException("wrong spatial dimension of LocalVerticesOut", "LocalVerticesOut");
            if (GlobalVerticesIn.GetLength(0) != LocalVerticesOut.GetLength(0))
                throw new ArgumentException("mismatch in number of vertices per cell.", "GlobalVerticesIn,LocalVerticesOut");

                        
            var Kref = this.Cells.GetRefElement(jCell);

            if (NewtonConvergence != null)
                ArrayTools.SetAll(NewtonConvergence, false);

#if NEWTON_DIAGNOSIS
            if (!printed_NEWTON_DIAGNOSIS_reminder) {
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("REMINDER: Newton algorithm diagnostic code is active;");
                Console.WriteLine("          serious performance impact expected;");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                printed_NEWTON_DIAGNOSIS_reminder = true;
            }

            var Points = new List<MultidimensionalArray>();
            var ResidualHistory = new List<double>();
#endif


            Cell Kj = this.Cells.GetCell(jCell);
            if (!Kj.Type.IsLinear()) {
                // ++++++++++++++
                // nonlinear cell
                // ++++++++++++++

                MultidimensionalArray RES = MultidimensionalArray.Create(1, N, D);
                MultidimensionalArray _RES = RES.ResizeShallow(N, D);

                double h = Kref.GetMaxDiameter();

                var _LocalVerticesOut = LocalVerticesOut;

                // Newton iteration start
                _LocalVerticesOut.Clear();

                // Newton iterations
                double[] OldResidualNorm = new double[N]; OldResidualNorm.SetAll(double.MaxValue); 
                double[] MimimumResidual_SoFar = new double[N]; MimimumResidual_SoFar.SetAll(double.MaxValue);
                double[] ResidualNorm = new double[N];
                int itercnt = 0;
                int[] DivergenceDetector = new int[N];
                int[] convergenceDetector = new int[N];

                while (true) {

                    // RES = GlobalVerticesIn - Trfao(x[i-1])
                    RES.Clear();
                    NodeSet __LocalVerticesOut = new NodeSet(Kref, _LocalVerticesOut);
                    this.TransformLocal2Global(__LocalVerticesOut, jCell, 1, RES, 0);
                    RES.Scale(-1.0);
                    _RES.Acc(1.0, GlobalVerticesIn);

#if DEBUG
                    if(_LocalVerticesOut.CheckForNanOrInf(true,true,false) > 0) {
                        throw new ArithmeticException("Newton algorithm divergence.");
                    }
#endif
                    int undecided = N;
                    for (int n = 0; n < N; n++) {
                        bool decicion = false;
                        ResidualNorm[n] = _RES.GetRow(n).L2Norm();
                        if (ResidualNorm[n] > OldResidualNorm[n]) {
                            DivergenceDetector[n]++;
                        } else {
                            DivergenceDetector[n] = 0;
                        }

                        if (ResidualNorm[n] >= MimimumResidual_SoFar[n])
                            convergenceDetector[n]++;
                        else
                            convergenceDetector[n] = 0;
                        MimimumResidual_SoFar[n] = Math.Min(MimimumResidual_SoFar[n], ResidualNorm[n]);

                        if (DivergenceDetector[n] > 10)
                            decicion = true;

                        if ((ResidualNorm[n] <= 1.0e-10 * h || convergenceDetector[n] > 8)) {
                            if (NewtonConvergence != null)
                                NewtonConvergence[n] = true;
                            decicion = true;
                        }

                        if (decicion)
                            undecided--;
                    }
                    

                    // inverse transformation matrix
                    //Kref.InverseJacobianOfTransformation(_LocalVerticesOut, __InverseMtx, 0, Kj.Type, Kj.TransformationParams);
                    MultidimensionalArray __InverseMtx = this.InverseJacobian.GetValue_Cell(__LocalVerticesOut, jCell, 1);



                    var __Inv = __InverseMtx.ExtractSubArrayShallow(0, -1, -1, -1);

                    // x[i] = x[i-1] + Inv*RES
                    _LocalVerticesOut.Multiply(1, __Inv, _RES, 1.0, "ki", "kij", "kj");
#if NEWTON_DIAGNOSIS
                    Points.Add(_LocalVerticesOut.CloneAs());
                    ResidualHistory.Add(ResidualNorm);
#endif

                    itercnt++;
                    if (itercnt > 100 || undecided == 0) {
#if NEWTON_DIAGNOSIS
                        FullMatrix Output = new FullMatrix(Points.Count, D*N);
                        for (int ii = 0; ii < Points.Count; ii++) {
                            var _Points = Points[ii];

                            for (int nn = 0; nn < N; nn++) {
                                for (int dd = 0; dd < D; dd++) {
                                    Output[ii, nn*D + dd] = _Points[nn, dd];
                                }
                            }
                        }
                        Output.ToTxtFile("C:\\tmp\\Newton-Problem.txt");
#endif
                        //throw new ArithmeticException("Newton algorithm divergence.");
                        return;
                    }
                }

            } else {
                // +++++++++++
                // linear cell
                // +++++++++++

                // linear cells always converge in one Newton step:
                if(NewtonConvergence != null)
                    NewtonConvergence.SetAll(true);
                                
                // affine Offset
                MultidimensionalArray offset = this.GlobalNodes.GetValue_Cell(Kref.Center, jCell, 1);

                // inverse transformation matrix
                MultidimensionalArray InverseMtx = this.InverseJacobian.GetValue_Cell(Kref.Center, jCell, 1);

                double[] vec = new double[this.SpatialDimension];

                // transform: loop over vertices...
                for (int n = 0; n < N; n++) {
                    for (int d = 0; d < D; d++) {
                        vec[d] = GlobalVerticesIn[n, d] - offset[0, 0, d];
                    }

                    for (int d = 0; d < D; d++) {
                        double s = 0;
                        for (int dd = 0; dd < D; dd++) {
                            s += InverseMtx[0,0, d, dd] * vec[dd];
                        }
                        LocalVerticesOut[n, d] = s;
                    }
                }
            }
        }

 
        /// <summary>
        /// Parallel setup, e.g. send and receive lists.
        /// - 'ínput data': cell-neighborship information obtained from <see cref="GridCommons.GetCellNeighbourship"/>.
        /// - 'output data': <see cref="CellData.CellNeighbours_global_tmp"/>, <see cref="CellData.CellNeighbours"/> and various entries of <see cref="Parallel"/>
        /// </summary>
        void ParallelSetup() {
            using (new FuncTrace()) {
                // global indices for all cells
                // ============================

                m_Grid.InitNumberOfCells();
                int Jglob = m_Grid.NumberOfCells;
                int J = m_Grid.NoOfUpdateCells;
                int J_BC = m_Grid.NoOfBcCells;
                var Part = m_Grid.CellPartitioning;
                var BcCellPart = m_Grid.BcCellPartitioning;
                int MyRank = this.MpiRank;

                // compute neighborship info 
                // =========================

                var CNglb = m_Grid.GetCellNeighbourship(true);
                Debug.Assert(CNglb.Length == (J + J_BC));
#if DEBUG
                for(int j = 0; j < J; j++) {
                    var CNglb_j = CNglb[j].ToArray();

                    for(int n1 = 0; n1 < CNglb_j.Length; n1++) {
                        for (int n2 = 0; n2 < CNglb_j.Length; n2++) {
                            if(n1 != n2) {
                                if((CNglb_j[n1].Neighbour_GlobalIndex == CNglb_j[n2].Neighbour_GlobalIndex) && (CNglb_j[n1].Neighbour_GlobalIndex >= 0 )  && (CNglb_j[n2].Neighbour_GlobalIndex >= 0 )) {
                                    long GlId0 = m_Grid.Cells[j].GlobalID;
                                    
                                    throw new ApplicationException("Fatal error in cell graph: edge between " + GlId0
                                            + " and some other cell is defined multiple times.");
                                }
                            }
                        }
                    }
                }
#endif


                // separate normal cells and boundary-condition -- cells
                // =====================================================

                m_Cells.CellNeighbours_global_tmp = new IEnumerable<GridCommons.Neighbour>[J];
                Array.Copy(CNglb, m_Cells.CellNeighbours_global_tmp, J);

                var BcCNglb = new IEnumerable<GridCommons.Neighbour>[J_BC];
                Array.Copy(CNglb, J, BcCNglb, 0, J_BC);

                var NeighGlobalIdx = m_Cells.CellNeighbours_global_tmp;

                // define External/ghost cells, sort them according to MPI process rank
                // ====================================================================

                // keys: MPI process rank 'p'
                // values: collection of ghost cells which belong to processor 'p'
                Dictionary<int, List<long>> ExternalCells = new Dictionary<int, List<long>>();

                for (int j = 0; j < J; j++) { // loop over cells
                    int Nj = NeighGlobalIdx[j].Count();

                    for (int n = 0; n < Nj; n++) { // loop over faces
                        var idx = NeighGlobalIdx[j].ElementAt(n).Neighbour_GlobalIndex;
                        Debug.Assert(idx >= 0);

                        if (idx >= Jglob) {
                            // link to a boundary cell - not relevant yet
                            continue;
                        }

                        if (!Part.IsInLocalRange(idx)) {
                            // found external cell

                            int rnk = Part.FindProcess(idx);

                            List<long> ext_rnk;
                            ExternalCells.TryGetValue(rnk, out ext_rnk);
                            if (ext_rnk == null) {
                                ext_rnk = new List<long>();
                                ExternalCells.Add(rnk, ext_rnk);
                            }

                            ext_rnk.Add(idx);
                        } else {
                            // nop
                        }
                    }
                }

                // sort list of external cells, remove duplicates
                int Jexternal = 0;
                foreach (var p in ExternalCells.Keys.ToArray()) {
                    var list = ExternalCells[p];

                    list.Sort();
                    List<long> list2 = new List<long>(list.Count);
                    long old_l = long.MinValue;
                    for (int i = 0; i < list.Count; i++) {
                        var l = list[i];
                        if (l != old_l)
                            list2.Add(l);
                        old_l = l;
                    }

                    ExternalCells[p] = list2;

                    Jexternal += list2.Count;
                }

                // record external cells, build receive lists
                // ==========================================
                Dictionary<long, int> Global2LocalIdx = new Dictionary<long, int>();
                m_Parallel.Global2LocalIdx = Global2LocalIdx;
                {
                    m_Parallel.GlobalIndicesExternalCells = new long[Jexternal];

                    m_Parallel.RcvCommListsInsertIndex = new int[MpiSize];
                    m_Parallel.RcvCommListsNoOfItems = new int[MpiSize];

                    int cnt = 0;
                    foreach (var kv in ExternalCells) {
                        int proc = kv.Key;
                        var list = kv.Value;

                        int L = list.Count;
                        m_Parallel.RcvCommListsNoOfItems[proc] = L;
                        m_Parallel.RcvCommListsInsertIndex[proc] = cnt + J;

                        for (int l = 0; l < L; l++)
                            Global2LocalIdx.Add(list[l], cnt + J + l); // here, the local indices for the external/ghost cells are determined;

                        list.CopyTo(m_Parallel.GlobalIndicesExternalCells, cnt);
                        cnt += L;
                    }

                    m_Parallel.ProcessesToReceiveFrom = ExternalCells.Keys.ToArray();
                }


                m_Cells.m_CellNeighbours = new int[J][];
                var ClNg = m_Cells.m_CellNeighbours;
                for (int j = 0; j < J; j++) {
                    int Nj = NeighGlobalIdx[j].Count();
                    ClNg[j] = new int[Nj];

                    int cnt = 0;
                    for (int n = 0; n < Nj; n++) {
                        var cn = NeighGlobalIdx[j].ElementAt(n);
                        long idx = cn.Neighbour_GlobalIndex;


                        if (idx >= 0 && idx < Jglob) {
                            if (!Part.IsInLocalRange(idx)) {
                                ClNg[j][cnt] = Global2LocalIdx[idx];
                            } else {
                                ClNg[j][cnt] = Part.TransformIndexToLocal((int)idx);
                            }
                            cnt++;
                        } else {
                            //ClNg[j][cnt] = int.MinValue;
                            //cnt++;
                        }
                    }

                    if(cnt < Nj) {
                        Array.Resize(ref ClNg[j], cnt);
                        Nj = cnt;
                    }


#if DEBUG
                    for(int n1 = 0; n1 < Nj; n1++) {
                        for(int n2 = 0; n2 < Nj; n2++) {
                            if(n1 != n2) {
                                if(ClNg[j][n1] >= 0)
                                    Debug.Assert(ClNg[j][n1] != ClNg[j][n2], "Same neighbor cell specified twice");
                            }
                        }
                    }
                    Debug.Assert(ClNg[j].Contains(j) == false, "cell seems to be its own neighbor.");
                    //Debug.Assert(ClNg[j].Where(i => i < 0).Count() == 0);

#endif
                }

                // build send lists
                // ================
                {
                    SerialisationMessenger sms = new SerialisationMessenger(MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    sms.SetCommPathsAndCommit(ExternalCells.Keys);

                    foreach (var msg in ExternalCells)
                        sms.Transmitt(msg.Key, msg.Value.ToArray());

                    m_Parallel.SendCommLists = new int[MpiSize][];

                    List<int> procsToSendTo = new List<int>();
                    int rcvRank;
                    long[] rcv;
                    while (sms.GetNext(out rcvRank, out rcv)) {
                        procsToSendTo.Add(rcvRank);


                        int I = rcv.Length;
                        int[] SendList = new int[I];
                        m_Parallel.SendCommLists[rcvRank] = SendList;

                        for (int i = 0; i < I; i++) {
                            Debug.Assert(Part.IsInLocalRange(rcv[i]));
                            SendList[i] = Part.TransformIndexToLocal((int)rcv[i]);
                        }
                    }

                    m_Parallel.ProcessesToSendTo = procsToSendTo.ToArray();

                    sms.Dispose();
                }

                // external cell data
                // ==================
                {
                    SerialisationMessenger sms = new SerialisationMessenger(MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    sms.SetCommPathsAndCommit(m_Parallel.ProcessesToSendTo);

                    foreach (var p in m_Parallel.ProcessesToSendTo) {
                        var sendList = m_Parallel.SendCommLists[p];
                        var L = sendList.Length;

                        Cell[] msg = new Cell[L];
                        for (int l = 0; l < L; l++)
                            msg[l] = m_Grid.Cells[sendList[l]];

                        sms.Transmitt(p, msg);
                    }

                    m_Parallel.ExternalCells = new Cell[Jexternal];

                    Cell[] rcvMsg;
                    int rcvRnk;
                    while (sms.GetNext(out rcvRnk, out rcvMsg)) {
                        Array.Copy(rcvMsg, 0, m_Parallel.ExternalCells, m_Parallel.RcvCommListsInsertIndex[rcvRnk] - J, rcvMsg.Length);
                    }
                }

                // Neighbour data
                // ==============
                {
                    Array.Resize(ref m_Cells.CellNeighbours_global_tmp, m_Cells.CellNeighbours_global_tmp.Length + Jexternal);
                    NeighGlobalIdx = m_Cells.CellNeighbours_global_tmp;

                    Dictionary<int, Helper> send_data = new Dictionary<int, Helper>();
                    foreach (int targ_proc in this.Parallel.ProcessesToSendTo) {
                        int[] SendList = this.Parallel.SendCommLists[targ_proc];

                        var buf = new GridCommons.Neighbour[SendList.Length][];
                        for (int i = 0; i < SendList.Length; i++) {
                            buf[i] = NeighGlobalIdx[SendList[i]].ToArray();
                        }

                        send_data.Add(targ_proc, new Helper() {
                            entries = buf
                        });
                    }


                    var rcv_data = SerialisationMessenger.ExchangeData(send_data, csMPI.Raw._COMM.WORLD);

                    foreach (var kv in rcv_data) {
                        int rcv_rank = kv.Key;
                        GridCommons.Neighbour[][] data = kv.Value.entries;

                        int L = data.Length;
                        Debug.Assert(L == this.Parallel.RcvCommListsNoOfItems[rcv_rank]);
                        int offset = this.Parallel.RcvCommListsInsertIndex[rcv_rank];
                        for (int l = 0; l < L; l++) {
                            NeighGlobalIdx[l + offset] = data[l];
                        }
                    }
                }

                // Exchange boundary-condition cells 
                // =================================
                {
                    Dictionary<int, BCElement> BcCells = new Dictionary<int, BCElement>();

                    Dictionary<int, Dictionary<int, BCElement>> SendBcCells = new Dictionary<int, Dictionary<int, BCElement>>();

                    int j_bc = 0;
                    int jCellGlob = Jglob + BcCellPart.i0;
                    foreach (var icn in BcCNglb) { // loop over all BcCell's on this processor
                        foreach (var cn in icn) {
                            long jVolCell = cn.Neighbour_GlobalIndex;

                            int R = Part.FindProcess((int)jVolCell); // the j_bc -- th boundary condition cell is required @ MPI process R

                            if (R == MyRank) {
                                BcCells.Add(jCellGlob, m_Grid.BcCells[j_bc]);
                            } else {
                                Dictionary<int, BCElement> _BcCells;
                                if (!SendBcCells.TryGetValue(R, out _BcCells)) {
                                    _BcCells = new Dictionary<int, BCElement>();
                                    SendBcCells.Add(R, _BcCells);
                                }

                                _BcCells.Add(jCellGlob, m_Grid.BcCells[j_bc]);
                            }
                        }

                        jCellGlob++;
                        j_bc++;
                    }

                    var rcvData = SerialisationMessenger.ExchangeData(SendBcCells, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    foreach (var T in rcvData.Values) {
                        foreach (var kv in T)
                            BcCells.Add(kv.Key, kv.Value);
                    }

                    m_BcCells_tmp = BcCells;
                }
            }
        }

        /// <summary>
        /// Temporary boundary condition cells during setup; 
        /// This dictionary (should) contain all boundary-condition cells
        /// which are required on the current MPI process.
        /// - key: global cell index of the boundary-condition cell. 
        /// - value: the boundary-condition cell.
        /// </summary>
        internal Dictionary<int, BCElement> m_BcCells_tmp;

        /// workaround for some mono bug in BinaryFormatter:
        /// It seems that the BinaryFormatter in mono has some problems
        /// (de)serializing staggered arrays directly. However, if the
        /// staggered array is packed into a class it seems to work.
        [Serializable]
        class Helper {
            public GridCommons.Neighbour[][] entries;
        }
                

        SubGrid[] m_RefElementSgrd;

        /// <summary>
        /// Subgrid of all cells which share the same reference element.
        /// </summary>
        /// <param name="iKref">
        /// index into <see cref="GridCommons.RefElements"/>.
        /// </param>
        public SubGrid GetRefElementSubGrid(int iKref) {
            MPICollectiveWatchDog.Watch();

            var KRefs = this.Grid.RefElements;
            if (m_RefElementSgrd == null) {
                m_RefElementSgrd = new SubGrid[KRefs.Length];
            }

            if (iKref < 0 || iKref >= KRefs.Length) {
                throw new ArgumentOutOfRangeException();
            }

            if (m_RefElementSgrd[iKref] == null) {
                m_RefElementSgrd[iKref] = new SubGrid(this.Cells.GetCells4Refelement(iKref).ToLogicalMask());
            }

            return m_RefElementSgrd[iKref];
        }

        SubGrid m_AffineLinearCells;

        /// <summary>
        /// subgrid that contains only cells 
        /// for which an affine-linear transformation to the reference element 
        /// can be defined.
        /// </summary>
        public SubGrid AffineLinearCells {
            get {
                MPICollectiveWatchDog.Watch();
                if (m_AffineLinearCells == null) {

                    if (m_Cells.ContainsNonlinearCell()) {
                        int J = this.Cells.NoOfLocalUpdatedCells;
                        BitArray msk = new BitArray(J);
                        for (int j = 0; j < J; j++)
                            msk[j] = (this.Cells.IsCellAffineLinear(j));

                        CellMask cm = new CellMask(this, msk);
                        m_AffineLinearCells = new SubGrid(cm);
                    } else {
                        m_AffineLinearCells = new SubGrid(CellMask.GetFullMask(this));
                    }
                }


                return m_AffineLinearCells;
            }
        }

        SubGrid m_NonlinearCells;

        /// <summary>
        /// the complement of <see cref="AffineLinearCells"/>.
        /// </summary>
        public SubGrid NonlinearCells {
            get {
                MPICollectiveWatchDog.Watch();
                if (m_NonlinearCells == null) {

                    if (m_Cells.ContainsNonlinearCell())
                        m_NonlinearCells = AffineLinearCells.Complement();
                    else
                        m_NonlinearCells = new SubGrid(CellMask.GetEmptyMask(this));
                }
                return m_NonlinearCells;
            }
        }

        EdgeMask m_BoundaryEdges;

        /// <summary>
        /// a quadrature execution mask, for edges, which contains all
        /// boundary edges, i.e. edges that bound only to one cell.
        /// </summary>
        public EdgeMask BoundaryEdges {
            get {
                if (m_BoundaryEdges == null)
                    DefineBoundaryMasks();

                return m_BoundaryEdges;
            }
        }

        SubGrid m_BoundaryCells;

        /// <summary>
        /// subgrid that contains only cells which lie on the boundary
        /// </summary>
        public SubGrid BoundaryCells {
            get {
                MPICollectiveWatchDog.Watch();
                if (m_BoundaryCells == null)
                    DefineBoundaryMasks();

                return m_BoundaryCells;
            }
        }

        SubGrid m_InnerCells;

        /// <summary>
        /// the complementary subgrid to <see cref="BoundaryCells"/>
        /// </summary>
        public SubGrid InnerCells {
            get {
                MPICollectiveWatchDog.Watch();
                if (m_InnerCells == null) {
                    m_InnerCells = new SubGrid(m_BoundaryCells.VolumeMask.Complement<CellMask>());
                }
                return m_InnerCells;
            }
        }

        Dictionary<RefElement, SubGrid> m_Subgrid4RefElement;


        /// <summary>
        /// for each cell type, the associated subgrid.
        /// </summary>
        /// <remarks>
        /// In contrast to <see cref="SubGrid"/>-objects, <see cref="CellMask"/>-objects have almost no 
        /// MPI-collective operations (with some exceptions).
        /// </remarks>
        public IDictionary<RefElement, SubGrid> Subgrid4RefElement {
            get {
                MPICollectiveWatchDog.Watch();
                if (m_Subgrid4RefElement == null) {
                    m_Subgrid4RefElement = new Dictionary<RefElement, SubGrid>(((Func<RefElement, RefElement, bool>)((a, b) => object.ReferenceEquals(a, b))).ToEqualityComparer());

                    var Krefs = this.Grid.RefElements;
                    for (int iKref = 0; iKref < Krefs.Length; iKref++) {
                        m_Subgrid4RefElement.Add(Krefs[iKref], this.GetRefElementSubGrid(iKref));
                    }
                }

                return m_Subgrid4RefElement;
            }
        }

        /*

        Dictionary<RefElement, CellMask> m_Cellmask4RefElement;


        /// <summary>
        /// for each cell type, the associated cell-mask.
        /// </summary>
        /// <remarks>
        /// In contrast to <see cref="SubGrid"/>-objects, <see cref="CellMask"/>-objects have almost no 
        /// MPI-collective operations (with some exceptions).
        /// </remarks>
        public IDictionary<RefElement, CellMask> Cellmask4RefElement {
            get {
                if (m_Cellmask4RefElement == null) {
                    m_Cellmask4RefElement = new Dictionary<RefElement, CellMask>(((Func<RefElement, RefElement, bool>)((a, b) => object.ReferenceEquals(a, b))).ToEqualityComparer());

                    var Krefs = this.Grid.RefElements;
                    for (int iKref = 0; iKref < Krefs.Length; iKref++) {
                        m_Cellmask4RefElement.Add(Krefs[iKref], this.GetRefElementCellMask(iKref));
                    }
                }

                return m_Cellmask4RefElement;
            }
        }
        */

        /// <summary>
        /// initializes <see cref="m_BoundaryCells"/>, <see cref="m_BoundaryEdges"/>;
        /// </summary>
        void DefineBoundaryMasks() {
            int J = m_Cells.NoOfLocalUpdatedCells;
            BitArray boundaryCells = new BitArray(J);

            int E = m_Edges.Count;
            BitArray boundaryEdges = new BitArray(E);
            BitArray boundaryCellsEdges = new BitArray(E);


            // loop over all Edges
            for (int e = 0; e < E; e++) {
                int Cel1 = m_Edges.CellIndices[e, 0];
                int Cel2 = m_Edges.CellIndices[e, 1];

                if (Cel2 < 0) {
                    // edge is located on the computational domain boundary
                    boundaryEdges[e] = true;

                    boundaryCells[Cel1] = true;
                }


                if (boundaryCells[Cel1] == true)
                    // edge belongs to a cell on the computational boundary
                    boundaryCellsEdges[e] = true;

                if (Cel2 >= 0 && Cel2 < J && boundaryCells[Cel2] == true)
                    // edge belongs to a cell on the computational boundary
                    boundaryCellsEdges[e] = true;
            }


            m_BoundaryEdges = new EdgeMask(this, boundaryEdges);
            m_BoundaryCells = new SubGrid(new CellMask(this, boundaryCells));
        }

        

        /// <summary>
        /// Gets the partitioning of cells over the MPI processes;
        /// </summary>
        public Partitioning CellPartitioning {
            get {
                return Grid.CellPartitioning;
            }
        }

      

        
        /// <summary>
        /// 
        /// </summary>
        private void Initialize_LocalCellIndexToEdges() {
            using (new FuncTrace()) {
                throw new NotImplementedException("todo: fk.");
                /*
                int EDG = Edges.GetLength(0);
                int J = NoOfLocalUpdatedCells;
                int E = m_Context.Grid.GridSimplex.NoOfEdges;

                LocalCellIndexToEdges = new int[J, E];
                ArrayTools.Set(LocalCellIndexToEdges, int.MaxValue);

                int setDings = 0;
                for (int edg = 0; edg < EDG; edg++) {

                    int sign = 1;
                    for (int k = 0; k < 2; k++) {
                        int j_cell = Edges[edg, k];
                        if (j_cell < 0)
                            continue;
                        if (j_cell >= J)
                            continue;

                        int e_cell = EdgeIndices[edg, k];

                        if (LocalCellIndexToEdges[j_cell, e_cell] != int.MaxValue)
                            // already set
                            throw new ApplicationException("algorithm error;");

                        LocalCellIndexToEdges[j_cell, e_cell] = sign * (edg + 1);

                        sign *= -1;
                        setDings++;
                    }
                }

                if (setDings != J * E)
                    // some entry remained uninitialized;
                    throw new ApplicationException("error in algorithm.");
                 */
            }
        }

        /// <summary>
        /// The spatial dimension of the grid (usually 1, 2 or 3).
        /// </summary>
        public int SpatialDimension {
            get {
                return m_Grid.SpatialDimension;
            }
        }

        private Permutation m_CurrentGlobalIdPermutation;

        /// <summary>
        /// the current GlobalID-permutation;
        /// Equal to <see cref="GridCommons.GetGlobalIDPermutation"/>;
        /// </summary>
        public Permutation CurrentGlobalIdPermutation {
            get {
                if (this.m_CurrentGlobalIdPermutation == null) {
                    this.m_CurrentGlobalIdPermutation = this.Grid.GetGlobalIDPermutation(false);
                }
                return m_CurrentGlobalIdPermutation;
            }
        }

        /// <summary>
        /// checks edge tag names (<see cref="GridCommons.EdgeTagNames"/>) for uniqueness
        /// </summary>
        private void CheckEdgeTagNames() {
            using (new FuncTrace()) {

                List<string> testColl = new List<string>();

                foreach (string EdgeTagName in m_Grid.EdgeTagNames.Values) {
                    if (testColl.Contains(EdgeTagName))
                        throw new ApplicationException("EdgeTagName '" + EdgeTagName + "' appears more than once.");
                    testColl.Add(EdgeTagName);
                }
            }
        }


        /// <summary>
        /// The bounding box of this part of the grid, which is stored on the local MPI process.
        /// </summary>
        public BoundingBox LocalBoundingBox {
            get {
                if (m_LocalBoundingBox == null) {
                    m_LocalBoundingBox = new BoundingBox(this.SpatialDimension);
                    int J = this.Cells.NoOfLocalUpdatedCells;
                    var CellBB = new BoundingBox(this.SpatialDimension);

                    for (int j = 0; j < J; j++) {
                        this.Cells.GetCellBoundingBox(j, CellBB);
                        m_LocalBoundingBox.AddBB(CellBB);
                    }

                }
                return m_LocalBoundingBox.CloneAs();
            }
        }

        BoundingBox m_LocalBoundingBox;


        /// <summary>
        /// the bounding box of the entire grid.
        /// </summary>
        public BoundingBox GlobalBoundingBox {
            get {
                if (m_GlobalBoundingBox == null) {
                    var _LocalBoundingBox = this.LocalBoundingBox;
                    int D = this.SpatialDimension;
                    double[] MinMax = new double[2 * D];
                    for(int d = 0; d < D; d++) {
                        MinMax[d] = _LocalBoundingBox.Min[d];
                        MinMax[d + D] = -_LocalBoundingBox.Max[d];
                    }
                    MinMax = MinMax.MPIMin(); // save collective MPI calls by min/max trick

                    for (int d = 0; d < D; d++) {
                        _LocalBoundingBox.Min[d] = MinMax[d];
                        _LocalBoundingBox.Max[d] = -MinMax[d + D];
                    }
                    m_GlobalBoundingBox = _LocalBoundingBox;
                }
                return m_GlobalBoundingBox.CloneAs();
            }
        }

        BoundingBox m_GlobalBoundingBox;

    }
}
