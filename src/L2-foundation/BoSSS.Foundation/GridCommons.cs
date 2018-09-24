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
using System.Reflection;
using System.Runtime.Serialization;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// Defines basic grid information which are used to create a
    /// <see cref="GridData"/> object.
    /// </summary>
    [Serializable]
    [DataContract]
    public partial class GridCommons : IGridInfo, ICloneable, IGrid {


        static ILog Logger = LogManager.GetLogger(typeof(GridCommons));

        /// <summary>
        /// Cells of the grid.
        /// </summary>
        /// <remarks>
        /// Cannot be implemented as property since NonSerialized only works
        /// for fields.
        /// </remarks>
        [NonSerialized]
        public Cell[] Cells;

        /// <summary>
        /// Optional elements that mark boundary conditions. Their global indices resp. global Id's 
        /// start after those of the <see cref="Cells"/>.
        /// </summary>
        [NonSerialized]
        public BCElement[] BcCells;

        /// <summary>
        /// creates a new grid object; no data of this object is set;
        /// A new Guid (<see cref="GridGuid"/>) is created;
        /// </summary>
        public GridCommons(RefElement[] RefElm, RefElement[] EdgeRefElm) {
            using (new FuncTrace() ) { 
                ilPSP.MPICollectiveWatchDog.Watch();
                m_GridGuid = Guid.NewGuid();
                m_GridGuid = m_GridGuid.MPIBroadcast(0);
                EdgeTagNames.Add(0, "inner edge");
                this.m_CreationTime = DateTime.Now;

                int D = RefElm.First().SpatialDimension;
                foreach (var Kref in RefElm) {
                    if (Kref.SpatialDimension != D)
                        throw new ArgumentException("All reference elements must have the same spatial dimension.");
                }

                for (int i = 0; i < RefElm.Length; i++) {
                    for (int j = i + 1; j < RefElm.Length; j++) {
                        if (RefElm[i].GetType() == RefElm[j].GetType()) {
                            throw new ArgumentException("Each cell reference element may occur only once: fount two '" + RefElm[i].GetType().Name + "'.", "RefElm");
                        }
                    }
                }
                for (int i = 0; i < EdgeRefElm.Length; i++) {
                    for (int j = i + 1; j < EdgeRefElm.Length; j++) {
                        if (EdgeRefElm[i].GetType() == EdgeRefElm[j].GetType()) {
                            throw new ArgumentException("Each edge reference element may occur only once: fount two '" + RefElm[i].GetType().Name + "'.");
                        }
                    }
                }

                if (D < 1)
                    throw new ArgumentException("Spatial dimension of volume cells must be 1 at least.");

                int De = EdgeRefElm.First().SpatialDimension;
                foreach (var KrefEdge in EdgeRefElm) {
                    if (KrefEdge.SpatialDimension != De)
                        throw new ArgumentException("All edge reference elements must have the same spatial dimension.");
                }

                if (D != (De + 1))
                    throw new ArgumentException("Mismatch between spatial dimension of edge and volume reference elements.");

                this.m_RefElements = RefElm.CloneAs();
                this.m_EdgeRefElements = EdgeRefElm.CloneAs();
            }
        }


        private void HackCheck() {

#if DEBUG
            // test that each the face reference element, 
            //     for each cell/volume reference element, 
            //         is contained in the list of edge reference elements, 
            // by reference-equality!

            RefElement[] Edge_KrefS = this.EdgeRefElements;
            foreach (var Cell_Kref in this.RefElements) {
                RefElement CellFaceKref = Cell_Kref.FaceRefElement;

                int FoundCount = 0;
                for (int i = 0; i < Edge_KrefS.Count(); i++) {
                    if (object.ReferenceEquals(CellFaceKref, Edge_KrefS.ElementAt(i)))
                        FoundCount++;
                }
                Debug.Assert(FoundCount == 1);

                //Debug.Assert(Edge_KrefS.Where(A => object.ReferenceEquals(Cell_Kref.FaceRefElement, A)).Count() == 1);
            }

#endif
        }


        private GridCommons() {
        }

        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        /// <param name="EdgeTagFunc"></param>
        public void DefineEdgeTags(Func<double[], byte> EdgeTagFunc) {
            
            int Jloc = this.Cells.Length;
            int minJloc = Jloc.MPIMin();
            if (minJloc <= 0) {
                // redist is necessary
                this.Redistribute(null, GridPartType.METIS, null);
            }

            
            var GrdDatTmp = new GridData(this);

            int D = SpatialDimension;

            double[] x = new double[D];
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.CreateWrapper(x, 1, D);

            int NoOfEdges = GrdDatTmp.Edges.Count;
            for (int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                if (GrdDatTmp.Edges.IsEdgeBoundaryEdge(iEdge)) {
                    int jCell = GrdDatTmp.Edges.CellIndices[iEdge, 0];
                    int iFace = GrdDatTmp.Edges.FaceIndices[iEdge, 0];

                    var Cell_j = this.Cells[jCell];
                    var KRef = GrdDatTmp.Cells.GetRefElement(jCell);

                    //var LocalVerticesIn = RefElm.FaceCenters.ExtractSubArrayShallow(new int[] { iFace, 0 }, new int[] { iFace, D - 1 });

                    GrdDatTmp.TransformLocal2Global(KRef.GetFaceCenter(iFace), GlobalVerticesOut, jCell);


                    byte et = EdgeTagFunc(x);
                    //if (et <= 0)
                    //    throw new ApplicationException("edge tag 0 is reserved for inner edges.");
                    if (et >= FIRST_PERIODIC_BC_TAG)
                        throw new ApplicationException("edge tags greater or equal to " + FIRST_PERIODIC_BC_TAG + " are reserved for periodic \"boundaries\"");
                    if (!m_EdgeTagNames.ContainsKey(et))
                        throw new ArgumentException("unable to find EdgeTagName for EdgeTag = " + et);

                    CellFaceTag CFT = new CellFaceTag() {
                        EdgeTag = et,
                        FaceIndex = iFace,
                        NeighCell_GlobalID = long.MinValue
                    };
                    CFT.AddToArray(ref this.Cells[jCell].CellFaceTags);
                }
            }
        }

        /// <summary>
        /// This is a mapping from each used <em>EdgeTag</em> to a string that
        /// provides a name and additional information about the EdgeTag. The
        /// intention for this member is to provide both, a name (e.g.
        /// 'Left wall') for different regions of the boundary as well as
        /// boundary condition type info (e.g. 'inlet' or 'wall' or 'outflow' ...).
        /// </summary>
        /// <remarks>
        /// The names have no impact on the application on this application
        /// layer (L2-layer of BoSSS). They may be used on a higher application
        /// layer; Usually, this member (as like mostly all other public
        /// variable of this class) should be initialized by grid generator
        /// programs (see <see cref="DefineEdgeTags"/>).
        /// </remarks>
        public IDictionary<byte, string> EdgeTagNames {
            get {
                return m_EdgeTagNames;
            }
        }

        /// <summary>
        /// <see cref="EdgeTagNames"/>
        /// </summary>
        [DataMember]
        SortedList<byte, string> m_EdgeTagNames = new SortedList<byte, string>();

        /// <summary>
        /// adds a predefined partitioning to the grid 
        /// </summary>
        /// <param name="id">
        /// The unique name of the partitioning
        /// </param>
        /// <param name="cellToRankMap">
        /// For each cell, the MPI rank which should own the cell;
        /// </param>
        public void AddPredefinedPartitioning(string id, int[] cellToRankMap) {
            if (cellToRankMap.Length != this.NoOfUpdateCells)
                throw new ArgumentException("wrong length", "part");

            GridPartitioningVector gr;
            gr.Guid = Guid.NewGuid().MPIBroadcast(0);
            gr.CellToRankMap = cellToRankMap;

            m_PredefinedGridPartitioning.Add(id, gr);
        }

        /// <summary>
        /// Adds a predefined partitioning to the grid. Here, the partitioning
        /// can be based on geometrical information in the form of cell centers
        /// </summary>
        /// <param name="id"></param>
        /// <param name="PartFunc">
        /// A function, which maps the center-of-gravity of each cell to a
        /// MPI-rank
        /// </param>
        public void AddPredefinedPartitioning(string id, Func<double[], int> PartFunc) {
            int J = this.Cells.Length;
            int[] part = new int[J];
            int D = this.SpatialDimension;

            double[] Center = new double[D];
            for (int j = 0; j < J; j++) {
                var Cell_j = this.Cells[j];

                int NN = Cell_j.TransformationParams.GetLength(0);
                for (int d = 0; d < D; d++) {
                    Center[d] = 0.0;

                    for (int nn = 0; nn < NN; nn++) {
                        Center[d] += Cell_j.TransformationParams[nn, d];
                    }
                    Center[d] *= (1.0 / ((double)NN));
                }

                part[j] = PartFunc(Center);
            }

            this.AddPredefinedPartitioning(id, part);
        }

        /// <summary>
        /// values of <see cref="m_PredefinedGridPartitioning"/>
        /// </summary>
        [Serializable]
        public struct GridPartitioningVector {

            /// <summary>
            /// Guid assigned to <see cref="CellToRankMap"/>
            /// </summary>
            public Guid Guid;

            /// <summary>
            /// partitioning vector
            /// </summary>
            [NonSerialized]
            public int[] CellToRankMap;
        }

        /// <summary>
        /// See <see cref="PredefinedGridPartitioning"/>.
        /// </summary>
        [DataMember]
        internal SortedList<string, GridPartitioningVector> m_PredefinedGridPartitioning =
            new SortedList<string, GridPartitioningVector>();

        /// <summary>
        /// - keys: a string 'id' which identifies the partitioning
        /// - entries: Vectors and their Guids that contain the partitioning for 'id'.
        /// </summary>
        public IDictionary<string, GridPartitioningVector> PredefinedGridPartitioning {
            get {
                return m_PredefinedGridPartitioning;
            }
        }

        /// <summary>
        /// see <see cref="GridGuid"/>;
        /// </summary>
        [DataMember]
        Guid m_GridGuid;

        /// <summary>
        /// Guid which identifies this grid;
        /// </summary>
        public Guid GridGuid {
            get {
                return m_GridGuid;
            }
        }

        /// <summary>
        /// This method enables the user to manually override
        /// the Guid of this grid object (<see cref="GridGuid"/>).
        /// This may be useful in some situations (e.g. the user wants to use 
        /// IO for fields, but he also wants to create a new grid object every
        /// time he starts the application, instead of loading the grid from
        /// disk) but should be used with care, because the IO may become
        /// confused.
        /// </summary>
        /// <param name="NewGridGuid">the new Guid</param>
        public void SetGridGuid(Guid NewGridGuid) {
            m_GridGuid = NewGridGuid;
        }

        /// <summary>
        /// see <see cref="BcCellsStorageGuid"/>.
        /// </summary>
        [DataMember]
        private Guid m_BcCellsStorageGuid = Guid.Empty;

        /// <summary>
        /// see <see cref="StorageGuid"/>.
        /// </summary>
        [DataMember]
        private Guid m_StorageGuid;

        /// <summary>
        /// a string to store some user-information about the grid;
        /// </summary>
        [DataMember]
        public string Description;

        /// <summary>
        /// returns the Guid of the vector in which
        /// <see cref="Cells"/> is stored in the database. (see <see cref="BoSSS.Foundation.IO.DatabaseDriver.SaveGrid"/>);
        /// </summary>
        public Guid StorageGuid {
            get {
                return m_StorageGuid;
            }
            internal set {
                m_StorageGuid = value;
            }
        }

        /// <summary>
        /// returns the Guid of the vector in which
        /// <see cref="BcCells"/> is stored in the database.
        /// </summary>
        public Guid BcCellsStorageGuid {
            get {
                return m_BcCellsStorageGuid;
            }
            internal set {
                m_BcCellsStorageGuid = value;
            }
        }

        /// <summary>
        /// encodes the beginning of the edge tag range which is dedicated
        /// to periodic boundary conditions.
        /// </summary>
        public const byte FIRST_PERIODIC_BC_TAG = 181;

        /// <summary>
        /// Periodic boundary conditions are treated by connecting an "outlet"
        /// with some "inlet". Beside the cell neighborship relations (see
        /// <see cref="GridData.CellData.CellNeighbours"/>) an linear
        /// transformation must be provided which maps an affine-linear
        /// manifold A (the "outlet") of dimension D-1 (D denotes the spatial
        /// dimension) to another affine-linear manifold B (the "inlet"). (Note
        /// that the terms "outlet" and "inlet" are exchangeable.)
        /// </summary>
        /// <param name="X1">
        /// D pairwise different vectors, each d-dimensional, that specify the
        /// affine-linear manifold A;
        /// </param>
        /// <param name="N1">
        /// normal onto manifold A
        /// </param>
        /// <param name="X2">
        /// the image (i.e. the result when the transformation is applied onto
        /// the vectors <paramref name="X1"/>) of the (unknown) 
        /// transformation in manifold B; 
        /// </param>
        /// <param name="N2">
        /// normal onto manifold B
        /// </param>
        /// <param name="PeriodicTrafo_Tag">
        /// The edge tag for the periodic transformation
        /// </param>
        public void ConstructPeriodicEdgeTrafo(double[][] X1, double[] N1, double[][] X2, double[] N2, out byte PeriodicTrafo_Tag) {

            // check for right usage
            // ---------------------

            //if (m_Context.m_Fields != null)
            //    throw new ApplicationException("it's not allowed to call this method after Context.Setup(..) has been called.");

            int D = SpatialDimension;

            if (X1.Length != D)
                throw new ArgumentException("must contain exactly " + D + " elements in " + D + "D.", "x");
            if (X2.Length != D)
                throw new ArgumentException("must contain exactly " + D + " elements in " + D + "D.", "y");
            for (int i = 0; i < D; i++) {
                if (X1[i].Length != D)
                    throw new ArgumentException("vectors must be " + D + "-dimensional.", "x");
                if (X2[i].Length != D)
                    throw new ArgumentException("vectors must be " + D + "-dimensional.", "y");
            }
            if (N1.Length != D)
                throw new ArgumentException("vectors must be " + D + "-dimensional.", "x");
            if (N2.Length != D)
                throw new ArgumentException("vectors must be " + D + "-dimensional.", "y");

            MultidimensionalArray preImage = MultidimensionalArray.Create(D + 1, D);
            MultidimensionalArray image = MultidimensionalArray.Create(D + 1, D);
            for (int i = 0; i < D; i++) {
                for (int d = 0; d < D; d++) {
                    preImage[i, d] = X1[i][d];
                    image[i, d] = X2[i][d];
                }
            }
            for (int d = 0; d < D; d++) {
                preImage[D, d] = N1[d] + X1[0][d];
                image[D, d] = N2[d] + X2[0][d];
            }
            AffineTrafo pet = AffineTrafo.FromPoints(preImage, image);


            int Tag = m_PeriodicTrafo.Count + FIRST_PERIODIC_BC_TAG;
            if (Tag >= (byte.MaxValue - 1)) // rem: tag=255 is reserved
                throw new ApplicationException("Can't handle more than " + (byte.MaxValue - FIRST_PERIODIC_BC_TAG + 1) + "periodic boundary conditions.");
            PeriodicTrafo_Tag = (byte)Tag;
            m_PeriodicTrafo.Add(pet);
        }

        /// <summary>
        /// list of transformations which describe how some edges should be
        /// transformed to other edges;
        /// </summary>
        /// <remarks>
        /// indices into this list are the
        /// <see cref="GridData.EdgeData.EdgeTags"/> minus
        /// <see cref="FIRST_PERIODIC_BC_TAG"/>
        /// </remarks>
        public IList<AffineTrafo> PeriodicTrafo {
            get {
                return m_PeriodicTrafo;
            }
        }

        /// <summary>
        /// <see cref="PeriodicTrafo"/>
        /// </summary>
        [DataMember]
        internal List<AffineTrafo> m_PeriodicTrafo = new List<AffineTrafo>();

        /// <summary>
        /// inverse mappings to <see cref="PeriodicTrafo"/>;
        /// </summary>
        public IList<AffineTrafo> InversePeriodicTrafo {
            get {
                if (m_InversePeriodicTrafo == null || m_InversePeriodicTrafo.Count != m_PeriodicTrafo.Count) {
                    m_InversePeriodicTrafo = m_PeriodicTrafo.Select(trafo => trafo.Invert()).ToList();
                }
                return m_InversePeriodicTrafo;
            }
        }

        /// <summary>
        /// <see cref="PeriodicTrafo"/>
        /// </summary>
        [NonSerialized]
        internal List<AffineTrafo> m_InversePeriodicTrafo;

        /// <summary>
        /// The number of locally updated (on this MPI process)
        /// cells (i.e. internal and border cells, but not external cells).
        /// </summary>
        public int NoOfUpdateCells {
            get {
                return Cells.GetLength(0);
            }
        }

        /// <summary>
        /// Number of boundary condition cells over all MPI processes.
        /// </summary>
        public int NoOfBcCells {
            get {
                int J_BC = this.BcCells != null ? this.BcCells.Length : 0;
                return J_BC;
            }
        }

        /// <summary>
        /// class type of entry of <see cref="m_RefElements"/>
        /// </summary>
        [DataMember]
        string[] m_ClassNameOfRefElement;

        /// <summary>
        /// class type of entry of <see cref="m_EdgeRefElements"/>
        /// </summary>
        [DataMember]
        string[] m_ClassNameOfEdgeRefElement;

        /// <summary>
        /// the total number of cells;
        /// Since this value is required for the loading of the cells (in order
        /// to compute an MPI-partition) this value is actually stored
        /// redundantly during serialization.
        /// a negative value indicates undefined;
        /// </summary>
        [DataMember]
        long m_TotalNumberOfCells = -1;

        /// <summary>
        /// The total number of boundary condition;
        /// Since this value is required for the loading of the cells (in order
        /// to compute an MPI-partition) this value is actually stored
        /// redundantly during serialization.
        /// a negative value indicates undefined;
        /// </summary>
        [DataMember]
        int m_TotalNumberOfBcCells = 0; // 0 should indicate "not initialized", since this member is newly introduced
        //                                 and there will be lots of grids with the default value 0



        internal void InitNumberOfCells() {
            m_TotalNumberOfCells = this.Cells.Length.MPISum();
            m_TotalNumberOfBcCells = (this.BcCells != null ? this.BcCells.Length : 0).MPISum() + 1;
        }

        /// <summary>
        /// The total number of cells in all MPI processes.
        /// Required by <see cref="BoSSS.Foundation.IO.DatabaseDriver.LoadGridData"/> to distribute the cells 
        /// over all given processors;
        /// </summary>
        public long NumberOfCells_l {
            get {
                if (m_TotalNumberOfCells < 0) {
                    if (this.Cells == null)
                        throw new ApplicationException("non-initialized member.");
                    else
                        m_TotalNumberOfCells = this.Cells.Length.MPISum();
                }

                return m_TotalNumberOfCells;
            }
        }

        /// <summary>
        /// The total number of boundary-condition cells (<see cref="BcCells"/>) in all MPI processes.
        /// </summary>
        public int NumberOfBcCells {
            get {
                if (m_TotalNumberOfBcCells <= 0) {
                    // only to support legacy grids!
                    // The member did not exist when the grid was created/saved, so it may have the default value.
                    // However, this is a potential deadlock.
                    m_TotalNumberOfBcCells = (this.BcCells != null ? this.BcCells.Length : 0).MPISum() + 1;
                    //throw new ApplicationException("non-initialized member.");
                }

                return m_TotalNumberOfBcCells - 1; // 0 == not initialized => deswegen intern Anzahl + 1 speichern
            }
        }

        /// <summary>
        /// This method initializes <see cref="m_ClassNameOfRefElement"/>
        /// from <see cref="m_RefElements"/>.
        /// </summary>
        /// <param name="context"></param>
        [OnSerializing]
        private void OnSerializing(StreamingContext context) {
            Init_ClassNameOfSimplex();

            if (m_TotalNumberOfCells < 0) {
                throw new ApplicationException("Un-initialized data member (m_TotalNumberOfCells).");
            }
            if (m_TotalNumberOfBcCells <= 0) {
                // this test may fail on legacy grids (before the introduction of the 'm_TotalNumberOfBcCells'-member).
                //throw new ApplicationException("Un-initialized data member (m_TotalNumberOfBcCells).");
            }
        }

        /// <summary>
        /// this method initializes <see cref="m_ClassNameOfRefElement"/>
        /// from <see cref="m_RefElements"/>.
        /// </summary>
        private void Init_ClassNameOfSimplex() {
            // --------------
            m_ClassNameOfRefElement = new string[m_RefElements.Length];
            for (int i = 0; i < m_RefElements.Length; i++)
                m_ClassNameOfRefElement[i] = m_RefElements[i].GetType().AssemblyQualifiedName;
            // --------------
            m_ClassNameOfEdgeRefElement = new string[m_EdgeRefElements.Length];
            for (int i = 0; i < m_EdgeRefElements.Length; i++)
                m_ClassNameOfEdgeRefElement[i] = m_EdgeRefElements[i].GetType().AssemblyQualifiedName;
            // -------------


            //if (m_TotalNumberOfCells < 0) {
            //    long dummy = NumberOfCells_l;
            //}
            //if (m_TotalNumberOfBcCells <= 0) {
            //    int dummy = NumberOfBcCells;
            //}

        }

        /// <summary>
        /// see <see cref="OnSerializing"/>;
        /// </summary>
        /// <param name="context"></param>
        [OnDeserialized]
        private void AfterDeserialisation(StreamingContext context) {
            Init_GridSimplices();
        }

        /// <summary>
        /// Workaround to map full qualified type names (<paramref name="fqn"/>) from the time before 
        /// the December-16-refactoring to the respective types.
        /// This method is necessary to deserialize the old data.
        /// </summary>
        Type GetRefElementType(string fqn) {
            Type type = Type.GetType(fqn);

            if (type == null) {
                // map full-qualified names before the December-16-refactoring to the respective type.
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                if (fqn.StartsWith("BoSSS.Foundation.Grid.Point")) {
                    type = typeof(BoSSS.Foundation.Grid.RefElements.Point);
                } else if (fqn.StartsWith("BoSSS.Foundation.Grid.Line")) {
                    type = typeof(BoSSS.Foundation.Grid.RefElements.Line);
                } else if (fqn.StartsWith("BoSSS.Foundation.Grid.Triangle")) {
                    type = typeof(BoSSS.Foundation.Grid.RefElements.Triangle);
                } else if (fqn.StartsWith("BoSSS.Foundation.Grid.Square")) {
                    type = typeof(BoSSS.Foundation.Grid.RefElements.Square);
                } else if (fqn.StartsWith("BoSSS.Foundation.Grid.Tetra")) {
                    type = typeof(BoSSS.Foundation.Grid.RefElements.Tetra);
                } else if (fqn.StartsWith("BoSSS.Foundation.Grid.Cube")) {
                    type = typeof(BoSSS.Foundation.Grid.RefElements.Cube);
                } else {
                    throw new NotImplementedException("unknown type: " + fqn);
                }
            }

            return type;
        }



        /// <summary>
        /// this method initializes <see cref="m_RefElements"/>
        /// from <see cref="m_ClassNameOfRefElement"/>.
        /// </summary>
        private void Init_GridSimplices() {
            m_RefElements = new RefElement[m_ClassNameOfRefElement.Length];
            for (int i = 0; i < m_RefElements.Length; i++) {
                Type type = GetRefElementType(m_ClassNameOfRefElement[i]);
                PropertyInfo inst = type.GetProperty("Instance", BindingFlags.Static | BindingFlags.Public);
                m_RefElements[i] = (RefElement)(inst.GetValue(null, null));
            }
            // ----------------
            m_EdgeRefElements = new RefElement[m_ClassNameOfEdgeRefElement.Length];
            for (int i = 0; i < m_EdgeRefElements.Length; i++) {
                Type type = GetRefElementType(m_ClassNameOfEdgeRefElement[i]);
                PropertyInfo inst = type.GetProperty("Instance", BindingFlags.Static | BindingFlags.Public);
                m_EdgeRefElements[i] = (RefElement)(inst.GetValue(null, null));
            }
            // -----------------

        }

        /// <summary>
        /// reference elements for this grid;
        /// </summary>
        [NonSerialized]
        private RefElement[] m_RefElements;

        /// <summary>
        /// reference elements for this grid;
        /// </summary>
        [NonSerialized]
        private RefElement[] m_EdgeRefElements;

        /// <summary>
        /// specific reference element for this grid.
        /// </summary>
        public RefElement GetRefElement(int i) {
            return m_RefElements[i];
        }

        [NonSerialized]
        Dictionary<CellType, RefElement> m_CellTypeToRefElemnt;

        /// <summary>
        /// Checks that each cell type <see cref="Element.Type"/>
        /// of all <see cref="Cells"/> and <see cref="BcCells"/> in this grid can be mapped 
        /// to a reference element (<see cref="RefElements"/>) resp. an edge reference element (<see cref="EdgeRefElements"/>) in this grid.
        /// </summary>
        public void CheckCellTypes() {
            m_CellTypeToRefElemnt = new Dictionary<CellType, RefElement>();

            CellType[][] ct4kref = this.m_RefElements.Select(Kref => Kref.SupportedCellTypes.ToArray()).ToArray();

            List<int> problems = new List<int>();

            int J = this.Cells.Length;
            for (int j = 0; j < J; j++) {
                var cell = this.Cells[j];
                if (m_CellTypeToRefElemnt.ContainsKey(cell.Type))
                    continue;


                int itest = -1;
                for (int iKref = 0; iKref < ct4kref.Length; iKref++) {
                    itest = Array.IndexOf(ct4kref[iKref], cell.Type);
                    if (itest >= 0) {
                        m_CellTypeToRefElemnt.Add(cell.Type, this.m_RefElements[iKref]);
                        break;
                    }
                }

                if (itest < 0)
                    problems.Add(j);
            }

            string ErrString = null;
            if (problems.Count > 0) {
                ErrString = string.Format("Unable to map cell type for {0} of {1} cells in this grid to the reference element. First problem for cell #{2}, GlobalId={3}, cell type is {4}.", problems.Count, J, problems[0], this.Cells[problems[0]].GlobalID, this.Cells[problems[0]].Type);
            }

            if (this.BcCells != null) {
                CellType[][] ct4ekref = this.m_EdgeRefElements.Select(Kref => Kref.SupportedCellTypes.ToArray()).ToArray();
                List<int> eproblems = new List<int>();

                int JE = this.BcCells.Length;
                for (int j = 0; j < JE; j++) {
                    var bcCell = this.BcCells[j];
                    if (m_CellTypeToRefElemnt.ContainsKey(bcCell.Type))
                        continue;

                    int itest = -1;
                    for (int iKref = 0; iKref < ct4ekref.Length; iKref++) {
                        itest = Array.IndexOf(ct4ekref[iKref], bcCell.Type);
                        if (itest >= 0) {
                            m_CellTypeToRefElemnt.Add(bcCell.Type, this.m_RefElements[iKref]);
                            break;
                        }
                    }

                    if (itest < 0)
                        eproblems.Add(j);
                }

                if (eproblems.Count > 0) {
                    string eErrString = string.Format("Unable to map boundary-condition cell type for {0} of {1} cells in this grid to the reference element. First problem for boundary cell #{2}, GlobalId={3}, cell type is {4}.", eproblems.Count, JE, eproblems[0], this.BcCells[eproblems[0]].GlobalID, this.BcCells[eproblems[0]].Type);
                    if (ErrString == null)
                        ErrString = eErrString;
                    else
                        ErrString = ErrString + " Additionaly: " + eErrString;
                }
            }

            if (ErrString != null) {
                throw new ApplicationException(ErrString);
            }
        }


        /// <summary>
        /// finds the reference element for a specific cell type.
        /// </summary>
        public RefElement GetRefElement(CellType ct) {
            if (m_CellTypeToRefElemnt == null) {
                CheckCellTypes();
            }

            if (!m_CellTypeToRefElemnt.ContainsKey(ct))
                throw new ArgumentException(string.Format("Cell type '{0}' is not used in this grid.", ct));
            return m_CellTypeToRefElemnt[ct];
        }

        /// <summary>
        /// all reference elements in the grid
        /// </summary>
        public RefElement[] RefElements {
            get {
                return m_RefElements;
            }
        }

        /// <summary>
        /// all reference elements for edges in the grid.
        /// </summary>
        public RefElement[] EdgeRefElements {
            get {
                return m_EdgeRefElements;
            }
        }

        /// <summary>
        /// the simplex dimension in the sense of measure-theory
        /// </summary>
        public int SpatialDimension {
            get {
                return m_RefElements[0].SpatialDimension;
            }
        }

        [NonSerialized]
        Partitioning m_CellPartitioning;

        /// <summary>
        /// Gets the partition of cells over the MPI processes;
        /// </summary>
        public Partitioning CellPartitioning {
            get {
                if (m_CellPartitioning == null) {
                    m_CellPartitioning = new Partitioning(this.Cells.Length);
                }
                return m_CellPartitioning;
            }
        }

        [NonSerialized]
        Partitioning m_BcCellPartitioning;

        /// <summary>
        /// Gets the partition of boundary-condition cells over the MPI
        /// processes;
        /// </summary>
        public Partitioning BcCellPartitioning {
            get {
                if (m_BcCellPartitioning == null) {
                    m_BcCellPartitioning = new Partitioning(this.NoOfBcCells);
                }
                return m_BcCellPartitioning;
            }
        }

        /// <summary>
        /// Returns a dictionary that maps GlobalId's (on the current MPI
        /// processor) to local cell indices.
        /// </summary>
        /// <returns>
        /// a mapping from the GlobalID's to local cell indices;<br/>
        /// Keys of the dictionary: GlobalId of a locally stored cell (also
        /// external cell);<br/>
        /// Values of the dictionary: local cell indices;<br/>
        /// Contains only internal cells.
        /// </returns>
        public IDictionary<long, int> GetGlobalId2CellIndexMap() {
            //if (m_GlobalId2CellIndexMap != null)
            //    throw new ApplicationException("internal error: dictionary has already been created.");

            SortedDictionary<long, int> _GlobalId2CellIndexMap = new SortedDictionary<long, int>();

            int J = NoOfUpdateCells;
            for (int j = 0; j < J; j++) {
                _GlobalId2CellIndexMap.Add(Cells[j].GlobalID, j);
            }

            return _GlobalId2CellIndexMap;
        }

        /// <summary>
        /// Returns the current GlobalID - permutation of this gird,
        /// i.e. a mapping from the global index to the global ID.
        /// </summary>
        public Permutation GetGlobalIDPermutation(bool IncludeBcCells) {

            ilPSP.MPICollectiveWatchDog.Watch();
            int J = NoOfUpdateCells;
            int J_BC = IncludeBcCells ? NoOfBcCells : 0;

            long NoOfCells = this.NumberOfCells_l;

            long[] gid = new long[J + J_BC];

            for (int j = 0; j < J; j++) {
                if (Cells[j].GlobalID < 0 || Cells[j].GlobalID >= NoOfCells)
                    throw new NotSupportedException("Found illegal GlobalID for cell;");

                gid[j] = Cells[j].GlobalID;
            }
            for (int j = 0; j < J_BC; j++) {
                if (BcCells[j].GlobalID < NoOfCells)
                    throw new NotSupportedException("Found illegal GlobalID for boundary-condition cell;");

                gid[j + J] = BcCells[j].GlobalID;
            }

            return new Permutation(gid, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// The inverse of <see cref="GetGlobalIDPermutation"/>,
        /// i.e. a mapping from the global ID to global indices.
        /// </summary>
        public Permutation GetInverseGlobalIDPermutation(bool IncludeBcCells) {
            return GetGlobalIDPermutation(IncludeBcCells).Invert();
        }

        /// <summary>
        /// checks <see cref="Cells"/> for the occurrence of NAN or INF values.
        /// </summary>
        public void CheckGridForNANorINF() {
            using (new FuncTrace()) {
                // loop over all cells:
                int J = NoOfUpdateCells;
                for (int j = 0; j < J; j++) {
                    Cell Cj = Cells[j];

                    for (int i = Cj.TransformationParams.GetLength(1) - 1; i >= 0; i--) {
                        for (int k = Cj.TransformationParams.GetLength(0) - 1; k >= 0; k--) {
                            double x = Cj.TransformationParams[k, i];
                            if (double.IsInfinity(x))
                                throw new ApplicationException("grid corrupted: contains INF - value.");
                            if (double.IsNaN(x))
                                throw new ApplicationException("grid corrupted: contains NAN - value.");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// MPI process rank (within world communicator)
        /// </summary>
        public int MyRank {
            get {
                int rank;
                MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out rank);
                return rank;
            }
        }

        /// <summary>
        /// MPI world communicator size 
        /// </summary>
        public int Size {
            get {
                int size;
                MPI.Wrappers.csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out size);
                return size;
            }
        }

        [NonSerialized]
        private Partitioning m_NodePartitioning;

        /// <summary>
        /// the partition of cell nodes among MPI processes
        /// </summary>
        public Partitioning NodePartitioning {
            get {
                if (m_NodePartitioning == null) {

                    int NodeIdMin_loc = int.MaxValue;
                    int NodeIdMax_loc = int.MinValue;
                    int J = this.NoOfUpdateCells;
                    for (int j = 0; j < J; j++) {
                        //var Kref = this.GetRefElement(this.Cells[j].MajorCellTypeIndex);


                        foreach (int nodeId in this.Cells[j].NodeIndices) {
                            NodeIdMax_loc = Math.Max(NodeIdMax_loc, nodeId);
                            NodeIdMin_loc = Math.Min(NodeIdMin_loc, nodeId);
                        }
                    }

                    int J_BC = this.NoOfBcCells;
                    for (int j = 0; j < J_BC; j++) {
                        var Bc = this.BcCells[j];

                        foreach (int nodeId in Bc.NodeIndices) {
                            NodeIdMax_loc = Math.Max(NodeIdMax_loc, nodeId);
                            NodeIdMin_loc = Math.Min(NodeIdMin_loc, nodeId);
                        }
                    }

                    int NodeIdMin, NodeIdMax;
                    {
                        var MaxMin = MPIExtensions.MPIMax(new int[] { NodeIdMax_loc, -NodeIdMin_loc });
                        NodeIdMax = MaxMin[0];
                        NodeIdMin = -MaxMin[1];
                    }

                    NodeIdMin = Math.Min(NodeIdMin, 0);
                    if (NodeIdMin < 0) {

                        throw new ApplicationException("Illegal node indexing: minimal node index is " + NodeIdMin + "(must be non-negative).");

                    }

                    int iNode_0 = (NodeIdMax + 1) * MyRank / Size;
                    int iNode_E = (NodeIdMax + 1) * (MyRank + 1) / Size;

                    m_NodePartitioning = new Partitioning(iNode_E - iNode_0);
                }
                return m_NodePartitioning;
            }
        }

        /// <summary>
        /// Checks whether the given <paramref name="nodes"/> are in
        /// monotonically increasing order.
        /// </summary>
        /// <param name="nodes"></param>
        public static void CheckMonotonicity(double[] nodes) {
            if (nodes.Length < 0)
                throw new ArgumentException("length of nodes array must be at least 2;");
            double dprev = nodes[0];
            for (int i = 1; i < nodes.Length; i++) {
                if (nodes[i] <= dprev)
                    throw new ArgumentException("nodes must be strictly ascending.");
                dprev = nodes[i];
            }
        }

        /// <summary>
        /// clone
        /// </summary>
        virtual public object Clone() {
            var PrivateCtor = this.GetType().GetConstructor(BindingFlags.NonPublic | BindingFlags.Public | BindingFlags.Instance, null, new Type[0], null);
            if (PrivateCtor == null)
                throw new NotSupportedException("Grid must provide an empty constructor.");

            GridCommons Ret = PrivateCtor.Invoke(new object[0]) as GridCommons;

            Ret.m_RefElements = this.m_RefElements.CloneAs();
            Ret.m_EdgeRefElements = this.m_EdgeRefElements.CloneAs();

            Ret.m_GridGuid = Guid.NewGuid();
            Ret.Description = this.Description == null ? null : this.Description.CloneAs();
            Ret.Cells = new Cell[this.NoOfUpdateCells];
            if (this.PeriodicTrafo != null)
                Ret.m_PeriodicTrafo = new List<AffineTrafo>(this.PeriodicTrafo.Select(Tr => Tr.CloneAs()));

            for (int j = 0; j < Ret.Cells.Length; j++) {
                Ret.Cells[j] = this.Cells[j].CloneAs();
            }
            Ret.m_EdgeTagNames = new SortedList<byte, string>();
            foreach (var kv in this.m_EdgeTagNames)
                Ret.m_EdgeTagNames.Add(kv.Key, kv.Value.CloneAs());

            Ret.HackCheck();
            return Ret;
        }

        /// <summary>
        /// Helper function, used for debugging: transforms the whole grid by an affine-linear transformation <paramref name="T"/>.
        /// </summary>
        /// <param name="T">some transformation</param>
        public GridCommons Transform(AffineTrafo T) {
            ilPSP.MPICollectiveWatchDog.Watch();

            var ret = this.CloneAs();

            int J = this.Cells.Length;
            for (int j = 0; j < J; j++) {
                ret.Cells[j].TransformationParams = T.Transform(this.Cells[j].TransformationParams);
            }

            if (ret.m_PeriodicTrafo != null) {
                var Tinv = T.Invert();

                for (int i = 0; i < ret.m_PeriodicTrafo.Count; i++) {
                    ret.m_PeriodicTrafo[i] = T * this.m_PeriodicTrafo[i] * Tinv;
                }
            }

            return ret;
        }

        static private int[][] CompressIndexRangeParallel(int[][] IDX) {
            using(new FuncTrace()) {
                // Pseudo-Parallel: collect all info on process 0 and do the serial job

                int Rank, size;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out Rank);

                //SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);



                if(Rank > 0) {
                    Dictionary<int, int[][]> Packet = new Dictionary<int, int[][]>();
                    Packet.Add(0, IDX);

                    var ret = SerialisationMessenger.ExchangeData(Packet, csMPI.Raw._COMM.WORLD);
                    Debug.Assert(ret.Count == 0);

                    Packet.Clear();
                    var R = SerialisationMessenger.ExchangeData(Packet, csMPI.Raw._COMM.WORLD);
                    Debug.Assert(R.Count == 1);

                    return R[0];

                } else {

                    Dictionary<int, int[][]> Packet = new Dictionary<int, int[][]>();
                    var other_IDX = SerialisationMessenger.ExchangeData(Packet, csMPI.Raw._COMM.WORLD);
                    Debug.Assert(other_IDX.Count == size - 1);

                    //int TotalNoofPackets = IDX.Length + other_IDX.Values.Sum(idx => idx.Length);
                    List<int[]> All = new List<int[]>();
                    All.AddRange(IDX);
                    for(int p = 1; p < size; p++)
                        All.AddRange(other_IDX[p]);

                    var Compressed_IDX = CompressIndexRange(All.ToArray());

                    Dictionary<int, int[][]> Return = new Dictionary<int, int[][]>();
                    int i0 = IDX.Length;
                    for(int p = 1; p < size; p++) {
                        int L = other_IDX[p].Length;
                        var send_p = Compressed_IDX.GetSubVector(i0, L);

                        Return.Add(p, send_p);
                        i0 += L;
                    }

                    SerialisationMessenger.ExchangeData(Return, csMPI.Raw._COMM.WORLD);

                    return Compressed_IDX.GetSubVector(0, IDX.Length);
                }
            }
        }

        /// <summary>
        /// Index compression on a single processor.
        /// </summary>
        /// <param name="_IDX">
        /// Input, a collection of numbers, structured into a staggered array.
        /// </param>
        /// <returns>
        /// A staggered array with equal number of entries as <paramref name="_IDX"/>.
        /// </returns>
        static private int[][] CompressIndexRange(int[][] _IDX) {

            // Pass 1: Bereich bestimmen
            // =========================

            int min = int.MaxValue;
            int max = int.MinValue;
            int Negcnt = 0;
            foreach (int[] IDX in _IDX) {
                int L = IDX.Length;
                for (int i = 0; i < L; i++) {
                    int h = IDX[i];
                    if (h < 0) {
                        Negcnt++;
                    } else {
                        min = Math.Min(h, min);
                        max = Math.Max(h, max);
                    }
                }
            }

            // Pass 2: 'verwendete' einträge markiren
            // ======================================

            BitArray Bussi = new BitArray(max - min + 1);
            foreach (int[] IDX in _IDX) {
                int L = IDX.Length;
                for (int i = 0; i < L; i++)
                    Bussi[IDX[i] - min] = true;
            }

            // Pass 3: 'kompaktifizierte' Indices vergeben
            // ===========================================
            int[] Mapping = new int[max - min + 1];
            Mapping.SetAll(int.MinValue);

            int cnt = 0;
            for (int i = min; i <= max; i++) {
                if (Bussi[i - min]) {
                    Mapping[i - min] = cnt;
                    cnt++;
                }
            }


            //foreach (int[] IDX in _IDX) {
            //    int L = IDX.Length;
            //    for (int i = 0; i < L; i++) {
            //        int h = IDX[i];
            //        if (h >= 0) {
            //            if (Bussi[h - min]) {
            //                Mapping[h - min] = cnt;
            //                cnt++;
            //            }
            //        }
            //    }
            //}

            // Pass 4: Indices transformieren
            // ==============================

            var _Ret = new int[_IDX.Length][];

            for (int k = 0; k < _Ret.Length; k++) {
                int[] IDX = _IDX[k];
                int L = IDX.Length;
                int[] Ret = new int[L];
                for (int i = 0; i < L; i++) {
                    int h = IDX[i];
                    if (h < 0) {
                        Ret[i] = h;
                    } else {
                        Ret[i] = Mapping[h - min];
                        Debug.Assert(Ret[i] >= 0);
                    }
                }
                _Ret[k] = Ret;
            }

            return _Ret;
        }

        /// <summary>
        /// This function ensures that the GlobalID's of cells
        /// (<see cref="Cells"/>) and boundary-condition cells
        /// (<see cref="BcCells"/>) 
        /// start at 0 and 
        /// occupy a continuous region of natural numbers.
        /// </summary>
        public void CompressGlobalID(IList<long> AdditionalGlobalIdsToTransform = default(long[])) {
            List<int> oldGlobalID = new List<int>();
            List<int> oldCellFaceTagIDs = new List<int>();

            int J = this.Cells.Length;


            for (int j = 0; j < J; j++) {
                var Cj = this.Cells[j];
                oldGlobalID.Add((int)Cj.GlobalID);

                if (Cj.CellFaceTags != null) {
                    oldCellFaceTagIDs.AddRange(Cj.CellFaceTags.Select(cft => (int)cft.NeighCell_GlobalID));
                }
                
            }
            int J_BC = this.NoOfBcCells;
            for (int j = 0; j < J_BC; j++) {
                var Bj = this.BcCells[j];

                oldGlobalID.Add((int)Bj.GlobalID);

                if (Bj.NeighCell_GlobalIDs != null) {
                    oldCellFaceTagIDs.AddRange(Bj.NeighCell_GlobalIDs.Select(f => (int)f));
                }
            }
            if(AdditionalGlobalIdsToTransform != null) {
                for(int i = 0; i < AdditionalGlobalIdsToTransform.Count; i++) {
                    oldCellFaceTagIDs.AddRange((int)(AdditionalGlobalIdsToTransform[i]));
                }
            }



            int[][] New = CompressIndexRangeParallel(new int[][] { oldGlobalID.ToArray(), oldCellFaceTagIDs.ToArray() });
            int[] newGlobalID = New[0];
            int[] newCellFaceTagIDs = New[1];
            Debug.Assert(oldGlobalID.Count == newGlobalID.Length);
            Debug.Assert(oldCellFaceTagIDs.Count == newCellFaceTagIDs.Length);

            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

            int c2 = 0;
            for (int j = 0; j < J; j++) {
                var Cj = this.Cells[j];
                Cj.GlobalID = newGlobalID[j];

                if (Cj.CellFaceTags != null) {
                    int L = Cj.CellFaceTags.Length;
                    for (int l = 0; l < L; l++) {
                        if (Cj.CellFaceTags[l].EdgeTag > 0 && Cj.CellFaceTags[l].EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG) {
                            Cj.CellFaceTags[l].NeighCell_GlobalID = long.MinValue;
                        } else {
                            Debug.Assert(Cj.CellFaceTags[l].NeighCell_GlobalID == oldCellFaceTagIDs[c2]);
                            Cj.CellFaceTags[l].NeighCell_GlobalID = newCellFaceTagIDs[c2];
                        }
                        c2++;
                    }
                }
                
            }
            for (int j = 0; j < J_BC; j++) {
                var Bj = this.BcCells[j];
                Bj.GlobalID = newGlobalID[j + J];

                if (Bj.NeighCell_GlobalIDs != null) {
                    int L = Bj.NeighCell_GlobalIDs.Length;

                    for (int l = 0; l < L; l++) {
                        Debug.Assert(Bj.NeighCell_GlobalIDs[l] == oldCellFaceTagIDs[c2]);
                        Bj.NeighCell_GlobalIDs[l] = newCellFaceTagIDs[c2];
                        c2++;
                    }
                }
            }
            if(AdditionalGlobalIdsToTransform != null) {
                for(int i = 0; i < AdditionalGlobalIdsToTransform.Count; i++) {
                    AdditionalGlobalIdsToTransform[i] = newCellFaceTagIDs[c2];
                    c2++;
                }
            }

            Debug.Assert(c2 == newCellFaceTagIDs.Length);
        }

        /// <summary>
        /// This function ensures that the node indices (see <see cref="Element.NodeIndices"/>)
        /// start at 0 and 
        /// occupy a continuous region of natural numbers.
        /// </summary>
        public void CompressNodeIndices() {
            using (new FuncTrace()) {
                List<int> oldNodeIdx = new List<int>();

                int J = this.Cells.Length;
                for (int j = 0; j < J; j++) {
                    var Cj = this.Cells[j];
                    oldNodeIdx.AddRange(Cj.NodeIndices);
                }
                int J_BC = this.NoOfBcCells;
                for (int j = 0; j < J_BC; j++) {
                    var Bj = this.BcCells[j];
                    oldNodeIdx.AddRange(Bj.NodeIndices);
                }


                int[][] newIdx = CompressIndexRangeParallel(new int[][] { oldNodeIdx.ToArray() });

                var newNodeIdx = newIdx[0];
                int cnt = 0;
                for (int j = 0; j < J; j++) {
                    var Cj = this.Cells[j];

                    int L = Cj.NodeIndices.Length;
                    for (int l = 0; l < L; l++) {
                        Debug.Assert(oldNodeIdx[cnt] == Cj.NodeIndices[l]);
                        Cj.NodeIndices[l] = newNodeIdx[cnt];
                        cnt++;
                    }
                }
                for (int j = 0; j < J_BC; j++) {
                    var Bj = this.BcCells[j];

                    int L = Bj.NodeIndices.Length;
                    for (int l = 0; l < L; l++) {
                        Debug.Assert(oldNodeIdx[cnt] == Bj.NodeIndices[l]);
                        Bj.NodeIndices[l] = newNodeIdx[cnt];
                        cnt++;
                    }
                }
                Debug.Assert(cnt == newNodeIdx.Length);
            }
        }

        /// <summary>
        /// Merges duplicate nodes and checks that unique nodes are geometrically at the same position.
        /// </summary>
        public void MergeAndCheckNodes() {
            using (var tr = new FuncTrace()) {
                if (this.CellPartitioning.MpiSize > 1)
                    throw new NotImplementedException("No MPI parallel support.");

                int D = this.SpatialDimension;
                int J = this.Cells.Length;

                // 1st pass: find maximum node index & alloc memory
                // ================================================
                int NoOfNodes = 0;
                for (int j = 0; j < J; j++) { // loop over cells
                    Cell Cl = this.Cells[j];

                    NoOfNodes = Math.Max(NoOfNodes, Cl.NodeIndices.Max());
                }
                NoOfNodes++;
                MultidimensionalArray GlobalNodes = MultidimensionalArray.Create(NoOfNodes, D);
                MultidimensionalArray hLocal = MultidimensionalArray.Create(NoOfNodes);
                int[] NodesIndices = new int[NoOfNodes];
                NodesIndices.SetAll(-1);


                // 2nd pass: get coordinates of all nodes
                // ======================================

                for (int j = 0; j < J; j++) { // loop over cells
                    Cell Cl = this.Cells[j];
                    RefElement Kref = this.GetRefElement(Cl.Type);

                    int[] NodesType = Kref.GetInterpolationNodes_NodeType(Cl.Type);
                    int[] VertexIdx = Kref.GetInterpolationNodes_EntityIndices(Cl.Type);
                    Debug.Assert(NodesType.Max() == D); // nodes of type D should be the Vertices of the element.

                    // find the index of the first and the last vertex node
                    int i0 = -1, iE = -1;
                    for (int i = 0; i < NodesType.Length; i++) {
                        if (NodesType[i] == D) {
                            if (i0 < 0)
                                i0 = i;
                            if (iE < 0)
                                iE = i;
                            else
                                iE++;
                        }
                    }
                    Debug.Assert(i0 >= 0);
                    Debug.Assert(iE == NodesType.Length - 1);
                    Debug.Assert(VertexIdx[i0] == 0);
                    Debug.Assert(VertexIdx[iE] == Kref.NoOfVertices - 1);

                    // get vertices
                    var Vtx = Cl.TransformationParams.ExtractSubArrayShallow(new int[] { i0, 0 }, new int[] { iE, D - 1 });
                    double h = Vtx.MindistBetweenRows();

                    for (int ivtx = 0; ivtx < Kref.NoOfVertices; ivtx++) {
                        int iNode = Cl.NodeIndices[ivtx];
                        double[] NodeCoord = Vtx.GetRow(VertexIdx[ivtx]);

                        if (NodesIndices[iNode] < 0) {
                            GlobalNodes.SetRow(iNode, NodeCoord);
                            hLocal[iNode] = h;
                            NodesIndices[iNode] = iNode;
                        } else {
                            hLocal[iNode] = Math.Min(hLocal[iNode], h);
                            double[] NodeCoord2 = GlobalNodes.GetRow(iNode);
                            double dist = GenericBlas.L2Dist(NodeCoord, NodeCoord2);

                            if (dist >= hLocal[iNode] * 10e-5)
                                throw new ArithmeticException(string.Format("Error in grid: cell {0}, global node index {1}, geometric distance is {2}.", j, iNode, dist));
                        }
                    }
                }

                // 3rd pass: elmininate duplicates
                // ===============================
                {
                    int NoOfElimNodes = 0;


                    BoundingBox bb = new BoundingBox(GlobalNodes);
                    bb.ExtendByFactor(0.005);
                    int[] Perm = new int[GlobalNodes.GetLength(0)];
                    PointLocalization locTree = new PointLocalization(GlobalNodes, bb, Perm);

                    // GlobalNodes.GetRow(Perm[iNode]) == locTree.Points.GetRow(iNode)

                    //int[] InvPerm = new int[Perm.Length];
                    //for(int i = 0; i < InvPerm.Length; i++) {
                    //    InvPerm[Perm[i]] = i;
                    //}

                    // eliminate duplicate points
                    using (new BlockTrace("duplicate Point elimination", tr)) {
                        for (int iNode = 0; iNode < NoOfNodes; iNode++) {
                            if (NodesIndices[iNode] < 0)
                                // unused index
                                continue;
                            if (NodesIndices[iNode] < iNode)
                                // already merged with some other node
                                continue;
                            Debug.Assert(NodesIndices[iNode] == iNode);

                            double[] NodeCoord = GlobalNodes.GetRow(iNode);
                            double eps = hLocal[iNode] * 1.0e-6;
                            Debug.Assert(eps > 0);

                            List<int> foundPoints = new List<int>();
                            locTree.FindNearPoints(foundPoints, eps, NodeCoord);
                            if (foundPoints.Count < 1) {
                                throw new ApplicationException("error in algorithm");
                            } else if (foundPoints.Count == 1) {
                                if (Perm[foundPoints[0]] != iNode)
                                    throw new ApplicationException("error in algorithm");
                            } else {
                                foreach (int kk in foundPoints) {
                                    int iOtherNode = Perm[kk];
                                    Debug.Assert(iOtherNode >= iNode);
                                    NodesIndices[iOtherNode] = iNode;
                                }
                                NoOfElimNodes += (foundPoints.Count - 1);
                            }
                        }
                    }

                    Console.WriteLine("Eliminated " + NoOfElimNodes + " nodes.");

                }


                // 4th pass: renumber the nodes in the cells
                // =========================================
                for (int j = 0; j < J; j++) { // loop over cells
                    Cell Cl = this.Cells[j];

                    for (int iVtx = 0; iVtx < Cl.NodeIndices.Length; iVtx++) {
                        Cl.NodeIndices[iVtx] = NodesIndices[Cl.NodeIndices[iVtx]];
                    }
                }

            }
        }


        /// <summary>
        /// Checks that the Jacobian determinant of all cells is positive.
        /// Elements with purly negative Jacobian determinant are fixed by a mirror operation.
        /// If the Jacobian determinant flips sign in one cell, this is an un-recoverable error.
        /// </summary>
        public void CheckAndFixJacobianDeterminat() {
            using (var tr = new FuncTrace()) {

                List<int> unrecoverAbleErrors = new List<int>();

                Dictionary<CellType, Tuple<int[], int[]>> Permutations = new Dictionary<CellType, Tuple<int[], int[]>>();

                int D = this.SpatialDimension;
                int flipCount = 0;
                int J = this.Cells.Length;
                for (int j = 0; j < J; j++) {
                    Cell Cl = this.Cells[j];
                    RefElement Kref = this.GetRefElement(Cl.Type);

                    PolynomialList[] Jpolys = Kref.GetInterpolationPolynomials1stDeriv(Cl.Type);

                    NodeSet nds = Kref.GetQuadratureRule(Kref.GetInterpolationDegree(Cl.Type) * 2 + 2).Nodes;

                    MultidimensionalArray[] DerivEval = new MultidimensionalArray[D];
                    for (int d = 0; d < D; d++) {
                        DerivEval[d] = Jpolys[d].Values.GetValues(nds);
                    }

                    MultidimensionalArray Jacobi = MultidimensionalArray.Create(nds.NoOfNodes, D, D);

                    for (int d1 = 0; d1 < D; d1++) {
                        Debug.Assert(Cl.TransformationParams.GetLength(0) == Jpolys[d1].Count);
                        MultidimensionalArray JacobiCol = Jacobi.ExtractSubArrayShallow(-1, -1, d1);
                        JacobiCol.Multiply(1.0, DerivEval[d1], Cl.TransformationParams, 0.0, "kd", "kn", "nd");
                    }

                    bool signPos = false;
                    bool signNeg = false;
                    bool zero = false;
                    for (int k = 0; k < nds.NoOfNodes; k++) {
                        double JacobiDet = Jacobi.ExtractSubArrayShallow(k, -1, -1).Determinant();

                        if (JacobiDet > 0)
                            signPos = true;
                        if (JacobiDet < 0)
                            signNeg = true;

                        if (JacobiDet.Abs() <= 1.0e-20)
                            zero = true;
                    }
                    if (signPos == signNeg || zero) {
                        unrecoverAbleErrors.Add(j);
                        continue;
                    }

                    if (signNeg) {
                        // flip it!

                        flipCount++;

                        Tuple<int[], int[]> perms;
                        if (!Permutations.TryGetValue(Cl.Type, out perms)) {
                            int[] pNds, pvtx;
                            MirrorPermutation(Cl.Type, out pNds, out pvtx);
                            perms = new Tuple<int[], int[]>(pNds, pvtx);
                            Permutations.Add(Cl.Type, perms);
                        }

                        int[] NodesPerm = perms.Item1;
                        int[] VertxPerm = perms.Item2;

                        if (Cl.NodeIndices != null) {
                            Debug.Assert(Cl.NodeIndices.Length == VertxPerm.Length);
                            int[] NewNodeIndices = new int[Cl.NodeIndices.Length];
                            for (int k = 0; k < NewNodeIndices.Length; k++) {
                                NewNodeIndices[k] = Cl.NodeIndices[VertxPerm[k]];
                            }
                            Cl.NodeIndices = NewNodeIndices;
                        }

                        {
                            Debug.Assert(Cl.TransformationParams.GetLength(0) == NodesPerm.Length);
                            MultidimensionalArray newTransformationParams = MultidimensionalArray.Create(Cl.TransformationParams.Lengths);
                            for (int k = 0; k < NodesPerm.Length; k++) {
                                for (int d = 0; d < D; d++) {
                                    newTransformationParams[k, d] = Cl.TransformationParams[NodesPerm[k], d];
                                }
                            }
                            Cl.TransformationParams = newTransformationParams;
                        }

                    }
                }

                tr.Info("Number of flipped elements: " + flipCount + " of " + J);

                if (unrecoverAbleErrors.Count > 0) {
                    throw new ArithmeticException("Found " + unrecoverAbleErrors.Count + " cells with issues on Jacobian determinat, which is zero or flips sign. First problem occured in cell " + unrecoverAbleErrors.First() + ".");
                }
            }
        }



        /// <summary>
        /// For some reference element, this computes the permutation of nodes under a mirror operation.
        /// The purpose of this operation is to fix cellc with negative Jacobian determinante by mirroring them,
        /// since the mirror operation flips the sign of the Jacobian determinat.
        /// The reference emenet must be centered around the origin, so that the mirror operation maps it onto itself.
        /// </summary>
        /// <param name="ct"></param>
        /// <param name="NodesPerm">Permutation of interpolation nodes, see <see cref="RefElement.GetInterpolationNodes"/>.</param>
        /// <param name="VertxPerm">Permutation of reference element vertices, see <see cref="RefElement.Vertices"/>.</param>
        void MirrorPermutation(CellType ct, out int[] NodesPerm, out int[] VertxPerm) {
            RefElement Kref = this.GetRefElement(ct);
            NodeSet _Nodes = Kref.GetInterpolationNodes(ct);
            NodeSet _Vertx = Kref.Vertices;

            int D = this.SpatialDimension;

            var Mirror = MultidimensionalArray.Create(D, D);
            Mirror.AccEye(1.0);
            Mirror[0, 0] = -1;

            MultidimensionalArray TrfNodes = MultidimensionalArray.Create(_Nodes.Lengths);
            MultidimensionalArray TrfVertx = MultidimensionalArray.Create(_Vertx.Lengths);

            TrfNodes.Multiply(1.0, Mirror, _Nodes, 0.0, "kd", "dl", "kl");
            TrfVertx.Multiply(1.0, Mirror, _Vertx, 0.0, "kd", "dl", "kl");

            double Thres = 1e-6;

            {
                int K = _Nodes.NoOfNodes;
                bool[] mark = new bool[K];
                NodesPerm = new int[K];
                for (int k = 0; k < K; k++) {
                    double[] pt = TrfNodes.GetRow(k);
                    int i = _Nodes.FindRow(pt, Thres);

                    if (i < 0)
                        throw new ArithmeticException("Unable to compute mirror operation for element " + ct + ".");
                    if (mark[i] == true)
                        throw new ArithmeticException("Unable to compute mirror operation for element " + ct + ".");
                    mark[i] = true;

                    NodesPerm[k] = i;
                }
            }

            {
                int K = _Vertx.NoOfNodes;
                bool[] mark = new bool[K];
                VertxPerm = new int[K];
                for (int k = 0; k < K; k++) {
                    double[] pt = TrfVertx.GetRow(k);
                    int i = _Vertx.FindRow(pt, Thres);

                    if (i < 0)
                        throw new ArithmeticException("Unable to compute mirror operation for element " + ct + ".");
                    if (mark[i] == true)
                        throw new ArithmeticException("Unable to compute mirror operation for element " + ct + ".");
                    mark[i] = true;

                    VertxPerm[k] = i;
                }
            }
        }




        #region IEquatable<IGridInfo> Members

        /// <summary>
        /// Compares this grid to another.
        /// </summary>
        /// <param name="other">The grid to compare to.</param>
        /// <returns>true, if the grid GUIDs are the same; false otherwise.</returns>
        public bool Equals(IO.IGridInfo other) {
            if (other == null) {
                return false;
            } else {
                return this.ID.Equals(other.ID);
            }
        }

        #endregion

        /// <summary>
        /// Equality.
        /// </summary>
        public override bool Equals(object obj) {
            return this.Equals(obj as IO.IGridInfo);
        }

        /// <summary>
        /// Number of cells if used as hash code.
        /// </summary>
        public override int GetHashCode() {
            return this.ID.GetHashCode();
        }

        /// <summary>
        /// Needed for a simplified creation of a hanging nodes grid.
        /// Note: stores only a BoundingBox <see cref="BoundingBox"/> and the 
        /// number of Cells in x and y direction 
        /// </summary>
        public struct GridBox {

            public int[] numOfCells;
            public BoundingBox boundingBox;

            public GridBox(double[] pt1, double[] pt2, int cellsInX, int cellsInY, int cellsInZ) :
                this() {
                if (cellsInZ != 0)
                    numOfCells = new int[] { cellsInX, cellsInY, cellsInZ };
                else
                    numOfCells = new int[] { cellsInX, cellsInY };
                boundingBox = new BoundingBox(pt1, pt2);
            }
            public GridBox(double[] pt1, double[] pt2, int cellsInEachDirection) :
                this(pt1, pt2, cellsInEachDirection, cellsInEachDirection, cellsInEachDirection) {
            }
            public GridBox(double[] pt1, double[] pt2, int cellsInX, int cellsInY) :
                this(pt1, pt2, cellsInX, cellsInY, 0) {
            }
            public GridBox(double pt1x, double pt1y, double pt2x, double pt2y, int cellsInX, int cellsInY) :
                this(new double[] { pt1x, pt1y }, new double[] { pt2x, pt2y }, cellsInX, cellsInY, 0) {
            }
            public GridBox(double pt1x, double pt1y, double pt2x, double pt2y, int cellsInEachDirection) :
                this(new double[] { pt1x, pt1y }, new double[] { pt2x, pt2y }, cellsInEachDirection, cellsInEachDirection, 0) {
            }
            public GridBox(double pt1x, double pt1y, double pt2x, double pt2y, int cellsInX, int cellsInY, int cellsInZ) :
               this(new double[] { pt1x, pt1y }, new double[] { pt2x, pt2y }, cellsInX, cellsInY, cellsInZ) {
            }
        }

    }
}

