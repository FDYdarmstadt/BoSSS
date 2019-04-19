using BoSSS.Foundation.IO;
using ilPSP;
using MPI.Wrappers;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Aggregation {

    /// <summary>
    /// A grid, which is composed throug the aggregation of cells in a parent grid (<see cref="ParentGrid"/>).
    /// </summary>
    [Serializable]
    [DataContract]
    public partial class AggregationGrid : IGridInfo, ICloneable, IGrid {

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="pGrid">
        /// Parent grid.
        /// </param>
        /// <param name="AggregationCells">
        /// Coarse cells which build up the fine cells.
        /// - 1st index: coarse (i.e. this) grid cell index
        /// - 2nd index: enumeration
        /// - content: local cell index into the parent grid <paramref name="pGrid"/>.
        /// </param>
        /// <param name="AggregationCellGids">
        /// GlobalID's for each of the <paramref name="AggregationCells"/>;  
        /// Optional, can be null - in this case, GlobalID's are chosen automatically;
        /// If not null, length must be equal to 1st length of <paramref name="AggregationCells"/>.
        /// </param>
        public AggregationGrid(IGrid pGrid, int[][] AggregationCells, long[] AggregationCellGids = null) {
            m_ParentGrid = pGrid;
            m_ParentGridID = pGrid.ID;
            m_GridGuid = Guid.NewGuid();
            if (AggregationCellGids != null) {
                if (AggregationCellGids.Length != AggregationCells.Length)
                    throw new ArgumentException("array length mismatch");
            }

            m_CellPartitioning = new Partitioning(AggregationCells.Length);

            var gdat = new AggregationGridData(this, AggregationCells);

            int J = AggregationCells.Length;
            if (AggregationCellGids == null) {
                AggregationCellGids = new long[J];
                int i0 = m_CellPartitioning.i0;
                for (int j = 0; j < J; j++)
                    AggregationCellGids[j] = j + i0;
            }

            long[] parentGid = pGrid.iGridData.CurrentGlobalIdPermutation.Values;
            AggCells = new AggCell[J];
            for (int j = 0; j < J; j++) {
                int[] AggCell = AggregationCells[j];
                int L = AggCell.Length;
                AggCells[j].GlobalID = AggregationCellGids[j];
                AggCells[j].PartGlobalId = new long[L];
                long[] _PartGlobalId = AggCells[j].PartGlobalId;
                for (int iPart = 0; iPart < L; iPart++) {
                    _PartGlobalId[iPart] = parentGid[AggCell[iPart]];
                }
            }


            m_GridData = gdat;
        }

        protected AggregationGrid()
        {

        }

        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        /// <param name="EdgeTagFunc"></param>
        public void DefineEdgeTags(Func<double[], byte> EdgeTagFunc)
        {
            int D = SpatialDimension;
            double[] x = new double[D];
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.CreateWrapper(x, 1, D);
            for (int iEdge = 0; iEdge < m_GridData.iGeomEdges.Count; ++iEdge)
            {
                if (m_GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge))
                {
                    int jCell = m_GridData.iGeomEdges.CellIndices[iEdge, 0];
                    int iFace = m_GridData.iGeomEdges.FaceIndices[iEdge, 0];
                    var KRef = m_GridData.iGeomCells.GetRefElement(jCell);

                    m_GridData.TransformLocal2Global(KRef.GetFaceCenter(iFace), GlobalVerticesOut, jCell);
                    byte et = EdgeTagFunc(x);

                    if (!EdgeTagNames.ContainsKey(et))
                        throw new ArgumentException("unable to find EdgeTagName for EdgeTag = " + et);
                    m_GridData.iGeomEdges.EdgeTags[iEdge] = et;
                }
            }
        }

        /// <summary>
        /// Structure to store one Aggregation cell
        /// </summary>
        [Serializable]
        [DataContract]
        public struct AggCell {
            /// <summary>
            /// ID of this cell
            /// </summary>
            [DataMember]
            public long GlobalID;

            /// <summary>
            /// ID's of cells in the parent grid (<see cref="ParentGrid"/>) which make this aggregate cell.
            /// </summary>
            [DataMember]
            public long[] PartGlobalId;
        }

        /// <summary>
        /// returns the Guid of the vector in which
        /// <see cref="AggCells"/> is stored in the database. (see <see cref="BoSSS.Foundation.IO.DatabaseDriver.SaveGrid"/>);
        /// </summary>
        public Guid AggCellStorageGuid {
            get {
                return m_StorageGuid;
            }
            internal set {
                m_StorageGuid = value;
            }
        }

        /// <summary>
        /// see <see cref="AggCellStorageGuid"/>.
        /// </summary>
        [DataMember]
        private Guid m_StorageGuid;


        [NonSerialized]
        [JsonIgnore]
        AggCell[] AggCells;

        public AggCell[] AggregationCells{ get => AggCells; set => AggCells = value; }

        /// <summary>
        /// a string to store some user-information about the grid;
        /// </summary>
        [DataMember]
        public string Description {
            get {
                return m_Description;
            }
            set {
                m_Description = value;
                if (Database != null) {
                    Database.Controller.SaveGridInfo(this);
                }
            }
        }

        [JsonProperty(PropertyName = "Aggregation_Description")]
        [DataMember]
        string m_Description;

        /// <summary>
        /// parent grid object
        /// </summary>
        public IGrid ParentGrid {
            get {
                return m_ParentGrid;
            }
        }

        [DataMember]
        IGrid m_ParentGrid;

        [DataMember]
        Guid m_ParentGridID;

        /// <summary>
        /// Guid (<see cref="IDatabaseEntityInfo{T}.ID"/> of parent grid <see cref="ParentGrid"/>
        /// </summary>
        public Guid ParentGridID {
            get {
                return m_ParentGridID;
            }
        }


        /// <summary>
        /// Either 1 for 1D, 2 for 2D or 3 for 3D.
        /// </summary>
        public int SpatialDimension {
            get {
                return m_ParentGrid.SpatialDimension;
            }
        }

        public IReadOnlyCollection<Guid> AllDataVectorIDs => throw new NotImplementedException();

        public IDictionary<byte, string> EdgeTagNames {
            get {
                return m_ParentGrid.EdgeTagNames;
            }
        }

        /// <summary>
        /// see <see cref="ID"/>;
        /// </summary>
        [DataMember]
        Guid m_GridGuid;

        /// <summary>
        /// Guid/Identification of this grid object in the database <see cref="Database"/>
        /// </summary>
        public Guid ID {
            get {
                return this.m_GridGuid;
            }
        }

        public DateTime CreationTime => throw new NotImplementedException();

        public DateTime WriteTime {
            get;
            set;
        }


        /// <summary>
        /// grid name: implementation of <see cref="IDatabaseEntityInfo{T}.Name"/>
        /// </summary>
        public string Name {
            get {
                return m_Name;
            }
            set {
                if (String.IsNullOrWhiteSpace(value) == false) {
                    m_Name = value.Trim();
                    if (Database != null) {
                        Database.Controller.SaveGridInfo(this);
                    }
                } else {
                    throw new Exception("New name of grid is invalid.");
                }
            }
        }

        [DataMember]
        private string m_Name;

        [DataMember]
        int m_NumberOfCells = -1;

        internal void InitNumberOfCells()
        {
            m_NumberOfCells = this.AggCells.Length.MPISum();
        }

        /// <summary>
        /// number of cells in the grid: implementation of <see cref="IGridInfo.NumberOfCells"/>
        /// </summary>
        public int NumberOfCells {
            get {
                if (m_NumberOfCells < 0) {
                    if (this.AggCells == null)
                        throw new ApplicationException("non-initialized member.");
                    else
                        m_NumberOfCells = this.AggCells.Length.MPISum();
                }

                return m_NumberOfCells; ;
            }
        }

        [NonSerialized]
        internal IO.IDatabaseInfo m_Database = null;

        /// <summary>
        /// Database which de-serialized this grid; implementation of <see cref="IDatabaseEntityInfo{t}.Database"/>
        /// </summary>
        public IDatabaseInfo Database {
            get {
                return m_Database;
            }
            set {
                if (value != null) {
                    m_Database = value;
                } else {
                    throw new ArgumentNullException();
                }
            }
        }

        [NonSerialized]
        AggregationGridData m_GridData;

        void InitGridData() {
            if (m_GridData != null)
                return;

            if(this.Size > 1) {
                Console.WriteLine("Warning: will probably not work in parallel");
            }

            // compute mapping: 
            // globalID -> Global Index for parent grid
            // =================================================
            var ParentGids = this.ParentGrid.iGridData.CurrentGlobalIdPermutation;
            var ParentIdx = ParentGids.Invert();

            // get parent grid indices for aggregation cells
            // =================================================
            int[][] AggIdx;
            {
                int J = this.AggCells.Length;
                AggIdx = new int[J][];

                int L = 0;
                for (int j = 0; j < J; j++) {
                    L += this.AggCells[j].PartGlobalId.Length;
                }

                long[] InputBuffer = new long[L];
                int l = 0;
                for (int j = 0; j < J; j++) {
                    long[] src = this.AggCells[j].PartGlobalId;
                    Array.Copy(src, 0, InputBuffer, l, src.Length);
                    l += src.Length;
                }

                long[] EvalBuffer = new long[L];
                ParentIdx.EvaluatePermutation(InputBuffer, EvalBuffer);

                int i0_Parent = ParentGrid.CellPartitioning.i0;
                int J_Parent = ParentGrid.iGridData.iLogicalCells.NoOfLocalUpdatedCells;

                l = 0;
                for (int j = 0; j < J; j++) {
                    int AG = this.AggCells[j].PartGlobalId.Length;
                    int[] AggIdx_j = new int[AG];
                    AggIdx[j] = AggIdx_j;

                    for (int k = 0; k < AG; k++) { // loop over parts of aggregation cell
                        Debug.Assert(EvalBuffer[l] >= 0);
                        Debug.Assert(EvalBuffer[l] < ParentGrid.NumberOfCells);
                        AggIdx_j[k] = (int)EvalBuffer[l] - i0_Parent;
                        Debug.Assert(AggIdx_j[k] >= 0);
                        Debug.Assert(AggIdx_j[k] < J_Parent);

                        l++;
                    }
                }
            }

            // create grid data
            // ===================
            m_GridData = new AggregationGridData(this, AggIdx);

        }


        /// <summary>
        /// returns the grid metrics, of type <see cref="AggregationGridData"/>
        /// </summary>
        public IGridData iGridData {
            get {
                if(m_GridData ==  null) {
                    InitGridData();
                }
                return m_GridData;
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

        /// <summary>
        /// Gets the partition of cells over the MPI processes;
        /// </summary>
        public Partitioning CellPartitioning {
            get {
                if (m_CellPartitioning == null) {
                    m_CellPartitioning = new Partitioning(this.AggCells.Length);
                }
                return m_CellPartitioning;
            }
        }

        [NonSerialized]
        Partitioning m_CellPartitioning;


        public object Clone() {
            throw new NotImplementedException();
        }

        public IGridInfo CopyFor(IDatabaseInfo targetDatabase) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Equality.
        /// </summary>
        public bool Equals(IGridInfo other) {
            if (other == null)
            {
                return false;
            }
            else
            {
                return this.ID == other.ID;
            }
        }

        public void Redistribute(IDatabaseDriver iom, GridPartType method, string PartOptions) {
            if (Size <= 1)
                return; // nothing to do

            throw new NotImplementedException();
        }

        public void RedistributeGrid(int[] part) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Returns the current GlobalID - permutation of this grid,
        /// i.e. a mapping from the global index to the global ID.
        /// </summary>
        public Foundation.Comm.Permutation CurrentGlobalIdPermutation()
        {
            ilPSP.MPICollectiveWatchDog.Watch();
            int NumberOfLocalMPICells = AggCells.Length;

            long NoOfCells = this.NumberOfCells;

            long[] gid = new long[NumberOfLocalMPICells];

            for (int j = 0; j < NumberOfLocalMPICells; j++)
            {
                if (AggCells[j].GlobalID < 0 || AggCells[j].GlobalID >= NoOfCells)
                    throw new NotSupportedException("Found illegal GlobalID for cell;");

                gid[j] = AggCells[j].GlobalID;
            }

            return new Foundation.Comm.Permutation(gid, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
        }
    }
}
