using System;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// An unstructured 3D-grid consisting of tetrahedral elements
    /// </summary>
    [Serializable]
    public sealed class UnstructuredTetraGrid : GridCommons {

        /// <summary>
        /// Constructs an empty grid. All public member variables must be
        /// initialized before this object can be used in BoSSS.
        /// </summary>
        public UnstructuredTetraGrid()
            : base() {
            m_GridSimplex = new Tetra();
        }

        /// <summary>
        /// The standard constructor. Checks the array dimensions and
        /// initializes the grid parameters.
        /// </summary>
        /// <param name="GlobalID">see <see cref="GridCommons.GlobalID"/></param>
        /// <param name="CellNeighbours">see <see cref="GridCommons.CellNeighbours"/>;</param>
        /// <param name="Vertices">see <see cref="GridCommons.Vertices"/>;</param>
        /// <param name="EdgeTags">see <see cref="GridCommons.EdgeTags"/>;</param>
        public UnstructuredTetraGrid(long[] GlobalID, double[, ,] Vertices, long[,] CellNeighbours, byte[,] EdgeTags)
            : this() {

            if (GlobalID.GetLength(0) != Vertices.GetLength(0))
                throw new ArgumentException("mismatch in dimension 0", "GlobalID,Vertices");
            if (GlobalID.GetLength(0) != CellNeighbours.GetLength(0))
                throw new ArgumentException("mismatch in dimension 0", "GlobalID,CellNeighbours");
            if (Vertices.GetLength(1) != m_GridSimplex.Vertices.GetLength(0))
                throw new ArgumentException("dimension 1; wrong number of vertices per element", "Vertices");
            if (Vertices.GetLength(2) != 3)
                throw new ArgumentException("dimension 2; wrong spatial dimension", "Vertices");
            if (CellNeighbours.GetLength(1) != 4)
                throw new ArgumentException("dimension 1; wrong number of neighbours", "CellNeighbours");

            base.GlobalID = GlobalID;
            base.Vertices = Vertices;
            base.CellNeighbours = CellNeighbours;
            base.EdgeTags = EdgeTags;
        }
    }
}
