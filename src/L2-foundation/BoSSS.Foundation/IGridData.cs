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

using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {
    public interface IGridData {

        /// <summary>
        /// MPI process rank (within world communicator)
        /// </summary>
        int MpiRank {
            get;
        }

        /// <summary>
        /// MPI world communicator size 
        /// </summary>
        int MpiSize {
            get;
        }

        /// <summary>
        /// Identification of the grid in the BoSSS database, 
        /// equal to the <see cref="BoSSS.Foundation.IO.IDatabaseEntityInfo{T}.ID"/>.
        /// </summary>
        Guid GridID {
            get;
        }

        /// <summary>
        /// This is a mapping from each used <em>EdgeTag</em>, (see <see cref="IGeometricalEdgeData.EdgeTags"/>) to a string that
        /// provides a name and additional information about the EdgeTag. The
        /// intention for this member is to provide both, a name (e.g.
        /// 'Left wall') for different regions of the boundary as well as
        /// boundary condition type info (e.g. 'inlet' or 'wall' or 'outflow' ...).
        /// </summary>
        IDictionary<byte, string> EdgeTagNames {
            get;
        }


        IGeometricalCellsData iGeomCells { get; }

        ILogicalCellData iLogicalCells { get; }

        /// <summary>
        /// Information about the vertices of the grid elements, see
        /// <see cref="IVertexData"/>.
        /// </summary>
        IVertexData iVertices {
            get;
        }

        /// <summary>
        /// metrics and operations which are associated to edges
        /// </summary>
        IGeometricalEdgeData iGeomEdges {
            get;
        }

        /// <summary>
        /// metrics and operations which are associated to edges
        /// </summary>
        ILogicalEdgeData iLogicalEdges {
            get;
        }

        /// <summary>
        /// see <see cref="Parallelization"/>
        /// </summary>
        IParallelization iParallel {
            get;
        }

        /// <summary>
        /// The spatial dimension of the grid (usually 1, 2 or 3).
        /// </summary>
        int SpatialDimension { get; }

        /// <summary>
        /// Evaluation of the DG basis on the grid
        /// </summary>
        Grid.BasisData ChefBasis {
            get;
        }

        /// <summary>
        /// Gets the partitioning of cells over the MPI processes;
        /// </summary>
        Partitioning CellPartitioning {
            get;
        }

        /// <summary>
        /// The global ID for each cell
        /// </summary>
        Comm.Permutation CurrentGlobalIdPermutation {
            get;
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
        void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, MultidimensionalArray GlobalVerticesOut, int jCell);


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
        void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, int j0, int Len, MultidimensionalArray GlobalVerticesOut, int OutArrayOffset);

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
        void TransformLocal2Global(MultidimensionalArray NS, int j0, int Len, MultidimensionalArray Nodesglob);


        /// <summary>
        /// Transforms vertices from the global coordinate system to the local coordinate systems 
        /// of geometrical cells <paramref name="j0"/> to <paramref name="j0"/>+<paramref name="Len"/>-1.
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
        void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int j0, int Len, int OutArrayOffset);

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
        /// Output; the vertices form <paramref name="GlobalVerticesIn"/>, transformed to local coordinates;
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
        void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int jCell, bool[] NewtonConvergence);



        /// <summary>
        /// Cached transformation of node sets to global coordinates.
        /// </summary>
        Caching.CacheLogic_CNs GlobalNodes {
            get;
        }

        /// <summary>
        /// Jacobian of transformation from reference to physical space, \f$ (\nabla T_j) \f$.
        /// </summary>
        Caching.CacheLogicImplBy_CNs Jacobian {
            get;
        }

        /// <summary>
        /// Adjungate of the Jacobian of the reference-to-physical coordinate transformation, 
        /// \f$ \mathrm{Adj}( \nabla T_j ) =  \determinant{ \nabla T_j } ( \nabla T_j )^{-1} \f$.
        /// </summary>
        Caching.CacheLogicImplBy_CNs AdjungateJacobian {
            get;
        }

        /// <summary>
        /// Inverse of the Jacobian of the reference-to-physical coordinate transformation, \f$ ( \nabla T_j )^{-1} \f$.
        /// </summary>
        Caching.CacheLogicImplBy_CNs InverseJacobian {
            get;
        }

        /// <summary>
        /// Determinant of the Jacobian of the reference-to-physical coordinate transformation, \f$ \determinant{ \nabla T_j } \f$.
        /// </summary>
        Caching.CacheLogicImplBy_CNs JacobianDeterminat {
            get;
        }

    }

    public interface IGeometricalCellsData {

        /// <summary>
        /// Number of geometrical cells.
        /// </summary>
        int Count {
            get;
        }

        /// <summary>
        /// Cell type for cell <paramref name="jCell"/>.
        /// </summary>
        CellType GetCellType(int jCell);

        /// <summary>
        /// For affine-linear cells,
        /// the absolute value of the (Jacobi) determinant of the 
        /// transformation from local cell coordinate system to global
        /// coordinate system.
        /// (see <see cref="EvaluateJacobian(NodeSet, int, int, MultidimensionalArray)"/>
        /// 1st index: local cell index;
        /// </summary>
        MultidimensionalArray JacobiDet {
            get;
        }

        /// <summary>
        /// Inverse matrices to <see cref="Transformation"/>
        /// </summary>
        /// <remarks>
        /// Indices are defined as follows:
        /// <list type="bullet">
        ///   <item>1st index: local cell index (locally updated and external cells);</item>
        ///   <item>2nd index: matrix row index;</item>
        ///   <item>3rd index: matrix column index;</item>
        /// </list>
        /// </remarks>
        MultidimensionalArray InverseTransformation {
            get;
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
        MultidimensionalArray Transformation {
            get;
        }


        /// <summary>
        /// See <see cref="CellInfo"/>.
        /// </summary>
        CellInfo[] InfoFlags {
            get;
        }

        
        /// <summary>
        /// Mapping from cells to vertices/nodes of the grid (stored in
        /// <see cref="IVertexData.Coordinates"/>) 
        /// - content: indices into <see cref="IVertexData.Coordinates"/> 
        /// - 1st index: local cell index (externals included) 
        /// - 2nd index: cell vertex index
        /// </summary>
        int[][] CellVertices {
            get;
        }

        /// <summary>
        /// Returns the reference element index (an index into <see cref="RefElements"/>) for cell <paramref name="jCell"/>.
        /// </summary>
        /// <param name="jCell">local cell index.</param>
        /// <returns>reference element index.</returns>
        int GetRefElementIndex(int jCell);

        /// <summary>
        /// Returns the reference element (one of <see cref="RefElements"/>) for cell <paramref name="j"/>.
        /// </summary>
        RefElement GetRefElement(int j);

        /// <summary>
        /// True, if the mapping from the reference/local coordinates to physical coordinates 
        /// </summary>
        /// <param name="jCell">local cell index.</param>
        /// <returns>reference element index.</returns>
        bool IsCellAffineLinear(int jCell);

        /// <summary>
        /// All reference elements for cells, see <see cref="GetRefElementIndex(int)"/> resp. <see cref="GetRefElement(int)"/>.
        /// </summary>
        RefElement[] RefElements {
            get;
        }

        /// <summary>
        /// Returns the volume (to be more exact: the
        /// <see cref="SpatialDimension"/> - dimensional measure) of the
        /// cell <paramref name="j"/>;
        /// </summary>
        /// <param name="j">local cell index</param>
        /// <returns></returns>
        double GetCellVolume(int j);

        /// <summary>
        /// Computes the bounding box of cell <paramref name="j"/>.
        /// </summary>
        /// <param name="j">local cell index.</param>
        /// <param name="bb">
        /// on exit, the bounding box of cell j.
        /// </param>
        void GetCellBoundingBox(int j, BoundingBox bb);


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
        int GetNoOfSimilarConsecutiveCells(CellInfo mask, int j0, int Lmax);

        /// <summary>
        /// The minimal Euclidean distance between two distinct vertices for each cell;
        /// (Can be used e.g. to compute the CFL number);
        /// - index: local cell index;
        /// </summary>
        MultidimensionalArray h_min {
            get;
        }

        /// <summary>
        /// The maximal Euclidean distance between two vertices for each cell;
        /// (Can be used e.g. to compute the CFL number);
        /// - index: local cell index;
        /// </summary>
        MultidimensionalArray h_max {
            get;
        }
        
        /// <summary>
        /// Cell-Mask of all geometric cells which share the same reference element.
        /// </summary>
        CellMask GetCells4Refelement(RefElement Kref);

        /// <summary>
        /// Mapping form geometrical to logical cells
        /// </summary>
        int[] GeomCell2LogicalCell {
            get;
        }

        /// <summary>
        /// polynomial interpolation degree of the Reference-to-Global coordinate transformation.
        /// </summary>
        /// <param name="jCell"></param>
        /// <returns></returns>
        int GetInterpolationDegree(int jCell);

        /// <summary>
        /// Center-of-gravity for the cell
        /// </summary>
        double[] GetCenter(int jCell);
    }

    public interface ILogicalCellData {

        /// <summary>
        /// Global identification (a cell index which remains constant under re-partitioning) for the <paramref name="j"/>-th cell.
        /// </summary>
        /// <param name="j">Local cell index.</param>
        /// <returns>Global Id.</returns>
        long GetGlobalID(int j);


        /// <summary>
        /// local indices of neighbor cells;
        ///  - 1st index: local cell index;
        ///  - 2nd index: enumeration
        /// </summary>
        int[][] CellNeighbours { get; }

        /// <summary>
        /// Number of locally updated cells - the cells which are computed on
        /// this processor (in contrast, see <see cref="NoOfExternalCells"/>);
        /// </summary>
        int NoOfLocalUpdatedCells {
            get;
        }

        /// <summary>
        /// Number of locally stored external cells - no computations are carried out for
        /// that cells, but their values are needed.
        /// </summary>
        int NoOfExternalCells {
            get;
        }

        /// <summary>
        /// <see cref="NoOfExternalCells"/> plus <see cref="NoOfLocalUpdatedCells"/>;
        /// </summary>
        int Count {
            get;
        }

        /// <summary>
        /// Mapping from logical cells to geometrical cells; only required if geometrical and logical cells are not identical.
        /// - 1st index: local (logical) cell index
        /// - 2nd index: enumeration of geometrical cells (parts)
        /// </summary>
        int[][] AggregateCellToParts {
            get;
        }

        /// <summary>
        /// Which edges (see <see cref="EdgeData.CellIndices"/>) bound to
        /// which cells? <br/>
        /// 1st index: local cell index <em>j</em>, only local updated<br/>
        /// 2nd index: collection, order is arbitrary; <br/>
        /// content <em>e</em>: 
        /// If <em>e</em> is positive, then cell <em>j</em> is the first
        /// (IN) cell of edge <em>e - 1</em>. Otherwise, if <em>e</em> is
        /// negative, then cell <em>j</em> is the second (OUT) cell of edge
        /// <em>-e - 1</em>.
        /// </summary>
        /// <remarks>
        /// Note: the second index does NOT correlate with the face index
        /// of the cell. This is because, in the case of hanging nodes, the
        /// number of edges that bound to one cell is not equal to the
        /// number of faces, i.e., more than one edge is associated with
        /// one face.
        /// </remarks>
        int[][] Cells2Edges {
            get;
        }


        /// <summary>
        /// True, if _all_ geometrical cells which make up this logical cell are affine-linear. 
        /// </summary>
        /// <param name="jCell">local cell index.</param>
        /// <returns>reference element index.</returns>
        bool IsCellAffineLinear(int jCell);

        /// <summary>
        /// Returns the volume (to be more exact: the
        /// <see cref="SpatialDimension"/> - dimensional measure) of the
        /// cell <paramref name="j"/>;
        /// </summary>
        /// <param name="j">local cell index</param>
        /// <returns></returns>
        double GetCellVolume(int j);

        /// <summary>
        /// Computes the bounding box of cell <paramref name="j"/>.
        /// </summary>
        /// <param name="j">local cell index.</param>
        /// <param name="bb">
        /// on exit, the bounding box of cell j.
        /// </param>
        void GetCellBoundingBox(int j, BoundingBox bb);

        
        /// <summary>
        /// polynomial interpolation degree of the Reference-to-Global coordinate transformation (maximum over all geometrical parts).
        /// </summary>
        /// <param name="jCell"></param>
        /// <returns></returns>
        int GetInterpolationDegree(int jCell);


        /// <summary>
        /// Center-of-gravity for the cell
        /// </summary>
        double[] GetCenter(int jCell);
    }

    public interface IVertexData {


        MultidimensionalArray Coordinates {
            get;
        }


        /// <summary>
        /// For each vertex, the local indices of the adjacent cells;
        ///  - 1st index: local vertex index
        ///  - 2nd index: collection
        /// </summary>
        int[][] VerticeToCell {
            get;
        }
    }



    public interface IGeometricalEdgeData {

        /// <summary>
        /// total number of all edges handled on this processor;
        /// </summary>
        int Count {
            get;
        }

        /// <summary>
        /// Transformations from edge coordinate system to local cell
        /// coordinate systems.
        /// </summary>
        IList<AffineTrafo> Edge2CellTrafos {
            get;
        }

        /// <summary>
        /// Square-root of the Gramian determinat for each transformation in <see cref="Edge2CellTrafos"/>, i.e.
        /// if \f$ \myMatrix{M} \f$ 
        /// is the matrix of the transformation, this number is 
        /// \f$ \sqrt{ \operatorname{det} ( \myMatrix{M}^T \cdot \myMatrix{M} ) } \f$.
        /// </summary>
        MultidimensionalArray Edge2CellTrafos_SqrtGramian {
            get;
        }

        /// <summary>
        /// For each edge-to-cell transformation, see <see cref="Edge2CellTrafos"/>,
        /// the index of the cell reference element, i.e. an index into <see cref="EdgeRefElements"/>.
        /// </summary>
        IList<int> Edge2CellTrafosRefElementIndices {
            get;
        }

        /// <summary>
        /// Edge-to-Cell - transformation index, i.e. index into <see cref="Edge2CellTrafos"/>;
        /// - 1st index: local edge index; resp. part index (for aggregate grids).
        /// - 2nd index: 0,1 first and second neighbor;
        /// </summary>
        int[,] Edge2CellTrafoIndex {
            get;
        }


        /// <summary>
        /// Edge Tags, index represents local edge index.
        /// </summary>
        byte[] EdgeTags {
            get;
        }

        /// <summary>
        /// Additional edge information, index represents local edge index.
        /// </summary>
        EdgeInfo[] Info {
            get;
        }




        /// <summary>
        /// Face index, where the numbering of faces is defined by the reference element, see e.g. <see cref="RefElement.FaceToVertexIndices"/>.
        /// - 1st index: local edge index; resp. part index (for aggregate grids).
        /// - 2nd index: 0 and 1 for first and second neighbor;
        /// </summary>
        byte[,] FaceIndices {
            get;
        }

        /// <summary>
        /// All reference elements for edges.
        /// </summary>
        RefElement[] EdgeRefElements {
            get;
        }

        /// <summary>
        /// For edge number <paramref name="e"/>, the index into
        /// <see cref="EdgeRefElements"/>.
        /// </summary>
        int GetRefElementIndex(int e);

        /// <summary>
        /// Some hack, used by <see cref="NodeSet.GetVolumeNodeSet"/>; 
        /// only effective (un-equal 0), if more than one grid is used in the application.
        /// </summary>
        int e2C_offet {
            get;
        }

        /// <summary>
        /// Local *geometric* cell indices of cells that belong to an edge.
        /// - 1st index: local edge index
        /// - 2nd index: 0,1 first and second neighbor;
        /// </summary>
        /// <remarks>
        /// Example: Let be <see cref="CellIndices"/>[i,0] = 123 and
        /// <see cref="CellIndices"/>[i,1] = 321; Then edge i is located
        /// on the intersection of (the closure of) cell 123 and cell 321;
        /// A negative cell index indicates that an edge is only subset of
        /// one cell (cells on the border of an domain). The negative cell
        /// index is always stored at the 2nd entry.
        /// </remarks>
        int[,] CellIndices {
            get;
        }

        /// <summary>
        /// Local *logical* cell indices of cells that belong to an edge.
        /// - 1st index: local edge index
        /// - 2nd index: 0,1 first and second neighbor;
        /// </summary>
        int[,] LogicalCellIndices {
            get;
        }

        /// <summary>
        /// For each edge that is affine-linear (i.e. <see cref="Info"/>[e]
        /// &amp; <see cref="EdgeInfo.EdgeIsAffineLinear"/> != 0), the
        /// square root of the Gram determinant; NaN for nonlinear edges.<br/>
        /// 1st index: local edge index;
        /// </summary>
        /// <remarks>
        /// Let be 
        /// \f[ 
        ///   \mathbb{R}^{D-1} 
        ///     \ni \vec{\xi} 
        ///       \mapsto
        ///         \gamma(\vec{\xi}) \in
        ///           \mathbb{R}^{D-1}
        /// \f]
        /// the mapping from the edge coordinate system to the physical coordinate system.
        /// Then the integral of \f$  f \f$  over the edge 
        /// \f$ \gamma(K_\textrm{ref}) \f$ 
        /// is given as 
        /// \f[ 
        ///   \int_{\vec{x} \in \gamma(K_\textrm{ref})} f(\vec{x}) \ \textrm{dS}
        ///   =
        ///   \int_{\xi \in K_\textrm{ref}} f(\gamma(\xi)) g(\vec{\xi}) \ \textrm{d} \vec{\xi}
        /// \f]
        /// with the square-root of the Gram determinant
        /// \f[ 
        ///   g(\vec{xi}) = \sqrt{ 
        ///      \textrm{det} ( (\partial \gamma)^T \cdot (\partial \gamma) )  
        ///   }.
        /// \f]
        /// If the transformation 
        /// \f$ \gamma \f$
        /// of the edge to the global coordinate system 
        /// is affine-linear, the Jacobian 
        /// \f$ \partial \gamma \f$
        /// is constant and 
        /// \f$ g \f$
        /// can be precomputed.
        /// (see Analysis 2, Königsberger, Springer-Verlag 2000, pp. 343)
        /// </remarks>
        MultidimensionalArray SqrtGramian {
            get;
        }

        /// <summary>
        /// Cached normals at nodes.
        /// </summary>
        Caching.EdgeNormalsCacheLogic_CNsFace NormalsCache {
            get;
        }

        /// <summary>
        /// true, if edge <paramref name="e"/> is affine-linear, false if not.
        /// </summary>
        bool IsEdgeAffineLinear(int e);

        /// <summary>
        /// Computes the normals on face <paramref name="iFace"/> in the
        /// volume coordinate system of a given cell
        /// <paramref name="jCell"/> at the given <paramref name="Nodes"/>
        /// and writes the result into <paramref name="NormalsOut"/>
        /// </summary>
        /// <param name="jCell">Cell index</param>
        /// <param name="iFace">Face index</param>
        /// <param name="Nodes">
        /// Evaluation nodes
        /// <list type="bullet">
        ///   <item>1st index: Node index</item>
        ///   <item>2nd index: Spatial dimension</item>
        /// </list>
        /// </param>
        /// <param name="NormalsOut">
        /// <list type="bullet">
        ///     <item>1st index: cell index</item>
        ///     <item>2nd index: Node index</item>
        ///     <item>3rd index: Spatial dimension</item>
        /// </list>
        /// </param>
        /// <param name="QuadMetric">
        /// A by-product: the integral transformation metric.
        /// </param>
        /// <param name="Offset">
        /// An offset into the first entry of <paramref name="NormalsOut"/> and <paramref name="QuadMetric"/>.
        /// </param>
        void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut, MultidimensionalArray QuadMetric, int Offset);


        /// <summary>
        /// Computes the normals on face <paramref name="iFace"/> in the
        /// volume coordinate system of a given cell
        /// <paramref name="jCell"/> at the given <paramref name="Nodes"/>
        /// and writes the result into <paramref name="NormalsOut"/>
        /// </summary>
        /// <param name="jCell">Cell index</param>
        /// <param name="iFace">Face index</param>
        /// <param name="Nodes">
        /// Evaluation nodes
        /// <list type="bullet">
        ///   <item>1st index: Node index</item>
        ///   <item>2nd index: Spatial dimension</item>
        /// </list>
        /// </param>
        /// <param name="NormalsOut">
        /// <list type="bullet">
        ///     <item>1st index: cell index</item>
        ///     <item>2nd index: Node index</item>
        ///     <item>3rd index: Spatial dimension</item>
        /// </list>
        /// </param>
        void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut);


        /// <summary>
        /// For each (edge) reference element, this method provides a
        /// mask containing all cells which are mapped from the specific
        /// reference element.
        /// </summary>
        /// <param name="Kref">
        /// Reference element for edges.
        /// </param>
        EdgeMask GetEdges4RefElement(RefElement Kref);
     
        
        /// <summary>
        /// returns the area (to be more exact: the (D-1) - dimensional measure) of the geometrical edge <paramref name="e"/>
        /// </summary>
        /// <param name="e">logical edge index</param>
        /// <returns></returns>
        double GetEdgeArea(int e);
    }


    /// <summary>
    /// 
    /// </summary>
    public interface ILogicalEdgeData {

        /// <summary>
        /// total number of all edges handled on this processor;
        /// </summary>
        int Count {
            get;
        }

        /// <summary>
        /// Local cell indices of cells that belong to an edge.
        /// - 1st index: local edge index
        /// - 2nd index: 0,1 first and second neighbor;
        /// </summary>
        /// <remarks>
        /// Example: Let be <see cref="CellIndices"/>[i,0] = 123 and
        /// <see cref="CellIndices"/>[i,1] = 321; Then edge i is located
        /// on the intersection of (the closure of) cell 123 and cell 321;
        /// A negative cell index indicates that an edge is only subset of
        /// one cell (cells on the border of an domain). The negative cell
        /// index is always stored at the 2nd entry.
        /// </remarks>
        int[,] CellIndices {
            get;
        }



        /// <summary>
        /// Only used for aggregation grids, where each edge can be triangulated.
        /// - 1st index: local edge index
        /// - 2nd index: enumeration of parts.
        /// </summary>
        int[][] EdgeToParts {
            get;
        }


        /// <summary>
        /// returns the area (to be more exact: the (D-1) - dimensional measure) of the logical edge <paramref name="e"/>,
        /// which is the sum of all geometrical parts
        /// </summary>
        /// <param name="e">logical edge index</param>
        /// <returns></returns>
        double GetEdgeArea(int e);

    }



    public interface IParallelization {


        /// <summary>
        /// Conversion of global cell indices to local cell indices, 
        /// i.e. the inverse of <see cref="GlobalIndicesExternalCells"/>. 
        ///  - keys: global indices of external/ghost cells 
        ///  - values: local indices of external/ghost cells
        /// </summary>
        Dictionary<long, int> Global2LocalIdx {
            get;
        }

        /// <summary>
        /// Global indices of external cells (local indices j in the range
        /// <see cref="ICellData.NoOfLocalUpdatedCells"/> &lt;= j &lt;
        /// <see cref="ICellData.NoOfCells"/>); Note that there is an index
        /// offset, so the entry at index 0 is the global index of cell at
        /// local index <see cref="CellData.NoOfLocalUpdatedCells"/>;
        /// See also <see cref="GlobalIndicesExternalCells"/>;
        /// </summary>
        long[] GlobalIndicesExternalCells {
            get;
        }


        /// <summary>
        /// list of processes (MPI ranks) which receive data from this process
        /// - index: enumeration
        /// - content: MPI process rank
        /// </summary>
        int[] ProcessesToSendTo {
            get;
        }

        /// <summary>
        /// List of processes (MPI ranks) which send data to this process.
        /// </summary>
        int[] ProcessesToReceiveFrom {
            get;
        }

        /// <summary>
        /// Local cell indices (only border cells) which must be send to
        /// other processors; for each processor, the communication list is
        /// stored in ascending order.
        /// </summary>
        /// <remarks>
        /// - 1st index: target processor; if the 'p'-th entry is null, there is no communication with processor 'p'
        /// - 2nd index: enumeration;
        /// - content: local cell index
        /// </remarks>
        int[][] SendCommLists {
            get;
        }

        /// <summary>
        /// For each process, a local cell index at which items received
        /// from other processes should be inserted;
        /// - index: MPI process rank of process from which data is received.
        /// - content: local cell index.
        /// </summary>
        int[] RcvCommListsInsertIndex {
            get;
        }

        /// <summary>
        /// For each process, the number of cells that are received from
        /// this process
        /// - index: MPI process rank of process from which data is received
        /// - content: number of cells
        /// </summary>
        int[] RcvCommListsNoOfItems {
            get;
        }
    }
}
