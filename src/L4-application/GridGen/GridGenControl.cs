using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.GridGen {

    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    [DataContract]
    public class GridGenControl : AppControl {


        /// <summary>
        /// Base-class for mesh blocks;
        /// </summary>
        [Serializable]
        [DataContract]
        abstract public class MeshBlock {

            /// <summary>
            /// override to create a specific mesh, e.g. a 3D Cartesian
            /// </summary>
            public abstract GridCommons CreateGrid();

        }

        /// <summary>
        /// a 3D Hex mesh block
        /// </summary>
        [Serializable]
        [DataContract]
        public class Cartesian3D : MeshBlock {

            /// <summary>
            /// Nodes of Hex-Mesh in x-Direction
            /// </summary>
            [DataMember]
            public double[] xNodes;
            
            /// <summary>
            /// Nodes of Hex-Mesh in y-Direction
            /// </summary>
            [DataMember]
            public double[] yNodes;

            /// <summary>
            /// Nodes of HGex-Mesh in z-Direction
            /// </summary>
            [DataMember]
            public double[] zNodes;

            /// <summary>
            /// 
            /// </summary>
            public CellType _CellType = CellType.Cube_Linear;

            /// <summary>
            /// 
            /// </summary>
            public bool periodicX = false;

            /// <summary>
            /// 
            /// </summary>
            public bool periodicY = false;

            /// <summary>
            /// 
            /// </summary>
            public bool periodicZ = false;

            /// <summary>
            /// 
            /// </summary>
            public BoundingBox[] CutOuts;


            /// <summary>
            /// instantiation of the mesh block
            /// </summary>
            public override GridCommons CreateGrid() {
                return Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes, _CellType, periodicX, periodicY, periodicZ, CutOuts);
            }
        }

        /// <summary>
        /// Collection of mesh blocks; if more than one block is present, the blocks patched together.
        /// </summary>
        public MeshBlock[] GridBlocks;


        /// <summary>
        /// Assignment of boundary conditions;
        /// These are pairs of bounding boxes and edge tag names; the latter ones determine the boundary conditions for the solver to use.
        /// Later boundary conditions 
        /// </summary>
        public List<(BoundingBox Region, string EdgeTagName)> BoundaryRegions = new List<(BoundingBox Region, string EdgeTagName)> ();


        /// <summary>
        /// Edge Tag names which are in any case added to the mesh.
        /// </summary>
        public string[] EdgeTagNamesToEnsure;

    }
}
