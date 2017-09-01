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
using System.Linq;
using System.Text;
using ilPSP;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation {

    /// <summary>
    /// Extends the functionality of the <see cref="MultidimensionalArray"/>-class to support caching.
    /// </summary>
    public class NodeSet : MultidimensionalArray {

        static object syncRoot = new Object();
        static int RefCounter = 123;

        /// <summary>
        /// Constructor: initializes this node set as a (non-shallow) clone of the array <paramref name="nds"/>.
        /// </summary>
        public NodeSet(RefElement r, MultidimensionalArray nds)
            : base(2) //
        {
        
            if(nds.Dimension != 2)
                throw new ArgumentException("Expecting 2D-array. 1st dim.: node index, 2nd dim.: spatial direction.");
            if(nds.GetLength(1) > 3)
                throw new ArgumentException("Spatial dimension is expected to be lower or equal to 3.");
            base.Allocate(nds.Lengths);
            base.Set(nds);
            base.LockForever();
            this.RefElement = r;
            lock(syncRoot) {
                this.Reference = RefCounter;
                if(RefCounter >= int.MaxValue)
                    throw new ApplicationException("NodeSet ref-counter overflow.");
                RefCounter++;
            }
        }

        /// <summary>
        /// Constructor: initializes this node set as a (non-shallow) clone of the array <paramref name="nds"/>.
        /// </summary>
        public NodeSet(RefElement r, double[,] nds)
            : base(2) //
        {

            if(nds.GetLength(1) > 3)
                throw new ArgumentException("Spatial dimension is expected to be lower or equal to 3.");
            base.Allocate(nds.GetLength(0), nds.GetLength(1));
            base.Set2DArray(nds);
            base.LockForever();
            this.RefElement = r;
            lock(syncRoot) {
                this.Reference = RefCounter;
                if(RefCounter >= int.MaxValue)
                    throw new ApplicationException("NodeSet ref-counter overflow.");
                RefCounter++;
            }
        }

        /// <summary>
        /// Constructor: initializes this node set
        /// containing only one point <paramref name="point"/>.
        /// </summary>
        public NodeSet(RefElement r, double[] point)
            : base(2) //
        {
            if(point.Length > 3)
                throw new ArgumentException("Spatial dimension is expected to be lower or equal to 3.");
            base.Allocate(1, point.Length);
            base.ExtractSubArrayShallow(0, -1).SetVector(point);
            base.LockForever();
            this.RefElement = r;
            lock(syncRoot) {
                this.Reference = RefCounter;
                if(RefCounter >= int.MaxValue)
                    throw new ApplicationException("NodeSet ref-counter overflow.");
                RefCounter++;
            }
        }


        /// <summary>
        /// Creates a node set with all nodes set to 0 -- these values can still be changed.
        /// Use this constructor with care!
        /// Note that this node set is in an invalid state until <see cref="MultidimensionalArray.LockForever"/> is called.
        /// </summary>
        public NodeSet(RefElement r, int NoOfNodes, int D)
            : base(2) //
        {
            base.Allocate(NoOfNodes, D);
            this.RefElement = r;
            lock(syncRoot) {
                this.Reference = RefCounter;
                if(RefCounter >= int.MaxValue)
                    throw new ApplicationException("NodeSet ref-counter overflow.");
                RefCounter++;
            }
        }


        /// <summary>
        /// unique refenence number for the node set.
        /// </summary>
        public int Reference {
            get;
            private set;
        }
        
        /// <summary>
        /// Clones this node set.
        /// </summary>
        new public NodeSet CloneAs() {
            if(!base.IsLocked)
                throw new NotSupportedException("NodeSet must be locked before first usage.");

            var R = new NodeSet(this.RefElement, this.GetLength(0), this.GetLength(1));
            R.Set(this);
            R.LockForever();
            return R;
        }

        /// <summary>
        /// The reference element for which this node set is valid; since cached values,
        /// like basis polynomials at nodes depend on the reference element,
        /// the node set needs to be associated with the reference element.
        /// </summary>
        public Grid.RefElements.RefElement RefElement {
            get;
            private set;
        }
        
        /// <summary>
        /// The index of <see cref="RefElement"/> within the 
        /// reference elements for cells in grid <paramref name="gridData"/>
        /// </summary>
        public int GetVolumeRefElementIndex(IGridData gridData) {
            if(!base.IsLocked)
                throw new NotSupportedException("NodeSet must be locked before first usage.");

            return gridData.iGeomCells.RefElements.IndexOf(this.RefElement, (A, B) => object.ReferenceEquals(A, B));
        }

        /// <summary>
        /// The index of <see cref="RefElement"/> within the 
        /// reference elements for edges in grid <paramref name="gridData"/>
        /// </summary>
        public int GetEdgeRefElementIndex(IGridData gridData) {
            if(!base.IsLocked)
                throw new NotSupportedException("NodeSet must be locked before first usage.");

            return gridData.iGeomEdges.EdgeRefElements.IndexOf(this.RefElement, (A, B) => object.ReferenceEquals(A, B));
        }

        /// <summary>
        /// Determines, for a given grid <paramref name="g"/>, whether 
        /// this node set describes local cell-coordinates or edge-coordinates.
        /// </summary>
        public NodeCoordinateSystem GetNodeCoordinateSystem(IGridData g) {
            if(!base.IsLocked)
                throw new NotSupportedException("NodeSet must be locked before first usage.");
            
            if(GetVolumeRefElementIndex(g) >= 0)
                return NodeCoordinateSystem.CellCoord;

            if(GetEdgeRefElementIndex(g) >= 0)
                return NodeCoordinateSystem.EdgeCoord;

            throw new ArgumentException("Unable to identify coordinate (either cell- or edge-reference) system for NodeSet.");
        }
                 
        /// <summary>
        /// Number of nodes; equal to the 0-th dimension of this array.
        /// </summary>
        public int NoOfNodes {
            get {
                if(!base.IsLocked)
                    throw new NotSupportedException("Node set must be locked before first usage.");
                Debug.Assert((base.GetLength(1) == this.RefElement.SpatialDimension) || (base.GetLength(1) == 1 && this.RefElement.SpatialDimension == 0), "Mismatch between number of spatial directions in node set and reference element.");
                return base.GetLength(0);
            }
        }

        /// <summary>
        /// Vector length of the nodes in this node set; equal to the 1-st dimension of this array.
        /// </summary>
        public int SpatialDimension {
            get {
                if(!base.IsLocked)
                    throw new NotSupportedException("Node set must be locked before first usage.");
                Debug.Assert((base.GetLength(1) == this.RefElement.SpatialDimension) || (base.GetLength(1) == 1 && this.RefElement.SpatialDimension == 0), "Mismatch between number of spatial directions in node set and reference element.");
                return base.GetLength(1);
            }
        }
        
        NodeSet[] VolumeNodeSets;

        /// <summary>
        /// If this is a node set defined in the edge coordinate system, this method provides 
        /// the nodes transformed to the cell coordinate system.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="Edge2CellTrafoIndex">
        /// The transformation index (from edge to cell coordinate system), i.e. an index into <see cref="GridData.EdgeData.Edge2CellTrafos"/>,
        /// see also <see cref="GridData.EdgeData.Edge2CellTrafoIndex"/>.
        /// </param>
        /// <returns></returns>
        public NodeSet GetVolumeNodeSet(IGridData g, int Edge2CellTrafoIndex) {
            if(!base.IsLocked)
                throw new NotSupportedException("NodeSet must be locked before first usage.");
            Debug.Assert((base.GetLength(1) == this.RefElement.SpatialDimension) || (base.GetLength(1) == 1 && this.RefElement.SpatialDimension == 0), "Mismatch between number of spatial directions in node set and reference element.");
                
#if DEBUG
            if(this.GetNodeCoordinateSystem(g) != NodeCoordinateSystem.EdgeCoord) {
                throw new NotSupportedException("Operation only supported for edge node sets.");
            }
#endif            
            int D = g.SpatialDimension;
            int NN = this.NoOfNodes;
            int idx = Edge2CellTrafoIndex + g.iGeomEdges.e2C_offet;

            if(VolumeNodeSets == null) {
                // alloc mem, if necessary
                // ++++++++++++++++++++++++
                VolumeNodeSets = new NodeSet[g.iGeomEdges.Edge2CellTrafos.Count + g.iGeomEdges.e2C_offet];
            }
            if(VolumeNodeSets.Length < (g.iGeomEdges.Edge2CellTrafos.Count + g.iGeomEdges.e2C_offet)) {
                // re-alloc mem, if necessary
                // ++++++++++++++++++++++++++
                var newVNS = new NodeSet[g.iGeomEdges.Edge2CellTrafos.Count + g.iGeomEdges.e2C_offet];
                Array.Copy(VolumeNodeSets, 0, newVNS, 0, VolumeNodeSets.Length);
                VolumeNodeSets = newVNS;
            }
            if(VolumeNodeSets[idx] == null) {
                // transform edge-nodes to cell-nodes, if necessary
                // ++++++++++++++++++++++++++++++++++++++++++++++++
                AffineTrafo Trafo = g.iGeomEdges.Edge2CellTrafos[Edge2CellTrafoIndex];
                int iKref = g.iGeomEdges.Edge2CellTrafosRefElementIndices[Edge2CellTrafoIndex];

                NodeSet volNS = new NodeSet(g.iGeomCells.RefElements[iKref], NN, D);
                Trafo.Transform(this, volNS);
                volNS.LockForever();

                VolumeNodeSets[idx] = volNS;
            }

            return VolumeNodeSets[idx];
        }
    }
}

