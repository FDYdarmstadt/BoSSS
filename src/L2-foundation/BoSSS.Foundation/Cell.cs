﻿/* =======================================================================
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
using System.Diagnostics;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// Base class for elements for the physical domain.
    /// </summary>
    [Serializable]
    abstract public class Element {

        /// <summary>
        /// Global identification number of the element,
        /// which remains constant under re-ordering.
        /// </summary>
        public long GlobalID;

        /// <summary>
        /// Parameters for the transformation from reference to physical domain.
        /// </summary>
        /// <remarks>
        /// <list type="bullet">
        ///    <item>
        ///      1st index: interpolation node in physical coordinates
        ///    </item>        
        ///    <item>
        ///      2nd index: spatial direction (0 is x, 1 is y, ...)
        ///    </item>
        /// </list>
        /// </remarks>
        public MultidimensionalArray TransformationParams;

        /// <summary>
        /// Global Node indices: usually, the connectivity information of cells
        /// is derived from this information. Alternatively, cell connectivity
        /// may be specified by <see cref="Cell.CellFaceTags"/>.
        /// </summary>
        public long[] NodeIndices;

        /// <summary>
        /// Element/Cell type;
        /// </summary>
        public CellType Type;

        /// <summary>
        /// 
        /// </summary>
        public override string ToString() {
            using(var stw = new System.IO.StringWriter()) {
                stw.WriteLine("Cell " + GlobalID + " (" + Type + "): ");
                if(NodeIndices != null) {
                    stw.Write("Node indices: ");
                    foreach (long ni in NodeIndices)
                        stw.Write(ni + " ");
                } else {
                    stw.Write("No Node Indices defined. ");
                }
                stw.Write("Transformation parameters: ");
                if (TransformationParams != null) {
                    stw.WriteLine();
                    IMatrixExtensions.SaveToStream(this.TransformationParams, stw);
                } else {
                    stw.WriteLine("NULL.");
                }
                return stw.ToString();
            }
        }

    }

    /// <summary>
    /// Represents one cell in physical domain
    /// </summary>
    [Serializable]
    public class Cell : Element, ICloneable {

        /// <summary>
        /// Alternative method to describe cell neighborship, can be used
        /// instead of <see cref="Element.NodeIndices"/>.
        /// Mandatory
        /// to specify periodic boundary conditions and connection between
        /// cells with hanging nodes. Should be null, or of length 0 if not
        /// used.
        /// </summary>
        public CellFaceTag[] CellFaceTags = null;

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var ret = new Cell();
            ret.GlobalID = this.GlobalID;
            ret.TransformationParams = this.TransformationParams.CloneAs();
            ret.NodeIndices = this.NodeIndices.CloneAs();
            ret.CellFaceTags = this.CellFaceTags == null ? null : this.CellFaceTags.CloneAs();
            Debug.Assert(typeof(CellFaceTag).IsValueType);
            ret.Type = this.Type;
            ret.CoarseningClusterID = CoarseningClusterID;
            ret.ParentCell = ParentCell != null ? ParentCell.CloneAs() : null;
            ret.RefinementLevel = this.RefinementLevel;
            ret.CoarseningClusterSize = this.CoarseningClusterSize;
            ret.CoarseningLeafIndex = this.CoarseningLeafIndex;
            return ret;
        }
        
                
        /// <summary>
        /// If this cell was created by adaptive mesh refinement, 
        /// the parent cell from which this cell was created from. This member is required to be able to coarsen the mesh again.
        /// </summary>
        public Cell ParentCell;

        /// <summary>
        /// How often the cell was refined.
        /// </summary>
        public int RefinementLevel;

        /// <summary>
        /// If this cell was created by adaptive mesh refinement, 
        /// the siblings in the coarsening cluster can be identified by this ID/Token.
        /// </summary>
        public int CoarseningClusterID;
        
        /// <summary>
        /// If this cell was created by adaptive mesh refinement, 
        /// the number of siblings in the coarsening cluster.
        /// </summary>
        public int CoarseningClusterSize;

        /// <summary>
        /// If this cell was created by adaptive mesh refinement, 
        /// the index within the coarsening cluster.
        /// </summary>
        public int CoarseningLeafIndex;

        /// <summary>
        /// 
        /// </summary>
        public override string ToString() {
            using (var stw = new System.IO.StringWriter()) {
                stw.Write(base.ToString());
                stw.WriteLine("Parent cell: " + (ParentCell != null));
                stw.WriteLine("Refinement Level: " + RefinementLevel);
                stw.WriteLine("CoarseningClusterID: " + CoarseningClusterID);
                stw.WriteLine("CoarseningClusterSize: " + CoarseningClusterSize);
                stw.WriteLine("CoarseningLeafIndex: " + CoarseningLeafIndex);

                if(CellFaceTags != null ) {
                    foreach (var e in CellFaceTags)
                        stw.WriteLine("Face Tag: " + e.ToString());
                } else {
                    stw.WriteLine("no cell face tags");
                }

                return stw.ToString();
            }
        }

    }

    /// <summary>
    /// Represents a boundary element of the physical domain (i.e., an element
    /// where boundary conditions are applied).
    /// </summary>
    [Serializable]
    public class BCElement : Element, ICloneable {

        /// <summary>
        /// The tag of this element
        /// </summary>
        public byte EdgeTag;

        /// <summary>
        /// Global ID's of neighboring cells.
        /// </summary>
        public long[] NeighCell_GlobalIDs;

        /// <summary>
        /// Bool indicating whether this a conformal element.
        /// </summary>
        public bool Conformal;

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var ret = new BCElement();
            ret.GlobalID = this.GlobalID;
            ret.TransformationParams = this.TransformationParams.CloneAs();
            ret.NodeIndices = this.NodeIndices.CloneAs();
            ret.EdgeTag = this.EdgeTag;
            ret.Type = this.Type;
            return ret;
        }
    }
}
