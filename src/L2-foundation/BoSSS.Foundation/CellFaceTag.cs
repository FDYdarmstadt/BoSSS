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
using System.Diagnostics;
using System.Runtime.Serialization;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// carries information which is 'tagged' onto the faces of cells
    /// </summary>
    [Serializable]
    [DataContract]
    public struct CellFaceTag {

        /// <summary>
        /// Tags (i.e. numbers) which define information for edge (connection
        /// between cells);
        /// </summary>
        /// <remarks>
        /// In BoSSS, different regions of the computational boundary (e.g.
        /// inlet, wall, moving wall,...) are identified by numbers (edge tags)
        /// which are assigned to the grid edges on the boundary. These numbers
        /// are small, in the region of 1 (including) to
        /// <see cref="GridCommons.FIRST_PERIODIC_BC_TAG"/>;
        /// in total, there are three different ranges of edge tags:
        /// <list type="bullet">
        ///    <item>
        ///      _Internal edges_:
        ///      A value 0 indicates an internal edge (this information is
        ///      redundant, however);
        ///    </item>
        ///    <item>
        ///      _Boundaries_:
        ///      Values between 1 (including) and
        ///      <see cref="GridCommons.FIRST_PERIODIC_BC_TAG"/> (excluding)
        ///      are reserved for user defined boundary conditions (by using
        ///      border-edge-flux functions);
        ///    </item>
        ///    <item>
        ///      _Periodic boundaries_:
        ///      Values between <see cref="GridCommons.FIRST_PERIODIC_BC_TAG"/>
        ///      (including) and 255 (including) are reserved for periodic
        ///      boundary conditions (which are, in contrast to user defined
        ///      boundary conditions, implemented on a lower level by
        ///      <see cref="GridCommons.GetCellNeighbourship"/> and the
        ///      <see cref="GridCommons.PeriodicTrafo"/>-transformations. So,
        ///      e.g. a tag of
        ///      <see cref="GridCommons.FIRST_PERIODIC_BC_TAG"/> + 3 indicates
        ///      a periodic boundary condition with transformation
        ///      <see cref="GridCommons.PeriodicTrafo"/>[3].
        ///     </item>
        /// </list>
        /// This entry may be null if the cell has only internal edges, i.e.
        /// all edge tags are 0.
        /// </remarks>
        public byte EdgeTag {
            set {
                if(value >= 255)
                    throw new ArgumentOutOfRangeException("255 is a reserved value");
                //bool signFlag = PeriodicInverse;
                //m_SignedEdgeTag = value;
                //PeriodicInverse = signFlag;

                if (m_SignedEdgeTag < 0)
                    throw new NotSupportedException();

                m_SignedEdgeTag = (value & 0xFF) | (m_SignedEdgeTag & ~0xFF);
            }
            get {
                int r;
                if (m_SignedEdgeTag < 0) {
                    r = Math.Abs(m_SignedEdgeTag);
                } else {
                    r = m_SignedEdgeTag & 0xFF;
                }
                Debug.Assert(r < 255, "in reserved value range");
                return (byte)r;
            }
        }

        const int PeriodicInverseMask = 0x4000000;
        const int EdgeMayBeEmptyMask = 0x2000000;

        /// <summary>
        /// Indicates on which 'side' of the periodic transformation this cell
        /// resides. If false, resides on the side associated with the original
        /// periodic transformation. If true, resides on the side associated
        /// with the inverse of this transformation.
        /// </summary>
        public bool PeriodicInverse {
            get {
                return (m_SignedEdgeTag & PeriodicInverseMask) != 0
                    || m_SignedEdgeTag < 0; // legacy stuff
            }
            set {
                if (m_SignedEdgeTag < 0)
                    throw new NotSupportedException();

                if (value)
                    m_SignedEdgeTag |= PeriodicInverseMask;
                else
                    m_SignedEdgeTag &= ~PeriodicInverseMask;
            }
        }

        /// <summary>
        /// Used for grid refinement operations; this flag indicates that an edge may be empty, and it will be tested geometrically.
        /// If it is empty, it will be ignored.
        /// </summary>
        public bool EdgeMayBeEmpty {
            get {
                if (m_SignedEdgeTag < 0)
                    return false;

                return (m_SignedEdgeTag & EdgeMayBeEmptyMask) != 0;
            }
            set {
                if (m_SignedEdgeTag < 0)
                    throw new NotSupportedException();

                if (value)
                    m_SignedEdgeTag |= EdgeMayBeEmptyMask;
                else
                    m_SignedEdgeTag &= ~EdgeMayBeEmptyMask;
            }
        }


        [DataMember]
        int m_SignedEdgeTag;

        /// <summary>
        /// Face index of the owner cell (*not* the cell referenced with <see cref="NeighCell_GlobalID"/>).
        /// </summary>
        [DataMember]
        public int FaceIndex;

        /// <summary>
        /// Two options, depending on  <see cref="EdgeTag"/>:
        /// - internal edges or periodic boundaries: the global identification of the neighbor cell, i.e. this entry is non-negative.
        /// - not used for cell face tags that represent boundary conditions: then, this entry should be negative.
        /// </summary>
        [DataMember]
        public long NeighCell_GlobalID;

        /// <summary>
        /// True, if the face vertices match 1-to-1 with the face vertices of
        /// some other cell. This information is optional, i.e. if this flag is
        /// set to false, the code will test if the connection to the neighbor
        /// is conformal. However, this causes some performance overhead and
        /// the result depends on certain tolerance thresholds. 
        /// </summary>
        [DataMember]
        public bool ConformalNeighborship;

       
    }
}
