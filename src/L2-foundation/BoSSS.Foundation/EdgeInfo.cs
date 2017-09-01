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
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {
    
    /// <summary>
    /// Additional information
    /// </summary>
    [Flags]
    public enum EdgeInfo {

        /// <summary>
        /// edge with nothing special
        /// </summary>
        Default = 0,

        /// <summary>
        /// masks the region where the reference element index for the
        /// edge is stored.
        /// </summary>
        EdgeSimplexIdxMask = 0x3, // binary: 0011, can store numbers 0 to 3;

        /// <summary>
        /// The edge is aligned to cell 1 (in-cell) in a non-conformal manner (hanging nodes).
        /// If this is _not_ set, the edge covers an entire face of the corresponding cell.
        /// </summary>
        Cell1_Nonconformal = 4, // binary: 0100

        /// <summary>
        /// the edge is aligned to cell 2 (out-cell) in a non-conformal manner (hanging nodes).
        /// If this is _not_ set, the edge covers an entire face of the corresponding cell.
        /// </summary>
        Cell2_Nonconformal = 8, // binary: 1000

        /// <summary>
        /// the transformation from the edge to the global coordinate
        /// system is affine-linear. In this case, also the normal
        /// vector is constant along the edge.
        /// </summary>
        EdgeIsAffineLinear = 0x10, // binary: 0001'0000

        /// <summary>
        /// edge is located on the boundary of the computational domain
        /// </summary>
        Boundary = 0x20, // binary: 0010'0000

        /// <summary>
        /// edge is locates on the border between two MPI processes
        /// </summary>
        Interprocess = 0x40, // binary: 0100'0000
        
        /// <summary>
        /// marks cells which are aggregated from smaller parts
        /// </summary>
        IsAggregate = 0x80, // binary: 1000'0000
        
        /// <summary>
        /// All flags on, implementing this ensures that the enum 
        /// is compiled using 32 bits.
        /// </summary>
        AllOn = unchecked((int)~0)
    }

}
