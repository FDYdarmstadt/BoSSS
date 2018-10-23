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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Control;

namespace CNS.Boundary {

    /// <summary>
    /// Delegate for all functions that can be evaluated to a boundary value.
    /// </summary>
    /// <param name="x">The global coordinate</param>
    /// <param name="time">The simulation time</param>
    /// <returns>The value at the boundary</returns>
    public delegate double BoundaryValue(double[] x, double time);

    /// <summary>
    /// Interface common to all implementations of a boundary condition map.
    /// Its main purpose is to provide a mapping between edge tags and the
    /// evaluation results of <see cref="BoundaryValue"/>s associated with this
    /// edge tag.
    /// </summary>
    public interface IBoundaryConditionMap {
        
        /// <summary>
        /// Mapping of edge tags to edge tag names
        /// </summary>
        IDictionary<byte, string> EdgeTagNames {
            get;
        }

        /// <summary>
        /// Retrieves the configured boundary condition for a given
        /// <paramref name="edgeTagName"/>.
        /// </summary>
        /// <param name="edgeTagName">The name of the edge</param>
        /// <returns>
        /// The boundary condition that has been assigned to edges with
        /// name <paramref name="edgeTagName"/>
        /// </returns>
        BoundaryCondition GetBoundaryCondition(string edgeTagName);

        /// <summary>
        /// Returns the value at the boundary value defined by
        /// <paramref name="EdgeTag"/> at global the coordinates
        /// <paramref name="x"/> at time <paramref name="time"/> for the inner
        /// values <paramref name="stateVector"/>
        /// </summary>
        /// <param name="EdgeTag">The edge tag of the boundary</param>
        /// <param name="time">The physical time</param>
        /// <param name="x">
        /// The global coordinates of the point at the boundary
        /// </param>
        /// <param name="normal">The outward unit normal vector</param>
        /// <param name="stateVector">The flow state inside the domain</param>
        /// <returns>The value at the boundary</returns>
        StateVector GetBoundaryState(byte EdgeTag, double time, double[] x, double[] normal, StateVector stateVector);

        /// <summary>
        /// Returns the <see cref="BoundaryCondition"/> associated to the given
        /// edge tag name
        /// </summary>
        /// <param name="edgeTageName">
        /// The name of the edge type
        /// </param>
        /// <returns>
        /// The respective boundary condition object.
        /// </returns>
        BoundaryCondition this[string edgeTageName] {
            get;
        }

        /// <summary>
        /// Returns the <see cref="BoundaryCondition"/> associated to the given
        /// edge tag
        /// </summary>
        /// <returns>
        /// The respective boundary condition object.
        /// </returns>
        /// <remarks>
        /// For use inside fluxes where only the EdgeTag is known
        /// </remarks>
        BoundaryCondition this[byte EdgeTag] {
            get;
        }
    }
}
