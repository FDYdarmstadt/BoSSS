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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Standard implementation of <see cref="IBoundaryConditionMap"/> which
    /// works with instances of <see cref="BoundaryConditionMap"/>.
    /// </summary>
    public class BoundaryConditionMap : BoSSS.Solution.Utils.BoundaryCondMap<CompressibleBcType>, IBoundaryConditionMap {

        private IGridData gridData;

        /// <summary>
        /// Cache for <see cref="ConditionMap"/>
        /// </summary>
        private Dictionary<string, BoundaryCondition> conditionMap;

        /// <summary>
        /// Mapping between edge tag names and associated boundary conditions
        /// </summary>
        protected Dictionary<string, BoundaryCondition> ConditionMap {
            get {
                if (conditionMap == null) {
                    conditionMap = new Dictionary<string, BoundaryCondition>();
                    foreach (byte edgeTag in gridData.iGeomEdges.EdgeTags) {
                        // Only process non-periodic boundary edges
                        if (edgeTag == 0 || edgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                            continue;
                        }

                        string edgeTagName = EdgeTagNames[edgeTag];
                        conditionMap[edgeTagName] = GetBoundaryCondition(edgeTagName);
                    }
                }

                return conditionMap;
            }
        }
        /*
        /// <summary>
        /// Configuration options
        /// </summary>
        protected readonly CNSControl control;
        */


        /// <summary>
        /// Constructs a new map by searching through all the edge tags
        /// (<see cref="GridData.EdgeData.EdgeTags"/> and instantiating sub classes
        /// of <see cref="BoundaryConditionMap"/> specific for the compressible
        /// Navier-Stokes solver depending on their edge tag names.
        /// </summary>
        /// <param name="gridData">The omnipresent grid data</param>
        /// <param name="control">Configuration options</param>
        public BoundaryConditionMap(IGridData gridData, AppControl control) : base(gridData, control.BoundaryValues) {
            this.gridData = gridData;
            this.control = control;
        }

        /// <summary>
        /// Searches the boundary config (<see cref="AppControl.BoundaryValues"/>)
        /// for a definition of a set of boundary values for an edge with name
        /// <paramref name="edgeTagName"/>.
        /// </summary>
        /// <param name="edgeTagName">
        /// The id of the boundary in question
        /// </param>
        /// <returns>The boundary values specified in the control file</returns>
        private AppControl.BoundaryValueCollection GetConfiguredBoundaryValues(string edgeTagName) {
            AppControl.BoundaryValueCollection boundaryValues = null;
            foreach (var tagConditionPair in control.BoundaryValues) {
                if (edgeTagName.StartsWith(tagConditionPair.Key, StringComparison.InvariantCultureIgnoreCase)) {
                    boundaryValues = tagConditionPair.Value;
                }
            }

            if (boundaryValues == null) {
                throw new ArgumentException(
                    "Unable to find a definition for boundary conditions for edges of type '" + edgeTagName + "'",
                    "edgeTagName");
            }

            return boundaryValues;
        }

        /// <summary>
        /// Finds the boundary evaluator (<see cref="AppControl.BoundaryValues"/>)
        /// associated with <paramref name="fieldName"/> for the condition
        /// <paramref name="boundaryValues"/>.
        /// </summary>
        /// <param name="boundaryValues">The boundary condition</param>
        /// <param name="fieldName">The id of the boundary in question</param>
        /// <returns>
        /// The evaluator stored in <see cref="AppControl.BoundaryValues"/>
        /// for a boundary with id <paramref name="fieldName"/>.
        /// </returns>
        private static Func<double[], double, double> GetBoundaryValueFunction(AppControl.BoundaryValueCollection boundaryValues, string fieldName) {
            Func<double[], double, double> boundaryValue = null;
            TryGetBoundaryValueFunction(boundaryValues, fieldName, out boundaryValue);

            if (boundaryValue == null) {
                throw new ArgumentException(
                    "Missing boundary specification for field '" + fieldName + "'",
                    "condition");
            }

            return boundaryValue;
        }

        /// <summary>
        /// Finds the boundary evaluator (<see cref="AppControl.BoundaryValues"/>)
        /// associated with <paramref name="fieldName"/> for the condition
        /// <paramref name="boundaryValues"/>.
        /// </summary>
        /// <param name="boundaryValues">
        /// The configured boundary values
        /// </param>
        /// <param name="fieldName">The id of the boundary in question</param>
        /// <param name="result">
        /// On exit: The evaluator stored in <see cref="AppControl.BoundaryValues"/>
        /// for a boundary with id <paramref name="fieldName"/> if it exists
        /// </param>
        /// <returns>
        /// True, if a boundary value has been found; false otherwise
        /// </returns>
        private static bool TryGetBoundaryValueFunction(AppControl.BoundaryValueCollection boundaryValues, string fieldName, out Func<double[], double, double> result) {
            return boundaryValues.Evaluators.TryGetValue(fieldName, out result);
        }

        /// <summary>
        /// Vector version of <see cref="GetBoundaryValueFunction"/> for the
        /// velocity components
        /// </summary>
        /// <param name="boundaryCondition">
        /// <see cref="GetBoundaryValueFunction"/>
        /// </param>
        /// <returns>
        /// A list of length
        /// <see cref="CNSEnvironment.NumberOfDimensions"/> of boundary
        /// values as returned by <see cref="GetBoundaryValueFunction"/>.
        /// </returns>
        private static Func<double[], double, double>[] GetVelocityBoundaryValueFunction(AppControl.BoundaryValueCollection boundaryCondition) {
            Func<double[], double, double>[] result = GetOptionalVelocityBoundaryValueFunction(boundaryCondition);
            
            if (result == null) {
                throw new Exception(String.Format(
                    "Missing definition of velocity component for boundary condition"));
            }

            return result;
        }

        /// <summary>
        /// Version of <see cref="GetVelocityBoundaryValueFunction"/> which
        /// supports optional vector fields.
        /// </summary>
        /// <param name="boundaryValues">
        /// <see cref="GetBoundaryValueFunction"/>
        /// </param>
        /// <returns>
        /// A list of length
        /// <see cref="CNSEnvironment.NumberOfDimensions"/> of boundary
        /// values as returned by <see cref="GetBoundaryValueFunction"/>.
        /// </returns>
        private static Func<double[], double, double>[] GetOptionalVelocityBoundaryValueFunction(AppControl.BoundaryValueCollection boundaryValues) {
            int numberOfDimensions = CNSEnvironment.NumberOfDimensions;

            // First check x-component only
            Func<double[], double, double> boundaryValue;
            bool found = TryGetBoundaryValueFunction(boundaryValues, "u0", out boundaryValue);

            // If x-component is missing, assume no velocity
            if (!found) {
                return null;
            }

            Func<double[], double, double>[] result = new Func<double[], double, double>[numberOfDimensions];
            result[0] = boundaryValue;

            // If x-component is present, _all_ other components have to be
            // there, too
            for (int i = 1; i < numberOfDimensions; i++) {
                result[i] = GetBoundaryValueFunction(boundaryValues, "u" + i);
            }

            return result;
        }

        #region IBoundaryConditionMap Members

        

        /// <summary>
        /// Retrieves the configured boundary condition for a given
        /// <paramref name="edgeTagName"/>.
        /// </summary>
        /// <param name="edgeTagName">The name of the edge</param>
        /// <returns>
        /// The boundary condition that has been assigned to edges with
        /// name <paramref name="edgeTagName"/>
        /// </returns>
        public virtual BoundaryCondition GetBoundaryCondition(string edgeTagName) {
            AppControl.BoundaryValueCollection boundaryValues =
                GetConfiguredBoundaryValues(edgeTagName);

            var key = boundaryValueMap.Keys.SingleOrDefault(tag => edgeTagName.StartsWith(tag, StringComparison.InvariantCultureIgnoreCase));
            if (key == null) {
                throw new NotImplementedException("Unknown edge tag name \"" + edgeTagName + "\"");
            } else {
                return boundaryValueMap[key](control, boundaryValues);
            }

        }

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
        public StateVector GetBoundaryState(byte EdgeTag, double time, double[] x, double[] normal, StateVector stateVector) {
            string edgeTagName = EdgeTagNames[EdgeTag];
            if (!ConditionMap.ContainsKey(edgeTagName)) {
                throw new ArgumentException("No boundary condition found for edge \"" + edgeTagName + "\"", "EdgeTag");
            }
            return ConditionMap[edgeTagName].GetBoundaryState(time, x, normal, stateVector);
        }

        /// <summary>
        /// <see cref="IBoundaryConditionMap"/>
        /// </summary>
        /// <param name="edgeTagName">
        /// <see cref="IBoundaryConditionMap"/>
        /// </param>
        /// <returns><see cref="IBoundaryConditionMap"/></returns>
        public BoundaryCondition this[string edgeTagName] {
            get {
                return ConditionMap[edgeTagName];
            }
        }

        /// <summary>
        /// <see cref="IBoundaryConditionMap"/>
        /// </summary>
        /// <param name="EdgeTag"></param>
        /// <returns></returns>
        public BoundaryCondition this[byte EdgeTag] {
            get {
                return this[EdgeTagNames[EdgeTag]];
            }
        }

        #endregion


        /*
        /// <summary>
        /// Mapping between edge tag names an the corresponding implementations
        /// of <see cref="BoundaryCondition"/>
        /// </summary>
        private static Dictionary<string, Func<CNSControl, AppControl.BoundaryValueCollection, BoundaryCondition>> boundaryValueMap =
            new Dictionary<string, Func<CNSControl, AppControl.BoundaryValueCollection, BoundaryCondition>>() {
                {
                    "adiabaticSlipWall",
                    (config, boundaryValues) => new AdiabaticSlipWall(
                        config,
                        GetOptionalVelocityBoundaryValueFunction(boundaryValues))
                },
                {
                    "symmetryPlane",
                    (config, boundaryValues) => new AdiabaticSlipWall(config)
                },
                {
                    "adiabaticWall",
                    (config, boundaryValues) => new AdiabaticWall(config)
                },
                {
                    "isothermalWall",
                    (config, boundaryValues) => new IsothermalWall(
                        config,
                        GetBoundaryValueFunction(boundaryValues, "T"),
                        GetOptionalVelocityBoundaryValueFunction(boundaryValues))
                },
                {
                    "subsonicInlet",
                    (config, boundaryValues) => new SubsonicInlet(
                        config,
                        GetBoundaryValueFunction(boundaryValues, "rho"),
                        GetVelocityBoundaryValueFunction(boundaryValues))
                },
                {
                    "subsonicPressureInlet",
                    (config, boundaryValues) => new SubsonicPressureInlet(
                        config,
                        GetBoundaryValueFunction(boundaryValues, "p0"),
                        GetBoundaryValueFunction(boundaryValues, "T0"))
                },
                {
                    "subsonicOutlet",
                    (config, boundaryValues) => new SubsonicOutlet(
                        config,
                        GetBoundaryValueFunction(boundaryValues, "p"))
                },
                {
                    "supersonicInlet",
                    (config, boundaryValues) => new SupersonicInlet(
                        config,
                        GetBoundaryValueFunction(boundaryValues, "rho"),
                        GetVelocityBoundaryValueFunction(boundaryValues),
                        GetBoundaryValueFunction(boundaryValues, "p"))
                },
                {
                    "supersonicOutlet",
                    (config, boundaryValues) => new SupersonicOutlet(config)
                }
            };

        */
    }
}
