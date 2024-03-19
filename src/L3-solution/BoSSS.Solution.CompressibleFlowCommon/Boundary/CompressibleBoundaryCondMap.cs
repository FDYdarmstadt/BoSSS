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
    /// works with instances of <see cref="CompressibleBoundaryCondMap"/>.
    /// </summary>
    public class CompressibleBoundaryCondMap : BoSSS.Solution.Utils.BoundaryCondMap<CompressibleBcType>, IBoundaryConditionMap {

        protected IGridData gridData;

        /// <summary>
        /// Mapping between edge tag names and associated boundary conditions
        /// - index: the edge tag, between 1 (including) and <see cref="GridCommons.FIRST_PERIODIC_BC_TAG"/> (excluding)
        /// - item: boundary condition for the respective edge tag
        /// </summary>
        public BoundaryCondition[] ConditionMap {
            get;
            private set;
        }

        /// <summary>
        /// Initialization of <see cref="ConditionMap"/>
        /// </summary>
        private void InitConditionMap() {
            ConditionMap = new BoundaryCondition[GridCommons.FIRST_PERIODIC_BC_TAG];

            foreach (byte edgeTag in gridData.EdgeTagNames.Keys) {
                // Only process non-periodic boundary edges
                if (edgeTag == 0 || edgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                    continue;
                }

                if (ConditionMap[edgeTag] != null) {
                    // already initialized
                    continue;
                }

                CompressibleBcType bcType = base.EdgeTag2Type[edgeTag];
                ConditionMap[edgeTag] = BoundaryConditionFactory(bcType, Material, edgeTag);
            }
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
        public virtual BoundaryCondition GetBoundaryCondition(string edgeTagName) {
            return ConditionMap[base.EdgeTagName2EdgeTag[edgeTagName]];
        }

        /// <summary>
        /// simulation properties
        /// </summary>
        protected MaterialProperty.Material Material {
            get;
            private set;
        }

        /// <summary>
        /// This calls for a better solution - sorry!
        /// </summary>
        static string[] bndFuncNames = new[] { "u0", "u1", "u2", "p", "T", "rho", "p0", "T0", "ringleb",
            "u0#A", "u1#A", "u2#A", "p#A", "T#A", "rho#A", "p0#A", "T0#A", "ringleb#A",
            "u0#B", "u1#B", "u2#B", "p#B", "T#B", "rho#B", "p0#B", "T0#B", "ringleb#B",            
            "u0#L", "u1#L", "u2#L", "p#L", "T#L", "rho#L", "p0#L", "T0#L", "ringleb#L",
            "u0#R", "u1#R", "u2#R", "p#R", "T#R", "rho#R", "p0#R", "T0#R", "ringleb#R"};

        /// <summary>
        /// Constructs a new map by searching through all the edge tags
        /// (<see cref="GridData.EdgeData.EdgeTags"/> and instantiating sub classes
        /// of <see cref="CompressibleBoundaryCondMap"/> specific for the compressible
        /// Navier-Stokes solver depending on their edge tag names.
        /// </summary>
        /// <param name="gridData">The omnipresent grid data</param>
        /// <param name="control">Configuration options</param>
        public CompressibleBoundaryCondMap(IGridData gridData, AppControl control, MaterialProperty.Material __material) : base(gridData, control.BoundaryValues, bndFuncNames) {
            this.gridData = gridData;
            this.Material = __material;
            InitConditionMap();
        }

        private Func<double[], double, double> GetBoundaryValueFunction(byte EdgeTag, string fieldName) {
            if (!base.bndFunction.ContainsKey(fieldName)) {
                throw new ArgumentException(
                    "Missing boundary specification for field '" + fieldName + "'",
                    "condition");
            }

            return base.bndFunction[fieldName][EdgeTag];
        }

        /// <summary>
        /// Vector version of <see cref="GetBoundaryValueFunction"/> for the
        /// velocity components
        /// </summary>
        private Func<double[], double, double>[] GetVelocityBoundaryValueFunction(byte EdgeTag, string hashPlusSpecies) {
            Func<double[], double, double>[] result = GetOptionalVelocityBoundaryValueFunction(EdgeTag, hashPlusSpecies);

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
        private Func<double[], double, double>[] GetOptionalVelocityBoundaryValueFunction(byte EdgeTag, string hashPlusSpecies) {
            int numberOfDimensions = gridData.SpatialDimension;

            // First check x-component only
            Func<double[], double, double> boundaryValue;
            bool found = base.bndFunction.ContainsKey("u0" + hashPlusSpecies);
            boundaryValue = base.bndFunction["u0" + hashPlusSpecies][EdgeTag];

            //bool found = TryGetBoundaryValueFunction(boundaryValues, "u0", out boundaryValue);

            // If x-component is missing, assume no velocity
            if (!found) {
                return null;
            }

            Func<double[], double, double>[] result = new Func<double[], double, double>[numberOfDimensions];
            result[0] = boundaryValue;

            // If x-component is present, _all_ other components have to be
            // there, too
            for (int i = 1; i < numberOfDimensions; i++) {
                result[i] = base.bndFunction["u" + i + hashPlusSpecies][EdgeTag];
            }

            return result;
        }

        /// <summary>
        /// Retrieves the configured boundary condition for a given
        /// <paramref name="edgeTag"/>.
        /// </summary>
        public virtual BoundaryCondition GetBoundaryCondition(byte edgeTag) {
            return ConditionMap[edgeTag];
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
            return ConditionMap[EdgeTag].GetBoundaryState(time, x, normal, stateVector);
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
                return GetBoundaryCondition(edgeTagName);
            }
        }

        /// <summary>
        /// <see cref="IBoundaryConditionMap"/>
        /// </summary>
        /// <param name="EdgeTag"></param>
        /// <returns></returns>
        public BoundaryCondition this[byte EdgeTag] {
            get {
                return GetBoundaryCondition(EdgeTag);
            }
        }

        /// <summary>
        /// Mapping between edge tag names an the corresponding implementations
        /// of <see cref="BoundaryCondition"/>
        /// </summary>
        protected BoundaryCondition BoundaryConditionFactory(CompressibleBcType bcType, MaterialProperty.Material material, byte EdgeTag, string speciesName = null) {
            //private static 
            //Dictionary<string, Func<MaterialProperty.Material, AppControl.BoundaryValueCollection, BoundaryCondition>> boundaryValueMap =
            //new Dictionary<string, Func<MaterialProperty.Material, AppControl.BoundaryValueCollection, BoundaryCondition>>() {

            string hashPlusSpecies = "";
            if (speciesName != null) {
                hashPlusSpecies = "#" + speciesName;
            }

            switch (bcType) {
                case CompressibleBcType.adiabaticSlipWall:
                    return new AdiabaticSlipWall(
                        material,
                        GetOptionalVelocityBoundaryValueFunction(EdgeTag, hashPlusSpecies));

                case CompressibleBcType.symmetryPlane:
                    return new AdiabaticSlipWall(material);

                case CompressibleBcType.adiabaticWall:
                    return new AdiabaticWall(material);

                case CompressibleBcType.isothermalWall:
                    return new IsothermalWall(
                        material,
                        GetBoundaryValueFunction(EdgeTag, "T" + hashPlusSpecies),
                        GetOptionalVelocityBoundaryValueFunction(EdgeTag, hashPlusSpecies));

                case CompressibleBcType.subsonicInlet:
                    return new SubsonicInlet(
                        material,
                        GetBoundaryValueFunction(EdgeTag, "rho" + hashPlusSpecies),
                        GetVelocityBoundaryValueFunction(EdgeTag, hashPlusSpecies));

                case CompressibleBcType.subsonicPressureInlet:
                    return new SubsonicPressureInlet(
                        material,
                        GetBoundaryValueFunction(EdgeTag, "p0" + hashPlusSpecies),
                        GetBoundaryValueFunction(EdgeTag, "T0" + hashPlusSpecies));

                case CompressibleBcType.subsonicOutlet:
                    return new SubsonicOutlet(
                        material,
                        GetBoundaryValueFunction(EdgeTag, "p" + hashPlusSpecies));

                case CompressibleBcType.supersonicInlet:
                    return new SupersonicInlet(
                        material,
                        GetBoundaryValueFunction(EdgeTag, "rho" + hashPlusSpecies),
                        GetVelocityBoundaryValueFunction(EdgeTag, hashPlusSpecies),
                        GetBoundaryValueFunction(EdgeTag, "p" + hashPlusSpecies));

                case CompressibleBcType.spaceTimeBoundary:
                    return new SpaceTimeBoundaryCondition(material);

                case CompressibleBcType.supersonicOutlet:
                    return new SupersonicOutlet(material);

                case CompressibleBcType.ringleb:
                    return new ExactRinglebBoundaryState(material);

                default:
                    throw new ArgumentException("unknown boundary type: " + bcType + hashPlusSpecies);
            }
        }
    }
}