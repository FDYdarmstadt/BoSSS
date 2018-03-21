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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Queries {

    /// <summary>
    /// Represents a query to be evaluates at some instant in time.
    /// </summary>
    /// <param name="application"></param>
    /// <param name="time">
    /// The physical time of the query
    /// </param>
    /// <returns>
    /// The query result
    /// </returns>
    public delegate double Query(IApplication<AppControl> application, double time);

    /// <summary>
    /// A library of useful standard queries
    /// </summary>
    public static class QueryLibrary {

        /// <summary>
        /// Computes the L2 norm of field <paramref name="fieldName"/> in the
        /// domain defined by <paramref name="maskFactory"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="maskFactory">
        /// A function that returns the domain of integration in terms of cells
        /// to be evaluated. If null, the full computational domain will be
        /// considered.
        /// </param>
        /// <returns>
        /// The L2 norm of <paramref name="fieldName"/>.
        /// </returns>
        public static Query L2Norm(string fieldName, Func<GridData, CellMask> maskFactory = null) {
            return (app, time) => GetField(app.IOFields, fieldName).L2Norm(GetMaskOrNull(maskFactory, app.GridData));
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to the given <paramref name="referenceField"/> in the
        /// domain defined by <paramref name="maskFactory"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="referenceField">
        /// A DG field containing the exact solution.
        /// </param>
        /// <param name="maskFactory">
        /// A function that returns the domain of integration in terms of cells
        /// to be evaluated. If null, the full computational domain will be
        /// considered.
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, DGField referenceField, Func<GridData, CellMask> maskFactory = null) {
            return delegate(IApplication<AppControl> app, double time) {
                DGField differenceField = referenceField - GetField(app.IOFields, fieldName);
                return differenceField.L2Norm(GetMaskOrNull(maskFactory, app.GridData));
            };
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to the given, time-dependent
        /// <paramref name="referenceSolution"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="referenceSolution">
        /// The reference solution
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, Func<double[], double, double> referenceSolution) {
            return (app, time) =>
                GetField(app.IOFields, fieldName).L2Error(referenceSolution.Vectorize(time));
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to the given, time-dependent
        /// <paramref name="referenceSolution"/> using a quadrature rule of
        /// order <paramref name="quadratureOrder"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="referenceSolution">
        /// The reference solution
        /// </param>
        /// <param name="quadratureOrder">
        /// The desired order of the integration rule
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, Func<double[], double, double> referenceSolution, int quadratureOrder) {
            return (app, time) =>
                GetField(app.IOFields, fieldName).L2Error(referenceSolution.Vectorize(time), quadratureOrder);
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to the given <paramref name="referenceSolution"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="referenceSolution">
        /// The reference solution
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, Func<double[], double> referenceSolution) {
            return (app, time) =>
                GetField(app.IOFields, fieldName).L2Error(referenceSolution.Vectorize());
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to the given <paramref name="referenceSolution"/> using a
        /// quadrature rule of order <paramref name="quadratureOrder"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="referenceSolution">
        /// The reference solution
        /// </param>
        /// <param name="quadratureOrder">
        /// The desired order of the integration rule
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, Func<double[], double> referenceSolution, int quadratureOrder) {
            return (app, time) =>
                GetField(app.IOFields, fieldName).L2Error(referenceSolution.Vectorize(), quadratureOrder);
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to a field with the same name stored in the time-step with
        /// <see cref="Guid"/> <paramref name="timestepGuid"/> within the
        /// current database.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="timestepGuid">
        /// The time-step id containing the reference field
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, Guid timestepGuid) {
            return delegate(IApplication<AppControl> app, double time) {
                DGField field = app.IOFields.Single(f => f.Identification == fieldName);
                DGField referenceField = GetStoredField(app.GridData, timestepGuid, fieldName);

                return L2Error(fieldName, referenceField)(app, time);
            };
        }

        /// <summary>
        /// Computes the L2 error of field <paramref name="fieldName"/> with
        /// respect to a field with the same name stored in the time-step with
        /// <see cref="TimestepNumber"/> <paramref name="timestepNumber"/> within the
        /// current database.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="timestepNumber">
        /// The time-step containing the reference field
        /// </param>
        /// <returns>
        /// The error of <paramref name="fieldName"/> in the L2 norm.
        /// </returns>
        public static Query L2Error(string fieldName, TimestepNumber timestepNumber) {
            return delegate (IApplication<AppControl> app, double time) {
                ITimestepInfo ts = app.CurrentSessionInfo.Timesteps.Single(t => t.TimeStepNumber.Equals(timestepNumber));
                DGField field = app.IOFields.Single(f => f.Identification == fieldName);
                DGField referenceField = GetStoredField(app.GridData, ts.ID, fieldName);

                return L2Error(fieldName, referenceField)(app, time);
            };
        }

        /// <summary>
        /// Evaluates the integral of field <paramref name="fieldName"/> over
        /// the domain defined by <paramref name="maskFactory"/>.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="maskFactory">
        /// A function that returns the domain of integration in terms of cells
        /// to be evaluated. If null, the full computational domain will be
        /// considered.
        /// </param>
        /// <returns>
        /// The integral of <paramref name="fieldName"/> over the cells defined
        /// by <paramref name="maskFactory"/>.
        /// </returns>
        public static Query Integral(string fieldName, Func<GridData, CellMask> maskFactory = null) {
            return (app, time) => GetField(app.IOFields, fieldName).IntegralOver(GetMaskOrNull(maskFactory, app.GridData));
        }

        /// <summary>
        /// Evaluates the field with name <paramref name="fieldName"/> at the
        /// given point.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <param name="point">
        /// The spatial coordinates
        /// </param>
        /// <returns>
        /// The value of field <paramref name="fieldName"/> at
        /// x=<paramref name="point"/>
        /// </returns>
        public static Query Probe(string fieldName, params double[] point) {
            return (app, time) => GetField(app.IOFields, fieldName).ProbeAt(point);
        }

        /// <summary>
        /// Evaluates the total number of DOF for field
        /// <paramref name="fieldName"/>. Useful for automated analysis of 
        /// parameter studies.
        /// </summary>
        /// <param name="fieldName">
        /// The name of the field to be evaluated
        /// </param>
        /// <returns>
        /// The total number of DOF.
        /// </returns>
        public static Query NumberOfDOF(string fieldName) {
            return delegate(IApplication<AppControl> app, double time) {
                long noOfCells = app.GridData.Grid.NumberOfCells;
                int degree = GetField(app.IOFields, fieldName).Basis.Degree;

                int DOFPerCell;
                switch (app.GridData.SpatialDimension) {
                    case 1:
                        DOFPerCell = degree + 1;
                        break;

                    case 2:
                        DOFPerCell = ((degree + 1) * (degree + 2)) / 2;
                        break;

                    case 3:
                        DOFPerCell = ((degree + 1) * (degree + 2) * (degree + 3)) / 6;
                        break;

                    default:
                        throw new Exception("Invalid number of dimensions");
                }

                return noOfCells * DOFPerCell;
            };
        }

        /// <summary>
        /// Returns the result of <paramref name="maskFactory"/>, if
        /// <paramref name="maskFactory"/> is not null. Otherwise, returns
        /// null.
        /// </summary>
        /// <param name="maskFactory"></param>
        /// <param name="gridData"></param>
        /// <returns></returns>
        private static CellMask GetMaskOrNull(Func<GridData, CellMask> maskFactory, GridData gridData) {
            return (maskFactory == null) ? null : maskFactory(gridData);
        }

        /// <summary>
        /// Retrieves the field with given <paramref name="fieldName"/>
        /// </summary>
        /// <param name="fields"></param>
        /// <param name="fieldName"></param>
        /// <returns></returns>
        private static DGField GetField(IEnumerable<DGField> fields, string fieldName) {
            return fields.Single(f => f.Identification == fieldName);
        }

        /// <summary>
        /// Loads the field with name <paramref name="fieldName"/> within the
        /// time-step with id <paramref name="timestepGuid"/> from the
        /// database.
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="timestepGuid"></param>
        /// <param name="fieldName"></param>
        /// <returns></returns>
        private static DGField GetStoredField(GridData gridData, Guid timestepGuid, string fieldName) {
            IDatabaseInfo database = gridData.Grid.Database;
            ITimestepInfo tsi = database.Controller.DBDriver.LoadTimestepInfo(
                timestepGuid, null, database);
            return database.Controller.DBDriver.LoadFields(tsi, gridData, new[] { fieldName }).Single();
        }
    }
}
