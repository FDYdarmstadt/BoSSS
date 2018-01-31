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

using BoSSS.Foundation.IO;
using BoSSS.Solution.Statistic;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Extension methods for <see cref="ITimestepInfo"/>
    /// </summary>
    public static class ITimestepInfoExtensions {

        /// <summary>
        /// Deletes a time-step from its database.
        /// </summary>
        /// <param name="timestep">The time-step to be deleted.</param>
        /// <param name="force">
        /// If true, the user will not be asked for confirmation
        /// </param>
        public static void Delete(this ITimestepInfo timestep, bool force = false) {
            bool sure = true;
            if (!force) {
                Console.WriteLine("Time-step: " + timestep.ToString());
                Console.Write("Do you really want to delete this time-step? [y/n]: ");
                string line = Console.ReadLine();

                sure = line.ToLower().Equals("y");
            }

            if (sure) {
                TimestepNumber number = timestep.TimeStepNumber;
                timestep.Database.Controller.DeleteTimestep(timestep);
                Console.WriteLine("Time-step " + number + " deleted.");
            } else {
                Console.WriteLine("Session delete canceled.");
            }
        }

        /// <summary>
        /// Deletes all entities inside a collection, provided that they are
        /// database entities that can be deleted, i.e. ISessionInfo and IGridInfo
        /// as of now.
        /// The user is not asked to confirm this action once, and not for every
        /// individual element in the collection.
        /// </summary>
        /// <param name="timesteps">The entities to be deleted.</param>
        public static void DeleteAll(this IEnumerable<ITimestepInfo> timesteps) {
            if (timesteps.IsNullOrEmpty()) {
                Console.WriteLine("Given collection is empty; nothing to delete");
                return;
            }

            Console.WriteLine("Time-steps to delete:");
            foreach (ITimestepInfo timestep in timesteps) {
                Console.WriteLine(timestep.ToString());
            }
            Console.Write("Do you really want to delete these time-steps? [y/n]: ");
            string line = Console.ReadLine();

            if (line.ToLower().Equals("y")) {
                foreach (ITimestepInfo timestep in timesteps) {
                    timestep.Database.Controller.DeleteTimestep(timestep);
                }
                Console.WriteLine("Time-steps successfully deleted.");
            } else {
                Console.WriteLine("Time-step delete canceled.");
            }
        }

        /// <summary>
        /// Convenience interface to create a
        /// <see cref="SessionExportInstruction"/> for a single time-step.
        /// </summary>
        /// <param name="timestep">
        /// The time-step to be exported
        /// </param>
        /// <returns>
        /// A new instance of <see cref="SessionExportInstruction"/> for the
        /// given <paramref name="timestep"/>
        /// </returns>
        public static SessionExportInstruction Export(this ITimestepInfo timestep) {
            return new SessionExportInstruction(timestep.Session).WithTimesteps(timestep.TimeStepNumber);
        }

        /// <summary>
        /// Convenience interface to create a
        /// <see cref="SessionExportInstruction"/> for a list of time-steps.
        /// </summary>
        /// <param name="timesteps">
        /// A non-empty list of time-steps to be exported
        /// </param>
        /// <returns>
        /// A new instance of <see cref="SessionExportInstruction"/> for the
        /// given <paramref name="timesteps"/>
        /// </returns>
        public static SessionExportInstruction Export(this IEnumerable<ITimestepInfo> timesteps) {
            if (timesteps.IsNullOrEmpty()) {
                throw new ArgumentException("Sequence must contain at least one time-step");
            }

            if (!timesteps.IsUniform(t => t.Session)) {
                throw new ArgumentException(
                    "All time-steps that are exported simultaneously must belong to the same session");
            }

            return new SessionExportInstruction(timesteps.First().Session).WithTimesteps(
                timesteps.Select(t => t.TimeStepNumber).ToArray());
        }

        /// <summary>
        /// Opens the directory where the export for the selected
        /// <paramref name="timestep"/> are stored in the explorer.
        /// </summary>
        /// <param name="timestep">
        /// The selected time-step.
        /// </param>
        /// <remarks>
        /// Obviously, this only works in Windows environments.
        /// </remarks>
        public static void OpenExportDirectory(this ITimestepInfo timestep) {
            Process.Start(Utils.GetExportDirectory(timestep.Session));
        }

        /// <summary>
        /// Finds the time-steps with the given
        /// <paramref name="timeStepNumber"/> within
        /// <paramref name="timesteps"/>
        /// </summary>
        /// <param name="timesteps">List to be searched</param>
        /// <param name="timeStepNumber">Number to be searched for</param>
        /// <returns>
        /// The single time-step within <paramref name="timesteps"/> with
        /// time-step number <paramref name="timeStepNumber"/>.
        /// </returns>
        public static ITimestepInfo Find(this IEnumerable<ITimestepInfo> timesteps, TimestepNumber timeStepNumber) {
            return timesteps.Single(t => t.TimeStepNumber.Equals(timeStepNumber));
        }

        /// <summary>
        /// Filters all elements of <paramref name="timesteps"/> which are
        /// sub-steps, i.e. which have a <see cref="TimestepNumber"/> with a
        /// length greater than 1.
        /// </summary>
        /// <param name="timesteps">
        /// The sequence of time-steps to be filtered
        /// </param>
        /// <returns>
        /// A sequence of full time-steps.
        /// </returns>
        public static IEnumerable<ITimestepInfo> WithoutSubSteps(this IEnumerable<ITimestepInfo> timesteps) {
            return timesteps.Where(t => t.TimeStepNumber.Length == 1);
        }

        /// <summary>
        /// Selects the previous time-step stored in the session.
        /// </summary>
        /// <param name="timestep">
        /// The current time-step
        /// </param>
        /// <returns>
        /// The last full time-step in the current session before
        /// <paramref name="timestep"/>.
        /// </returns>
        public static ITimestepInfo Previous(this ITimestepInfo timestep) {
            return timestep.Session.Timesteps.WithoutSubSteps().TakeWhile(t => t != timestep).Last();
        }

        /// <summary>
        /// Selects the next time-step stored in the session.
        /// </summary>
        /// <param name="timestep">
        /// The current time-step
        /// </param>
        /// <returns>
        /// The first full time-step in the current session after
        /// <paramref name="timestep"/>.
        /// </returns>
        public static ITimestepInfo Next(this ITimestepInfo timestep) {
            return timestep.Session.Timesteps.WithoutSubSteps().SkipWhile(t => t != timestep).Second();
        }

        /// <summary>
        /// Determines the mean time-step size between <paramref name="timestep"/>
        /// and <see cref="ITimestepInfoExtensions.Previous"/>
        /// </summary>
        /// <param name="timestep">
        /// The considered time-step
        /// </param>
        /// <returns>
        /// The average time-step in the given based on
        /// <see cref="ITimestepInfo.PhysicalTime"/>.
        /// </returns>
        public static double GetTimeStepSize(this ITimestepInfo timestep) {
#if DEBUG
            Console.WriteLine("Warning: Computation of time-step size uses averaging at the moment. Ask Björn for details");
#endif
            ITimestepInfo previous = timestep.Previous();
            int noOfSteps = timestep.TimeStepNumber.MajorNumber - previous.TimeStepNumber.MajorNumber;
            return (timestep.PhysicalTime - previous.PhysicalTime) / noOfSteps;
        }

        /// <summary>
        /// Converts a list of time-steps into a <see cref="Plot2Ddata"/>
        /// </summary>
        /// <param name="timesteps">
        /// A list of time-steps
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="ySelector">
        /// Selector for the relevant data at the data points.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and <paramref name="ySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> xSelector, Func<ITimestepInfo, double> ySelector) {
            List<double> xValues = new List<double>(timesteps.Count());
            List<double> yValues = new List<double>(xValues.Count);

            foreach (var timestep in timesteps) {
                try {
                    xValues.Add(xSelector(timestep));
                    yValues.Add(ySelector(timestep));
                } catch (Exception e) {
                    Console.WriteLine(
                        "Extracting requested information from time-step"
                            + " {0} of session {1} failed with message "
                            + " '{2}'. Proceeding without this data point",
                        timestep.TimeStepNumber,
                        timestep.Session.ID,
                        e.Message);
                }
            }

            return new Plot2Ddata(xValues, yValues);
        }

        /// <summary>
        /// Converts a list of time-steps into a <see cref="Plot2Ddata"/> while
        /// grouping the results using <paramref name="groupKeySelector"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A list of time-steps
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="ySelector">
        /// Selector for the relevant data at the data points.
        /// </param>
        /// <param name="groupKeySelector">
        /// A function defining a group id for each key-value pair
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and <paramref name="ySelector"/>,
        /// grouped by means of <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> xSelector, Func<ITimestepInfo, double> ySelector, Func<ITimestepInfo, string> groupKeySelector) {
            var OrderedTimeSteps = timesteps.
                Select(t => new KeyValuePair<double, ITimestepInfo>(xSelector(t), t)).
                OrderBy(p => p.Key).
                ToArray();

            // Process each item individually so we can swallow exceptions
            var data = new List<KeyValuePair<string, double[][]>>(timesteps.Count());
            foreach (var group in OrderedTimeSteps.GroupBy(p => groupKeySelector(p.Value))) {
                var xyPairs = new List<KeyValuePair<double, double>>(group.Count());
                foreach (var pair in group) {
                    try {
                        xyPairs.Add(new KeyValuePair<double, double>(pair.Key, ySelector(pair.Value)));
                    } catch (Exception e) {
                        Console.WriteLine(
                            "Extracting requested information from time-step"
                                + " {0} of session {1} failed with message "
                                + " '{2}'. Proceeding without this data point",
                            pair.Value.TimeStepNumber,
                            pair.Value.Session.ID,
                            e.Message);
                    }
                }

                // Group might contain only entries with failed conversions
                if (xyPairs.Count > 0) {
                    data.Add(new KeyValuePair<string, double[][]>(
                        group.Key,
                        new double[][] {
                        xyPairs.Select(p => p.Key).ToArray(),
                        xyPairs.Select(p => p.Value).ToArray()
                    }));
                }
            }

            return new Plot2Ddata(data.ToArray());
        }

        /// <summary>
        /// Converts a list of time-steps into a <see cref="Plot2Ddata"/> based on
        /// the results a query named <paramref name="queryName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A list of time-steps
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="queryName">
        /// Name of the query whose results will be used as error measure.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and the results of the query named
        /// <paramref name="queryName"/>.
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> xSelector, string queryName) {
            return timesteps.ToDataSet(
                xSelector,
                t => t.Session.QueryResults()[queryName]);
        }

        /// <summary>
        /// Converts a list of time-steps into a <see cref="Plot2Ddata"/> while
        /// grouping the results using the DG degree of the field identified
        /// by <paramref name="groupingFieldName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A list of time-steps
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="ySelector">
        /// Selector for the relevant data at the data points.
        /// </param>
        /// <param name="groupingFieldName">
        /// A name of a DG field present in all <paramref name="timesteps"/>
        /// whose DG degree will be used as a grouping function.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and <paramref name="ySelector"/>,
        /// grouped by means of the DG degree of the field named
        /// <paramref name="groupingFieldName"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> xSelector, Func<ITimestepInfo, double> ySelector, string groupingFieldName) {
            return timesteps.ToDataSet(
                xSelector,
                ySelector,
                t => t.Fields.Find(groupingFieldName).Basis.Degree.ToString());
        }

        /// <summary>
        /// Converts a list of time-steps into a <see cref="Plot2Ddata"/> based on
        /// the results of a query named <paramref name="queryName"/>, while
        /// grouping the results using <paramref name="groupKeySelector"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A list of time-steps
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="queryName">
        /// Name of the query whose results will be used as error measure.
        /// </param>
        /// <param name="groupKeySelector">
        /// A function defining a group id for each key-value pair
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and the results of the query named
        /// <paramref name="queryName"/>, grouped by means of
        /// <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> xSelector, string queryName, Func<ITimestepInfo, string> groupKeySelector) {
            return timesteps.ToDataSet(
                xSelector,
                t => t.Session.QueryResults()[queryName],
                groupKeySelector);
        }

        /// <summary>
        /// Converts a list of time-steps into a <see cref="Plot2Ddata"/> based on
        /// the results a query named <paramref name="queryName"/>, while
        /// grouping the results using the DG degree of the field identified
        /// by <paramref name="groupingFieldName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A list of time-steps
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="queryName">
        /// Name of the query whose results will be used as error measure.
        /// </param>
        /// <param name="groupingFieldName">
        /// A name of a DG field present in all <paramref name="timesteps"/>
        /// whose DG degree will be used as a grouping function.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and the results of the query named
        /// <paramref name="queryName"/>, grouped by means of the DG degree of
        /// the field named <paramref name="groupingFieldName"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> xSelector, string queryName, string groupingFieldName) {
            return timesteps.ToDataSet(
                xSelector,
                t => t.Session.QueryResults()[queryName],
                groupingFieldName);
        }

        /// <summary>
        /// Extracts information about the grid convergence of the given
        /// <paramref name="timesteps"/> with respect to the given
        /// <paramref name="errorFunctional"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a grid convergence study. The
        /// DG degree of the relevant fields should be the same for all
        /// time-steps.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of
        /// <paramref name="errorFunctional"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> errorFunctional) {
            return timesteps.ToDataSet(t => t.Grid.GetMeshSize(), errorFunctional).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts information about the grid convergence of the given
        /// <paramref name="timesteps"/> with respect to the given
        /// <paramref name="errorFunctional"/>, while grouping the results
        /// using <paramref name="groupKeySelector"/>
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a grid convergence study. The
        /// DG degree of the relevant fields should be the same for all
        /// time-steps.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <param name="groupKeySelector">
        /// A grouping function.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of
        /// <paramref name="errorFunctional"/>, grouped by means of
        /// <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> errorFunctional, Func<ITimestepInfo, string> groupKeySelector) {
            return timesteps.ToDataSet(t => t.Grid.GetMeshSize(), errorFunctional, groupKeySelector).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts information about the grid convergence of the given
        /// <paramref name="timesteps"/> with respect to the given
        /// <paramref name="errorFunctional"/> and groups the results according
        /// to the DG degree of a field named
        /// <paramref name="groupingFieldName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a grid convergence study.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <param name="groupingFieldName">
        /// A name of a relevant DG field which is present in all
        /// <paramref name="timesteps"/> and whose DG degree will be used as a
        /// group key.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/>, the
        /// values are determined via the logarithm of
        /// <paramref name="errorFunctional"/> and the results are grouped with
        /// respect to the DG degree of a field named
        /// <paramref name="groupingFieldName"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> errorFunctional, string groupingFieldName) {
            return timesteps.ToDataSet(t => t.Grid.GetMeshSize(), errorFunctional, groupingFieldName).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts information about the grid convergence of the given
        /// <paramref name="timesteps"/> with respect to the results of the
        /// query named <paramref name="queryName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a grid convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used an error measure.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of the results of a query
        /// named <paramref name="queryName"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, string queryName) {
            return timesteps.ToDataSet(t => t.Grid.GetMeshSize(), queryName).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts information about the grid convergence of the given
        /// <paramref name="timesteps"/> with respect to the results of the
        /// query named <paramref name="queryName"/>, while
        /// grouping the results using the DG degree of the field identified
        /// by <paramref name="groupingFieldName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a grid convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used as an error measure.
        /// </param>
        /// <param name="groupingFieldName">
        /// A name of a DG field present in all <paramref name="timesteps"/>
        /// whose DG degree will be used as a grouping function.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/>, the
        /// values are determined via the logarithm of the results of a query
        /// named <paramref name="queryName"/> and the results are grouped with
        /// respect to the DG degree of a field named
        /// <paramref name="groupingFieldName"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, string queryName, string groupingFieldName) {
            return timesteps.ToDataSet(t => t.Grid.GetMeshSize(), queryName, groupingFieldName).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts information about the grid convergence of the given
        /// <paramref name="timesteps"/> with respect to the results of the
        /// query named <paramref name="queryName"/>, while grouping the results
        /// using <paramref name="groupKeySelector"/>
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a grid convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used as an error measure.
        /// </param>
        /// <param name="groupKeySelector">
        /// A grouping function.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of the results of a query
        /// named <paramref name="queryName"/>, grouped by means of
        /// <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, string queryName, Func<ITimestepInfo, string> groupKeySelector) {
            return timesteps.ToDataSet(t => t.Grid.GetMeshSize(), queryName, groupKeySelector).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts information about the time convergence of the given 
        /// <paramref name="timesteps"/> with respect to the given
        /// <paramref name="errorFunctional"/>
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a time convergence study.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="ITimestepInfoExtensions.GetTimeStepSize"/>
        /// and the values are determined via the logarithm of
        /// <paramref name="errorFunctional"/>.
        /// </returns>
        public static Plot2Ddata ToTimeConvergenceData(this IEnumerable<ITimestepInfo> timesteps, Func<ITimestepInfo, double> errorFunctional) {
            return timesteps.
                ToDataSet(t => t.GetTimeStepSize(), errorFunctional, t => t.GridID.ToString()).
                WithLogX().
                WithLogY();
        }

        /// <summary>
        /// Extracts information about the time convergence of the given 
        /// <paramref name="timesteps"/> with respect to the results of the
        /// query named <paramref name="queryName"/>
        /// </summary>
        /// <param name="timesteps">
        /// A set of time-steps representing a time convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used as an error measure.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="ITimestepInfoExtensions.GetTimeStepSize"/>
        /// and the values are determined via the logarithm of the results of a
        /// query named <paramref name="queryName"/>
        /// </returns>
        public static Plot2Ddata ToTimeConvergenceData(this IEnumerable<ITimestepInfo> timesteps, string queryName) {
            return timesteps.
                ToDataSet(
                    t => t.GetTimeStepSize(),
                    queryName,
                    t => t.Session.KeysAndQueries["ExplicitScheme"].ToString() + " " 
                        + t.Session.KeysAndQueries["ExplicitOrder"].ToString()).
                WithLogX().
                WithLogY();
        }

        /// <summary>
        /// Estimates the errors of the fields with name
        /// <paramref name="fieldName"/> in the given
        /// <paramref name="timesteps"/> by computing the errors with respect
        /// to the solution on the finest corresponding grid by making use of
        /// <see cref="BoSSS.Solution.Statistic.DGFieldComparison.ComputeErrors"/>.
        /// The result is then grouped according to the polynomial degree of field
        /// <paramref name="fieldName"/>.
        /// </summary>
        /// <param name="timesteps">
        /// The time-steps containing the fields whose errors should be
        /// estimated.
        /// </param>
        /// <param name="fieldName">
        /// The name of the DG field whose error should be estimated.
        /// </param>
        /// <returns>
        /// A data set containing information about the grid resolution and the
        /// corresponding errors with respect to the finest corresponding grid,
        /// grouped by the polynomial degree. Obviously, the time-step
        /// associated with the finest grid for each polynomial degree has an
        /// estimated error of zero (by definition) and is thus excluded from
        /// the result.
        /// </returns>
        public static Plot2Ddata ToEstimatedGridConvergenceData(this IEnumerable<ITimestepInfo> timesteps, string fieldName) {
            Dictionary<string, double[][]> dataGroups = new Dictionary<string, double[][]>();
            foreach (var group in timesteps.GroupBy(t => t.Fields.Find(fieldName).Basis.Degree)) {
                double[] resolution;
                Dictionary<string, double[]> errors;
                DGFieldComparison.ComputeErrors(
                    new string[] { fieldName },
                    group.ToArray(),
                    out resolution,
                    out errors);

                Debug.Assert(errors.ContainsKey(fieldName));
                Debug.Assert(errors[fieldName].Length == resolution.Length);

                double[][] resolutionsAndErrors = new double[2][] { resolution, errors[fieldName] };
                dataGroups.Add(group.Key.ToString(), resolutionsAndErrors);
            }

            return new Plot2Ddata(dataGroups.ToArray()).WithLogX().WithLogY();
        }

        /// <summary>
        /// Extracts time-steps with a time-step number greater than or equal
        /// to <paramref name="number"/> from the given list of
        /// <paramref name="sortedTimesteps"/>.
        /// </summary>
        /// <param name="sortedTimesteps">
        /// A list of time-steps, in ascending order of their time-step number
        /// </param>
        /// <param name="number">
        /// The minimum time-step number to be selected.
        /// </param>
        /// <returns>
        /// The time-steps with a time-step number greater than or equal
        /// to <paramref name="number"/>.
        /// </returns>
        public static IEnumerable<ITimestepInfo> From(this IEnumerable<ITimestepInfo> sortedTimesteps, TimestepNumber number) {
            var enumerator = sortedTimesteps.GetEnumerator();
            while (enumerator.MoveNext()) {
                if (enumerator.Current.TimeStepNumber.CompareTo(number) >= 0) {
                    yield return enumerator.Current;
                    break;
                }
            }

            while (enumerator.MoveNext()) {
                yield return enumerator.Current;
            }
        }

        /// <summary>
        /// Extracts time-steps with a time-step number smaller than or equal
        /// to <paramref name="number"/> from the given list of
        /// <paramref name="sortedTimesteps"/>.
        /// </summary>
        /// <param name="sortedTimesteps">
        /// A list of time-steps, in ascending order of their time-step number
        /// </param>
        /// <param name="number">
        /// The minimum time-step number to be selected.
        /// </param>
        /// <returns>
        /// The time-steps with a time-step number greater than or equal
        /// to <paramref name="number"/>.
        /// </returns>
        public static IEnumerable<ITimestepInfo> To(this IEnumerable<ITimestepInfo> sortedTimesteps, TimestepNumber number) {
            var enumerator = sortedTimesteps.GetEnumerator();
            while (enumerator.MoveNext() && enumerator.Current.TimeStepNumber.CompareTo(number) <= 0) {
                yield return enumerator.Current;
            }
        }
    }
}
