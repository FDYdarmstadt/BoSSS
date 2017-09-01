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

using BoSSS.Foundation;
using BoSSS.Solution.Control;
using ilPSP.Tracing;
using log4net;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace BoSSS.Solution.Queries {

    /// <summary>
    /// Handles the creation and evaluation of queries (see <see cref="Query"/>)
    /// which can be used to evaluate certain quantities of interest after the
    /// BoSSS run has finished
    /// </summary>
    public class QueryHandler {

        static ILog Logger = LogManager.GetLogger(typeof(QueryHandler));

        /// <summary>
        /// The currently running application.
        /// </summary>
        protected readonly IApplication<AppControl> Application;

        /// <summary>
        /// A mapping between a query id and its associated query (new format).
        /// </summary>
        public Dictionary<string, Query> QueryMap {
            get;
            private set;
        }

        /// <summary>
        /// A mapping between a query id and the result of its evaluation
        /// (after <see cref="EvaluateQueries"/> has been called).
        /// </summary>
        public Dictionary<string, object> QueryResults {
            get;
            private set;
        }

        /// <summary>
        /// Constructs a new <see cref="QueryHandler"/>.
        /// </summary>
        /// <param name="app"></param>
        public QueryHandler(IApplication<AppControl> app) {
            this.Application = app;
            QueryMap = new Dictionary<string, Query>();
            QueryResults = new Dictionary<string, object>();
        }

        /// <summary>
        /// Evaluates all registered queries and collects
        /// the results into <see cref="QueryResults"/>.
        /// </summary>
        /// <param name="fields">
        /// DG fields required for the evaluation of the queries.
        /// </param>
        /// <param name="time">
        /// Physical simulation time, required for time dependent evaluation of the queries
        /// </param>
        /// <returns>
        /// A mapping between query ids and their
        /// respective evaluation results.
        /// </returns>
        public void EvaluateQueries(IEnumerable<DGField> fields, double time) {
            using (new FuncTrace()) {
                foreach (var idQueryPair in QueryMap) {
                    QueryResults.Add(idQueryPair.Key, idQueryPair.Value(Application, time));
                }
            }
        }

        /// <summary>
        /// Writes the given query results to the given log file.
        /// </summary>
        /// <param name="resultsMap">
        /// The query results.
        /// </param>
        /// <param name="log">
        /// A log file.
        /// </param>
        static public void LogQueryResults(Dictionary<string, object> resultsMap, TextWriter log) {
            // a 'table'-style 
            var SortedKeys = new List<string>(resultsMap.Keys);
            SortedKeys.Sort();

            for (int i = 0; i < SortedKeys.Count; i++) {
                log.Write(SortedKeys[i]);
                if (i < SortedKeys.Count - 1)
                    log.Write("\t");
            }
            log.WriteLine();

            for (int i = 0; i < SortedKeys.Count; i++) {
                var val = resultsMap[SortedKeys[i]];
                if (val is double) {
                    log.Write(((double)val).ToString("E16", NumberFormatInfo.InvariantInfo));
                } else {
                    log.Write(val.ToString());
                }
                if (i < SortedKeys.Count - 1)
                    log.Write("\t");
            }
            log.WriteLine();
        }
        
        /// <summary>
        /// Adds the given query (new format) to the list of queries to be
        /// evaluated.
        /// </summary>
        /// <param name="id">
        /// The unique id of the query
        /// </param>
        /// <param name="query">
        /// The query function to be evaluated
        /// </param>
        public void AddQuery(string id, Query query) {
            if (QueryMap.ContainsKey(id)) {
                throw new Exception("A query with id '" + id + "' already exists");
            }
            QueryMap.Add(id, query);
        }

        public void ValueQuery(string name, double value, bool OverWriteIfExistent = false) {
            Query query = (app, t) => value;
            if (QueryMap.ContainsKey(name)) {
                if (!OverWriteIfExistent) {
                    throw new Exception("A query with id '" + name + "' already exists");
                } else {
                    QueryMap[name] = query;
                }
            } else {
                QueryMap.Add(name, query);
            }
        }
    }
}
