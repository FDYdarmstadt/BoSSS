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
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Standard implementation of a database info object.
    /// </summary>
    public class DatabaseInfo : IDatabaseInfo {

        /// <summary>
        /// Stores the path
        /// </summary>
        /// <param name="path">Path to the database</param>
        public DatabaseInfo(string path) {
            this.Path = path;
            if (path == null) {
                Controller = NullDatabaseController.Instance;
            } else {
                Controller = new DatabaseController(this);
            }
        }

        /// <summary>
        /// Full path to the base directory of the database.
        /// </summary>
        public string Path {
            get;
            private set;
        }

        /// <summary>
        /// Provides functionality to copy/move/delete info objects stored in
        /// the database
        /// </summary>
        public IDatabaseController Controller {
            get;
            private set;
        }

        /// <summary>
        /// Returns a string representation of this database.
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            // TO DO: THIS REQUIRES IO OPERATIONS => BAD!
            //return "{ Session Count = " + Controller.Sessions.Count()
            //    + "; Grid Count = " + Controller.Grids.Count()
            //    + "; Path = " + Path + " }";

            // Temporary workaround
            return "{ Session Count = " + ((DatabaseController)Controller).SessionCount
                + "; Grid Count = " + ((DatabaseController)Controller).GridCount
                + "; Path = " + Path + " }";
        }

        /// <summary>
        /// The sessions of this database.
        /// </summary>
        public IEnumerable<ISessionInfo> Sessions {
            get {
                //Stopwatch stw = new Stopwatch();

                //Console.WriteLine("aquire sessions...");
                //stw.Reset();
                //stw.Start();
                var allsessions = Controller.Sessions;
                //stw.Stop();
                //Console.WriteLine("done. " + stw.ElapsedMilliseconds);

                //Console.WriteLine("sorting sessions...");
                //stw.Reset();
                //stw.Start();
                var R = allsessions.OrderByDescending(s => s.WriteTime);
                //stw.Stop();
                //Console.WriteLine("done. " + stw.ElapsedMilliseconds);

                return R;
            }
        }
        
        /// <summary>
        /// The grids of this database.
        /// </summary>
        public IEnumerable<IGridInfo> Grids {
            get {
                return Controller.Grids.OrderByDescending(g => g.WriteTime);
            }
        }

        /// <summary>
        /// Sessions sorted according to projects, see <see cref="ISessionInfo.ProjectName"/>.
        /// </summary>
        public IDictionary<string, IEnumerable<ISessionInfo>> Projects {
            get {
                Dictionary<string, IEnumerable<ISessionInfo>> R = new Dictionary<string, IEnumerable<ISessionInfo>>();

                foreach(var s in this.Sessions) {
                    string PrjNmn = s.ProjectName;
                    if(PrjNmn == null || PrjNmn.Length <= 0)
                        PrjNmn = "__unknown_project__";

                    IEnumerable<ISessionInfo> ProjectSessions;
                    if(!R.TryGetValue(PrjNmn, out ProjectSessions)) {
                        ProjectSessions = new List<ISessionInfo>();
                        R.Add(PrjNmn, ProjectSessions);
                    }
                    Debug.Assert(ProjectSessions != null);

                    ((List<ISessionInfo>)ProjectSessions).Add(s);
                }

                
                return R;
            }
        }
    }
}
