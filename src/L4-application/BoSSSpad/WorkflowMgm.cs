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
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP.Utils;
using System.Data;
using System.Reflection;
using System.Threading;
using ilPSP;
using BoSSS.Solution.Control;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Workflow management.
    /// </summary>
    public partial class WorkflowMgm {

        /// <summary>
        /// Not intended for user interaction.
        /// </summary>
        internal WorkflowMgm() { }

        string m_CurrentProject;

        /// <summary>
        /// Name of the current project;
        /// </summary>
        public string CurrentProject {
            get {
                return m_CurrentProject;
            }
        }

        /// <summary>
        /// Clears/Invalidates all cached data.
        /// </summary>
        public void InvalidateCaches() {
            m_Sessions = null;
            m_Grids = null;
        }


        TimeSpan m_UpdatePeriod = new TimeSpan(0, 5, 0);

        /// <summary>
        /// Data like <see cref="Sessions"/> is cached for performance reasons; after this time span is elapsed, the data is re-read from disk.
        /// </summary>
        public TimeSpan UpdatePeriod {
            get {
                return m_UpdatePeriod;
            }
            set {
                m_UpdatePeriod = value;
            }
        }


        /// <summary>
        /// Defines the name of the current project;
        /// </summary>
        public void Init(string ProjectName) {
            if ((m_CurrentProject == null) || (!m_CurrentProject.Equals(ProjectName)))
                InvalidateCaches();
            m_CurrentProject = ProjectName;
            Console.WriteLine("Project name is set to '{0}'.", ProjectName);
        }

        DateTime m_Sessions_CacheTime;
        ISessionInfo[] m_Sessions;

        /// <summary>
        /// A list of all sessions in the current project.
        /// </summary>
        public ISessionInfo[] Sessions {
            get {
                if (CurrentProject.IsEmptyOrWhite()) {
                    Console.WriteLine("Workflow management not initialized yet - call Init(...)!");
                    return new ISessionInfo[0];
                }

                if (m_Sessions == null || ((DateTime.Now - m_Sessions_CacheTime) > UpdatePeriod)) {

                    List<ISessionInfo> ret = new List<ISessionInfo>();

                    if (InteractiveShell.databases != null) {
                        foreach (var db in InteractiveShell.databases) {
                            var SS = db.Sessions.Where(delegate( ISessionInfo si) {
#if DEBUG 
                                return si.ProjectName.Equals(this.CurrentProject);
#else
                                Guid g = Guid.Empty;
                                try {
                                    g = si.ID;
                                    return si.ProjectName.Equals(this.CurrentProject);
                                } catch(Exception e) {
                                    Console.WriteLine("Warning: " + e.Message + " reading session " + g + ".");
                                    return false;
                                }
#endif
                            });
                            ret.AddRange(SS);
                        }
                    }

                    m_Sessions = ret.ToArray();
                    m_Sessions_CacheTime = DateTime.Now;
                }

                return m_Sessions;
            }
        }

        /// <summary>
        /// A list of all tags in all sessions.
        /// </summary>
        public string[] Tags {
            get {
                if (CurrentProject.IsEmptyOrWhite()) {
                    Console.WriteLine("Workflow management not initialized yet - call Init(...)!");
                    return new string[0];
                }

                HashSet<string> r = new HashSet<string>();

                foreach (var s in this.Sessions) {
                    r.AddRange(s.Tags);
                }

                return r.ToArray();
            }
        }

        IGridInfo[] m_Grids;
        DateTime m_Grids_CacheTime;

        /// <summary>
        /// A list of all grids which are used in the current project.
        /// </summary>
        public IGridInfo[] Grids {
            get {
                if (CurrentProject.IsEmptyOrWhite()) {
                    Console.WriteLine("Workflow management not initialized yet - call Init(...)!");
                    return new IGridInfo[0];
                }

                if (m_Grids == null || ((DateTime.Now - m_Sessions_CacheTime) > UpdatePeriod)) {
                    HashSet<IGridInfo> grids = new HashSet<IGridInfo>(
                    new ilPSP.FuncEqualityComparer<IGridInfo>((g1, g2) => g1.ID.Equals(g2.ID), g => g.ID.GetHashCode()));

                    foreach (var s in this.Sessions) {
                        Console.Write("Session " + s.ID + " ... ");
                        grids.AddRange(s.GetGrids());
                        Console.WriteLine(" done.");
                    }

                    m_Grids = grids.ToArray();
                    m_Grids_CacheTime = DateTime.Now;
                }

                return m_Grids;
            }
        }

        /// <summary>
        /// The keys and queries <see cref="ISessionInfo.KeysAndQueries"/> of all sessions in the 
        /// project (see <see cref="Sessions"/>) in one table.
        /// </summary>
        public DataTable SessionTable {
            get {
                var adiColi = AdditionalSessionTableColums.Select(kv => new Tuple<string, Func<ISessionInfo, object>>(kv.Key, kv.Value)).ToArray();
                return this.Sessions.GetSessionTable(adiColi);
            }
        }

        Dictionary<string, Func<ISessionInfo, object>> m_AdditionalSessionTableColums = new Dictionary<string, Func<ISessionInfo, object>>();

        /// <summary>
        /// Custom, user-defined columns for the session table (<see cref="SessionTable"/>).
        /// - keys: column name
        /// - values: functions which map the session info to a column value.
        /// </summary>
        public IDictionary<string, Func<ISessionInfo, object>> AdditionalSessionTableColums {
            get {
                return m_AdditionalSessionTableColums;
            }
        }

        



        Dictionary<string, Job> m_AllJobs = new Dictionary<string, Job>();

        /// <summary>
        /// Lists all compute jobs which are currently known by the work flow management system.
        /// - key: job name
        /// - item 
        /// </summary>
        public IDictionary<string, Job> AllJobs {
            get {
                return m_AllJobs;
            }
        }


        /// <summary>
        /// Blocks until all jobs in <see cref="AllJobs"/> ale either <see cref="JobStatus.Failed"/>
        /// or <see cref="JobStatus.FinishedSuccessful"/>.
        /// </summary>
        /// <param name="TimeOutSeconds">
        /// If positive, this method should terminate at latest after approximately this time period.
        /// </param>
        /// <param name="PollingIntervallSeconds">
        /// Seconds to wait before checking the jobs status again; should be in the order of seconds, not to overload the IO.
        /// </param>
        public void BlockUntilAllJobsTerminate(double TimeOutSeconds = -1, double PollingIntervallSeconds = 2) {
            DateTime start = DateTime.Now;
            while(true) {
                Thread.Sleep((int)PollingIntervallSeconds);

                if(TimeOutSeconds > 0) {
                    double RuntimeSoFar = (DateTime.Now - start).TotalSeconds;
                    if(RuntimeSoFar > TimeOutSeconds) {
                        Console.WriteLine("Timeout.");
                        return;
                    }
                }

                bool terminate = true;
                foreach(var J in this.AllJobs) {
                    var s = J.Value.Status;
                    if(s!= JobStatus.Failed && s != JobStatus.FinishedSuccessful) {
                        terminate = false;
                        break;
                    }
                }

                if (terminate) {
                    Console.WriteLine("All jobs finished.");
                    m_Sessions = null;
                    return;
                }

                Thread.Sleep(1000);
            }
            
        }

        List<Tuple<AppControl, int>> RegisteredControls = new List<Tuple<AppControl, int>>();


        /// <summary>
        /// Records the control object <paramref name="C"/> in an internal list, for its entire lifetime,
        /// and provides an index for it. 
        /// </summary>
        public int RegisterControl(AppControl C) {
            int max = 0;
            foreach (var t in RegisteredControls) {
                if (object.ReferenceEquals(t.Item1, C))
                    return t.Item2;
                max = Math.Max(t.Item2, max);
            }

            RegisteredControls.Add(Tuple.Create(C, max + 1));
            return max + 1;
        }


    }


    /*
    public static class MetaJobManager {

        static Dictionary<string, BatchProcessorClient> m_Computers;

        /// <summary>
        /// 
        /// </summary>
        static public Dictionary<string, BatchProcessorClient> Computers {
            get {
                if (m_Computers == null) {
                    m_Computers = new Dictionary<string, BatchProcessorClient>();
                }

                return m_Computers;
            }
        }

    }
    */
}



