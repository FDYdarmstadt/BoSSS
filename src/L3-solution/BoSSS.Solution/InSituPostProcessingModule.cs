using BoSSS.Foundation.XDG;
using ilPSP;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {
    
    /// <summary>
    /// Base-class for in-situ post-processing, i.e. post-processing task 
    /// which are performed while the simulation is running.
    /// Objects of this type can be added to the list <see cref="Control.AppControl.PostprocessingModules"/>,
    /// and the respective post-processing (implemented in the derived class) will be executed.
    /// </summary>
    [Serializable]
    public abstract class InSituPostProcessingModule : IDisposable { 

    

        /// <summary>
        /// Name of the logfile; 
        /// </summary>
        [JsonIgnore]
        protected abstract string LogFileName {
            get;
        }

        /// <summary>
        /// if the log is a CSV-table, the character which should be used for column separation
        /// </summary>
        virtual protected char ColumnSeperator => '\t';



        /// <summary>
        /// Override to write a header for the logging file
        /// </summary>
        virtual protected void WriteHeader(TextWriter textWriter) {

        }

        [JsonIgnore]
        [NonSerialized]
        bool m_Setupdone = false;

        /// <summary>
        /// 
        /// </summary>
        virtual public void Setup(IApplication solverMain) {
            
            if(m_Setupdone)
                throw new NotSupportedException("Setup Routine may only be called once.");
            m_Setupdone = true;

            this.SolverMain = solverMain;
            if (solverMain.MPIRank == 0) {
                if (LogFileName != null) {
                    if (solverMain.CurrentSessionInfo.ID != null && !solverMain.CurrentSessionInfo.ID.Equals(Guid.Empty)) {
                        if (solverMain.MPIRank == 0) {
                            this.Log = SolverMain.DatabaseDriver.FsDriver.GetNewLog(LogFileName, solverMain.CurrentSessionInfo.ID);
                        } else {
                            this.Log = new StreamWriter(Stream.Null);
                        }
                    } else {
                        this.Log = new StreamWriter(
                            new FileStream(LogFileName + ".txt", FileMode.Create, FileAccess.Write, FileShare.Read));
                    }
                }
                WriteHeader(this.Log);
            }
        }

        /// <summary>
        /// reference to solver application class
        /// </summary>
        [JsonIgnore]
        protected IApplication SolverMain {
            get;
            private set;
        }

        /// <summary>
        /// reference to application control
        /// </summary>
        [JsonIgnore]
        protected Control.AppControl Control {
            get {
                return SolverMain.ControlBase;
            }
        }

        /// <summary>
        /// helper for casting <see cref="Control"/>
        /// </summary>
        protected void GetControl<C>(out C catControl) where C : Control.AppControl {
            catControl = (C)this.Control;
        }


        /// <summary>
        /// Access to the level set tracker 
        /// </summary>
        [JsonIgnore]
        protected LevelSetTracker LsTrk {
            get {
                return SolverMain.LsTrk;
            }
        }

        /// <summary>
        /// Main implementation point for the post-processing routine after each time-step
        /// </summary>
        abstract protected void PerformTimestepPostProcessing(int iTimestep, double PhysTime);

        /// <summary>
        /// Driver routine for the application to call the post-processing
        /// </summary>
        public void DriverTimestepPostProcessing(int iTimestep, double PhysTime) {
            if(iTimestep % LogPeriod == 0) {
                PerformTimestepPostProcessing(iTimestep, PhysTime);

                if(Log != null && ColumnCounter > 0) {
                    Log.WriteLine();
                    Log.Flush();
                    ColumnCounter = 0;
                }
            }
        }


        /// <summary>
        /// finalization
        /// </summary>
        public void Dispose() {
            if(Log != null) {
                Log.Flush();
                Log.Close();
                Log.Dispose();
                Log = null;
            }
        }

        int ColumnCounter = 0;

        /// <summary>
        /// writes an integer number to the current log
        /// </summary>
        protected void AppendToLog(int I) {
            if (this.SolverMain.MPIRank != 0) return;
            if (ColumnCounter > 0)
                Log.Write(ColumnSeperator);
            Log.Write(I);
            ColumnCounter++;
        }

        /// <summary>
        /// writes a string to the current log
        /// </summary>
        protected void AppendToLog(string s) {
            if (this.SolverMain.MPIRank != 0) return;
            if (ColumnCounter > 0)
                Log.Write(ColumnSeperator);
            Log.Write(s);
            ColumnCounter++;
        }

        /// <summary>
        /// writes a string line to the current log
        /// </summary>
        protected void AppendLineToLog(string s) {
            if (this.SolverMain.MPIRank != 0) return;
              Log.Write(s + "\n");
            ColumnCounter++;
        }

        /// <summary>
        /// appends a floating-point number to the current log
        /// </summary>
        protected void AppendToLog(double d) {
            if (this.SolverMain.MPIRank != 0) return;
            if (ColumnCounter > 0)
                Log.Write(ColumnSeperator);
            Log.Write(d.ToStringDot());
            ColumnCounter++;
        }

        /// <summary>
        /// writes an entire collection of floating-point values to the log
        /// </summary>
        protected void AppendToLog(IEnumerable<double> values) {
            if (this.SolverMain.MPIRank != 0) return;
            foreach (var d in values)
                AppendToLog(d);
        }

        /*
        /// <summary>
        /// writes an entire collection of values to the log
        /// </summary>
        protected void AppendToLog(System.Runtime.CompilerServices.ITuple t) {
            for(int i = 0; i < t.Length; i++) {
                object o = t[i];
                if(o is double d) {
                    AppendToLog(d);
                } else if(o is int k) {
                    AppendToLog(k);
                } else {
                    AppendToLog(o.ToString());
                }
            }
        }
        */

        /// <summary>
        /// Logfile 
        /// </summary>
        [JsonIgnore]
        protected TextWriter Log {
            get;
            private set;
        }


        /// <summary>
        /// <see cref="Application{T}.QueryResultTable"/>
        /// </summary>
        protected QueryResultTable QueryResultTable {
            get {
                return SolverMain.QueryResultTable;
            }
        }


        /// <summary>
        /// <see cref="Application{T}.QueryResultTable"/>
        /// </summary>
        protected Queries.QueryHandler QueryHandler {
            get {
                return SolverMain.QueryHandler;
            }
        }

        /// <summary>
        /// Period between two logged timesteps (1 is logging every timetesp, 2 is every second...)
        /// </summary>
        [DataMember]
        public int LogPeriod = 1;

        /// <summary>
        /// Pre-solver or post-solver (0: both, 1: only pre-solver and 2: only post-solver )
        /// </summary>
        [DataMember]
        public int SolverStage = 0;
    }
}
