using BoSSS.Foundation.XDG;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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
    //public abstract class InSituPostProcessingModule<S,C> : IDisposable 
    //    where C : Control.AppControl, new()
    //    where S : Application<C> //
    

        /// <summary>
        /// Name of the logfile; 
        /// </summary>
        [JsonIgnore]
        protected abstract string LogFileName {
            get;
        }



        /// <summary>
        /// Override to write a header for the logging file
        /// </summary>
        virtual protected void WriteHeader(TextWriter textWriter) {

        }

        /// <summary>
        /// 
        /// </summary>
        public void Setup(IApplication solverMain) {
            this.SolverMain = solverMain;
            if(LogFileName != null) {
                if(solverMain.CurrentSessionInfo.ID != null && !solverMain.CurrentSessionInfo.ID.Equals(Guid.Empty)) {
                    this.Log = SolverMain.DatabaseDriver.FsDriver.GetNewLog(LogFileName, solverMain.CurrentSessionInfo.ID);
                } else {
                    this.Log = new StreamWriter(LogFileName + ".txt");
                }
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
        /// Main implementation point for the post-processing routine
        /// </summary>
        abstract public void PerformPostProcessing(int iTimestep, double PhysTime, double dt);


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


        /// <summary>
        /// Logfile 
        /// </summary>
        [JsonIgnore]
        protected TextWriter Log {
            get;
            private set;
        }




    }
}
