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
//#define TEST
using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using log4net;
using MPI.Wrappers;

namespace ilPSP.Tracing {

    /// <summary>
    /// This module contains methods to log trace information to the trace files. 
    /// </summary>
    static public class Tracer {

        /// <summary>
        /// a list of all name-spaces for which <see cref="FuncTrace"/> should perform tracing/logging;
        /// </summary>
        internal static string[] m_NamespacesToLog = new string[0];

        /// <summary>
        /// a list of all name-spaces for which <see cref="FuncTrace"/> should perform tracing/logging;
        /// </summary>
        public static string[] NamespacesToLog {
            get {
                return ((string[])(m_NamespacesToLog.Clone()));
            }
            set {
                //Console.Write("Resetting logging namespaces: ");
                //if(value == null || value.Length <= 0) {
                //     dbg_launch();
                //    Console.WriteLine("NIX2LOG.");
                //} else {
                //    foreach(string s in value)
                //        Console.Write($"<{s}> ");
                //    Console.WriteLine();
                //}
                var NameSpaceList = value;
                if (NameSpaceList == null)
                    throw new ArgumentNullException();
                m_NamespacesToLog = NameSpaceList;
            }
        }


        /// <summary>
        /// overrides <see cref="NamespacesToLog"/> to trace everything; 
        /// This is intended to be used only by the test runner
        /// </summary>
        public static bool NamespacesToLog_EverythingOverrideTestRunner = false;


        /// <summary>
        /// Must be initialized to write memory-tracing
        /// </summary>
        public static TextWriter MemtraceFile;// = new System.IO.StreamWriter("memory." + ilPSP.Environment.MPIEnv.MPI_Rank + "of" + ilPSP.Environment.MPIEnv.MPI_Size + ".csv");


        internal static System.Collections.Generic.LinkedList<string> MemtraceFileTemp = new System.Collections.Generic.LinkedList<string>();
        
        


        /// <summary>
        /// 
        /// </summary>
        public static MemoryInstrumentationLevel MemoryInstrumentationLevel {
            get;
            set;
        } = MemoryInstrumentationLevel.None;


        /// <summary>
        /// Explicit switch for turning cone instrumentation on/off; this is useful when to much overhead is caused by instrumentation.
        /// </summary>
        public static bool InstrumentationSwitch {
            get {
                return m_InstrumentationSwitch;
            }
            set {
                m_InstrumentationSwitch = value;
            }
        }

        static bool m_InstrumentationSwitch = true;

        static Tracer() {
            _Root = new MethodCallRecord(null, "root_frame");
            Current = _Root;
            TotalTime = new Stopwatch();
            TotalTime.Reset();
            TotalTime.Start();
        }

        /// <summary>
        /// the root of the call tree; the runtime (see <see cref="MethodCallRecord.m_TicksSpentInMethod"/>) of the root object is equal to
        /// the overall time spend in the application so far.
        /// </summary>
        static public MethodCallRecord Root {
            get {
                TotalTime.Stop();
                _Root.m_TicksSpentInMethod = TotalTime.Elapsed.Ticks;
                TotalTime.Start();
//#if TEST
//                Console.WriteLine("memory measuring activated. Use this only for Debugging / Testing. This will have an impact on performance.");
//#endif
                return _Root;
            }
        }


        static private Stopwatch TotalTime;

        static private MethodCallRecord _Root;

        /// <summary>
        /// The record corresponding to the current function.
        /// </summary>
        public static MethodCallRecord Current {
            get;
            private set;
        }

        static private long GetMPITicks() {
            return ((MPI.Wrappers.IMPIdriver_wTimeTracer)MPI.Wrappers.csMPI.Raw).TicksSpent;
        }

        

        private static readonly object padlock = new object();


        internal static int Push_MethodCallRecord(string _name, out MethodCallRecord mcr) {
            Debug.Assert(InstrumentationSwitch == true);
            

            //if (Tracer.Current != null) {
            lock(padlock) {
                if(!Tracer.Current.Calls.TryGetValue(_name, out mcr)) {
                    mcr = new MethodCallRecord(Tracer.Current, _name);
                    Tracer.Current.Calls.Add(_name, mcr);
                }
            }
            Tracer.Current = mcr;
            mcr.CallCount++;
            mcr.m_TicksSpentinBlocking = -GetMPITicks();
            //mcr.m_Memory = -GetMemory();
            //} else {
            //    Debug.Assert(Tracer.Root == null);
            //    var mcr = new MethodCallRecord(Tracer.Current, _name);
            //    Tracer.Root = mcr;
            //    Tracer.Current = mcr;
            //}

            return mcr.Depth;
        }

        internal static int Pop_MethodCallrecord(long ElapsedTicks, long Memory_increase, long PeakMemory_increase, out MethodCallRecord mcr) {
            Debug.Assert(InstrumentationSwitch == true, "instrumentation switch off!");


            Debug.Assert(!object.ReferenceEquals(Current, _Root), "root frame cannot be popped");
            //if(!object.ReferenceEquals(Current, _Root) == false) {
            //    Console.Error.WriteLine("root frame cannot be popped");
            //    throw new Exception("root frame cannot be popped");
            //}
            Tracer.Current.m_TicksSpentInMethod += ElapsedTicks;
            Tracer.Current.m_TicksSpentinBlocking += GetMPITicks();
            Tracer.Current.m_MemoryIncrease = Math.Max(Tracer.Current.m_MemoryIncrease, Memory_increase);
            Tracer.Current.m_PeakMemoryIncrease = Math.Max(Tracer.Current.m_PeakMemoryIncrease, PeakMemory_increase);

            //fails for some reason on lichtenberg:
            //Debug.Assert(ElapsedTicks > Tracer.Current.m_TicksSpentinBlocking, $"ticks are fucked up: elapsed = {ElapsedTicks}, blocking = {Tracer.Current.m_TicksSpentinBlocking}");

            mcr = Tracer.Current;
            Tracer.Current = Tracer.Current.ParrentCall;
            return Tracer.Current.Depth;
        }


        internal static MethodCallRecord LogDummyblock(long ticks, string _name) {
            Debug.Assert(InstrumentationSwitch == true);

            MethodCallRecord mcr;
            if (!Tracer.Current.Calls.TryGetValue(_name, out mcr)) {
                mcr = new MethodCallRecord(Tracer.Current, _name);
                //mcr.IgnoreForExclusive = true;
                Tracer.Current.Calls.Add(_name, mcr);
            }
            mcr.CallCount++;
            //Debug.Assert(mcr.IgnoreForExclusive == true);
            mcr.m_TicksSpentInMethod += ticks;

            return mcr;
        }
    }

    /// <summary>
    /// Levels of memory instrumentation during func tracing
    /// </summary>
    public enum MemoryInstrumentationLevel {

        /// <summary>
        /// Off
        /// </summary>
        None = 0,

        /// <summary>
        /// only memory managed by garbage collector, causing almost no overhead
        /// </summary>
        OnlyGcTotalMemory = 1,


        /// <summary>
        /// Garbage collector memory and un-manged private memory;
        /// Very expensive instrumentation option, slows down the application by a factor of two to three!!!
        /// </summary>
        GcAndPrivateMemory = 2

    }


    /// <summary>
    /// baseclass for runtime measurement of functions and blocks
    /// </summary>
    /// <example>
    /// void SomeFunction() {
    ///     using(var tr = new FuncTrace) {
    ///         // method-body
    ///     }
    /// }
    /// </example>
    abstract public class Tmeas : IDisposable {

        public static int Memtrace_lineCount = 0;

        static long PreviousLineMem = 0;


        /// <summary>
        /// Hack to write <see cref="Info"/> also to cout.
        /// </summary>
        public bool InfoToConsole = false;

        /// <summary>
        /// Enables the stdout output for all ranks
        /// (normally, suppressed for all ranks but 0)
        /// </summary>
        public void StdoutOnAllRanks() {
            m_StdoutOnlyOnRank0Bkup = ilPSP.Environment.StdoutOnlyOnRank0;
            ilPSP.Environment.StdoutOnlyOnRank0 = true;
        }

        /// <summary>
        /// logger to write the enter/leave -- messages to;
        /// </summary>
        internal protected ILog m_Logger = null;

        /// <summary>
        /// logger to write the enter/leave -- messages to;
        /// </summary>
        virtual public ILog Logger {
            get {
                return m_Logger;
            }
        }

      

        /// <summary>
        /// current process to trace memory
        /// </summary>
        static protected long GetPrivateMemory() {
            if(!Tracer.InstrumentationSwitch)
                return 0;

            
            // expensive: 
            switch (Tracer.MemoryInstrumentationLevel) {
                case MemoryInstrumentationLevel.GcAndPrivateMemory:
                    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    Console.WriteLine("WARNING: `Tracer.MemoryInstrumentationLevel` set to " + MemoryInstrumentationLevel.GcAndPrivateMemory);
                    Console.WriteLine("This gives the most accurate memory allocation report, but it is a");
                    Console.WriteLine("very expensive instrumentation option, slows down the application by a factor of two to three!!!");
                    Console.WriteLine("Should not be used for production runs, but in order to identify memory problems.");
                    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    var p = Process.GetCurrentProcess(); // process object must be fresh, otherwise old data
                    return p.WorkingSet64;

                case MemoryInstrumentationLevel.OnlyGcTotalMemory:
                    long virt = GC.GetTotalMemory(false);
                    return virt;

                case MemoryInstrumentationLevel.None:
                    return 0;

                default:
                    throw new NotImplementedException();
                    //var p = Process.GetCurrentProcess();
                    //return p.WorkingSet64;
            }

        }

        /// <summary>
        /// Currently allocated memory in bytes;
        /// </summary>
        /// <returns></returns>
        public long GetMemory() {
            return GetPrivateMemory();
        }

        /// <summary>
        /// Currently allocated memory Megabytes
        /// </summary>
        /// <returns></returns>
        public double GetMemoryMB() {
            return (double)GetPrivateMemory() /(1024.0 * 1024.0);
        }
                       

        /// <summary>
        /// the internal stopwatch
        /// </summary>
        protected Stopwatch Watch;

        long PeakWorkingSet_onEntry;
        long WorkingSet_onEntry;

        long PeakWorkingSet_onExit;
        long WorkingSet_onExit;

        /// <summary>
        /// ctor
        /// </summary>
        protected Tmeas() {
            //startTicks = Watch.ElapsedTicks;
            Watch = new Stopwatch();
            Watch.Start();

            {
                WorkingSet_onEntry = GetPrivateMemory();
                PeakWorkingSet_onEntry = 0;
            }

            
        }

        bool? m_StdoutOnlyOnRank0Bkup;

        /// <summary>
        /// logs an 'inclusive' block;
        /// </summary>
        public MethodCallRecord LogDummyblock(long ticks, string name) {
            if(!Tracer.InstrumentationSwitch)
                return new MethodCallRecord(null, "dummy");
            else 
                return Tracer.LogDummyblock(ticks, name);
        }


        /// <summary>
        /// time elapsed after stopping.
        /// </summary>
        public TimeSpan Duration {
            get;
            private set;
        }


        /// <summary>
        /// Increase of the working set memory during this method
        /// </summary>
        public long AllocatedMem {
            get {
                return Math.Max(0, WorkingSet_onExit - WorkingSet_onEntry);
            }
        }

        /// <summary>
        /// Detected increase of peak memory during execution
        /// </summary>
        public long PeakMem {
            get {
                return Math.Max(0, Math.Max(PeakWorkingSet_onExit - PeakWorkingSet_onEntry, AllocatedMem));
            }
        }

        /// <summary>
        /// (selective) Info - message
        /// </summary>
        /// <param name="o"></param>
        public void Info(object o) {
            if (DoLogging)
                m_Logger.Info(o);

            if(InfoToConsole)
                Console.WriteLine(o);
        }


        /// <summary>
        /// (selective) error - message
        /// </summary>
        /// <param name="o"></param>
        public void Error(object o) {
            if (DoLogging)
                m_Logger.Error(o);

            Console.Error.WriteLine(o);
        }

        /// <summary>
        /// (selective) Warning - message
        /// </summary>
        /// <param name="o"></param>
        public void Warning(object o) {
            if (DoLogging)
                m_Logger.Warn(o);

            Console.WriteLine(o);
        }



        /// <summary>
        /// writes information about system memory usage to trace file;
        /// This seems to have a severe performance impact on server OS, therefore deactivated (fk,21dec20)
        /// </summary>
        public void LogMemoryStat() {

            if(!Tracer.InstrumentationSwitch)
                return;

            Process myself = Process.GetCurrentProcess();

            {
                string s = "MEMORY STAT.: garbage collector memory: ";
                try {
                    long virt = GC.GetTotalMemory(false) / (1024 * 1024);
                    s += (virt + " Meg");
                } catch(Exception e) {
                    s += e.GetType().Name + ": " + e.Message;
                }
                Info(s);
            }

            {
                string s = "MEMORY STAT.: working set memory: ";
                try {
                    long virt = myself.WorkingSet64 / (1024 * 1024);
                    s += (virt + " Meg");
                } catch(Exception e) {
                    s += e.GetType().Name + ": " + e.Message;
                }
                Info(s);
            }
            {
                string s = "MEMORY STAT.: peak working set memory: ";
                try {
                    long virt = myself.PeakWorkingSet64 / (1024 * 1024);
                    s += (virt + " Meg");
                } catch(Exception e) {
                    s += e.GetType().Name + ": " + e.Message;
                }
                Info(s);
            }

            /* also not very interesting: 
            {
                string s = "MEMORY STAT.: private memory: ";
                try {
                    long virt = myself.PrivateMemorySize64 / (1024 * 1024);
                    s += (virt + " Meg");
                } catch(Exception e) {
                    s += e.GetType().Name + ": " + e.Message;
                }
                Console.WriteLine(s);
                Info(s);
            }
            

            /*
             * virtual memory seems to be useless information; 
             * some incredibly high value
            {
                string s = "MEMORY STAT.: peak virtual memory: ";
                try {
                    long virt = myself.PeakVirtualMemorySize64 / (1024 * 1024);
                    s += (virt + " Meg");
                } catch(Exception e) {
                    s += e.GetType().Name + ": " + e.Message;
                }
                Console.WriteLine(s);
                Info(s);
            }
            {
                string s = "MEMORY STAT.: virtual memory: ";
                try {
                    long virt = myself.VirtualMemorySize64 / (1024 * 1024);
                    s += (virt + " Meg");
                } catch(Exception e) {
                    s += e.GetType().Name + ": " + e.Message;
                }
                Console.WriteLine(s);
                Info(s);
            }
            */
        }

        static internal ILog s_FtLogger = LogManager.GetLogger(typeof(FuncTrace));

        
        bool m_DoLogging;

        /// <summary>
        /// true, if the timing measurement from this object is logged;
        /// this behavior is controlled by <see cref="Tracer.NamespacesToLog"/>
        /// </summary>
        public bool DoLogging {
            get {
                return m_DoLogging;
            }
            protected set {
                m_DoLogging = value;
            }
        }

        /// <summary>
        /// name of the time tracing block
        /// </summary>
        protected string _name;

        /// <summary>
        /// Message when the measurement starts.
        /// </summary>
        protected void EnterMessage(string elo, string _name) {
            int newDepth = Tracer.Push_MethodCallRecord(_name, out var mcr);
            this._name = _name;

            if (DoLogging) {
                s_FtLogger.Info(elo + _name + " new stack depth = " + newDepth);
            }

            
            LogMemtrace(WorkingSet_onEntry, ">", mcr);
        }


        /// <summary>
        /// Message when the measurement is finished.
        /// </summary>
        protected void LeaveLog() {
            if(!Tracer.InstrumentationSwitch)
                return;

            int newDepht = Tracer.Pop_MethodCallrecord(this.Duration.Ticks, this.AllocatedMem, this.PeakMem, out var mcr);

            if (DoLogging) {
                
                string time = this.Duration.TotalSeconds.ToString(NumberFormatInfo.InvariantInfo);
                string str = string.Format("LEAVING {0} ({1} sec, return to stack depth = {2})", _name, time, newDepht);

                try {
                    s_FtLogger.Info(str);
                } catch (Exception nre) {
                    Console.Error.WriteLine("ERRROR (logging): " + nre.Message);
                    Console.Error.WriteLine(nre.StackTrace);
                }
            }

            LogMemtrace(WorkingSet_onExit, "<", mcr);
        }


        void LogMemtrace(long WorkingSet, string Ch, MethodCallRecord mcr) {
            string PrintMeg(long bytes) {
                double megs = (double)bytes / (1024.0 * 1024.0);
                return Math.Round(megs).ToString();
            }


            if (Tracer.MemoryInstrumentationLevel != MemoryInstrumentationLevel.None) {

                string Line;
                using (var Memtrace = new StringWriter()) {


                    Memtrace.Write(Memtrace_lineCount); // col 0: line number
                    Memtrace_lineCount++;
                    Memtrace.Write("\t");
                    Memtrace.Write(WorkingSet); // col 1: mem in Bytes
                    long diff = WorkingSet - PreviousLineMem;
                    PreviousLineMem = WorkingSet;
                    Memtrace.Write("\t");
                    Memtrace.Write($"{PrintMeg(WorkingSet)}"); // col 2: mem in megs
                    Memtrace.Write("\t");
                    Memtrace.Write($"{diff}"); // col 3: diff-mem in bytes
                    Memtrace.Write("\t");
                    Memtrace.Write($"{PrintMeg(diff)}"); // col 4: diff-mem in megs
                    Memtrace.Write("\t");

                    /*
                    long gcTot = WorkingSet;//.MPISum();
                    long gcMax = WorkingSet;//.MPIMax();
                    long gcMin = WorkingSet;//.MPIMin();

                    Memtrace.Write($"{gcTot}");
                    Memtrace.Write("\t");
                    Memtrace.Write($"{gcMax}");
                    Memtrace.Write("\t");
                    Memtrace.Write($"{gcMin}");
                    Memtrace.Write("\t");


                    long mpiDiff = gcTot - PreviousLineMem_mpi;
                    PreviousLineMem_mpi = gcTot;
                    Memtrace.Write($"{PrintMeg(mpiDiff)}");
                    Memtrace.Write("\t");
                    */

                    using (var stw = new StringWriter()) {
                        stw.Write(Ch);
                        var _mcr = mcr;
                        while (_mcr != null) {
                            stw.Write(_mcr.Name);


                            _mcr = _mcr.ParrentCall;
                            if (_mcr != null)
                                stw.Write(">");
                        }

                        Memtrace.Write(stw.ToString());
                    }

                    Line = Memtrace.ToString();
                }



                if (Tracer.MemtraceFile != null) {
                    while (Tracer.MemtraceFileTemp.Count > 0) {
                            Tracer.MemtraceFile.WriteLine(Tracer.MemtraceFileTemp.First.Value);
                            Tracer.MemtraceFileTemp.RemoveFirst();
                    }
                    Tracer.MemtraceFile.WriteLine(Line);
                } else {
                    Tracer.MemtraceFileTemp.AddLast(Line);
                    while (Tracer.MemtraceFileTemp.Count > 500) {
                        Tracer.MemtraceFileTemp.RemoveFirst();
                    }
                }
            }
        }


#region IDisposable Members

        
        /// <summary>
        /// stops the measurement
        /// </summary>
        virtual public void Dispose() {
            //this.DurationTicks = Watch.ElapsedTicks - startTicks;
            Watch.Stop();
            this.Duration = Watch.Elapsed;
            {
                WorkingSet_onExit = GetPrivateMemory();
                PeakWorkingSet_onExit = 0;
            }

            if (m_StdoutOnlyOnRank0Bkup != null)
                ilPSP.Environment.StdoutOnlyOnRank0 = m_StdoutOnlyOnRank0Bkup.Value;

            LeaveLog();
        }

#endregion
    }

    /// <summary>
    /// measures and logs the runtime of a function
    /// </summary>
    public class FuncTrace : Tmeas {


        /// <summary>
        /// ctor: logs the 'enter' - message
        /// </summary>
        public FuncTrace() : base() {
            if(!Tracer.InstrumentationSwitch)
                return;

            string _name;
            Type callingType;
            {
                StackFrame fr = new StackFrame(1, true);

                System.Reflection.MethodBase m = fr.GetMethod();
                _name = m.DeclaringType.FullName + "." + m.Name;
                callingType = m.DeclaringType;
            }

            if (Tracer.NamespacesToLog_EverythingOverrideTestRunner) {
                DoLogging = true;
            } else {
                for (int i = Tracer.m_NamespacesToLog.Length - 1; i >= 0; i--) {
                    if (_name.StartsWith(Tracer.m_NamespacesToLog[i])) {
                        DoLogging = true;
                        break;
                    }
                }
            }

            m_Logger = LogManager.GetLogger(callingType);
            base.EnterMessage("ENTERING ", _name);
        }

        /// <summary>
        /// ctor: logs the 'enter' - message
        /// </summary>
        public FuncTrace(string UserName) : base() {
            if(!Tracer.InstrumentationSwitch)
                return;

            string _name = UserName;

            Type callingType = null;
            string filtername;
            {
                StackFrame fr = new StackFrame(1, true);

                System.Reflection.MethodBase m = fr.GetMethod();
                callingType = m.DeclaringType;
                filtername = callingType.FullName;
            }

            if (Tracer.NamespacesToLog_EverythingOverrideTestRunner) {
                DoLogging = true;
            } else {
                for (int i = Tracer.m_NamespacesToLog.Length - 1; i >= 0; i--) {
                    if (filtername.StartsWith(Tracer.m_NamespacesToLog[i])) {
                        DoLogging = true;
                        break;
                    }
                }
            }

            m_Logger = LogManager.GetLogger(callingType);
            base.EnterMessage("ENTERING ", _name);
        }


        /// <summary>
        /// dtor: logs the 'leave' - message
        /// </summary>
        public override void Dispose() {
            base.Dispose();
        }

    

     

        

        /*
        /// <summary>
        /// logs an 'inclusive' block (see <see cref="MethodCallRecord.IgnoreForExclusive"/> );
        /// </summary>
        public MethodCallRecord LogDummyblock(long ticks, string name) {
            if(!Tracer.InstrumentationSwitch)
                return new MethodCallRecord(null, "dummy");
            else 
                return Tracer.LogDummyblock(ticks, name);
        }
        */
    }

    /// <summary>
    /// measures and logs the runtime of some 
    /// </summary>
    public class BlockTrace : Tmeas {

        readonly FuncTrace _f;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="Title">
        /// a title under which the block appears in the logfile
        /// </param>
        /// <param name="f">
        /// tracing of function which contains the block
        /// </param>
        /// <param name="timeToCout">
        /// if true, the time elapsed in the block will be printed to console 
        /// </param>
        public BlockTrace(string Title, FuncTrace f, bool timeToCout = false) {
            if(!Tracer.InstrumentationSwitch)
                return;
            base.DoLogging = f.DoLogging;
            string _name = Title;
            _f = f;
            m_Logger = _f.m_Logger;
            m_timeToCout = timeToCout;
            base.EnterMessage("BLKENTER ", _name);
        }

        readonly bool m_timeToCout = false;

        /*
        /// <summary>
        /// stops the watch
        /// </summary>
        public override void Dispose() {
            base.Dispose();
            if(!Tracer.InstrumentationSwitch)
                return;
            int newDepth = Tracer.Pop_MethodCallrecord(base.Duration.Ticks, this.AllocatedMem, this.PeakMem);

            if (_f.DoLogging) {
                FuncTrace.s_FtLogger.Info("LEAVING " + _name + " ("
                    + base.Duration.TotalSeconds.ToString(NumberFormatInfo.InvariantInfo)
                    + " sec, return to stack depth = " + newDepth + ")");
            }
        }
        */

        /// <summary>
        /// 
        /// </summary>
        public override void Dispose() {
            base.Dispose();

            if(m_timeToCout) {
                Console.WriteLine(base._name + ": " + base.Watch.Elapsed.TotalSeconds + " sec.");
            }
        }
    }
}
