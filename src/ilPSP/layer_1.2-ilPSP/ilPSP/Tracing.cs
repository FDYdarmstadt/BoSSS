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

namespace ilPSP.Tracing {

    /// <summary>
    /// This module contains methods to log trace information to the trace files. 
    /// </summary>
    static public class Tracer {

        //static ILog Logger = LogManager.GetLogger(typeof(Tracer));



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
                //    Debugger.Launch();
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
        /// Explicit switch for turning cone instrumentation on/off; this is useful when to much overhead is caused by instrumentation.
        /// </summary>
        public static bool InstrumentationSwitch = true;


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
#if TEST
                Console.WriteLine("memory measuring activated. Use this only for Debugging / Testing. This will have an impact on performance.");
#endif
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

        /*
        static private long GetMemory() {
            long mem = 0;
            Process myself = Process.GetCurrentProcess();
#if TEST
            {
                try {
                    //mem = myself.WorkingSet64 / (1024 * 1024);
                    mem = myself.PrivateMemorySize64 / (1024 * 1024);
                    //mem = GC.GetTotalMemory(false) / (1024 * 1024);
                } catch (Exception e) {
                    mem = 0;
                }
            }
#endif
            return mem;
        }
        */


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
            Debug.Assert(InstrumentationSwitch == true);

            Debug.Assert(!object.ReferenceEquals(Current, _Root), "root frame cannot be popped");
            Tracer.Current.m_TicksSpentInMethod += ElapsedTicks;
            Tracer.Current.m_TicksSpentinBlocking += GetMPITicks();
            Tracer.Current.m_MemoryIncrease = Math.Max(Tracer.Current.m_MemoryIncrease, Memory_increase);
            Tracer.Current.m_PeakMemoryIncrease = Math.Max(Tracer.Current.m_PeakMemoryIncrease, PeakMemory_increase);
            Debug.Assert(ElapsedTicks > Tracer.Current.m_TicksSpentinBlocking);
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
        /// Must be initialized to write memory-tracing
        /// </summary>
        public static TextWriter Memtrace;// = new System.IO.StreamWriter("memory_nocache." + ilPSP.Environment.MPIEnv.MPI_Rank + ".csv");

        public static int Memtrace_lineCount = 0;
        
        static long PreviousLineMem = 0;


        /// <summary>
        /// Very expensive instrumentation option, slows down by a factor of two to three!!!
        /// </summary>
        public static bool LogPrivateMem = true;

        /// <summary>
        /// current process to trace memory
        /// </summary>
        static protected long GetPrivateMemory() {


            // expensive: 
            if(LogPrivateMem) {
                var p = Process.GetCurrentProcess(); // process object must be fresh, otherwise old data
                //long q = p.PeakWorkingSet64;
                return p.WorkingSet64;
            } else {
                // almost no overhead, does not include un-managed memory.
                long virt = GC.GetTotalMemory(false);
                return virt;
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
                       

        Stopwatch Watch;

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
            m_Logger.Info(o);
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

        string _name;

        /// <summary>
        /// Message when the measurement starts.
        /// </summary>
        protected void EnterMessage(string elo, string _name) {
            int newDepth = Tracer.Push_MethodCallRecord(_name, out var mcr);
            this._name = _name;

            if (DoLogging) {
                s_FtLogger.Info(elo + _name + " new stack depth = " + newDepth);
            }

            /*if(Memtrace != null) {
                Memtrace.Write(Memtrace_lineCount); 
                Memtrace_lineCount++;
                Memtrace.Write("\t");
                Memtrace.Write(WorkingSet_onEntry);
                long diff = WorkingSet_onEntry - PreviousLineMem;
                PreviousLineMem = WorkingSet_onEntry;
                Memtrace.Write("\t");
                Memtrace.Write($"{Math.Round((double)WorkingSet_onEntry / (1024.0 * 1024.0))}");
                Memtrace.Write("\t");
                Memtrace.Write($"{Math.Round((double)diff / (1024.0 * 1024.0))}");
                Memtrace.Write("\t");
                Memtrace.Write(WorkingSet_onEntry);
                Memtrace.Write("\t");

                var _mcr = mcr;
                while(_mcr != null) {
                    Memtrace.Write(">");
                    Memtrace.Write(_mcr.Name);
                    _mcr = _mcr.ParrentCall;
                }
                Memtrace.WriteLine();
            }*/
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
            if(Memtrace != null) {
                Memtrace.Write(Memtrace_lineCount); 
                Memtrace_lineCount++;
                Memtrace.Write("\t");
                Memtrace.Write(WorkingSet);
                long diff = WorkingSet - PreviousLineMem;
                PreviousLineMem = WorkingSet;
                Memtrace.Write("\t");
                Memtrace.Write($"{Math.Round((double)WorkingSet / (1024.0 * 1024.0))}");
                Memtrace.Write("\t");
                Memtrace.Write($"{Math.Round((double)diff / (1024.0 * 1024.0))}");
                Memtrace.Write("\t");

                Memtrace.Write(Ch);
                var _mcr = mcr;
                while(_mcr != null) {
                    Memtrace.Write(_mcr.Name);
                    _mcr = _mcr.ParrentCall;
                    if(_mcr != null)
                        Memtrace.Write(">");
                }
                Memtrace.WriteLine();
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

            LeaveLog();
        }

#endregion
    }

    /// <summary>
    /// measures and logs the runtime of a function
    /// </summary>
    public class FuncTrace : Tmeas {

        

        /// <summary>
        /// Hack to write <see cref="Info"/> also to cout.
        /// </summary>
        public bool InfoToConsole = false;

        /// <summary>
        /// ctor: logs the 'enter' - message
        /// </summary>
        public FuncTrace() : base() {
            if(!Tracer.InstrumentationSwitch)
                return;

            string _name;
            Type callingType = null;
            {
                StackFrame fr = new StackFrame(1, true);

                System.Reflection.MethodBase m = fr.GetMethod();
                _name = m.DeclaringType.FullName + "." + m.Name;
                callingType = m.DeclaringType;
            }
            
            for (int i = Tracer.m_NamespacesToLog.Length - 1; i >= 0; i--) {
                if (_name.StartsWith(Tracer.m_NamespacesToLog[i])) {
                    DoLogging = true;
                    break;
                }
            }

            m_Logger = LogManager.GetLogger(callingType);
            base.EnterMessage("ENTERING ", _name);
        }

        // <summary>
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
            
            for (int i = Tracer.m_NamespacesToLog.Length - 1; i >= 0; i--) {
                if (filtername.StartsWith(Tracer.m_NamespacesToLog[i])) {
                    DoLogging = true;
                    break;
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

        /// <summary>
        /// (selective) Info - message
        /// </summary>
        /// <param name="o"></param>
        public void Info(object o) {
            if (base.DoLogging)
                m_Logger.Info(o);

            if(InfoToConsole)
                Console.WriteLine(o);
        }


        /// <summary>
        /// (selective) error - message
        /// </summary>
        /// <param name="o"></param>
        public void Error(object o) {
            if (base.DoLogging)
                m_Logger.Error(o);

            Console.Error.WriteLine(o);
        }

        /// <summary>
        /// (selective) Warning - message
        /// </summary>
        /// <param name="o"></param>
        public void Warning(object o) {
            if (base.DoLogging)
                m_Logger.Warn(o);

            Console.WriteLine(o);
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

        FuncTrace _f;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="Title">
        /// a title under which the block appears in the logfile
        /// </param>
        /// <param name="f">
        /// tracing of function which contains the block
        /// </param>
        public BlockTrace(string Title, FuncTrace f) {
            if(!Tracer.InstrumentationSwitch)
                return;
            string _name = Title;
            _f = f;
            m_Logger = _f.m_Logger;
            
            base.EnterMessage("BLKENTER ", _name);
        }

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
    }
}
