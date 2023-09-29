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

using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text.RegularExpressions;

namespace ilPSP.Tracing {

    /// <summary>
    /// Accumulator for the total time spend in some traced method;
    /// </summary>
    [Serializable]
    [DataContract]
    public class MethodCallRecord {

        /// <summary>
        /// ctor.
        /// </summary>
        public MethodCallRecord(MethodCallRecord parrent, string n) {
            ParrentCall = parrent;
            this.Name = n;
        }

        /// <summary>
        /// Serializes this object and its childs into a JSON string.
        /// </summary>
        public string Serialize() {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };

            using (var wrt = new StringWriter()) {
                using (JsonWriter writer = new JsonTextWriter(wrt)) {  // Alternative: binary writer: JsonTextWriter
                    formatter.Serialize(writer, this);
                }

                return wrt.ToString();
            }

        }

        /// <summary>
        /// De-Serializes from a JSON string.
        /// </summary>
        public static MethodCallRecord Deserialize(string JsonString) {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };

            using (var rd = new StringReader(JsonString)) {
                MethodCallRecord mcr;
                using (JsonReader reader = new JsonTextReader(rd)) {  // Alternative: binary writer: JsonTextWriter
                    mcr = formatter.Deserialize<MethodCallRecord>(reader);
                }

                mcr.FixParrentRecursive();

                return mcr;
            }
        }

        void FixParrentRecursive() {
            foreach(var ch in this.Calls.Values) {
                ch.ParrentCall = this;
                ch.FixParrentRecursive();
            }
        }

        /// <summary>
        /// Sets the time spend in respective method (see <see cref="TicksSpentInMethod"/>) and in all child calls to zero. Can be configured to reset call count as well.
        /// </summary>
        public void ResetRecursive(bool EliminateKilledTimeFromParrent = true, bool EliminateCallcount = false) {

            if (EliminateKilledTimeFromParrent) {
                long toRemove = this.TicksExclusive;

                for (MethodCallRecord p = this.ParrentCall; p != null; p = p.ParrentCall) {

                    p.m_TicksSpentInMethod -= toRemove;
                }

            }
            this.m_TicksSpentInMethod -= this.TicksExclusive;

            foreach (var c in Calls.Values) {
                c.ResetRecursive(EliminateCallcount: EliminateCallcount);
            }
            if (EliminateCallcount) {
                this.CallCount = 0;
            }
        }


        /// <summary>
        /// method name
        /// </summary>
        [DataMember]
        public string Name {
            get;
            private set;
        }

        /*
        /// <summary>
        /// to prevent that 'dummy' blocks get accumulated in
        /// <see cref="TicksSpentInChildCalls"/>
        /// </summary>
        public bool IgnoreForExclusive = false;
        */

        /// <summary>
        /// the method which called this method
        /// </summary>
        [JsonIgnore]
        public MethodCallRecord ParrentCall;

        /// <summary>
        /// Depth in the call stack
        /// </summary>
        public int Depth {
            get {
                if(ParrentCall == null)
                    return 1;
                else
                    return ParrentCall.Depth + 1;
            }
        }

        /// <summary>
        /// all traces sub-calls
        /// </summary>
        [DataMember]
        public Dictionary<string, MethodCallRecord> Calls = new Dictionary<string, MethodCallRecord>();

        /// <summary>
        /// Adds a new entry to <see cref="Calls"/> (if <paramref name="_name"/> does not exist), otherwise 
        /// the existing entry is modified.
        /// </summary>
        public MethodCallRecord AddSubCall(string _name, long ElapsedTicks) {

            MethodCallRecord mcr;
            if (!this.Calls.TryGetValue(_name, out mcr)) {
                mcr = new MethodCallRecord(this, _name);
                this.Calls.Add(_name, mcr);
            }
            mcr.CallCount++;
            mcr.m_TicksSpentInMethod += ElapsedTicks;

            return mcr;
        }


        /// <summary>
        /// Accumulated time in method.
        /// </summary>
        [DataMember]
        internal long m_TicksSpentInMethod = 0;

        /// <summary>
        /// Accumulated time in method.
        /// </summary>
        [JsonIgnore]
        public long TicksSpentInMethod {
            get {
                return m_TicksSpentInMethod;
            }
        }

        [DataMember]
        internal long m_MemoryIncrease = 0;

        [DataMember]
        internal long m_PeakMemoryIncrease = 0;

        /// <summary>
        /// memory allocation of call without child calls
        /// </summary>
        [JsonIgnore]
        public long ExclusiveMemoryIncrease {
            get {
                long childs = 0;
                foreach (var c in Calls.Values) {
                    childs += c.m_MemoryIncrease;
                }
                return m_MemoryIncrease - childs;
            }
        }

        /// <summary>
        /// memory allocation of call including all child calls
        /// </summary>
        [JsonIgnore]
        public long InclusiveMemoryIncrease {
            get {
                return m_MemoryIncrease;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        internal long m_TicksSpentinBlocking = 0;

        [JsonIgnore]
        public TimeSpan TimeSpentInMPIBlocking {
            get {
                return new TimeSpan(m_TicksSpentinBlocking);
            }
        }

        /// <summary>
        /// Ticks spent in blocking MPI routines.
        /// ticks are exclusive the child calls.
        /// </summary>
        [JsonIgnore]
        public long ExclusiveBlockingTicks {
            get {
                long childs = 0;
                foreach (var c in Calls.Values) {
                    childs += c.m_TicksSpentinBlocking;
                }
                return m_TicksSpentinBlocking - childs;
            }
        }

        /// <summary>
        /// Time spent in blocking MPI routines exclusive the child calls
        /// </summary>
        public TimeSpan ExclusiveBlockingTime {
            get {
                return new TimeSpan(ExclusiveBlockingTicks);
            }
        }

        /// <summary>
        /// Accumulated, total time spend in method (inclusive child calls).
        /// </summary>
        [JsonIgnore]
        public TimeSpan TimeSpentInMethod {
            get {
                return new TimeSpan(m_TicksSpentInMethod);
            }
        }


        /// <summary>
        /// time spent in traced child calls
        /// </summary>
        [JsonIgnore]
        public long TicksSpentInChildCalls {
            get {
                long r = 0;
                foreach (var c in Calls.Values) {
                    //if (!c.IgnoreForExclusive) {
                    r += c.m_TicksSpentInMethod;
                    //}
                }
                return r;
            }
        }

        /// <summary>
        /// Accumulated, total time spend in method (inclusive child calls).
        /// </summary>
        [JsonIgnore]
        public TimeSpan TimeSpentInChildCalls {
            get {
                return new TimeSpan(TicksSpentInChildCalls);
            }
        }

        /// <summary>
        /// Time spent the method itself, exclusive of child calls.
        /// </summary>
        [JsonIgnore]
        public long TicksExclusive {
            get {
                return this.m_TicksSpentInMethod - TicksSpentInChildCalls;
            }
        }

        /// <summary>
        /// Time spent the method itself, exclusive of child calls.
        /// </summary>
        [JsonIgnore]
        public TimeSpan TimeExclusive {
            get {
                return new TimeSpan(TicksExclusive);
            }
        }

        /// <summary>
        /// time spend in this method, relative to the <see cref="Root"/>;
        /// </summary>
        [JsonIgnore]
        public double TimeFractionOfRoot {
            get {
                double rT = (double)(this.m_TicksSpentInMethod);
                double rR = (double)(this.Root.m_TicksSpentInMethod);
                return rT / rR;
            }
        }

        /// <summary>
        /// time spend in this method, relative to the <see cref="Root"/>;
        /// </summary>
        [JsonIgnore]
        public double ExclusiveTimeFractionOfRoot {
            get {
                double rT = (double)(this.TicksExclusive);
                double rR = (double)(this.Root.m_TicksSpentInMethod);
                return rT / rR;
            }
        }

        /// <summary>
        /// how often the method was called
        /// </summary>
        [DataMember]
        public int CallCount = 0;

     
        /// <summary>
        /// root of the <see cref="MethodCallRecord"/>-tree, equal to <see cref="Tracer.Root"/>
        /// </summary>
        [JsonIgnore]
        public MethodCallRecord Root {
            get {
                if (this.ParrentCall == null) {
                    //Debug.Assert(object.ReferenceEquals(this,Tracer._Root));
                    return this;
                } else {
                    return this.ParrentCall.Root;
                }

            }
        }

        /// <summary>
        /// creates an report for this method;
        /// </summary>
        public MiniReport GetMiniReport(MethodCallRecord altRoot = null) {
            if (altRoot == null)
                altRoot = this.Root;

            return new MiniReport() {
                M = this,
                RelativeRoot = altRoot
            };
        }

        /// <summary>
        /// summarizes some statistics about a method call
        /// <see cref="MiniReport.M"/>, created by <see cref="GetMiniReport"/>;
        /// </summary>
        public class MiniReport {

            /// <summary>
            /// Constructs an empty report
            /// </summary>
            internal MiniReport() {
            }

            /// <summary>
            /// see <see cref="RuntimeFraction"/>;
            /// </summary>
            public MethodCallRecord RelativeRoot {
                get;
                internal set;
            }

            /// <summary>
            /// the creator of this object;
            /// </summary>
            public MethodCallRecord M {
                get;
                internal set;
            }

            /// <summary>
            /// time spent in method traced by <see cref="M"/>
            /// </summary>
            TimeSpan TimeSpentInMethod {
                get {
                    return new TimeSpan(this.M.m_TicksSpentInMethod);
                }
            }

            /// <summary>
            /// the relative runtime with respect to <see cref="RelativeRoot"/>;
            /// </summary>
            public double RuntimeFraction {
                get {
                    double rT = (double)(this.M.m_TicksSpentInMethod);
                    double rR = (double)(this.RelativeRoot.m_TicksSpentInMethod);
                    return rT / rR;
                }
            }

            /// <summary>
            /// a report;
            /// </summary>
            override public string ToString() {
                using(var stw = new StringWriter()) {
                    stw.Write(string.Format(
                        "'{0}': {1} calls, {2:0.###E-00} seconds, {3:F3} %  of '{4}', {5} MB allocated, {6} MB allocated exclusive, called by: ",
                        M.Name, // 0
                        M.CallCount, // 1
                        TimeSpentInMethod.TotalSeconds, //  2
                        RuntimeFraction * 100, // 3 
                        RelativeRoot.Name, // 4
                        (double)(M.InclusiveMemoryIncrease) / (1024.0 * 1024.0), //  5
                        (double)(M.ExclusiveMemoryIncrease) / (1024.0 * 1024.0))); // 6

                    if(M.ParrentCall == null) {
                        stw.Write("NOBODY");
                    } else {
                        var mcr = M.ParrentCall;
                        while(mcr != null) {
                            stw.Write(mcr.Name);
                            mcr = mcr.ParrentCall;
                            if(mcr != null)
                                stw.Write(">");
                            else
                                stw.Write(">X");
                        }


                    }
                    return stw.ToString();
                }
            }
        }

        /// <summary>
        /// similar to <see cref="FindChildren(string)"/>, but returns exactly on child or null;
        /// </summary>
        public MethodCallRecord FindChild(string wildcard) {
            return FindChild(wildcard.WildcardToRegex());
        }

        /// <summary>
        /// similar to <see cref="FindChildren(Regex)"/>, but returns exactly on child or null;
        /// </summary>
        public MethodCallRecord FindChild(Regex r) {
            var childs = FindChildren(r);
            if (childs.Count() > 1) {
                throw new ArgumentException("unable to find unique selection");
            }
            if (childs.Count() <= 0)
                return null;
            return FindChildren(r).Single();
        }

        /// <summary>
        /// collects recursively all child calls (see <see cref="Calls"/>)
        /// whose name (<see cref="Name"/>) matches the wildcard string
        /// <paramref name="wildcard"/>;
        /// </summary>
        public IEnumerable<MethodCallRecord> FindChildren(string wildcard) {
            return FindChildren(wildcard.WildcardToRegex());
        }

        /// <summary>
        /// collects recursively all child calls (see <see cref="Calls"/>)
        /// whose name (<see cref="Name"/>) matches <paramref name="n"/>;
        /// </summary>
        public IEnumerable<MethodCallRecord> FindChildren(Regex n) {
            List<MethodCallRecord> ret = new List<MethodCallRecord>();
            this.FindChildrenRec(ret, n);
            return ret;
        }

        /// <summary>
        /// Recursive implementation of <see cref="FindChildren(Regex)"/>
        /// </summary>
        /// <param name="L"></param>
        /// <param name="n"></param>
        private void FindChildrenRec(List<MethodCallRecord> L, Regex n) {
            if (n.IsMatch(this.Name))
                L.Add(this);
            foreach (var c in Calls.Values) {
                c.FindChildrenRec(L, n);
            }
        }

        /// <summary>
        /// flattens the call tree and creates a collective description by
        /// grouping all calls with the same name
        /// </summary>
        public IEnumerable<CollectionReport> CompleteCollectiveReport() {
            var R = new Dictionary<string, List<MethodCallRecord>>();
            CollectiveReportRec(R);
            return R.Values.Select(l => new CollectionReport(l));
        }



        /// <summary>
        /// Recursive implementation of <see cref="CompleteCollectiveReport"/>
        /// </summary>
        /// <param name="R"></param>
        private void CollectiveReportRec(Dictionary<string, List<MethodCallRecord>> R) {
            List<MethodCallRecord> H;
            if (!R.TryGetValue(this.Name, out H)) {
                H = new List<MethodCallRecord>();
                R.Add(this.Name, H);
            }
            H.Add(this);

            foreach (MethodCallRecord c in this.Calls.Values)
                c.CollectiveReportRec(R);
        }


        /// <summary>
        /// flattens this entry and all childs recusively into a list
        /// </summary>
        public IEnumerable<MethodCallRecord> Flatten() {
            var Ret = new List<MethodCallRecord>();
            FlattenRecursive(Ret);
            return Ret;
        }


        private void FlattenRecursive(List<MethodCallRecord> r) {
            r.Add(this);
            foreach(var mcr in this.Calls.Values)
                mcr.FlattenRecursive(r);
        }

        /// <summary>
        /// 
        /// </summary>
        public override string ToString() {
            var m = this.GetMiniReport();
            return m.ToString();
        }

    }

    /// <summary>
    /// Provides summarized statistics on a range of method calls with potentially
    /// different parent's (see <see cref="MethodCallRecord.ParrentCall"/>).
    /// </summary>
    public class CollectionReport {

        /// <summary>
        /// ctor.
        /// </summary>
        public CollectionReport(ICollection<MethodCallRecord> ac) {
            AllCalls = ac.ToArray();
        }

        private MethodCallRecord[] AllCalls;

        /// <summary>
        /// name of the method
        /// </summary>
        public string Name {
            get {
                return AllCalls.First().Name;
            }
        }

        /// <summary>
        /// Total number of calls
        /// </summary>
        public int CallCount {
            get {
                return AllCalls.Sum(a => a.CallCount);
            }
        }

        /// <summary>
        /// total ticks spend in the whole method
        /// </summary>
        public long TicksSpentInMethod {
            get {
                return AllCalls.Sum(a => a.m_TicksSpentInMethod);
            }
        }

        /// <summary>
        /// ticks spent in the method itself (not inside traced child calls)
        /// </summary>
        public long ExclusiveTicks {
            get {
                return AllCalls.Sum(a => a.TicksExclusive);
            }
        }

        /// <summary>
        /// ticks spent in traced child calls
        /// </summary>
        public long TicksSpentInChildCalls {
            get {
                return AllCalls.Sum(a => a.TicksSpentInChildCalls);
            }
        }

        /// <summary>
        /// ticks spent in blocking MPI routines within traced child calls
        /// </summary>
        public long ExclusiveBlockingTicks {
            get {
                return AllCalls.Sum(a => a.ExclusiveBlockingTicks);
            }
        }

        /// <summary>
        /// Memory spent in traced child calls
        /// </summary>
        public long ExclusiveMemoryIncrease {
            get {
                return AllCalls.Sum(a => a.ExclusiveMemoryIncrease);
            }
        }

        /// <summary>
        /// Memory spent in traced child calls
        /// </summary>
        public long InclusiveMemoryIncrease {
            get {
                return AllCalls.Sum(a => a.InclusiveMemoryIncrease);
            }
        }

        /// <summary>
        /// time spend in this method, relative to the
        /// note: not usable in case of recursive calls
        /// <see cref="MethodCallRecord.Root"/>;
        /// </summary>
        public double TimeFractionOfRoot {
            get {
                var Root = this.AllCalls[0].Root;

                double rT = (double)(this.TicksSpentInMethod);
                double rR = (double)(Root.m_TicksSpentInMethod);
                return rT / rR;
            }
        }

        /// <summary>
        /// time spent in this method, relative to the
        /// <see cref="MethodCallRecord.Root"/>;
        /// </summary>
        public double ExclusiveTimeFractionOfRoot {
            get {
                var Root = this.AllCalls[0].Root;
                double rT = (double)(this.ExclusiveTicks);
                double rR = (double)(Root.m_TicksSpentInMethod);
                return rT / rR;
            }
        }

        /// <summary>
        /// formats a runtime report
        /// </summary>
        override public string ToString() {
            List<string> calledBy = new List<string>();
            List<int> calledByCnt = new List<int>();
            foreach (var mcr in this.AllCalls) {
                if (mcr.ParrentCall != null) {
                    if (!calledBy.Contains(mcr.ParrentCall.Name)) {
                        calledBy.Add(mcr.ParrentCall.Name);
                        calledByCnt.Add(mcr.CallCount);
                    }
                }
            }

            using(StringWriter stw = new StringWriter()) {


                stw.Write($"'{this.Name}': {this.CallCount} calls, ");
                stw.Write($"{(this.ExclusiveTimeFractionOfRoot * 100):F3}% / {(new TimeSpan(this.ExclusiveTicks)).TotalSeconds:0.##E-00} sec. runtime exclusive, ");
                stw.Write($"{(this.TimeFractionOfRoot * 100):F3}% / {(new TimeSpan(this.TicksSpentInMethod)).TotalSeconds:0.##E-00} sec. runtime inclusive child calls, ");
                stw.Write($"{(double)(this.InclusiveMemoryIncrease) / (1024.0 * 1024.0):0.##E-00}/{(double)(this.ExclusiveMemoryIncrease) / (1024.0 * 1024.0):0.##E-00} MB allocated incl./excl., ");
                stw.Write("called by: ");
                        
                stw.Write("'");
                for(int i = 0; i < calledBy.Count; i++) {
                    stw.Write(calledBy[i]);
                    stw.Write("*");
                    stw.Write(calledByCnt[i]);
                    if(i < calledBy.Count - 1)
                        stw.Write(",");
                    else
                        stw.Write("'");
                }

                return stw.ToString();
            }
        }
    }
}
