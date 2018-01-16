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
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace ilPSP.Tracing {

    /// <summary>
    /// Accumulator for the total time spend in some traced method;
    /// </summary>
    [Serializable]
    public class MethodCallRecord {

        /// <summary>
        /// ctor.
        /// </summary>
        public MethodCallRecord(MethodCallRecord parrent, string n) {
            ParrentCall = parrent;
            this.Name = n;
        }

        /// <summary>
        /// method name
        /// </summary>
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
        public MethodCallRecord ParrentCall;

        /// <summary>
        /// all traces sub-calls
        /// </summary>
        public Dictionary<string, MethodCallRecord> Calls =
            new Dictionary<string, MethodCallRecord>();

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
        internal long m_TicksSpentInMethod = 0;

        /// <summary>
        /// Accumulated time in method.
        /// </summary>
        public long TicksSpentInMethod {
            get {
                return m_TicksSpentInMethod;
            }
        }


        /// <summary>
        /// Accumulated, total time spend in method (inclusive child calls).
        /// </summary>
        public TimeSpan TimeSpentInMethod {
            get {
                return new TimeSpan(m_TicksSpentInMethod);
            }
        }


        /// <summary>
        /// time spent in traced child calls
        /// </summary>
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
        public TimeSpan TimeSpentInChildCalls {
            get {
                return new TimeSpan(TicksSpentInChildCalls);
            }
        }

        /// <summary>
        /// Time spent the method itself, exclusive of child calls.
        /// </summary>
        public long TicksExclusive {
            get {
                return this.m_TicksSpentInMethod - TicksSpentInChildCalls;
            }
        }

        /// <summary>
        /// Time spent the method itself, exclusive of child calls.
        /// </summary>
        public TimeSpan TimeExclusive {
            get {
                return new TimeSpan(TicksExclusive);
            }
        }

        /// <summary>
        /// time spend in this method, relative to the <see cref="Root"/>;
        /// </summary>
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
        public int CallCount = 0;

        private static Regex WildcardToRegex(string pattern) {
            return new Regex("^" + Regex.Escape(pattern).
            Replace("\\*", ".*").
            Replace("\\?", ".") + "$");
        }

        /// <summary>
        /// root of the <see cref="MethodCallRecord"/>-tree, equal to <see cref="Tracer.Root"/>
        /// </summary>
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
                return string.Format(
                    "'{0}': {1} calls, {2:0.###E-00} seconds, {3:F3} %  of '{4}', called by '{5}'",
                    M.Name,
                    M.CallCount,
                    TimeSpentInMethod.TotalSeconds,
                    RuntimeFraction * 100,
                    RelativeRoot.Name,
                    M.ParrentCall != null ? M.ParrentCall.Name : "nobody");
            }
        }

        /// <summary>
        /// similar to <see cref="FindChildren(string)"/>, but returns exactly on child or null;
        /// </summary>
        public MethodCallRecord FindChild(string wildcard) {
            return FindChild(WildcardToRegex(wildcard));
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
            return FindChildren(WildcardToRegex(wildcard));
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
        /// time spend in this method, relative to the
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
            foreach (var mcr in this.AllCalls) {
                if (mcr.ParrentCall != null) {
                    if (!calledBy.Contains(mcr.ParrentCall.Name)) {
                        calledBy.Add(mcr.ParrentCall.Name);
                    }
                }
            }

            StringWriter stw = new StringWriter();


            stw.Write(string.Format(
                "'{0}': {1} calls, {2:F3}% / {4:0.##E-00} sec. runtime exclusive, {3:F3}% / {5:0.##E-00} sec. runtime inclusive child calls, called by: ",
                    this.Name,
                    this.CallCount,
                    this.ExclusiveTimeFractionOfRoot * 100,
                    this.TimeFractionOfRoot * 100,
                    (new TimeSpan(this.ExclusiveTicks)).TotalSeconds,
                    (new TimeSpan(this.TicksSpentInMethod)).TotalSeconds));

            stw.Write("'");
            for (int i = 0; i < calledBy.Count; i++) {
                stw.Write(calledBy[i]);
                if (i < calledBy.Count - 1)
                    stw.Write(",");
                else
                    stw.Write("'");
            }

            return stw.ToString();
        }
    }
}
