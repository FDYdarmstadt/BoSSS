using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using ilPSP.Utils;
//using static BoSSS.Solution.Gnuplot.Plot2Ddata;
using NUnit.Framework;


namespace BoSSS.Application.BoSSSpad {


    /// <summary>
    /// Combines all `memory.*.txt`-files in a session directory into a single time-line
    /// </summary>
    public class SessionMemtrace : MemtraceAnalyze {

        /// <summary>
        /// Combines all `memory.*.txt`-files in a session directory into a single time-line
        /// </summary>
        /// <param name="MaxLines">
        /// Limit to the number of lines to be analyzed; the algorithm might get to expensive, otherwise for long runs
        /// </param>
        /// <param name="LineOffset">
        /// Lines to skip at the beginning
        /// </param>
        /// <param name="SessionDirectory"></param>
        public SessionMemtrace(DirectoryInfo SessionDirectory, int LineOffset = 0, int MaxLines = 10000) {
            string FilePattern = $"memory.*.txt";

            FileInfo[] AllFiles = SessionDirectory.GetFiles(FilePattern);

            MPIsize = AllFiles.Length;

            myRecord[][] AllContent = new myRecord[AllFiles.Length][];
            for (int i = 0; i < AllFiles.Length; i++) {
                var all = myRecord.FromFile(AllFiles[i].FullName);
                if (all.Length - LineOffset > MaxLines)
                    Console.WriteLine("Note: truncating memory trace ...");
                AllContent[i] = all.Skip(LineOffset).Take(MaxLines).ToArray();
            }

            CombinedRanksS = LongestCommonSubsequence(AllContent, false);
            base.Timeline = CombinedRanksS;
        }


        internal myRecord[] CombinedRanksS;
    
    
        /// <summary>
        /// Number of cores for the run
        /// </summary>
        public int MPIsize {
            get;
        }

        /// <summary>
        /// maximal memory over all MPI ranks over time
        /// </summary>
        public long[] MaximumMem {
            get {
                var r = new long[CombinedRanksS.Length];
                for (int i = 0; i < r.Length; i++) {
                    CombinedRanksS[i].GetMPIMem(out _, out long max, out _);
                    r[i] = max;
                }
                return r;
            }
        }

        /// <summary>
        /// maximal memory over all MPI ranks over time in megabytes
        /// </summary>
        public double[] MaximumMemMeg {
            get {
                return MaximumMem.Select(noOfBytes => noOfBytes / (1024.0 * 1024.0)).ToArray();
            }
        }

        /// <summary>
        /// average memory over all MPI ranks over time
        /// </summary>
        public long[] AverageMem {
            get {
                var r = new long[CombinedRanksS.Length];
                for (int i = 0; i < r.Length; i++) {
                    CombinedRanksS[i].GetMPIMem(out _, out _, out long avg);
                    r[i] = avg;
                }
                return r;
            }
        }

        /// <summary>
        /// average memory over all MPI ranks over time in megabytes
        /// </summary>
        public double[] AverageMemMeg {
            get {
                return AverageMem.Select(noOfBytes => noOfBytes / (1024.0 * 1024.0)).ToArray();
            }
        }

        /// <summary>
        /// minimal memory over all MPI ranks over time
        /// </summary>
        public long[] MinimumMem {
            get {
                var r = new long[CombinedRanksS.Length];
                for (int i = 0; i < r.Length; i++) {
                    CombinedRanksS[i].GetMPIMem(out var min, out _, out _);
                    r[i] = min;
                }
                return r;
            }
        }

        /// <summary>
        /// minimal memory over all MPI ranks over time in megabyte
        /// </summary>
        public double[] MinimumMemMegs {
            get {
                return MinimumMem.Select(noOfBytes => noOfBytes / (1024.0 * 1024.0)).ToArray();
            }
        }

        /// <summary>
        /// total memory (aka. sum) over all MPI ranks over time
        /// </summary>
        public long[] TotalMem {
            get {
                var r = new long[CombinedRanksS.Length];
                for (int i = 0; i < r.Length; i++) {
                    var sum = CombinedRanksS[i].GetMPIMem(out _, out _, out _);
                    r[i] = sum;
                }
                return r;
            }
        }


        /// <summary>
        /// total memory (aka. sum) over all MPI ranks over time in megabyte
        /// </summary>
        public double[] TotalMemMegs {
            get {
                return TotalMem.Select(noOfBytes => noOfBytes / (1024.0 * 1024.0)).ToArray();
            }
        }

        /// <summary>
        /// Reports the largest memory-allocating routines in descending order
        /// </summary>
        public (int TimelineIndex, double Megs, string Name)[] ReportLargestAllocators() {
            return ReportLargest(+1);
        }

        private (int TimelineIndex, double Megs, string Name)[] ReportLargest(double sign) {
            var ret = new List<(int TimelineIndex, double Megs, string Name)>();

            var _TotalMemMegs = TotalMemMegs;
            var timelNames = base.GetTimeLine();
            double Scale = ((int.MaxValue / 16) / Math.Max(_TotalMemMegs.Max(), BLAS.MachineEps));
            Scale *= sign;

            // calc allocation difference
            double PrevMegs = 0.0;
            for (int iLine = 0; iLine < base.NoOfTimeEntries; iLine++) {
                ret.Add((iLine, _TotalMemMegs[iLine] - PrevMegs, timelNames[iLine]));
                PrevMegs = _TotalMemMegs[iLine];
            }


            // sort
            int ComparerFunc(ValueTuple<int, double, string> A, ValueTuple<int, double, string> B) {

                double Megs_A = A.Item2;
                double Megs_B = B.Item2;

                return (int)Math.Round(Scale * (Megs_B - Megs_A));
            }
            ret.Sort(FuncComparerExtensions.ToComparer<ValueTuple<int, double, string>>(ComparerFunc));

            // return
            return ret.ToArray();
        }


        /// <summary>
        /// Reports the largest memory-freeing routines in descending order
        /// </summary>
        public (int TimelineIndex, double Megs, string Name)[] ReportLargestDeallocators() {

            return ReportLargest(-1);
        }

        /// <summary>
        /// Plotting of 
        /// minimum, average and maximum memory allocations over all MPI ranks over time
        /// </summary>
        public Plot2Ddata GetMinAvgMaxMemPlot() {
            var ret = new Plot2Ddata();
            int L = NoOfTimeEntries;


            ret.AddDataGroup(new Plot2Ddata.XYvalues(
                $"Min Mem [MegB] at {MPIsize} cores",
                L.ForLoop(i => (double)i),
                MinimumMemMegs),
                new PlotFormat(Style: Styles.Lines, lineColor: LineColors.Blue));

            ret.AddDataGroup(new Plot2Ddata.XYvalues(
                $"Max Mem [MegB] at {MPIsize} cores",
                L.ForLoop(i => (double)i),
                MaximumMemMeg),
                new PlotFormat(Style: Styles.Lines, lineColor: LineColors.Red));

            ret.AddDataGroup(new Plot2Ddata.XYvalues(
                $"Avg Mem [MegB] at {MPIsize} cores",
                L.ForLoop(i => (double)i),
                AverageMemMeg),
                new PlotFormat(Style: Styles.Lines, lineColor: LineColors.Black));

            ret.Title = "Memory of session ";


            return ret;
        }

        /// <summary>
        /// total memory (aka. sum) over all MPI ranks over time
        /// </summary>
        public Plot2Ddata GetMPItotalMemory() {
            
            int L = NoOfTimeEntries;

            var ret = new Plot2Ddata();

            ret.AddDataGroup(new Plot2Ddata.XYvalues(
                $"Tot Mem [MegB] at {MPIsize} cores",
                L.ForLoop(i => (double)i),
                TotalMemMegs));

            ret.Title = "Total memory of session ";


            return ret;
        }

        /// <summary>
        /// writes a CSV
        /// </summary>
        public void WriteCombinedRanksOutfile(string path) {
            using (var stw = new StreamWriter(path)) {

                long PrevLine = 0;
                foreach (var e in this.CombinedRanksS) {
                    long TotMem = e.GetMPIMem(out long min, out long max, out long avg);

                    stw.Write(TotMem);
                    stw.Write("\t");
                    stw.Write(min);
                    stw.Write("\t");
                    stw.Write(max);
                    stw.Write("\t");
                    stw.Write(avg);
                    stw.Write("\t");

                    long diff = TotMem - PrevLine;
                    PrevLine = TotMem;
                    stw.Write(diff);
                    stw.Write("\t");

                    stw.Write(Math.Round(diff / 1024.0 / 1024.0));
                    stw.Write("\t");


                    stw.Write(e.Name);
                    stw.WriteLine();
                }


                stw.Flush();
                stw.Close();
            }
        }

    }

    /// <summary>
    /// Combines all `memory.*.txt`-files in multiple session directory into a single time-line
    /// </summary>
    public class SessionsComparisonMemtrace : MemtraceAnalyze {

        /// <summary>
        /// Combines all `memory.*.txt`-files in multiple session directory into a single time-line
        /// </summary>
        public SessionsComparisonMemtrace(DirectoryInfo[] SessionDirectories, int LineOffset = 0, int MaxLines = 10000) {

            SessionMemtrace[] singleRuns = new SessionMemtrace[SessionDirectories.Length];
            for (int i = 0; i < singleRuns.Length; i++)
                singleRuns[i] = new SessionMemtrace(SessionDirectories[i], LineOffset, MaxLines);


            CombinesRuns = LongestCommonSubsequence(singleRuns.Select(run => run.CombinedRanksS).ToArray(), true);
            base.Timeline = CombinesRuns;
        }

        myRecord[] CombinesRuns;


        /// <summary>
        /// MPI size (number of cores) for each run in this log
        /// </summary>
        public int[] MPIsizeOfRuns {
            get {
                List<int> sz = new List<int>();

                var f0 = CombinesRuns.First();

                sz.Add(f0.MPISize);

                foreach(var p in f0.peerRun) {
                    sz.Add(p.MPISize);
                }

                sz.Sort();
                return sz.ToArray();
            }
        }

        /// <summary>
        /// The total memory, in bytes, for each MPI size in <see cref="MPIsizeOfRuns"/>
        /// </summary>
        public (int MPISz, long[] TotalMem)[] GetTotalMemory() {

            int[] sizes = this.MPIsizeOfRuns;
            int NoOfLines = CombinesRuns.Length;

            (int MPISz, long[] TotalMem)[] ret = sizes.Select(sz => (sz, new long[NoOfLines])).ToArray();




            List<myRecord> line = new List<myRecord> (); 
            for(int iLine = 0; iLine < NoOfLines; iLine++) {
                line.Clear();
                line.Add(CombinesRuns[iLine]);
                line.AddRange(CombinesRuns[iLine].peerRun);

                bool[] marker = new bool[line.Count];

                // search for `line`-entry with matching rank; 
                // be very careful to do this in a stable fashion, i.e. in the case of two equal MPI sizes, don't mess up the sequence

                for(int k = 0; k < marker.Length; k++) {
                    int sz_k = sizes[k];
                    for(int i = 0; i < marker.Length; i++) {
                        if (marker[i])
                            continue;
                        if (line[i].MPISize == sz_k) {
                            Debug.Assert(line[i].MPISize == sz_k);
                            ret[i].TotalMem[iLine] = line[i].GetMPIMem(out long _, out long _, out long _);
                        }
                    }
                }


            }


            return ret;
        }


        /// <summary>
        /// The total memory, in megabytes, for each MPI size in <see cref="MPIsizeOfRuns"/>
        /// </summary>
        public (int MPISz, double[] TotalMem)[] GetTotalMemoryMegs() {
            var a = GetTotalMemory();

            return a.Select(tt => (
                tt.MPISz, 
                tt.TotalMem.Select(noOfBytes => noOfBytes / (1024.0 * 1024.0)).ToArray())).ToArray();
        }


        /// <summary>
        /// Reports the largest differences in memory allocation between the multiple runs
        /// </summary>
        public (int TimelineIndex, double Imbalance, double[] AllocMegs, string Name)[] ReportLargestAllocatorImbalance() {

            var ret = new List<(int TimelineIndex, double Imbalance, double[] AllocMegs, string Name)>();

            var _TotalMemMegs = GetTotalMemoryMegs();
            var timelNames = base.GetTimeLine();
            var MPIsizes = MPIsizeOfRuns;
            int NoOfRuns = MPIsizes.Length; ;

            // calc allocation difference
            double[] PrevMegs = new double[NoOfRuns];
            double maxMegs = BLAS.MachineEps;
            for (int iLine = 0; iLine < base.NoOfTimeEntries; iLine++) {
                double[] allocMegs_iLine = new double[NoOfRuns];
                for(int iRun = 0; iRun < NoOfRuns; iRun++) {
                    double megs = _TotalMemMegs[iRun].TotalMem[iLine];
                    allocMegs_iLine[iRun] = megs - PrevMegs[iRun];
                    PrevMegs[iRun] = megs;
                    maxMegs = Math.Max(maxMegs, megs);
                }

                double ImBalance = allocMegs_iLine.Max() - allocMegs_iLine.Min();


                ret.Add((iLine, ImBalance, allocMegs_iLine, timelNames[iLine]));
            }

            double Scale = (int.MaxValue / 16) / maxMegs;
            // sort
            int ComparerFunc(ValueTuple<int, double, double[], string> A, ValueTuple<int, double, double[], string> B) {

                double Megs_A = A.Item2;
                double Megs_B = B.Item2;

                return (int)Math.Round(Scale * (Megs_B - Megs_A));
            }
            ret.Sort(FuncComparerExtensions.ToComparer<ValueTuple<int, double, double[], string>>(ComparerFunc));

            // return
            return ret.ToArray();
        }

        /*
        static void WriteCombinedRunsOutfile(string path, IEnumerable<myRecord> CombList) {
            using (var stw = new StreamWriter(path)) {

                
                //int[] MPIsizes =   //allruns.Select(list => list.ElementAt(0).MPISize).ToArray();
                int NoOfRuns = 1 + CombList.ElementAt(0).peerRun.Count;


                foreach (var e in CombList) {

                    for (int i = 0; i < NoOfRuns; i++) {
                        var e_run = i == 0 ? e : e.peerRun[i - 1];
                        long TotMem = e_run.GetMPIMem(out long _, out long _, out long _);

                        if (e.Name != e_run.Name)
                            throw new ApplicationException();

                        stw.Write(TotMem);
                        stw.Write("\t");
                    }


                    stw.Write(e.Name);
                    stw.WriteLine();
                }


                stw.Flush();
                stw.Close();
            }
        }
        */

    }



    /// <summary>
    /// IO for memory tracing files in session directory;
    /// </summary>
    public abstract class MemtraceAnalyze {
        
        /*
        /// <summary>
        /// Combines all `memory.*.txt`-files in a session directory into a single time-line
        /// </summary>
        /// <param name="SessionDirectory"></param>
        public MemtraceAnalyze(DirectoryInfo SessionDirectory) {

            string FilePattern = $"memory.*.csv";
           
            FileInfo[] AllFiles = SessionDirectory.GetFiles(FilePattern);

            myRecord[][] AllContent = new myRecord[AllFiles.Length][];
            for (int i = 0; i < AllFiles.Length; i++) {
                AllContent[i] = myRecord.FromFile(AllFiles[i].FullName);
            }

            CombinedRanksS = LongestCommonSubsequence(AllContent, false);

        }


        myRecord[] CombinedRanksS;
        */

        /// <summary>
        /// internal ctor
        /// </summary>
        protected MemtraceAnalyze() {

        }

        /// <summary>
        /// memory tracing entries over time
        /// </summary>
        internal protected myRecord[] Timeline;






        /// <summary>
        /// Time-line, i.e. the name/call-stack of memory allocation
        /// </summary>
        public string[] GetTimeLine() {
            return Timeline.Select(e => e.Name).ToArray();
        }

        /// <summary>
        /// No of Lines, i.e. records over time
        /// </summary>
        public int NoOfTimeEntries {
            get {
                return Timeline.Length;
            }
        }


        /// <summary>
        /// memory tracing record
        /// </summary>
        internal protected class myRecord : IComparable<myRecord>, IEquatable<myRecord>, ICloneable {

            /// <summary>
            /// parsing from single line in text file
            /// </summary>
            public static myRecord ParseString(string s) {
                //line No | Mem in Bytes | mem in megs | diff-mem in bytes | diff-mem in megs | Call Stack
                //13130   440481432       420     31280   0       <BoSSS.Foundation.Quadrature.Quadrature`2.Execute>BoSSS.Solution.AdvancedSolvers.AggregationGridBasis.BuildInjector_Lv1>Aggregation_basis_init>BoSSS.Solution.Application`1.RunSolverMode

                string[] Parts = s.Split('\t');
                if (Parts.Length != 6)
                    throw new IOException();

                return new myRecord() {
                    line = int.Parse(Parts[0]),
                    Mem = long.Parse(Parts[1]),
                    //MemDiff = long.Parse(Parts[3]),
                    Name = Parts[5]
                };

            }

            /// <summary>
            /// 
            /// </summary>
            static public myRecord[] FromFile(string fname) {
                using (var f = new StreamReader(fname)) {
                    var ret = new List<myRecord>();
                    for (string s = f.ReadLine(); s != null; s = f.ReadLine()) {
                        try {
                            ret.Add(myRecord.ParseString(s));
                        } catch (IOException) {

                        }
                    }

                    return ret.ToArray();
                }
            }


            /// <summary>
            /// 
            /// </summary>
            public int CompareTo(myRecord other) {
                return this.Name.CompareTo(other.Name);
            }

            /// <summary>
            /// The call stack at the point of time at which the memory usage was measured
            /// </summary>
            public string Name;

            /// <summary>
            /// 
            /// </summary>
            public int line;

            /// <summary>
            /// memory in bytes, measured at the respective point of time
            /// </summary>
            public long Mem;

            /*
            /// <summary>
            /// allocation with respect to the previous line
            /// </summary>
            public long MemDiff;
            */

            List<myRecord> mpiBrothers = new List<myRecord>();

            /// <summary>
            /// Adds another MPI rank
            /// </summary>
            public void AddMpiBrother(myRecord b) {
                if (!b.Name.Equals(this.Name))
                    throw new ArgumentException();
                mpiBrothers.Add(b);
            }

            /// <summary>
            /// true if this is a combination of multiple MPI entries.
            /// </summary>
            public bool isMPICombined {
                get {
                    return mpiBrothers.Count > 0;
                }
            }

            /// <summary>
            /// Number of MPI cores 
            /// </summary>
            public int MPISize {
                get {
                    return mpiBrothers.Count + 1;
                }
            }

            /// <summary>
            /// Entries (of other runs or ranks or runs) which have been paired with this one
            /// </summary>
            public List<myRecord> peerRun = new List<myRecord>();

            /// <summary>
            /// Adds another MPI rank
            /// </summary>
            public void AddPeerRun(myRecord b) {
                if (!b.Name.Equals(this.Name))
                    throw new ArgumentException();
                peerRun.Add(b);
            }


            /// <summary>
            /// 
            /// </summary>
            public long GetMPIMem(out long min, out long max, out long avg) {
                long sum = this.Mem;
                min = sum;
                max = sum;

                foreach (var b in mpiBrothers) {
                    sum += b.Mem;
                    min = Math.Min(min, b.Mem);
                    max = Math.Max(max, b.Mem);
                }

                avg = sum / (mpiBrothers.Count + 1);

                return sum;
            }


            /// <summary>
            /// 
            /// </summary>
            public override bool Equals(object obj) {
                var otherRec = obj as myRecord;
                if (otherRec == null)
                    return false;
                return this.Name.Equals(otherRec.Name);
            }

            /// <summary>
            /// 
            /// </summary>
            public bool Equals(myRecord other) {
                return this.Name.Equals(other.Name);
            }

            /// <summary>
            /// 
            /// </summary>
            public override int GetHashCode() {
                return this.line;
            }

            /// <summary>
            /// non-shallow copy
            /// </summary>
            public object Clone() {
                var r = new myRecord() {
                    line = this.line,
                    Mem = this.Mem,
                    Name = this.Name,
                };

                foreach (var p in this.peerRun) {
                    r.peerRun.Add(p.CloneAs());
                }

                foreach (var bro in this.mpiBrothers) {
                    r.mpiBrothers.Add(bro.CloneAs());
                }

                return r;
            }
        }


        internal static int[,] LongestCommonSubsequence(myRecord[] s1, myRecord[] s2) {
            int[,] c = new int[s1.Length + 1, s2.Length + 1];
            for (int i = 1; i <= s1.Length; i++)
                for (int j = 1; j <= s2.Length; j++) {
                    if (s1[i - 1].Equals(s2[j - 1]))
                        c[i, j] = c[i - 1, j - 1] + 1;
                    else
                        c[i, j] = c[i - 1, j] > c[i, j - 1] ? c[i - 1, j] : c[i, j - 1];
                }
            return c;
        }


        internal static myRecord[] Cat(myRecord[] a, myRecord[] b) {
            var ret = new myRecord[a.Length + b.Length];
            Array.Copy(a, 0, ret, 0, a.Length);
            Array.Copy(b, 0, ret, a.Length, b.Length);
            return ret;
        }

        /// <summary>
        /// Reference verions from Wikipedia;
        /// Causes stack overflow for larger data sets
        /// </summary>
        internal static myRecord[] BackTrack(int[,] c, myRecord[] s1, myRecord[] s2, int i, int j, bool RunOrRankCombine) {
            if (i == 0 || j == 0) {
                return new myRecord[0];
            } else if (s1[i - 1].Equals(s2[j - 1])) {
                if (!RunOrRankCombine)
                    s1[i - 1].AddMpiBrother(s2[j - 1]); // 
                else
                    s1[i - 1].AddPeerRun(s2[j - 1]); // 

                return Cat(BackTrack(c, s1, s2, i - 1, j - 1, RunOrRankCombine), new[] { s1[i - 1] });
            } else if (c[i, j - 1] > c[i - 1, j]) {
                return BackTrack(c, s1, s2, i, j - 1, RunOrRankCombine);
            } else {
                return BackTrack(c, s1, s2, i - 1, j, RunOrRankCombine);
            }
  
        }

        /// <summary>
        /// Non-recursive version
        /// </summary>
        internal static myRecord[] BackTrackNonRecursive(int[,] matrix, myRecord[] s1, myRecord[] s2, bool RunOrRankCombine) {
            int i = s1.Length;// self.characters.count
            int j = s2.Length;// other.characters.count
           

            var lcs = new List<myRecord>();

            while (i > 0 && j > 0) {
                if (matrix[i,j] == matrix[i,j - 1]) {
                    // Indicates propagation without change: no new char was added to lcs.
                    j -= 1;
                } else if (matrix[i,j] == matrix[i - 1,j]) {
                    // Indicates propagation without change: no new char was added to lcs.
                    i -= 1;
                } else {
                    // Value on the left and above are different than current cell.
                    // This means 1 was added to lcs length.

                    Debug.Assert(s1[i - 1].Equals(s2[j - 1]));
                    if (!RunOrRankCombine)
                        s1[i - 1].AddMpiBrother(s2[j - 1]); // 
                    else
                        s1[i - 1].AddPeerRun(s2[j - 1]); // 

                    i -= 1;
                    j -= 1;

                    lcs.Add(s1[i]);
                }
            }

            lcs.Reverse();

            return lcs.ToArray(); 
        }

        internal static myRecord[] LongestCommonSubsequence(IEnumerable<myRecord> a, IEnumerable<myRecord> b, bool RunOrRankCombine) {
            myRecord[] _a = a.ToArray();
            myRecord[] _b = b.ToArray();


            myRecord[] _aC = null;
            myRecord[] _bC = null;

            if(_a.Length*_b.Length < 400*400) {
                _aC = _a.Select(e => e.CloneAs()).ToArray();
                _bC = _b.Select(e => e.CloneAs()).ToArray();
            }

            var aa = LongestCommonSubsequence(_a.ToArray(), _b.ToArray());


            var r = BackTrackNonRecursive(aa, _a, _b, RunOrRankCombine);

            if (_aC != null && _bC != null) {
                var r1 = BackTrack(aa, _aC, _bC, _aC.Length, _bC.Length, RunOrRankCombine);
                if (!r1.ListEquals(r))
                    throw new ApplicationException("mismatch between backtracking implementations");
            }

            return r;
        }

        internal static myRecord[] LongestCommonSubsequence(IEnumerable<myRecord>[] lists, bool RunOrRankCombine) {
            myRecord[] ret = lists[0].ToArray();
            for (int i = 1; i < lists.Length; i++) {
                //Console.Write($"combining {i} of {lists.Length} ...");
                ret = LongestCommonSubsequence(ret, lists[i], RunOrRankCombine);
                //Console.WriteLine("done.");
            }
            return ret;
        }




    }
}
