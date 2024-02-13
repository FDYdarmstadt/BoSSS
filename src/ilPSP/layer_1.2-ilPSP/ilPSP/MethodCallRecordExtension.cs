using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.Tracing;
using ilPSP.Utils;

namespace ilPSP {
    public static class MethodCallRecordExtension {
        
        
        public static void PrintMostExpensiveCalls(this MethodCallRecord mcr, int count) {

            GetMostExpensiveCalls(Console.Out, mcr, count);
            Console.Out.Flush();

        }

        public static void PrintMostExpensiveBlocking(this MethodCallRecord mcr, int count) {

            GetMostExpensiveBlocking(Console.Out, mcr, count);
            Console.Out.Flush();

        }

        public static void PrintMostExpensiveMemory(this MethodCallRecord mcr, int count) {
            GetMostExpensiveMemory(Console.Out, mcr, count);
            Console.Out.Flush();
        }
       

        /// <summary>
        /// Writes a summary on 
        /// most expensive calls
        /// (sum over all calls, i.e. no distinction by parent) 
        /// </summary>
        public static void GetMostExpensiveCalls(TextWriter wrt, MethodCallRecord R, int cnt = 0) {
            int i = 1;
            var mostExpensive = R.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveTicks);
            foreach(var cr in mostExpensive) {
                wrt.Write("Rank " + i + ": \t");
                wrt.Write($"{(cr.ExclusiveTimeFractionOfRoot * 100):F3}%\t{(new TimeSpan(cr.ExclusiveTicks)).TotalSeconds:0.##E-00}\t");
                wrt.WriteLine(cr.ToString());
                if(i == cnt) return;
                i++;
            }
        }

        /// <summary>
        /// Writes a summary on 
        /// most expensive calls
        /// (sum over all calls, i.e. no distinction by parent) 
        /// </summary>
        public static void GetMostMemoryConsumingCalls(TextWriter wrt, MethodCallRecord R, int cnt = 0) {
            int i = 1;
            var mostExpensive = R.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveMemoryIncrease);
            foreach(var cr in mostExpensive) {
                wrt.Write("Rank " + i + ": \t");
                wrt.Write($"{Math.Round((double)cr.ExclusiveMemoryIncrease / (1024.0 * 1024.0))}\t");
                wrt.WriteLine(cr.ToString());
                if(i == cnt) return;
                i++;
            }
        }
        
        /// <summary>
        /// Writes a detailed summary on 
        /// most expensive calls
        /// (distinction by parent calls) 
        /// </summary>
        public static void GetMostExpensiveCallsDetails(TextWriter wrt, MethodCallRecord R, int cnt = 0) {
            int i = 1;
            var mostExpensive = R.Flatten().OrderByDescending(cr => cr.ExclusiveTimeFractionOfRoot);
            foreach(MethodCallRecord cr in mostExpensive) {
                wrt.Write("Rank " + i + ": \t");
                wrt.Write($"{(cr.ExclusiveTimeFractionOfRoot * 100):F3}%\t{cr.TimeExclusive.TotalSeconds:0.##E-00}\t");
                wrt.WriteLine(cr.ToString());
                if(i == cnt) return;
                i++;
            }
        }

        /// <summary>
        /// Writes a detailed summary on 
        /// most expensive calls
        /// (distinction by parent calls) 
        /// </summary>
        public static void GetMostMemoryConsumingCallsDetails(TextWriter wrt, MethodCallRecord R, int cnt = 0) {
            int i = 1;
            var mostExpensive = R.Flatten().OrderByDescending(cr => cr.ExclusiveMemoryIncrease);
            foreach(var cr in mostExpensive) {
                wrt.Write("Rank " + i + ": \t");
                wrt.Write($"{Math.Round((double)cr.ExclusiveMemoryIncrease / (1024.0 * 1024.0))}\t");
                
                wrt.WriteLine(cr.ToString());
                if(i == cnt) return;
                i++;
            }
        }


        private struct Stats {
            public Stats(double[] times) {
                m_Max = times.Max();
                m_Min = times.Min();
                m_Average = times.Sum() / times.Length;
                m_Imbal = m_Max - m_Min;
                Debug.Assert(m_Imbal >= 0);
                Debug.Assert(m_Average >= 0);
            }
            private double m_Min;
            private double m_Max;
            private double m_Average;
            private double m_Imbal;

            public double Min {
                get { return m_Min; }
            }
            public double Max {
                get { return m_Max; }
            }
            public double Average {
                get { return m_Average; }
            }
            public double Imbalance {
                get { return m_Imbal; }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public static Dictionary<string, (double RelInbalance, double Imbalance, int CallCount)> GetFuncImbalance(MethodCallRecord[] mcrs) {
            return GetImbalance(mcrs, s => s.TimeExclusive.TotalSeconds);
        }

        /// <summary>
        /// 
        /// </summary>
        public static Dictionary<string, (double RelInbalance, double Imbalance, int CallCount)> GetMPIImbalance(MethodCallRecord[] mcrs) {
            return GetImbalance(mcrs, s => s.ExclusiveBlockingTime.TotalSeconds);
        }

        /// <summary>
        /// Analyzes imbalances across MPI cores
        /// </summary>
        /// <param name="mcrs">Instrumentation of a run;</param>
        /// <param name="ProperyEval">
        /// evaluates a numeric property (r.g. some measured runtime, <see cref="MethodCallRecord.ExclusiveBlockingTime"/>) which is then compared across several MPI cores
        /// </param>
        /// <returns>
        /// - dictionary key: method name
        /// - value 1 (relative imbalance): imbalance / average
        /// - value 2 (imbalance): difference between maximum and minimum value across all MPI cores
        /// - value 3: how often the function was called
        /// </returns>
        public static Dictionary<string, (double RelInbalance, double Imbalance, int CallCount)> GetImbalance(MethodCallRecord[] mcrs, Func<MethodCallRecord, double> ProperyEval) {
            var kv = new Dictionary<string, Stats>();
            var methodImblance = new Dictionary<string, (double, double, int)>();
            List<string> method_names = new List<string>();

            mcrs[0].CompleteCollectiveReport().ForEach(r => method_names.Add(r.Name));
            double[] rootTimes = new double[mcrs.Length];

            for(int j = 0; j < mcrs.Length; j++) {
                rootTimes[j] = mcrs[j].TimeSpentInMethod.TotalSeconds;
            }

            var rootStat = new Stats(rootTimes);

            foreach(string method in method_names) {
                double[] times = new double[mcrs.Length];
                int cnt = 0;

                for(int j = 0; j < times.Length; j++) {
                    mcrs[j].FindChildren(method).ForEach(s => times[j] += ProperyEval(s));
                }

                mcrs[0].FindChildren(method).ForEach(s => cnt += s.CallCount);

                var TStats = new Stats(times);
                kv.Add(method, TStats);
                methodImblance.Add(method, (TStats.Imbalance / rootStat.Average, TStats.Imbalance, cnt));
            }
            return methodImblance;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="wrt"></param>
        /// <param name="R"></param>
        /// <param name="printcnt"></param>
        private static void GetMostExpensiveBlocking(TextWriter wrt, MethodCallRecord R, int printcnt = 0) {
            int i = 1;
            var mostExpensive = R.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveBlockingTicks);
            foreach(var kv in mostExpensive) {
                wrt.Write("#" + i + ": ");
                wrt.WriteLine(string.Format(
                "'{0}': {1} calls, {2:0.##E-00} sec. exclusive runtime",
                    kv.Name,
                    kv.CallCount,
                    new TimeSpan(kv.ExclusiveBlockingTicks).TotalSeconds));
                if(i == printcnt) return;
                i++;
            }
            Console.Out.Flush();
        }

        private static void GetMostExpensiveMemory(TextWriter wrt, MethodCallRecord R, int printcnt = 0) {
            int i = 1;
            var mostExpensive = R.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveMemoryIncrease);
            foreach(var kv in mostExpensive) {
                wrt.Write("#" + i + ": ");
                wrt.WriteLine(string.Format(
                "'{0}': {1} calls, {2} MB, exclusive memory",
                    kv.Name,
                    kv.CallCount,
                    kv.ExclusiveMemoryIncrease));
                if(i == printcnt) return;
                i++;
            }
            Console.Out.Flush();
        }
    }
}
