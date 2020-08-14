using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.Tracing;

namespace ilPSP
{
    public static class MethodCallRecordExtension
    {
        public static void PrintMostExpensiveCalls(this MethodCallRecord mcr, int count) {

            GetMostExpensiveCalls(Console.Out, mcr, count);
            Console.Out.Flush();

        }

        public static void GetMostExpensiveCalls(TextWriter wrt, MethodCallRecord R, int cnt = 0) {
            int i = 1;
            var mostExpensive = R.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveTicks);
            foreach (var cr in mostExpensive) {
                wrt.Write("Rank " + i + ": ");
                wrt.WriteLine(cr.ToString());
                if (i == cnt) return;
                i++;
            }
        }

        private struct Stats
        {
            public Stats(double[] times) {
                m_Max=times.Max();
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

        public static Dictionary<string, Tuple<double, double, int>> GetProfilingStats(MethodCallRecord[] mcr) {
            var kv = new Dictionary<string, Stats>();
            var methodImblance = new Dictionary<string, Tuple<double, double, int>>();
            List<string> method_names = new List<string>();

            mcr[0].CompleteCollectiveReport().ForEach(r => method_names.Add(r.Name));
            double[] rootTimes = new double[mcr.Length];

            for (int j = 0; j < mcr.Length; j++) {
                rootTimes[j] = mcr[j].TimeSpentInMethod.TotalSeconds;
            }

            var rootStat = new Stats(rootTimes);

            foreach (string method in method_names) {
                double[] times = new double[mcr.Length];

                for (int j = 0; j < times.Length; j++) {
                    mcr[j].FindChildren(method).ForEach(s => times[j]=s.TimeExclusive.TotalSeconds);
                }

                var TStats = new Stats(times);
                kv.Add(method, TStats);
                methodImblance.Add(method, new Tuple<double, double, int>(TStats.Imbalance / rootStat.Average, TStats.Imbalance, mcr[0].CallCount));
            }
            return methodImblance;
        }
    }
}
