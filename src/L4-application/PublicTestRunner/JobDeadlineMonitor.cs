using BoSSS.Application.BoSSSpad;
using MathNet.Numerics.Optimization;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;

namespace PublicTestRunner {
    public class JobTimeEntry {
        public string name;
        public double avgSeconds;
        public double dueMargin;
        public DateTime lastupdate;

        public JobTimeEntry(string name, double avgSeconds) {
            this.name = name;
            this.avgSeconds = Math.Max(60, avgSeconds); 
            UpdateMargin();
            this.lastupdate = DateTime.Now;
        }

        internal void UpdateMargin() {
            var margin = this.avgSeconds*0.1;
            margin = Math.Max(60 * 5, margin); // minimum: 5 minutes, to avoid problems with very quick tests.
            margin = Math.Ceiling(margin / 60.0) * 60.0; // round up to full minutes
            this.dueMargin = margin;
        }
    }

    public class JobDeadlineMonitor {
        private string path;
        private bool shouldUpdateTimes;
        private Dictionary<string, JobTimeEntry> overview;

        public JobDeadlineMonitor(string DirectoryOffset, bool shouldUpdateTimes = false) {

            this.path = Path.Combine(DirectoryOffset, "TimeRecords.json");
            this.shouldUpdateTimes = shouldUpdateTimes;
            this.overview = new Dictionary<string, JobTimeEntry>();
            var filepath = path;
            if ( File.Exists(filepath) ) {
                var content = File.ReadAllText(filepath);
                this.overview = JsonConvert.DeserializeObject<Dictionary<string, JobTimeEntry>>(content);
            } else {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine($"No Benchmark file detected for {filepath}");
                Console.WriteLine("Continuing without it, be sure to manually update missing entries!");
                Console.ForegroundColor = ConsoleColor.White;
            }

            //this.shouldUpdateTimes = true;
            //foreach ( var entry in this.overview.Values )
            //    entry.UpdateMargin();
            //this.Save();
        }

        private string Trim(string name, string PrefixToTrim) {
            if ( name.StartsWith(PrefixToTrim) ) {
                name = name.Substring(PrefixToTrim.Length);
                name = name.TrimStart('-');
                return name;
            } else {
                return name;
            }
        }

        public bool JobExists(Job job, string PrefixToTrim) {
            return this.JobExists(Trim(job.Name, PrefixToTrim));
        }

        private bool JobExists(string name) {
            return this.overview.ContainsKey(name);
        }

        public bool Overdue(Job job, string PrefixToTrim) {
            TimeSpan span = job.LatestDeployment.RunTime;
            return this.Overdue(Trim(job.Name, PrefixToTrim), span.TotalSeconds);
        }

        private bool Overdue(string name, double currentSeconds) {
            if ( this.overview.TryGetValue(name, out var result) ) {
                return (currentSeconds > result.avgSeconds + result.dueMargin);
            }
            return false;
        }

        public void UpdateEntry(Job job, string PrefixToTrim) {
            TimeSpan span = job.LatestDeployment.RunTime;
            this.UpdateEntry(Trim(job.Name, PrefixToTrim), span.TotalSeconds);
        }

        private void UpdateEntry(string name, double seconds) {
            this.overview[name] = new JobTimeEntry(name, seconds);
        }

        public double GetAvgTime(Job job, string PrefixToTrim) {
            return GetAvgTime(Trim(job.Name, PrefixToTrim));
        }

        private double GetAvgTime(string name) {
            if ( this.overview.TryGetValue(name, out var time) ) {
                return time.avgSeconds;
            }
            return 0;
        }

        //public double GetPercentTime(Job job) {
        //    return this.GetPercentTime(job.Name, job.LatestDeployment.RunTime.TotalSeconds);
        //}

        public double GetPercentTime(string name, double currentSeconds) {
            return currentSeconds / this.GetAvgTime(name);
        }

        public void Save() {
            var content = JsonConvert.SerializeObject(this.overview, Formatting.Indented);
            var filepath = this.shouldUpdateTimes ? "TimeRecords.json" : "TimeRecordsTmp.json";
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"Saving Timings to: {Path.GetFullPath(filepath)}");
            Console.ForegroundColor = ConsoleColor.White;
            File.WriteAllText(filepath, content);
            try {
                File.Copy(filepath, @"C:\tmp\TimeRecordings\TimeRecords-" + DateTime.Now.ToFileTimeUtc() + ".json");
            } catch (Exception e) {
                Console.Error.WriteLine("Creating side-copy of TimeRecords: " + e);
            }
        }
    }
}
