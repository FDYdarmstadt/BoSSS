using BoSSS.Application.BoSSSpad;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;

namespace PublicTestRunner {
    internal class JobTimeEntry {
        public string name;
        public double avgSeconds;
        public double dueMargin;
        public DateTime lastupdate;

        public JobTimeEntry(string name, double avgSeconds) {
            this.name = name;
            this.avgSeconds = avgSeconds;
            var margin = 1 / Math.Log10(avgSeconds + 1);
            this.dueMargin = margin;
            this.lastupdate = DateTime.Now;
        }
    }

    internal class JobDeadlineMonitor {
        private string path;
        private bool shouldUpdateTimes;
        private Dictionary<string, JobTimeEntry> overview;

        public JobDeadlineMonitor(string path, bool shouldUpdateTimes = false) {

            this.path = path;
            this.shouldUpdateTimes = shouldUpdateTimes;
            this.overview = new Dictionary<string, JobTimeEntry>();
            var filepath = "TimeRecords.json";
            if ( File.Exists(filepath) ) {
                var content = File.ReadAllText(filepath);
                this.overview = JsonConvert.DeserializeObject<Dictionary<string, JobTimeEntry>>(content);
            } else {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine($"No Benchmark file detected for {filepath}");
                Console.WriteLine("Continuing without it, be sure to manually update missing entries!");
                Console.ForegroundColor = ConsoleColor.White;
            }
        }

        public bool JobExists(Job job) {
            return this.JobExists(job.Name);
        }

        public bool JobExists(string name) {
            return this.overview.ContainsKey(name);
        }

        public bool Overdue(Job job) {
            TimeSpan span = job.LatestDeployment.RunTime;
            return this.Overdue(job.Name, span.TotalSeconds);
        }

        public bool Overdue(string name, double currentSeconds) {
            if ( this.overview.TryGetValue(name, out var result) ) {
                return (currentSeconds - result.dueMargin * result.avgSeconds) < result.avgSeconds;
            }
            return false;
        }

        public void UpdateEntry(Job job) {
            TimeSpan span = job.LatestDeployment.RunTime;
            this.UpdateEntry(job.Name, span.TotalSeconds);
        }

        public void UpdateEntry(string name, double seconds) {
            this.overview[name] = new JobTimeEntry(name, seconds);
        }

        public double GetAvgTime(Job job) {
            return GetAvgTime(job.Name);
        }

        public double GetAvgTime(string name) {
            if ( this.overview.TryGetValue(name, out var time) ) {
                return time.avgSeconds;
            }
            return 0;
        }

        public double GetPercentTime(Job job) {
            return this.GetPercentTime(job.Name, job.LatestDeployment.RunTime.TotalSeconds);
        }

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
        }
    }
}
