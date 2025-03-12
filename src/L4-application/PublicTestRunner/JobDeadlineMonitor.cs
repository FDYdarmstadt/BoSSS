using BoSSS.Application.BoSSSpad;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;

namespace PublicTestRunner {
    internal class JobTimeEntry {
        public string name;
        public double avgSeconds;
        public double dueMargin;
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
            var margin = 1 / Math.Log10(seconds + 1);
            this.overview.Add(name, new JobTimeEntry { name = name, avgSeconds = seconds, dueMargin = margin });
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
            // var content = JsonConvert.SerializeObject(this.overview, Formatting.Indented);
            // var filepath = this.shouldUpdateTimes ? "TimeRecords.json" : "TimeRecordsTmp.json";
            // Console.ForegroundColor = ConsoleColor.Red;
            // Console.WriteLine($"Saving Timings to: {Path.GetFullPath(filepath)}");
            // Console.ForegroundColor = ConsoleColor.White;
            // File.WriteAllText(filepath, content);
            using ( XmlWriter writer = XmlWriter.Create("TestResults.xml") ) {
                writer.WriteStartDocument();
                writer.WriteStartElement("testsuites");
                writer.WriteAttributeString("id", DateTime.Now.ToString());
                writer.WriteAttributeString("name", "Testsuites");
                writer.WriteAttributeString("tests", "250");
                writer.WriteAttributeString("failures", "10");
                writer.WriteAttributeString("time", "400.80");
                // Start children testsuite
                for ( int i = 0; i < 2; i++ ) {
                    writer.WriteStartElement("testsuite");
                    writer.WriteAttributeString("id", DateTime.Now.ToString());
                    writer.WriteAttributeString("name", "Testsuite" + i);
                    writer.WriteAttributeString("tests", "150");
                    writer.WriteAttributeString("failures", "5");
                    writer.WriteAttributeString("time", "200.80");
                    // Start Testcase
                    for ( int j = 0; j < 10; j++ ) {
                        writer.WriteStartElement("testcase");
                        writer.WriteAttributeString("id", DateTime.Now.ToString());
                        writer.WriteAttributeString("name", "Testcase" + j);
                        writer.WriteAttributeString("time", "10");
                        // Failure
                        if ( j % 2 == 0 ) {
                            writer.WriteStartElement("failure");
                            writer.WriteAttributeString("message", "Well it failed");
                            writer.WriteAttributeString("type", "Profiling");
                            writer.WriteEndElement();

                            if ( j % 4 == 0 ) {
                                writer.WriteStartElement("failure");
                                writer.WriteAttributeString("message", "two times ?");
                                writer.WriteAttributeString("type", "Profiling");
                                writer.WriteEndElement();// end failure
                            }
                        }

                        writer.WriteEndElement(); // end testcase
                    }
                    writer.WriteEndElement(); // End Testsuite1
                }
                writer.WriteEndElement(); // end testsuites
                writer.WriteEndDocument();
            }
        }
    }
}
