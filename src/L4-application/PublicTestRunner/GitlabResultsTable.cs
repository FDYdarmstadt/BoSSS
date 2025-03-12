using BoSSS.Application.BoSSSpad;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml;

namespace PublicTestRunner {
    internal class GitlabResultsTable {
        public static void WriteResultsTable(List<Job> jobs, string testsuite) {
            var accumulatedTime = jobs.Select(j=>j.LatestDeployment.RunTime.TotalSeconds).Sum();
            var noTests = jobs.Count;
            var noFailed = jobs.Where(j => j.LatestDeployment.Status != JobStatus.FinishedSuccessful).Count();

            using ( XmlWriter writer = XmlWriter.Create("TestResults.xml") ) {
                writer.WriteStartDocument();
                writer.WriteStartElement("testsuites");
                writer.WriteAttributeString("id", testsuite + DateTime.Now);
                writer.WriteAttributeString("name", testsuite);
                writer.WriteAttributeString("tests", noTests.ToString());
                writer.WriteAttributeString("failures", noFailed.ToString());
                writer.WriteAttributeString("time", accumulatedTime.ToString());
                {
                    writer.WriteStartElement("testsuite");
                    writer.WriteAttributeString("id", testsuite + DateTime.Now);
                    writer.WriteAttributeString("name", testsuite);
                    writer.WriteAttributeString("tests", noTests.ToString());
                    writer.WriteAttributeString("failures", noFailed.ToString());
                    writer.WriteAttributeString("time", accumulatedTime.ToString());

                    foreach ( Job job in jobs ) {
                        writer.WriteStartElement("testcase");
                        writer.WriteAttributeString("id", job.LatestDeployment.BatchProcessorIdentifierToken);
                        writer.WriteAttributeString("name", job.Name);
                        writer.WriteAttributeString("time", job.LatestDeployment.RunTime.TotalSeconds.ToString());
                        // Failure
                        if ( job.LatestDeployment.Status != JobStatus.FinishedSuccessful ) {
                            writer.WriteStartElement("failure");
                            writer.WriteAttributeString("message", job.LatestDeployment.DeploymentDirectory.ToString());
                            writer.WriteAttributeString("type", job.LatestDeployment.Status.ToString());
                            writer.WriteElementString("system-out", job.Stdout);
                            writer.WriteElementString("system-err", job.Stderr);
                            writer.WriteEndElement();
                        }
                        writer.WriteEndElement(); // end testcase
                    }
                    writer.WriteEndElement(); // End testsuite
                }
                writer.WriteEndElement(); // End testsuites
                writer.WriteEndDocument();
            }
        }
    }
}
