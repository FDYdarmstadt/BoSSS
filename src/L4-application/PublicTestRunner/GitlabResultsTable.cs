using BoSSS.Application.BoSSSpad;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;

namespace PublicTestRunner {
    internal class GitlabResultsTable {
        public static void WriteResultsTable(List<Job> jobs, string testsuite) {

            bool IsValidXmlChar(char c) {
                // XML 1.0 valid characters:
                // #x9 | #xA | #xD | [#x20-#xD7FF] | [#xE000-#xFFFD] | [#x10000-#x10FFFF]
                return c == 0x9 || c == 0xA || c == 0xD ||
                    (c >= 0x20 && c <= 0xD7FF) ||
                    (c >= 0xE000 && c <= 0xFFFD);
            }


            string SanitizeXmlString(string text) {
                if ( string.IsNullOrEmpty(text) )
                    return text;

                var sb = new StringBuilder(text.Length);
                foreach ( char c in text ) {
                    if ( IsValidXmlChar(c) )
                        sb.Append(c);
                    // Optionally replace with space or other character:
                    // else sb.Append(' ');
                }
                return sb.ToString();
            }



            try {
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
                                writer.WriteElementString("system-out", SanitizeXmlString(job.Stdout));
                                writer.WriteElementString("system-err", SanitizeXmlString(job.Stderr));
                                writer.WriteEndElement();
                            }
                            writer.WriteEndElement(); // end testcase
                        }
                        writer.WriteEndElement(); // End testsuite
                    }
                    writer.WriteEndElement(); // End testsuites
                    writer.WriteEndDocument();
                }
            } catch ( Exception e ) {
                Console.Error.WriteLine("Exception writing gitlab table: " + e);
            }
        }
    }
}
