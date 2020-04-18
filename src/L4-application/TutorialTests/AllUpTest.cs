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

using ilPSP;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.IO;
using System.Threading;

namespace BoSSS.Application.TutorialTests {

    /// <summary>
    /// Runs all the worksheets contained in the BoSSS handbook.
    /// </summary>
    [TestFixture]
    static public class AllUpTest {

        /// <summary>
        /// Finalization.
        /// </summary>
        static public void OneTimeTearDown(bool killBatch) {
            //csMPI.Raw.mpiFinalize();

            if (killBatch) {
                // try to terminate batch processor, if still running:
                int timeoucount = 0;
                while (MiniBatchProcessor.Server.IsRunning) {
                    Console.WriteLine("Terminating MiniBatchProcessor...");
                    MiniBatchProcessor.Server.SendTerminationSignal(TimeOutInSeconds: -1);
                    Thread.Sleep(10000);

                    timeoucount++;
                    if (timeoucount > 100) {
                        Assert.Fail("Unable to kill MiniBatchProcessor - server");
                    }
                }
                Console.WriteLine("MiniBatchProcessor terminated.");
            }
        }

        /// <summary>
        /// Init.
        /// </summary>
        //[OneTimeSetUp]
        static public bool OneTimeSetUp() {
            //BoSSS.Solution.Application.InitMPI(new string[0]);

            return MiniBatchProcessor.Server.StartIfNotRunning(RunExternal: false);

            //string preExistingDb = BoSSS.Application.BoSSSpad.InteractiveShell.GetDefaultDatabaseDir();
            //if (Directory.Exists(preExistingDb)) {
            //    //preExistingDb.Delete(true);
            //    Directory.Delete(preExistingDb, true);
            //}
        }

        //static string DirectoryOffset = Path.Combine("..", "..", "..", "..", "..", "doc", "handbook");
        internal static string DirectoryOffset = "";

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("quickStartCNS/IsentropicVortex.tex")]
        [Test]
        static public void Run__IsentropicVortex() {
            RunWorksheet("IsentropicVortex.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("MetaJobManager/MetaJobManager.tex")]
        [Test]
        static public void Run__MetaJobManager() {
            RunWorksheet("MetaJobManager.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("GridGeneration/GridGeneration.tex")]
        [Test]
        static public void Run__GridGeneration() {
            RunWorksheet("GridGeneration.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("quickStartIBM/channel.tex")]
        [Test]
        static public void Run__channel() {
            RunWorksheet("channel.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("shortTutorialMatlab/tutorialMatlab.tex")]
        [Test]
        static public void Run__tutorialMatlab() {
            RunWorksheet("tutorialMatlab.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial2/uebung2tutorial.tex")]
        [Test]
        static public void Run__uebung2tutorial() {
            RunWorksheet("uebung2tutorial.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial4/tutorial4.tex")]
        //[Test]
        static public void Run__tutorial4() {
            RunWorksheet("tutorial4.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial5/uebung5tutorial.tex")]
        //[Test]
        static public void Run__uebung5tutorial() {
            RunWorksheet("uebung5tutorial.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial6/tutorial6.tex")]
        //[Test]
        static public void Run__tutorial6() {
            RunWorksheet("tutorial6.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial9-SIP/sip.tex")]
        //[Test]
        static public void Run__sip() {
            RunWorksheet("sip.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial10-PoissonSystem/Poisson.tex")]
        [Test]
        static public void Run__Poisson() {
            RunWorksheet("Poisson.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack( "tutorial11-Stokes/StokesEq.tex")]
        [Test]
        static public void Run__StokesEq() {
            RunWorksheet("StokesEq.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack( "CsharpAndBoSSSpad/CsharpAndBoSSSpad.tex")]
        [Test]
        static public void Run__CsharpAndBoSSSpad() {
            RunWorksheet("CsharpAndBoSSSpad.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("convergenceStudyTutorial/convStudy.tex")]
        //[Test]
        static public void Run__convStudy() {
            RunWorksheet("convStudy.tex");
        }

        /// <summary>
        /// Runs some worksheet contained in the BoSSS handbook.
        /// </summary>
        static public void RunWorksheet(string TexFileName) {

            // run test:
            string FullTexName = Path.Combine(DirectoryOffset, TexFileName);
            Assert.IsTrue(File.Exists(FullTexName), "unable to find TeX source: " + FullTexName);

            bool iStartedThisShit = OneTimeSetUp();

            int ErrCount = BoSSS.Application.BoSSSpad.BoSSSpadMain.Main(new string[] { "--texbatch", FullTexName });

            Console.WriteLine("TutorialTests.exe: finished '{0}', error count is {1}.", FullTexName, ErrCount);

            Assert.LessOrEqual(ErrCount, 0, "Found " + ErrCount + " errors in worksheet: " + FullTexName + " (negative numbers may indicate file-not-found, etc.).");
            Assert.IsTrue(ErrCount >= 0, "Fatal return code: " + ErrCount + " in worksheet: " + FullTexName + " (negative numbers may indicate file-not-found, etc.).");

            // shutting down the local mini batch processor:
            OneTimeTearDown(iStartedThisShit);

            //foreach(var db in BoSSS.Application.BoSSSpad.InteractiveShell.databases) {
            //    db.
            //}
        }

    }
}
