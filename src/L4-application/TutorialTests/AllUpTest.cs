﻿/* =======================================================================
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
using System.Collections.Generic;
using System.IO;
using System.Linq;
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
            
            if (killBatch) {
                Console.WriteLine("Must ... finish ... ...  MiniBatchProcessor ... ");
                Console.Out.Flush();

                // try to terminate batch processor, if still running:
                int timeoucount = 0;
                while (MiniBatchProcessor.Server.GetIsRunning(null)) {
                    Console.WriteLine("Terminating MiniBatchProcessor...");
                    MiniBatchProcessor.Server.SendTerminationSignal(TimeOutInSeconds: -1);
                    if(timeoucount > 0)
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
        static public bool OneTimeSetUp() {
            /*
            Console.WriteLine("OneTimeSetup: starting 'MiniBatchProcessor'...");
            bool r = MiniBatchProcessor.Server.StartIfNotRunning(RunExternal: false, Reset: true);
            if(r)
                Console.WriteLine("started within this process.");
            else
                Console.WriteLine("already running.");
            
            return r;
            */
            return false;
        }

        internal static string DirectoryOffset = "";

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("quickStartCNS/IsentropicVortex.ipynb")]
        [Test]
        static public void Run__IsentropicVortex() {
            RunWorksheet("quickStartCNS/IsentropicVortex.ipynb");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("InitialValues/InitialValues.tex")]
        [Test]
        static public void Run__InitialValues() {
            RunWorksheet("InitialValues/InitialValues.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("BoundaryAndInitialData/BoundaryAndInitialData.ipynb")]
        [Test]
        static public void Run__BoundaryAndInitialData() {
            // --test=BoSSS.Application.TutorialTests.AllUpTest.Run__BoundaryAndInitialData
            RunWorksheet("BoundaryAndInitialData/BoundaryAndInitialData.ipynb");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("MetaJobManager/MetaJobManager.ipynb")]
        [Test]
        static public void Run__MetaJobManager() {
            RunWorksheet("MetaJobManager/MetaJobManager.ipynb");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("GridGeneration/GridGeneration.tex")]
        [Test]
        static public void Run__GridGeneration() {
            RunWorksheet("GridGeneration/GridGeneration.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("quickStartIBM/channel.tex")]
        [Test]
        static public void Run__channel() {
            RunWorksheet("quickStartIBM/channel.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("shortTutorialMatlab/tutorialMatlab.tex")]
        [Test]
        static public void Run__tutorialMatlab() {
            RunWorksheet("shortTutorialMatlab/tutorialMatlab.tex");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("ue2Basics/ue2Basics.tex")]
        [Test]
        static public void Run__ue2Basics() {
            RunWorksheet("ue2Basics/ue2Basics.tex");
        }

#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("SpatialOperatorNexpTimeInt/SpatialOperatorNexpTimeInt.ipynb")]
        [Test]
        static public void Run__SpatialOperatorNexpTimeInt() {
            RunWorksheet("SpatialOperatorNexpTimeInt/SpatialOperatorNexpTimeInt.ipynb");
        }
#endif

#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("ue5NumFluxConv/ue5NumFluxConv.ipynb")]
        [Test]
        static public void Run__ue5NumFluxConv() {
            RunWorksheet("ue5NumFluxConv/ue5NumFluxConv.ipynb");
        }
#endif

#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("ue6ScalarConvStability/ue6ScalarConvStability.ipynb")]
        [Test]
        static public void Run__ue6ScalarConvStability() {
            RunWorksheet("ue6ScalarConvStability/ue6ScalarConvStability.ipynb");
        }
#endif

#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial9-SIP/sip.ipynb")]
        [Test]
        static public void Run__sip() {
            RunWorksheet("tutorial9-SIP/sip.ipynb");
        }
#endif

#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("tutorial10-PoissonSystem/Poisson.ipynb")]
        [Test]
        static public void Run__Poisson() {
            // --test=BoSSS.Application.TutorialTests.AllUpTest.Run__Poisson
            RunWorksheet("tutorial10-PoissonSystem/Poisson.ipynb");
        }
#endif

#if !DEBUG

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack( "tutorial11-Stokes/StokesEq.ipynb")]
        [Test]
        static public void Run__StokesEq() {
            RunWorksheet("tutorial11-Stokes/StokesEq.ipynb");
        }
#endif
 
#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("CsharpAndBoSSSpad/CsharpAndBoSSSpad.ipynb")]
        [Test]
        static public void Run__CsharpAndBoSSSpad() {
            // --test=BoSSS.Application.TutorialTests.AllUpTest.Run__CsharpAndBoSSSpad
            RunWorksheet("CsharpAndBoSSSpad/CsharpAndBoSSSpad.ipynb");
        }
#endif

#if !DEBUG
        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("convergenceStudyTutorial/convStudy.ipynb")]
        [Test]
        static public void Run__convStudy() {
            RunWorksheet("convergenceStudyTutorial/convStudy.ipynb");
        }
#endif

        /// <summary>
        /// Runs some worksheet contained in the BoSSS handbook.
        /// </summary>
        static public void RunWorksheet(string TexPartialPath) {

            // locate script
            string TexFileName = TexPartialPath.Split(new[] { '/', '\\' }, StringSplitOptions.RemoveEmptyEntries).Last();
            string FullTexName;
            if (!File.Exists(TexFileName)) {
                FullTexName = LocateFile(TexPartialPath).Single();
            } else {
                FullTexName = TexFileName;
            }

            Assert.IsTrue(File.Exists(FullTexName), "unable to find source: " + FullTexName);

            // start the minibatchprocessor which is used internally
            bool iStartedThisShit = OneTimeSetUp();

            BoSSSpad.Job.UndocumentedSuperHack = true;
            BoSSSpad.ReadEvalPrintLoop.WriteFullExceptionInfo = true;
            
            try {
                // run test:
                string mode;
                if(Path.GetExtension(FullTexName).Equals(".tex", StringComparison.InvariantCultureIgnoreCase))
                    mode = "--texbatch";
                else
                    mode = "--JupyterBatch";

                
                int ErrCount = BoSSS.Application.BoSSSpad.BoSSSpadMain.Main(new string[] { mode, FullTexName });

                Console.WriteLine("TutorialTests.exe: finished '{0}', error count is {1}.", FullTexName, ErrCount);
                Assert.LessOrEqual(ErrCount, 0, "Found " + ErrCount + " errors in worksheet: " + FullTexName + " (negative numbers may indicate file-not-found, etc.).");
                Assert.IsTrue(ErrCount >= 0, "Fatal return code: " + ErrCount + " in worksheet: " + FullTexName + " (negative numbers may indicate file-not-found, etc.).");
            } finally {
                // shutting down the local mini batch processor:
                OneTimeTearDown(iStartedThisShit);
            }
        }


        static string[] LocateFile(string PartialPath) {
            DirectoryInfo repoRoot;
            if(!DirectoryOffset.IsEmptyOrWhite())
                repoRoot = new DirectoryInfo(DirectoryOffset);
            else
                repoRoot = new DirectoryInfo(Directory.GetCurrentDirectory());

            // if we get here, we probably have access to the repository root directory.
            string[] r = LocateFileRecursive("", repoRoot, PartialPath);
            if (r == null || r.Length <= 0) {
                throw new IOException("unable to find file '" + PartialPath + "'");
            }

            return r;
        }


        static string[] LocateFileRecursive(string RelPath, DirectoryInfo absPath, string SomeFileName) {
            List<string> ret = new List<string>();

            string _SomeFileName = "*" + SomeFileName;
            
            foreach (var f in absPath.GetFiles()) {
                string RelName = RelPath + f.Name;

                if (RelName.EndsWith(SomeFileName))
                    ret.Add(f.FullName);
                else if (SomeFileName.WildcardMatch(RelName))
                    ret.Add(f.FullName);
                else if (_SomeFileName.WildcardMatch(RelName))
                    ret.Add(f.FullName);

            }

            foreach (var d in absPath.GetDirectories()) {
                ret.AddRange(LocateFileRecursive(RelPath + d.Name + "/", d, SomeFileName));
            }


            return ret.ToArray();
        }

    }
}
