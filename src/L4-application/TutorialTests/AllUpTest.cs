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

namespace BoSSS.Application.TutorialTests {

    /// <summary>
    /// Runs all the worksheets contained in the BoSSS handbook.
    /// </summary>
    [TestFixture]
    static public class AllUpTest {

        /// <summary>
        /// MPI finalization.
        /// </summary>
        [TestFixtureTearDown]
        static public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// MPI init.
        /// </summary>
        [TestFixtureSetUp]
        static public void TestFixtureSetUp() {
            BoSSS.Solution.Application.InitMPI(new string[0]);

             if (System.Environment.MachineName.ToLowerInvariant().EndsWith("rennmaschin")
                //|| System.Environment.MachineName.ToLowerInvariant().Contains("jenkins")
                ) {
                // This is Florians Laptop;
                // he is to poor to afford MATLAB, so he uses OCTAVE
                BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;
                BatchmodeConnector.MatlabExecuteable = "C:\\cygwin64\\bin\\bash.exe";
            } 
        }

        static string DirectoryOffset = Path.Combine("..", "..", "..", "..", "..", "doc", "handbook");

        /// <summary>
        /// Runs all the worksheets contained in the BoSSS handbook.
        /// </summary>
        [Test]
        static public void RunWorksheets([Values(
            "quickStartCNS/IsentropicVortex.tex",
            "MetaJobManager/MetaJobManager.tex",
            "GridGeneration/GridGeneration.tex",
            "quickStartIBM/channel.tex",
            "shortTutorialMatlab/tutorialMatlab.tex",
            // ----
            "tutorial2/uebung2tutorial.tex",
            "tutorial4/tutorial4.tex",
            "tutorial5/uebung5tutorial.tex",
            "tutorial6/tutorial6.tex",
            "tutorial9-SIP/sip.tex",
            // ---
            "tutorial10-PoissonSystem/Poisson.tex",
            "tutorial11-Stokes/StokesEq.tex",
            "CsharpAndBoSSSpad/CsharpAndBoSSSpad.tex"//   */
            "ParameterStudy\ParameterStudy.tex"
            )] string TexFileName) {

            string FullTexName = Path.Combine(DirectoryOffset, TexFileName);
            Assert.IsTrue(File.Exists(FullTexName), "unable to find TeX source: " + FullTexName);

            int ErrCount = BoSSS.Application.BoSSSpad.BoSSSpadMain.Main(new string[] { "--texbatch", FullTexName });

            Console.WriteLine("TutorialTests.exe: finished '{0}', error count is {1}.", FullTexName, ErrCount);

            Assert.LessOrEqual(ErrCount, 0, "Found " + ErrCount + " errors in worksheet: " + FullTexName + " (negative numbers may indicate file-not-found, etc.).");
            Assert.IsTrue(ErrCount >= 0, "Fatal return code: " + ErrCount + " in worksheet: " + FullTexName + " (negative numbers may indicate file-not-found, etc.).");
        }

    }
}
