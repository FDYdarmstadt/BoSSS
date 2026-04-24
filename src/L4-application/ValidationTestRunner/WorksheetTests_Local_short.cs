using BoSSS.Application.BoSSSpad;
using BoSSS.Application.LsTest;
using BoSSS.Application.TutorialTests;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using FreeXNSE;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using NUnit.Framework.Constraints;
using PublicTestRunner;
using SAIDT;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;
using XESTSF;

namespace ValidationTestRunner {

   

    /// <summary>
    /// NUnit entry point for each example worksheet which represents a short-running validation test;
    /// </summary>
    /// <remarks>
    /// - short running rests are fully re-computed every time
    /// - All these tests here are intended to be run at the local MS windows HPC cluster (aka. FDYcluster) at Chair of Fluid Dynamics (FDY)
    /// </remarks>
    [TestFixture]
    [NUnitNumThreads(1)]
    static public class WorksheetTests_Local_short {

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("BoundaryAndInitialData/BoundaryAndInitialData.ipynb")]
        [Test]
        static public void Run__BoundaryAndInitialData() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__BoundaryAndInitialData
            Mutex JupyterMutex = new Mutex(false, "BoundaryAndInitialData");
            try {
                JupyterMutex.WaitOne();

                WorkflowMgmTestUtils.DeleteDatabase("Demo_BoundaryAndInitialData");
                WorkflowMgmTestUtils.DeleteDeployments("Demo_BoundaryAndInitialData*");
                ValidationTestRunnerMain.RunWorksheet("BoundaryAndInitialData/BoundaryAndInitialData.ipynb");
            } finally {
                JupyterMutex.ReleaseMutex();
            }
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("MetaJobManager/MetaJobManager.ipynb")]
        [Test]
        static public void Run__MetaJobManager() {
            //--test=ValidationTestRunner.WorksheetTests_Local_short.Run__MetaJobManager
            Mutex JupyterMutex = new Mutex(false, "MetaJobManager_Tutorial");
            try {
                // 0.) Clean leftover of previous runs
                // -----------------------------------
                JupyterMutex.WaitOne();
                WorkflowMgmTestUtils.DeleteDatabase("MetaJobManager_Tutorial");
                WorkflowMgmTestUtils.DeleteDeployments("MetaJobManager_Tutorial*");

                var allDepl = WorkflowMgmTestUtils.GetAllDeployments("MetaJobManager_Tutorial*");
                Assert.Zero(allDepl.Length, $"expecting 0 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));

                // 1.) First run
                // -----------------------------------

                ValidationTestRunnerMain.RunWorksheet("MetaJobManager/MetaJobManager.ipynb");

                BoSSSshell.WorkflowMgm.Init("MetaJobManager_Tutorial");

                Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed");
                Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, "MetaJobManager_Tutorial", "un-expected project name");
                allDepl = WorkflowMgmTestUtils.GetAllDeployments("MetaJobManager_Tutorial*");
                Assert.AreEqual(allDepl.Length, 1, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));
                Console.WriteLine("successfully finished first run.");

                /*
                // 2.) Persistence test: Test that, if the Worksheet is already executed, the job will **not** be re-submitted
                // --------------------------------------------------------------------------------------------------------------------------
                // execute the worksheet a second time!
                ValidationTestRunnerMain.RunWorksheet("MetaJobManager/MetaJobManager.ipynb");

                allDepl = WorkflowMgmTestUtils.GetAllDeployments("MetaJobManager_Tutorial*");
                Assert.AreEqual(allDepl.Length, 1, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));
                BoSSSshell.WorkflowMgm.ResetSessionsCache();
                Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed twice");
                Console.WriteLine("successfully finished second run.");

                // 3.) delete deployments and test a third time:
                // -----------------------------------------------
                WorkflowMgmTestUtils.DeleteDeployments("MetaJobManager_Tutorial*");
                ValidationTestRunnerMain.RunWorksheet("MetaJobManager/MetaJobManager.ipynb");

                BoSSSshell.WorkflowMgm.ResetSessionsCache();
                Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed twice");
                allDepl = WorkflowMgmTestUtils.GetAllDeployments("MetaJobManager_Tutorial*");
                Assert.Zero(allDepl.Length, $"expecting 0 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));
                Console.WriteLine("successfully finished third run.");
                */


            } finally {
                JupyterMutex.ReleaseMutex();
            }
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("GridGeneration/GridGeneration.ipynb")]
        [Test]
        static public void Run__GridGeneration() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__GridGeneration
            ValidationTestRunnerMain.RunWorksheet("GridGeneration/GridGeneration.ipynb");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("quickStartIBM/channel.ipynb")]
        [Test]
        static public void Run__channel() {
            ValidationTestRunnerMain.RunWorksheet("quickStartIBM/channel.ipynb");
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("convergenceStudyTutorial/convStudy.ipynb")]
        [Test]
        static public void Run__convStudy() {
            Mutex JupyterMutex = new Mutex(false, "ConvStudyTutorial");
            try {
                JupyterMutex.WaitOne();
                WorkflowMgmTestUtils.DeleteDatabase("ConvStudyTutorial");
                WorkflowMgmTestUtils.DeleteDeployments("ConvStudyTutorial*");
                ValidationTestRunnerMain.RunWorksheet("convergenceStudyTutorial/convStudy.ipynb");
            } finally {
                JupyterMutex.ReleaseMutex();
            }
        }


    }
}
