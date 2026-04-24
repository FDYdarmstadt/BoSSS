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
    /// NUnit entry point for each example worksheet which represents a long-term validation test
    /// </summary>
    /// <remarks>
    /// All these tests here are intended to be run at the Lichtenberg HPC;
    /// The worksheets itself run on a server at FDY, only compute jobs are send to Lichtenberg.
    /// 
    /// On this local FDY server, 
    /// we are currently using a <see cref="MiniBatchProcessorClient"/> in order to execute the worksheets.
    /// (Note that the tests/worksheets are actually executed on our local server;
    /// the compute jobs which these worksheets run, are deployed to Lichtenberg.)
    /// In order to be able to run tests/worksheets in parallel, we use the 
    /// `runjobmanager` option for the test runner, which submits to the local mini batch processor.
    /// 
    /// This is heavily hard-coded to our test server environment and subject to change at some point.
    /// </remarks>
    [TestFixture]
    static public class WorksheetTests_Lichtenberg {

        /// <summary>
        /// Linear solver performance:
        /// - Steady-State XDG Stokes problem (water droplet in air)
        /// - one MPI core
        /// </summary>
        [NUnitFileToCopyHack(
            "handbook/apdx-MPISolverPerformance/unified/ParLinslvPerf_ConstPoisson.ipynb",
            "handbook/apdx-MPISolverPerformance/unified/ParLinslvPerf_XdgPoisson.ipynb",
            "handbook/apdx-MPISolverPerformance/unified/ParLinslvPerf_BottiPietroStokes3D.ipynb",
            "handbook/apdx-MPISolverPerformance/unified/ParLinslvPerf_XdgStokes.ipynb",
            "handbook/apdx-MPISolverPerformance/unified/ParLinslvPerf_Evaluation.ipynb")]
        [Test]
        static public void Run__LinslvPerfPar() {
            // --test=ValidationTestRunner.WorksheetTests_Lichtenberg.Run__LinslvPerfPar

            string PROJECT_NAME = System.Environment.GetEnvironmentVariable("LinslvPerfPar") ?? "LinslvPerfPar"; // this allows to modify the project name for testing purposes

            Console.WriteLine("skipping Run__LinslvPerfSer");
            return;

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                PROJECT_NAME,
                new TimeSpan(days: 40, hours: 0, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_ConstPoisson.ipynb", allowErrors: true);
            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_XdgPoisson.ipynb", allowErrors: true);
            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_BottiPietroStokes3D.ipynb", allowErrors: true);
            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_XdgStokes.ipynb", allowErrors: true);

            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_Evaluation.ipynb", allowErrors: false);

            // if we reach this point, its a success - backup the datebase
            ValidationTestRunnerMain.BackupDatabase(PROJECT_NAME);
        }

    }



    /// <summary>
    /// NUnit entry point for each example worksheet which validates that 
    /// the basic functionality of the Lichtenberg HPC somehow works
    /// (i.e. access granted, paths correct, submission to slurm works)
    /// </summary>
    /// <remarks>
    /// All these tests here are intended to be run at the Lichtenberg HPC;
    /// The worksheets itself run on a server at FDY, only compute jobs are send to Lichtenberg.
    /// </remarks>
    [TestFixture]
    [NUnitNumThreads(1)]

    static public class WorksheetTests_LichtenbergItselfWorking {

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("ue2Basics/ue2Basics.ipynb")]
        [Test]
        static public void Run__ue2Basics() {
            //--test=ValidationTestRunner.WorksheetTests_LichtenbergItselfWorking.Run__ue2Basics
            ValidationTestRunnerMain.RunWorksheet("ue2Basics/ue2Basics.ipynb");

            try {
                File.Move("ue2Basics.ipynb", "ue2Basics-Lichtenberg.html");
            } catch(Exception e) {
                Console.Error.WriteLine($"File copy exception: {e.GetType()} : {e.Message}");
            }
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("MetaJobManager/MetaJobManager.ipynb")]
        [Test]
        static public void Run__MetaJobManager() {
            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "MetaJobManager_Tutorial",
                new TimeSpan(days: 5, hours: 1, minutes: 0, seconds: 1));


            Console.WriteLine("Testing MetaJobManager @ Lichtenberg...");
            ValidationTestRunnerMain.RunWorksheet("MetaJobManager/MetaJobManager.ipynb");
            Console.WriteLine("Finished MetaJobManager @ Lichtenberg.");

            // if we reach this point, its a success - backup the datebase
            ValidationTestRunnerMain.BackupDatabase("MetaJobManager_Tutorial");
        }
    }
}
