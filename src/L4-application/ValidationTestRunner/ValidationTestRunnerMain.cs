using ilPSP.Utils;
using PublicTestRunner;
using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using ilPSP;
using BoSSS.Application.BoSSSpad;
using System.Collections.Generic;
using System.Diagnostics;
using SAIDT;
using BoSSS.Application.TutorialTests;
using System.Threading;

namespace ValidationTestRunner {

    /// <summary>
    /// Extends the public tests to some which are only available in the internal par of BoSSS
    /// </summary>
    class ValidationTests : ITestTypeProvider {

        public Type[] FullTest {
            get {
                var ret = new Type[] {
                    typeof(ValidationTestRunnerMain),
                    typeof(SAIDTMain) // required to have the SAIDT binary available
                };
                return ret;
            }
        }

        public Type[] ReleaseOnlyTests {
            get {
                var ret = new Type[]{

                };
                return ret;
            }
        }

        /// <summary>
        /// MPI, DEBUG and RELEASE
        /// </summary>
        public (Type type, int NoOfProcs)[] MpiFullTests {
            get {
                var ret = new (Type, int)[] { };
                return ret;
            }
        }

        /// <summary>
        /// MPI, only RELEASE
        /// </summary>
        public (Type type, int NoOfProcs)[] MpiReleaseOnlyTests {
            get {
                var ret = new (Type, int)[] { };
                // e.g.
                //(typeof(XDGShock.Program), 4).AddToArray(ref ret); // 2nd entry is number of MPI cores
                //(typeof(XDGShock_MPITests.XDGShock_MPITestsMain), 4).AddToArray(ref ret);
                return ret;
            }
        }

        public DirectoryInfo GetRepositoryBaseDir() {
            DirectoryInfo repoRoot;
            try {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                repoRoot = dir.Parent.Parent.Parent.Parent.Parent.Parent.Parent;

                var publicDir = repoRoot.GetDirectories("public").SingleOrDefault();
                if(publicDir == null)
                    return null;

                var src = publicDir.GetDirectories("src").SingleOrDefault();
                var libs = publicDir.GetDirectories("libs").SingleOrDefault();
                var doc = publicDir.GetDirectories("doc").SingleOrDefault();

                if(src == null || !src.Exists)
                    return null;
                //throw new Exception();
                if(libs == null || !libs.Exists)
                    return null;
                //throw new Exception();
                if(doc == null || !doc.Exists)
                    return null;
                //throw new Exception();

            } catch(Exception) {
                return null;
                //throw new IOException("Unable to find repository root. 'runjobmanger' must be invoked from its default location within the BoSSS git repository.");
            }

            return repoRoot;
        }

        virtual public bool CopyManagedAssembliesCentrally => false;

        virtual public int RetryCount => 2;

        virtual public bool DeleteSuccessfulTestFiles => false;
    }


    /// <summary>
    /// NUnit entry point for each example worksheet which represents a long-term validation test
    /// </summary>
    /// <remarks>
    /// - long-term tests are typicalle executed from some backup database; therefore, the file `BOSSS_RUNTESTFROMBACKUP.txt` must be present in the local dir
    /// - All these tests here are intended to be run at the local MS windows HPC cluster (aka. FDYcluster) at Chair of Fluid Dynamics (FDY)
    /// </remarks>
    [TestFixture]
    [NUnitNumThreads(1)]
    static public class WorksheetTests_Local_long {


        /// <summary>
        /// XDG-IST Solver, 
        /// publication results for: Vandergrift, Kummer: An extended discontinuous Galerkin shock tracking method, https://onlinelibrary.wiley.com/doi/full/10.1002/fld.5293
        /// </summary>
        [NUnitFileToCopyHack("ShockFitting/Studies/ConvergenceStudy/ConvergenceStudy_BowShock_HPC.ipynb", "ShockFitting/Studies/ConvergenceStudy/bosss_db_levelSets.zip", "ShockFitting/Studies/ConvergenceStudy/BowShockPoints.txt", "ShockFitting/Studies/ConvergenceStudy/ConvergenceStudy_BowShock_PostProcessing.ipynb")]
        [Test]
        static public void Run__XDGIST_BowShock()
        {
            // delete the database if it is more than 75 days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "XESF_BowShock_ConvStudy",
                "XESF_BowShock_ConvStudy",
                "DELETE_XDGISTBowShock",
                new TimeSpan(days: 75, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("ShockFitting/Studies/ConvergenceStudy/ConvergenceStudy_BowShock_HPC.ipynb");
            ValidationTestRunnerMain.RunWorksheet("ShockFitting/Studies/ConvergenceStudy/ConvergenceStudy_BowShock_PostProcessing.ipynb");

            Console.WriteLine("XDGISTBowShock @ FDYcluster");
        }

        /// <summary>
        /// XDG-IST Solver, 
        /// thesis results for: Vandergrift: Implicit Discontinuous Galerkin Shock Tracking Methods for Compressible Flows with Shocks (2024)
        /// </summary>
        [NUnitFileToCopyHack("ShockFitting/Studies/ConvergenceStudy/AcousticWave1D_ConvergenceStudy.ipynb", "ShockFitting/Studies/ConvergenceStudy/AcousticWave1D_ConvergenceStudy_PostProcessing.ipynb")]
        [Test]
        static public void Run__XDGIST_1DShockAcoustic() {

            // delete the database if it is more than 25 days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "XESTSF_ShockAcousticInteraction1D_ConvergenceStudy",
                "XESTSF_ShockAcousticInteraction1D_ConvergenceStudy",
                "DELETE_XESTSFShockAcousticInteraction1D",
                new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("ShockFitting/Studies/ConvergenceStudy/AcousticWave1D_ConvergenceStudy.ipynb");
            ValidationTestRunnerMain.RunWorksheet("ShockFitting/Studies/ConvergenceStudy/AcousticWave1D_ConvergenceStudy_PostProcessing.ipynb");

            Console.WriteLine("XDGIST1DShockAcoustic @ FDYcluster");
        }

        /// <summary>
        /// CNS Solver, 
        /// thesis results for: Vandergrift: Implicit Discontinuous Galerkin Shock Tracking Methods for Compressible Flows with Shocks (2024)
        /// </summary>
        [NUnitFileToCopyHack("ShockFitting/Studies/ConvergenceStudy/CNSAcousticWave1DHPC_ConvStudy.ipynb", "ShockFitting/Studies/ConvergenceStudy/CNSAcousticWave1DHPC_ConvStudy_PostProcessing.ipynb")]
        [Test]
        static public void Run__CNS_1DShockAcoustic()
        {
            // delete the database if it is more than 25 days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "CNS_AcousticWave1D_ConvStudy",
                "CNS_AcousticWave1D_ConvStudy",
                "DELETE_CNSShockAcousticInteraction1D",
                new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("ShockFitting/Studies/ConvergenceStudy/CNSAcousticWave1DHPC_ConvStudy.ipynb");
            ValidationTestRunnerMain.RunWorksheet("ShockFitting/Studies/ConvergenceStudy/CNSAcousticWave1DHPC_ConvStudy_PostProcessing.ipynb");

            Console.WriteLine("CNS1DShockAcoustic @ FDYcluster");
        }

        /// <summary>
        /// Rheology Solver, 
        /// publication results for: Kikker, Kummer, Oberlack: A fully coupled high-order discontinuous Galerkin solver for viscoelastic fluid flow, https://onlinelibrary.wiley.com/doi/10.1002/fld.4950
        /// </summary>
        [NUnitFileToCopyHack("rheology/ConfinedCylinder_ConvergenceStudy.ipynb", "rheology/ConfinedCylinder_ConvergenceStudy_Postprocessing.ipynb", "rheology/mesh_karman_OriginalBox_MEDIUM_*_half.msh")]
        [Test]
        static public void Run__RheologyConfinedCylinder() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__RheologyConfinedCylinder

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "ConfinedCylinder_ConvergenceStudy",
                "ConfinedCylinder_ConvergenceStudy*",
                "DELETE_RheologyConfinedCylinder",
                new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("rheology/ConfinedCylinder_ConvergenceStudy.ipynb");
            ValidationTestRunnerMain.RunWorksheet("rheology/ConfinedCylinder_ConvergenceStudy_Postprocessing.ipynb");

            Console.WriteLine("RheologyConfinedCylinder @ FDYcluster");
        }


        /// <summary>
        /// Hagen-Poiseulle flow (aka. pipe flow) for the helical symmetric solver
        /// Maintainer: Schahin Akbari
        /// </summary>
        [NUnitFileToCopyHack("HelicalSymmetricSolver/HagenPoiseulle.ipynb", "HelicalSymmetricSolver/Post_Processing_HagenPoiseulle.ipynb")]
        [Test]
        static public void Run__Helical_HagenPoiseulle() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__Helical_HagenPoiseulle

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "Helical_HagenPoiseulle",
                "Helical_HagenPoiseulle*",
                "DELETE_Helical_HagenPoiseulle",
                new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("HelicalSymmetricSolver/HagenPoiseulle.ipynb");
            ValidationTestRunnerMain.RunWorksheet("HelicalSymmetricSolver/Post_Processing_HagenPoiseulle.ipynb");

            Console.WriteLine("Helical_HagenPoiseulle @ FDYcluster");
        }


        /// <summary>
        /// Centrifugal flow (aka. centrifugal flow) for the helical symmetric solver
        /// Maintainer: Schahin Akbari
        /// </summary>
        [NUnitFileToCopyHack("HelicalSymmetricSolver/Centrifugal.ipynb", "HelicalSymmetricSolver/Post_Processing_Centrifugal.ipynb")]
        [Test]
        static public void Run__Helical_Centrifugal() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__Helical_Centrifugal

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "Helical_Centrifugal",
                "Helical_Centrifugal*",
                "DELETE_Helical_Centrifugal",
                new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("HelicalSymmetricSolver/Centrifugal.ipynb");
            ValidationTestRunnerMain.RunWorksheet("HelicalSymmetricSolver/Post_Processing_Centrifugal.ipynb");

            Console.WriteLine("Helical_Centrifugal @ FDYcluster");
        }


        /// <summary>
        /// Contact Line at heated wall,
        /// Maintainer: Matthias Rieckmann
        /// </summary>
        [NUnitFileToCopyHack("XNSFE_Solver/HeatedWall_Validation/HeatedWallSimple_VerificationFastMarching.ipynb", "XNSFE_Solver/HeatedWall_Validation/*.json")]
        [Test]
        static public void Run__HeatedWallSimple() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__HeatedWallSimple
            string really = System.Environment.GetEnvironmentVariable("RUN_HEATEDWALLSIMPLE");
            if (really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping Run__HeatedWallSimple ");
                return;
            } else {
                Console.WriteLine("RUN_HEATEDWALLSIMPLE = " + really);
            }
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "HeatedWallSimple_VerificationFastMarching",
                "HeatedWall_Simple*",
                "delete_HeatedWallSimple",
                new TimeSpan(days: 10, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("XNSFE_Solver/HeatedWall_Validation/HeatedWallSimple_VerificationFastMarching.ipynb");
            Console.WriteLine("HeatedWallSimple @ FDYcluster");
        }

        /// <summary>
        /// Contact Line at heated wall,
        /// Maintainer: Matthias Rieckmann
        /// </summary>
        [NUnitFileToCopyHack("XNSFE_Solver/HeatedWall_Validation/HeatedWallConvergenceValidation_*.ipynb", "XNSFE_Solver/HeatedWall_Validation/HeatedWall_Validation.zip")]
        [Test]
        static public void Run__HeatedWallConvergence() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__HeatedWallConvergence

            string really = System.Environment.GetEnvironmentVariable("RUN_HEATEDWALLCONVERGENCE");
            if (really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping Run__HeatedWallConvergence ");
                return;
            } else {
                Console.WriteLine("RUN_HEATEDWALLCONVERGENCE = " + really);
            }
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "HeatedWallSimple_VerificationFastMarching",
                "HeatedWall_Simple*",
                "delete_HeatedWallConvergence",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("XNSFE_Solver/HeatedWall_Validation/HeatedWallConvergenceValidation_Controls.ipynb");
            ValidationTestRunnerMain.RunWorksheet("XNSFE_Solver/HeatedWall_Validation/HeatedWallConvergenceValidation_Postprocessing.ipynb");
            ValidationTestRunnerMain.RunWorksheet("XNSFE_Solver/HeatedWall_Validation/HeatedWallConvergenceValidation_Comparison.ipynb");
            Console.WriteLine("HeatedWallConvergence @ FDYcluster");
        }

        /// <summary>
        /// Printing Nip Stokes Simulations
        /// </summary>
        [NUnitFileToCopyHack("PrintingNip/*.ipynb", "PrintingNip/*.sh", "PrintingNip/*.tex", "PrintingNip/*.txt")]
        [Test]
        static public void Run__PrintingNip() {

            // additional artifacts in ./PrintingNip/Figures (raw gnuplot figures), ./PrintingNip/Output (preview pdf with figures), ./PrintingNip/Files (datatables in csv format) can be generated when also running "Part0" and "Part6". This is not part of testing!
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__PrintingNip

            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "PrintingNip_Part1",
                "PrintingNip_Part1*",
                "delete_PRINTINGNIP",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 0));
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "PrintingNip_Part2",
                "PrintingNip_Part2*",
                "delete_PRINTINGNIP",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 0));
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "PrintingNip_Part3",
                "PrintingNip_Part3*",
                "delete_PRINTINGNIP",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 0));
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "PrintingNip_Part4",
                "PrintingNip_Part4*",
                "delete_PRINTINGNIP",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 0));

            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part1_PrintingNip_Correlation_Run.ipynb");
            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part1_PrintingNip_Correlation_Evaluate.ipynb");

            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part2_PrintingNip_ConstantStagnationPoint_Run.ipynb");
            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part2_PrintingNip_ConstantStagnationPoint_Evaluate.ipynb");

            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part3_PrintingNip_SimulateExperiment_Run.ipynb");
            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part3_PrintingNip_SimulateExperiment_Evaluate.ipynb");

            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part4_PrintingNip_SimulateUXMap_Run.ipynb");

            ValidationTestRunnerMain.RunWorksheet("PrintingNip/Part5_PrintingNip_Validation.ipynb");

            Console.WriteLine("PrintingNip @ FDYcluster");
        }

        /// <summary>
        /// Test of the Low-Mach solver;
        /// </summary>
        [NUnitFileToCopyHack("LowMach/HeatedCouetteFlow/CouetteTemperatureDifference_ConvStudy.ipynb"
            )]
        [Test]
        static public void Run__CouetteTemperatureDifference_ConvStudy() {
            Console.WriteLine("HeatedCouette @ FDYcluster");

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "CouetteTemperatureDifference_ConvStudy",
                "CouetteTemperatureDifference_ConvStudyca*",
                "delete_CouetteTemperatureDifference_ConvStudy",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedCouetteFlow/CouetteTemperatureDifference_ConvStudy.ipynb");

        }


        /// <summary>
        /// Test of the Low-Mach solver;
        /// Publication Results from: Gutierrez and Kummer, 2021, A fully coupled high-order Discontinuous Galerkin method for diffusion flames in a low-Mach number framework
        /// </summary>
        [NUnitFileToCopyHack("LowMach/HeatedSquareCavity/HeatedCavity_RaSweep.ipynb",
                             "LowMach/HeatedSquareCavity/HeatedCavity_RaSweepPostProc.ipynb",
                              "LowMach/HeatedSquareCavity/*.txt" //
            )]
        [Test]
        static public void Run__HeatedCavityRayleighSweep() {
            //--test=ValidationTestRunner.WorksheetTests_Local.Run__HeatedCavityRayleighSweep 
            Console.WriteLine("HeatedCavity_RaSweep @ FDYcluster");

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "HeatedCavity_RayleighSweepStudy",
                "HeatedCavity_RayleighSweepStudy*",
                "delete_HeatedCavityRayleighSweep",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedSquareCavity/HeatedCavity_RaSweep.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedSquareCavity/HeatedCavity_RaSweepPostProc.ipynb");

        }
        /// <summary>
        /// Test of the Low-Mach solver;
        /// Publication Results from: Gutierrez and Kummer, 2021, A fully coupled high-order Discontinuous Galerkin method for diffusion flames in a low-Mach number framework
        /// </summary>
        [NUnitFileToCopyHack("LowMach/HeatedSquareCavity/HeatedCavity_ConvStudy.ipynb", "LowMach/HeatedSquareCavity/HeatedCavity_ConvStudyPostProc.ipynb")]
        [Test]
        static public void Run__HeatedCavityConvergenceStudy() {
            Console.WriteLine("HeatedCavity_ConvStudy @ FDYcluster");


            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "HeatedCavity_ConvergenceStudy",
                "HeatedCavity_ConvergenceStudy*",
                "delete_HeatedCavityConvergenceStudy",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedSquareCavity/HeatedCavity_ConvStudy.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedSquareCavity/HeatedCavity_ConvStudyPostProc.ipynb");

        }


        /// <summary>
        /// Test of the Low-Mach solver;
        /// Publication Results from: Gutierrez and Kummer, 2021, A fully coupled high-order Discontinuous Galerkin method for diffusion flames in a low-Mach number framework
        /// </summary>
        [NUnitFileToCopyHack("LowMach/HeatedSquareCavity/HeatedCavity_NusseltStudy.ipynb", "LowMach/HeatedSquareCavity/HeatedCavity_NusseltStudyPostProc.ipynb")]
        [Test]
        static public void Run__HeatedCavityNusseltStudy() {
            //ValidationTestRunner.WorksheetTests_Local.Run__HeatedCavityNusseltStudy
            Console.WriteLine("HeatedCavityNusseltStudy @ FDYcluster");
            //System.Diagnostics.Debugger.Launch();
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "HeatedCavity_NusseltStudy",
                "HeatedCavity_NusseltStudy*",
                "delete_HeatedCavityNusseltStudy",
                new TimeSpan(days: 30, hours: 0, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedSquareCavity/HeatedCavity_NusseltStudy.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedSquareCavity/HeatedCavity_NusseltStudyPostProc.ipynb");

        }

        /// <summary>
        /// Steady state Counter diffusion flame calculation using the Low-Mach solver
        /// Publication Results from: Gutierrez and Kummer, 2021, A fully coupled high-order Discontinuous Galerkin method for diffusion flames in a low-Mach number framework
        /// </summary>
        [NUnitFileToCopyHack(
            "LowMach/DiffusionFlames/CounterDiffusionFlame/CounterFlowFlame_Calculations.ipynb",
            "LowMach/DiffusionFlames/CounterDiffusionFlame/CounterFlowFlame_PostProc.ipynb",
            "LowMach/DiffusionFlames/CounterDiffusionFlame/ML*.txt")]
        [Test]
        static public void Run__CounterDiffusionFlame() {
            Console.WriteLine("CounterFlowFlame @ FDYcluster");
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "CounterFlowFlame_MF_FullComparison",
                "CounterFlowFlame_MF_FullComparison*",
                "delete_CounterDiffusionFlame",
                new TimeSpan(days: 30, hours: 0, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("LowMach/DiffusionFlames/CounterDiffusionFlame/CounterFlowFlame_Calculations.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LowMach/DiffusionFlames/CounterDiffusionFlame/CounterFlowFlame_PostProc.ipynb");
        }





        [NUnitFileToCopyHack("LowMach/DiffusionFlames/CoFlowDiffusionFlame/CoFlowFlame_Calculations.ipynb")]
        [Test]
        static public void Run__CoFlowDiffusionFlame() {
            Console.WriteLine("CoFlowDiffusionFlame @ FDYcluster");

            string really = System.Environment.GetEnvironmentVariable("RUN_COFLOWFLAME");
            if (really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping Run__CoFlowDiffusionFlame ");
                return;
            } else {
                Console.WriteLine("RUN_COFLOWFLAME = " + really);
            }
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "CounterFlowFlame",
                "CounterFlowFlame*",
                "delete_CoFlowDiffusionFlame",
                new TimeSpan(days: 30, hours: 0, minutes: 0, seconds: 1));


            ValidationTestRunnerMain.RunWorksheet("LowMach/DiffusionFlames/CoFlowDiffusionFlame/CoFlowFlame_Calculations.ipynb");

        }


        [NUnitFileToCopyHack("LowMach/HeatedBackwardFacingStep/HeatedBackwardFacingStep.ipynb", "LowMach/HeatedBackwardFacingStep/HeatedBackwardFacingStep_PostProc.ipynb", "LowMach/HeatedBackwardFacingStep/*.txt")]

        [Test]
        static public void Run__HeatedBackwardFacingStep() {


            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "HeatedBackwardFacingStep",
                "HeatedBackwardFacingStep*",
                "delete_HeatedBackwardFacingStep",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedBackwardFacingStep/HeatedBackwardFacingStep.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LowMach/HeatedBackwardFacingStep/HeatedBackwardFacingStep_PostProc.ipynb");

        }





        /// <summary>
        /// Convergence study of a pseudo 1-D configuration for a diffusion flame calculation using the Low-Mach solver
        /// Publication Results from: Gutierrez and Kummer, 2021, A fully coupled high-order Discontinuous Galerkin method for diffusion flames in a low-Mach number framework
        /// </summary>
        [NUnitFileToCopyHack(
            "LowMach/DiffusionFlames/ChamberedDiffusionFlame/ChamberFlame_ConvStudy_Calculations.ipynb",
            "LowMach/DiffusionFlames/ChamberedDiffusionFlame/ChamberFlame_ConvStudy_PostProc.ipynb" //
            )]
        [Test]
        static public void Run__DiffusionFlameConvergenceStudy() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__DiffusionFlameConvergenceStudy
            Console.WriteLine("Convergence study of diffusion flame @ FDYcluster");

            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "DiffFlameConvergenceStudy",
                "DiffFlameConvergenceStudy*",
                "delete_DiffusionFlameConvergenceStudy",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("LowMach/DiffusionFlames/ChamberedDiffusionFlame/ChamberFlame_ConvStudy_Calculations.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LowMach/DiffusionFlames/ChamberedDiffusionFlame/ChamberFlame_ConvStudy_PostProc.ipynb");
        }


        /*
        /// <summary> 
        /// Dummy to test the basic functionality of the validation test runner
        /// </summary>
        [NUnitFileToCopyHack("ue2Basics/ue2Basics.ipynb")]
        [Test]
        static public void Run__ue2Basics() {
            ValidationTestRunnerMain.RunWorksheet("ue2Basics/ue2Basics.ipynb");

            try {
                File.Move("ue2Basics.ipynb", "ue2Basics-LocalValidation.html");
            } catch(Exception e) {
                Console.Error.WriteLine($"File copy exception: {e.GetType()} : {e.Message}");
            }
        }
        */


        /// <summary> 
        /// 3D oscillating droplet, using the two-phase solver;
        /// simulations for the DACH-Cooperation with TU Graz (Prof. Brenn)
        /// </summary>
        [NUnitFileToCopyHack(
            "Oscillating-Droplet/Droplet3D.ipynb",
            "Oscillating-Droplet/data/InitialValues/m*/surfaceDrop*.txt",
            "Oscillating-Droplet/data/InitialValues/m*/radialVel*.txt",
            "Oscillating-Droplet/data/InitialValues/m*/polarVel*.txt")]
        [Test]
        static public void Run__Droplet3D() {

            string really = null; // System.Environment.GetEnvironmentVariable("RUN_DROPLET");
            if (really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping Run__Droplet3D ");
                return;
            } else {
                Console.WriteLine("RUN_DROPLET = " + really);
            }

            Console.WriteLine("Lets go...");

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "OscillatingDroplet3D",
                "OscillatingDroplet3D*",
                "delete_Droplet3D",
                new TimeSpan(days: 60, hours: 1, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("Oscillating-Droplet/Droplet3D.ipynb");

            ValidationTestRunnerMain.RunWorksheet("Oscillating-Droplet/Droplet3D.ipynb");
        }


        /// <summary> 
        /// combustion of a droplet,
        /// from the project of Juan Gutierrez
        /// </summary>
        [NUnitFileToCopyHack("examples/CombustingDroplet/CombustingDroplet.ipynb")]
        [Test]
        static public void Run__CombustingDroplet() {
            //ValidationTestRunner.WorksheetTests_Local.Run__CombustingDroplet
            string really = System.Environment.GetEnvironmentVariable("RUN_COMBDROPLET");
            if (really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping Run__CombustingDroplet ");
                return;
            } else {
                Console.WriteLine("RUN_COMBDROPLET = " + really);
            }

            Console.WriteLine("Lets go...");

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "CombustingDroplet",
                "CombustingDroplet*",
                "delete_CombustingDroplet",
                new TimeSpan(days: 60, hours: 1, minutes: 0, seconds: 1));
            ValidationTestRunnerMain.RunWorksheet("CombustingDroplet.ipynb");
        }

        /// <summary> 
        /// 3D oscillating droplet, using the two-phase solver;
        /// simulations for the DACH-Cooperation with TU Graz (Prof. Brenn)
        /// </summary>
        [NUnitFileToCopyHack("Oscillating-Droplet/Droplet3D-FirstPeriodStudy.ipynb", "Oscillating-Droplet/data/InitialValues/m*/surfaceDrop*.txt")]
        [Test]
        static public void Run__Droplet3D_FirstPeriodStudy() {

            string really = null; // System.Environment.GetEnvironmentVariable("RUN_DROPLET_FIRSTPERIOD");
            if (really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping Run__Droplet3D_FirstPeriodStudy ");
                return;
            } else {
                Console.WriteLine("RUN_DROPLET_FIRSTPERIOD = " + really);
            }

            Console.WriteLine("Lets go...");

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "OscillatingDroplet3D_FirstPeriodStudy",
                "OscillatingDroplet3D_FirstPeriodStudy*",
                "delete_Droplet3D_FirstPeriodStudy",
                new TimeSpan(days: 30, hours: 1, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("Oscillating-Droplet/Droplet3D-FirstPeriodStudy.ipynb");
        }


        /*
        /// <summary>
        /// Linear solver performance:
        /// - constant coefficient Poisson problem (only DG, no XDG)
        /// - one MPI core
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-NodeSolverPerformance/PoissonConstCoeff/LinslvPerf_ConstPoissonMpi1.ipynb", "handbook/apdx-NodeSolverPerformance/PoissonConstCoeff/LinslvPerf_ConstPoissonMpi1-Pt2.ipynb")]
        [Test]
        static public void Run__PoissonPerformance() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__PoissonPerformance
            

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "LinslvPerf_ConstPoissonMpi1", 
                "LinslvPerf_ConstPoissonMpi1*", new TimeSpan(days: 10, hours: 1, minutes: 0, seconds: 1)); 

            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_ConstPoissonMpi1.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_ConstPoissonMpi1-Pt2.ipynb");
        }


        /// <summary>
        /// Linear solver performance:
        /// - XDG Poisson problem (1:1000 diffusion factor)
        /// - one MPI core
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-NodeSolverPerformance/XDGPoisson/LinslvPerf_XdgPoissonSer.ipynb", "handbook/apdx-NodeSolverPerformance/XDGPoisson/LinslvPerf_XdgPoissonSer-Pt2.ipynb")]
       
        [Test]
        static public void Run__XdgPoissonPerformance() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__XdgPoissonPerformance
            

            //string really = System.Environment.GetEnvironmentVariable("RUN_XDGPOISSONPERFORMANCE");
            //if(really.IsEmptyOrWhite()) {
            //    Console.WriteLine("skipping Run__XdgPoissonPerformance ");
            //    return;
            //} else {
            //    Console.WriteLine("RUN_XDGPOISSONPERFORMANCE = " + really);
            //}

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "LinslvPerf_XdgPoissonSer", 
                "LinslvPerf_XdgPoissonSer*", new TimeSpan(days: 10, hours: 1, minutes: 0, seconds: 1)); 

            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_XdgPoissonSer.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_XdgPoissonSer-Pt2.ipynb");

        }

        /*
        /// <summary>
        /// Linear solver performance:
        /// - Steady-State Stokes problem 
        /// - one MPI core
        /// </summary>
        /// <remarks>
        /// Benchmark is originally proposed in 
        /// 'p‑Multilevel Preconditioners for HHO Discretizations of the Stokes Equations with Static Condensation',
        /// by L. Botti and D. Di Pietro (https://doi.org/10.1007/s42967-021-00142-5)
        /// </remarks>
        [NUnitFileToCopyHack("handbook/apdx-NodeSolverPerformance/BottiPietroStokes/LinslvPerf_BottiPietroStokes2D.ipynb")]

        [Test]
        static public void Run__BottiPietroStokes2DPerformance() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__BottiPietroStokes2DPerformance

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "LinslvPerf_BottiPietroStokes2D",
                "LinslvPerf_BottiPietroStokes2D*", new TimeSpan(days: 10, hours: 1, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_BottiPietroStokes2D.ipynb");
        }

        /*
        /// <summary>
        /// Linear solver performance:
        /// - Steady-State XDG Stokes problem (water droplet in air)
        /// - one MPI core
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-NodeSolverPerformance/XDGStokes/LinslvPerf_XdgStokes.ipynb")]

        [Test]
        static public void Run__XdgStokesPerformance() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__XdgStokesPerformance

            //string really = System.Environment.GetEnvironmentVariable("RUN_XDGSTOKESPERFORMANCE");
            //if(really.IsEmptyOrWhite()) {
            //    Console.WriteLine("skipping Run__XdgStokesPerformance ");
            //    return;
            //} else {
            //    Console.WriteLine("RUN_XDGSTOKESPERFORMANCE = " + really);
            //}

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                "LinslvPerf_XdgStokes",
                "LinslvPerf_XdgStokes*", new TimeSpan(days: 10, hours: 1, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_XdgStokes.ipynb");
        }
        */


        /// <summary>
        /// Linear solver performance:
        /// - Steady-State XDG Stokes problem (water droplet in air)
        /// - one MPI core
        /// </summary>
        [NUnitFileToCopyHack(
            "handbook/apdx-NodeSolverPerformance/unified/LinslvPerf_ConstPoissonMpi1.ipynb",
            "handbook/apdx-NodeSolverPerformance/unified/LinslvPerf_XdgPoissonSer.ipynb",
            "handbook/apdx-NodeSolverPerformance/unified/LinslvPerf_BottiPietroStokes2D.ipynb",
            "handbook/apdx-NodeSolverPerformance/unified/LinslvPerf_BottiPietroStokes3D.ipynb",
            "handbook/apdx-NodeSolverPerformance/unified/LinslvPerf_XdgStokes.ipynb",
            "handbook/apdx-NodeSolverPerformance/unified/LinslvPerf_Evaluation.ipynb")]

        [Test]
        static public void Run__LinslvPerfSer() {
            // --test=ValidationTestRunner.WorksheetTests_Local.Run__LinslvPerfSer

            string PROJECT_NAME = System.Environment.GetEnvironmentVariable("LinslvPerfSer") ?? "LinslvPerfSer"; // this allows to modify the project name for testing purposes


            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                PROJECT_NAME,
                $"{PROJECT_NAME}*",
                "delete_LinslvPerfSer",
                new TimeSpan(days: 60, hours: 0, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_ConstPoissonMpi1.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_XdgPoissonSer.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_BottiPietroStokes2D.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_BottiPietroStokes3D.ipynb");
            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_XdgStokes.ipynb");

            ValidationTestRunnerMain.RunWorksheet("LinslvPerf_Evaluation.ipynb");
        }


        // <summary>
        // Serial Shock Fitting Solver:
        // - SAIDT - space-time Scalar Advection in 1D
        /// - BUIDT - space-time Burgers Equation in 1D
        // - XESF  - Inviscid Euler Equation in 2D
        // </summary>
        //[NUnitFileToCopyHack(
        //    "examples/ShockFitting/SAIDT/SAIDT_Validation.ipynb"
        //    //"../internal/src/private-seb/Notebooks/BUIDT/BUIDT_Validation.ipynb",
        //    //"../internal/src/private-seb/Notebooks/XESF/XESF_Validation.ipynb",
        //    //"../internal/src/private-seb/Notebooks/ValidationEvaluation.ipynb"
        //    )]
        //[Test] DEACTIVATED: TEST TAKES MORE THAN 8 HOURS
        //static public void Run__ShockFittingSer() {
        //    // --test=ValidationTestRunner.WorksheetTests_Local.Run__LinslvPerfSer

        //    string PROJECT_NAME = System.Environment.GetEnvironmentVariable("ShockFitting") ?? "ShockFitting"; // this allows to modify the project name for testing purposes


        //    // delete the database if it is more than XX days old;
        //    // this will cause a re-execution of all computations
        //    // otherwise, i.e. if the database is not deleted, sessions from the database 
        //    ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
        //        PROJECT_NAME,
        //        $"{PROJECT_NAME}*",
        //        "delete_ShockFitting",
        //        new TimeSpan(days: 40, hours: 0, minutes: 0, seconds: 1));

        //    ValidationTestRunnerMain.RunWorksheet("SAIDT_Validation.ipynb");
        //    //ValidationTestRunnerMain.RunWorksheet("BUIDT_Validation.ipynb");
        //    //ValidationTestRunnerMain.RunWorksheet("XESF_Validation.ipynb");
        //    //ValidationTestRunnerMain.RunWorksheet("ValidationEvaluation.ipynb");
        //}

        /// <summary> 
        /// Testing of memory scaling.
        /// </summary>
        [NUnitFileToCopyHack("memprofile/memprofile.ipynb")]
        [Test]
        static public void Run__memprofile() {

            const string PROJECT_NAME = "memprofile";

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                PROJECT_NAME,
                $"{PROJECT_NAME}*",
                "delete_memprofile",
                new TimeSpan(days: 10, hours: 0, minutes: 0, seconds: 1));

            NotebookRunner.DeleteDeployments("memprofile*");

            ValidationTestRunnerMain.RunWorksheet("memprofile/memprofile.ipynb");

        }
    }



    /// <summary>
    /// NUnit entry point for each example worksheet which represents a short-running validation test;
    /// </summary>
    /// <remarks>
    /// - short running rests are fully re-computed every timem
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

                NotebookRunner.DeleteDatabase("Demo_BoundaryAndInitialData");
                NotebookRunner.DeleteDeployments("Demo_BoundaryAndInitialData*");
                ValidationTestRunnerMain.RunWorksheet("BoundaryAndInitialData/BoundaryAndInitialData.ipynb");
            } finally {
                JupyterMutex.ReleaseMutex();
            }
        }

        /// <summary> Testing of respective worksheet. </summary>
        [NUnitFileToCopyHack("MetaJobManager/MetaJobManager.ipynb")]
        [Test]
        static public void Run__MetaJobManager() {
            //--test=ValidationTestRunner.WorksheetTests_Local.Run__MetaJobManager
            Mutex JupyterMutex = new Mutex(false, "MetaJobManager_Tutorial");
            try {
                JupyterMutex.WaitOne();
                NotebookRunner.DeleteDatabase("MetaJobManager_Tutorial");
                NotebookRunner.DeleteDeployments("MetaJobManager_Tutorial*");
                ValidationTestRunnerMain.RunWorksheet("MetaJobManager/MetaJobManager.ipynb");
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
                NotebookRunner.DeleteDatabase("ConvStudyTutorial");
                NotebookRunner.DeleteDeployments("ConvStudyTutorial*");
                ValidationTestRunnerMain.RunWorksheet("convergenceStudyTutorial/convStudy.ipynb");
            } finally {
                JupyterMutex.ReleaseMutex();
            }
        }
    }

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





        /*

        #region prallel performance
        /// <summary>
        /// - strong scaling for XDG Poisson problem (1:1000 diffusion factor)
        /// - core sweep: 4,8,16,32,64,128,256
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-MPISolverPerformance/XdgPoisson/*.ipynb")]
        [Test]
        static public void Run__Strong_XdgPoisson() {
            if(!IsActive("RUN_PARALLELPERFORMANCE")) return;
            string dbname = "XdgPoisson_Strong";
            System.Environment.SetEnvironmentVariable("DATABASE_NAME", dbname);
            DelteDBWhenOld(dbname);
            RunSheets("strong","XP");
        }

        /*
        /// <summary>
        /// - weak scaling for XDG Poisson problem (1:1000 diffusion factor)
        /// - core sweep: 4,8,16,32,64,128,256
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-MPISolverPerformance/XdgPoisson/*.ipynb")]
        [Test]
        static public void Run__Weak_XdgPoisson() {
            if(!IsActive("RUN_PARALLELPERFORMANCE")) return;
            string dbname = "XdgPoisson_Weak";
            System.Environment.SetEnvironmentVariable("DATABASE_NAME", dbname);
            DelteDBWhenOld(dbname);
            RunSheets("weak","XP");
        }

        /// <summary>
        /// - strong scaling for (X)DG Stokes problem (rotating Sphere with inflow)
        /// - core sweep: 4,8,16,32,64,128,256
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-MPISolverPerformance/IBMrotSphere/*.ipynb")]
        [Test]
        static public void Run__Strong_IBMrotSphere() {
            if(!IsActive("RUN_PARALLELPERFORMANCE")) return;
            string dbname = "DGrotSphere_Strong";
            System.Environment.SetEnvironmentVariable("DATABASE_NAME", dbname);
            DelteDBWhenOld(dbname);
            RunSheets("strong","IrS");
        }

        /// <summary>
        /// - weak scaling for (X)DG Stokes problem (rotating Sphere with inflow)
        /// - core sweep: 4,8,16,32,64,128,256
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-MPISolverPerformance/IBMrotSphere/*.ipynb")]
        [Test]
        static public void Run__Weak_IBMrotSphere() {
            if(!IsActive("RUN_PARALLELPERFORMANCE")) return;
            string dbname = "DGrotSphere_Weak";
            System.Environment.SetEnvironmentVariable("DATABASE_NAME", dbname);
            DelteDBWhenOld(dbname);
            RunSheets("weak","IrS");
        }

        /// <summary>
        /// - strong scaling for (X)DG NSE problem (rotating Cube with inflow)
        /// - core sweep: 4,8,16,32,64,128,256
        /// </summary>
        [NUnitFileToCopyHack("handbook/apdx-MPISolverPerformance/IBMrotCube/*.ipynb")]
        [Test]
        static public void Run__Strong_IBMrotCube() {
            if(!IsActive("RUN_PARALLELPERFORMANCE")) return;
            string dbname = "DGrotCube_Strong";
            System.Environment.SetEnvironmentVariable("DATABASE_NAME", dbname); 
            DelteDBWhenOld(dbname);
            RunSheets("strong","IrC");
        }
        

        static private string GetJobName() {
            var envout = new EnvironmentVariableTarget();
            System.Environment.GetEnvironmentVariable("JOB_NAME", envout);
            return Convert.ToString(envout);
        }

        static private void DelteDBWhenOld(string dbname) {
            string fullname = GetJobName() + "_" + dbname;
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                fullname,
                fullname + "*", 
                "delete_" + dbname,
                new TimeSpan(days: 0, hours: 0, minutes: 0, seconds: 1));
        }

        static private void RunSheets(string scalingtype, string suffix) {
            var sheets2execute = new List<string>();
            if(!IsActive("PROCESSING_ONLY")) {
                sheets2execute.Add($"Part1_0-Calculations_{scalingtype}_{suffix}");
            }
            // add more files to this list, missing ones will be skipped
            sheets2execute.AddRange(new string[]{
                $"Part1_1-Export_JSON_{suffix}",
                $"Part2_0-scaling_{suffix}",
                $"Part3_0-profiling_{suffix}"
            });

            foreach(string sname in sheets2execute) {
                string filename = sname + ".ipynb";
                if(!System.IO.File.Exists(filename)) {
                    throw new FileNotFoundException($"{filename} is non existent!");
                }
                try {
                    ValidationTestRunnerMain.RunWorksheet(filename);
                } catch (Exception ex) {
                    Console.WriteLine(ex.Message);
                } finally {
                    System.IO.File.Move(sname + ".html", sname + $"_{scalingtype}" + ".html", true);
                }
            }
        }

        static private bool IsActive(string EnvVar) {
            string really = System.Environment.GetEnvironmentVariable(EnvVar);
            if(really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping test:"+ EnvVar + "is empty or white");
                return false;
            } else {
                Console.WriteLine(EnvVar + " = " + really);
                return true;
            }
        }
        #endregion
        */



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

            // delete the database if it is more than XX days old;
            // this will cause a re-execution of all computations
            // otherwise, i.e. if the database is not deleted, sessions from the database 
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
                PROJECT_NAME,
                $"{PROJECT_NAME}*",
                "delete_LinslvPerfPar",
                new TimeSpan(days: 40, hours: 0, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_ConstPoisson.ipynb", allowErrors: true);
            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_XdgPoisson.ipynb", allowErrors: true);
            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_BottiPietroStokes3D.ipynb", allowErrors: true);
            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_XdgStokes.ipynb", allowErrors: true);

            ValidationTestRunnerMain.RunWorksheet("ParLinslvPerf_Evaluation.ipynb", allowErrors: false);
        }

        /*
        /// <summary>
                /// </summary>
        [OneTimeSetUp]
        static public void StartMinibatch() {
            MiniBatchProcessor.Server.StartIfNotRunning(RunExternal: true);
        }

        /// <summary>
        /// counterpart to <see cref="StartMinibatch"/>
        /// </summary>
        [OneTimeTearDown]
        static public void StopMiniBatch() {
            
        }
        */
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
                "MetaJobManager_Tutorial*",
                "delete_MetaJobManager",
                new TimeSpan(days: 5, hours: 1, minutes: 0, seconds: 1));


            Console.WriteLine("Testing MetaJobManager @ Lichtenberg...");
            ValidationTestRunnerMain.RunWorksheet("MetaJobManager/MetaJobManager.ipynb");
            Console.WriteLine("Finished MetaJobManager @ Lichtenberg.");
        }


    }


    static class ValidationTestRunnerMain {

        /// <summary>
        /// Deletes a database <paramref name="Directory"/> if it older than specified by <paramref name="DeletionAge"/>
        /// - this will cause a re-execution of all computations specified in the test.
        /// - otherwise, i.e. if the database is not deleted, existing sessions from the database may be taken if they match 
        ///    the respective compute job (<see cref="Job"/>)
        /// 
        /// Note: the database must be located beneath the <see cref="BatchProcessorClient.AllowedDatabasesPaths"/>
        /// of the <see cref="BoSSSshell.GetDefaultQueue"/>.
        /// </summary>
        public static void DeleteDatabaseAndDeploymentsWhenOld(string Directory, string DeployMents, string EnforceDeletionEnvVar, TimeSpan DeletionAge) {
            bool runfromBackup = File.Exists("BOSSS_RUNTESTFROMBACKUP.txt");
            if(runfromBackup) {
                Console.WriteLine("Run-From-Backup mode activated ('BOSSS_RUNTESTFROMBACKUP.txt' file exists); not messing with any database or deployment directories.");
                return;
            }
            bool masterEnforceDeletion = File.Exists("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt");
            if(runfromBackup) {
                Console.WriteLine("Run-From-Backup mode activated ('BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt' file exists); enforcing deletion.");
            }

            var nau = DateTime.Now;

            EnforceDeletionEnvVar = EnforceDeletionEnvVar.ToUpperInvariant();
            

            bool enforce = false;
            if(EnforceDeletionEnvVar != null) {
                string really = System.Environment.GetEnvironmentVariable(EnforceDeletionEnvVar) ?? "";
                string really2 = System.Environment.GetEnvironmentVariable("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER") ?? "";
                if ((really.IsNonEmpty() && really.Trim() != "0")
                    || (really2.IsNonEmpty() && really2.Trim() != "0")
                    || masterEnforceDeletion) {

                    Console.WriteLine($"Enforcing deletion of old runs, since environment variables {EnforceDeletionEnvVar}/BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER are set to {really}/{really2} (not 0) or file 'BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt' exists!");
                    enforce = true;

                }
            }

            //Console.WriteLine("Note Deletion deactivated !!!!!!!!!!!!!!!!!!!!!!!!!");
            //return; /*

            List<string> DirsToIgnore = new List<string>();
            bool delete = false; // deletion is decided based on the age of the database folder
            foreach(var q in BoSSSshell.ExecutionQueues) {

                Console.WriteLine("Checking queue: " + q);

                // move old databases to backup
                // ============================

                {
                    //var dirsToDelete = new HashSet<DirectoryInfo>();
                    //bool delete = false;
                    foreach(var allowedPath in q.AllowedDatabasesPaths) {
                        var localBaseDir = new DirectoryInfo(allowedPath.LocalMountPath);

                        var dbDirs = localBaseDir.GetDirectories(Directory, SearchOption.TopDirectoryOnly);
                        foreach(var db in dbDirs) {
                            var age = nau - db.CreationTime;
                            if((age > DeletionAge || enforce) && db.Exists) {
                                delete = true;
                                string NewName = Path.Combine(Path.GetDirectoryName(db.FullName), "bkup-" + db.CreationTime.ToString("yyyyMMMdd_HHmmss") + "." + db.Name);
                                Console.WriteLine("Renaming existing database: " + db.FullName + " -> " + NewName);
                                db.MoveTo(NewName);
                                DirsToIgnore.Add(NewName);
                            }
                        }
                    }
                }
            }

            if(delete) {
                Console.WriteLine("deletion of old deployments in execution queues");
                foreach(var q in BoSSSshell.ExecutionQueues) {
                    // delete old deployment
                    // =====================

                    {
                        Console.WriteLine("searching for deployments in: " + q.DeploymentBaseDirectory);

                        var deplDirs = (new DirectoryInfo(q.DeploymentBaseDirectory)).GetDirectories(DeployMents, SearchOption.TopDirectoryOnly);
                        foreach(var d in deplDirs) {
                            if(DirsToIgnore.Any(keepDir => keepDir.Equals(d.FullName)))
                                continue; // don't delete the db

                            if(delete && d.Exists) {
                                Console.WriteLine("Deleting deployment dir: " + d.FullName);
                                d.Delete(true);
                            }
                        }
                    }


                    Console.Out.Flush();
                }
            } else {
                Console.WriteLine("No deletion of old deployments");
            }
            //*/
        }

        /// <summary>
        /// Runs some worksheet contained in the BoSSS handbook.
        /// </summary>
        static public void RunWorksheet(string NotebookPartialPath, bool allowErrors = false) {
            using(new BoSSS.Application.TutorialTests.NotebookRunner(NotebookPartialPath, "", allowErrors)) { }
        }

        static int Main(string[] args) {
            PublicTestRunner.PublicTestRunnerMain.TimeOutSec = 24 * 3600 * 10; // 10 days
            try {
                return PublicTestRunner.PublicTestRunnerMain._Main(args, new ValidationTests());
            } catch(Exception e) {
                // note: this seemingly useless try-catch is here since our test runner server (FDYGITRUNNER)
                // seems to silently fail on all exceptions thrown after MPI init.
                Console.WriteLine("Got some exception: " + e);
                return -667;
            }
        }
    }
}
