using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using ilPSP.LinSolvers.PETSc;
using NUnit.Framework;
using System;

namespace AdvancedSolverTests.SolverChooser {
    [TestFixture]
    public static class ConfigTest {

        [Test]
        [Ignore("Fails due to timeout", Until = "2025-04-10 00:00:00Z")]
        public static void TestLinearSolverConfigurations() {
            // --test=AdvancedSolverTests.SolverChooser.ConfigTest.TestLinearSolverConfigurations
            ////Arrange --- configs
            //var ACS = new AppControlSolver();
            //var lconfig = ACS.LinearSolver;
            //NonLinearSolverConfig nlconfig = ACS.NonLinearSolver; // is not of interest in this test, but we have to set this ...
            //lconfig.verbose = true;
            //lconfig.NoOfMultigridLevels = 3;
            //lconfig.TargetBlockSize = 10;


            Console.WriteLine("[ConfigTest] Building XDG test matrix and multigrid hierarchy...");

            //Arrange --- Multigrid stuff
            using ( var O = Utils.CreateTestMGOperator(Resolution: 10, DGOrder: 3) ) {
                Console.WriteLine("[ConfigTest] Matrix and multigrid hierarchy are ready.");
                AggregationGridData[] seq = O.MGSeq;
                var MGO = O.MGOp;
                var changeofbasisis = Utils.GetAllMGConfig(MGO);
                var agggridbasisis = Utils.GetAllAggGridBasis(MGO);

                //Arrange --- get available lincodes
                var lincodes = (LinearSolverCode[])Enum.GetValues(typeof(LinearSolverCode));

                //Act and Assert

                for(int i = 0; i < lincodes.Length; i++) {
                    LinearSolverCode code = lincodes[i];
                    CollectiveCheckpoint($"before {code}", 900000 + 2 * i);
                    var config = code.GetConfig();
                    ISolverSmootherTemplate solver = null;
                    try {
                        if(code == LinearSolverCode.SchurPrecond || code == LinearSolverCode.SchurComplementForStokes) {
                            Console.WriteLine($"Skipping {LinearSolverCode.SchurPrecond}: the Schur Complement requires a (D+1)x(D+1) - system, an the test system here is only 2x2.");
                            continue;
                        }

                        if(code == LinearSolverCode.direct_petsc_cholesky) {
                            Console.WriteLine("Skipping direct_petsc_cholesky: the generic XDG test matrix is not guaranteed SPD.");
                            continue;
                        }

                        Console.WriteLine($"[ConfigTest] [{i + 1}/{lincodes.Length}] Initializing {code}...");
                        solver = config.CreateInstance(O.MGOp);
                        Assert.IsNotNull(solver);
                        Console.WriteLine($"[ConfigTest] [{i + 1}/{lincodes.Length}] {code} initialized.");
                    } catch(Exception error) when(IsOkWhenFails(config, error)) {
                        Console.WriteLine($"Skipping {code}: optional PETSc direct backend is unavailable ({error.Message})");
                    } finally {
                        (solver as IDisposable)?.Dispose();
                        CollectiveCheckpoint($"after {code}", 900001 + 2 * i);
                    }
                }
            }
        }

        /// <summary>
        /// Makes transitions between solver configurations collective.
        /// </summary>
        static void CollectiveCheckpoint(string label, int token) {
            //Console.WriteLine($"[ConfigTest] MPI checkpoint: {label}.");
            ilPSP.MPICollectiveWatchDog.WatchAtRelease(
                MPI.Wrappers.csMPI.Raw._COMM.WORLD,
                WaitTimeSeconds: 15.0,
                quitAppIfExpired: true,
                token: token);
        }

        static bool IsOkWhenFails(ISolverFactory config, Exception error) {
            if(!(config is PETScDirectSolverConfig petscConfig))
                return false;
            if(error is NotSupportedException)
                return true;

            return error is ApplicationException
                // PETSC_ERR_SUP (92) is returned when an optional PETSc
                // factorization package was not compiled into the selected
                // native PETSc variant.
                && error.Message.Contains("PETSc error code 92");
        }

        [Test]
        public static void TestNonLinearSolverConfigurations() {

            //Arrange --- get test multigrid operator stuff
            using ( var O = Utils.CreateTestMGOperator(Resolution: 10, DGOrder: 3) ) {
                var MGO = O.MGOp;
                var map = MGO.Mapping;
                var changeofbasisis = Utils.GetAllMGConfig(MGO);
                var agggridbasisis = Utils.GetAllAggGridBasis(MGO);

                //Arrange --- get nonlinear codes available
                var nonlincodes = (NonLinearSolverCode[])Enum.GetValues(typeof(NonLinearSolverCode));
                NonlinearSolver NLsolver = null;

                //Arrange --- get test linear Solver to set in NLsolver
                LinearSolverCode[] LinTestcandidates = { LinearSolverCode.direct_pardiso, LinearSolverCode.exp_gmres_levelpmg, LinearSolverCode.exp_Kcycle_schwarz, LinearSolverCode.exp_Kcycle_schwarz_CoarseMesh, LinearSolverCode.exp_Kcycle_schwarz_PerProcess }; // in order to test the GMRES variants of the NL solver

                //Act and Assert
                foreach ( var lincode in LinTestcandidates ) {
                    Assert.DoesNotThrow(() => lincode.GetConfig().CreateInstance(MGO), "", null);

                    //lconfig.SolverCode = lincode;
                    //TestDelegate nldlg = () => SF.GenerateNonLin(out NLsolver, out LinSolver, null, agggridbasisis, changeofbasisis, seq);
                    //SF.Clear();

                    foreach ( NonLinearSolverCode nlcode in nonlincodes ) {
                        var nlconfig = new NonLinearSolverConfig();
                        nlconfig.SolverCode = nlcode;
                        nlconfig.verbose = true;

                        var SF = new SolverFactory(nlconfig, lincode.GetConfig());
                        SF.GenerateNonLin(out NLsolver, null, agggridbasisis, changeofbasisis);

                        Assert.IsNotNull(NLsolver);
                    }
                }
            }
        }
    }
}
