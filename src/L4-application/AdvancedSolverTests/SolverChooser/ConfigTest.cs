using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AdvancedSolverTests.SolverChooser
{
    [TestFixture]
    public static class ConfigTest
    {

        [Test]
        public static void TestLinearSolverConfigurations() {
            //Arrange --- configs
            var ACS = new AppControlSolver();
            LinearSolverConfig lconfig = ACS.LinearSolver;
            NonLinearSolverConfig nlconfig = ACS.NonLinearSolver; // is not of interest in this test, but we have to set this ...
            lconfig.verbose = true;
            lconfig.NoOfMultigridLevels = 3;
            lconfig.TargetBlockSize = 10;
            var SF = new SolverFactory(nlconfig, lconfig);
            
            //Arrange --- Multigrid stuff
            AggregationGridData[] seq;
            var MGO = Utils.CreateTestMGOperator(out seq, Resolution: 10);
            var changeofbasisis = Utils.GetAllMGConfig(MGO);
            var agggridbasisis = Utils.GetAllAggGridBasis(MGO);

            //Arrange --- get available lincodes
            var lincodes = (LinearSolverCode[])Enum.GetValues(typeof(LinearSolverCode));
            ISolverSmootherTemplate LinSolver = null;
            TestDelegate lindlg = () => SF.GenerateLinear(out LinSolver, agggridbasisis, changeofbasisis);

            //Act and Assert
            foreach (LinearSolverCode code in lincodes) {
                SF.Clear();
                lconfig.SolverCode = code;
                if(code==LinearSolverCode.selfmade)
                    SF.Selfmade_linsolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };
                Assert.DoesNotThrow(lindlg, "", null);
                Assert.IsNotNull(LinSolver);
            }
        }

        [Test]
        public static void TestNonLinearSolverConfigurations() {
            //Arrange --- set configs
            var ACS = new AppControlSolver();
            NonLinearSolverConfig nlconfig = ACS.NonLinearSolver;
            LinearSolverConfig lconfig = ACS.LinearSolver;
            lconfig.verbose = true;
            lconfig.NoOfMultigridLevels = 3;
            lconfig.TargetBlockSize = 10;
            nlconfig.verbose = true;
            var SF = new SolverFactory(nlconfig, lconfig);

            //Arrange --- get test multigrid operator stuff
            AggregationGridData[] seq;
            var MGO = Utils.CreateTestMGOperator(out seq, Resolution: 10);
            var map = MGO.Mapping;
            var changeofbasisis = Utils.GetAllMGConfig(MGO);
            var agggridbasisis = Utils.GetAllAggGridBasis(MGO);

            //Arrange --- get nonlinear codes available
            var nonlincodes = (NonLinearSolverCode[])Enum.GetValues(typeof(NonLinearSolverCode));
            NonlinearSolver NLsolver = null;

            //Arrange --- get test linear Solver to set in NLsolver
            ISolverSmootherTemplate LinSolver = null;
            LinearSolverCode[] LinTestcandidates = { LinearSolverCode.classic_pardiso, LinearSolverCode.exp_gmres_levelpmg }; // in order to test the GMRES variants of the NL solver

            //Act and Assert
            foreach (var lincode in LinTestcandidates) {
                lconfig.SolverCode = lincode;
                TestDelegate nldlg = () => SF.GenerateNonLin(out NLsolver, out LinSolver, null, agggridbasisis, changeofbasisis, seq);
                SF.Clear();
                foreach (NonLinearSolverCode nlcode in nonlincodes) {
                    nlconfig.SolverCode = nlcode;
                    if (nlconfig.SolverCode == NonLinearSolverCode.selfmade) {
                        SF.Selfmade_nonlinsolver = new Newton(null, agggridbasisis, changeofbasisis);
                    }
                    Assert.DoesNotThrow(nldlg, "", null);
                    Assert.IsNotNull(NLsolver);
                }
                Console.WriteLine("====");
            }
        }

    }
}
