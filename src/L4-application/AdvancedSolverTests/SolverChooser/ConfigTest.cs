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
            // --test=AdvancedSolverTests.SolverChooser.ConfigTest.TestLinearSolverConfigurations
            ////Arrange --- configs
            //var ACS = new AppControlSolver();
            //var lconfig = ACS.LinearSolver;
            //NonLinearSolverConfig nlconfig = ACS.NonLinearSolver; // is not of interest in this test, but we have to set this ...
            //lconfig.verbose = true;
            //lconfig.NoOfMultigridLevels = 3;
            //lconfig.TargetBlockSize = 10;


            //Arrange --- Multigrid stuff
            using(var O = Utils.CreateTestMGOperator(Resolution: 10, DGOrder: 3)) {
                AggregationGridData[] seq = O.MGSeq;
                var MGO = O.MGOp;
                var changeofbasisis = Utils.GetAllMGConfig(MGO);
                var agggridbasisis = Utils.GetAllAggGridBasis(MGO);

                //Arrange --- get available lincodes
                var lincodes = (LinearSolverCode[])Enum.GetValues(typeof(LinearSolverCode));

                //Act and Assert
                foreach(LinearSolverCode code in lincodes) {
                    

                    Assert.DoesNotThrow(() => code.GetConfig().CreateInstance(O.MGOp), "", null);
                    Assert.IsNotNull(code.GetConfig().CreateInstance(O.MGOp));
                }
            }
        }

        [Test]
        public static void TestNonLinearSolverConfigurations() {

            //Arrange --- get test multigrid operator stuff
            using(var O = Utils.CreateTestMGOperator(Resolution: 10, DGOrder: 3)) {
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
                foreach(var lincode in LinTestcandidates) {
                    Assert.DoesNotThrow(() => lincode.GetConfig().CreateInstance(MGO), "", null);
                    
                    //lconfig.SolverCode = lincode;
                    //TestDelegate nldlg = () => SF.GenerateNonLin(out NLsolver, out LinSolver, null, agggridbasisis, changeofbasisis, seq);
                    //SF.Clear();
                    
                    foreach(NonLinearSolverCode nlcode in nonlincodes) {
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
