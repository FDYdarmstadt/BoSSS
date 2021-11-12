using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.Tests
{
    /// <summary>
    /// 
    /// </summary>
    [TestFixture]
    public class TwoPhaseCouplingTests
    {
        
        [Test]
        public static void SolidBallInChannel([Values(2)] int p = 2,
                                    [Values(16)] int res = 16) {

            var C = ZwoLevelSetSolver.ControlFiles.HardCodedControl.AcceleratedBallInChannel(p, res);
            C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 0.0;
            //C.NonLinearSolver.MinSolverIterations = 200;
            C.NonLinearSolver.MaxSolverIterations = 20;

            using(var q = new ZLS()) {
                q.Init(C);
                q.RunSolverMode();
                Assert.IsTrue(q.LastSolverSuccess, "Nonlinear solver did not converge.");
            }
        }

        [Test]
        public static void FluidBallInChannel([Values(2)] int p = 2,
                                    [Values(16)] int res = 16) {
            var C = ZwoLevelSetSolver.ControlFiles.HardCodedControl.AcceleratedFluidBallInChannel(p, res);
            C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 0.0;
            //C.NonLinearSolver.MinSolverIterations = 200;
            C.NonLinearSolver.MaxSolverIterations = 20;

            using(var q = new ZLS()) {
                q.Init(C);
                q.RunSolverMode();
                Assert.IsTrue(q.LastSolverSuccess, "Nonlinear solver did not converge.");
            }
        }

        [Test]
        public static void ExtensionSolidBallInChannel([Values(2)] int p = 2,
                                    [Values(16)] int res = 16) {
            var C = ZwoLevelSetSolver.ControlFiles.HardCodedControl.AcceleratedBallInChannel(p, res);
            C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 0.0;
            C.DisplacementExtension = true;
            //C.NonLinearSolver.MinSolverIterations = 200;
            C.NonLinearSolver.MaxSolverIterations = 20;

            using(var q = new ZLS()) {
                q.Init(C);
                q.RunSolverMode();
                Assert.IsTrue(q.LastSolverSuccess, "Nonlinear solver did not converge.");
            }
        }
    }
}
