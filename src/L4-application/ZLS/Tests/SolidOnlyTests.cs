using BoSSS.Solution.Statistic;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Tests {

    /// <summary>
    /// 
    /// </summary>
    [TestFixture]
    static public class SolidOnlyTests {


        [Test]
        public static void RotationConvergenceTest([Values(2, 3, 4)] int p = 2
            ) {
            double dt = 1.0e200;
            // --test=ZwoLevelSetSolver.Tests.SolidOnlyTests.RotationConvergenceTest

            BoSSS.Solution.Application.DeleteOldPlotFiles();
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();

            // Displacement-Divergence
            //ZLS_Control.DisplacementDegOffset = 0;
            //ZLS.displacementViscosity = 1.0;
            SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 0.0; // Newton divergence when not 0.0
            SolidPhase.NavierCauchy.EulerAlamansiPenalty = +1.0; // Newton divergence when negative...
            SolidPhase.Continuity.ContinuityInDisplacement = true;
            SolidPhase.Continuity.ContinuityStabilization = true; // seems to be required

            //// Velocity-Divergence
            //ZLS_Control.DisplacementDegOffset = 0;
            //ZLS.displacementViscosity = 0.0;
            //SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 1.0;
            //SolidPhase.NavierCauchy.EulerAlamansiPenalty = 0.0;
            //SolidPhase.Continuity.ContinuityInDisplacement = false;
            //SolidPhase.Continuity.ContinuityStabilization = false; // does not help

            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 8));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 16));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 32));
            if(p <= 2)
                controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 64));

            foreach(var c in controlFiles) {
                Assert.IsTrue(c.SkipSolveAndEvaluateResidual == false);
                c.dtFixed = dt;
                Assert.IsTrue(c.TimesteppingMode == BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady);
                c.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
                c.NonLinearSolver.ConvergenceCriterion = 0.0; // as accurate as possible
                c.NonLinearSolver.MaxSolverIterations = 100; 


            }

            controlFiles.SolverConvergenceTest_Experimental("SolidSolverConvP" + p,
                (VariableNames.DisplacementX, p - 1.5, NormType.L2_embedded),
                (VariableNames.DisplacementY, p - 1.5, NormType.L2_embedded),
                (BoSSS.Solution.NSECommon.VariableNames.Pressure, p - 1.5, NormType.L2noMean_embedded)
                // in a steady-state setting, the velocity is 0.0; therefore, we cannot measure convergence
                //(BoSSS.Solution.NSECommon.VariableNames.VelocityX, p - 1.2, NormType.L2_embedded),
                //(BoSSS.Solution.NSECommon.VariableNames.VelocityY, p - 1.2, NormType.L2_embedded)
                );

        }

        /// <summary>
        /// Here to fix something that works somehow, eventually.
        /// Just one working test, for a positive feeling!
        /// </summary>
        [Test]
        public static void RunSolver([Values(2)] int p = 2,
                                    [Values(16)] int res = 16) {
            // --test=ZwoLevelSetSolver.Tests.SolidOnlyTests.RunSolver
            
            var C = ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, res);
            C.SkipSolveAndEvaluateResidual = false;
            C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 0.0;
            //C.NonLinearSolver.MinSolverIterations = 200;
            C.NonLinearSolver.MaxSolverIterations = 20;

            // Displacement - Divergence
            ZLS.displacementViscosity = 0.0;
            SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 0.0; // Newton divergence when not 0.0
            SolidPhase.NavierCauchy.EulerAlamansiPenalty = +1.0; // Newton divergence when negative...
            SolidPhase.Continuity.ContinuityInDisplacement = true;
            SolidPhase.Continuity.ContinuityStabilization = true;

            ////Velocity - Divergence : this is shit
            //ZLS.displacementViscosity = 0.0;
            //SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 1.0; // must be positive for Newton convergence when Zero
            //SolidPhase.NavierCauchy.EulerAlamansiPenalty = 0.0; // seems to be not required; works for +1, 0, -1
            //SolidPhase.Continuity.ContinuityInDisplacement = false;
            //SolidPhase.Continuity.ContinuityStabilization = false; // seems to have no benefit for condition number


            //C.dtFixed = 0.1;
            //C.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
            

            using(var q = new ZLS()) {
                q.Init(C);
                q.RunSolverMode();
                Assert.IsTrue(q.LastSolverSuccess, "Nonlinear solver did not converge.");
            }
        }

    }
}
