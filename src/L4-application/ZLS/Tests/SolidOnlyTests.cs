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
        public static void RotationConvergenceTest([Values(2,3,4)] int p = 2
            ) {
            double dt = 1.0e200;
            BoSSS.Solution.Control.NonLinearSolverCode slvCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            // --test=ZwoLevelSetSolver.Tests.SolidOnlyTests.RotationConvergenceTest

            BoSSS.Solution.Application.DeleteOldPlotFiles();
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();

            // Displacement-Divergence
            //ZLS_Control.DisplacementDegOffset = 0;
            ZLS.displacementViscosity = 1.0;
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

                //c.NonLinearSolver.SolverCode = slvCode;
                //if(dt >= 1e10)
                //    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
                //else {
                //    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
                //    c.dtFixed = dt;
                //}


            }

            controlFiles.SolverConvergenceTest_Experimental("SolisSolverConvP" + p,
                (VariableNames.DisplacementX, p - 1.5, NormType.L2_embedded),
                (VariableNames.DisplacementY, p - 1.5, NormType.L2_embedded),
                (BoSSS.Solution.NSECommon.VariableNames.Pressure, p - 1.5, NormType.L2_embedded),
                (BoSSS.Solution.NSECommon.VariableNames.VelocityX, p - 1.2, NormType.L2_embedded),
                (BoSSS.Solution.NSECommon.VariableNames.VelocityY, p - 1.2, NormType.L2_embedded));

        }



    }
}
