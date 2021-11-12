using BoSSS.Solution.AdvancedSolvers.Testing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.Tests;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Solution.Statistic;
using NUnit.Framework;

namespace ZwoLevelSetSolver {
    class ZLSmain {

        static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();

            //RunSolver(args);
            //ConditionNumberScaling();
            Tests.SolidOnlyTests.RotationConvergenceTest(2);

            BoSSS.Solution.Application.FinalizeMPI();
        }

        static void RunSolver(string[] args) {
            ZLS._Main(args, false, delegate () {
                //Control file from runtime via args
                var p = new ZLS();
                return p;
            });
        }

        static void ConditionNumberScaling(int p = 2,
            double dt = 1.0e200,
            BoSSS.Solution.Control.NonLinearSolverCode slvCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton
            ) {
            
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();

            // Displacement-Divergence
            //ZLS_Control.DisplacementDegOffset = 0;
            
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
            //controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 64));

            foreach(var c in controlFiles) {
                c.ArtificialViscosity = 1.0;
                c.SkipSolveAndEvaluateResidual = true;
                c.NonLinearSolver.SolverCode = slvCode;
                if(dt >= 1e10)
                    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
                else {
                    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
                    c.dtFixed = dt;
                }

                
            }
            ConditionNumberScalingTest.Perform(controlFiles, true);
        }
        
    }

}

