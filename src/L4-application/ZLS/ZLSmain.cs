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
            //BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();

            RunSolver(args);
            //ConditionNumberScaling();
            //Tests.SolidOnlyTests.RotationConvergenceTest(2);
            //ZwoLevelSetSolver.Tests.SolidOnlyTests.RotationConvergenceTest(2);

            //ParameterSweep();
            //BoSSS.Solution.Application.FinalizeMPI();
        }

        static void RunSolver(string[] args) {

            //SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 0.0; // Newton divergence when not 0.0
            //SolidPhase.NavierCauchy.EulerAlamansiPenalty = +1.0; // Newton divergence when negative...
            //SolidPhase.Continuity.ContinuityInDisplacement = true;
            //SolidPhase.Continuity.ContinuityStabilization = true; // seems to be required

            //var C = ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(2, 32);
            //C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            //C.NonLinearSolver.ConvergenceCriterion = 0.0; // as accurate as possible
            //C.NonLinearSolver.MaxSolverIterations = 100;


            //using (var q = new ZLS()) {
            //    q.Init(C);
            //    q.RunSolverMode();
            //    //q.OperatorAnalysis();
            //}
            //Assert.AreEqual(true, false, "Remove me");
                        
            
            ZLS._Main(args, false, delegate () {
                //Control file from runtime via args
                var p = new ZLS();
                return p;
            });
            //*/
        }

        static void ConditionNumberScaling(int p = 2,
            double dt = 1.0e200,
            BoSSS.Solution.Control.NonLinearSolverCode slvCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton
            ) {
            
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
            //controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 64));

            foreach(var c in controlFiles) {
                c.SkipSolveAndEvaluateResidual = true;
                c.NonLinearSolver.SolverCode = slvCode;
                if(dt >= 1e10)
                    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
                else {
                    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
                    c.dtFixed = dt;
                }

                
            }
            ConditionNumberScalingTest.Perform(controlFiles, new ConditionNumberScalingTest.Config() { plot = true, title= "ConditionNumberScaling" });
        }



        static void ConvergenceNumberScaling(int p = 2,
            double dt = 1.0e200,
            BoSSS.Solution.Control.NonLinearSolverCode slvCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton
            ) {
            
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
            //controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 64));

            foreach(var c in controlFiles) {
                c.SkipSolveAndEvaluateResidual = true;
                c.NonLinearSolver.SolverCode = slvCode;
                if(dt >= 1e10)
                    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
                else {
                    c.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
                    c.dtFixed = dt;
                }

                
            }
            ConditionNumberScalingTest.Perform(controlFiles, new ConditionNumberScalingTest.Config() { plot = true, title = "ConvergenceNumberScaling" });
        }

          /*
        class Case {
            public Dictionary<string, object> Keys;
            public Dictionary<string, double> Results;


            public void AppendToTable() {

            }
        }

        static void ParameterSweep() {

            int NoOfCases = 0;
            for(int sweep = 0; sweep < 2; sweep++) {
                int cnt = 0;


                foreach(int DisplOffset in new int[] { -2, -1, 0, 1 }) { // 1: DG-degree offset of displacements
                    ZLS_Control.DisplacementDegOffset = DisplOffset;
                    foreach(double av in new double[] { 0, 0.001,  1 }) { // 2: Displacement artificial viscosity
                        ZLS.displacementViscosity = av;
                        foreach(double dt in new double[] { 1, 10, 1e300 }) { // 3: Timestep size
                                                                                         // kkkk
                            foreach(bool contiInDispl in new bool[] { false, true }) { // 4: Divergence-freenes of displacements or velocities?
                                SolidPhase.Continuity.ContinuityInDisplacement = contiInDispl;
                                foreach(bool contiStab in new bool[] { false, true }) { // 5: Stabilization on/off
                                    SolidPhase.Continuity.ContinuityStabilization = contiStab;
                                    foreach(var slvCode in new[] { BoSSS.Solution.Control.NonLinearSolverCode.Newton, BoSSS.Solution.Control.NonLinearSolverCode.Picard }) { // 6: newton vs picard
                                                                                                                                                                             //
                                        foreach(double penalty2 in new double[] { -1, 0, 1 }) { // additional displacement penalty
                                            SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = penalty2;

                                            foreach(double penalty3 in new double[] { -1.3, 0, 1.3 }) { // penalty in Euler-Almansi stress term
                                                SolidPhase.NavierCauchy.EulerAlamansiPenalty = penalty3;
                                                
                                                cnt++;

                                                if(sweep == 0) {
                                                    NoOfCases++;
                                                } else {
                                                    Console.WriteLine($"Running case {cnt} of {NoOfCases}...");
                                                }

                                            }

                                        }
                                    }
                                }

                            }


                        }
                    }
                }
            }

        }
        */
    }

}

