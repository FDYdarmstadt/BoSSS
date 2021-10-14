using BoSSS.Solution.AdvancedSolvers.Testing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.Tests;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;


namespace ZwoLevelSetSolver {
    class ZLSmain {

        static void Main(string[] args) {
            RunSolver(args);
            //ConditionNumberScaling();

            //ParameterSweep();
        }

        static void RunSolver(string[] args) {
            /*
            BoSSS.Solution.Application.InitMPI();

            const int res = 8;

            var C0 = BoSSS.Application.SipPoisson.SipHardcodedControl.Square(res, res, 2);
            using(var r = new BoSSS.Application.SipPoisson.SipPoissonMain()) {
                r.Init(C0);
                r.RunSolverMode();


                var MTX = BoSSS.Application.SipPoisson.SipPoissonMain.LastMatrix;
                long J = r.GridData.CellPartitioning.TotalLength;
                MTX.SaveToTextFileSparse($"Mtx_ipPoisson-J{J}.txt");
                double condNo = MTX.condest();
                Console.WriteLine($"Matlab condition number estimate {J} cells: " + condNo);

                
            }

            var C = ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(2, res);
            C.SkipSolveAndEvaluateResidual = true;
            C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            using(var q = new ZLS()) {
                q.Init(C);
                q.RunSolverMode();
                q.OperatorAnalysis();
            }

            BoSSS.Solution.Application.FinalizeMPI();
            */


            BoSSS.Solution.Application.InitMPI();

            const int res = 4;

            var C = ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(2, res);
            C.SkipSolveAndEvaluateResidual = false;
            C.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            C.NonLinearSolver.MinSolverIterations = 200;
            C.NonLinearSolver.MaxSolverIterations = 220;

            ZLS_Control.DisplacementDegOffset = 0;
            ZLS.displacementViscosity = 0.0;
            SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 0.0;
            SolidPhase.NavierCauchy.EulerAlamansiPenalty = +1.0; // divergence when negative...
            SolidPhase.Continuity.ContinuityInDisplacement = true;
            SolidPhase.Continuity.ContinuityStabilization = true;


            //ZLS_Control.DisplacementDegOffset = 0;
            //ZLS.displacementViscosity = 0.0;
            //SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 1.3;
            //SolidPhase.NavierCauchy.EulerAlamansiPenalty = +1.0; // divergence when negative...
            //SolidPhase.Continuity.ContinuityInDisplacement = false;
            //SolidPhase.Continuity.ContinuityStabilization = false;


            //C.dtFixed = 1.0;
            C.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            using(var q = new ZLS()) {
                q.Init(C);
                q.RunSolverMode();
                //q.OperatorAnalysis();
            }

            BoSSS.Solution.Application.FinalizeMPI();

            /*
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
            ZLS.InitMPI();
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();

            ZLS_Control.DisplacementDegOffset = 0;
            ZLS.displacementViscosity = 0.0;
            SolidPhase.DisplacementEvolution.onlyPenaltyPenalty = 0.0;
            SolidPhase.NavierCauchy.EulerAlamansiPenalty = 1.0;
            SolidPhase.Continuity.ContinuityInDisplacement = true;
            SolidPhase.Continuity.ContinuityStabilization = true;

            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 8));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 16));
            //controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 32));
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
            ConditionNumberScalingTest.Perform(controlFiles, true);
        }

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
    }

}

