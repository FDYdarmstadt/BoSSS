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
            //RunSolver(args);
            ConditionNumberScaling();
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
            ZLS.displacementViscosity = 1.0;
            SolidPhase.Continuity.ContinuityInDisplacement = false;
            SolidPhase.Continuity.ContinuityStabilization = false;
            
       
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 8));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 16));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 32));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 64));

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

        static void ParameterSweep() {
            foreach(int DisplOffset  in new int[] {-2, -1, 0, 1}) { // 1: DG-degree offset of displacements
                ZLS_Control.DisplacementDegOffset = DisplOffset;
                foreach(double av in new double[] {0, 0.001, 0.01, 0.1, 1}) { // 2: Displacement artificial viscosity
                    ZLS.displacementViscosity = av;
                    foreach(double dt in new double[] {1, 10, 100, 1000, 1e300}) { // 3: Timestep size
                        // kkkk
                        foreach(bool contiInDispl in new bool[] { false, true }) { // 4: Divergence-freenes of displacements or velocities?
                            SolidPhase.Continuity.ContinuityInDisplacement = contiInDispl;
                            foreach(bool contiStab in new bool[] { false, true }) { // 5: Stabilization on/off
                                SolidPhase.Continuity.ContinuityStabilization = contiStab; 
                                foreach(var slvCode in new[] { BoSSS.Solution.Control.NonLinearSolverCode.Newton, BoSSS.Solution.Control.NonLinearSolverCode.Picard }) { // 6: newton vs picard
                                    //
                                }
                            }

                        }


                    }
                }
            }
        }

    }

}

