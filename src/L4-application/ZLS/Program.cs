using BoSSS.Solution.AdvancedSolvers.Testing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.Tests;

namespace ZwoLevelSetSolver {
    class Program {

        static void Main(string[] args) {
            //RunSolver(args);
            ConditionNumberScaling();
        }

        static void RunSolver(string[] args) {
            ZLS._Main(args, false, delegate () {
                //Control file from runtime via args
                var p = new ZLS();
                return p;
            });
        }

        static void ConditionNumberScaling() {
            ZLS.InitMPI();
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();
            int p = 4;
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 4));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 8));
            controlFiles.Add(ZwoLevelSetSolver.ControlFiles.Vortex.SteadyVortex(p, 16));
            ConditionNumberScalingTest.Perform(controlFiles, true);
        }
    }
}
