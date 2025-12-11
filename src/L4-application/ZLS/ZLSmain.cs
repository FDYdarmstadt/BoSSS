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
using System.Diagnostics;

namespace ZwoLevelSetSolver {
    class ZLSmain {

        static void Main(string[] args) {
            //BoSSS.Solution.Application.InitMPI(num_threads:1);
            //BoSSS.Solution.Application.DeleteOldPlotFiles();

            RunSolver(args);
            //ConditionNumberScaling();
            //Tests.SolidOnlyTests.RotationConvergenceTest(2);
            //Tests.ThreePhaseTests.ThreePhaseContactLine_TensionBalance(true);

            //BoSSS.Solution.Application.FinalizeMPI();
        }

        static void RunSolver(string[] args) {
            //Debugger.Launch();
            ZLS._Main(args, false, delegate () {
                //Control file from runtime via args
                var p = new ZLS();
                return p;
            });
        }
    }

}

