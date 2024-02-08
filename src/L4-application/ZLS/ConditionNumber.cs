using BoSSS.Solution.AdvancedSolvers.Testing;
using CommandLine;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.ControlFiles;
using static BoSSS.Solution.AdvancedSolvers.Testing.ConditionNumberScalingTest;

namespace ZwoLevelSetSolver {
    class ConditionNumber {
        public static void ConditionNumberScaling() {
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();
            //controlFiles.Add(HardCodedControl.Test_Convergence(2,8));
            //controlFiles.Add(HardCodedControl.Test_Convergence(2,16));
            //controlFiles.Add(HardCodedControl.Test_Convergence(2,32));
            controlFiles.Add(ConvergenceTests.RycroftPaper(2, 8));
            controlFiles.Add(ConvergenceTests.RycroftPaper(2, 16));
            controlFiles.Add(ConvergenceTests.RycroftPaper(2, 32));


            // StencilCondNo-bndyUncut
            var conf = new ConditionNumberScalingTest.Config() {
                plot = true,
                title = "ZwoLevelSetSolver-CondScaleing",
                ComputeGlobalCondNo = false
            };

            // 
            conf.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_bndyUncut] = (XAxisDesignation.Grid_1Dres, 0.5, -0.5);


            ConditionNumberScalingTest.Perform(controlFiles, conf);
        }
    }
}
