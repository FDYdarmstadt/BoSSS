using BoSSS.Solution.AdvancedSolvers.Testing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.ControlFiles;

namespace ZwoLevelSetSolver
{
    class ConditionNumber
    {
        public static void ConditionNumberScaling()
        {
            List<ZLS_Control> controlFiles = new List<ZLS_Control>();
            //controlFiles.Add(HardCodedControl.Test_Convergence(2,8));
            //controlFiles.Add(HardCodedControl.Test_Convergence(2,16));
            //controlFiles.Add(HardCodedControl.Test_Convergence(2,32));
            controlFiles.Add(ConvergenceTests.RycroftPaper(2, 8));
            controlFiles.Add(ConvergenceTests.RycroftPaper(2, 16));
            controlFiles.Add(ConvergenceTests.RycroftPaper(2, 32));
            ConditionNumberScalingTest.Perform(controlFiles);
        }
    }
}
