using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Dynamic configuration of multi-level-solvers;
    /// 
    /// </summary>
    public class AllSolveConfig : IterativeSolverConfig {
        public override string Name => throw new NotImplementedException();

        public override string Shortname => "AllSolve";

        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level)
        {
            throw new NotImplementedException();
        }
    }
}
