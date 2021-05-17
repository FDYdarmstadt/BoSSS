using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// control object for <see cref="XNSFE"/>
    /// </summary>
    public class XNSFE_Control : XNSE_Control {


        public XNSFE_Control() : base() {
            this.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton; // XNSFE should always use Newton solver
        }

        /// <summary>
        /// type for solver factory
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XNSFE);
        }        
    }
}
