using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IntersectingLevelSetTest
{
    class PlotControl : AppControl
    {

        /// <summary>
        /// Ctor.
        /// </summary>
        public PlotControl()
        {
            base.savetodb = false;
            base.ImmediatePlotPeriod = 1;
            base.SuperSampling = 2;
        }


        /// <summary>
        /// nix supported.
        /// </summary>
        public override Type GetSolverType()
        {
            throw new NotSupportedException();
        }

    }
}
