using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StokesHelical_Ak {
    public partial class HelicalMain : Application<HelicalControl> {

        public double cond_Full_0Vars;
        public double cond_Full_12Vars;
        public double cond_Full_012Vars;
        public double cond_Inner_0Vars;
        public double cond_Inner_12Vars;
        public double cond_Inner_012Vars;
        public double cond_OneCell_0Vars;
        public double cond_OneCell_12Vars;
        public double cond_OneCell_012Vars;
        public double cond_Full_AllVars;

    }
}
