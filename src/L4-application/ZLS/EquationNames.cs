using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver {
    public static class EquationNames {

        public static string DisplacementEvolutionX = "DisplacementEvolutionX";

        public static string DisplacementEvolutionY = "DisplacementEvolutionY";

        public static string DisplacementEvolutionZ = "DisplacementEvolutionZ";

        static public string DisplacementEvolutionComponent(int d) {
            switch(d) {
                case 0: return DisplacementEvolutionX;
                case 1: return DisplacementEvolutionY;
                case 2: return DisplacementEvolutionZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }
    }
}
