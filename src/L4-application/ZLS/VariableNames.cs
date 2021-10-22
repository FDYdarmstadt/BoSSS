using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver {
    public static class VariableNames {

        public static string SolidLevelSetCG = BoSSS.Solution.NSECommon.VariableNames.LevelSetCGidx(1);

        public static string SolidLevelSetDG = BoSSS.Solution.NSECommon.VariableNames.LevelSetDGidx(1);

        public static string SolidCurvature = "SolidCurvature";

        public static string DisplacementX = "DisplacementX";

        public static string DisplacementY = "DisplacementY";

        public static string DisplacementZ = "DisplacementZ";

        public static string[] DisplacementVector(int D) {
            if(D == 2)
                return new string[] { DisplacementX, DisplacementY };
            else if(D == 3)
                return new string[] { DisplacementX, DisplacementY, DisplacementZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        public static string DisplacementX0 = "DisplacementX0";

        public static string DisplacementY0 = "DisplacementY0";

        public static string DisplacementZ0 = "DisplacementZ0";

        public static string[] Displacement0Vector(int D) {
            if(D == 2)
                return new string[] { DisplacementX0, DisplacementY0 };
            else if(D == 3)
                return new string[] { DisplacementX0, DisplacementY0, DisplacementZ0 };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        public static string TensionX = "TensionX";

        public static string TensionY = "TensionY";
        
        public static string TensionZ = "TensionZ";

        public static string[] TensionVector(int D) {
            if(D == 2)
                return new string[] { TensionX, TensionY };
            else if(D == 3)
                return new string[] { TensionX, TensionY, TensionZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

    }
}
