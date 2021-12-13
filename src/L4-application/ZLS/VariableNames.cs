using System;

namespace ZwoLevelSetSolver
{
    public static class VariableNames
    {
        public static string SolidLevelSetCG = BoSSS.Solution.NSECommon.VariableNames.LevelSetCGidx(1);

        public static string SolidLevelSetDG = BoSSS.Solution.NSECommon.VariableNames.LevelSetDGidx(1);

        public static string SolidCurvature = "SolidCurvature";

        public static string DisplacementX = "DisplacementX";

        public static string DisplacementY = "DisplacementY";

        public static string DisplacementZ = "DisplacementZ";

        public static string[] DisplacementVector(int D)
        {
            if (D == 2)
                return new string[] { DisplacementX, DisplacementY };
            else if (D == 3)
                return new string[] { DisplacementX, DisplacementY, DisplacementZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        public static string DisplacementComponent(int d)
        {
            switch (d){
                case 0:
                    return DisplacementX;

                case 1:
                    return DisplacementY;

                case 2:
                    return DisplacementZ;

                default:
                    throw new NotSupportedException();
            }
        }

        public static string DisplacementX0 = "DisplacementX0";

        public static string DisplacementY0 = "DisplacementY0";

        public static string DisplacementZ0 = "DisplacementZ0";

        public static string[] Displacement0Vector(int D)
        {
            if (D == 2)
                return new string[] { DisplacementX0, DisplacementY0 };
            else if (D == 3)
                return new string[] { DisplacementX0, DisplacementY0, DisplacementZ0 };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        public static string DisplacementLaplaceX = "DisplacementLaplaceX";

        public static string DisplacementLaplaceY = "DisplacementLaplaceY";

        public static string[] DisplacementLaplaceVector(int D) {
            if(D == 2)
                return new string[] { DisplacementLaplaceX, DisplacementLaplaceY };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");

        }
    }
}