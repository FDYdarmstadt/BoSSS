using BoSSS.Foundation.XDG.Quadrature.HMF;

namespace BoSSS.Foundation.XDG.Quadrature.BruteForce
{
    struct CombinedID
    {
        public int LevSet0;
        public JumpTypes Jmp0;
        public int LevSet1;
        public JumpTypes Jmp1;

        public bool Equals(CombinedID otherID)
        {
            if ((LevSet0 == otherID.LevSet0)
                && (Jmp0 == otherID.Jmp0)
                && (LevSet1 == otherID.LevSet1)
                && (Jmp1 == otherID.Jmp1))
            {
                return true;
            }
            else if ((LevSet0 == otherID.LevSet1)
                && (Jmp0 == otherID.Jmp1)
                && (LevSet1 == otherID.LevSet0)
                && (Jmp1 == otherID.Jmp0))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
}
