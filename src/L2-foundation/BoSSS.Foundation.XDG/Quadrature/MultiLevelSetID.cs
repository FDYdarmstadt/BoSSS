using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using System.Threading;
using BoSSS.Platform;
using System.ComponentModel;
using System.Collections;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Foundation.XDG.Quadrature
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
