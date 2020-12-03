using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {

    class EvapMicroRegion : Coefficient {
        public override IList<string> CoefficientsNames => new string[] { "EvapMicroRegion" };
        public override DelCoefficientFactory Factory => EvapMicroRegionFactory;

        public EvapMicroRegion() {
        }

        (string, object)[] EvapMicroRegionFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            
            BitArray EvapMicroRegion = lstrk.GridDat.GetBoundaryCells().GetBitMask();
            EvapMicroRegion.SetAll(false);
            var Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], EvapMicroRegion);
            
            return Ret;
        }       
    }
}
