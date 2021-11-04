using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Application.XNSEC {



    /// <summary>
    /// Reynolds number for homotopy 
    /// </summary>
    public class ReynoldsNumber : Coefficient {
        double m_reynolds;
        public override IList<string> CoefficientsNames => new string[] { "Reynolds" };
        public override DelCoefficientFactory Factory => ReynoldsNumberInit;

        public ReynoldsNumber(XNSEC_OperatorConfiguration config) {            
            m_reynolds = config.Reynolds;
        }

        (string, object)[] ReynoldsNumberInit(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {            

            var Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], m_reynolds);

            return Ret;
        }
    }
}
