using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// EvapMicroRegion, provides a BitArray to tag cells for a Evaporation Model
    /// </summary>
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

    /// <summary>
    /// provides some (scalar) predefined MassFlux
    /// </summary>
    class PrescribedMassFlux : Coefficient {
        XNSFE_OperatorConfiguration config;
        public override IList<string> CoefficientsNames => new string[] { "PrescribedMassFlux" };
        public override DelCoefficientFactory Factory => PrescribedMassFluxFactory;

        public PrescribedMassFlux(XNSFE_OperatorConfiguration config) {
            this.config = config;
        }

        (string, object)[] PrescribedMassFluxFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

            double[] dummyX = new double[] { 0.0, 0.0 };
            double PrescribedMassFlux = config.prescribedMassflux(dummyX, time);

            var Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], PrescribedMassFlux);

            return Ret;
        }
    }
}
