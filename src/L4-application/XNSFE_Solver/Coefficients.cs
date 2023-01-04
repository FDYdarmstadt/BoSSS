using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Application.XNSFE_Solver {

    /// <summary>
    /// EvapMicroRegion, provides a BitArray to tag cells for a Evaporation Model
    /// </summary>
    public class EvapMicroRegion : Coefficient {
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
    /// provides SlipLengths, MultidimensionalArray with slip length for each cell
    /// </summary>
    public class SlipLengths : Coefficient {
        DoNotTouchParameters dntParams;
        ThermalParameters thermParams;

        /// <summary>
        /// SlipLength coefficient for Generalized Navier-Stokes Boundary conditions
        /// </summary>
        public override IList<string> CoefficientsNames => new string[] { "ThermalSlipLengths" };

        /// <summary>
        /// Creates Slip length array, depending on NavierSlip_Localization and NavierSlip_SlipLength
        /// </summary>
        public override DelCoefficientFactory Factory => SlipLengthsFactory;

        /// <summary>
        /// provides SlipLengths, MultidimensionalArray with slip length for each cell
        /// </summary>
        /// <param name="config"></param>
        /// <param name="degU"></param>
        public SlipLengths(IXHeat_Configuration config) {
            dntParams = config.getDntParams;
            thermParams = config.getThermParams;
        }

        (string, object)[] SlipLengthsFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

            CellMask SlipArea;
            switch (dntParams.ThermalSlip_Localization) {
                case ThermalSlip_Localization.Bulk: {
                        SlipArea = lstrk.GridDat.BoundaryCells.VolumeMask;
                        break;
                    }
                case ThermalSlip_Localization.ContactLine: {
                        throw new NotImplementedException();
                        break;
                    }
                case ThermalSlip_Localization.Nearband: {
                        SlipArea = lstrk.GridDat.BoundaryCells.VolumeMask.Intersect(lstrk.Regions.GetNearFieldMask(lstrk.NearRegionWidth));
                        break;
                    }
                case ThermalSlip_Localization.Prescribed: {
                        throw new NotImplementedException();
                    }
                case ThermalSlip_Localization.Everywhere: {
                        SlipArea = CellMask.GetFullMask(lstrk.GridDat);
                        break;
                    }
                default:
                    throw new ArgumentException();
            }

            MultidimensionalArray SlipLengths;
            SlipLengths = lstrk.GridDat.Cells.h_min.CloneAs();
            SlipLengths.Clear();
            //SlipLengths.AccConstant(-1.0);
            if (SlipArea != null) {
                foreach (int i in SlipArea.ItemEnum) {                    
                    SlipLengths[i] = thermParams.sliplength;                        
                }
            }

            var Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], SlipLengths);

            return Ret;
        }

    }

    /// <summary>
    /// provides some (scalar) predefined MassFlux
    /// </summary>
    public class PrescribedMassFlux : Coefficient {
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
