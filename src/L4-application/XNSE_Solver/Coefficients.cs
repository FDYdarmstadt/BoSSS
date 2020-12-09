using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.XSolver;

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

    /// <summary>
    /// provides SlipLengths, MultidimensionalArray with slip length for each cell
    /// </summary>
    class SlipLengths : Coefficient {
        XNSFE_OperatorConfiguration config;
        int degU;
        public override IList<string> CoefficientsNames => new string[] { "SlipLengths" };
        public override DelCoefficientFactory Factory => PrescribedMassFluxFactory;

        public SlipLengths(XNSFE_OperatorConfiguration config, int degU) {
            this.config = config;
            this.degU = degU;
        }

        (string, object)[] PrescribedMassFluxFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

            CellMask SlipArea;
            switch (config.dntParams.GNBC_Localization) {
                case NavierSlip_Localization.Bulk: {
                    SlipArea = lstrk.GridDat.BoundaryCells.VolumeMask;
                    break;
                }
                case NavierSlip_Localization.ContactLine: {
                    SlipArea = null;
                    break;
                }
                case NavierSlip_Localization.Nearband: {
                    SlipArea = lstrk.GridDat.BoundaryCells.VolumeMask.Intersect(lstrk.Regions.GetNearFieldMask(lstrk.NearRegionWidth));
                    break;
                }
                case NavierSlip_Localization.Prescribed: {
                    throw new NotImplementedException();
                }
                default:
                    throw new ArgumentException();
            }

            MultidimensionalArray SlipLengths;
            SlipLengths = lstrk.GridDat.Cells.h_min.CloneAs();
            SlipLengths.Clear();
            //SlipLengths.AccConstant(-1.0);
            if (SlipArea != null) {
                foreach (Chunk cnk in SlipArea) {
                    for (int i = cnk.i0; i < cnk.JE; i++) {
                        switch (config.dntParams.GNBC_SlipLength) {
                            case NavierSlip_SlipLength.hmin_DG: {                                
                                SlipLengths[i] = lstrk.GridDat.Cells.h_min[i] / (degU + 1);
                                break;
                            }
                            case NavierSlip_SlipLength.hmin_Grid: {
                                SlipLengths[i] = SlipLengths[i] = lstrk.GridDat.Cells.h_min[i];
                                break;
                            }
                            case NavierSlip_SlipLength.Prescribed_SlipLength: {
                                SlipLengths[i] = config.physParams.sliplength;
                                break;
                            }
                            case NavierSlip_SlipLength.Prescribed_Beta: {
                                SlipLengths[i] = -1.0;
                                break;
                            }
                        }
                    }
                }
            }

            var Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], SlipLengths);

            return Ret;
        }

    }
}
