﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP;

namespace BoSSS.Solution.XNSECommon {
    /// <summary>
    /// provides SlipLengths, MultidimensionalArray with slip length for each cell
    /// </summary>
    public class SlipLengths : Coefficient {
        DoNotTouchParameters dntParams;
        PhysicalParameters physParams;
        int degU;

        /// <summary>
        /// SlipLength coefficient for Generalized Navier-Stokes Boundary conditions
        /// </summary>
        public override IList<string> CoefficientsNames => new string[] { "SlipLengths" };
        
        /// <summary>
        /// Creates Slip length array, depending on NavierSlip_Localization and NavierSlip_SlipLength
        /// </summary>
        public override DelCoefficientFactory Factory => PrescribedMassFluxFactory;

        /// <summary>
        /// Slip length coefficients for Generalized Navier-Stokes Boundary conditions.
        /// </summary>
        /// <param name="config"></param>
        /// <param name="degU"></param>
        public SlipLengths(ISolver_Configuration config, int degU) {
            dntParams = config.getDntParams;
            physParams = config.getPhysParams;
            this.degU = degU;
        }

        (string, object)[] PrescribedMassFluxFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

            CellMask SlipArea;
            switch (dntParams.GNBC_Localization) {
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
                foreach (int i in SlipArea.ItemEnum) {
                    switch (dntParams.GNBC_SlipLength) {
                        case NavierSlip_SlipLength.hmin_DG: {
                            SlipLengths[i] = lstrk.GridDat.Cells.h_min[i] / (degU + 1);
                            break;
                        }
                        case NavierSlip_SlipLength.hmin_Grid: {
                            SlipLengths[i] = SlipLengths[i] = lstrk.GridDat.Cells.h_min[i];
                            break;
                        }
                        case NavierSlip_SlipLength.Prescribed_SlipLength: {
                            SlipLengths[i] = physParams.sliplength;
                            break;
                        }
                        case NavierSlip_SlipLength.Prescribed_Beta: {
                            SlipLengths[i] = -1.0;
                            break;
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
