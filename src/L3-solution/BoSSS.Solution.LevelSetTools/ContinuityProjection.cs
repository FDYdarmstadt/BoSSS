/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution.LevelSetTools {

    /// <summary>
    /// Options for enforcing the continuity of the level-set function
    /// </summary>
    public enum ContinuityProjectionOption {

        /// <summary>
        /// Do not perform continuity projection
        /// </summary>
        None,

        /// <summary>
        /// projection on a spectral finite element field
        /// </summary>
        SpecFEM,

        /// <summary>
        /// L2-projection with continuity constraints at inner cell boundaries
        /// </summary>
        ContinuousDG
    }


    /// <summary>
    /// Projects a DG Field onto another DGField of higher order
    /// </summary>
    public class ContinuityProjection {
        /// <summary>
        /// Selects algorithm variant based on the <paramref name="Option"/>
        /// </summary>
        /// <param name="DGLevelSet">Discontinuous Field</param>
        /// <param name="gridData"></param>
        /// <param name="Option">Choice of algorithm</param>
        /// <param name="SmoothedLevelSet">The Continuous Field</param>
        public ContinuityProjection(SinglePhaseField DGLevelSet, Foundation.Grid.Classic.GridData gridData, ContinuityProjectionOption Option) {
            int k = DGLevelSet.Basis.Degree + 1;
            myOption = Option;
            switch (Option) {
                case ContinuityProjectionOption.SpecFEM: {
                        var ContinuousLevelSetBasis = new SpecFemBasis(gridData, k);
                        MyProjection = new ContinuityProjectionSpecFem(ContinuousLevelSetBasis);
                        break;
                    }
                case ContinuityProjectionOption.ContinuousDG: {
                        var ContinuousLevelSetDGBasis = new Basis(gridData, k);
                        MyProjection = new ContinuityProjectionCDG(ContinuousLevelSetDGBasis);
                        break;
                    }
                case ContinuityProjectionOption.None: {
                        MyProjection = new NoProjection();
                        break;
                    }
                default:
                    throw new ArgumentException();
            }

        }

        /// <summary>
        /// Creates the LevelSet field for the Tracker based on the Option selected
        /// </summary>
        /// <param name="DGLevelSet">the unfiltered Level-set</param>
        /// <param name="gridData"></param>
        /// <param name="Option"></param>
        /// <returns>The filtered Level-Set Field </returns>
        public static LevelSet CreateField(SinglePhaseField DGLevelSet, Foundation.Grid.Classic.GridData gridData, ContinuityProjectionOption Option) {
            int k = DGLevelSet.Basis.Degree + 1;
            switch (Option) {
                case ContinuityProjectionOption.SpecFEM: {
                        var ContinuousLevelSetBasis = new SpecFemBasis(gridData, k);
                        return new LevelSet(ContinuousLevelSetBasis.ContainingDGBasis, "Phi");
                    }
                case ContinuityProjectionOption.ContinuousDG: {
                        var ContinuousLevelSetDGBasis = new Basis(gridData, k);
                        return new LevelSet(ContinuousLevelSetDGBasis, "Phi");
                    }
                case ContinuityProjectionOption.None: {
                        Console.WriteLine("WARNING: No additional enforcement of the level-set continuity!");
                        LevelSet SmoothedLevelSet = new LevelSet(DGLevelSet.Basis, "Phi");
                        DGLevelSet = SmoothedLevelSet;
                        return SmoothedLevelSet;
                    }
                default:
                    throw new ArgumentException();
            }
        }


        IContinuityProjection MyProjection;

        ContinuityProjectionOption myOption;

        /// <summary>
        /// Makes <paramref name="DGLevelSet"/> a continuous function <paramref name="LevelSet"/>,
        /// according to the Option from the initialization
        /// </summary>
        /// <param name="DGLevelSet"></param>
        /// <param name="LevelSet"></param>
        /// <param name="Domain"></param>
        /// <param name="PosMask"></param>
        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain, CellMask PosMask, bool setFarFieldConstant = true) {
            MyProjection.MakeContinuous(DGLevelSet, LevelSet, Domain);

            if (myOption != ContinuityProjectionOption.None && setFarFieldConstant) {
                SetFarField(LevelSet, Domain, PosMask);
            }

        }

        /// <summary>
        /// Sets the positive Far-field of the level-set to +1 and the negative side to -1 
        /// </summary>
        /// <param name="LevelSet"></param>
        /// <param name="Domain"></param>
        /// <param name="PosMask"></param>
        public void SetFarField(SinglePhaseField LevelSet, CellMask Domain, CellMask PosMask) {
            // set cells outside narrow band to +/- 1
            var Pos = PosMask.Except(Domain);
            var Neg = PosMask.Complement().Except(Domain);

            LevelSet.Clear(Pos);
            LevelSet.AccConstant(1, Pos);

            LevelSet.Clear(Neg);
            LevelSet.AccConstant(-1, Neg);
        }

    }



   
    interface IContinuityProjection {
        void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain);
    }

    // ===============================
    // The actual implementations
    // ===============================

    ///<summary>
    /// Smoothing based on SpecFem
    /// => Actually ContinuousFunctionSpace
    ///</summary>
    class ContinuityProjectionSpecFem : IContinuityProjection  {

        public ContinuityProjectionSpecFem(SpecFemBasis myBasis) {
            FEMLevSet = new SpecFemField(myBasis);
        }
        SpecFemField FEMLevSet;

        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain) {
            FEMLevSet.ProjectDGField(1.0, DGLevelSet, Domain);
            LevelSet.Clear();
            FEMLevSet.AccToDGField(1.0, LevelSet);
        }
    }

    ///<summary>
    /// Smoothing based on ContinuousDGField 
    /// => Lagrange-Multiplier Approach
    ///</summary>
    class ContinuityProjectionCDG : IContinuityProjection {

        public ContinuityProjectionCDG(Basis myBasis) {
            CDGField = new ContinuousDGField(myBasis);
        }
        ContinuousDGField CDGField;

        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain) {
            CDGField.ProjectDGField(1.0, DGLevelSet, Domain);
            LevelSet.Clear();
            CDGField.AccToDGField(1.0, LevelSet);
        }
    }

    ///<summary>
    /// Does nothing => no enforcement of continuity  
    ///</summary>
    class NoProjection : IContinuityProjection {

        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain) {

            if (ReferenceEquals(DGLevelSet, LevelSet)) {
                // Do nothing
            }
            else{
                LevelSet.Clear();
                LevelSet.AccLaidBack(1.0, DGLevelSet);
            }
            
        }
    }

}
