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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// refinement around the iso-contour MixtureFraction = z_stoichiomeric
    ///
    /// </summary>
    [Serializable]
    public class AMR_onFlameSheet : AMRLevelIndicatorWithLevelset {

        /// <summary>
        /// Empty constructor for serialization
        /// </summary>
        private AMR_onFlameSheet() { }

        public AMR_onFlameSheet(double zSt, int maxRefinementLevelval) {
            m_zSt = zSt;
            maxRefinementLevel = maxRefinementLevelval;
        }


        /// <summary>
        /// The value of the mixture fraction at stoichiometric conditions
        /// </summary>
        [DataMember]
        public double m_zSt;

        [DataMember]
        private DGField MixtureFraction {
            get {
                return (XDGField)((XNSEC)this.SolverMain).IOFields.Where(f => f.Identification == VariableNames.MixtureFraction).Single();
                //return (XDGField)((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MixtureFraction).Single();
            }
        }

        public override int[] DesiredCellChanges() {
           

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            Cell[] cells = GridData.Grid.Cells;


            for(int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;
                double localMinVal, localMaxVal;
                MixtureFraction.GetExtremalValuesInCell(out localMinVal, out localMaxVal, j);
                bool cutcell = false; ;
                try {
                    double lvlsetlocalMinVal, lvlsetlocalMaxVal;
                    var lvlsetImmersedBoundary = (LevelSet)this.LsTrk.LevelSets[1];
                    lvlsetImmersedBoundary.GetExtremalValuesInCell(out lvlsetlocalMinVal, out lvlsetlocalMaxVal, j);
                    cutcell = (lvlsetlocalMaxVal > 0) && (0 > lvlsetlocalMinVal) ? true : false;
                } catch {

                }

                if ((localMaxVal > m_zSt) && (m_zSt > localMinVal) && currentLevel < maxRefinementLevel && !cutcell) {
                    //if((Math.Abs(zMeanValue - m_zSt) < zBand) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                } else if(currentLevel >0) {
                    levels[j] = -1;
                } else {
                    //
                }
                
                //else if ( (localMaxVal < m_zSt || m_zSt < localMinVal) && currentLevel > 0) {
                //    levels[j] = -1;
                //} 

                //else if(!(Math.Abs(zMeanValue - m_zSt) < zBand) && currentLevel > 0) {
                //    levels[j] = -1;
                //}
            }
            Console.WriteLine("Number of refined cells around the flame sheet (z = zst)" + levels.Where(val => val > 0).Sum() + ".\n");
            Console.WriteLine("Number of coarsened cells around the flame sheet (z = zst)" + levels.Where(val => val < 0).Sum() + ".\n");


            //}
            return levels;
        }

        ///// <summary>
        ///// the level-set which represents the fluid-ImmersedBoundary interface
        ///// </summary>
        //protected LevelSet LevSetImmersedBoundary {
        //    get {
        //        return this.LsTrk.LevelSets[1] != null ? (LevelSet)(this.LsTrk.LevelSets[1]) : null;
        //    }
        //}

    }


    /// <summary>
    /// Calculate the reaction rates field and refine in areas with high values
    /// </summary>
    [Serializable]
    public class AMR_onReactiveZones : AMRLevelIndicatorWithLevelset {


        [DataMember]
        public double refinementTreshhold;

        public AMR_onReactiveZones(int maxRefinementLevelval, double _refinementTreshhold) {
            maxRefinementLevel = maxRefinementLevelval;
            refinementTreshhold = _refinementTreshhold;
        }

        /// <summary>
        /// The value of the mixture fraction at stoichiometric conditions
        /// </summary>
        [DataMember]
        public double m_zSt;

     
        [DataMember]
        private XDGField ReactionRate {
            get {
                return (XDGField)((XNSEC)this.SolverMain).Parameters.Where(f => f.Identification == "kReact").Single();
            }
        }

        public override int[] DesiredCellChanges() {

            var solver = ((XNSEC)this.SolverMain);

            solver.XOperator.InvokeParameterUpdate(0.0, solver.CurrentState.ToArray(), solver.Parameters.ToArray()); ;
     


            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            //Tecplot.PlotFields(new DGField[] { kReact }, "Kreact", 0.0, 3);
            double globalMinVal, globalMaxVal;

            ReactionRate.GetExtremalValues(out globalMinVal, out globalMaxVal);

            Cell[] cells = GridData.Grid.Cells;
            for(int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;

                double localMinVal, localMaxVal;
                ReactionRate.GetExtremalValuesInCell(out localMinVal, out localMaxVal, j);
                if ((localMaxVal / globalMaxVal > refinementTreshhold) && currentLevel < maxRefinementLevel && localMaxVal > 1.0) {
                    levels[j] = 1;
                }
                //else if (localMinVal < 0) {
                //    levels[j] = 1;
                //}
                else if ((localMaxVal / globalMaxVal <= refinementTreshhold)) {
                    levels[j] = -1;
                }

            }
            Console.WriteLine("Number of refined cells around areas with high reaction rates: " + levels.Sum() + ".\n");
            return levels;
        }
    }
}