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
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// refinement of cells based on ALL residuals
    /// </summary>
    [Serializable]
    public class AMR_BasedOnResiduals : AMRLevelIndicatorWithLevelset {
 
  
        public AMR_BasedOnResiduals() {
        }

        public AMR_BasedOnResiduals(int _maxRefinementLevelval, double _tresholdValue) {
            maxRefinementLevel = _maxRefinementLevelval;
            tresholdValue = _tresholdValue;
        }

        public double tresholdValue;


        protected IList<DGField> FieldForRefinement {
            get {
                return ((XNSEC)this.SolverMain).CurrentResidualVector.Fields;
                //return ((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == FieldName).Single();
            }
        }

        public int delayCounter = 0;
        public override int[] DesiredCellChanges() {
            //==========================================
            //Refinement based on residuals
            //==========================================
            int FirstAMR = 5;
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            foreach(var field in FieldForRefinement) {

                Cell[] cells = GridData.Grid.Cells;

                double globalMinVal, globalMaxVal;

                field.GetExtremalValues(out globalMinVal, out globalMaxVal);
                
                double maxGlobalAbsoluteValue = Math.Max(Math.Abs(globalMinVal), Math.Abs(globalMaxVal));


                for(int j = 0; j < J; j++) {
                    double localMinVal, localMaxVal;
                    field.GetExtremalValuesInCell(out localMinVal, out localMaxVal, j);

                    double maxAbsoluteValue = Math.Max(Math.Abs(localMinVal), Math.Abs(localMaxVal));



                    int currentLevel = cells[j].RefinementLevel;
                    if( maxAbsoluteValue/maxGlobalAbsoluteValue >0.1 &&
                        maxAbsoluteValue > 1e-2 && 
                        currentLevel < maxRefinementLevel
                        //&& delayCounter > FirstAMR
                        ) {
                        levels[j] = 1;
                    }

                    //if(currentLevel > 0 && maxAbsoluteValue < 1e-2) {
                    //    levels[j] = -1;
                    //}

                    
                }
            }

            return levels;
        }
    }




    /// <summary>
    /// refinement of cells based on a residuals
    /// </summary>
    [Serializable]
    public class AMR_BasedOnResidual_GRADIENT : AMRLevelIndicatorWithLevelset {


        public AMR_BasedOnResidual_GRADIENT() {
        }

        public AMR_BasedOnResidual_GRADIENT(int _maxRefinementLevelval, double _tresh) {
            maxRefinementLevel = _maxRefinementLevelval;
            this.tresh = _tresh;
        }

        protected IList<DGField> FieldForRefinement {
            get {
                return ((XNSEC)this.SolverMain).CurrentResidualVector.Fields;
                //return ((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == FieldName).Single();
            }
        }
        public  double tresh = 1;
        public int delayCounter = 0;
        public override int[] DesiredCellChanges() {
            //==========================================
            //Refinement based on residuals
            //==========================================
            int FirstAMR = 3;
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            int D = GridData.SpatialDimension;

            foreach(var field in FieldForRefinement) {

                Cell[] cells = GridData.Grid.Cells;


                // Calculate magnitude gradient of variable
                VectorField<XDGField> GradientOfField = new VectorField<XDGField>(D.ForLoop(d => new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Temperature].Degree), "Grad[" + d + "]")));
                GradientOfField.Gradient(1.0, field);

                var MagnitudeGradient = new SinglePhaseField(new Basis(base.GridData, 0), "Magnitude_Grad");
                MagnitudeGradient.Clear();
                MagnitudeGradient.ProjectFunction(1.0,
                    (ilPSP.Vector X, double[] U, int jCell) => Math.Sqrt(U[0].Pow2() + U[1].Pow2()),
                    new Foundation.Quadrature.CellQuadratureScheme(),
                    GradientOfField.ToArray());



                double globalMinVal, globalMaxVal;
                MagnitudeGradient.GetExtremalValues(out globalMinVal, out globalMaxVal);
                double maxGlobalAbsoluteValue = Math.Max(Math.Abs(globalMinVal), Math.Abs(globalMaxVal));


                for(int j = 0; j < J; j++) {
                    double localMinVal, localMaxVal;
                    MagnitudeGradient.GetExtremalValuesInCell(out localMinVal, out localMaxVal, j);

                    double maxAbsoluteValue = Math.Max(Math.Abs(localMinVal), Math.Abs(localMaxVal));

                    int currentLevel = cells[j].RefinementLevel;
                    if(maxAbsoluteValue > tresh &&
                        currentLevel < maxRefinementLevel &&
                        delayCounter > FirstAMR &&
                        field.Identification == VariableNames.Temperature
                        ) {
                        levels[j] = 1;
                    }

                    //if(currentLevel > 0 && maxAbsoluteValue < 1e-2) {
                    //    levels[j] = -1;
                    //}


                }
            }
            Console.WriteLine("delayCounter " + delayCounter);
            delayCounter++;
            return levels;
        }
    }





}