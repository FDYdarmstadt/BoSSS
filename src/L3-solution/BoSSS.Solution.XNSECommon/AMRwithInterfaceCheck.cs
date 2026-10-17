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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Tracing;
using NUnit.Framework.Internal;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.XNSECommon {


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class AMRcurvaturebased : LevelSetTools.SolverWithLevelSetUpdater.AMRLevelIndicatorWithLevelset {

        public int levelSet = 0; // level set this Indicator should be active on

        public double radiusMaxThreshold = 1.0;
        public double curvSignThreshold = 0.1;

        public bool printErrors = false;

        public override int[] DesiredCellChanges() {
            using(var tr = new FuncTrace("AMRcurvaturebased")) {
                tr.InfoToConsole = true;

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];

                CellMask CCmask = this.LsTrk.Regions.GetCutCellMask4LevSet(levelSet);
                SinglePhaseField curvField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "CurvatureByBonnet_AbsoluteValue");
                SinglePhaseField radiusField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "equivalentRadius_AbsoluteValue");
                SinglePhaseField refValueField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "referenceValue_AbsoluteValue");
                SinglePhaseField curvSignField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "CurvatureSign");

                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for(int j = 0; j < J; j++) {
                    int currentLevel = cells[j].RefinementLevel;
                    if(!CCmask.Contains(j)) {
                        if(currentLevel > 0) {
                            levels[j] = -1;
                            cellsToCoarse++;
                        }
                        continue;
                    }

                    var result = EvaluateCurvatureInCutCell(j); //, curvTotField, curvTotFlxField);
                    double radius = 1.0 / result.curv.Abs();
                    double cellSize = GridData.Cells.h_min[j];

                    double refValue = radius / cellSize;
                    double refSign = result.sign.Abs();

                    if(printErrors) {
                        curvField.SetMeanValue(j, result.curv.Abs());
                        radiusField.SetMeanValue(j, radius);
                        refValueField.SetMeanValue(j, refValue);
                        curvSignField.SetMeanValue(j, result.sign);
                    }

                    // check for curvature threshold
                    if(refValue < radiusMaxThreshold && currentLevel < maxRefinementLevel) {
                        tr.Info($"cell {j}: ref value below max radius threshold: {refValue} (threshold = {radiusMaxThreshold})");
                        levels[j] = +1;
                        cellsToRefine++;
                    }
                    //else if(refValue > 2.0 * radiusMaxThreshold && currentLevel > 0) {
                    //    tr.Info($"cell {j}: ref value above coarsen radius threshold: {refValue} (threshold = {2 * radiusMaxThreshold})");
                    //    levels[j] = -1;
                    //    cellsToCoarse++;
                    //}

                    // check for curvature sign threshold (excluding small cut cells)
                    double cellVol = GridData.Cells.GetCellVolume(j);
                    var ccVolumes = LsTrk.GetXDGSpaceMetrics().CutCellMetrics.CutCellVolumes;
                    double ccVA = ccVolumes[LsTrk.GetSpeciesId("A")][j];
                    double ccVB = ccVolumes[LsTrk.GetSpeciesId("B")][j];
                    double minCutCellVol = ccVA < ccVB ? ccVA : ccVB;
                    if(levels[j] == 0 && refSign < curvSignThreshold && currentLevel < maxRefinementLevel && minCutCellVol / cellVol > 0.1) {
                        tr.Info($"cell {j}: ref value below curv sign threshold: {refSign} (threshold = {curvSignThreshold})");
                        levels[j] = +1;
                        cellsToRefine++;
                    }
                    //else if(refSign > 2.0 * radiusMaxThreshold && currentLevel > 0) {
                    //    levels[j] = currentLevel - 1;
                    //    cellsToCoarse++;
                    //}
                }

                if(printErrors) {
                    Tecplot.Tecplot.PlotFields(new DGField[] { curvField, radiusField, refValueField, curvSignField, (DGField)LsTrk.LevelSets[levelSet] }, "AMRcurvature", 0.0, 3);
                }

                return levels;
            }

        }


        private (double curv, double sign) EvaluateCurvatureInCutCell(int jCell) {

            //int order = 0;
            //if(LsTrk.GetCachedOrders().Count() > 0) {
            //    order = LsTrk.GetCachedOrders().Max();
            //} else {
            //    order = 2 * ((DGField)LsTrk.LevelSets[levelSet]).Basis.Degree + 1;
            //}

            var metrics = LsTrk.GetXDGSpaceMetrics();
            //XQuadSchemeHelper schemeHelper = LsTrk.GetXDGSpaceMetrics().XQuadSchemeHelper;

            GridData grdDat = LsTrk.GridDat;
            BitArray cellIntDomBA = new BitArray(grdDat.Cells.NoOfLocalUpdatedCells);
            cellIntDomBA[jCell] = true;
            CellMask cellIntDom = new CellMask(grdDat, cellIntDomBA); //, MaskType.Geometrical);

            double InterfaceArea = 0.0;
            double curvAlongInterface = 0.0;
            double sign = 0.0;

            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                metrics.XQuadSchemeHelper.GetLevelSetQuadScheme(0, cellIntDom).Compile(LsTrk.GridDat, metrics.CutCellQuadOrder),
                delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                    int qN = QR.NoOfNodes;

                    var curv = MultidimensionalArray.Create(length, qN);
                    ((LevelSet)LsTrk.LevelSets[0]).EvaluateTotalCurvature(i0, length, QR.Nodes, curv);
                    // curv.Scale(0.5);    // mean curvature

                    for(int i = 0; i < length; i++) {
                        for(int qn = 0; qn < qN; qn++) {
                            EvalResult[i, qn, 0] = 1.0;
                            EvalResult[i, qn, 1] += curv[i, qn];
                            if(curv[i, qn] >= 0.0) {
                                sign += 1.0;
                            } else {
                                sign -= 1.0;
                            }
                        }
                    }
                    sign /= qN;
                },
                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < length; i++) {
                        InterfaceArea += ResultsOfIntegration[i, 0];
                        curvAlongInterface += ResultsOfIntegration[i, 1];
                    }
                }
            ).Execute();

            return (curvAlongInterface / InterfaceArea, sign);
        }


        public override bool Equals(object obj) {
            if(!base.Equals(obj))
                return false;
            var other = obj as AMRcurvaturebased;
            if(other == null)
                return false;
            if(other.levelSet != this.levelSet)
                return false;
            if(other.radiusMaxThreshold != this.radiusMaxThreshold)
                return false;
            if(other.curvSignThreshold != this.curvSignThreshold)
                return false;
            return true;
        }


    }


    /// <summary>
    /// checks the Gauss theorem and Stokes theorem on cut cells and refines if a given error threshold is exceeded
    /// </summary>
    [Serializable]
    public class AMRwithIntegralTheoremChecks : LevelSetTools.SolverWithLevelSetUpdater.AMRLevelIndicatorWithLevelset {

        public int levelSet = 0; // level set this Indicator should be active on

        public double errThreshold_Gauss = 1.0e-10;
        public double errThreshold_Stokes = 1.0e-10;

        public bool printErrors = false;

        public override int[] DesiredCellChanges() {
            using(var tr = new FuncTrace("AMRcheckIntegralsTheorems")) {
                tr.InfoToConsole = false;

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];

                var testField = GetTestField(GridData);

                CellMask CCmask = this.LsTrk.Regions.GetCutCellMask4LevSet(levelSet);
                SinglePhaseField gaussErrorFieldA = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "gaussError#A");
                SinglePhaseField gaussErrorFieldB = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "gaussError#B");
                SinglePhaseField stokesErrorField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "stokesError");

                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for(int j = 0; j < J; j++) {
                    int currentLevel = cells[j].RefinementLevel;
                    if(!CCmask.Contains(j))
                        continue;

                    int order = 0;
                    if(LsTrk.GetCachedOrders().Count() > 0) {
                        order = LsTrk.GetCachedOrders().Max();
                    } else {
                        order = 2 * ((DGField)LsTrk.LevelSets[levelSet]).Basis.Degree + 1;
                    }

                    var resultGaussA = XNSEUtils.CheckGaussInCutCell(j, LsTrk, LsTrk.GetSpeciesId("A"), testField, order);
                    var resultGaussB = XNSEUtils.CheckGaussInCutCell(j, LsTrk, LsTrk.GetSpeciesId("B"), testField, order);
                    var resultStokes = XNSEUtils.CheckStokesForCell(j, LsTrk, testField, order);


                    if(printErrors) {
                        gaussErrorFieldA.SetMeanValue(j, resultGaussA.error.Abs());
                        gaussErrorFieldB.SetMeanValue(j, resultGaussB.error.Abs());
                        stokesErrorField.SetMeanValue(j, resultStokes.error.Abs());
                    }

                    bool gaussViolated = resultGaussA.error.Abs() > errThreshold_Gauss || resultGaussB.error.Abs() > errThreshold_Gauss;
                    if(gaussViolated) {
                        tr.Info($"Gauss theorem error in cell {j} above threshhold {errThreshold_Gauss}: {resultGaussA.error} (A), {resultGaussB.error} (B)");
                    }

                    bool stokesViolated = resultStokes.error.Abs() > errThreshold_Stokes;
                    if(stokesViolated) {
                        tr.Info($"Stokes theorem error in cell {j} above threshhold {errThreshold_Stokes}: {resultStokes.error} ");
                    }

                    if((gaussViolated || stokesViolated) && currentLevel < maxRefinementLevel) {
                        levels[j] = 1;
                        cellsToRefine++;
                    } else if(!gaussViolated && !stokesViolated && currentLevel > 0) {
                        levels[j] = -1;
                        cellsToCoarse++;
                    }
                }
                tr.Info($"done");

                if(printErrors) {
                    Tecplot.Tecplot.PlotFields(new DGField[] { gaussErrorFieldA, gaussErrorFieldB, stokesErrorField, (DGField)LsTrk.LevelSets[levelSet] }, "AMRintegralErrors", 0.0, 3);
                }

                return levels;
            }

        }


        public VectorField<SinglePhaseField> GetTestField(GridData grdDat) {

            Basis ScalarBasis = new Basis(grdDat, 0);

            SinglePhaseField testFieldX = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fx = (X => 1.0);
            testFieldX.ProjectField(fx);
            SinglePhaseField testFieldY = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fy = (X => 1.0);
            testFieldY.ProjectField(fy);
            SinglePhaseField testFieldZ = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fz = (X => 1.0);
            testFieldZ.ProjectField(fz);

            if(grdDat.SpatialDimension == 2) {
                return new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY });
            } else {
                return new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY, testFieldZ });
            }

        }


        public override bool Equals(object obj) {
            if(!base.Equals(obj))
                return false;
            var other = obj as AMRwithIntegralTheoremChecks;
            if(other == null)
                return false;
            if(other.levelSet != this.levelSet)
                return false;
            if(other.errThreshold_Gauss != this.errThreshold_Gauss)
                return false;
            if(other.errThreshold_Stokes != this.errThreshold_Stokes)
                return false;
            return true;
        }


    }




    /// <summary>
    /// checks the Gauss theorem on cut cells and refines if a given error threshold is exceeded
    /// </summary>
    [Serializable]
    public class AMRwithInterfaceCheck : LevelSetTools.SolverWithLevelSetUpdater.AMRLevelIndicatorWithLevelset {

        public int levelSet = 0; // level set this Indicator should be active on

        public double errThreshold = 1.0e-10;

        public override int[] DesiredCellChanges() {
            using (var tr = new FuncTrace("AMRcheckGauss")) {
                tr.InfoToConsole = true;

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];

                var testField = GetTestField(GridData);

                CellMask CCmask = this.LsTrk.Regions.GetCutCellMask4LevSet(levelSet);

                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for (int j = 0; j < J; j++) {
                    int currentLevel = cells[j].RefinementLevel;
                    if (!CCmask.Contains(j))
                        continue;

                    var result = XNSEUtils.CheckGaussInCutCell(j, LsTrk, LsTrk.GetSpeciesId("A"), testField, LsTrk.GetCachedOrders().Max());
  
                    bool gaussViolated = result.error.Abs() > errThreshold;
                    if (gaussViolated) {
                        tr.Info($"Gauss theorem error in cell {j} above threshhold {errThreshold}: {result.error} ");
                    }
                    if (gaussViolated && currentLevel < maxRefinementLevel) {
                        levels[j] = 1;
                        cellsToRefine++;
                    } else if (!gaussViolated && currentLevel > 0) {
                        levels[j] = -1;
                        cellsToCoarse++;
                    }
                }

            return levels;
            }
            
        }


        public VectorField<SinglePhaseField> GetTestField(GridData grdDat) {

            Basis ScalarBasis = new Basis(grdDat, 0);

            SinglePhaseField testFieldX = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fx = (X => 1.0);
            testFieldX.ProjectField(fx);
            SinglePhaseField testFieldY = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fy = (X => 1.0);
            testFieldY.ProjectField(fy);
            SinglePhaseField testFieldZ = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fz = (X => 1.0);
            testFieldZ.ProjectField(fz);

            if (grdDat.SpatialDimension == 2) {
                return new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY });
            } else {
                return new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY, testFieldZ });
            }

        }


        public override bool Equals(object obj) {
            if (!base.Equals(obj))
                return false;
            var other = obj as AMRwithInterfaceCheck;
            if (other == null)
                return false;
            if (other.levelSet != this.levelSet)
                return false;
            if (other.errThreshold != this.errThreshold)
                return false;
            return true;
        }


    }


}
