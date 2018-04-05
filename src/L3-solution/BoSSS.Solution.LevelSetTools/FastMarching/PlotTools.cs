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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {

    /// <summary>
    /// Plots a single step of the ConstructExtension Algo.
    /// In each step a cell is updated. The Textfile saves the coordinates of each cell center that has been updated.       
    /// </summary>
    class PlotStepper {

        DGField[] Fields;
        int counter_Plotstep;
        int timestepNo;
        int totalSteps;
        GridData GridData;
        Dictionary<int, int> cellStepNumber = new Dictionary<int, int>();

        public PlotStepper() {
            this.counter_Plotstep = 0;
            this.totalSteps = 0;
        }

        public void setup(DGField[] fields, int TimestepNo) {
            GridData = (GridData)(fields.First().GridDat);
            this.Fields = fields;
            totalSteps += counter_Plotstep;
            counter_Plotstep = 0;
            timestepNo = TimestepNo;
        }

        public void plotstep(BitArray Accepted_Mutable) {
            //Create Subgridmask and save to textfile. This can be portrayed in VisIt via a scatter plot. 
            CellMask subgridMask = new CellMask(this.GridData, Accepted_Mutable);
            mapStepNumber(Accepted_Mutable);
            subgridMask.SaveToTextFile("reinit - " + timestepNo + "." + counter_Plotstep + ".txt", true, infoStepNumber);
            //Plot all fields 
            Tecplot.Tecplot.PlotFields(this.Fields, "reinit-" + timestepNo + "." + counter_Plotstep, counter_Plotstep + totalSteps, 3);
            ++counter_Plotstep;

        }

        public double infoStepNumber(double[] x, int jCellLog, int jCelCeom) {
            return cellStepNumber[jCellLog];
        }

        public void mapStepNumber(BitArray Accepted_Mutable) {
            //find new entry and map a number 
            for (int jCell = 0; jCell < Accepted_Mutable.Length; jCell++) {
                if (!cellStepNumber.ContainsKey(jCell) && Accepted_Mutable[jCell] == true) {
                    cellStepNumber[jCell] = counter_Plotstep;
                }
            }
        }

    }

    /*

    List<Tuple<int,double[]>> DependencyVisualizationArrows = new List<Tuple<int, double[]>>();

    void plotDependencyArrow(int cnt, int j1, int j2) {
        double[] pt1 = this.GridDat.Cells.CellCenter.GetRow(j1);
        double[] pt2 = this.GridDat.Cells.CellCenter.GetRow(j2);

        double len = GenericBlas.L2Dist(pt1, pt2);

        double[] ortho = new double[] { pt2[1] - pt1[1], pt1[0] - pt2[0] };
        ortho.ScaleV(0.01 / len);

        int PointsPerStay = 100; // one line ('stay') will be approximated with this number of points.
        double delta = 1.0 / (PointsPerStay - 1);


        double[] stPoint1 = pt1.CloneAs();
        stPoint1.AccV(+1.0, ortho);
        double[] stPoint2 = pt1.CloneAs();
        stPoint2.AccV(-1.0, ortho);

        for(int i = 0; i < PointsPerStay; i++) {
            double s = delta * i;
            double[] pta = new double[2];
            double[] ptb = new double[2];
            pta.AccV(s, stPoint1);
            pta.AccV(1 - s, pt2);
            DependencyVisualizationArrows.Add(new Tuple<int, double[]>(cnt, pta));
            ptb.AccV(s, stPoint2);
            ptb.AccV(1 - s, pt2);
            DependencyVisualizationArrows.Add(new Tuple<int, double[]>(cnt, ptb));
        }
    }

    void PlottAlot(string FileName) {
        using(var outStr = new System.IO.StreamWriter(FileName)) {

            foreach(var tt in DependencyVisualizationArrows) {
                outStr.Write(tt.Item1);
                outStr.Write('\t');

                for(int d = 0; d < 2; d++) {
                    outStr.Write(tt.Item2[d].ToStringDot());
                    if(d < (2 - 1))
                        outStr.Write('\t');
                }

                outStr.WriteLine();
            }

            outStr.Flush();
            outStr.Close();
        }
    }

    void Ploti(CellMask Recalc, CellMask Accept, CellMask Trial, LevelSet Phi, VectorField<SinglePhaseField> Phi_gradient, SinglePhaseField optEik, int cnt) {
        var LSaccept = new LevelSet(this.LevelSetBasis, "Accepted");
        var LStrial = new LevelSet(this.LevelSetBasis, "Trial");
        var LSrecalc = new LevelSet(this.LevelSetBasis, "Recalc");

        LSaccept.AccConstant(-1.0);
        LStrial.AccConstant(-1.0);
        LSrecalc.AccConstant(-1.0);

        LSaccept.Clear(Accept);
        LStrial.Clear(Trial);
        LSrecalc.Clear(Recalc);


        LSaccept.Acc(1.0, Phi, Accept);
        LStrial.Acc(1.0, Phi, Trial);
        LSrecalc.Acc(1.0, Phi, Recalc);

        Tecplot.Tecplot.PlotFields(new DGField[] { Phi, LSaccept, LStrial, LSrecalc, Phi_gradient[0], Phi_gradient[1], optEik }, "reinit-" + cnt, "reinit-" + cnt, cnt, 3);
    }

     */
}
