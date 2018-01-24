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
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.IO;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.Utils;
using ilPSP;
using System.Collections.Generic;
using System.Diagnostics;

namespace ALTSTests {
    /// <summary>
    /// Test program for the implementation of A-LTS.
    /// A scalar transport equation is being used for testing.
    /// In total 1 test run is carried out:
    /// - 1 run with A-LTS featuring a change in the cell metric
    /// to trigger the dynamic (re-)clustering.
    /// </summary>
    class Program : Application {
        static void Main(string[] args) {
            Application._Main(args, true, delegate () {
                Program p = new Program();
                p.m_GridPartitioningType = BoSSS.Foundation.Grid.GridPartType.none;
                return p;
            });
        }

        // Settings
        int dgDegree = 0;

        double inflowDirichletValue = 0.5;
        double rightValue = 1;
        double pod = 0.45;

        //int numOfCellsX = 3;
        int numOfCellsX = 4;
        int numOfCellsY = 1;
        double xMin = 0;
        double xMax = 1;
        double yMin = 0;
        double yMax = 1;
        double xVel = 1;
        double yVel = 0;

        DGField c;
        VectorField<SinglePhaseField> Velocity;
        SinglePhaseField viscosity;
        ExplicitEuler timeStepper;
        SpatialOperator diffOp;

        internal int ABOrder = 1;
        internal int numOfSubgrids = 4;
        internal double energyNorm;

        double dtFixed = 1e-5;
        double endTime = 10e-5;
        private const double ChangeMetricTime = 4e-5;

        protected override GridCommons CreateOrLoadGrid() {
            GridCommons grd;

            double[] xnodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
            double[] ynodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
            grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicX: false, periodicY: false);

            return grd;
        }

        protected override void CreateFields() {
            //Debugger.Launch();

            Basis cBasis = new Basis(this.GridData, dgDegree);
            c = new SinglePhaseField(cBasis, "c");
            Velocity = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, new Basis(this.GridData, dgDegree), "Velocity", SinglePhaseField.Factory);
            viscosity = new SinglePhaseField(new Basis(GridData, 2), "viscosity");

            m_IOFields.Add(c);
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            diffOp = new SpatialOperator(
                new string[] { "c" },
                new string[] { "viscosity", "VelocityX", "VelocityY" },
                new string[] { "codom1" },
                QuadOrderFunc.Linear());
            diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux2D(inflowDirichletValue));
            diffOp.Commit();

            CoordinateMapping coordMap;
            coordMap = new CoordinateMapping(viscosity, Velocity[0], Velocity[1]);

            // 3 sub-grids
            MultidimensionalArray metricOne = MultidimensionalArray.Create(numOfCellsX);
            MultidimensionalArray metricTwo = MultidimensionalArray.Create(numOfCellsX);

            // 3 cells
            //metricOne[0] = 2;
            //metricOne[1] = 1;
            //metricOne[2] = 0.5;

            //metricTwo[0] = 1;
            //metricTwo[1] = 0.5;
            //metricTwo[2] = 2;

            // 4 cells
            metricOne[0] = 2;
            metricOne[1] = 1;
            metricOne[2] = 0.5;
            metricOne[3] = 0.25;

            metricTwo[0] = 0.5;
            metricTwo[1] = 2;
            metricTwo[2] = 0.25;
            metricTwo[3] = 1;

            CustomTimestepConstraint = new SurrogateConstraint(GridData, dtFixed, dtFixed, double.MaxValue, endTime, metricOne, metricTwo);

            timeStepper = new AdamsBashforthLTS(
                diffOp,
                new CoordinateMapping(c),
                coordMap,
                order: ABOrder,
                numOfClusters: this.numOfSubgrids,
                timeStepConstraints: new List<TimeStepConstraint>() { CustomTimestepConstraint },
                fluxCorrection: false,
                reclusteringInterval: 1);

            // Sub-grid visualization
            //AdamsBashforthLTS timeStepper2 = timeStepper as AdamsBashforthLTS;
            //timeStepper2.SubGridField.Identification = "clusterLTS";
            //m_IOFields.Add(timeStepper2.SubGridField);
            //timeStepper = timeStepper2;
        }

        private SurrogateConstraint CustomTimestepConstraint;

        private class SurrogateConstraint : TimeStepConstraint {

            private MultidimensionalArray MetricOne;

            private MultidimensionalArray MetricTwo;

            internal bool FirstMetricActive = true;

            public SurrogateConstraint(GridData gridData, double dtMin, double dtMax, double dtFraction, double EndTime, MultidimensionalArray MetricOne, MultidimensionalArray MetricTwo) :
                base(gridData, dtMin, dtMax, dtFraction, EndTime) {

                this.MetricOne = MetricOne;
                this.MetricTwo = MetricTwo;
            }

            public override double GetLocalStepSize(int i0, int Length) {
                int i0global = (int)gridData.iLogicalCells.GetGlobalID(i0);

                MultidimensionalArray currentMetric;
                if (FirstMetricActive) {
                    currentMetric = MetricOne;
                } else {
                    currentMetric = MetricTwo;
                }

                Debug.Assert(Length == 1);

                return currentMetric[i0global];
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot plt1 = new Tecplot(GridData, true, false, (uint)superSampling);
            plt1.PlotFields("ALTSTests_" + timestepNo, physTime, m_IOFields);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                base.NoOfTimesteps = int.MaxValue;
                base.EndTime = endTime;

                // Set time step size
                if (dt <= 0)
                    dt = dtFixed;

                Console.Write("Timestep " + TimestepNo + " ...");

                if (phystime >= ChangeMetricTime) {
                    CustomTimestepConstraint.FirstMetricActive = false;
                }

                timeStepper.Perform(dt);

                Console.WriteLine("finished");

                // Print L2 Norm at the end of the simulation run
                if (EndTime - phystime - dt < 1E-10) {
                    energyNorm = c.L2Norm();
                    Console.WriteLine("Energy norm = {0:0.000000000000000E-00}", energyNorm);
                }

                return dt;
            }
        }

        protected override void SetInitial() {
            c.ProjectField(delegate (double x, double y) {
                if (x <= pod)
                    return inflowDirichletValue;
                else
                    return rightValue;
            });

            Velocity[0].ProjectField((_2D)((x, y) => xVel));
            Velocity[1].ProjectField((_2D)((x, y) => yVel));

            double energyNormInit = c.L2Norm();
            Console.WriteLine("Energy norm = {0:0.000000000000000E-00}", energyNormInit);
        }

    }
}
