/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software.
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled or executed, partly or as a whole, without an explicit
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
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
            Application._Main(args, true, "", delegate () {
                return new Program();
            });
        }

        // Settings
        int dgDegree = 0;

        double inflowDirichletValue = 0.5;
        double rightValue = 1;
        double pod = 0.45;

        int numOfCellsX = 3;
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
        internal int numOfSubgrids = 3;
        internal double energyNorm;

        protected override GridCommons CreateOrLoadGrid() {
            GridCommons grd;

            double[] xnodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
            double[] ynodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
            grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicX: false, periodicY: false);

            return grd;
        }

        protected override void CreateFields() {
            Basis cBasis = new Basis(this.GridData, dgDegree);
            c = new SinglePhaseField(cBasis, "c");
            Velocity = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, new Basis(this.GridData, dgDegree), "Velocity", SinglePhaseField.Factory);
            viscosity = new SinglePhaseField(new Basis(GridData, 2), "viscosity");

            m_IOFields.Add(c);
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
            diffOp = new SpatialOperator(
                new string[] { "c" },
                new string[] { "viscosity", "VelocityX", "VelocityY" },
                new string[] { "codom1" },
                QuadOrderFunc.Linear());
            diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux2D(inflowDirichletValue));
            diffOp.Commit();

            CoordinateMapping coordMap;
            coordMap = new CoordinateMapping(viscosity, Velocity[0], Velocity[1]);

            MultidimensionalArray metricOne = MultidimensionalArray.Create(numOfCellsX);
            MultidimensionalArray metricTwo = MultidimensionalArray.Create(numOfCellsX);

            // 3 sub-grids
            metricOne[0] = 2;
            metricOne[1] = 1;
            metricOne[2] = 0.5;

            metricTwo[0] = 1;
            metricTwo[1] = 0.5;
            metricTwo[2] = 2;

            timeStepper = new AdamsBashforthLTS(diffOp, new CoordinateMapping(c), coordMap, order: ABOrder, numOfSubgrids: this.numOfSubgrids, fluxCorrection: false, reclusteringInterval: 1);

            AdamsBashforthLTS timeStepper2 = timeStepper as AdamsBashforthLTS;
            timeStepper = timeStepper2;

            // Set cell metrics for dynamic clustering
            timeStepper2.MetricOne = metricOne;
            timeStepper2.MetricTwo = metricTwo;
            timeStepper2.ChangeMetricTime = 4e-5;
            timeStepper2.IsNUnitTest = true;

            // Enable sub-grid visualization
            timeStepper2.SgrdField.Identification = "cluster";
            m_IOFields.Add(timeStepper2.SgrdField);
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot plt1 = new Tecplot(GridData, true, false, (uint)superSampling);
            plt1.PlotFields("ALTSTests_" + timestepNo, physTime, m_IOFields);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                base.NoOfTimesteps = int.MaxValue;
                base.EndTime = 10e-5;

                // Set time step size
                if (dt <= 0)
                    dt = 1e-5;

                Console.Write("Timestep " + TimestepNo + " ...");

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
