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
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;

namespace LTS_NUnit {
    /// <summary>
    /// NUnit test project for the implementation of LTS:
    /// The scalar transport equation of the BoSSS Tutorial is used to test LTS.
    /// </summary>
    class Program : Application {
        static void Main(string[] args) {
            Application._Main(args, true, "", delegate () {
                return new Program();
            });
        }

        private int a1 = 10;
        private int a2 = 20;
        private int a3 = 10;
        private int pol_order = 6;

        internal double L2error;
        ExplicitEuler timeStepper;
        DGField u;

        internal int ABorder;
        internal double dt_input;
        internal bool LTS;
        internal bool ALTS;
        internal int numOfSubgrids;

        //For testing
        //internal int ABorder = 2;
        //internal double dt_input = 2E-3 / 4;
        //internal bool LTS = false;
        //internal bool ALTS = true;
        //internal int numOfSubgrids = 1;

        protected override GridCommons CreateOrLoadGrid() {
            double[] xnodes1 = GenericBlas.Linspace(-1, 0, a1 + 1);
            double[] xnodes2 = GenericBlas.Linspace(0, 0.5, a2 + 1);
            double[] xnodes3 = GenericBlas.Linspace(0.5, 1, a3 + 1);

            double[] xComplete = new double[xnodes1.Length + xnodes2.Length + xnodes3.Length - 2];
            for (int i = 0; i < xnodes1.Length; i++) {
                xComplete[i] = xnodes1[i];
            }
            for (int i = 1; i < xnodes2.Length; i++) {
                xComplete[i + xnodes1.Length - 1] = xnodes2[i];
            }
            for (int i = 1; i < xnodes3.Length; i++) {
                xComplete[i + xnodes1.Length + xnodes2.Length - 2] = xnodes3[i];
            }

            double[] ynodes = GenericBlas.Linspace(0, 0.1, 2);
            GridCommons grd = Grid2D.Cartesian2DGrid(xComplete, ynodes, CellType.Square_Linear, true, false);
            grd.Description = "2D cartesian grid 35x1 cells";

            return grd;
        }

        protected override void CreateFields() {
            Basis uBasis = new Basis(this.GridData, pol_order);
            u = new SinglePhaseField(uBasis, "u");
            m_IOFields.Add(u);
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
            SpatialOperator diffOp = new SpatialOperator(1, 0, 1, QuadOrderFunc.MaxDegTimesTwo(), "u", "codom1");
            diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux());
            diffOp.Commit();

            if (LTS) {
                AdamsBashforthLTS ltsTimeStepper = new AdamsBashforthLTS(diffOp, new CoordinateMapping(u), null, ABorder, numOfSubgrids, fluxCorrection: true, reclusteringInterval: 0);
                ltsTimeStepper.SgrdField.Identification = "clusterLTS";
                m_IOFields.Add(ltsTimeStepper.SgrdField);
                timeStepper = ltsTimeStepper;
            } else if (ALTS) {
                AdamsBashforthLTS ltsTimeStepper = new AdamsBashforthLTS(diffOp, new CoordinateMapping(u), null, ABorder, numOfSubgrids, fluxCorrection: false, reclusteringInterval: 1);
                ltsTimeStepper.SgrdField.Identification = "clusterLTS";
                m_IOFields.Add(ltsTimeStepper.SgrdField);
                timeStepper = ltsTimeStepper;
            } else {
                timeStepper = new AdamsBashforth(diffOp, new CoordinateMapping(u), null, ABorder);
                //timeStepper = new RungeKutta(RungeKutta.RungeKuttaScheme.Heun, diffOp, new CoordinateMapping(u),null);
                //timeStepper = RungeKutta.Factory(ABorder, diffOp, new CoordinateMapping(u));
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(m_IOFields, "transport_" + timestepNo, physTime, superSampling);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                if (dt <= 0) {
                    NoOfTimesteps = 10000;
                    EndTime = 2;
                    dt = dt_input;
                    //if (TimestepNo < 3)
                    //    dt /= 3;
                    if (EndTime - phystime < dt)
                        dt = EndTime - phystime;
                }
                if (TimestepNo % 100 == 0)
                    Console.Write("Timestep " + TimestepNo + " ...");
                timeStepper.Perform(dt);
                if (TimestepNo % 100 == 0)
                    Console.WriteLine("finished");

                // Plot and print L2 Error Norm
                if (EndTime - phystime - dt < 1E-10) {
                    L2error = u.L2Error(Jump, 30);
                    Console.WriteLine("L2Error:" + L2error);
                }

                // TEST
                //Console.Write("Timestep " + TimestepNo + " ...");

                //timeStepper.Perform(dt);

                //Console.WriteLine("finished");
                //L2error = u.L2Error(Jump, 30);
                //Console.WriteLine("L2Error:" + L2error);

                //if (TimestepNo == NoOfTimesteps)
                //    Console.ReadKey();

                return dt;
            }
        }

        protected override void SetInitial() {
            u.ProjectField(Jump);
            double error = u.L2Error(Jump, 30);
            Console.WriteLine("L2Error:" + error);
        }

        void Jump(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];
                output[i] = Math.Sin(Math.PI * x);
            }
        }
    }
}
