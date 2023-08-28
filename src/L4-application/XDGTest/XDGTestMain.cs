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
using System.Linq;
using System.Text;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using NUnit.Framework;

namespace BoSSS.Application.XDGTest {

    /// <summary>
    /// 
    /// </summary>
    public class XDGTestControl : BoSSS.Solution.Control.AppControl {

        /// <summary>
        /// 
        /// </summary>
        public XDGTestControl() {
            //base.MemoryInstrumentationLevel = MemoryInstrumentationLevel.None;
            SetDGdegree(2);
        }


        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            base.FieldOptions.Clear();
            base.AddFieldOption("Pressure", p, Solution.Control.FieldOpts.SaveToDBOpt.TRUE);
            base.AddFieldOption("Phi", Math.Max(2, p), Solution.Control.FieldOpts.SaveToDBOpt.TRUE);
        }
    }


    /// <summary>
    /// Some basic tests for the XDG framework, e.g. DG projection and extrapolation under level-set movement
    /// </summary>
    class XDGTestMain : BoSSS.Solution.Application<XDGTestControl> {
        static void Main(string[] args) {
            InitMPI();
            //DeleteOldPlotFiles();
            //VariousTests.MultipleTrackerUpdateCalls(1);
            //UnitTest.AllUp();
            UnitTest.RestartTest();
            FinalizeMPI();

            /*
            _Main(args, false, delegate () {
                XDGTestMain p = new XDGTestMain();
                //p.m_GridPartitioningType = BoSSS.Foundation.Grid.GridPartType.none;
                return p;
            });
            */
        }

        /// <summary>
        /// pressure computed in projection step
        /// </summary>
        XDGField Pressure;

        /// <summary>
        /// the level set
        /// </summary>
        LevelSet LevSet;


        /// <summary>
        /// marks all cells which contain species A
        /// </summary>
        SinglePhaseField Amarker;

        /// <summary>
        /// marks all cells which contain species B
        /// </summary>
        SinglePhaseField Bmarker;

        /// <summary>
        /// marks cut, near and far cells.
        /// </summary>
        SinglePhaseField LevelSetDistancce;


        static public double PressureExactA(double[] X) {
            return 2 + 0.3 * X[0] * X[1];
        }
        static public double PressureExactB(double[] X) {
            return 1 - X[1].Pow2();
        }

        static public double Phi(double[] X, double x0) {
            return ((X[0] - 0.83 - x0) / 0.8).Pow2() + (X[1] / 0.8).Pow2() - 1.0;
        }
        static public double Phi0(double[] X) {
            return Phi(X, 0.0);
        }


        protected override void CreateFields() {
            this.LevSet = new LevelSet(new Basis(this.GridData, Control.FieldOptions["Phi"].Degree), "Phi");
            this.LsTrk = new LevelSetTracker((GridData) this.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, LevSet);
            Pressure = new XDGField(new XDGBasis(this.LsTrk, Control.FieldOptions["Pressure"].Degree), "Pressure");
            IOFields.Add(this.LevSet);
            IOFields.Add(this.Pressure);
            Amarker = new SinglePhaseField(new Basis(this.GridData, 0), "SpeciesA");
            Bmarker = new SinglePhaseField(new Basis(this.GridData, 0), "SpeciesB");
            LevelSetDistancce = new SinglePhaseField(new Basis(this.GridData, 0), "LevelSetDistancce");
        }

        protected override int BurstSave => 2;

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
        }

        internal double AutoExtrapolationErr = double.MinValue;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            dt = Control.dtFixed;

          
            Console.WriteLine("Timestep #{0}, dt = {1} ...", TimestepNo, dt);

            // advance level-set
            // -----------------

            double x0 = 0.17 * (phystime + dt);
            this.Pressure.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            this.LevSet.ProjectField((_2D)((x, y) => Phi(new[] { x, y }, x0)));
            this.LsTrk.PushStacks();
            this.LsTrk.UpdateTracker(phystime + dt);


            // update markers
            // --------------
            this.Amarker.Clear();
            this.Bmarker.Clear();
            this.Amarker.AccConstant(1.0, this.LsTrk.Regions.GetSpeciesMask("A"));
            this.Bmarker.AccConstant(1.0, this.LsTrk.Regions.GetSpeciesMask("B"));

            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for(int j = 0; j < J; j++) {
                double dist = this.LsTrk.Regions.GetLevelSetDistance(0, j);
                this.LevelSetDistancce.SetMeanValue(j, dist);
            }

            // test the auto-extrapolation
            // ---------------------------

            var RefPressure = new XDGField(this.Pressure.Basis);
            RefPressure.GetSpeciesShadowField("A").ProjectField(PressureExactA);
            RefPressure.GetSpeciesShadowField("B").ProjectField(PressureExactB);
            RefPressure.Acc(-1.0, Pressure);

            AutoExtrapolationErr = RefPressure.L2Norm();
            Console.WriteLine("Error of extrapolation: " + AutoExtrapolationErr);
            Assert.LessOrEqual(AutoExtrapolationErr, 1.0e-8, "Error after auto-extrapolation to high.");


            Console.WriteLine("done.");

            return dt;
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { this.Pressure, this.LevSet, this.Amarker, this.Bmarker, this.LevelSetDistancce }, "XNSE_prj" + timestepNo, physTime, superSampling);
        }


    }
}
