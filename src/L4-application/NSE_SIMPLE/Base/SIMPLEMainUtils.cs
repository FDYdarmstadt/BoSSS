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
using System.IO;
using System.Globalization;
using BoSSS.Platform;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using NSE_SIMPLE.LowMach;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Some utilities for special test cases.
    /// </summary>
    public partial class NSE_SIMPLEMain {

        #region OrrSommerfeld

        /// <summary>
        /// Write quadrature nodes to text file for Matlab.
        /// </summary>
        void WriteQuadNodesOrrSommerfeld() {
            StreamWriter Nodes = new StreamWriter("\\\\fdyprime\\redirected\\klein\\Documents\\BoSSS\\Orr-Sommerfeld\\Nodes32x80.txt");

            this.WorkingSet.Velocity.Current[0].ProjectField(delegate(MultidimensionalArray input, MultidimensionalArray output) {
                for (int i = 0; i < input.GetLength(0); i++) {
                    double x = input[i, 0];
                    double y = input[i, 1];

                    Nodes.Write(x.ToString("E", System.Globalization.CultureInfo.InvariantCulture).PadLeft(14));
                    Nodes.Write("\t" + y.ToString("E", System.Globalization.CultureInfo.InvariantCulture).PadLeft(14));
                    Nodes.Write("\n");
                }
            });

            Nodes.Close();

            Console.WriteLine("Write QuadNodes finished.");
        }

        /// <summary>
        /// Initialize Orr-Sommerfeld
        /// </summary>
        void InitOrrSommerfeld() {
            string[] fileLines = File.ReadAllLines("\\\\fdyprime\\redirected\\klein\\Documents\\BoSSS\\Orr-Sommerfeld\\uv_init_32x80.txt");
            double[,] uv_init = new double[fileLines.Length, 2];

            for (int i = 0; i < fileLines.Length; i++) {
                string[] uv_split = fileLines[i].Split('\t');
                uv_init[i, 0] = Convert.ToDouble(uv_split[0], CultureInfo.InvariantCulture);
                uv_init[i, 1] = Convert.ToDouble(uv_split[1], CultureInfo.InvariantCulture);
            }

            int i0 = 0;
            WorkingSet.Velocity.Current[0].ProjectField(delegate(MultidimensionalArray input, MultidimensionalArray output) {
                for (int i = 0; i < input.GetLength(0); i++) {
                    output[i] = uv_init[i + i0, 0];
                }
                i0 += input.GetLength(0);
            });

            i0 = 0;
            WorkingSet.Velocity.Current[1].ProjectField(delegate(MultidimensionalArray input, MultidimensionalArray output) {
                for (int i = 0; i < input.GetLength(0); i++) {
                    output[i] = uv_init[i + i0, 1];
                }
                i0 += input.GetLength(0);
            });

            Console.WriteLine("Initialized OrrSommerfeld.");
        }

        // Some variables for OrrSommerfeld
        SinglePhaseField Energy;
        TextWriter Log_Energy;
        double EnergyIntegral;
        double Energy_t0;
        double omega;

        /// <summary>
        /// Initialize logging of energy for Orr-Sommerfeld problem
        /// </summary>
        void InitLogEnergyOrrSommerfeld() {
            Energy = new SinglePhaseField(new Basis(base.gridData, 20));

            if (base.MPIRank == 0) {
                Log_Energy = base.DatabaseDriver.FsDriver.GetNewLog("PerturbationEnergy", CurrentSessionInfo.ID);
            }

            Energy.Clear();
            Energy.ProjectFunction(
                1.0,
                (X, U, cell) => (1.0 - X[1] * X[1] - U[0]) * (1.0 - X[1] * X[1] - U[0]) + U[1] * U[1],
                null,
                WorkingSet.Velocity.Current[0],
                WorkingSet.Velocity.Current[1]);

            Energy_t0 = Energy.IntegralOver(null);

            if (base.MPIRank == 0) {
                Log_Energy.Write("Time\tPerturbationEnergy\tOmega");
                Log_Energy.WriteLine();
                Log_Energy.Write("{0}\t{1}\t{2}",
                    (0.0).ToString("0.000", NumberFormatInfo.InvariantInfo),
                    Energy_t0.ToString("0.0000000000E+00", NumberFormatInfo.InvariantInfo),
                    (2.234976E-03).ToString("0.000000E+00", NumberFormatInfo.InvariantInfo));
                Log_Energy.Flush();
            }

            if (base.MPIRank == 0) {
                Console.WriteLine("E(t0) = " + Energy_t0);
            }
        }

        /// <summary>
        /// Log energy for Orr-Sommerfeld problem.
        /// </summary>
        void LogEnergyOrrSommerfeld(int TimestepNo, double phystime, double dt) {
            Energy.Clear();
            Energy.ProjectFunction(
                1.0,
                (X, U, cell) => (1.0 - X[1] * X[1] - U[0]) * (1.0 - X[1] * X[1] - U[0]) + U[1] * U[1],
                null,
                WorkingSet.Velocity.Current[0],
                WorkingSet.Velocity.Current[1]);

            EnergyIntegral = Energy.IntegralOver(null);
            omega = 1 / (2 * (phystime + dt)) * Math.Log(EnergyIntegral / Energy_t0);

            if (base.MPIRank == 0) {
                Log_Energy.WriteLine();
                Log_Energy.Write("{0}\t{1}\t{2}",
                    (phystime + dt).ToString("0.000", NumberFormatInfo.InvariantInfo),
                    EnergyIntegral.ToString("0.0000000000E+00", NumberFormatInfo.InvariantInfo),
                    omega.ToString("0.000000E+00", NumberFormatInfo.InvariantInfo));
                Log_Energy.Flush();
                if (TimestepNo == base.Control.NoOfTimesteps) {
                    Log_Energy.Close();
                }
            }
        }

        #endregion

        #region TaylorVortex

        /// <summary>
        /// Initialize TaylorVortex with time steps from analytical solution.
        /// </summary>
        void TaylorVortexHack() {
            WorkingSet.Velocity.IncreaseHistoryLength(SolverConf.Control.TimeOrder);

            // time step calculation
            double dt = base.Control.Endtime / base.Control.NoOfTimesteps;

            // init velocity
            for (int j = SolverConf.Control.TimeOrder - 1; j >= 0; j--) {
                WorkingSet.Velocity.Current[0].ProjectField(
                    ((_2D)((x, y) =>
                        (-Math.Cos(Math.PI * x) * Math.Sin(Math.PI * y) * Math.Exp(-2.0 * Math.PI * Math.PI * (-(double)j * dt) / 100.0)))).Vectorize());
                WorkingSet.Velocity.Current[1].ProjectField(
                    ((_2D)((x, y) =>
                        (Math.Sin(Math.PI * x) * Math.Cos(Math.PI * y) * Math.Exp(-2.0 * Math.PI * Math.PI * (-(double)j * dt) / 100.0)))).Vectorize());

                // push start value(s) to history
                WorkingSet.Velocity.Push();
            }

            // init pressure
            WorkingSet.Pressure.ProjectField(
                    ((_2D)((x, y) => (-0.25 * (Math.Cos(2.0 * Math.PI * x) + Math.Cos(2.0 * Math.PI * y))))).Vectorize());
        }

        #endregion

        #region DragAndLift

        Force DragAndLift;
        TextWriter Log_DragAndLift;

        /// <summary>
        /// Drag and lift for square cylinder.
        /// </summary>
        void CalculateDragAndLift(double phystime) {

            if (DragAndLift == null) {
                DragAndLift = new Force(this, SolverConf.Control.EdgeTagsDragAndLift);

                if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                    Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("DragAndLift", CurrentSessionInfo.ID);
                    string DragAndLiftHeader = "PhysTime\t xForce\t yForce\t zForce \t xMoment \t yMoment \t zMoment";
                    Console.WriteLine(DragAndLiftHeader);
                    Log_DragAndLift.WriteLine(DragAndLiftHeader);
                }
            }

            DragAndLift.CalculateForce();

            if ((base.MPIRank == 0) && (Log_DragAndLift != null)) {
                string DragAndLiftData = (phystime + "\t" + DragAndLift.XForce + "\t" + DragAndLift.YForce + "\t" + DragAndLift.ZForce + "\t" + DragAndLift.XMoment + "\t" + DragAndLift.YMoment + "\t" + DragAndLift.ZMoment);
                Console.WriteLine(DragAndLiftData);
                Log_DragAndLift.WriteLine(DragAndLiftData);

                Log_DragAndLift.Flush();
            }

            //var fs = File.Open("C:\\tmp\\DragAndLift.txt", FileMode.Append, FileAccess.Write);
            //StreamWriter sw = new StreamWriter(fs);
            //sw.WriteLine(phystime + "\t" + DragAndLift.XForce + "\t" + DragAndLift.YForce);
            //sw.Close();
            //fs.Close();
        }

        #endregion

        #region NusseltNumber

        NusseltNumber NusseltNum;
        TextWriter LogNusselt;

        /// <summary>
        /// Calculate Nusselt number for Low-Mach number flows.
        /// </summary>
        /// <param name="Timestep"></param>
        /// <param name="GridDat"></param>
        /// <param name="Temperature"></param>
        /// <param name="SolverConf"></param>
        void CalculateNusselt(TimestepNumber Timestep, IGridData GridDat, SinglePhaseField Temperature, SIMPLEControl SolverConf) {

            // Initialize calculation of Nusselt number.
            if (NusseltNum == null) {
                LowMachSIMPLEControl lowMachConf = SolverConf as LowMachSIMPLEControl;
                NusseltNum = new NusseltNumber(GridDat, Temperature, lowMachConf.EoS, lowMachConf.EdgeTagsNusselt);

                if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                    LogNusselt = base.DatabaseDriver.FsDriver.GetNewLog("Nusselt", CurrentSessionInfo.ID);
                    LogNusselt.Write("Timestep");
                    for (int bc = 0; bc < NusseltNum.Nusselt.Length; bc++)
                        LogNusselt.Write("\t" + lowMachConf.EdgeTagsNusselt[bc]);
                    LogNusselt.WriteLine();
                }
            }

            // Calculate Nusselt number
            NusseltNum.CalculateNusseltNumber();

            // Write result to text file
            if ((base.MPIRank == 0) && (LogNusselt != null)) {
                LogNusselt.Write(Timestep.ToString());
                for (int bc = 0; bc < NusseltNum.Nusselt.Length; bc++)
                    LogNusselt.Write("\t" + NusseltNum.Nusselt[bc].ToString("0.0000000000E+00", NumberFormatInfo.InvariantInfo));
                LogNusselt.WriteLine();
                LogNusselt.Flush();
            }
        }

        #endregion NusseltNumber

        #region TemperatureHistoryTallCavity

        TextWriter LogTemperatureHistory;

        /// <summary>
        /// Log temperature history for tall cavity.
        /// </summary>
        /// <param name="time"></param>        
        /// <param name="Temperature"></param>        
        void LogTemperature(double time, SinglePhaseField Temperature) {

            // Initialize
            if (LogTemperatureHistory == null) {
                if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                    LogTemperatureHistory = base.DatabaseDriver.FsDriver.GetNewLog("TemperatureHistory", CurrentSessionInfo.ID);
                    LogTemperatureHistory.Write("Time");
                    LogTemperatureHistory.Write("\t" + "Temperature");
                    LogTemperatureHistory.WriteLine();
                }
            }

            double TemperatureValue = Temperature.ProbeAt(new double[] { 0.181, 7.37 });

            // Write result to text file
            if (base.MPIRank == 0) {
                LogTemperatureHistory.Write(time);
                LogTemperatureHistory.Write("\t" + TemperatureValue.ToString("0.0000000000E+00", NumberFormatInfo.InvariantInfo));
                LogTemperatureHistory.WriteLine();
                LogTemperatureHistory.Flush();
            }
        }

        #endregion TemperatureHistoryTallCavity
    }
}
