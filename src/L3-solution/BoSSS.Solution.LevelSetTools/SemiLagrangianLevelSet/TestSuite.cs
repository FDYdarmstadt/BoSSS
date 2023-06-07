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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.IO;

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.Timestepping;
using System.Diagnostics;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTool.SemiLagrangianLevelSet;

namespace BoSSS.Application.SemiLagrangianLevelSetTestSuite
{
    /// <summary>
    /// Test application for level-set movement.
    /// </summary>
    public class TestSuite : Solution.Application
    {
        static void Main(string[] args)
        {

            _Main(args, true,
                delegate () {
                    // dbg_launch();
                    return new TestSuite();
                });
        }

        //private TestCase Case = new MovingCircle(LagrangianMode.Marker, 0.1, 100);
        private TestCase Case = new TheKartoffel(LagrangianMode.Marker);

        /// <summary>
        /// Level set function
        /// </summary>

        /// <summary>
        /// A level set tracker for different analysis such as calculating the areas/volumes of the negative (denoted by NLS) and positive (denoted by PLS) parts of the level set function
        /// </summary>

        /// <summary>
        /// velocity
        /// </summary>

        /// <summary>
        /// Initialize Fields based on <see cref="LevelSetControl"/>
        /// </summary>
        protected override void CreateFields()
        {
            Case.CreateFields(this.GridData);
        }

        protected override IGrid CreateOrLoadGrid()
        {
            return Case.CreateGrid();
        }

        /// <summary>
        /// Initialize <see cref="MovementAlgorithm"/> and ReInit
        /// </summary>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L)
        {
        }

        /// <summary>
        /// Performing one timestep time marching
        /// </summary>
        /// <param name="TimestepNo">Major Timestep number for physical Timesteps</param>
        /// <param name="phystime"> Physical time </param>
        /// <param name="dt"> TimestepSize </param>
        /// <returns></returns>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt)
        {
            dt = Case.dt;
            Console.WriteLine("Performing Timestep #{0} - Endtime: {1}, dt {2} ", TimestepNo, phystime + dt, dt);
            Case.TimestepWork(phystime, dt,TimestepNo);
            PlotCurrentState(phystime, TimestepNo, 3);

            phystime += dt;
            Console.WriteLine();
            return dt;
        }

        /// <summary>
        /// Set Intial Values and perform first SpecFemSmoothing, if set in Control File
        /// </summary>
        protected override void SetInitial(double t)
        {
            base.NoOfTimesteps = Case.NbrofTimesteps;
            Case.InitialWork();
            PlotCurrentState(0.0, 0, 3);
        }

        /// <summary>
        /// Write Out
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="timestepNo"></param>
        /// <param name="superSampling"></param>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling)
        {
            Case.Plot(timestepNo, physTime);
        }
    }
}
