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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using ilPSP;
using System.Diagnostics;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Tests whether the XNSE solver (<see cref="XNSE_SolverMain"/>) also works MPI-parallel for larger cases (8 MPI cores)
    /// </summary>
    [TestFixture]
    public static class XNSE_Solver_LargeMPItest {

        [Test]
        static public void ParallelRotatingSphere() {
            var C = XNSE_Solver_MPItest.Rotating_Sphere(k: 1, Res: 20, SpaceDim: 3, useAMR: true, useLoadBal: true);
            //C.TracingNamespaces = "*";

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
            //Assert.IsTrue(false, "this should fail!");
        }

        [Test]
        public static void ParallelRotatingTilted3DTorus() {
            var C = HardcodedControl.RotatingTiltedXRigid(k: 1, Res: 20, SpaceDim: 3, AMR: true, AMRLevel: 1, TiltAngle: Math.PI / 4, SolverOn: true);

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        public static void RotatingTilted3DTorusAgg0() {
            var C = HardcodedControl.RotatingTilted3DTorusAgg0();

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }
        
        /// <summary>
        /// Initiates all the test cases
        /// </summary>
        static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            ParallelRotatingSphere();
            ParallelRotatingTilted3DTorus();
            BoSSS.Solution.Application.FinalizeMPI();
        }




    }
}
