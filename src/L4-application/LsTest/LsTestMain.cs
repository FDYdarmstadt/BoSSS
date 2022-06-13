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

using BoSSS.Solution.LevelSetTools;

namespace BoSSS.Application.LsTest {
    partial class SolverWithLevelSetUpdaterTestCenter {

        static void Main(string[] args) {

            string ArgName = "BOSSS_ARG_0";
            string ArgValue = System.Environment.GetEnvironmentVariable(ArgName);
            if (ArgValue == null && args.Length == 0) {
                InitMPI();
                //LevelSetUnitTests.LevelSetCubeProjectionConvergenceTest(2, LevelSetEvolution.StokesExtension, BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting, false);
                //LevelSetUnitTests.LevelSetCircleProjectionConvergenceTest(2, LevelSetEvolution.StokesExtension, BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting, false, true);
                //LevelSetUnitTests.LevelSetZalasakDiscConvergenceTest(3, LevelSetEvolution.StokesExtension, BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting, false);
                //LevelSetUnitTests.LevelSetZalesakDisc(2, 0, LevelSetEvolution.StokesExtension, BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting, false);
                LevelSetUnitTests.LevelSetSwirlingFlow(2, 0, LevelSetEvolution.PrescribedVelocity, BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting, false);

                FinalizeMPI();
            }

            {
                SolverWithLevelSetUpdaterTestCenter._Main(args, false, delegate () {
                    var p = new SolverWithLevelSetUpdaterTestCenter();
                    return p;
                });
            }
        }

    }
}
