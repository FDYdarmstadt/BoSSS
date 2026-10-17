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
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP;


namespace BoSSS.Application.XNSE_Solver.DongBoundaryConditionTests {


    public class XNSE_DongTest : XNSE_TestSolver<XNSE_Test_Control> {

        // ===========
        //  Main file
        // ===========
        static void Main(string[] args) {
            XNSE_DongTest._Main(args, false, delegate () {
                var p = new XNSE_DongTest();
                return p;
            });
        }


        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            dt = base.RunSolverOneStep(TimestepNo, phystime, dt);

            double l2normVelX = this.Velocity[0].L2Error((ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {
                KovasznayFlowSolutions.KovasznayFlow_u.EvaluateV(nodes, 0.0, results);
            });
            Console.WriteLine($"l2 error - VelocityX = {l2normVelX}");

            double l2normVelY = this.Velocity[1].L2Error((ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {
                KovasznayFlowSolutions.KovasznayFlow_v.EvaluateV(nodes, 0.0, results);
            });
            Console.WriteLine($"l2 error - VelocityY = {l2normVelY}");

            //double[] PrefPoint = new double[] { -0.5, 0.0 };
            //double Prefvalue = this.Pressure.ProbeAt(PrefPoint);
            //double PsolValue = KovasznayFlowSolutions.KovasznayFlow_p.Evaluate(PrefPoint, 0.0);
            //double Pdiff = Prefvalue - PsolValue;
            //Console.WriteLine($"PrefValue = {Prefvalue}; PsolValue = {PsolValue}; Pdiff = {Pdiff}");
            //this.Pressure.AccConstant(Pdiff);
            double l2normP = this.Pressure.L2Error((ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {
                KovasznayFlowSolutions.KovasznayFlow_p.EvaluateV(nodes, 0.0, results);
            });
            Console.WriteLine($"l2 error - Pressure = {l2normP}");

            return dt;
        }
    }


    public class XNSE_TestSolver<T> : XNSE<T> where T : XNSE_Test_Control, new() {

        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystem(D, opFactory, lsUpdater);

            if (this.Control.UseManufacturedComps) {
                if (D != 2)
                    throw new NotSupportedException("Test setup only for 2D case");

                XNSE_OperatorConfiguration config = new XNSE_OperatorConfiguration(this.Control);

                Console.WriteLine("add manufactured solution components for Kovasznay flow with Dong B.C.");
                for (int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new ManufacturedSolutionComp(d, boundaryMap, config, "A"));
                    opFactory.AddEquation(new ManufacturedSolutionComp(d, boundaryMap, config, "B"));
                }
            }
        }

    }


}
