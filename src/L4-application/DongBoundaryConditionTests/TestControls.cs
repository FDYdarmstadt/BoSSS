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
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.LevelSetTools;
using System.Runtime.Serialization;


namespace BoSSS.Application.XNSE_Solver.DongBoundaryConditionTests {

    [DataContract]
    public class XNSE_Test_Control : XNSE_Control {

        public override Type GetSolverType() {
            return typeof(XNSE_TestSolver<XNSE_Test_Control>);
        }

        [DataMember]
        public bool UseManufacturedComps = false;
    }


    public static class DongBoundaryConditionTests {

        public static XNSE_Test_Control KovasznayFlow_SteadyState(int k = 3, int Res = 4) {

            XNSE_Test_Control C = new XNSE_Test_Control();

            C.SetDGdegree(k);

            // physical parameters
            double rhoSpc = 1.0;
            C.PhysicalParameters.rho_A = rhoSpc;
            double muSpc = 1.0 / 40.0;
            C.PhysicalParameters.mu_A = muSpc;

            C.PhysicalParameters.IncludeConvection = true;

            C.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(-0.5, 0.1, (1 * Res) + 1);
                double[] yNodes = GenericBlas.Linspace(-0.5, 0.5, (2 * Res) + 1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicY: true);

                grd.DefineEdgeTags(delegate (double[] X) {
                    string ret = null;
                    if ((X[0] + 0.5).Abs() <= 1e-8)
                        ret = IncompressibleBcType.Velocity_Inlet.ToString();
                    if ((X[0] - 0.1).Abs() <= 1e-8)
                        ret = IncompressibleBcType.Velocity_Inlet.ToString();
                    //if ((X[1] + 0.5).Abs() <= 1e-8)
                    //    ret = IncompressibleBcType.Velocity_Inlet.ToString();
                    //if ((X[1] - 1.5).Abs() <= 1e-8)
                    //    ret = IncompressibleBcType.Velocity_Inlet.ToString();
                    return ret;
                });

                return grd;

            };

            // initial conditions
            C.AddInitialValue("VelocityX#A", KovasznayFlowSolutions.KovasznayFlow_u);
            C.AddInitialValue("VelocityY#A", KovasznayFlowSolutions.KovasznayFlow_v);
            C.AddInitialValue("Pressure#A", KovasznayFlowSolutions.KovasznayFlow_p);

            // boundary conditions
            C.AddBoundaryValue(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityX#A", KovasznayFlowSolutions.KovasznayFlow_u);
            C.AddBoundaryValue(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityY#A", KovasznayFlowSolutions.KovasznayFlow_v);

            //C.AddBoundaryValue(IncompressibleBcType.Dong_OutFlow.ToString(), "VelocityX#A", KovasznayFlow_u);
            //C.AddBoundaryValue(IncompressibleBcType.Dong_OutFlow.ToString(), "VelocityY#A", KovasznayFlow_v);

            C.SkipSolveAndEvaluateResidual = true;

            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.Globalization = Newton.GlobalizationOption.LineSearch;
            C.NonLinearSolver.MaxSolverIterations = 50;

            // C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            return C;

        }

    }

}
