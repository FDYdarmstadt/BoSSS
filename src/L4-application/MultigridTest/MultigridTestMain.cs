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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Linq;
using System.Runtime.Serialization;
using BoSSS.Platform;
using System.Collections;
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.Statistic;

namespace BoSSS.Application.MultigridTest {

    /// <summary>
    /// Les main class
    /// </summary>
    public static class MultigridMain {

        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            TestProgramCurvedMG.Init();
            TestProgramLinearMG.Init();

            BoSSS.Application.MultigridTest.TestProgramCurvedMG.RestrictionOfSystemOpTest();
            //BoSSS.Application.MultigridTest.TestProgramLinearMG.RestrictionOfSystemOpTest();
            //BoSSS.Application.MultigridTest.TestProgramLinearMG.XDG_Norm(0, 0.1, 1, MultigridOperator.Mode.SymPart_DiagBlockEquilib);
            //BoSSS.Application.MultigridTest.TestProgramLinearMG.XDG_MatrixPolynomialRestAndPrlgTest_2(0, 0.3d);
            //BoSSS.Application.MultigridTest.TestProgram.XDG_MatrixPolynomialRestAndPrlgTest(1, 0.0d, 1);
            //BoSSS.Application.MultigridTest.TestProgramLinearMG.XDG_ProlongationTest(3, 0.0d, 1, MultigridOperator.Mode.IdMass);
            //TestProgram.XDG_ProlongationTest(0, 0.3, 1, MultigridOperator.Mode.IdMass);
            //BoSSS.Application.MultigridTest.TestProgram.XDG_MatrixPolynomialRestAndPrlgTest(0, 0.0d, 0);
            //BoSSS.Application.MultigridTest.TestProgram.XDG_MatrixPolynomialRestAndPrlgTest(2, 0.0d, 0);
            //BoSSS.Application.MultigridTest.TestProgram.XDG_MatrixPolynomialRestAndPrlgTest(3, 0.0d, 1);
            //BoSSS.Application.MultigridTest.TestProgramCurvedMG.XDG_MatrixPolynomialRestAndPrlgTest_2(0, 0.0);
            //BoSSS.Application.MultigridTest.TestProgramLinearMG.XDG_MatrixPolynomialRestAndPrlgTest(0, 0.0d, 0);
            foreach (int w in new int[] { 0 }) {
                for (int p = 1; p <= 1; p++) {
                    //TestProgram.ProlongationTest(p);
                    //TestProgram.PolynomialRestAndPrlgTest(p);
                    //XDG_MatrixPolynomialRestAndPrlgTest_2
                    //TestProgram.RestictionMatrixTest(p);
                    //TestProgram.XDG_PolynomialRestAndPrlgTest(p, 0.3, w);
                    //TestProgram.XDG_MatrixPolynomialRestAndPrlgTest_2(p, 0.3);
                    //TestProgram.XDG_MatrixPolynomialRestAndPrlgTest_2(p, 0.3, w, MultigridOperator.Mode.IdMass);
                    //TestProgram.XDG_ProlongationTest(0, 0.0, w, MultigridOperator.Mode.IdMass);
                }
            }


            BoSSS.Solution.Application.InitMPI();
        }





        

  
    }
}

