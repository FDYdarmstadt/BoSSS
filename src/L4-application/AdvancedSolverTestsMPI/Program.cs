using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using BoSSS.Solution;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;
using AdvancedSolverTests.SubBlocking;
using System.Diagnostics;

namespace AdvancedSolverTests {

    public class MPITests {



        public static void Main() {
            BoSSS.Solution.Application.InitMPI();
            //Debugger.Launch();
            Test();
            BoSSS.Solution.Application.FinalizeMPI();
        }

        public static void Test() {
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubBlockExtraction(XDGusage.none, 2, MatrixShape.diagonal_var_spec, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var_spec, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.GetExternalRowsTest(XDGusage.all, 2, 4);
            AdvancedSolverTestsMPI.Solver.SchwarzForCoarseMeshTest.TestInit(XDGusage.all, 2, 4, -1, 1);
            //AdvancedSolverTests.Solver.addSchwarzTest.RunTest();
            //AdvancedSolverTests.SubBlocking.ExternalTests.GetExternalRowsTest(XDGusage.none, 2, 4);

        }

    }
}