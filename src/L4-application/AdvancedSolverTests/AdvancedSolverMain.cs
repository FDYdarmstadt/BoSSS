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

    public class AdvancedSolverMain {



        public static void Main() {
            BoSSS.Solution.Application.InitMPI();
            Test();
            BoSSS.Solution.Application.FinalizeMPI();
        }

        public static void Test() {
            //Console.WriteLine("wer hat den output eingestellt: " + ilPSP.Environment.StdoutOnlyOnRank0);
            //AdvancedSolverTests.SubBlocking.LocalTests.MapConsistencyTest(XDGusage.none, 2);
            //AdvancedSolverTests.SubBlocking.LocalTests.SubMatrixExtractionWithCoupling(XDGusage.none, 2, MatrixShape.diagonal);
            //AdvancedSolverTests.SubBlocking.LocalTests.SubMatrixExtractionWithCoupling(XDGusage.all, 2, MatrixShape.diagonal);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubBlockExtraction(XDGusage.none, 2, MatrixShape.diagonal_var_spec, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var_spec, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var, 4);
            //AdvancedSolverTests.SubBlocking.ExternalTests.GetExternalRowsTest(XDGusage.all, 2, 4);
            AdvancedSolverTests.Solver.addSchwarzTest.RunTest();
        }

    }
}
