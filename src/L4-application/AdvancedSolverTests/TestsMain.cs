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

    public class TestsMain {



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
            //AdvancedSolverTests.Solver.mklILU.CompareFactorization();
            //AdvancedSolverTests.SolverChooser.ConfigTest.TestLinearSolverConfigurations();
            //AdvancedSolverTests.Solver.mklILUtest.CompareFactorization();
            //AdvancedSolverTests.SubBlocking.LocalTests.SubSelection(SelectionType.all_combined);
            //Debugger.Launch();
            AdvancedSolverTests.SolverChooser.ConfigTest.TestNonLinearSolverConfigurations();
        }

    }
}
