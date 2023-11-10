using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.Diagnostics;
using IntersectingQuadrature.TensorAnalysis;

namespace BoSSS.Application.DatabaseTests {
    public class DatabaseTestsProgram {
        public static void Main() {
            BoSSS.Solution.Application.InitMPI();
            PrintActiveThreads();
#if DEBUG
            StartDebugger();
#endif
            RunTest();
            BoSSS.Solution.Application.FinalizeMPI();
        }

        static void RunTest() {
            //var tst = new MiscTests();
            //tst.Init();
            //tst.GridEquivalenceTest();
            //tst.CleanUp();

            var t2 = new DBDriverTests();
            t2.Init();
            t2.TestCopySession();
            t2.CleanUp();
        }

        static void PrintActiveThreads() {
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);

            ilPSP.Environment.StdoutOnlyOnRank0 = false;
            Console.WriteLine("Hello from " + rank + " of " + size + ".");
        }

        static void StartDebugger() {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            if(size > 1) {
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
                if(rank == 1)
                     Debugger.Launch();
            }
        }
    }
}
