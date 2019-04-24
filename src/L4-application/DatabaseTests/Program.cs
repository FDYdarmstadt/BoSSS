using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.Diagnostics;

namespace BoSSS.Application.DatabaseTests
{
    class Program
    {
        public static void Main()
        {
            MPITest.InitOnce();
            PrintActiveThreads();
#if DEBUG
            StartDebugger();
#endif
            RunTest();
            MPITest.TearDown();
        }

        static void RunTest()
        {
            var tst = new DBDriverTests();
            tst.Init();
            tst.TestRenameGrid();
            tst.CleanUp();
        }

        static void PrintActiveThreads()
        {
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);

            ilPSP.Environment.StdoutOnlyOnRank0 = false;
            Console.WriteLine("Hello from " + rank + " of " + size + ".");
        }

        static void StartDebugger()
        {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            if(size > 1)
            {
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
                if (rank == 1)
                    Debugger.Launch();
            }
        }
    }
}
