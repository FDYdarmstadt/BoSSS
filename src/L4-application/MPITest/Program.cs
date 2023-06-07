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
using System.Globalization;
using System.Linq;
using System.Threading;
using ilPSP;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;

namespace MPITest {

    /// <summary>
    /// Testing if MPI and some other base routines are working
    /// </summary>
    [TestFixture]
    public class Program {

        public static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            Test();
            MPI.Wrappers.csMPI.Raw.mpiFinalize(); 
        }

   
        [Test]
        public static void Test() {
            TestAllreduce();
            TestLapack();
        }

        private static void TestAllreduce() {
            int rank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);

            int size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

            ilPSP.Environment.StdoutOnlyOnRank0 = false;
            Console.WriteLine("Hello from " + rank + " of " + size + ".");



            int result = 0;
            unsafe {
                csMPI.Raw.Allreduce(
                    (IntPtr)(&rank),
                    (IntPtr)(&result),
                    1,
                    csMPI.Raw._DATATYPE.INT,
                    csMPI.Raw._OP.SUM,
                    csMPI.Raw._COMM.WORLD);
            }

            int expectedResult = Enumerable.Range(0, size).Sum();
            Assert.IsTrue(result == expectedResult);
        }

        public static void TestLapack() {


            Random rnd = new Random();
            for(int n = 2; n < 50; n++) {

                var SomeMtx = MultidimensionalArray.Create(n, n);
                for(int i = 0; i < n; i++) {
                    for(int j = 0; j < n; j++) {
                        double r = rnd.NextDouble();
                        SomeMtx[i, j] += r;
                        SomeMtx[j, i] += r;
                    }
                }


                var eigs = SomeMtx.EigsSymm();

                var (EigVal, EigVec) = SomeMtx.EigenspaceSymm();
                
                var D = MultidimensionalArray.CreateDiagMtx(EigVal);

                var Test = EigVec.GEMM(D, EigVec.TransposeTo());
                Test.Acc(-1.0, SomeMtx);
                double err = Test.InfNorm();
                Assert.LessOrEqual(err, 1e-8);
            }

            

        }
    }
}
