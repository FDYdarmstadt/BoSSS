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
using System.Diagnostics;
using NUnit.Framework;

namespace ilPSP {

    /// <summary>
    /// Various tests for <see cref="MultidimensionalArray"/>, especially for the 
    /// arithmetical functions.
    /// </summary>
    [TestFixture]
    public class MultidimensionalArray_Tests {

        /// <summary>
        /// for direct execution...
        /// </summary>
        public static void Main() {
            //MultiplyTest1();
            for (int i = 0; i < 3; i++) { // due to JIT-compilation, we don't get realistic runtimes at the first time
                Console.WriteLine("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
                Console.WriteLine("RUN #" + i);
                Console.WriteLine("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
                
                MultiplyTest0();
                MultiplyTest1();
                MultiplyTest2();
                MultiplyTest3();

                MultiplyTrafoTest0();

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("");
            }

        }



        /// <summary>
        /// test for <see cref="MultidimensionalArray.Multiply(double,MultidimensionalArray,MultidimensionalArray,double,string,string,string)"/>
        /// </summary>
        [Test]
        public static void MultiplyTest0() { // one summation index 
            foreach(int C in new int[] { 2, 60 }) { // test unrolling and standard path

                int I = 120;
                int J = 21;
                int K = 40;


                MultidimensionalArray A = MultidimensionalArray.Create(I, C, K);
                MultidimensionalArray B = MultidimensionalArray.Create(K, C, J);

                Console.WriteLine("number of operands in A: " + A.Length);
                Console.WriteLine("number of operands in B: " + B.Length);

                MultidimensionalArray ResTst1 = MultidimensionalArray.Create(J, K, I);
                MultidimensionalArray ResTst2 = MultidimensionalArray.Create(J, K, I);
                MultidimensionalArray ResChck = MultidimensionalArray.Create(J, K, I);

                // fill operands with random values
                Random rnd = new Random();
                A.ApplyAll(x => rnd.NextDouble());
                B.ApplyAll(x => rnd.NextDouble());

                double alpha = 0.99;
                double beta = 0;

                Stopwatch sw = new Stopwatch();


                // tensorized multiplication:
                sw.Reset();
                sw.Start();
                ResTst1.Multiply(alpha, A, B, beta, "jki", "ick", "kcj");
                ResTst2.Multiply(alpha, B, A, beta, "jki", "kcj", "ick");
                sw.Stop();

                Console.WriteLine("runtime of tensorized multiplication: " + sw.ElapsedMilliseconds + " millisec.");

                // old-fashioned equivalent with loops:
                double errsum = 0;
                sw.Reset();
                sw.Start();
                for(int i = 0; i < I; i++) {
                    for(int j = 0; j < J; j++) {
                        for(int k = 0; k < K; k++) {

                            // summation:
                            double sum = 0;
                            for(int c = 0; c < C; c++)
                                sum += A[i, c, k] * B[k, c, j];

                            ResChck[j, k, i] = sum * alpha + ResChck[j, k, i] * beta;

                            errsum += Math.Abs(ResTst1[j, k, i] - ResChck[j, k, i]);
                            errsum += Math.Abs(ResTst2[j, k, i] - ResChck[j, k, i]);
                        }
                    }
                }
                sw.Stop();

                Console.WriteLine("runtime of loop multiplication: " + sw.ElapsedMilliseconds + " millisec.");
                Console.WriteLine("total error: " + errsum);

                double thres = 1.0e-13;
                Assert.IsTrue(errsum < thres);
            }
        }


        /// <summary>
        /// test for <see cref="MultidimensionalArray.Multiply(double,MultidimensionalArray,MultidimensionalArray,double,string,string,string)"/>
        /// </summary>
        [Test]
        public static void MultiplyTest1() { // two summation indices (k,c) and (k) occurs TWICE in second array like in matrix trace operation.
            int I = 120;
            int J = 21;
            int K = 40;
            int C = 60;

            MultidimensionalArray A = MultidimensionalArray.Create(I, C, K);
            MultidimensionalArray B = MultidimensionalArray.Create(K, C, J, K);

            Console.WriteLine("number of operands in A: " + A.Length);
            Console.WriteLine("number of operands in B: " + B.Length);

            MultidimensionalArray ResTst1 = MultidimensionalArray.Create(J, I);
            MultidimensionalArray ResTst2 = MultidimensionalArray.Create(J, I);
            MultidimensionalArray ResChck = MultidimensionalArray.Create(J, I);

            // fill operands with random values
            Random rnd = new Random();
            A.ApplyAll(x => rnd.NextDouble());
            B.ApplyAll(x => rnd.NextDouble());

            double alpha = 0.99;
            double beta = 0;

            Stopwatch sw = new Stopwatch();


            // tensorized multiplication:
            sw.Reset();
            sw.Start();
            ResTst1.Multiply(alpha, A, B, beta, "ji", "ick", "kcjk");
            sw.Stop();
            ResTst2.Multiply(alpha, B, A, beta, "ji", "kcjk", "ick");
            

            Console.WriteLine("runtime of tensorized multiplication: " + sw.ElapsedMilliseconds + " millisec.");

            // old-fashioned equivalent with loops:
            double errsum = 0;
            sw.Reset();
            sw.Start();
            for (int i = 0; i < I; i++) {
                for (int j = 0; j < J; j++) {
                    
                    // summation:
                    double sum = 0;
                    for (int c = 0; c < C; c++)
                        for (int k = 0; k < K; k++)
                            sum += A[i, c, k] * B[k, c, j, k];

                    ResChck[j, i] = sum * alpha + ResChck[j, i] * beta;

                    //errsum += Math.Abs(ResTst1[j, i] - ResChck[j, i]);
                    errsum += Math.Abs(ResTst2[j, i] - ResChck[j, i]);
                }
            }
            sw.Stop();

            Console.WriteLine("runtime of loop multiplication: " + sw.ElapsedMilliseconds + " millisec.");
            Console.WriteLine("total error: " + errsum);

            double thres = 1.0e-13;
            Assert.IsTrue(errsum < thres);
        
        }

        /// <summary>
        /// test for <see cref="MultidimensionalArray.Multiply(double,MultidimensionalArray,MultidimensionalArray,double,string,string,string)"/>
        /// </summary>
        [Test]
        public static void MultiplyTest2() { // no summation, only tenzorization
            int I = 125;
            int K = 21;
            int M = 43;
            int N = 63;

            MultidimensionalArray A = MultidimensionalArray.Create(I, K, M);
            MultidimensionalArray B = MultidimensionalArray.Create(I, K, N);

            Console.WriteLine("number of operands in A: " + A.Length);
            Console.WriteLine("number of operands in B: " + B.Length);

            MultidimensionalArray ResTst1 = MultidimensionalArray.Create(I, K, M, N);
            MultidimensionalArray ResTst2 = MultidimensionalArray.Create(I, K, M, N);
            MultidimensionalArray ResChck = MultidimensionalArray.Create(I, K, M, N);

            // fill operands with random values
            Random rnd = new Random();
            A.ApplyAll(x => rnd.NextDouble());
            B.ApplyAll(x => rnd.NextDouble());


            // tensorized multiplication:
            Stopwatch TenMult = new Stopwatch();
            TenMult.Start();
            ResTst1.Multiply(1.0, A, B, 0.0, "ikmn", "ikm", "ikn");
            TenMult.Stop();
            ResTst2.Multiply(1.0, B, A, 0.0, "ikmn", "ikn", "ikm");
            Console.WriteLine("runtime of tensorized multiplication: " + TenMult.ElapsedMilliseconds + " millisec.");

            // comparison code
            Stopwatch RefMult = new Stopwatch();
            RefMult.Start();
            
            double errSum = 0;
            for (int i = 0; i < I; i++) {
                for (int k = 0; k < K; k++) {
                    for (int n = 0; n < N; n++) {
                        for (int m = 0; m < M; m++) {
                            ResChck[i, k, m, n] = A[i, k, m] * B[i, k, n];
                            errSum += Math.Abs(ResChck[i, k, m, n] - ResTst1[i, k, m, n]);
                            errSum += Math.Abs(ResChck[i, k, m, n] - ResTst2[i, k, m, n]);
                        }
                    }

                }

            }
            RefMult.Stop();
            Console.WriteLine("runtime of loop multiplication: " + RefMult.ElapsedMilliseconds + " millisec.");

            Console.WriteLine("total error: " + errSum);

            double thres = 1.0e-13;
            Assert.IsTrue(errSum < thres);
        }


        /// <summary>
        /// test for <see cref="MultidimensionalArray.Multiply(double,MultidimensionalArray,MultidimensionalArray,double,string,string,string)"/>
        /// </summary>
        [Test]
        public static void MultiplyTest3() { // two summation indices (k,r) 
            foreach(int K in new int[] { 2, 21 }) { // test unrolling and standard path
                int I = 12;
                int M = 43;
                int N = 63;
                int R = 21;

                MultidimensionalArray A = MultidimensionalArray.Create(I, R, K, M);
                MultidimensionalArray B = MultidimensionalArray.Create(I, K, N, R);

                Console.WriteLine("number of operands in A: " + A.Length);
                Console.WriteLine("number of operands in B: " + B.Length);

                MultidimensionalArray ResTst1 = MultidimensionalArray.Create(I, M, N);
                MultidimensionalArray ResTst2 = MultidimensionalArray.Create(I, M, N);
                MultidimensionalArray ResChck = MultidimensionalArray.Create(I, M, N);

                // fill operands with random values
                Random rnd = new Random();
                A.ApplyAll(x => rnd.NextDouble());
                B.ApplyAll(x => rnd.NextDouble());
                ResTst1.ApplyAll(x => rnd.NextDouble());
                ResChck.Set(ResTst1);
                ResTst2.Set(ResTst1);


                double alpha = 0.67;
                double beta = 1.3;


                // tensorized multiplication:
                Stopwatch TenMult = new Stopwatch();
                TenMult.Start();
                ResTst1.Multiply(alpha, A, B, beta, "imn", "irkm", "iknr");
                TenMult.Stop();
                ResTst2.Multiply(alpha, B, A, beta, "imn", "iknr", "irkm");
                Console.WriteLine("runtime of tensorized multiplication: " + TenMult.ElapsedMilliseconds + " millisec.");

                // comparison code
                Stopwatch RefMult = new Stopwatch();
                RefMult.Start();

                double errSum = 0;
                for(int i = 0; i < I; i++) {
                    for(int n = 0; n < N; n++) {
                        for(int m = 0; m < M; m++) {

                            // summation:
                            double sum = 0;
                            for(int r = 0; r < R; r++)
                                for(int k = 0; k < K; k++)
                                    sum += A[i, r, k, m] * B[i, k, n, r];

                            ResChck[i, m, n] = sum * alpha + ResChck[i, m, n] * beta;

                            errSum += Math.Abs(ResTst1[i, m, n] - ResChck[i, m, n]);
                            errSum += Math.Abs(ResTst2[i, m, n] - ResChck[i, m, n]);

                        }

                    }

                }
                RefMult.Stop();
                Console.WriteLine("runtime of loop multiplication: " + RefMult.ElapsedMilliseconds + " millisec.");

                Console.WriteLine("total error: " + errSum);

                double thres = 1.0e-6;
                Assert.IsTrue(errSum < thres);
            }
        }

        /// <summary>
        /// test for <see cref="MultidimensionalArray.Multiply(double,MultidimensionalArray,MultidimensionalArray,double,string,string,string)"/>
        /// </summary>
        [Test]
        public static void MultiplyTrafoTest0() { // two summation indices (k,r) 
            foreach(int K in new int[] { 210 }) { // test unrolling and standard path
                int I = 120;
                int M = 43;

                MultidimensionalArray A = MultidimensionalArray.Create(I, K, 2 * M);
                MultidimensionalArray B = MultidimensionalArray.Create(2 * M, K); ;

                Console.WriteLine("number of operands in A: " + A.Length);
                Console.WriteLine("number of operands in B: " + B.Length);

                int[] mTrafo = new int[M];

                MultidimensionalArray ResTst1 = MultidimensionalArray.Create(I, M);
                MultidimensionalArray ResTst2 = MultidimensionalArray.Create(I, M);
                MultidimensionalArray ResChck = MultidimensionalArray.Create(I, M);

                // fill operands with random values
                Random rnd = new Random();
                A.ApplyAll(x => rnd.NextDouble());
                B.ApplyAll(x => rnd.NextDouble());

                ResTst1.ApplyAll(x => rnd.NextDouble());
                ResChck.Set(ResTst1);
                ResTst2.Set(ResTst1);

                for(int m = 0; m < M; m++) {
                    mTrafo[m] = rnd.Next(2 * M);
                    //mTrafo[m] = m;
                    Debug.Assert(mTrafo[m] < 2 * M);
                }


                double alpha = 0.67;
                double beta = 1.3;

                var mp1 = MultidimensionalArray.MultiplyProgram.Compile("im", "ikT(m)", "T(m)k", true);
                var mp2 = MultidimensionalArray.MultiplyProgram.Compile("im", "T(m)k", "ikT(m)", true);


                // tensorized multiplication:
                Stopwatch TenMult = new Stopwatch();
                TenMult.Start();
                ResTst1.Multiply(alpha, A, B, beta, ref mp1, mTrafo);
                TenMult.Stop();
                ResTst2.Multiply(alpha, B, A, beta, ref mp2, mTrafo);
                Console.WriteLine("runtime of tensorized multiplication: " + TenMult.ElapsedMilliseconds + " millisec.");

                // comparison code
                Stopwatch RefMult = new Stopwatch();
                RefMult.Start();

                double errSum = 0;
                for(int i = 0; i < I; i++) {

                    for(int m = 0; m < M; m++) {

                        int m_trf = mTrafo[m];

                        // summation:
                        double sum = 0;
                        for(int k = 0; k < K; k++)
                            sum += A[i, k, m_trf] * B[m_trf, k];

                        ResChck[i, m] = sum * alpha + ResChck[i, m] * beta;

                        errSum += Math.Abs(ResTst1[i, m] - ResChck[i, m]);
                        errSum += Math.Abs(ResTst2[i, m] - ResChck[i, m]);

                    }

                }

                
                RefMult.Stop();
                Console.WriteLine("runtime of loop multiplication: " + RefMult.ElapsedMilliseconds + " millisec.");

                Console.WriteLine("total error: " + errSum);

                double thres = 1.0e-6;
                Assert.IsTrue(errSum < thres);
            }
        }
    }
}
