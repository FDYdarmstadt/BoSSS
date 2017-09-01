using System;
using ilPSP;
using ilPSP.Threading;
using System.Runtime.InteropServices;

namespace ThreadTest {
    unsafe class ThreadTest {


        [DllImport("Kernel32.dll")]
        static extern bool QueryPerformanceCounter(out long lpPerformanceCount);


        void Loop(int i0, int iE) {
            for (int i = i0; i < iE; i++) {
                double acc = 0;
                for (int j = 0; j < N; j++)
                    acc += Mtx[i, j] * A[j];
                B[i] = acc;
            }
        }

        double* _pA;
        double* _pB;
        double* _pMtx;


        unsafe long UnsafeLoop(int i0, int iE) {
            double* pA = _pA;
            double* pB = _pB;
            double* pMtx = _pMtx;

            long st = DateTime.Now.Ticks;
            int _N = N;

            for (int i = i0; i < iE; i++) {
                double acc = 0;
                for (int j = 0; j < _N; j++) {
                    acc += pMtx[i * _N + j] * pA[j];
                }
                pB[i] = acc;
            }
            long en = DateTime.Now.Ticks;


            return (en - st);
        }

        void UnsafeRoundRobin(int rnk, int sz) {
            double* pA = _pA;
            double* pB = _pB;
            double* pMtx = _pMtx;
                      

            for (int i = rnk; i < N; i += sz) {
                double acc = 0;
                for (int j = 0; j < N; j++)
                    acc += pMtx[i * N + j] * pA[j];
                    //acc += pMtx[i * N + j] * Math.Pow(-1, j + i*N);
                pB[i] = acc;
            }
        }

        int N = 5000;
        double[] A;
        double[] B;

        double[,] Mtx;

        static void Main(string[] args) {
            //ilPSP.Enviroment.Bootstrap(args, out dummy);
            (new ThreadTest()).Run();
        }

        void Run() {



            N = 5000;
            A = new double[N];
            B = new double[N];

            Mtx = new double[N, N];
            Random rnd = new Random(0);
            for (int i = 0; i < N; i++) {
                A[i] = rnd.NextDouble();
                for (int j = 0; j < N; j++) {
                    Mtx[i, j] = rnd.NextDouble();
                }
            }

            {
                _pA = (double*)Marshal.AllocHGlobal(sizeof(double) * N);
                _pB = (double*)Marshal.AllocHGlobal(sizeof(double) * N);
                _pMtx = (double*)Marshal.AllocHGlobal(sizeof(double) * N * N);
            }

            long[] parRunTime = new long[System.Environment.ProcessorCount];
            {

                for( int kase = 1; kase <= 3; kase++) {
                    Console.WriteLine("=================================================");
                    switch (kase) {
                        case 1: Console.WriteLine("Managed FOR"); break;
                        case 2: Console.WriteLine("Unsafe Round-Robin FOR"); break;
                        case 3: Console.WriteLine("Unsafe FOR"); break;
                    }
                    Console.WriteLine("=================================================");


                    for (int run = 0; run < 2; run++) {
                        for (int numThr = 0; numThr < parRunTime.Length; numThr++) {

                            Paralleism.NumThreads = numThr + 1;
                            long st;// = DateTime.Now.Ticks;
                            QueryPerformanceCounter(out st);

                            //GCHandle pinA = GCHandle.Alloc(A, GCHandleType.Pinned);
                            //GCHandle pinB = GCHandle.Alloc(B, GCHandleType.Pinned);
                            //GCHandle pinMtx = GCHandle.Alloc(Mtx, GCHandleType.Pinned);
                            //_pA = (double*)Marshal.UnsafeAddrOfPinnedArrayElement(A, 0);
                            //_pB = (double*)Marshal.UnsafeAddrOfPinnedArrayElement(B, 0);
                            //_pMtx = (double*)Marshal.UnsafeAddrOfPinnedArrayElement(Mtx, 0);

                            long t = 0;
                            switch(kase) {
                                case 1: Paralleism.For(0, N, this.Loop); break;
                                case 2: Paralleism.Run(this.UnsafeRoundRobin); break;
                                case 3: t = Paralleism.ReduceFor<long>(0, N, this.UnsafeLoop,
                                    delegate(ref long tot, ref long loc) {
                                        tot += loc;   
                                    }); break;
                            }

                            //pinA.Free();
                            //pinB.Free();
                            //pinMtx.Free();

                            //Paralleism.For(0, N, this.Loop);

                            long en;// = DateTime.Now.Ticks; ;
                            QueryPerformanceCounter(out en);

                            //if (kase != 3) {
                            parRunTime[numThr] = en - st;
                            //} else {
                            //    parRunTime[numThr] = t / (numThr + 1);
                            //}

                            if (run > 0) {
                                Console.Write((numThr + 1));
                                Console.Write(" threads: ");
                                Console.Write(parRunTime[numThr]);
                                Console.Write(" msec; ");

                                if (numThr > 0) {
                                    double speedup = 1.0 / ((double)parRunTime[numThr] / (double)parRunTime[0]);
                                    Console.Write("speedup = ");
                                    Console.Write(speedup);
                                }
                                Console.WriteLine();
                            }
                        }
                    }

                    Console.WriteLine("--------------------");
                }
            }






        }
    }
}
