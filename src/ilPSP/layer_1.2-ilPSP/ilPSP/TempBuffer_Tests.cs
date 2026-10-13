using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;

namespace ilPSP {

    /// <summary>
    /// Tests for <see cref="TempBuffer"/>: thread-local slots, shared-pool fall-back, pool growth and misuse detection.
    /// </summary>
    [TestFixture]
    public static class TempBuffer_Tests {

        /// <summary>
        /// More buffers than thread-local slots are held simultaneously: all handles positive and distinct,
        /// all buffers distinct and large enough; released in non-LIFO order and re-acquired.
        /// </summary>
        [Test]
        public static void NestedBuffersAreDistinct() {
            int N = TempBuffer.NO_OF_THREADLOCAL_SLOTS + 4;
            for(int pass = 0; pass < 3; pass++) {
                var handles = new int[N];
                var bufs = new double[N][];
                for(int i = 0; i < N; i++) {
                    int size = 10 + 7 * i + pass;
                    bufs[i] = TempBuffer.GetTempBuffer(out handles[i], size);
                    Assert.Greater(handles[i], 0, "handles must be positive (callers use 0 as 'no buffer')");
                    Assert.GreaterOrEqual(bufs[i].Length, size);
                    for(int j = 0; j < size; j++) bufs[i][j] = i;
                }
                Assert.AreEqual(N, handles.Distinct().Count(), "handles must be distinct");
                for(int i = 0; i < N; i++)
                    for(int i2 = i + 1; i2 < N; i2++)
                        Assert.IsFalse(ReferenceEquals(bufs[i], bufs[i2]), "buffers must be distinct");
                for(int i = 0; i < N; i++)
                    for(int j = 0; j < 10 + 7 * i + pass; j++)
                        Assert.AreEqual(i, bufs[i][j], "buffer content overwritten");

                // free in non-LIFO order
                foreach(int i in Enumerable.Range(0, N).OrderBy(i => (i * 5) % N))
                    TempBuffer.FreeTempBuffer(handles[i]);
            }
        }

        /// <summary>
        /// Large requests and the exhausted thread-local slots are served by the shared pool;
        /// the pool grows beyond its initial size instead of throwing.
        /// </summary>
        [Test]
        public static void PoolFallbackAndGrowth() {
            // large request: bypasses thread-local slots
            var big = TempBuffer.GetTempBuffer(out int hBig, TempBuffer.MAX_THREADLOCAL_SIZE + 1);
            Assert.Greater(hBig, 0);
            Assert.GreaterOrEqual(big.Length, TempBuffer.MAX_THREADLOCAL_SIZE + 1);
            TempBuffer.FreeTempBuffer(hBig);

            // many small buffers: more than the initial pool size (64*NumThreads)
            int N = 64 * Math.Max(1, ilPSP.Environment.NumThreads) + 50;
            var handles = new List<int>();
            var bufs = new List<double[]>();
            try {
                for(int i = 0; i < N; i++) {
                    bufs.Add(TempBuffer.GetTempBuffer(out int h, 3));
                    handles.Add(h);
                }
                Assert.AreEqual(N, handles.Distinct().Count());
                Assert.IsTrue(handles.All(h => h > 0));
                Assert.AreEqual(N, bufs.Distinct(ReferenceEqualityComparer.Instance).Count());
            } finally {
                foreach(int h in handles)
                    TempBuffer.FreeTempBuffer(h);
            }
        }

        /// <summary>
        /// Freeing a thread-local buffer twice, or on a different thread, is detected.
        /// </summary>
        [Test]
        public static void MisuseIsDetected() {
            TempBuffer.GetTempBuffer(out int h, 17);
            TempBuffer.FreeTempBuffer(h);
            Assert.Throws<ApplicationException>(() => TempBuffer.FreeTempBuffer(h), "double free not detected");

            TempBuffer.GetTempBuffer(out int h2, 17);
            Exception crossThread = null;
            var th = new Thread(() => {
                try {
                    TempBuffer.FreeTempBuffer(h2);
                } catch(Exception e) {
                    crossThread = e;
                }
            });
            th.Start();
            th.Join();
            Assert.IsInstanceOf<ApplicationException>(crossThread, "free on a different thread not detected");
            TempBuffer.FreeTempBuffer(h2); // correct free on the owning thread
        }

        /// <summary>
        /// Concurrent use from many threads: each thread fills its buffers with a thread/iteration specific
        /// value and checks that no other holder has written into them.
        /// </summary>
        [Test]
        public static void ConcurrentUse() {
            int errors = 0;
            Parallel.For(0, 64, new ParallelOptions { MaxDegreeOfParallelism = 8 }, delegate (int task) {
                var rnd = new Random(task);
                for(int it = 0; it < 2000; it++) {
                    int n = 1 + rnd.Next(5);
                    var hs = new int[n];
                    var bs = new double[n][];
                    var sz = new int[n];
                    for(int i = 0; i < n; i++) {
                        sz[i] = 1 + rnd.Next(rnd.Next(10) == 0 ? 5000 : 200);
                        bs[i] = TempBuffer.GetTempBuffer(out hs[i], sz[i]);
                        double v = task * 1e6 + it * 10 + i;
                        for(int j = 0; j < sz[i]; j++) bs[i][j] = v;
                    }
                    for(int i = 0; i < n; i++) {
                        double v = task * 1e6 + it * 10 + i;
                        for(int j = 0; j < sz[i]; j++)
                            if(bs[i][j] != v) { Interlocked.Increment(ref errors); break; }
                    }
                    for(int i = n - 1; i >= 0; i--)
                        TempBuffer.FreeTempBuffer(hs[i]);
                }
            });
            Assert.AreEqual(0, errors, "buffers were shared between concurrent holders");
        }

        sealed class ReferenceEqualityComparer : IEqualityComparer<double[]> {
            public static readonly ReferenceEqualityComparer Instance = new ReferenceEqualityComparer();
            public bool Equals(double[] x, double[] y) => ReferenceEquals(x, y);
            public int GetHashCode(double[] obj) => System.Runtime.CompilerServices.RuntimeHelpers.GetHashCode(obj);
        }
    }
}
