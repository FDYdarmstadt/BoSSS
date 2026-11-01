using System;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace BoSSS.Application.BoSSSpad {
    /// <summary>Regression checks for associating MPI memory traces with run labels.</summary>
    [TestFixture]
    public static class MemtraceAnalysisTests {
        /// <summary>Shuffled and equal-sized runs retain the memory values of their own session.</summary>
        [Test]
        public static void MemoryValuesFollowSortedMPISizes() {
            // Include shuffled sizes and repeated sizes with different memory values.
            // Equal-sized runs must retain their original relative order.
            foreach (int[] sizes in new[] { new[] { 3, 1, 2 }, new[] { 2, 1, 2 } }) {
                string root = Path.Combine(Directory.GetCurrentDirectory(), "memtrace-test-" + Guid.NewGuid());
                Directory.CreateDirectory(root);
                try {
                    var directories = sizes.Select((size, run) => {
                        var dir = Directory.CreateDirectory(Path.Combine(root, "run" + run));
                        for (int rank = 0; rank < size; rank++) {
                            long bytes = (run + 1L) * 1024 * 1024;
                            File.WriteAllLines(Path.Combine(dir.FullName, $"memory.{rank}.txt"), new[] {
                                $"0\t{bytes}\t0\t0\t0\tcheckpoint-a",
                                $"1\t{2 * bytes}\t0\t0\t0\tcheckpoint-b"
                            });
                        }
                        return dir;
                    }).ToArray();
                    var actual = new SessionsComparisonMemtrace(directories).GetTotalMemory();
                    var order = Enumerable.Range(0, sizes.Length).OrderBy(i => sizes[i]).ToArray();
                    Assert.AreEqual(order.Length, actual.Length);
                    for (int i = 0; i < order.Length; i++) {
                        int run = order[i];
                        long bytes = sizes[run] * (run + 1L) * 1024 * 1024;
                        Assert.AreEqual(sizes[run], actual[i].MPISz);
                        CollectionAssert.AreEqual(new[] { bytes, 2 * bytes }, actual[i].TotalMem);
                    }
                } finally {
                    Directory.Delete(root, recursive: true);
                }
            }
        }
    }
}
