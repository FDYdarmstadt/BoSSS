using System;
using System.IO;
using System.Runtime.InteropServices;
using BoSSS.Foundation.Grid;
using NUnit.Framework;

namespace BoSSS.Application.ExternalBinding {
    /// <summary>Checks the managed degree selection independently of native field storage.</summary>
    [TestFixture, NonParallelizable]
    public static class FieldDegreeTests {
        static Initializer initializer;
        static OpenFOAMGrid grid;

        [OneTimeSetUp]
        public static void Init() {
            initializer = new Initializer();
            initializer.BoSSSInitialize();
            grid = GridImportTestSmall.GenerateFOAMGrid();
        }

        [OneTimeTearDown]
        public static void Cleanup() {
            grid = null;
            initializer.BoSSSFinalize();
        }

        sealed class TemporaryCase : IDisposable {
            readonly string originalDirectory = Directory.GetCurrentDirectory();
            readonly string directory = Path.Combine(Path.GetTempPath(), "BoSSS-field-degree-" + Guid.NewGuid());
            public TemporaryCase(string json = null) {
                Directory.CreateDirectory(Path.Combine(directory, "bosss"));
                if (json != null) File.WriteAllText(Path.Combine(directory, "bosss", "config.json"), json);
                Directory.SetCurrentDirectory(directory);
            }
            public void Dispose() {
                Directory.SetCurrentDirectory(originalDirectory);
                Directory.Delete(directory, true);
            }
        }

        [TestCase(1, 4)]
        [TestCase(2, 10)]
        [TestCase(3, 20)]
        public static void ManagedDegreeOverridesNativeStorage(int managedDegree, int modes) {
            using (new TemporaryCase("{\"DgDgrees\":{\"c\":" + managedDegree + "}}")) {
                IntPtr name = Marshal.StringToCoTaskMemUTF8("c");
                try {
                    // Exercise the same pointer/name constructor used by generated native adapters.
                    var field = new OpenFoamDGField(grid, 2, 1, name);
                    Assert.AreEqual(managedDegree, field.Fields[0].Basis.Degree);
                    Assert.AreEqual(modes, field.Fields[0].Basis.Length);
                    Assert.AreEqual("c", field.Fields[0].Identification);
                    var matrix = new OpenFoamMatrix(grid, field);
                    Assert.AreEqual(managedDegree, matrix.ColMap.BasisS[0].Degree);

                    field.SetDGcoordinate(0, 0, modes - 1, 0.125);
                    // OpenFOAM reads only its ten native slots. Its view must not resize
                    // the managed field or erase modes that never cross the boundary.
                    for (int n = 0; n < 10; n++) {
                        double value = field.GetDGcoordinate(0, 0, n);
                        if (n < modes) field.SetDGcoordinate(0, 0, n, value);
                        else Assert.AreEqual(0, value);
                    }
                    field.SetMean(0, 0, 0.375);
                    Assert.AreEqual(0.125, field.GetDGcoordinate(0, 0, modes - 1));
                    Assert.AreEqual(0.375, field.GetMean(0, 0), 1e-12);
                    Assert.AreEqual(modes, field.Fields[0].Basis.Length);
                } finally { Marshal.FreeCoTaskMem(name); }
            }
        }

        [Test]
        public static void MissingConfigurationUsesNativeDegree() {
            using (new TemporaryCase()) {
                var field = new OpenFoamDGField(grid, 2, 1, "c");
                Assert.AreEqual(2, field.Fields[0].Basis.Degree);
            }
        }

        [Test]
        public static void UnmatchedAndUnnamedFieldsKeepNativeDegree() {
            using (new TemporaryCase("{\"DgDgrees\":{\"c\":3}}")) {
                var velocity = new OpenFoamDGField(grid, 2, 3, "u");
                var unnamed = new OpenFoamDGField(grid, 2, 1);
                Assert.AreEqual(2, velocity.Fields[0].Basis.Degree);
                Assert.AreEqual(2, unnamed.Fields[0].Basis.Degree);
            }
        }

        [TestCase("{broken")]
        [TestCase("{\"DgDgrees\":{\"c\":-1}}")]
        [TestCase("{\"DgDgrees\":{\"c\":\"3\"}}")]
        public static void InvalidOverrideReportsCaseFile(string json) {
            using (new TemporaryCase(json)) {
                var error = Assert.Throws<InvalidDataException>(() => new OpenFoamDGField(grid, 2, 1, "c"));
                Assert.That(error.Message, Does.Contain("config.json"));
                Assert.That(error.Message, Does.Contain("'c'"));
            }
        }
    }
}
