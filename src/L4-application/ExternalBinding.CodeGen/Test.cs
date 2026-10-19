using NUnit.Framework;
using System;
using System.IO;
using System.Linq;

namespace BoSSS.Application.ExternalBinding.CodeGen {
    [TestFixture]
    public static class Test {
        /// <summary>Regeneration is deterministic, preserves unrelated files and emits both sides.</summary>
        [Test]
        public static void RunGenerator() {
            string directory = Path.Combine(Path.GetTempPath(), "BoSSS-codegen-" + Guid.NewGuid());
            Directory.CreateDirectory(directory);
            try {
                string sentinel = Path.Combine(directory, "user.cpp");
                File.WriteAllText(sentinel, "preserve me");
                CodeGenMain.Generate(directory, directory);
                string managed = File.ReadAllText(Path.Combine(directory, "NativeExports.g.cs"));
                Assert.That(managed, Does.Contain("UnmanagedCallersOnly"));
                Assert.That(managed, Does.Contain("FixedOperators_CahnHilliard"));
                var snapshot = Directory.GetFiles(directory).ToDictionary(Path.GetFileName, File.ReadAllText);
                CodeGenMain.Generate(directory, directory);
                foreach (var file in snapshot)
                    Assert.AreEqual(file.Value, File.ReadAllText(Path.Combine(directory, file.Key)), file.Key);
                Assert.AreEqual("preserve me", File.ReadAllText(sentinel));
                Assert.That(File.ReadAllText(Path.Combine(directory, "Globals.cpp")), Does.Contain("get_hostfxr_path"));
            } finally { Directory.Delete(directory, true); }
        }
    }
}
