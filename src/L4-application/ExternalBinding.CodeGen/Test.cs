using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {

    /// <summary>
    /// Unit/regression tests
    /// </summary>
    [TestFixture]
    static public class Test {

        /// <summary>
        /// Executes the Code generator to see if it runs without exception.
        /// </summary>
        [Test]
        public static void RunGenerator() {
            CodeGenMain.Main(new[] { "." });
        }
    }
}
