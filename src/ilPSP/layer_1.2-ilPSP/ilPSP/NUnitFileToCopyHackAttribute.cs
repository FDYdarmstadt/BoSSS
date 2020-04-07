using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP {

    /// <summary>
    /// This is a super-ugly hack to mark data files and directories
    /// within the source code repository
    /// that must be copied for the NUnit tests.
    /// </summary>
    public class NUnitFileToCopyHackAttribute : Attribute {

        /// <summary>
        /// Wäh
        /// </summary>
        public NUnitFileToCopyHackAttribute(params string[] SomeFileNames) {
            this.SomeFileNames = SomeFileNames.CloneAs();
        }

        /// <summary>
        /// partial path or file name within the source repo
        /// </summary>
        public string[] SomeFileNames;
    }
}
