using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {
    class CodeGenHeaderFile : CodeGenBase {

        public const string HeaderFileSuffix = ".h";


        protected override string FilenameExt {
            get {
                return HeaderFileSuffix;
            }
        }
    }
}
