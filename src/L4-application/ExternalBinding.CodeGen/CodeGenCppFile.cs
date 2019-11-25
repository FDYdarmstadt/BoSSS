using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {


    class CodeGenCppFile : CodeFileBase {

        public CodeGenCppFile() {
            IncludeDirectives.Add("#include <stdlib.h>");
            IncludeDirectives.Add("#include <assert.h>");
            IncludeDirectives.Add("#define _CRT_SECURE_NO_WARNINGS");
            IncludeDirectives.Add("#include <string.h>");
            IncludeDirectives.Add("");
            IncludeDirectives.Add("#include <mono/metadata/mono-config.h>");
            IncludeDirectives.Add("#include <mono/jit/jit.h>");
            IncludeDirectives.Add("#include <mono/metadata/assembly.h>");
            IncludeDirectives.Add("#include <mono/metadata/debug-helpers.h>");
            IncludeDirectives.Add("");
            IncludeDirectives.Add("#include \"Globals.h\"");
            
        }



        public const string CppFileSuffix = ".cpp";

        protected override string FilenameExt {
            get {
                return CppFileSuffix;
            }
        }
    }
}
