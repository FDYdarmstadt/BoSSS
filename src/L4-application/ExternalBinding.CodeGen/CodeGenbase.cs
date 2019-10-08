using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {




    abstract class CodeGenBase {

        public string FileName;
        
        public List<string> IncludeDirectives = new List<string>();

        public List<BracedSection> MainCode = new List<BracedSection>();

    }
}
