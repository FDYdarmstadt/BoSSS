using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {

    abstract class CodeFileBase {

        protected static DateTime CreationTime = DateTime.Now;


        public string FileName;
        
        public List<string> IncludeDirectives = new List<string>();

        

        public List<BracedSection> MainCode = new List<BracedSection>();
        
        protected abstract string FilenameExt {
            get;
        }

        public string[] GetCode() {
            var R = new List<string>();
            R.AddRange(IncludeDirectives);
            foreach(var b in MainCode) {
                R.AddRange(b.GetCode());
            }

            return R.ToArray();
        }

        public override string ToString() {
            using(var stw = new StringWriter()) {

                stw.WriteLine("// ##########################################################################");
                stw.WriteLine("// This is AUTO-GENERATED code created by the ");
                stw.WriteLine("// BoSSS External Language Binding code generator.");
                stw.WriteLine("// **Any manual changes are over-written if the code-generator is executed.**");
                stw.WriteLine("// Creation Date: " + CreationTime);
                stw.WriteLine("// ##########################################################################");
                

                foreach(var l in GetCode()) {
                    stw.WriteLine(l);
                }
                return stw.ToString();
            }
        }

        public void WriteFile(string __Directory) {
            if (!Directory.Exists(__Directory))
                throw new IOException("Output directory is not existent.");

            string FullName = Path.Combine(__Directory, FileName + FilenameExt);


            File.WriteAllText(FullName, this.ToString());
        }
    }
}
