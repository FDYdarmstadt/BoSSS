﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {




    abstract class CodeGenBase {

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
                foreach(var l in GetCode()) {
                    stw.WriteLine(l);
                }
                return stw.ToString();
            }
        }

        public void WriteFile() {
            File.WriteAllText(FileName + FilenameExt, this.ToString());
        }
    }
}
