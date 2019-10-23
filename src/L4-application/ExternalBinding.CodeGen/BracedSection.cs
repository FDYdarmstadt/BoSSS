using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {
    class BracedSection {

        public bool NoBraces = false;
        public bool ClosingSemicolon = false;


        public List<string> OutsideCode = new List<string>();

        public void AddOutside(string s) {
            OutsideCode.Add(s);
        }
        public void AddOutside(string s, params object[] os) {
            OutsideCode.Add(string.Format(s, os));
        }
        

        public List<object> Children = new List<object>();

        public void AddInner(string s) {
            Children.Add(s);
        }
        public void AddInner(string s, params object[] os) {
            Children.Add(string.Format(s, os));
        }
        
        
        public string[] GetCode() {
            var ret = new List<string>();

            ret.AddRange(OutsideCode);

            if(!NoBraces)
                ret.Add("{");

            foreach(object s in Children) {
                if(s is string line) {
                    ret.Add(line);
                } else if(s is BracedSection child) {
                    ret.AddRange(child.GetCode());
                } else {
                    throw new NotImplementedException();
                }
            }
            
            if(!NoBraces)
                ret.Add("}" + (ClosingSemicolon ? ";" : ""));
            else {
                if (ClosingSemicolon)
                    throw new NotSupportedException("illegal config - no braces, but semicolon.");
            }

            return ret.ToArray();
        }



        public override string ToString() {
            using(var stw = new StringWriter()) {
                foreach(var l in GetCode()) {
                    stw.WriteLine(l);
                }
                return stw.ToString();
            }
        }

    }
}
