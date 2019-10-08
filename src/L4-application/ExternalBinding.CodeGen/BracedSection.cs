using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding.CodeGen {
    class BracedSection {

        public bool NoBraces = false;

        public List<string> OutsideCode = new List<string>();

        public List<object> Children = new List<object>();
        
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
                ret.Add("}");

            return ret.ToArray();
        }

    }
}
