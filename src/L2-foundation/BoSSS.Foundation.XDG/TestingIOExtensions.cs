using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {


    /// <summary>
    /// Another extension class 
    /// </summary>
    public static class TestingIOExtensions {

        /// <summary>
        /// Adds an XDG field.
        /// </summary>
        public static void AddDGField(this TestingIO t, XDGField f) {
            var trk = f.Basis.Tracker;

            foreach(string spc in trk.SpeciesNames) {
                var fs = f.GetSpeciesShadowField(spc);
                t.AddDGField(fs);
            }
        }

    }
}
