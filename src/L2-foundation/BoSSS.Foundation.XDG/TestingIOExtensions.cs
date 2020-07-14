using BoSSS.Foundation.Grid;
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
                SinglePhaseField fsFatClone = new SinglePhaseField(fs.Basis, fs.Identification);
                CellMask msk = trk.Regions.GetSpeciesMask(spc);

                fsFatClone.Acc(1.0, fs, msk);

                t.AddDGField(fsFatClone);
            }
        }
        
        /// <summary>
        /// Adds an XDG field.
        /// </summary>
        public static XDGField LocalError(this TestingIO t, XDGField f) {
            var trk = f.Basis.Tracker;

            var ErrLoc = f.CloneAs();
            ErrLoc.Clear();

            foreach(string spc in trk.SpeciesNames) {
                var fs = ErrLoc.GetSpeciesShadowField(spc);
                t.OverwriteDGField(fs);
            }

            ErrLoc.Scale(-1);
            ErrLoc.Acc(1.0, f);
            
            ErrLoc.Identification = "Error-" + ErrLoc.Identification;
            return ErrLoc;
        }

        /// <summary>
        /// Overwrites the memory of a XDG field with the reference data 
        /// </summary>
        static public void OverwriteDGField(this TestingIO t, XDGField f) {
             var trk = f.Basis.Tracker;

            foreach(string spc in trk.SpeciesNames) {
                var fs = f.GetSpeciesShadowField(spc);
                t.OverwriteDGField(fs);
            }
        }

    }
}
