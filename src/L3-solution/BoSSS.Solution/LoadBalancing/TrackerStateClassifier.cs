using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.LoadBalancing {


    /// <summary>
    /// (base-class only)
    /// Classification of cells for load balancing according to the state of the level-set-tracker
    /// </summary>
    abstract public class TrackerStateClassifier : ICellClassifier {

        /// <summary>
        /// A species which is designated to be void, i.e. causes negligible compute load, and won't be counted, etc.
        /// Set to null or an invalid species name, if all species should be considered
        /// </summary>
        public string VoidSpecies = "C";

        public abstract int[] ClassifyCells(IApplication program);

        abstract public object Clone();
    }
}
