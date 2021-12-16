using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Control {

    /// <summary>
    /// A constant value (i.e. constant in space and time) for initial values/boundary values
    /// </summary>
    [Serializable]
    [DataContract]
    public class ConstantValue : IBoundaryAndInitialData {

        /// <summary>
        /// Ctor
        /// </summary>
        public ConstantValue(double val) {
            DaValue = val;
        }

        /// <summary>
        /// The specified value
        /// </summary>
        public double DaValue {
            get;
            set;
        }

        /// <summary>
        /// 
        /// </summary>
        public double Evaluate(double[] X, double t) {
            return DaValue;
        }

        /// <summary>
        /// 
        /// </summary>
        public void Evaluate(MultidimensionalArray input, double time, MultidimensionalArray output) {
            output.SetAll(DaValue);
        }
    }
}
