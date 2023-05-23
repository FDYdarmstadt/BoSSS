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
        [DataMember]
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
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            output.SetAll(DaValue);
        }

        /// <summary>
        /// true, if the specified values are approximately equal
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as ConstantValue;
            if(other == null)
                return false;

            return this.DaValue.ApproxEqual(other.DaValue);
        }

        /// <summary>
        /// 
        /// </summary>
        public override int GetHashCode() {
            return (int)(DaValue * 12378.1234);
        }

    }
}
