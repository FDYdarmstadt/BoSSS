using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.CompressibleFlowCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNS {

    /// <summary>
    /// Represents auxiliary variables that are updated via some user-defined
    /// update function
    /// </summary>
    public class DerivedVariable : Variable {

        /// <summary>
        /// The update function to be invoked after each (sub-)time-step
        /// </summary>
        public Action<DGField, CellMask, IProgram<CNSControl>> UpdateFunction;

        /// <summary>
        /// See <see cref="Variable"/> and <see cref="UpdateFunction"/>
        /// </summary>
        /// <param name="name"></param>
        /// <param name="type"></param>
        /// <param name="updateFunction"></param>
        public DerivedVariable(string name, VariableTypes type, Action<DGField, CellMask, IProgram<CNSControl>> updateFunction)
            : base(name, type) {
            this.UpdateFunction = updateFunction;
        }
    }
}
