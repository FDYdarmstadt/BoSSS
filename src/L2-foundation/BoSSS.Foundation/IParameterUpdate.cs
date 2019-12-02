using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation {

    /// <summary>
    /// Update of parameter fields (e.g. when computing finite difference Jacobian, see e.g. <see cref="SpatialOperator.GetFDJacobianBuilder"/>).
    /// </summary>
    /// <param name="DomainVar">
    /// Input fields, current state of domain variables
    /// </param>
    /// <param name="ParameterVar">
    /// Output fields, updated states of parameter fields
    /// </param>
    public delegate void DelParameterUpdate(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar);


    public interface IParameterUpdate {


        /// <summary>
        /// Allocates DG fields to store parameters
        /// </summary>
        DGField[] AllocateParameters(IEnumerable<DGField> DomainVar);
        
        /// <summary>
        /// Update of parameter fields (e.g. when computing finite difference Jacobian, see e.g. <see cref="SpatialOperator.GetFDJacobianBuilder"/>).
        /// </summary>
        /// <param name="DomainVar">
        /// Input fields, current state of domain variables
        /// </param>
        /// <param name="ParameterVar">
        /// Output fields, updated states of parameter fields
        /// </param>
        void DelParameterUpdate(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar);
    }
}
