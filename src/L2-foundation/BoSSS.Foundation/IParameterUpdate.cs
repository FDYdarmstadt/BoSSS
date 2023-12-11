using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation {

    /// <summary>
    /// Update of parameter fields (e.g. when computing finite difference Jacobian, see e.g. <see cref="DifferentialOperator.GetFDJacobianBuilder"/>),
    /// event used in <see cref="IDifferentialOperator.ParameterUpdates"/>.
    /// </summary>
    /// <param name="DomainVarFields">
    /// Input fields, current state of domain variables.
    /// - key: variable name as defined by the operator (<see cref="IDifferentialOperator.DomainVar"/>)
    /// - value: DG field to store the respective state
    /// </param>
    /// <param name="ParameterVarFields">
    /// Output fields, updated states of parameter fields
    /// - key: variable name as defined by the operator (<see cref="IDifferentialOperator.ParameterVar"/>)
    /// - value: DG field to store the respective state
    /// </param>
    /// <param name="phystime">
    /// physical timestamp
    /// </param>
    /// <remarks>
    /// Note:
    /// 1. Alternatively, equation components which implement <see cref="IParameterHandling"/> can be used.
    /// 2. For an <em>event</em>, such as <see cref="IDifferentialOperator.ParameterUpdates"/> multiple handlers can be added so it is not necessary to 
    ///    put the update of all parameter fields for the operator into one big piece of spaghetti code.
    ///    Hence, it can be split among different handlers.
    /// </remarks>
    public delegate void DelPartialParameterUpdate(double phystime, IReadOnlyDictionary<string,DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields);


    /// <summary>
    /// Factory for the allocation of storage for storing the parameters for this component
    /// (alternatively, <see cref="IParameterHandling.MyParameterAlloc"/> can be used.)
    /// </summary>
    /// <param name="DomainVarFields"></param>
    /// <returns>
    /// a list of pairs, containing:
    /// - parameter name: must match one name in <see cref="IDifferentialOperator.ParameterVar"/>
    /// - a DG field to store the respective parameter
    /// </returns>
    public delegate (string ParameterName, DGField ParamField)[] DelParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields);


   
}
