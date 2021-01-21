using System.Collections.Generic;

namespace BoSSS.Foundation.XDG.OperatorFactory {
    
    
    /// <summary>
    /// A parameter (for a spatial operator) is some field data (<see cref="ISpatialOperator.ParameterVar"/>),
    /// upon which the PDE depends on; It is neither a residual nor some field to solve for.
    /// A typical example would be a non-constant diffusion coefficient.
    /// 
    /// This class (or factory) provides means for updating and generating parameter fields 
    /// </summary>
    /// <remarks>
    /// This class is a member of the operator factory framework, 
    /// which provides a driver framework (to specify entire equations at once)
    /// upon the fine-grained specification for the spatial operator (where equations are build up from individual components).
    /// </remarks>
    public abstract class ParameterS {
        
        /// <summary>
        /// To correlate the parameters handled by this class
        /// with the entire parameter list of the spatial operator (<see cref="ISpatialOperator.ParameterVar"/>)
        /// </summary>
        public abstract IList<string> ParameterNames { get; }

        /// <summary>
        /// Allocation of DG fields; 
        /// </summary>
        public abstract DelParameterFactory Factory { get; }

        /// <summary>
        /// Update of DG fields;
        /// </summary>
        public DelPartialParameterUpdate Update;
    }

    class ParameterList {
        List<ParameterS> parameters;

        public ParameterList(int capacity = 10) {
            parameters = new List<ParameterS>(capacity);
        }

        public void AddParameter(ParameterS parameter) {
            parameters.Add(parameter);
        }

        public ICollection<DelParameterFactory> Factories(IList<string> names) {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelParameterFactory> parameterFactories = new LinkedList<DelParameterFactory>();
            //Find parameters and remove all found parameters from list;

            while(nameList.Count > 0) {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for(int i = 0; i < parameters.Count; ++i) {
                    ParameterS parameter = parameters[i];
                    if(parameter.ParameterNames.Contains(name)) {
                        if(parameter.Factory != null) {
                            parameterFactories.AddLast(parameter.Factory);
                        }
                        foreach(string otherParamName in parameter.ParameterNames) {
                            nameList.Remove(otherParamName);
                        }
                        break;
                    }
                }
            }
            return parameterFactories;
        }

        public ICollection<DelPartialParameterUpdate> ParameterUpdates(IList<string> names) {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelPartialParameterUpdate> parameterUpdates = new LinkedList<DelPartialParameterUpdate>();

            //Find parameters and remove all found parameters from list;
            while(nameList.Count > 0) {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for(int i = 0; i < parameters.Count; ++i) {
                    ParameterS parameter = parameters[i];
                    if(parameter.ParameterNames.Contains(name)) {
                        if(parameter.Update != null) {
                            parameterUpdates.AddLast(parameter.Update);
                        }
                        foreach(string otherParamName in parameter.ParameterNames) {
                            nameList.Remove(otherParamName);
                        }
                        break;
                    }
                }
            }
            return parameterUpdates;
        }
    }
}
