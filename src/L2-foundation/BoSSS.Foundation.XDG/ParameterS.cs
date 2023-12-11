using System.Collections.Generic;

namespace BoSSS.Foundation.XDG.OperatorFactory {
    
    
    /// <summary>
    /// A parameter (for a spatial operator) is some field data (<see cref="IDifferentialOperator.ParameterVar"/>),
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
        /// with the entire parameter list of the spatial operator (<see cref="IDifferentialOperator.ParameterVar"/>)
        /// </summary>
        public abstract IList<string> ParameterNames { get; }

        /// <summary>
        /// Allocation of DG fields; 
        /// </summary>
        public abstract DelParameterFactory Factory { get; }

        /// <summary>
        /// Update of DG fields;
        /// </summary>
        virtual public DelPartialParameterUpdate Update { 
            get { return null; }
        }
        
    }

    
}
