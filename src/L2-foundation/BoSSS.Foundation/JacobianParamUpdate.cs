using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Foundation {

    /// <summary>
    /// Utility to define parameter variables for Jacobian operators (<see cref="SpatialOperator.GetJacobiOperator"/>)
    /// </summary>
    public class JacobianParamUpdate : IParameterUpdate {

        /// <summary>
        /// Not intended for direct user interaction, but for use by <see cref="SpatialOperator.GetJacobiOperator"/>
        /// </summary>
        public JacobianParamUpdate(IEnumerable<string> __DomainVar, IEnumerable<string> __ParamterVar, List<IEquationComponent> comps, Func<IEquationComponent,TermActivationFlags> extractTaf, int SpatialDimension, DelParameterUpdate __originalUpdate) {
            Components.AddRange(comps);
            originalUpdate = __originalUpdate;
            DomainVar = __DomainVar.ToArray();
            var ParameterVar = __ParamterVar != null ? __ParamterVar.ToArray() : new string[0];
            FindJacobianParams(SpatialDimension, ParameterVar, extractTaf);
        }

        
        /// <summary>
        /// Implementation of constructor functionality.
        /// </summary>
        private void FindJacobianParams(int SpatialDimension, string[] ParameterVar, Func<IEquationComponent,TermActivationFlags> extractTaf) {
            string[] DomVar = DomainVar.ToArray();
            bool[] DomVarAsParam = new bool[DomVar.Length];
            bool[] DomVarDerivAsParam = new bool[DomVar.Length];
            NoOfOrgParams = ParameterVar.Length;

            foreach (var eq in Components) {

                if (!(eq is ISupportsJacobianComponent))
                    throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));

                if (eq.ArgumentOrdering != null) {

                    // collect term activation flags
                    TermActivationFlags termActivationFlags = extractTaf(eq);
                    
                    foreach (string var in eq.ArgumentOrdering) {
                        int iVar = DomVar.IndexOf(var, (a, b) => a.Equals(b));

                        if ((termActivationFlags & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)
                            DomVarAsParam[iVar] |= true;
                        if ((termActivationFlags & (TermActivationFlags.GradUxV | TermActivationFlags.GradUxGradV)) != 0)
                            DomVarDerivAsParam[iVar] |= true;
                    }
                }
            }
            

            string[] newParamVar = new string[0];
            DomainToParam = new int[DomVar.Length];
            DomainDerivToParam = new int[DomVar.Length, SpatialDimension];
            DomainToParam.SetAll(-1234);
            DomainDerivToParam.SetAll(-1235);

            // insert parameters from the original operator at the beginning
            newParamVar = newParamVar.Cat(ParameterVar);

            // parameters for field values in the middle 
            for (int iVar = 0; iVar < DomVar.Length; iVar++) {
                if (DomVarAsParam[iVar]) {
                    DomainToParam[iVar] = newParamVar.Length;
                    newParamVar = newParamVar.Cat(DomVar[iVar] + "_lin");
                }
            }

            // after this, the gradients
            for (int iVar = 0; iVar < DomVar.Length; iVar++) {
                if (DomVarDerivAsParam[iVar]) {
                    for (int d = 0; d < SpatialDimension; d++) {
                        DomainDerivToParam[iVar, d] = newParamVar.Length;
                        newParamVar = newParamVar.Cat(DomVar[iVar] + "_lin_d[" + d + "]");
                    }
                }
            }

            JacobianParameterVars = newParamVar;
        }

        /// <summary>
        /// List of parameter variables required for the Jacobian. 
        /// </summary>
        public string[] JacobianParameterVars {
            get;
            private set;
        }

        DelParameterUpdate originalUpdate;

        int NoOfOrgParams;

        List<IEquationComponent> Components = new List<IEquationComponent>();

        string[] DomainVar;

        int[] DomainToParam;

        int[,] DomainDerivToParam;

        /// <summary>
        /// Parameter update function, 
        /// compatible with <see cref="DelParameterUpdate"/>.
        /// </summary>
        virtual public void ParameterUpdate(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar) {
            DGField[] __DomainVar = DomainVar.ToArray();
            DGField[] __ParameterVar = ParameterVar.ToArray();

            if(originalUpdate != null) {
                originalUpdate(DomainVar, __ParameterVar.GetSubVector(0, NoOfOrgParams));
            }


            int D = __DomainVar[0].GridDat.SpatialDimension;
            if (D != DomainDerivToParam.GetLength(1))
                throw new ApplicationException("spatial dimension mismatch.");

            if (__DomainVar.Length != DomainToParam.Length)
                throw new ApplicationException("mismatch in number of domain variables");

            for (int i = 0; i < __DomainVar.Length; i++) {
                int iDest = DomainToParam[i];
                if (iDest < 0)
                    continue;

                DGField src = __DomainVar[i];
                DGField dst = __ParameterVar[iDest];

                if (object.ReferenceEquals(src, dst))
                    continue;

                dst.Clear();
                dst.Acc(1.0, src);
            }


            for (int i = 0; i < __DomainVar.Length; i++) {
                for (int d = 0; d < D; d++) {

                    int iDest = DomainDerivToParam[i, d];
                    if (iDest < 0)
                        continue;

                    DGField src = __DomainVar[i];
                    DGField dst = __ParameterVar[iDest];

                    if (object.ReferenceEquals(src, dst))
                        throw new ApplicationException("Separate DG field must be allocated to store Parameters.");

                    dst.Clear();
                    dst.Derivative(1.0, src, d);
                }
            }
        }

        /// <summary>
        /// creates clones of the domain fields to store parameter fields
        /// </summary>
        virtual public DGField[] AllocateParameters(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar) {
            DGField[] ret = new DGField[this.JacobianParameterVars.Length];
            DGField[] __DomainVar = DomainVar.ToArray();
            DGField[] __ParameterVar = ParameterVar != null ? ParameterVar.ToArray() : new DGField[0];
            
            if (DomainVar.Count() != DomainToParam.Length)
                throw new ApplicationException("mismatch in number of domain variables");
            if (__ParameterVar.Count() != NoOfOrgParams)
                throw new ApplicationException("mismatch in number of parameter variables");
            int GAMMA = DomainToParam.Length;
            int D = DomainVar.First().GridDat.SpatialDimension;
            if (D != DomainDerivToParam.GetLength(1))
                throw new ApplicationException("spatial dimension mismatch.");

            for(int i = 0; i < __ParameterVar.Length; i++) {
                ret[i] = __ParameterVar[i];
            }
            
            for (int i = 0; i < __DomainVar.Length; i++) {
                int iDest = DomainToParam[i];
                if (iDest < 0)
                    continue;

                DGField src = __DomainVar[i];
                DGField dst = src;

                ret[iDest] = dst;
            }


            for (int i = 0; i < __DomainVar.Length; i++) {
                for (int d = 0; d < D; d++) {

                    int iDest = DomainDerivToParam[i, d];
                    if (iDest < 0)
                        continue;

                    DGField src = __DomainVar[i];
                    DGField dst = src.CloneAs();
                    dst.Identification = dst.Identification + "_d[" + d + "]";

                    ret[iDest] = dst;
                }
            }

            return ret;
        }
    }
}
