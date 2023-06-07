using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation {

    /// <summary>
    /// Extension functions for <see cref="ISpatialOperator"/>
    /// </summary>
    static public class ISpatialOparatorExtensions {


        /// <summary>
        /// Involves all <see cref="ISpatialOperator.ParameterFactories"/> events 
        /// and all <see cref="IParameterHandling.MyParameterAlloc"/> methods in the operators equation components 
        /// in order to allocate operator storage.
        /// </summary>
        public static DGField[] InvokeParameterFactory(this ISpatialOperator op, IEnumerable<DGField> __DomainFields) {

            if (!op.IsCommitted)
                throw new NotSupportedException("not allowed before commit.");
            var _DomainFields = __DomainFields.ToArray();
            if (_DomainFields.Length != op.DomainVar.Count) {
                string fl_domNames = _DomainFields.Select(f => f?.Identification ?? "NULL").ToConcatString("[", ",", "]");
                string op_domNames = op.DomainVar.ToConcatString("[", ",", "]");
                throw new ArgumentException($"Mismatch in number of domain variables: provided domain fields {fl_domNames}, specified by operator {op_domNames}.");
            }

            int NoOfParams = op.ParameterVar.Count;
            DGField[] ret = new DGField[NoOfParams];

            var DomainVarsDict = new Dictionary<string, DGField>();
            for (int iVar = 0; iVar < _DomainFields.Length; iVar++) {
                DomainVarsDict.Add(op.DomainVar[iVar], _DomainFields[iVar]);
            }

            // invoke factories by the operator
            if (op.ParameterFactories != null) {
                foreach (DelParameterFactory Factory in op.ParameterFactories) {
                    var ttt = Factory(DomainVarsDict);

                    foreach (var tt in ttt) {
                        int idx = op.ParameterVar.IndexOf(tt.ParameterName);
                        if (idx < 0)
                            throw new Exception($"Illegal parameter name {tt.ParameterName} provided by parameter factory -- not in the operator parameter list.");
                        ret[idx] = tt.ParamField;
                    }
                }
            }

            bool ComponentFactoryUseful(IParameterHandling _ph, out int[] targidx) {
                var phParams = _ph.ParameterOrdering;
                if (phParams == null) {
                    targidx = new int[0];
                    return false;
                }

                targidx = new int[phParams.Count];
                int c = 0;
                bool useful = false;
                foreach(string phParamName in phParams) {
                    int idx = op.ParameterVar.IndexOf(phParamName);
                    if (idx < 0)
                        throw new ApplicationException("should not happen if operator is committed and verified");
                    targidx[c] = idx;
                    if (ret[idx] == null)
                        useful = true; // some parameter has not been allocated yet
                    c++;
                }
                return useful;
            }

            DGField[] GetArguments(IEquationComponent c) {
                var cArgs = c.ArgumentOrdering;
                if (cArgs == null)
                    return new DGField[0];
                DGField[] rr = new DGField[cArgs.Count];
                for(int i = 0; i < rr.Length; i++) {
                    rr[i] = DomainVarsDict[cArgs[i]];
                }
                return rr;
            }

            // invoke factories in equation components
            foreach (string codName in op.CodomainVar) {
                foreach (IEquationComponent comp in op.EquationComponents[codName]) {
                    if (comp is IParameterHandling ph) {
                        if(ComponentFactoryUseful(ph, out int[] targIdx)) {
                            DGField[] newParams = ph.MyParameterAlloc(GetArguments(ph));
                            if (newParams.Length != targIdx.Length)
                                throw new NotSupportedException($"Illegal implementation of parameter allocation for {ph.GetType().Name}: the length of the array returned by parameter allocation must be equal to number of parameters.");
                        
                        
                            for(int i = 0; i < newParams.Length; i++) {
                                if (newParams[i] != null)
                                    ret[targIdx[i]] = newParams[i];
                            }
                        }


                    }
                }
            }

            // return 
            return ret;
        }

        /// <summary>
        /// Involves all <see cref="ISpatialOperator.ParameterUpdates"/> events 
        /// and all <see cref="IParameterHandling.MyParameterUpdate"/> methods in the operators equation components 
        /// in order to allocate operator storage.
        /// </summary>
        public static void InvokeParameterUpdate(this ISpatialOperator op, double time, DGField[] __DomainFields, DGField[] __ParameterFields) {
            using(new FuncTrace()) {
                if(!op.IsCommitted)
                    throw new NotSupportedException("not allowed before commit.");
                if(__DomainFields.Length != op.DomainVar.Count)
                    throw new ArgumentException("Mismatch in number of domain variables.");
                if(__ParameterFields.Length != op.ParameterVar.Count)
                    throw new ArgumentException("Mismatch in number of parameter variables.");

                int NoOfParams = op.ParameterVar.Count;
                DGField[] ret = new DGField[NoOfParams];

                var DomainVarsDict = new Dictionary<string, DGField>();
                for(int iVar = 0; iVar < __DomainFields.Length; iVar++) {
                    DomainVarsDict.Add(op.DomainVar[iVar], __DomainFields[iVar]);
                }

                var ParameterVarsDict = new Dictionary<string, DGField>();
                for(int iVar = 0; iVar < __ParameterFields.Length; iVar++) {
                    ParameterVarsDict.Add(op.ParameterVar[iVar], __ParameterFields[iVar]);
                }

                // invoke update functions in the operator
                foreach(var PartialParameterUpdate in op.ParameterUpdates) {
                    PartialParameterUpdate(time, DomainVarsDict, ParameterVarsDict);
                }

                // invoke update functions in the equation components
                bool[] ParameterIsUpdated = new bool[NoOfParams]; // for the equation components, we mark parameters that are already updated, to avoid doing the update multiple times

                bool ComponentUpdateUseful(IParameterHandling _ph, out DGField[] _ph_args, out DGField[] _ph_params, out int[] targIdx) {
                    var phParams = _ph.ParameterOrdering;
                    targIdx = new int[phParams != null ? phParams.Count : 0];

                    _ph_params = new DGField[phParams != null ? phParams.Count : 0];
                    bool useful = false;
                    if(_ph_params.Length > 0) {
                        int c = 0;
                        foreach(string phParamName in phParams) {
                            int idx = op.ParameterVar.IndexOf(phParamName);
                            if(idx < 0)
                                throw new ApplicationException("should not happen if operator is committed and verified");
                            _ph_params[c] = ParameterVarsDict[phParamName];
                            targIdx[c] = idx;

                            if(ParameterIsUpdated[idx] == false) {
                                useful = true;
                            }

                            if(ret[idx] == null)
                                useful = true; // some parameter has not been allocated yet
                            c++;
                        }
                    }

                    var phArgs = _ph.ArgumentOrdering;
                    _ph_args = new DGField[phArgs != null ? phArgs.Count : 0];
                    if(_ph_args.Length > 0) {
                        int c = 0;
                        foreach(string phArgName in phArgs) {
                            _ph_args[c] = DomainVarsDict[phArgName];
                            c++;
                        }
                    }
                    return useful;
                }


                foreach(string codName in op.CodomainVar) {
                    foreach(IEquationComponent comp in op.EquationComponents[codName]) {
                        if(comp is IParameterHandling ph) {
                            if(ComponentUpdateUseful(ph, out var ph_argFields, out var ph_paramFields, out int[] targIdx)) {
                                ph.MyParameterUpdate(ph_argFields, ph_paramFields);

                                for(int i = 0; i < ph_paramFields.Length; i++) {
                                    if(ph_paramFields[i] != null) // hack: if the 'MyParameterUpdate' sets some entry to null, it signals that it did not updated the parameter variable
                                        ParameterIsUpdated[targIdx[i]] = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
