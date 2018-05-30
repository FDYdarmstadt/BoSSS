/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;

using BoSSS.Platform;
using ilPSP.Utils;
using System.Diagnostics;
using ilPSP;


namespace BoSSS.Foundation.Quadrature.FluxQuadCommon {


    /// <summary>
    /// used by the spatial operator (see <see cref="BoSSS.Foundation.SpatialOperator.GetArgMapping{T}"/>)
    /// to give a collection of all equation components 
    /// of a certain type (<typeparamref name="T"/>)
    /// </summary>
    public class EquationComponentArgMapping<T> where T : IEquationComponent {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="DiffOp"></param>
        /// <param name="CoDomVarName">
        /// the name of the variable in the codomain (<see cref="SpatialOperator.CodomainVar"/>-member
        /// of <paramref name="DiffOp"/>, for which this object should be defined;
        /// </param>
        /// <param name="_fieldList">
        /// list of domain variable names
        /// </param>
        /// <param name="_fieldList2">
        /// optional (i.e. can be null)
        /// list of parameter variable names;
        /// will be concatenated with <paramref name="_fieldList"/>
        /// </param>
        /// <param name="F">optional filter, can be null.</param>
        /// <param name="vectorizer">
        /// Function for the vectorization of the evaluation of <paramref name="F"/>
        /// </param>
        internal EquationComponentArgMapping(SpatialOperator DiffOp, string CoDomVarName,
            IList<string> _fieldList, IList<string> _fieldList2, Func<T, bool> F, Func<IEquationComponent, IEquationComponent> vectorizer) {
            m_CoDomVarName = CoDomVarName;
            //m_DomainFields = DomainMapping.Fields;


            // =================
            // concat field list
            // =================
            IList<string> fieldList;
            if (_fieldList2 != null) {
                fieldList = Enumerable.Concat(_fieldList, _fieldList2).ToList();
            } else {
                fieldList = _fieldList;
            }

            // =========================================
            // collect all equation components of type T
            // =========================================

            ICollection<IEquationComponent> eqCompS = DiffOp.EquationComponents[CoDomVarName];

            List<T> AllComponentsofMyType = new List<T>();

            foreach (IEquationComponent eqComp in eqCompS) {
                if (eqComp is T) {
                    T _eqComp = (T)eqComp;

                    if (F == null || F(_eqComp)) {
                        AllComponentsofMyType.Add(_eqComp);
                    }
                } else if (vectorizer != null) {
                    IEquationComponent VeqComp = vectorizer(eqComp);

                    if (VeqComp != null && VeqComp is T) {
                        T _VeqComp = (T)VeqComp;

                        if (F == null || F(_VeqComp)) {
                            AllComponentsofMyType.Add(_VeqComp);
                        }
                    }
                }
            }

            m_AllComponentsOfMyType = AllComponentsofMyType.ToArray();

            // ======================
            // build argument mapping
            // ======================

            int argMapMaxCount = 0;
            foreach (var w in m_AllComponentsOfMyType) {
                int c = w.ArgumentOrdering.Count;
                if (w.ParameterOrdering != null)
                    c += w.ParameterOrdering.Count;
                argMapMaxCount = Math.Max(argMapMaxCount, c);
            }
            AllToSub = new int[m_AllComponentsOfMyType.Length, argMapMaxCount];
            NoOfArguments = new int[m_AllComponentsOfMyType.Length];
            NoOfParameters = new int[m_AllComponentsOfMyType.Length];
            //ArrayTools.SetAll(AllToSub, int.MinValue);
            AllToSub.SetAll(int.MinValue);

            //IList<string> fieldList = DiffOp.DomainVar;
            for (int i = 0; i < AllToSub.GetLength(0); i++) {
                // arguments...
                IList<string> argMap = m_AllComponentsOfMyType[i].ArgumentOrdering;
                NoOfArguments[i] = argMap.Count;
                for (int j = 0; j < NoOfArguments[i]; j++) {
                    int ifound = fieldList.IndexOf(argMap[j]);
                    AllToSub[i, j] = ifound;
                    Debug.Assert(AllToSub[i, j] >= 0);
                }

                // parameters...
                IList<string> paramMap = m_AllComponentsOfMyType[i].ParameterOrdering;
                if (paramMap != null) {
                    NoOfParameters[i] = paramMap.Count;

                    for (int j = 0; j < NoOfParameters[i]; j++) {
                        AllToSub[i, j + NoOfArguments[i]] = fieldList.IndexOf(paramMap[j]);
                        Debug.Assert(AllToSub[i, j + NoOfArguments[i]] >= 0);
                    }
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public T[] m_AllComponentsOfMyType;

        /// <summary>
        /// 
        /// </summary>
        string m_CoDomVarName;

        /// <summary>
        /// reorders the values of the domain variables (<paramref name="FieldValues"/>)
        /// to comply with the argument and paramter ordering of component <paramref name="comp"/>.
        /// </summary>
        public MultidimensionalArray[] MapArguments(MultidimensionalArray[] FieldValues, T comp, bool OnlyArgs = false) {
            int CompInd = Array.IndexOf<T>(m_AllComponentsOfMyType, comp);

            int cnt = NoOfArguments[CompInd] + (OnlyArgs ? 0 : NoOfParameters[CompInd]);
            //for (cnt = 0; cnt < AllToSub.GetLength(1); cnt++)
            //    if (AllToSub[CompInd, cnt] < 0)
            //        break;
            MultidimensionalArray[] ret = new MultidimensionalArray[cnt];

            for (int i = 0; i < cnt; i++)
                ret[i] = FieldValues[AllToSub[CompInd, i]];

            return ret;
        }

        /// <summary>
        /// content: the index of a DG - field (global field index within the integrator)<br/>
        /// <list type="bullet">
        ///   <item> 1st index: equation component index; the order of the equation
        ///          components correlates, resp. is induced by <see cref="m_AllComponentsOfMyType"/>;
        ///   </item>
        ///   <item> 2nd index: a field index, where the order of the fields is defined by 
        ///          <see cref="IEquationComponent.ArgumentOrdering"/>.
        ///          (equation-specific field index)
        ///   </item>
        /// </list>
        /// </summary>
        public int[,] AllToSub;

        ///// <summary>
        ///// The number of arguments (<see cref="IEquationComponent.ArgumentOrdering"/>.Count),
        ///// for each equation component.
        ///// 1st index: equation component index; the order of the equation
        /////          components is defined by <see cref="m_AllComponentsOfMyType"/>;
        ///// </summary>
        //internal int[] NumberOfCodomArgs;

        /// <summary>
        /// content: number of arguments for each equation component; 
        /// (<see cref="IEquationComponent.ArgumentOrdering"/>.Count) <br/>
        /// index: equation component index; the order of the equation
        /// components is defined by <see cref="m_AllComponentsOfMyType"/>;
        /// </summary>
        public int[] NoOfArguments;

        /// <summary>
        /// content: number of parameters for each equation component; 
        /// (<see cref="IEquationComponent.ParameterOrdering"/>.Count) <br/>
        /// index: equation component index; the order of the equation
        /// components is defined by <see cref="m_AllComponentsOfMyType"/>;
        /// </summary>
        public int[] NoOfParameters;

    }

    /// <summary>
    /// Extensions
    /// </summary>
    public static class EquationComponentArgMapping_Extensions {

        /// <summary>
        /// Utility function to determine whether it is required to evaluate some DG field resp its gradient.
        /// </summary>
        public static void DetermineReqFields<T>(this EquationComponentArgMapping<T>[] ecm, bool[] ValueRequired,  Func<T, bool> det) where T : IEquationComponent {
            int Gamma = ecm.Length;

            for (int gamma = 0; gamma < Gamma; gamma++) {
                for (int iComp = 0; iComp < ecm[gamma].m_AllComponentsOfMyType.Length; iComp++) {
                    int NoArgs = ecm[gamma].NoOfArguments[iComp];

                    T comp = ecm[gamma].m_AllComponentsOfMyType[iComp];

                    for (int na = 0; na < NoArgs; na++) {
                        int iDomField = ecm[gamma].AllToSub[iComp, na];

                        if (det(comp))
                            ValueRequired[iDomField] = true;
                    }
                }
            }
        }

    
        /// <summary>
        /// Creating performance stopwatches and registering them among <see cref="IQuadrature.CustomTimers"/> (and related).
        /// </summary>
        static public Stopwatch[][] InitStopWatches<T>(this EquationComponentArgMapping<T>[] fluxes, int FluxEvalPt, IQuadrature quad)
            where T : IEquationComponent //
        {
            Debug.Assert(quad.CustomTimers.Length == quad.CustomTimers_Names.Length);
            Debug.Assert(quad.CustomTimers.Length == quad.CustomTimers_RootPointer.Length);
            Debug.Assert(Array.IndexOf(quad.CustomTimers_Names, "Flux-Eval") == FluxEvalPt);

            Stopwatch[][] FluxesStopwatch = new Stopwatch[fluxes.Length][];
            List<string> names = new List<string>();
            List<Stopwatch> watches = new List<Stopwatch>();
            for (int i = 0; i < fluxes.Length; i++) {

                FluxesStopwatch[i] = new Stopwatch[fluxes[i].m_AllComponentsOfMyType.Length];
                int j = 0;
                foreach (object obj in fluxes[i].m_AllComponentsOfMyType) {
                    names.Add("Flux-Eval:" + obj.GetType().Name + ":" + obj.ToString());
                    FluxesStopwatch[i][j] = new Stopwatch();
                    watches.Add(FluxesStopwatch[i][j]);
                    j++;
                }
            }
            quad.CustomTimers = quad.CustomTimers.Cat(watches);
            quad.CustomTimers_Names = quad.CustomTimers_Names.Cat(names);
            int[] Pts = new int[names.Count];
            ArrayTools.SetAll(Pts, FluxEvalPt);
            quad.CustomTimers_RootPointer = quad.CustomTimers_RootPointer.Cat(Pts);

            return FluxesStopwatch;
        }

    }
}
