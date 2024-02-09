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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;
using ilPSP;
using BoSSS.Platform;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature.Linear;

namespace BoSSS.Foundation.Quadrature.NonLin {

    /// <summary>
    /// Common implementations for edge- and volume-quadrature of 
    /// 'N'onlinear 'E'quation 'C'omponents;
    /// </summary>
    abstract internal class NECQuadratureCommon {
        
        /// <summary>
        /// the differential operator
        /// </summary>
        readonly protected DifferentialOperator Operator;

        readonly protected IGridData m_GrdDat;

        /// <summary>
        /// the gird on which this quadrature object operates on.
        /// </summary>
        public IGridData GridDat {
            get {
                return m_GrdDat;
            }
        }


        

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="context">the context which Bjoern loves so much</param>
        /// <param name="DiffOp">the spatial operator</param>
        /// <param name="_DomainFields">
        /// the mapping for the DG fields (variables) in the domain of the differential operator <paramref name="DiffOp"/>;
        /// </param>
        /// <param name="CodomainMapping">
        /// the mapping for the DG fields (variables) in the codomain of the differential operator <paramref name="DiffOp"/>;
        /// </param>
        /// <param name="ParamFields">
        /// the mapping for the DG fields (variables) in the parameter list of the differential operator <paramref name="DiffOp"/>;
        /// </param>
        protected NECQuadratureCommon(IGridData context,
                                      DifferentialOperator DiffOp,
                                      IEnumerable<DGField> _DomainFields,
                                      IEnumerable<DGField> ParamFields,
                                      UnsetteledCoordinateMapping CodomainMapping) {
            // ---------------
            // check arguments
            // ---------------
            m_GrdDat = context;
            m_CodomainMapping = CodomainMapping;
            this.m_DomainFields = _DomainFields.ToArray();
            if (ParamFields != null && ParamFields.Count() > 0) {
                // concatenate parameters to domain mapping
                List<DGField> dom = _DomainFields.ToList(), param = ParamFields.ToList();
                DGField[] fld = new DGField[dom.Count + param.Count];
                int __i;
                for (__i = 0; __i < dom.Count; __i++)
                    fld[__i] = dom[__i];
                for (int j = 0; j < param.Count; j++)
                    fld[j + __i] = param[j];
                _DomainFields = fld;
            }

            m_CodomainBasisS = m_CodomainMapping.BasisS.ToArray();
            m_DomainAndParamFields = _DomainFields.ToArray();
            

            Operator = DiffOp;

            if ((DiffOp.DomainVar.Count + DiffOp.ParameterVar.Count) != m_DomainAndParamFields.Length) {
                string extMsg;
                extMsg = "[DiffOp domain and parameter vars: ";
                for (int ii = 0; ii < DiffOp.DomainVar.Count; ii++) {
                    extMsg += (DiffOp.DomainVar[ii]);
                    if(ii < DiffOp.DomainVar.Count - 1)
                        extMsg += ", ";
                }
                extMsg += "; ";
                for (int ii = 0; ii < DiffOp.ParameterVar.Count; ii++) {
                    extMsg += (DiffOp.ParameterVar[ii]);
                    if (ii < DiffOp.ParameterVar.Count - 1)
                        extMsg += ", ";
                }
                extMsg += "\n";
                extMsg += "Domain/Parameter mapping vars: ";
                for (int ii = 0; ii < m_DomainAndParamFields.Length; ii++) {
                    extMsg += (m_DomainAndParamFields[ii].Identification);
                    if (ii < m_DomainAndParamFields.Length - 1)
                        extMsg += ", ";
                }
                extMsg += "]";
                throw new ArgumentException("mismatch between number of domain variables: " + extMsg, "DomainMapping,DiffOp");
            }
            if (DiffOp.CodomainVar.Count != CodomainMapping.BasisS.Count)
                throw new ArgumentException("mismatch between number of codomain variables", "CodomainMapping,DiffOp");


            IList<Basis> CoDomBasisS = m_CodomainMapping.BasisS;
            m_NoOfTestFunctions = new int[CoDomBasisS.Count];
            m_MyMap = new int[CoDomBasisS.Count];

            int i = 0;
            int c = 0;
            foreach (Basis b in CoDomBasisS) {
                m_NoOfTestFunctions[i] = b.Length;
                m_MyMap[i] = c;
                c += m_NoOfTestFunctions[i];
                i++;
            }
            
            //m_MaxCodBasis = m_CodomainBasisS.ElementAtMax(bs => bs.Degree);
            //foreach(Basis bs in m_CodomainBasisS) {
            //    if(!bs.IsSubBasis(m_MaxCodBasis))
            //        throw new NotImplementedException();
            //}



        }


        protected abstract class ThreadLocals {

            readonly protected NECQuadratureCommon m_owner;
            protected int m_iThread;

            protected ThreadLocals(int iThread, NECQuadratureCommon _owner, IQuadrature q) {
                m_iThread = iThread;
                m_owner = _owner;


                // ------------------------
                // sort equation components
                // ------------------------
                m_NonlinFluxes = EquationComponentArgMapping<INonlinearFlux>.GetArgMapping(m_owner.Operator, true);
                m_NonlinFluxesEx = EquationComponentArgMapping<INonlinearFluxEx>.GetArgMapping(m_owner.Operator, true);

                if (q == null)
                    return;
            }



            /// <summary>
            /// array index: equation index
            /// </summary>
            internal EquationComponentArgMapping<INonlinearFlux>[] m_NonlinFluxes;

            /// <summary>
            /// array index: equation index
            /// </summary>
            internal EquationComponentArgMapping<INonlinearFluxEx>[] m_NonlinFluxesEx;

            /// <summary>
            /// true, if this integrator is responsible for any component
            /// </summary>
            public virtual bool IsNonEmpty {
                get {
                    return m_NonlinFluxes.IsNonEmpty() ||
                        m_NonlinFluxesEx.IsNonEmpty();
                }
            }

            /// <summary>
            /// array index: equation index
            /// </summary>
            protected Stopwatch[][] m_NonlinFluxesWatches;

            /// <summary>
            /// array index: equation index
            /// </summary>
            protected Stopwatch[][] m_NonlinFluxesExWatches;

        }



        /// <summary>
        /// <see cref="Time"/>
        /// </summary>
        protected double m_Time = Double.NaN;

        /// <summary>
        /// time value passed to flux functions
        /// </summary>
        public double Time {
            get {
                return m_Time;
            }
            set {
                m_Time = value;
            }
        }


        /// <summary>
        /// Index offset for each codomain variable
        /// </summary>
        readonly protected int[] m_MyMap;


        /// <summary>
        /// DG coordinate mapping for the codomain (output) of this
        /// quadrature object;
        /// index mapping for <see cref="m_Output"/>;
        /// </summary>
        readonly public UnsetteledCoordinateMapping m_CodomainMapping;

        /// <summary>
        /// The basis of the codomain variables - these are used as test functions;
        /// equal to the <see cref="BoSSS.Foundation.UnsetteledCoordinateMapping.BasisS"/>-member
        /// of <see cref="m_CodomainMapping"/>;
        /// </summary>
        readonly internal Basis[] m_CodomainBasisS;


        /// <summary>
        /// domain fields AND parameters (for the evaluation, there is no real difference between a domain field and a parameter),
        /// i.e., correlates with the concatenation of <see cref="IDifferentialOperator.DomainVar"/> and <see cref="IDifferentialOperator.ParameterVar"/>
        /// </summary>
        readonly internal DGField[] m_DomainAndParamFields;

        /// <summary>
        /// domain fields AND parameters (for the evaluation, there is no real difference between a domain field and a parameter),
        /// i.e., correlates with the concatenation of <see cref="IDifferentialOperator.DomainVar"/>
        /// </summary>
        readonly internal DGField[] m_DomainFields;

        /// <summary>
        /// Number of test functions;
        /// The test functions are the basis functions of the fields in the codomain;
        /// corresponds with <see cref="m_CodomainBasisS"/>;
        /// </summary>
        protected int[] m_NoOfTestFunctions;

        /// <summary>
        /// output array for the results of the integration;
        /// indices into it should be computed with <see cref="m_CodomainMapping"/>;
        /// The results of the integration/quadrature are multiplied by <see cref="m_alpha"/>
        /// and accumulated here;
        /// </summary>
        public IList<double> m_Output = null;

        /// <summary>
        /// Scaling for the result of the integration, see <see cref="m_Output"/>;
        /// </summary>
        public double m_alpha = 1.0;

        /// <summary>
        /// The performer of the quadrature.
        /// </summary>
        protected IQuadrature m_Quad;

        /// <summary>
        /// to be set by ctor of derivative class
        /// </summary>
        protected bool IsNonEmpty;

        /// <summary>
        /// executes the quadrature
        /// </summary>
        public void Execute() {
            if(IsNonEmpty == false)
                return;
            m_Quad.Execute();
        }

    }
}
