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

namespace BoSSS.Foundation.Quadrature.NonLin {

    /// <summary>
    /// Common implementations for edge- and volume-quadrature of 
    /// 'N'onlinear 'E'quation 'C'omponents;
    /// </summary>
    abstract internal class NECQuadratureCommon {
        
        /// <summary>
        /// the differential operator
        /// </summary>
        SpatialOperator m_DifferentialOperator;

        private IGridData m_GrdDat;

        /// <summary>
        /// the gird on which this quadrature object operates on.
        /// </summary>
        public IGridData GridDat {
            get {
                return m_GrdDat;
            }
        }

        internal static void DetermineReqFields<T>(bool[] ValueRequired, EquationComponentArgMapping<T>[] ecm, Func<T, bool> det) where T : IEquationComponent {
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
                                      SpatialOperator DiffOp,
                                      IList<DGField> _DomainFields,
                                      IList<DGField> ParamFields,
                                      UnsetteledCoordinateMapping CodomainMapping) {
            // ---------------
            // check arguments
            // ---------------
            m_GrdDat = context;
            m_CodomainMapping = CodomainMapping;
            if (ParamFields != null && ParamFields.Count > 0) {
                // concatenate parameters to domain mapping
                IList<DGField> dom = _DomainFields, param = ParamFields;
                DGField[] fld = new DGField[dom.Count + param.Count];
                int __i;
                for (__i = 0; __i < dom.Count; __i++)
                    fld[__i] = dom[__i];
                for (int j = 0; j < param.Count; j++)
                    fld[j + __i] = param[j];
                _DomainFields = fld;
            }

            m_CodomainBasisS = m_CodomainMapping.BasisS.ToArray();
            m_DomainFields = _DomainFields.ToArray();
            
            _DomainFields = null;

            m_DifferentialOperator = DiffOp;

            if ((DiffOp.DomainVar.Count + DiffOp.ParameterVar.Count) != m_DomainFields.Length) {
                string extMsg;
                extMsg = "[DiffOp domain and parameter vars: ";
                for (int ii = 0; ii < DiffOp.DomainVar.Count; ii++)
                    extMsg += (DiffOp.DomainVar[ii] + ", ");
                extMsg += "; ";
                for (int ii = 0; ii < DiffOp.ParameterVar.Count; ii++)
                    extMsg += (DiffOp.ParameterVar[ii] + ", ");
                extMsg += "\n";
                extMsg += "Domain/Parameter mapping vars: ";
                for (int ii = 0; ii < m_DomainFields.Length; ii++)
                    extMsg += (m_DomainFields[ii].Identification + ", ");
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


            // ------------------------
            // sort equation components
            // ------------------------
            m_NonlinFluxes = DiffOp.GetArgMapping<INonlinearFlux>(true);
            m_NonlinFluxesEx = DiffOp.GetArgMapping<INonlinearFluxEx>(true);
        }

        


        /// <summary>
        /// array index: equation index
        /// </summary>
        protected EquationComponentArgMapping<INonlinearFlux>[] m_NonlinFluxes;

        /// <summary>
        /// array index: equation index
        /// </summary>
        protected EquationComponentArgMapping<INonlinearFluxEx>[] m_NonlinFluxesEx;
        
        /// <summary>
        /// array index: equation index
        /// </summary>
        protected Stopwatch[][] m_NonlinFluxesWatches;

        /// <summary>
        /// array index: equation index
        /// </summary>
        protected Stopwatch[][] m_NonlinFluxesExWatches;


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
        /// used by <see cref="MyMap"/>;
        /// </summary>
        int[] m_MyMap;

        /// <summary>
        /// internal index mapping;
        /// </summary>
        /// <param name="variableIndex">internal field index; <see cref="m_CodomainBasisS"/>;</param>
        /// <param name="CoordInd">coordinate index, basis function index;</param>
        /// <returns></returns>
        protected int MyMap(int variableIndex, int CoordInd) {
            return m_MyMap[variableIndex] + CoordInd;
        }

        /// <summary>
        /// DG coordinate mapping for the codomain (output) of this
        /// quadrature object;
        /// index mapping for <see cref="m_Output"/>;
        /// </summary>
        public UnsetteledCoordinateMapping m_CodomainMapping;

        /// <summary>
        /// The basis of the codomain variables - these are used as test functions;
        /// equal to the <see cref="BoSSS.Foundation.UnsetteledCoordinateMapping.BasisS"/>-member
        /// of <see cref="m_CodomainMapping"/>;
        /// </summary>
        protected Basis[] m_CodomainBasisS;

        ///// <summary>
        ///// Largest basis, which includes all other basises found in <see cref="m_CodomainBasisS"/>.
        ///// </summary>
        //protected Basis m_MaxCodBasis;

        /// <summary>
        /// equal to the <see cref="BoSSS.Foundation.CoordinateMapping.Fields"/>-member
        /// of <see cref="m_CodomainMapping"/>
        /// </summary>
        protected DGField[] m_DomainFields;

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
        /// executes the quadrature
        /// </summary>
        public void Execute() {
            m_Quad.Execute();
        }

    }
}
