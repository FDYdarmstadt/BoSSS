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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.XDG {
    
    /// <summary>
    /// integrator for coupling components, linear components;
    /// </summary>
    internal class LECQuadratureLevelSet<M,V> : BoSSS.Foundation.Quadrature.CellQuadrature  
        where M : IMutableMatrix
        where V : IList<double> //
    {
        
        UnsetteledCoordinateMapping m_RowMap;
        UnsetteledCoordinateMapping m_ColMap;

        DGField[] m_Parameters;

        bool AffineOnly {
            get {
                return this.OperatorMatrix == null;
            }
        }

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// index of the level set to evaluate
        /// </summary>
        int m_LevSetIdx;

        LevelSetTracker m_lsTrk;

        /// <summary>
        /// Negative and positive (with respect to level-set) species.
        /// </summary>
        Tuple<SpeciesId,SpeciesId> m_SpeciesPair;

        /// <summary>
        /// Negative species/Species A
        /// </summary>
        SpeciesId SpeciesA {
            get {
                return m_SpeciesPair.Item1;
            }
        }

        /// <summary>
        /// Positive species/Species B
        /// </summary>
        SpeciesId SpeciesB {
            get {
                return m_SpeciesPair.Item2;
            }
        }

        

                
        /// <summary>
        /// ctor.
        /// </summary>
        internal LECQuadratureLevelSet(IGridData context,
                                     XSpatialOperator DiffOp,
                                     M Matrix, V OffsetVec,
                                     UnsetteledCoordinateMapping RowMap, IList<DGField> ParamsMap, UnsetteledCoordinateMapping ColMap,
                                     LevelSetTracker lsTrk, int _iLevSet, 
                                     Tuple<SpeciesId,SpeciesId> SpeciesPair,
                                     //Tuple<CoefficientSet,CoefficientSet> NegPosCoeff,
                                     ICompositeQuadRule<QuadRule> domAndRule) //
            : base(new int[] { RowMap.MaxTotalNoOfCoordinatesPerCell, 1 + ((Matrix == null) ? 0 : ColMap.MaxTotalNoOfCoordinatesPerCell) },
                 context,
                 domAndRule) //
        {

            // ------------------
            // init custom timers
            // ------------------

            base.CustomTimers = new Stopwatch[] { new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch() };
            base.CustomTimers_Names = new string[] { "Flux-Eval", "Basis-Eval", "Loops", "ParametersAndNormals" };
            base.CustomTimers_RootPointer = new int[4];
            ArrayTools.SetAll(base.CustomTimers_RootPointer, -1);

            // -----------------------------------
            // set members / check ctor parameters
            // -----------------------------------

            m_RowMap = RowMap;
            m_ColMap = ColMap;
            m_Parameters = (ParamsMap != null) ? ParamsMap.ToArray() : new DGField[0];
                       
            if (m_RowMap.BasisS.Count != DiffOp.CodomainVar.Count)
                throw new ArgumentException("mismatch between number of codomain variables in spatial operator and row mapping");
            if (m_ColMap.BasisS.Count != DiffOp.DomainVar.Count)
                throw new ArgumentException("mismatch between number of domain variables in spatial operator and column mapping");
            if (m_Parameters.Length != DiffOp.ParameterVar.Count)
                throw new ArgumentException("mismatch between number of parameter variables in spatial operator and given parameters");

            int Gamma = m_RowMap.BasisS.Count;
            m_lsTrk = lsTrk;

            if (Matrix != null && (Matrix.RowPartitioning.LocalLength != RowMap.LocalLength))
                throw new ArgumentException("mismatch between matrix number of rows and row mapping.");
            if (Matrix != null && (Matrix.ColPartition.LocalLength != ColMap.LocalLength))
                throw new ArgumentException("mismatch between matrix number of columns and column mapping.");

            this.m_LevSetIdx = _iLevSet;
            this.m_SpeciesPair = SpeciesPair;

            this.OperatorMatrix = Matrix;
            this.OperatorAffine = OffsetVec;

            // ------------------------
            // sort equation components
            // ------------------------

            //m_BiLinForms = DiffOp.GetArgMapping<IBilinearForm>(true, Compfilter<IBilinearForm>);
            //m_2ndDerivFlux = DiffOp.GetArgMapping<ILinear2ndDerivativeCouplingFlux>(true, Compfilter<ILinear2ndDerivativeCouplingFlux>);

            m_LsForm_UxV = DiffOp.GetArgMapping<ILevelSetForm_UxV>(true,
               eq => ((eq.LevelSetTerms & TermActivationFlags.UxV) != 0) && Compfilter(eq),
               eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_GradUxV = DiffOp.GetArgMapping<ILevelSetForm_GradUxV>(true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.GradUxV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_UxGradV = DiffOp.GetArgMapping<ILevelSetForm_UxGradV>(true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.UxGradV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_GradUxGradV = DiffOp.GetArgMapping<ILevelSetForm_GradUxGradV>(true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.GradUxGradV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_V = DiffOp.GetArgMapping<ILevelSetForm_V>(true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.V) != 0 && Compfilter(eq)),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_GradV = DiffOp.GetArgMapping<ILevelSetForm_GradV>(true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.GradV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);

            this.m_LsForm_UxV_Watches = this.m_LsForm_UxV.InitStopWatches(0, this);
            this.m_LsForm_GradUxV_Watches = this.m_LsForm_GradUxV.InitStopWatches(0, this);
            this.m_LsForm_UxGradV_Watches = this.m_LsForm_UxGradV.InitStopWatches(0, this);
            this.m_LsForm_GradUxGradV_Watches = this.m_LsForm_GradUxGradV.InitStopWatches(0, this);
            this.m_LsForm_V_Watches = this.m_LsForm_V.InitStopWatches(0, this);
            this.m_LsForm_GradV_Watches = this.m_LsForm_GradV.InitStopWatches(0, this);

            // -----
            // alloc
            // -----

            AllocEmpty(m_LsForm_UxV, out Koeff_UxV, out Sum_Koeff_UxV, 5, false);
            AllocEmpty(m_LsForm_GradUxV, out Koeff_NablaUxV, out Sum_Koeff_NablaUxV, 6, false);
            AllocEmpty(m_LsForm_UxGradV, out Koeff_UxNablaV, out Sum_Koeff_UxNablaV, 6, false);
            AllocEmpty(m_LsForm_GradUxGradV, out Koeff_NablaUxNablaV, out Sum_Koeff_NablaUxNablaV, 7, false);

            AllocEmpty(m_LsForm_V, out Koeff_V, out Sum_Koeff_V, 3, true);
            AllocEmpty(m_LsForm_GradV, out Koeff_NablaV, out Sum_Koeff_NablaV, 4, true);
        }




        private bool Compfilter(IEquationComponent c) {

            ILevelSetForm b = (ILevelSetForm)c;
            
            if (this.m_LevSetIdx != b.LevelSetIndex)
                // component is not relevant for this level-set
                return false;

            if ((this.SpeciesA == b.NegativeSpecies && this.SpeciesB == b.PositiveSpecies)
                || (this.SpeciesB == b.NegativeSpecies && this.SpeciesA == b.PositiveSpecies)) {

            } else {
                // component is not relevant for this level-set
                return false;
            }

            // filter passed
            return true;
        }
        
        void AllocEmpty<T>(EquationComponentArgMapping<T>[] components,
            out MultidimensionalArray[][] compBuffer, out MultidimensionalArray[,] sumBuffer, int Dim, bool affine)
            where T : IEquationComponent {

            int Gamma = this.m_RowMap.BasisS.Count; // no of codom vars/test functions
            if (Gamma != components.Length)
                throw new ApplicationException("internal error");

            int Delta = this.m_ColMap.BasisS.Count; // no of domain vars/trial functions

            sumBuffer = new MultidimensionalArray[Gamma, affine ? 1 : Delta];

            compBuffer = new MultidimensionalArray[Gamma][];
            for (int gamma = 0; gamma < Gamma; gamma++) { // loop over codom vars/test functions...
                var bg = components[gamma];
                int NoComp = bg.m_AllComponentsOfMyType.Length;

                compBuffer[gamma] = new MultidimensionalArray[NoComp];

                bool[] UsedMarker = new bool[Delta];

                int iComp = -1;
                foreach (var b in bg.m_AllComponentsOfMyType) { // loop over equation components...
                    iComp++;
                    compBuffer[gamma][iComp] = new MultidimensionalArray(Dim);
                    //usedAtAll = true;

                    int NoOfArgs = b.ArgumentOrdering.Count;
                    if (affine) {
                        UsedMarker[0] = true;
                    } else {
                        for (int j = 0; j < NoOfArgs; j++) {
                            int targ = bg.AllToSub[iComp, j];
                            UsedMarker[targ] = true;
                        }
                    }
                }

                for(int delta = 0; delta < (affine ? 1 : Delta); delta++) {
                    if (UsedMarker[delta])
                        sumBuffer[gamma, delta] = new MultidimensionalArray(Dim - (affine ? 0 : 1));
                }
            }
        }

        

        EquationComponentArgMapping<ILevelSetForm_UxV>[] m_LsForm_UxV;
        EquationComponentArgMapping<ILevelSetForm_GradUxV>[] m_LsForm_GradUxV;
        EquationComponentArgMapping<ILevelSetForm_UxGradV>[] m_LsForm_UxGradV;
        EquationComponentArgMapping<ILevelSetForm_GradUxGradV>[] m_LsForm_GradUxGradV;
        EquationComponentArgMapping<ILevelSetForm_V>[] m_LsForm_V;
        EquationComponentArgMapping<ILevelSetForm_GradV>[] m_LsForm_GradV;

        Stopwatch[][] m_LsForm_UxV_Watches;
        Stopwatch[][] m_LsForm_GradUxV_Watches;
        Stopwatch[][] m_LsForm_UxGradV_Watches;
        Stopwatch[][] m_LsForm_GradUxGradV_Watches;
        Stopwatch[][] m_LsForm_V_Watches;
        Stopwatch[][] m_LsForm_GradV_Watches;

        /// <summary>
        /// values of parameter fields, positive side of level Set
        /// <list type="bullet">
        ///   <item>1st index: parameter index</item>
        ///   <item>2nd index: cell index, plus some offset</item>
        ///   <item>3rd index: node index</item>
        /// </list>
        /// </summary>
        MultidimensionalArray m_ParamFieldValuesPos;
        
        /// <summary>
        /// values of parameter fields, positive side of level Set
        /// <list type="bullet">
        ///   <item>1st index: parameter index</item>
        ///   <item>2nd index: cell index, plus some offset</item>
        ///   <item>3rd index: node index</item>
        /// </list>
        /// </summary>
        MultidimensionalArray m_ParamFieldValuesNeg;

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // buffers for 'IBilinearForm'
        // 1st index: codomain variable/test function
        // 2nd index: equation component
        MultidimensionalArray[][] Koeff_UxV;
        MultidimensionalArray[][] Koeff_NablaUxV;
        MultidimensionalArray[][] Koeff_UxNablaV;
        MultidimensionalArray[][] Koeff_NablaUxNablaV;
        MultidimensionalArray[][] Koeff_V;
        MultidimensionalArray[][] Koeff_NablaV;
        // ----------------
        // 1st index: codomain
        // 2nd index: domain 
        MultidimensionalArray[,] Sum_Koeff_UxV;
        MultidimensionalArray[,] Sum_Koeff_NablaUxV;
        MultidimensionalArray[,] Sum_Koeff_UxNablaV;
        MultidimensionalArray[,] Sum_Koeff_NablaUxNablaV;
        MultidimensionalArray[,] Sum_Koeff_V; // second index is always 0
        MultidimensionalArray[,] Sum_Koeff_NablaV;  // second index is always 0
        // ------------------------------------------------------------------


        static T[] ToArray<T>(params T[] obj) {
            return obj;
        }
        
        /// <summary>
        /// memory allocation.
        /// </summary>
        protected override void AllocateBuffers(int Nitm, NodeSet ruleNodes) {
            int Nnod = ruleNodes.GetLength(0);
            base.AllocateBuffers(Nitm, ruleNodes);
            if (m_Parameters.Length > 0) {
                m_ParamFieldValuesPos = MultidimensionalArray.Create(m_Parameters.Length, Nitm, Nnod);
                m_ParamFieldValuesNeg = MultidimensionalArray.Create(m_Parameters.Length, Nitm, Nnod);
            }
            int D = gridData.SpatialDimension;

            int _Delta = m_ColMap.BasisS.Count;

            Alloc(m_LsForm_UxV, Koeff_UxV, Sum_Koeff_UxV, 2, Nitm, Nnod, _Delta, 2, 2);
            Alloc(m_LsForm_GradUxV, Koeff_NablaUxV, Sum_Koeff_NablaUxV, 2, Nitm, Nnod, _Delta, 2, 2, D);
            Alloc(m_LsForm_UxGradV, Koeff_UxNablaV, Sum_Koeff_UxNablaV, 2, Nitm, Nnod, _Delta, 2, 2, D);
            Alloc(m_LsForm_GradUxGradV, Koeff_NablaUxNablaV, Sum_Koeff_NablaUxNablaV, 2, Nitm, Nnod, _Delta, 2, 2, D, D);

            Alloc(m_LsForm_V, Koeff_V, Sum_Koeff_V, -1, Nitm, Nnod, 2);
            Alloc(m_LsForm_GradV, Koeff_NablaV, Sum_Koeff_NablaV, -1, Nitm, Nnod, 2, D);
        }

        void Alloc<T>(EquationComponentArgMapping<T>[] components,
            MultidimensionalArray[][] compBuffer, MultidimensionalArray[,] sumBuffer,
            int i_of_Delta,
            params int[] Lengths) where T : IEquationComponent {

            int Gamma = this.m_RowMap.BasisS.Count;
            Debug.Assert(Gamma == components.Length);
            int Delta = this.m_ColMap.BasisS.Count;
            if (i_of_Delta >= 0)
                Debug.Assert(Lengths[i_of_Delta] == Delta);

            int[] SumbufLengths;
            if (i_of_Delta >= 0) {
                SumbufLengths = new int[Lengths.Length - 1];
                int ii = 0;
                for (int i = 0; i < Lengths.Length; i++) {
                    if (i == i_of_Delta)
                        continue;
                    SumbufLengths[ii] = Lengths[i];
                    ii++;
                }
            } else {
                SumbufLengths = Lengths;
            }

            for (int gamma = 0; gamma < Gamma; gamma++) {
                var bg = components[gamma];
                
                int iComp = -1;
                foreach (var b in bg.m_AllComponentsOfMyType) {
                    iComp++;

                    if (i_of_Delta >= 0)
                        Lengths[i_of_Delta] = b.ArgumentOrdering.Count;

                    compBuffer[gamma][iComp].Allocate(Lengths);

                    Debug.Assert(sumBuffer != null, "oops!");
                }

                for(int delta = 0; delta < (i_of_Delta >= 0 ? Delta : 1); delta++) {
                    if (sumBuffer[gamma, delta] != null)
                        sumBuffer[gamma, delta].Allocate(SumbufLengths);
                }
            }
        }

        protected override void Evaluate(int i0, int Len, QuadRule QR, MultidimensionalArray EvalResult) {
            NodeSet QuadNodes = QR.Nodes;
            int D = gridData.SpatialDimension;
            int NoOfNodes = QuadNodes.NoOfNodes;
            int DELTA = m_ColMap.BasisS.Count; // DELTA: number of basis variables, number of domain variables
            int GAMMA = m_RowMap.BasisS.Count;  // GAMMA: number of codom variables


            // Evaluate Parameter fields
            // -------------------------
            base.CustomTimers[3].Start();
            for (int i = 0; i < m_Parameters.Length; i++) {
                var bufPos = m_ParamFieldValuesPos.ExtractSubArrayShallow(i, -1, -1);
                var bufNeg = m_ParamFieldValuesNeg.ExtractSubArrayShallow(i, -1, -1);
                DGField m_Field = m_Parameters[i];
                if (m_Field != null) {
                    if (m_Field is XDGField) {
                        // jump in parameter i at level-set: separate evaluation for both sides
                        var xfi = m_Field as XDGField;

                        //xfi.GetSpeciesShadowField(this.SpeciesA).Evaluate(i0, Len, QuadNodes, bufNeg);
                        //xfi.GetSpeciesShadowField(this.SpeciesB).Evaluate(i0, Len, QuadNodes, bufPos);
                        xfi.GetSpeciesShadowField(m_lsTrk.GetSpeciesIdFromSign(-1)).Evaluate(i0, Len, QuadNodes, bufNeg);
                        xfi.GetSpeciesShadowField(m_lsTrk.GetSpeciesIdFromSign(+1)).Evaluate(i0, Len, QuadNodes, bufPos);

                    } else {
                        // no jump at level set: positive and negative limit of parameter i are equal
                        m_Parameters[i].Evaluate(i0, Len, QuadNodes, bufPos);
                        bufNeg.Set(bufPos);
                    }
                } else {
                    bufPos.Clear();
                    bufNeg.Clear();
                }
            }

            // Evaluate level sets and normals
            // -------------------------------
            var NoOfLevSets = m_lsTrk.LevelSets.Count;
            MultidimensionalArray Normals = m_lsTrk.DataHistories[m_LevSetIdx].Current.GetLevelSetNormals(QuadNodes, i0, Len);
            base.CustomTimers[3].Stop();

            // Evaluate basis and test functions
            // ---------------------------------

            bool[] ReqU = new bool[DELTA];
            bool[] ReqGradU = new bool[DELTA];
            bool[] ReqV = new bool[GAMMA];
            bool[] ReqGradV = new bool[GAMMA];

            for(int gamma = 0; gamma < GAMMA; gamma++) {
                for(int delta = 0; delta < DELTA; delta++) {
                    if (Sum_Koeff_UxV[gamma, delta] != null) {
                        ReqU[delta] = true;
                        ReqV[gamma] = true;
                    }
                    if (Sum_Koeff_NablaUxV[gamma, delta] != null) {
                        ReqGradU[delta] = true;
                        ReqV[gamma] = true;
                    }
                    if (Sum_Koeff_UxNablaV[gamma, delta] != null) {
                        ReqU[delta] = true;
                        ReqGradV[gamma] = true;
                    }
                    if (Sum_Koeff_NablaUxNablaV[gamma, delta] != null) {
                        ReqGradU[delta] = true;
                        ReqGradV[gamma] = true;
                    }
                }
                if (Sum_Koeff_V[gamma, 0] != null) {
                    ReqV[gamma] = true;
                }
                if (Sum_Koeff_NablaV[gamma, 0] != null) {
                    ReqGradV[gamma] = true;
                    //ReqV[gamma] = true;
                }
            }
            
            MultidimensionalArray[] BasisValues; //           index: domain variable/trial function
            MultidimensionalArray[] BasisGradientValues; //   index: domain variable/trial function
            MultidimensionalArray[] TestValues; //            index: codom variable/test function
            MultidimensionalArray[] TestGradientValues; //    index: codom variable/test function
            int[,] sectionsBasis;
            int[,] sectionsTest;
            base.CustomTimers[1].Start();
            EvalBasis(i0, Len, this.m_RowMap.BasisS, ReqV, ReqGradV, out TestValues, out TestGradientValues, out sectionsTest, QuadNodes);
            EvalBasis(i0, Len, this.m_ColMap.BasisS, ReqU, ReqGradU, out BasisValues, out BasisGradientValues, out sectionsBasis, QuadNodes);
            base.CustomTimers[1].Stop();

            // compute offsets into matrix blocks
            // ----------------------------------

            
            int[] offsetDom = new int[DELTA];
            CompOffsets(i0, Len, offsetDom, m_ColMap);

            int[] offsetCod = new int[GAMMA];
            CompOffsets(i0, Len, offsetCod, m_RowMap);


            // Transform nodes to global coordinates
            // -------------------------------------

            MultidimensionalArray NodesGlobal = gridData.GlobalNodes.GetValue_Cell(QuadNodes, i0, Len);

            // Evaluate Integral components
            // ----------------------------

            var _inParams = new BoSSS.Foundation.XDG.LevSetIntParams();
            _inParams.i0 = i0;

            
            
            // loop over codomain variables ...
            for (int gamma = 0; gamma < GAMMA; gamma++) {

                // prepare parameters
                // - - - - - - - - - 

                // set Normal's
                _inParams.Normal = Normals;
                // set Nodes Global
                _inParams.X = NodesGlobal;
                _inParams.time = this.time;
                _inParams.LsTrk = this.m_lsTrk;
                // set length scales


                // clear summation buffers
                // - - - - - - - - - - - -
                for(int delta = 0; delta < DELTA; delta++) {
                    if (Sum_Koeff_UxV[gamma, delta] != null)
                        Sum_Koeff_UxV[gamma, delta].Clear();
                    if (Sum_Koeff_NablaUxV[gamma, delta] != null)
                        Sum_Koeff_NablaUxV[gamma, delta].Clear();
                    if (Sum_Koeff_UxNablaV[gamma, delta] != null)
                        Sum_Koeff_UxNablaV[gamma, delta].Clear();
                    if (Sum_Koeff_NablaUxNablaV[gamma, delta] != null)
                        Sum_Koeff_NablaUxNablaV[gamma, delta].Clear();
                }
                if (Sum_Koeff_NablaV[gamma, 0] != null)
                    Sum_Koeff_NablaV[gamma, 0].Clear();
                if (Sum_Koeff_V[gamma, 0] != null)
                    Sum_Koeff_V[gamma, 0].Clear();


                // Evaluate Bilin. forms
                // - - - - - - - - - - -
                {
                    EvalComponent(_inParams, gamma, this.m_LsForm_UxV[gamma], this.m_LsForm_UxV_Watches[gamma], 
                        Koeff_UxV, Sum_Koeff_UxV , 2,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_UxV _comp, int _gamma, int i, LevSetIntParams inp) {
                            _comp.LevelSetForm_UxV(_inParams, Koeff_UxV[_gamma][i]);
                        });
                }
                {
                    EvalComponent(_inParams, gamma, this.m_LsForm_GradUxV[gamma], this.m_LsForm_GradUxV_Watches[gamma],
                        Koeff_NablaUxV, Sum_Koeff_NablaUxV, 2,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_GradUxV _comp, int _gamma, int i, LevSetIntParams inp) {
                            _comp.LevelSetForm_GradUxV(_inParams, Koeff_NablaUxV[_gamma][i]);
                        });
                }
                {
                    EvalComponent(_inParams, gamma, this.m_LsForm_UxGradV[gamma], this.m_LsForm_UxGradV_Watches[gamma],
                        Koeff_UxNablaV, Sum_Koeff_UxNablaV, 2,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_UxGradV _comp, int _gamma, int i, LevSetIntParams inp) {
                            _comp.LevelSetForm_UxGradV(_inParams, Koeff_UxNablaV[_gamma][i]);
                        });
                }
                {
                    EvalComponent(_inParams, gamma, this.m_LsForm_GradUxGradV[gamma], this.m_LsForm_GradUxGradV_Watches[gamma],
                        Koeff_NablaUxNablaV, Sum_Koeff_NablaUxNablaV, 2,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_GradUxGradV _comp, int _gamma, int i, LevSetIntParams inp) {
                            _comp.LevelSetForm_GradUxGradV(_inParams, Koeff_NablaUxNablaV[_gamma][i]);
                        });
                }

                {
                    EvalComponent(_inParams, gamma, this.m_LsForm_V[gamma], this.m_LsForm_V_Watches[gamma],
                        Koeff_V, Sum_Koeff_V, -1,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_V _comp, int _gamma, int i, LevSetIntParams inp) {
                            _comp.LevelSetForm_V(_inParams, Koeff_V[_gamma][i]);
                        });
                }
                {
                    EvalComponent(_inParams, gamma, this.m_LsForm_GradV[gamma], this.m_LsForm_GradV_Watches[gamma],
                        Koeff_NablaV, Sum_Koeff_NablaV, -1,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_GradV _comp, int _gamma, int i, LevSetIntParams inp) {
                            _comp.LevelSetForm_GradV(_inParams, Koeff_NablaV[_gamma][i]);
                        });
                }
            }
            

            // Summation Loops: multiply with test and trial functions
            // -------------------------------------------------------

            for (int gamma = 0; gamma < GAMMA; gamma++) {
                // Evaluate Integrand
                // - - - - - - - - - 

                var TestVal = TestValues[gamma];
                var TestGradVal = TestGradientValues[gamma];
                int N;
                if (TestVal != null)
                    N = TestVal.GetLength(2);
                else if (TestGradVal != null)
                    N = TestGradVal.GetLength(2);
                else
                    N = 0;

                base.CustomTimers[2].Start();

                // loop over domain variables/basis variables ...
                for (int delta = 0; (delta < DELTA) && !AffineOnly; delta++) {

                    var BasisVal = BasisValues[delta];
                    var BasisGradVal = BasisGradientValues[delta];
                    int M;
                    if (BasisVal != null)
                        M = BasisVal.GetLength(2);
                    else if (BasisGradVal != null)
                        M = BasisGradVal.GetLength(2);
                    else
                        M = 0;


                    // Das Schleifenmonster: ---------------------------------------------
                    for (int cr = 0; cr < 2; cr++) {  // 
                        for (int cc = 0; cc < 2; cc++) {

                            int[] extr0 = new int[] { 0, 0, sectionsTest[gamma, cr] * N + offsetCod[gamma], sectionsBasis[delta, cc] * M + 1 + offsetDom[delta] };
                            int[] extrE = new int[] { Len - 1, NoOfNodes - 1, extr0[2] + N - 1, extr0[3] + M - 1 };
                            var SubRes = EvalResult.ExtractSubArrayShallow(extr0, extrE);

                            if (Sum_Koeff_UxV[gamma, delta] != null) {
                                var Sum_Koeff_UxV_CrCc = Sum_Koeff_UxV[gamma, delta].ExtractSubArrayShallow(-1, -1, cr, cc);
                                SubRes.Multiply(1.0, Sum_Koeff_UxV_CrCc, BasisVal, TestVal, 1.0, "jknm", "jk", "jkm", "jkn");
                            }

                            if (Sum_Koeff_NablaUxV[gamma, delta] != null) {
                                var Sum_Koeff_NablaUxV_CrCc = Sum_Koeff_NablaUxV[gamma, delta].ExtractSubArrayShallow(-1, -1, cr, cc, -1);
                                SubRes.Multiply(1.0, Sum_Koeff_NablaUxV_CrCc, BasisGradVal, TestVal, 1.0, "jknm", "jkd", "jkmd", "jkn");
                            }
                            if (Sum_Koeff_UxNablaV[gamma, delta] != null) {
                                var Sum_Koeff_UxNablaV_CrCc = Sum_Koeff_UxNablaV[gamma, delta].ExtractSubArrayShallow(-1, -1, cr, cc, -1);
                                SubRes.Multiply(1.0, Sum_Koeff_UxNablaV_CrCc, BasisVal, TestGradVal, 1.0, "jknm", "jkd", "jkm", "jknd");
                            }
                            if (Sum_Koeff_NablaUxNablaV[gamma, delta] != null) {
                                var Sum_Koeff_NablaUxNablaV_CrCc = Sum_Koeff_NablaUxNablaV[gamma, delta].ExtractSubArrayShallow(-1, -1, delta, cr, cc, -1, -1);
                                SubRes.Multiply(1.0, Sum_Koeff_NablaUxNablaV_CrCc, BasisGradVal, TestGradVal, 1.0, "jknm", "jkde", "jkmd", "jkne");
                            }
                        }
                    }
                    // --------------------------------------- 

                } // end loop delta

                // affine offset
                {

                    // Das Schleifenmonster: ---------------------------------------------
                    for (int cr = 0; cr < 2; cr++) { // 
                        int[] extr0 = new int[] { 0, 0, sectionsTest[gamma, cr] * N + offsetCod[gamma], 0 };
                        int[] extrE = new int[] { Len - 1, NoOfNodes - 1, extr0[2] + N - 1, -1 };
                        var SubRes = EvalResult.ExtractSubArrayShallow(extr0, extrE);

                        if (Sum_Koeff_V[gamma, 0] != null) {
                            var Sum_Koeff_V_Cr = Sum_Koeff_V[gamma, 0].ExtractSubArrayShallow(-1, -1, cr);
                            SubRes.Multiply(1.0, Sum_Koeff_V_Cr, TestVal, 1.0, "jkn", "jk", "jkn");
                        }
                        if (Sum_Koeff_NablaV[gamma, 0] != null) {
                            var Sum_Koeff_NablaV_Cr = Sum_Koeff_NablaV[gamma, 0].ExtractSubArrayShallow(-1, -1, cr, -1);
                            SubRes.Multiply(1.0, Sum_Koeff_NablaV_Cr, TestGradVal, 1.0, "jkn", "jkd", "jknd");
                        }
                    }
                    // ---------------------------------------
                }

                base.CustomTimers[2].Stop();
            }
        }

        static internal void CompOffsets(int i0, int L, int[] offset, UnsetteledCoordinateMapping Map) {
            int DELTA = offset.Length;
            Debug.Assert(DELTA == Map.BasisS.Count);
            
            int _o0 = Map.LocalUniqueCoordinateIndex(0, i0, 0);
            for (int delta = 0; delta < offset.Length; delta++)
                offset[delta] = Map.LocalUniqueCoordinateIndex(delta, i0, 0) - _o0;

#if DEBUG
            int[] Ns = new int[DELTA];
            for (int delta = 0; delta < DELTA; delta++)
                Ns[delta] = Map.BasisS[delta].GetLength(i0);
            for (int i = 1; i < L; i++) {
                for (int delta = 0; delta < DELTA; delta++)
                    Debug.Assert(Ns[delta] == Map.BasisS[delta].GetLength(i0+i));
                
            }
#endif
        }

        delegate void CallComponent<T>(T comp, int gamma, int i, LevSetIntParams _inParams);

        static private void EvalComponent<T>(LevSetIntParams _inParams,
            int gamma, EquationComponentArgMapping<T> bf, Stopwatch[] timers,
            MultidimensionalArray[][] argsPerComp, MultidimensionalArray[,] argsSum,
            int componentIdx,
            MultidimensionalArray ParamFieldValuesPos, MultidimensionalArray ParamFieldValuesNeg,
            int DELTA,
            Stopwatch timer,
            CallComponent<T> ComponentFunc) where T : ILevelSetForm {
            timer.Start();
            
            
            for (int i = 0; i < bf.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components
                var comp = bf.m_AllComponentsOfMyType[i];

                //LengthScales.TryGetValue(comp.NegativeSpecies, out _inParams.NegCellLengthScale);
                //LengthScales.TryGetValue(comp.PositiveSpecies, out _inParams.PosCellLengthScale);
                
                argsPerComp[gamma][i].Clear();

                int NoOfArgs = bf.NoOfArguments[i];
                Debug.Assert(NoOfArgs == comp.ArgumentOrdering.Count);
                int NoOfParams = bf.NoOfParameters[i];
                Debug.Assert(NoOfParams == ((comp.ParameterOrdering != null) ? comp.ParameterOrdering.Count : 0));


                // map parameters
                _inParams.ParamsPos = new MultidimensionalArray[NoOfParams];
                _inParams.ParamsNeg = new MultidimensionalArray[NoOfParams];
                for (int c = 0; c < NoOfParams; c++) {
                    int targ = bf.AllToSub[i, c + NoOfArgs] - DELTA;
                    Debug.Assert(targ >= 0);
                    _inParams.ParamsPos[c] = ParamFieldValuesPos.ExtractSubArrayShallow(targ, -1, -1);
                    _inParams.ParamsNeg[c] = ParamFieldValuesNeg.ExtractSubArrayShallow(targ, -1, -1);
                }

                // evaluate equation components
                timers[i].Start();
                ComponentFunc(comp, gamma, i, _inParams);
                timers[i].Stop();
#if DEBUG
                argsPerComp[gamma][i].CheckForNanOrInf();
#endif

                // sum up bilinear forms:
                {
                    MultidimensionalArray Summand = argsPerComp[gamma][i];
                    
                    if (componentIdx >= 0) {
                        for (int c = 0; c < NoOfArgs; c++) { // loop over arguments of equation component
                            int targ = bf.AllToSub[i, c];
                            int[] selSum = new int[Summand.Dimension];
                            selSum.SetAll(-1);
                            selSum[componentIdx] = c;

                            MultidimensionalArray Accu = argsSum[gamma, targ];

                            //int[] selAccu = new int[Accu.Dimension];
                            //selAccu.SetAll(-1);
                            //selAccu[componentIdx] = targ;
#if DEBUG
                            Summand.ExtractSubArrayShallow(selSum).CheckForNanOrInf();
#endif

                            Accu.Acc(1.0, Summand.ExtractSubArrayShallow(selSum));
                        }
                    } else {
                        // affin
#if DEBUG
                        Summand.CheckForNanOrInf();
#endif
                        MultidimensionalArray Accu = argsSum[gamma, 0];
                        Accu.Acc(1.0, Summand);
                    }
                    
                }
                

            }
            timer.Stop();
        }

        /// <summary>
        /// Evaluation of DG basis and DG basis gradients
        /// </summary>
        /// <param name="i0">1st cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <param name="basisEnum">set of basis objects</param>
        /// <param name="ReqB">true, if the gradient of the basis is required</param>
        /// <param name="ReqGrad">true, if the gradient of the basis is required</param>
        /// <param name="TestValues"></param>
        /// <param name="TestGradientValues"></param>
        /// <param name="sectionsTest"></param>
        /// <param name="Nodes">input, nodes where the basis should be evaluated</param>
        static internal void EvalBasis(int i0, int Len, IList<Basis> basisEnum, bool[] ReqB, bool[] ReqGrad, out MultidimensionalArray[] TestValues, out MultidimensionalArray[] TestGradientValues, out int[,] sectionsTest, NodeSet Nodes) {
            
            TestGradientValues = new MultidimensionalArray[basisEnum.Count];
            TestValues = new MultidimensionalArray[basisEnum.Count];
            sectionsTest = new int[TestValues.Length, 2];

            Debug.Assert(basisEnum.Count == ReqB.Length);
            Debug.Assert(basisEnum.Count == ReqGrad.Length);


            int i = 0;
            foreach (var basis in basisEnum) {
                if (basis is XDGBasis) {
                    var Xbasis = ((XDGBasis)basis);

                    TestValues[i] = ReqB[i] ? Xbasis.NonX_Basis.CellEval(Nodes, i0, Len) : null;
                    TestGradientValues[i] = ReqGrad[i] ? Xbasis.NonX_Basis.CellEvalGradient(Nodes, i0, Len) : null;

                    sectionsTest[i, 0] = 0;
                    sectionsTest[i, 1] = 1;
                } else {
                    TestValues[i] = ReqB[i] ? basis.CellEval(Nodes, i0, Len) : null;
                    TestGradientValues[i] = ReqGrad[i] ? basis.CellEvalGradient(Nodes, i0, Len) : null;

                    sectionsTest[i, 0] = 0;
                    sectionsTest[i, 1] = 0;
                }
                i++;
            }
        }

        /// <summary>
        /// result of the integration
        /// </summary>
        private M OperatorMatrix;

        /// <summary>
        /// affine part of da operator.
        /// </summary>
        private V OperatorAffine;


        /// <summary>
        /// writes the dammed result of the integration to the sparse matrix
        /// </summary>
        protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {

            int Mmax = m_ColMap.MaxTotalNoOfCoordinatesPerCell;
            int Nmax = m_RowMap.MaxTotalNoOfCoordinatesPerCell;

            int[] _i0 = new int[] { int.MaxValue, 0, 1 };
            int[] _iE = new int[] { -1, Nmax - 1, _i0[2] + Mmax - 1 };

            int[] _i0aff = new int[] { int.MaxValue, 0, 0 };
            int[] _iEaff = new int[] { -1, Nmax - 1, -1 };

            bool saveMtx = this.OperatorMatrix != null;
            bool saveAff = this.OperatorAffine != null;

            // loop over cells...
            for (int i = 0; i < Length; i++) {
                int jCell = i + i0;
                int Row0 = m_RowMap.LocalUniqueCoordinateIndex(0, jCell, 0);
                int Row0_g = m_RowMap.i0 + Row0;

                if (saveMtx) {
                    _i0[0] = i;
                    _iE[0] = -1;

                    var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(_i0, _iE);

                    int Col0_g = m_ColMap.GlobalUniqueCoordinateIndex(0, jCell, 0);

                    
                    //for (int r = BlockRes.NoOfRows - 1; r >= 0; r--)
                    //    for (int c = BlockRes.NoOfCols - 1; c >= 0; c--)
                    //        OperatorMatrix[Row0_g + r, Col0_g + c] += BlockRes[r, c];
                    OperatorMatrix.AccBlock(Row0_g, Col0_g, 1.0, BlockRes);
                }

                if (saveAff) {
                    _i0aff[0] = i;
                    _iEaff[0] = -1;

                    var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(_i0aff, _iEaff);

                    for (int r = BlockRes.GetLength(0) - 1; r >= 0; r--)
                        OperatorAffine[Row0 + r] += BlockRes[r];
                }
            }
        }
    }
}
