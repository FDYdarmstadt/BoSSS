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

        ConventionalDGField[] m_ParametersA;
        ConventionalDGField[] m_ParametersB;

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
        internal int m_LevSetIdx {
            get;
            private set;
        }

        
        /// <summary>
        /// les tracker
        /// </summary>
        LevelSetTracker m_lsTrk;

        /// <summary>
        /// index into <see cref="LevelSetTracker.RegionsHistory"/>, etc.
        /// </summary>
        int m_LsTrkHistoryIndex;

        /*
        /// <summary>
        /// Normals, etc.
        /// </summary>
        LevelSetTracker.LevelSetData m_lsData;

        /// <summary>
        /// cut-cell masks, etc.
        /// </summary>
        LevelSetTracker.LevelSetRegions m_lsRegions;
        */

        /// <summary>
        /// Negative and positive (with respect to level-set) species.
        /// </summary>
        Tuple<SpeciesId,SpeciesId> m_SpeciesPair;

        /// <summary>
        /// Negative species/Species A
        /// </summary>
        public SpeciesId SpeciesA {
            get {
                return m_SpeciesPair.Item1;
            }
        }

        /// <summary>
        /// Positive species/Species B
        /// </summary>
        public SpeciesId SpeciesB {
            get {
                return m_SpeciesPair.Item2;
            }
        }

        public override Quadrature<QuadRule, CellMask> CloneForThreadParallelization() {
            return new LECQuadratureLevelSet<M, V>(
                this.gridData, m_DiffOp,
                this.OperatorMatrix, this.OperatorAffine,
                m_RowMap, m_ParametersA, m_ColMap,
                this.m_lsTrk, this.m_LevSetIdx, this.m_LsTrkHistoryIndex,
                this.m_SpeciesPair,
                base.m_compositeRule) {
                m_ParametersA = this.m_ParametersA,
                m_ParametersB = this.m_ParametersB
            };
        }

        XDifferentialOperatorMk2 m_DiffOp;
        /// <summary>
        /// ctor.
        /// </summary>
        internal LECQuadratureLevelSet(IGridData context,
                                     XDifferentialOperatorMk2 DiffOp,
                                     M Matrix, V OffsetVec,
                                     UnsetteledCoordinateMapping RowMap, IList<DGField> ParamsMap, UnsetteledCoordinateMapping ColMap,
                                     LevelSetTracker lsTrk, int _iLevSet, int TrackerHistoryIndex,
                                     Tuple<SpeciesId,SpeciesId> SpeciesPair,
                                     //Tuple<CoefficientSet,CoefficientSet> NegPosCoeff,
                                     ICompositeQuadRule<QuadRule> domAndRule) //
            : base(new int[] { RowMap.GetNonXBasisLengths(0).Sum()*2, 1 + ((Matrix == null) ? 0 : ColMap.GetNonXBasisLengths(0).Sum()*2) }, // we always integrate over species in pairs (neg + pos), so we need to alloc mem only 2 species
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
            m_DiffOp = DiffOp;
            m_RowMap = RowMap;
            m_ColMap = ColMap;

            int Gamma = m_RowMap.BasisS.Count;
            m_lsTrk = lsTrk;
            m_LsTrkHistoryIndex = TrackerHistoryIndex;

            if (Matrix != null && (Matrix.RowPartitioning.LocalLength != RowMap.LocalLength))
                throw new ArgumentException("mismatch between matrix number of rows and row mapping.");
            if (Matrix != null && (Matrix.ColPartition.LocalLength != ColMap.LocalLength))
                throw new ArgumentException("mismatch between matrix number of columns and column mapping.");

            this.m_LevSetIdx = _iLevSet;
            this.m_SpeciesPair = SpeciesPair;

            this.OperatorMatrix = Matrix;
            this.OperatorAffine = OffsetVec;

            var _Parameters = (ParamsMap != null) ? ParamsMap.ToArray() : new DGField[0];
            m_ParametersA = new ConventionalDGField[_Parameters.Length];
            m_ParametersB = new ConventionalDGField[_Parameters.Length];
            for(int i = 0; i < _Parameters.Length; i++) {
                var f = _Parameters[i];
                if(f == null) {
                    m_ParametersA[i] = null;
                    m_ParametersB[i] = null;
                } else if(f is XDGField xf) {
                    m_ParametersA[i] = xf.GetSpeciesShadowField(this.SpeciesA);
                    m_ParametersB[i] = xf.GetSpeciesShadowField(this.SpeciesB);
                } else if(f is ConventionalDGField cf) {
                    m_ParametersA[i] = cf;
                    m_ParametersB[i] = null;
                } else {
                    throw new NotImplementedException("missing implementation for " + f.GetType().Name);
                }
            }

            if (m_RowMap.BasisS.Count != DiffOp.CodomainVar.Count)
                throw new ArgumentException("mismatch between number of codomain variables in spatial operator and row mapping");
            if (m_ColMap.BasisS.Count != DiffOp.DomainVar.Count)
                throw new ArgumentException("mismatch between number of domain variables in spatial operator and column mapping");
            if (_Parameters.Length != DiffOp.ParameterVar.Count)
                throw new ArgumentException("mismatch between number of parameter variables in spatial operator and given parameters");



            //RowXbSw = m_RowMap.XorNonXbasis().Select(b => b ? 1 : 0).ToArray();
            //ColXbSw = m_ColMap.XorNonXbasis().Select(b => b ? 1 : 0).ToArray();


            TestNegativeAndPositiveSpecies(domAndRule, m_lsTrk, m_LsTrkHistoryIndex, SpeciesA, SpeciesB, m_LevSetIdx);

            // ------------------------
            // sort equation components
            // ------------------------

            //m_BiLinForms = DiffOp.GetArgMapping<IBilinearForm>(true, Compfilter<IBilinearForm>);
            //m_2ndDerivFlux = DiffOp.GetArgMapping<ILinear2ndDerivativeCouplingFlux>(true, Compfilter<ILinear2ndDerivativeCouplingFlux>);

            m_LsForm_UxV = EquationComponentArgMapping<ILevelSetForm_UxV>.GetArgMapping(DiffOp, true,
               eq => ((eq.LevelSetTerms & TermActivationFlags.UxV) != 0) && Compfilter(eq),
               eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_GradUxV = EquationComponentArgMapping<ILevelSetForm_GradUxV>.GetArgMapping(DiffOp, true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.GradUxV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_UxGradV = EquationComponentArgMapping<ILevelSetForm_UxGradV>.GetArgMapping(DiffOp, true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.UxGradV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_GradUxGradV = EquationComponentArgMapping<ILevelSetForm_GradUxGradV>.GetArgMapping(DiffOp, true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.GradUxGradV) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_V = EquationComponentArgMapping<ILevelSetForm_V>.GetArgMapping(DiffOp, true,
                eq => ((eq.LevelSetTerms & TermActivationFlags.V) != 0 && Compfilter(eq)),
                eq => (eq is ILevelSetForm) ? new LinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_LsForm_GradV = EquationComponentArgMapping<ILevelSetForm_GradV>.GetArgMapping(DiffOp, true,
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

            if ((this.m_lsTrk.GetSpeciesName(this.SpeciesA) == b.NegativeSpecies && this.m_lsTrk.GetSpeciesName(this.SpeciesB) == b.PositiveSpecies)
                || (this.m_lsTrk.GetSpeciesName(this.SpeciesB) == b.NegativeSpecies && this.m_lsTrk.GetSpeciesName(this.SpeciesA) == b.PositiveSpecies)) {

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
            Debug.Assert(m_ParametersA.Length == m_ParametersB.Length);
            if (m_ParametersA.Length > 0) {
                m_ParamFieldValuesPos = MultidimensionalArray.Create(m_ParametersB.Length, Nitm, Nnod);
                m_ParamFieldValuesNeg = MultidimensionalArray.Create(m_ParametersA.Length, Nitm, Nnod);
            }
            int D = gridData.SpatialDimension;

            int _Delta = m_ColMap.BasisS.Count;

            Alloc(m_LsForm_UxV, Koeff_UxV, Sum_Koeff_UxV, 4, Nitm, Nnod, 2, 2, _Delta);
            Alloc(m_LsForm_GradUxV, Koeff_NablaUxV, Sum_Koeff_NablaUxV, 4, Nitm, Nnod, 2, 2, _Delta, D);
            Alloc(m_LsForm_UxGradV, Koeff_UxNablaV, Sum_Koeff_UxNablaV, 4, Nitm, Nnod, 2, 2, _Delta, D);
            Alloc(m_LsForm_GradUxGradV, Koeff_NablaUxNablaV, Sum_Koeff_NablaUxNablaV, 4, Nitm, Nnod, 2, 2, _Delta, D, D);

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


            TestNegativeAndPositiveSpecies(i0, Len, m_lsTrk, m_LsTrkHistoryIndex, this.SpeciesA, this.SpeciesB,  this.m_LevSetIdx);




            // Evaluate Parameter fields
            // -------------------------
            base.CustomTimers[3].Start();
            for (int i = 0; i < m_ParametersA.Length; i++) {
                var bufPos = m_ParamFieldValuesPos.ExtractSubArrayShallow(i, -1, -1);
                var bufNeg = m_ParamFieldValuesNeg.ExtractSubArrayShallow(i, -1, -1);
                if (m_ParametersA[i] != null || m_ParametersB[i] != null) {
                    if (m_ParametersB[i] != null) {
                        // jump in parameter i at level-set: separate evaluation for both sides
                        m_ParametersA[i].Evaluate(i0, Len, QuadNodes, bufNeg);
                        m_ParametersB[i].Evaluate(i0, Len, QuadNodes, bufPos);

                    } else {
                        // no jump at level set: positive and negative limit of parameter i are equal
                        m_ParametersA[i].Evaluate(i0, Len, QuadNodes, bufPos);
                        bufNeg.Set(bufPos);
                    }
                } else {
                    bufPos.Clear();
                    bufNeg.Clear();
                }
            }

            // Evaluate level sets and normals
            // -------------------------------
            MultidimensionalArray Normals = m_lsTrk.DataHistories[m_LevSetIdx][m_LsTrkHistoryIndex].GetLevelSetNormals(QuadNodes, i0, Len);
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
            base.CustomTimers[1].Start();
            EvalBasis(i0, Len, this.m_RowMap.BasisS, ReqV, ReqGradV, out TestValues, out TestGradientValues, QuadNodes);
            EvalBasis(i0, Len, this.m_ColMap.BasisS, ReqU, ReqGradU, out BasisValues, out BasisGradientValues, QuadNodes);
            base.CustomTimers[1].Stop();

            // compute offsets into matrix blocks
            // ----------------------------------

            
            CompOffsets(i0, Len, out int[] offsetDom, out int[] NnonxDom, m_ColMap);
            CompOffsets(i0, Len, out int[] offsetCod, out int[] NnonxCod, m_RowMap);


            // Transform nodes to global coordinates
            // -------------------------------------

            MultidimensionalArray NodesGlobal = gridData.GlobalNodes.GetValue_Cell(QuadNodes, i0, Len);

            // Evaluate Integral components
            // ----------------------------

            //var _inParams = new BoSSS.Foundation.XDG.LevSetIntParams();
            EdgeFormParams _inParams = default(EdgeFormParams);
            _inParams.e0 = i0;
            _inParams.Len = Len;
            
            
            // loop over codomain variables ...
            for (int gamma = 0; gamma < GAMMA; gamma++) {


                // prepare parameters
                // - - - - - - - - - 

                // set Normal's
                _inParams.Normals = Normals;
                // set Nodes Global
                _inParams.Nodes = NodesGlobal;
                _inParams.time = this.time;
                _inParams.GridDat = this.GridDat;

                Debug.Assert(i0 == _inParams.e0);
                Debug.Assert(Len == _inParams.Len);

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
                    EvalComponent(ref _inParams, gamma, this.m_LsForm_UxV[gamma], this.m_LsForm_UxV_Watches[gamma], 
                        Koeff_UxV, Sum_Koeff_UxV, 4,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_UxV _comp, int _gamma, int i, ref EdgeFormParams inp) {
                            _comp.InternalEdge_UxV(ref inp, Koeff_UxV[_gamma][i]);
                        });
                }
                {
                    EvalComponent(ref _inParams, gamma, this.m_LsForm_GradUxV[gamma], this.m_LsForm_GradUxV_Watches[gamma],
                        Koeff_NablaUxV, Sum_Koeff_NablaUxV, 4,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_GradUxV _comp, int _gamma, int i, ref EdgeFormParams inp) {
                            _comp.InternalEdge_GradUxV(ref inp, Koeff_NablaUxV[_gamma][i]);
                        });
                }
                {
                    EvalComponent(ref _inParams, gamma, this.m_LsForm_UxGradV[gamma], this.m_LsForm_UxGradV_Watches[gamma],
                        Koeff_UxNablaV, Sum_Koeff_UxNablaV, 4,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_UxGradV _comp, int _gamma, int i, ref EdgeFormParams inp) {
                            _comp.InternalEdge_UxGradV(ref inp, Koeff_UxNablaV[_gamma][i]);
                        });
                }
                {
                    EvalComponent(ref _inParams, gamma, this.m_LsForm_GradUxGradV[gamma], this.m_LsForm_GradUxGradV_Watches[gamma],
                        Koeff_NablaUxNablaV, Sum_Koeff_NablaUxNablaV, 4,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_GradUxGradV _comp, int _gamma, int i, ref EdgeFormParams inp) {
                            _comp.InternalEdge_GradUxGradV(ref inp, Koeff_NablaUxNablaV[_gamma][i]);
                        });
                }

                {
                    EvalComponent(ref _inParams, gamma, this.m_LsForm_V[gamma], this.m_LsForm_V_Watches[gamma],
                        Koeff_V, Sum_Koeff_V, -1,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_V _comp, int _gamma, int i, ref EdgeFormParams inp) {
                            _comp.InternalEdge_V(ref inp, Koeff_V[_gamma][i]);
                        });
                }
                {
                    EvalComponent(ref _inParams, gamma, this.m_LsForm_GradV[gamma], this.m_LsForm_GradV_Watches[gamma],
                        Koeff_NablaV, Sum_Koeff_NablaV, -1,
                        m_ParamFieldValuesPos, m_ParamFieldValuesNeg,
                        DELTA,
                        base.CustomTimers[0],
                        delegate (ILevelSetForm_GradV _comp, int _gamma, int i, ref EdgeFormParams inp) {
                            _comp.InternalEdge_GradV(ref inp, Koeff_NablaV[_gamma][i]);
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

                if(N <= 0)
                    continue; // nothing to do for codomain/row variable 'gamma'

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

                    Debug.Assert(NnonxCod[gamma] == N);

                    // Das Schleifenmonster: ---------------------------------------------
                    for (int cr = 0; cr < 2; cr++) {  // loop over neg/pos species, row...
                        for (int cc = 0; cc < 2; cc++) {  // loop over neg/pos species, column...

                            if(M <= 0)
                                continue; // nothing to do for domain/column variable 'delta'
                            Debug.Assert(NnonxDom[delta] == M);
                                                        
                            int[] extr0 = new int[] { 0, 0, 
                                cr * N + offsetCod[gamma], // row
                                cc * M + 1 + offsetDom[delta] }; // col
                                //sectionsTest[gamma, cr] * N + offsetCod[gamma], // row
                                //sectionsBasis[delta, cc] * M + 1 + offsetDom[delta] }; // col
                            int[] extrE = new int[] { Len - 1, NoOfNodes - 1, 
                                extr0[2] + N - 1, // row
                                extr0[3] + M - 1 }; // col
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
                                var Sum_Koeff_NablaUxNablaV_CrCc = Sum_Koeff_NablaUxNablaV[gamma, delta].ExtractSubArrayShallow(-1, -1, cr, cc, -1, -1);
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
                        //int[] extr0 = new int[] { 0, 0, sectionsTest[gamma, cr] * N + offsetCod[gamma], 0 };
                        int[] extr0 = new int[] { 0, 0, cr * N + offsetCod[gamma], 0 };
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

        /// <summary>
        /// Test for the correct configuration of positive and negative species.
        /// Should be pretty quick, so we can do this also in Release.
        /// </summary>
        internal static void TestNegativeAndPositiveSpecies(ICompositeQuadRule<QuadRule> rule, LevelSetTracker lsTrk, int HistoryIndex, SpeciesId spcNeg, SpeciesId spcPos, int iLevSet) {
            foreach(var crp in rule) {

                TestNegativeAndPositiveSpecies(crp.Chunk.i0, crp.Chunk.Len, lsTrk, HistoryIndex, spcNeg, spcPos, iLevSet);
            }
        }


        /// <summary>
        /// Test for the correct configuration of positive and negative species.
        /// Should be pretty quick, so we can do this also in Release.
        /// </summary>
        internal static void TestNegativeAndPositiveSpecies(int i0, int Len, LevelSetTracker lsTrk, int HistoryIndex, SpeciesId spcNeg, SpeciesId spcPos, int iLevSet) {
            //var negSignCodes = lsTrk.GetLevelSetSignCodes(spcNeg);
            //var posSignCodes = lsTrk.GetLevelSetSignCodes(spcPos);

            var Regions = lsTrk.RegionsHistory[HistoryIndex];

            for(int j = i0; j < (i0 + Len); j++) {
                Regions.GetSpeciesIndex(spcNeg, j);
                Regions.GetSpeciesIndex(spcPos, j);

                //ushort RegionCode = m_lsTrk.Regions.m_LevSetRegions[j];
                LevelsetCellSignCode csc = Regions.GetCellSignCode(j);

                if(!(csc.GetSign(iLevSet) == LevelsetSign.Both))
                    throw new ApplicationException("Seem to perform level-set integration in a non-cut cell.");


                var cscNeg = csc; cscNeg.SetSign(iLevSet, LevelsetSign.Negative);
                if(lsTrk.ContainesSpecies(spcNeg, cscNeg) == false)
                    throw new ApplicationException("Pos/Neg species mishmash."); // for negative sign, cell MUST contain negative species
                if(lsTrk.ContainesSpecies(spcPos, cscNeg) == true)
                    throw new ApplicationException("Pos/Neg species mishmash."); // for negative sign, cell should NOT contain positive species

                var cscPos = csc; cscPos.SetSign(iLevSet, LevelsetSign.Positive);
                if(lsTrk.ContainesSpecies(spcPos, cscPos) == false)
                    throw new ApplicationException("Pos/Neg species mishmash."); // for positive sign, cell MUST contain positive species
                if(lsTrk.ContainesSpecies(spcNeg, cscPos) == true)
                    throw new ApplicationException("Pos/Neg species mishmash."); // for positive sign, cell should NOT contain negative species
            }
        }

        static internal void CompOffsets(int i0, int L, out int[] offset, out int[] Nnonx, UnsetteledCoordinateMapping Map) {
            int DELTA = Map.BasisS.Count;
            Nnonx = Map.GetNonXBasisLengths(i0);
            Debug.Assert(DELTA == Nnonx.Length);
            offset = new int[DELTA];
            for(int delta = 1; delta < DELTA; delta++) {
                offset[delta] = offset[delta - 1] + Nnonx[delta - 1] * 2;
            }

#if DEBUG
            
            for (int i = 1; i < L; i++) {
                int[] __Nnonx = Map.GetNonXBasisLengths(i0 + i);
                Debug.Assert(ArrayTools.ListEquals(Nnonx, __Nnonx));
            }
#endif
        }

        delegate void CallComponent<T>(T comp, int gamma, int i, ref EdgeFormParams _inParams);

        static private void EvalComponent<T>(ref EdgeFormParams _inParams,
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
                _inParams.ParameterVars_OUT = new MultidimensionalArray[NoOfParams];
                _inParams.ParameterVars_IN = new MultidimensionalArray[NoOfParams];
                for (int c = 0; c < NoOfParams; c++) {
                    int targ = bf.AllToSub[i, c + NoOfArgs] - DELTA;
                    Debug.Assert(targ >= 0);
                    _inParams.ParameterVars_OUT[c] = ParamFieldValuesPos.ExtractSubArrayShallow(targ, -1, -1);
                    _inParams.ParameterVars_IN[c] = ParamFieldValuesNeg.ExtractSubArrayShallow(targ, -1, -1);
                }

                // evaluate equation components
                timers[i].Start();
                ComponentFunc(comp, gamma, i, ref _inParams);
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
        /// <param name="Nodes">input, nodes where the basis should be evaluated</param>
        static internal void EvalBasis(int i0, int Len, IList<Basis> basisEnum, bool[] ReqB, bool[] ReqGrad, out MultidimensionalArray[] TestValues, out MultidimensionalArray[] TestGradientValues, NodeSet Nodes) {
            
            TestGradientValues = new MultidimensionalArray[basisEnum.Count];
            TestValues = new MultidimensionalArray[basisEnum.Count];
            //sectionsTest = new int[TestValues.Length, 2];

            Debug.Assert(basisEnum.Count == ReqB.Length);
            Debug.Assert(basisEnum.Count == ReqGrad.Length);


            int i = 0;
            foreach (var basis in basisEnum) {
                if (basis is XDGBasis Xbasis) {
                    
                    TestValues[i] = ReqB[i] ? Xbasis.NonX_Basis.CellEval(Nodes, i0, Len) : null;
                    TestGradientValues[i] = ReqGrad[i] ? Xbasis.NonX_Basis.CellEvalGradient(Nodes, i0, Len) : null;

                    //sectionsTest[i, 0] = 0;
                    //sectionsTest[i, 1] = 1;
                } else {
                    TestValues[i] = ReqB[i] ? basis.CellEval(Nodes, i0, Len) : null;
                    TestGradientValues[i] = ReqGrad[i] ? basis.CellEvalGradient(Nodes, i0, Len) : null;

                    //sectionsTest[i, 0] = 0;
                    //sectionsTest[i, 1] = 0;
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
        /// scaling factor for accumulation
        /// </summary>
        internal double m_alpha = 1.0;

        /// <summary>
        /// writes the dammed result of the integration to the sparse matrix
        /// </summary>
        protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {

            int DELTA = m_ColMap.BasisS.Count; // DELTA: number of domain/column/trial variables
            int GAMMA = m_RowMap.BasisS.Count;  // GAMMA: number of codom/row/test variables

            CompOffsets(i0, Length, out int[] offsetCol, out int[] ColNonxN, m_ColMap);
            CompOffsets(i0, Length, out int[] offsetRow, out int[] RowNonxN, m_RowMap);
            

            SpeciesId[] spcS = new[] { this.SpeciesA, this.SpeciesB };

            int[] _i0 = new int[3];// { int.MaxValue, 0, 1 };
            int[] _iE = new int[3];// { -1, Nmax - 1, _i0[2] + Mmax - 1 };

            int[] _i0aff = new int[3];// { int.MaxValue, 0, 0 };
            int[] _iEaff = new int[3];// { -1, Nmax - 1, -1 };

            bool saveMtx = this.OperatorMatrix != null;
            bool saveAff = this.OperatorAffine != null;

            double a = m_alpha;

            // loop over cells...
            for (int i = 0; i < Length; i++) {
                int jCell = i + i0;
#if DEBUG
                Debug.Assert(RowNonxN.ListEquals(m_RowMap.GetNonXBasisLengths(jCell)));
                Debug.Assert(ColNonxN.ListEquals(m_ColMap.GetNonXBasisLengths(jCell)));
#endif

                _i0[0] = i;
                _iE[0] = -1;

                _i0aff[0] = i;
                _iEaff[0] = -1;
                _i0aff[2] = 0;
                _iEaff[2] = -1;

                for(int gamma = 0; gamma < GAMMA; gamma++) { // loop over rows...
                    for(int cr = 0; cr < 2; cr++) { // loop over neg/pos species row...
                        SpeciesId rowSpc = spcS[cr];
                        int Row0 = m_RowMap.LocalUniqueCoordinateIndex(m_lsTrk, gamma, jCell, rowSpc, 0);
                        long Row0_g = m_RowMap.i0 + Row0;

                        //_i0aff[1] = offsetRow[gamma] + RowXbSw[gamma] * RowNonxN[gamma] * cr; // the 'RowXbSw' is 0 for non-xdg, so both species will be added
                        _i0aff[1] = offsetRow[gamma] + RowNonxN[gamma] * cr;
                        _iEaff[1] = _i0aff[1] + RowNonxN[gamma] - 1;


                        if(saveMtx) {
                            _i0[1] = _i0aff[1];
                            _iE[1] = _iEaff[1];

                            for(int delta = 0; delta < DELTA; delta++) {
                                for(int cc = 0; cc < 2; cc++) {
                                    SpeciesId colSpc = spcS[cc];
                                    int Col0 = m_ColMap.LocalUniqueCoordinateIndex(m_lsTrk, delta, jCell, colSpc, 0);
                                    long Col0_g = m_ColMap.i0 + Col0;

                                    //_i0[2] = offsetCol[delta] + ColXbSw[delta] * ColNonxN[delta] * cc + 1;
                                    _i0[2] = offsetCol[delta] + ColNonxN[delta] * cc + 1;
                                    _iE[2] = _i0[2] + ColNonxN[delta] - 1;

                                    var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(_i0, _iE);

                                    OperatorMatrix.AccBlock(Row0_g, Col0_g, a, BlockRes);
                                }
                            }
                        }

                        if(saveAff) {
                            var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(_i0aff, _iEaff);

                            for(int r = BlockRes.GetLength(0) - 1; r >= 0; r--)
                                OperatorAffine[Row0 + r] += BlockRes[r]*a;
                        }
                    }
                }
            }
        }

        
    }
}
