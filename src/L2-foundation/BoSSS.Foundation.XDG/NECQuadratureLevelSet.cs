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
    internal class NECQuadratureLevelSet<V> : BoSSS.Foundation.Quadrature.CellQuadrature  
        where V : IList<double> {

        /// <summary>
        /// Result of quadrature, Output storage
        /// </summary>
        public V ResultVector;
        
        /// <summary>
        /// Mapping of DG-DOF into <see cref="ResultVector"/>
        /// </summary>
        UnsetteledCoordinateMapping m_CodomainMap;
                
        /// <summary>
        /// All Domain and Parameter fields (domain first (0 to <see cref="DELTA"/>-1), then parameters) are defined by the operator.
        /// </summary>
        DGField[] m_DomainAndParamFields;

        /// <summary>
        /// Number of domain fields
        /// </summary>
        int DELTA; 
        
        /// <summary>
        /// Switch, whether the evaluation (the field value) of <see cref="m_DomainAndParamFields"/> is actually necessary.
        /// </summary>
        bool[] m_ValueRequired;

        /// <summary>
        /// Switch, whether the evaluation (the gradient value) of <see cref="m_DomainAndParamFields"/> is actually necessary.
        /// </summary>
        bool[] m_GradientRequired;
        
        /// <summary>
        /// ye good old level set tracker
        /// </summary>
        LevelSetTracker m_lsTrk;
        
        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// index of the level set to evaluate
        /// </summary>
        int m_LevSetIdx;

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
        internal NECQuadratureLevelSet(IGridData context,
                                     XSpatialOperator DiffOp,
                                     V __ResultVector,
                                     IList<DGField> __DomainFields,
                                     IList<DGField> __Parameters,
                                     UnsetteledCoordinateMapping CodomainMap,
                                     LevelSetTracker lsTrk, int _iLevSet, Tuple<SpeciesId, SpeciesId> SpeciesPair,
                                     ICompositeQuadRule<QuadRule> domAndRule) //
            : base(new int[] { CodomainMap.NoOfCoordinatesPerCell },
                 context,
                 domAndRule) //
        {



            // -----------------------------------
            // set members / check ctor parameters
            // -----------------------------------
            m_lsTrk = lsTrk;
            this.m_LevSetIdx = _iLevSet;
            this.m_SpeciesPair = SpeciesPair;
            this.ResultVector = __ResultVector;
            m_CodomainMap = CodomainMap;
            var _Parameters = (__Parameters != null) ? __Parameters.ToArray() : new DGField[0];

            if(__DomainFields.Count != DiffOp.DomainVar.Count)
                throw new ArgumentException("mismatch between number of domain variables in spatial operator and given domain variables");
            if(_Parameters.Length != DiffOp.ParameterVar.Count)
                throw new ArgumentException("mismatch between number of parameter variables in spatial operator and given parameters");
            if(m_CodomainMap.NoOfVariables != DiffOp.CodomainVar.Count)
                throw new ArgumentException("mismatch between number of codomain variables in spatial operator and given codomain mapping");

            m_DomainAndParamFields = ArrayTools.Cat(__DomainFields, _Parameters);
            this.DELTA = __DomainFields.Count;

            // ------------------------
            // sort equation components
            // ------------------------

            int Gamma = m_CodomainMap.NoOfVariables;

            m_NonlinLsForm_V = DiffOp.GetArgMapping<INonlinLevelSetForm_V>(true,
               eq => ((eq.LevelSetTerms & (TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradUxV)) != 0) && Compfilter(eq),
               eq => (eq is ILevelSetForm) ? new NonlinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);
            m_NonlinLsForm_GradV = DiffOp.GetArgMapping<INonlinLevelSetForm_GradV>(true,
                eq => ((eq.LevelSetTerms & (TermActivationFlags.GradV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxGradV)) != 0) && Compfilter(eq),
                eq => (eq is ILevelSetForm) ? new NonlinearLevelSetFormVectorizer((ILevelSetForm)eq) : null);


            m_ValueRequired = new bool[m_DomainAndParamFields.Length];
            m_GradientRequired = new bool[m_DomainAndParamFields.Length];

            m_NonlinLsForm_V.DetermineReqFields(m_GradientRequired,
                comp => ((comp.LevelSetTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0));
            m_NonlinLsForm_GradV.DetermineReqFields(m_GradientRequired,
                comp => ((comp.LevelSetTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0));
            m_NonlinLsForm_V.DetermineReqFields(m_ValueRequired,
                comp => ((comp.LevelSetTerms & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0));
            m_NonlinLsForm_GradV.DetermineReqFields(m_ValueRequired,
                comp => ((comp.LevelSetTerms & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0));

            for(int i = __DomainFields.Count; i < m_DomainAndParamFields.Length; i++) {
                m_ValueRequired[i] = true; // parameters are always required!
            }


            // -----
            // alloc
            // -----

            Koeff_V = new MultidimensionalArray[Gamma];
            Koeff_GradV = new MultidimensionalArray[Gamma];
            for(int gamma = 0; gamma < Gamma; gamma++) {
                Koeff_V[gamma] = new MultidimensionalArray(3);
                Koeff_GradV[gamma] = new MultidimensionalArray(4);
            }

            int L = m_DomainAndParamFields.Length;
            m_FieldValuesPos = new MultidimensionalArray[L];
            m_FieldValuesNeg = new MultidimensionalArray[L];
            m_FieldGradientValuesPos = new MultidimensionalArray[L];
            m_FieldGradientValuesNeg = new MultidimensionalArray[L];

            for(int l = 0; l < L; l++) {
                if(m_ValueRequired[l]) {
                    m_FieldValuesPos[l] = new MultidimensionalArray(2);
                    m_FieldValuesNeg[l] = new MultidimensionalArray(2);
                }

                if(m_GradientRequired[l]) {
                    m_FieldGradientValuesPos[l] = new MultidimensionalArray(3);
                    m_FieldGradientValuesNeg[l] = new MultidimensionalArray(3);
                }
            }


            // ------------------
            // init custom timers
            // ------------------

            base.CustomTimers = new Stopwatch[] { new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch() };
            base.CustomTimers_Names = new string[] { "Flux-Eval", "Basis-Eval", "Loops", "ParametersAndNormals", "Field-Eval" };
            base.CustomTimers_RootPointer = new int[5];
            ArrayTools.SetAll(base.CustomTimers_RootPointer, -1);

            this.m_NonlinLsForm_V_Watches = this.m_NonlinLsForm_V.InitStopWatches(0, this);
            this.m_NonlinLsForm_GradV_Watches = this.m_NonlinLsForm_GradV.InitStopWatches(0, this);

            Flux_Eval = base.CustomTimers[0];
            Basis_Eval = base.CustomTimers[1];
            Loops = base.CustomTimers[2];
            ParametersAndNormals = base.CustomTimers[3];
            Field_Eval = base.CustomTimers[4];
        }
        

       private bool Compfilter(IEquationComponent c) {

            ILevelSetForm b = (ILevelSetForm)c;
            
            if (this.m_LevSetIdx != b.LevelSetIndex)
                // component is not relevant for this level-set
                return false;

            if (!(this.SpeciesA == b.NegativeSpecies && this.SpeciesB == b.PositiveSpecies))
                // component is not relevant for this level-set
                return false;

            // filter passed
            return true;
        }
        
        Stopwatch Flux_Eval;
        Stopwatch Field_Eval;
        Stopwatch Basis_Eval;
        Stopwatch Loops;
        Stopwatch ParametersAndNormals;

        EquationComponentArgMapping<INonlinLevelSetForm_V>[] m_NonlinLsForm_V;
        EquationComponentArgMapping<INonlinLevelSetForm_GradV>[] m_NonlinLsForm_GradV;

        Stopwatch[][] m_NonlinLsForm_V_Watches;
        Stopwatch[][] m_NonlinLsForm_GradV_Watches;

        /// <summary>
        /// values of domain & parameter fields, positive side of level Set 
        /// <list type="bullet">
        ///   <item>1st index: parameter index, correlates with <see cref="m_DomainAndParamFields"/></item>
        ///   <item>2nd index: cell index, plus some offset</item>
        ///   <item>3rd index: node index</item>
        /// </list>
        /// </summary>
        MultidimensionalArray[] m_FieldValuesPos;
        
        /// <summary>
        /// values of domain & parameter fields, positive side of level Set
        /// <list type="bullet">
        ///   <item>1st index: parameter index, correlates with <see cref="m_DomainAndParamFields"/></item>
        ///   <item>2nd index: cell index, plus some offset</item>
        ///   <item>3rd index: node index</item>
        /// </list>
        /// </summary>
        MultidimensionalArray[] m_FieldValuesNeg;


        /// <summary>
        /// values of domain & parameter fields, positive side of level Set 
        /// <list type="bullet">
        ///   <item>1st index: parameter index, correlates with <see cref="m_DomainAndParamFields"/></item>
        ///   <item>2nd index: cell index, plus some offset</item>
        ///   <item>3rd index: node index</item>
        ///   <item>4th index: spatial direction</item>
        /// </list>
        /// </summary>
        MultidimensionalArray[] m_FieldGradientValuesPos;
        
        /// <summary>
        /// values of domain & parameter fields, positive side of level Set
        /// <list type="bullet">
        ///   <item>1st index: parameter index, correlates with <see cref="m_DomainAndParamFields"/></item>
        ///   <item>2nd index: cell index, plus some offset</item>
        ///   <item>3rd index: node index</item>
        ///   <item>4th index: spatial direction</item>
        /// </list>
        /// </summary>
        MultidimensionalArray[] m_FieldGradientValuesNeg;




        /// <summary>
        /// result buffers for <see cref="m_NonlinLsForm_V"/>,
        /// 1st index: codomain variable/test function
        /// </summary>
        MultidimensionalArray[] Koeff_V;

        /// <summary>
        /// result buffers for <see cref="m_NonlinLsForm_GradV"/>, 
        /// 1st index: codomain variable/test function
        /// </summary>
        MultidimensionalArray[] Koeff_GradV;

        
        /// <summary>
        /// memory allocation.
        /// </summary>
        protected override void AllocateBuffers(int Nitm, NodeSet ruleNodes) {
            int Nnod = ruleNodes.GetLength(0);
            int D = this.GridDat.SpatialDimension;
            base.AllocateBuffers(Nitm, ruleNodes);

            int L = m_DomainAndParamFields.Length;
            for(int l = 0; l < L; l++) {
                Debug.Assert((m_FieldValuesPos[l] != null) == (m_FieldValuesNeg[l] != null));
                if(m_FieldValuesPos[l] != null) {
                    m_FieldValuesPos[l].Allocate(Nitm, Nnod);
                    m_FieldValuesNeg[l].Allocate(Nitm, Nnod);
                }

                Debug.Assert((m_FieldGradientValuesPos[l] != null) == (m_FieldGradientValuesNeg[l] != null));
                if(m_FieldGradientValuesPos[l] != null) {
                    m_FieldGradientValuesPos[l].Allocate(Nitm, Nnod, D);
                    m_FieldGradientValuesNeg[l].Allocate(Nitm, Nnod, D);
                }
            }
        }


        protected override void Evaluate(int i0, int Len, QuadRule QR, MultidimensionalArray EvalResult) {
            NodeSet QuadNodes = QR.Nodes;
            int D = gridData.SpatialDimension;
            int NoOfNodes = QuadNodes.NoOfNodes;
            int GAMMA = m_CodomainMap.NoOfVariables;  // GAMMA: number of codom variables
           

            // Evaluate Domain & Parameter fields
            // --------------------------------

            Field_Eval.Start();

            for (int i = 0; i < m_DomainAndParamFields.Length; i++) {

                if(m_ValueRequired[i]) {
                    DGField _Field = m_DomainAndParamFields[i];
                    if(_Field != null) {
                        if (_Field is XDGField) {
                            // jump in parameter i at level-set: separate evaluation for both sides
                            var _xField = _Field as XDGField;

                            _xField.GetSpeciesShadowField(this.SpeciesA).Evaluate(i0, Len, QuadNodes, m_FieldValuesNeg[i]);
                            _xField.GetSpeciesShadowField(this.SpeciesB).Evaluate(i0, Len, QuadNodes, m_FieldValuesPos[i]);

                    } else {
                        // no jump at level set: positive and negative limit of parameter i are equal
                        _Field.Evaluate(i0, Len, QuadNodes, m_FieldValuesPos[i]);
                        m_FieldValuesNeg[i].Set(m_FieldValuesPos[i]);
                    }
                    } else {
                        m_FieldValuesPos[i].Clear();
                        m_FieldValuesNeg[i].Clear();
                    }
                }

                if(m_GradientRequired[i]) {
                    DGField _Field = m_DomainAndParamFields[i];
                    if(_Field != null) {
                        if (_Field is XDGField) {
                            // jump in parameter i at level-set: separate evaluation for both sides
                            var _xField = _Field as XDGField;

                            _xField.GetSpeciesShadowField(this.SpeciesA).EvaluateGradient(i0, Len, QuadNodes, m_FieldGradientValuesNeg[i]);
                            _xField.GetSpeciesShadowField(this.SpeciesB).EvaluateGradient(i0, Len, QuadNodes, m_FieldValuesPos[i]);

                    } else {
                        // no jump at level set: positive and negative limit of parameter i are equal
                        _Field.EvaluateGradient(i0, Len, QuadNodes, m_FieldGradientValuesPos[i]);
                        m_FieldGradientValuesNeg[i].Set(m_FieldGradientValuesPos[i]);
                    }
                    } else {
                        m_FieldGradientValuesPos[i].Clear();
                        m_FieldGradientValuesNeg[i].Clear();
                    }
                }
                
            }

            Field_Eval.Stop();

            // Evaluate level sets and normals
            // -------------------------------

            ParametersAndNormals.Start();
            var NoOfLevSets = m_lsTrk.LevelSets.Count;
            MultidimensionalArray Normals = m_lsTrk.DataHistories[m_LevSetIdx].Current.GetLevelSetNormals(QuadNodes, i0, Len);
            MultidimensionalArray NodesGlobal = gridData.GlobalNodes.GetValue_Cell(QuadNodes, i0, Len);
            ParametersAndNormals.Stop();

            // Evaluate basis and test functions
            // ---------------------------------

            bool[] ReqV = new bool[GAMMA];
            bool[] ReqGradV = new bool[GAMMA];

            for(int gamma = 0; gamma < GAMMA; gamma++) {
                if (Koeff_V[gamma] != null) {
                    ReqV[gamma] = true;
                }
                if (Koeff_GradV[gamma] != null) {
                    ReqGradV[gamma] = true;
                }
            }
            
            MultidimensionalArray[] TestValues; //            index: codom variable/test function
            MultidimensionalArray[] TestGradientValues; //    index: codom variable/test function
            int[,] sectionsTest;
            Basis_Eval.Start();
            LECQuadratureLevelSet<IMutableMatrixEx,double[]>
                .EvalBasis(i0, Len, this.m_CodomainMap.BasisS, ReqV, ReqGradV, out TestValues, out TestGradientValues, out sectionsTest, QuadNodes);
            Basis_Eval.Stop();


            // Evaluate Integral components
            // ----------------------------

                       
            // loop over codomain variables ...
            for (int gamma = 0; gamma < GAMMA; gamma++) {

                // prepare parameters
                // - - - - - - - - - 

                // set Normal's
                LevSetIntParams _inParams = new LevSetIntParams();
                _inParams.Normal = Normals;
                // set Nodes Global
                _inParams.X = NodesGlobal;
                _inParams.time = this.time;
                _inParams.LsTrk = this.m_lsTrk;
                // set length scales


                // clear summation buffers
                // - - - - - - - - - - - -

                if (Koeff_V[gamma] != null)
                    Koeff_V[gamma].Clear();
                if (Koeff_GradV[gamma] != null)
                    Koeff_GradV[gamma].Clear();


                // Evaluate Bilin. forms
                // - - - - - - - - - - -
                
                {
                    EvalComponent(_inParams, gamma, this.m_NonlinLsForm_V[gamma], this.m_NonlinLsForm_V_Watches[gamma],
                        Koeff_V[gamma], 
                        m_FieldValuesPos, m_FieldValuesNeg, m_FieldGradientValuesPos, m_FieldGradientValuesNeg,
                        DELTA,
                        Flux_Eval,
                        delegate (INonlinLevelSetForm_V _comp, LevSetIntParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray SumBuf) {
                            _comp.LevelSetForm_V(_inParams, uA, uB, Grad_uA, Grad_uB, SumBuf);
                        });
                }
                {
                    EvalComponent(_inParams, gamma, this.m_NonlinLsForm_GradV[gamma], this.m_NonlinLsForm_GradV_Watches[gamma],
                        Koeff_GradV[gamma], 
                        m_FieldValuesPos, m_FieldValuesNeg, m_FieldGradientValuesPos, m_FieldGradientValuesNeg,
                        DELTA,
                        Flux_Eval,
                        delegate (INonlinLevelSetForm_GradV _comp, LevSetIntParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray SumBuf) {
                            _comp.LevelSetForm_GradV(_inParams, uA, uB, Grad_uA, Grad_uB, SumBuf);
                        });
                }
            }            

            // Summation Loops: multiply with test and trial functions
            // -------------------------------------------------------

            int[] offsetCod = new int[GAMMA];
            LECQuadratureLevelSet<IMutableMatrixEx,double[]>.
                CompOffsets(i0, Len, offsetCod, m_CodomainMap);

            for(int gamma = 0; gamma < GAMMA; gamma++) {
                // Evaluate Integrand
                // - - - - - - - - - 

                var TestVal = TestValues[gamma];
                var TestGradVal = TestGradientValues[gamma];
                int N;
                if(TestVal != null)
                    N = TestVal.GetLength(2);
                else if(TestGradVal != null)
                    N = TestGradVal.GetLength(2);
                else
                    N = 0;

                Loops.Start();



                // affine offset


                for(int cr = 0; cr < 2; cr++) { // loop over negative/positive species
                    int[] extr0 = new int[] { 0, 0, sectionsTest[gamma, cr] * N + offsetCod[gamma] };
                    int[] extrE = new int[] { Len - 1, NoOfNodes - 1, extr0[2] + N - 1 };
                    var SubRes = EvalResult.ExtractSubArrayShallow(extr0, extrE);

                    if(Koeff_V[gamma] != null) {
                        var Sum_Koeff_V_Cr = Koeff_V[gamma].ExtractSubArrayShallow(-1, -1, cr);
                        SubRes.Multiply(1.0, Sum_Koeff_V_Cr, TestVal, 1.0, "jkn", "jk", "jkn");
                    }
                    if(Koeff_GradV[gamma] != null) {
                        var Sum_Koeff_NablaV_Cr = Koeff_GradV[gamma].ExtractSubArrayShallow(-1, -1, cr, -1);
                        SubRes.Multiply(1.0, Sum_Koeff_NablaV_Cr, TestGradVal, 1.0, "jkn", "jkd", "jknd");
                    }
                }


                Loops.Stop();
            }
        }




        static private void EvalComponent<T>(LevSetIntParams _inParams,
            int gamma, EquationComponentArgMapping<T> bf, Stopwatch[] timers,
            MultidimensionalArray SumBuf,
            MultidimensionalArray[] FieldValuesPos, MultidimensionalArray[] FieldValuesNeg, MultidimensionalArray[] FieldGradientValuesPos, MultidimensionalArray[] FieldGradientValuesNeg,
            int DELTA,
            Stopwatch timer,
            Action<T, LevSetIntParams, MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray> ComponentFunc) where T : ILevelSetForm {
            timer.Start();


            for(int i = 0; i < bf.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components
                var comp = bf.m_AllComponentsOfMyType[i];


                int NoOfArgs = bf.NoOfArguments[i];
                Debug.Assert(NoOfArgs == comp.ArgumentOrdering.Count);
                int NoOfParams = bf.NoOfParameters[i];
                Debug.Assert(NoOfParams == ((comp.ParameterOrdering != null) ? comp.ParameterOrdering.Count : 0));

                // map arguments
                var uA = bf.MapArguments(FieldValuesNeg, comp, true);
                var uB = bf.MapArguments(FieldValuesPos, comp, true);
                var Grad_uA = bf.MapArguments(FieldGradientValuesNeg, comp, true);
                var Grad_uB = bf.MapArguments(FieldGradientValuesPos, comp, true);

                // map parameters
                _inParams.ParamsPos = new MultidimensionalArray[NoOfParams];
                _inParams.ParamsNeg = new MultidimensionalArray[NoOfParams];
                for(int c = 0; c < NoOfParams; c++) {
                    int targ = bf.AllToSub[i, c + NoOfArgs];
                    Debug.Assert(targ >= 0);
                    _inParams.ParamsPos[c] = FieldValuesPos[targ];
                    _inParams.ParamsNeg[c] = FieldValuesNeg[targ];
                }

                // evaluate equation components
                timers[i].Start();
                ComponentFunc(comp, _inParams, uA, uB, Grad_uA, Grad_uB, SumBuf);
                timers[i].Stop();
#if DEBUG
                SumBuf.CheckForNanOrInf();
#endif

            }
            timer.Stop();
        }
       


        /// <summary>
        /// writes the dammed result of the integration to the sparse matrix
        /// </summary>
        protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {

            int Nmax = m_CodomainMap.MaxTotalNoOfCoordinatesPerCell;

            //int[] _i0 = new int[] { int.MaxValue, 0, 1 };
            //int[] _iE = new int[] { -1, Nmax - 1, _i0[2] + Mmax - 1 };

            int[] _i0aff = new int[] { int.MaxValue, 0 };
            int[] _iEaff = new int[] { -1, Nmax - 1 };


            // loop over cells...
            for(int i = 0; i < Length; i++) {
                int jCell = i + i0;
                int Row0 = m_CodomainMap.LocalUniqueCoordinateIndex(0, jCell, 0);
                int Row0_g = m_CodomainMap.i0 + Row0;

                _i0aff[0] = i;
                _iEaff[0] = -1;

                var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(_i0aff, _iEaff);

                for(int r = BlockRes.GetLength(0) - 1; r >= 0; r--)
                    ResultVector[Row0 + r] += BlockRes[r];

            }
        }
    }
}
