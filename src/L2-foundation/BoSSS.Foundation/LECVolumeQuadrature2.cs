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
using System.Text;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Foundation.Quadrature.Linear {
    class LECVolumeQuadrature2<M, V>
        where M : IMutableMatrix
        where V : IList<double> //
    {
        public LECVolumeQuadrature2(SpatialOperator op) {
            Operator = op;
            m_VolumeForm_UxV = op.GetArgMapping<IVolumeForm_UxV>(true, eq => ((eq.VolTerms & TermActivationFlags.UxV) != 0), eq => (eq is IVolumeForm) ? new LinearVolumeFormVectorizer((IVolumeForm)eq) : null);
            m_VolumeForm_UxGradV = op.GetArgMapping<IVolumeForm_UxGradV>(true, eq => ((eq.VolTerms & TermActivationFlags.UxGradV) != 0), eq => (eq is IVolumeForm) ? new LinearVolumeFormVectorizer((IVolumeForm)eq) : null);
            m_VolumeForm_GradUxV = op.GetArgMapping<IVolumeForm_GradUxV>(true, eq => ((eq.VolTerms & TermActivationFlags.GradUxV) != 0), eq => (eq is IVolumeForm) ? new LinearVolumeFormVectorizer((IVolumeForm)eq) : null);
            m_VolumeForm_GradUxGradV = op.GetArgMapping<IVolumeForm_GradUxGradV>(true, eq => ((eq.VolTerms & TermActivationFlags.GradUxGradV) != 0), eq => (eq is IVolumeForm) ? new LinearVolumeFormVectorizer((IVolumeForm)eq) : null);
            m_VolumeSource_V = op.GetArgMapping<IVolumeSource_V>(true, eq => ((eq.VolTerms & TermActivationFlags.V) != 0), eq => (eq is IVolumeForm) ? new LinearVolumeFormVectorizer((IVolumeForm)eq) : null);
            m_VolumeSource_GradV = op.GetArgMapping<IVolumeSource_GradV>(true, eq => ((eq.VolTerms & TermActivationFlags.GradV) != 0), eq => (eq is IVolumeForm) ? new LinearVolumeFormVectorizer((IVolumeForm)eq) : null);
        }

        SpatialOperator Operator;


        /// <summary>
        /// number of codomain variables
        /// </summary>
        int GAMMA {
            get {
                return Operator.CodomainVar.Count;
            }
        }

        /// <summary>
        /// number of domain variables
        /// </summary>
        int DELTA {
            get {
                return Operator.DomainVar.Count;
            }
        }

        /// <summary>
        /// index: test function/codomain variable
        /// </summary>
        EquationComponentArgMapping<IVolumeSource_GradV>[] m_VolumeSource_GradV;

        /// <summary>
        /// index: test function/codomain variable
        /// </summary>
        EquationComponentArgMapping<IVolumeSource_V>[] m_VolumeSource_V;


        /// <summary>
        /// index: test function/codomain variable
        /// </summary>
        EquationComponentArgMapping<IVolumeForm_GradUxGradV>[] m_VolumeForm_GradUxGradV;

        /// <summary>
        /// index: test function/codomain variable
        /// </summary>
        EquationComponentArgMapping<IVolumeForm_GradUxV>[] m_VolumeForm_GradUxV;

        /// <summary>
        /// index: test function/codomain variable
        /// </summary>
        EquationComponentArgMapping<IVolumeForm_UxGradV>[] m_VolumeForm_UxGradV;

        /// <summary>
        /// index: test function/codomain variable
        /// </summary>
        EquationComponentArgMapping<IVolumeForm_UxV>[] m_VolumeForm_UxV;


        Stopwatch[][] m_VolumeForm_UxV_Watches;
        Stopwatch[][] m_VolumeForm_GradUxV_Watches;
        Stopwatch[][] m_VolumeForm_UxGradV_Watches;
        Stopwatch[][] m_VolumeForm_GradUxGradV_Watches;
        Stopwatch[][] m_VolumeSource_V_Watches;
        Stopwatch[][] m_VolumeSource_GradV_Watches;


        /// <summary>
        /// Performs the integration of edge components.
        /// </summary>
        /// <param name="domNrule"></param>
        /// <param name="RowMap"></param>
        /// <param name="ParamsMap"></param>
        /// <param name="ColMap">
        /// Domain variables resp. trial variables resp. matrix column variables
        /// </param>
        /// <param name="Matrix">output accumulator: matrix where the linear part of the operator is stored.</param>
        /// <param name="Vector">output accumulator: vector where the affine part of the operator is stored.</param>
        /// <param name="time"></param>
        public void Execute(ICompositeQuadRule<QuadRule> domNrule,
            UnsetteledCoordinateMapping RowMap, IList<DGField> ParamsMap, UnsetteledCoordinateMapping ColMap,
            M Matrix, V Vector, double time) {       
            if(RowMap.BasisS.Count != GAMMA)
                throw new ArgumentException("Mismatch in number of codomain (rew. row-variables, resp. test-variables) variables.", "RowMap");
            if(ColMap.BasisS.Count != DELTA)
                throw new ArgumentException("Mismatch in number of domain (rew. column-variables, resp. trial-variables) variables.", "ColMap");

            m_GridDat = RowMap.GridDat;

            if (!object.ReferenceEquals(m_GridDat, ColMap.GridDat))
                throw new ArgumentException();

            m_ParameterFields = ParamsMap == null ? new DGField[0] : ParamsMap.ToArray();
            
            m_RowL = RowMap.BasisS.Select(b => b.Length).ToArray();
            m_ColL = ColMap.BasisS.Select(b => b.Length).ToArray();
            m_Vfunctions = RowMap.BasisS.ToArray();
            m_Ufunctions = ColMap.BasisS.ToArray();
            m_RowMap = RowMap;
            m_ColMap = ColMap;
            
            m_time = time;
            m_Matrix = Matrix;
            m_Vector = Vector;
            if(m_Matrix != null) {
                if(Matrix.RowPartitioning.LocalLength != RowMap.LocalLength)
                    throw new ArgumentException("Matrix mismatch: local number of rows", "Matrix");
                if(Matrix.ColPartition.LocalLength != ColMap.LocalLength)
                    throw new ArgumentException("Matrix mismatch: number of columns", "Matrix");
            }
            if(m_Vector != null) {
                if(m_Vector.Count != RowMap.LocalLength) {
                    throw new ArgumentException("affine vector mismatch: local number of rows", "AffineVector");
                }
            }

            var q = CellQuadrature.GetQuadrature2(new int[] { m_RowL.Sum(), (this.LinearRequired ? m_ColL.Sum() : 0) + (this.AffineRequired ? 1 : 0) },
                m_GridDat, domNrule,
                this.EvaluateEx,
                this.SaveIntegrationResults,
                _AllocateBuffers:this.AllocateBuffers);

            q.CustomTimers = new Stopwatch[] { new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch() };
            q.CustomTimers_Names = new string[] { "Flux-Eval", "Basis-Eval", "Loops", "ParametersAndNormals", "Flux-Trafo" };
            q.CustomTimers_RootPointer = new int[q.CustomTimers_Names.Length];
            ArrayTools.SetAll(q.CustomTimers_RootPointer, -1);
            this.Flux_Eval = q.CustomTimers[0];
            this.Basis_Eval = q.CustomTimers[1];
            this.Loops = q.CustomTimers[2];
            this.ParametersAndNormals = q.CustomTimers[3];
            this.FluxTrafo = q.CustomTimers[4];

            Debug.Assert(Array.IndexOf(q.CustomTimers_Names, "Flux-Eval") == 0);
            this.m_VolumeForm_UxV_Watches = this.m_VolumeForm_UxV.InitStopWatches(0, q);
            this.m_VolumeForm_UxGradV_Watches = this.m_VolumeForm_UxGradV.InitStopWatches(0, q);
            this.m_VolumeForm_GradUxV_Watches = this.m_VolumeForm_GradUxV.InitStopWatches(0, q);
            this.m_VolumeForm_GradUxGradV_Watches = this.m_VolumeForm_GradUxGradV.InitStopWatches(0, q);

            this.m_VolumeSource_V_Watches = this.m_VolumeSource_V.InitStopWatches(0, q);
            this.m_VolumeSource_GradV_Watches = this.m_VolumeSource_GradV.InitStopWatches(0, q);

            q.Execute();

        }

        double m_time;

        Stopwatch Flux_Eval;
        Stopwatch Basis_Eval;
        Stopwatch Loops;
        Stopwatch ParametersAndNormals;
        Stopwatch FluxTrafo;

        UnsetteledCoordinateMapping m_RowMap;
        UnsetteledCoordinateMapping m_ColMap;

        /// <summary>
        /// parameter fields
        /// </summary>
        DGField[] m_ParameterFields;

        /// <summary>
        /// length of codomain basis per cell
        /// </summary>
        int[] m_RowL;

        /// <summary>
        /// 
        /// </summary>
        Basis[] m_Vfunctions;

        /// <summary>
        /// length of domain basis per cell
        /// </summary>
        int[] m_ColL;

        /// <summary>
        /// 
        /// </summary>
        Basis[] m_Ufunctions;



        IGridData m_GridDat;

        /// <summary>
        /// Matrix where the linear part of the operator is saved.
        /// </summary>
        M m_Matrix;

        /// <summary>
        /// Vector where the affine part of the operator is saved.
        /// </summary>
        V m_Vector;

        /// <summary>
        /// true, if integration of <see cref="m_Matrix"/> is required.
        /// </summary>
        bool LinearRequired {
            get {
                if(m_Matrix == null)
                    return false;
                for(int gamma = this.GAMMA - 1; gamma >= 0; gamma--) {
                    if(m_VolumeForm_UxV[gamma].m_AllComponentsOfMyType.Length > 0)
                        return true;
                    if(m_VolumeForm_UxGradV[gamma].m_AllComponentsOfMyType.Length > 0)
                        return true;
                    if(m_VolumeForm_GradUxV[gamma].m_AllComponentsOfMyType.Length > 0)
                        return true;
                    if(m_VolumeForm_GradUxGradV[gamma].m_AllComponentsOfMyType.Length > 0)
                        return true;
                }
                return false;
            }
        }

        /// <summary>
        /// true, if integration of <see cref="m_Vector"/> is required.
        /// </summary>
        bool AffineRequired {
            get {
                if(m_Vector == null)
                    return false;
                for(int gamma = this.GAMMA - 1; gamma >= 0; gamma--) {
                    if(m_VolumeSource_V[gamma].m_AllComponentsOfMyType.Length > 0)
                        return true;
                    if(m_VolumeSource_GradV[gamma].m_AllComponentsOfMyType.Length > 0)
                        return true;
                }
                return false;
            }
        }


        
        MultidimensionalArray[] m_ParameterFieldsValues;
        
        // 1st index (staggered): test/codomain variable
        // 2nd index (staggered): equation component
        MultidimensionalArray[][] m_GradUxGradVComponentBuffer;
        MultidimensionalArray[][] m_UxGradVComponentBuffer;
        MultidimensionalArray[][] m_GradUxVComponentBuffer;
        MultidimensionalArray[][] m_UxVComponentBuffer;
        //MultidimensionalArray[][] m_VComponentBuffer;
        //MultidimensionalArray[][] m_GradVComponentBuffer;

        // 1st index: test/codomain variable
        // 2nd index: trial/domain component
        MultidimensionalArray[,] m_GradUxGradVSumBuffer;
        MultidimensionalArray[,] m_UxGradVSumBuffer;
        MultidimensionalArray[,] m_GradUxVSumBuffer;
        MultidimensionalArray[,] m_UxVSumBuffer;
        MultidimensionalArray[,] m_VSumBuffer;     // affine part: second index only 0
        MultidimensionalArray[,] m_GradVSumBuffer; // affine part: second index only 0


        // 1st index: test/codomain variable
        // 2nd index: trial/domain component
        MultidimensionalArray[,] m_Trf_GradUxGradVSumBuffer;
        MultidimensionalArray[,] m_Trf_UxGradVSumBuffer;
        MultidimensionalArray[,] m_Trf_GradUxVSumBuffer;
        MultidimensionalArray[,] m_Trf_UxVSumBuffer;
        MultidimensionalArray[,] m_Trf_VSumBuffer;     // affine part: second index only 0
        MultidimensionalArray[,] m_Trf_GradVSumBuffer; // affine part: second index only 0


        void TestFunctionRequired(bool[] Req, MultidimensionalArray[,] SumBuf, ref int MaxDeg, ref int MaxN) {
            Debug.Assert(Req.Length == SumBuf.GetLength(0));

            for(int i = SumBuf.GetLength(0) - 1; i >= 0; i--) {
                for(int j = SumBuf.GetLength(1) - 1; j >= 0; j--) {
                    if(SumBuf[i, j] != null) {
                        MaxDeg = Math.Max(this.m_Vfunctions[i].Degree, MaxDeg);
                        MaxN = Math.Max(MaxN, this.m_RowL[i]);
                        Req[i] = true;
                        break;
                    }
                }
            }
        }

        void TrialFunctionRequired(bool[] Req, MultidimensionalArray[,] SumBuf, bool LinReq, ref int MaxDeg, ref int MaxN) {
            Debug.Assert(Req.Length == SumBuf.GetLength(1));
            if(LinReq) {
                for(int i = SumBuf.GetLength(1) - 1; i >= 0; i--) {
                    for(int j = SumBuf.GetLength(0) - 1; j >= 0; j--) {
                        if(SumBuf[j, i] != null) {
                            MaxDeg = Math.Max(this.m_Ufunctions[i].Degree, MaxDeg);
                            MaxN = Math.Max(MaxN, this.m_ColL[i]);
                            Req[i] = true;
                            break;
                        }
                    }
                }
            } else {
                Req.SetAll(false);
            }
        }

        protected void AllocateBuffers(int NoOfItems, MultidimensionalArray rule) {
            int NoOfNodes = rule.GetLength(0);
            int D = m_GridDat.SpatialDimension;

            if (m_ParameterFieldsValues == null) {
                m_ParameterFieldsValues = m_ParameterFields.Select(f => new MultidimensionalArray(2)).ToArray();
            }

            foreach (var pfv in m_ParameterFieldsValues) {
                if (pfv.GetLength(0) != NoOfItems || pfv.GetLength(1) != NoOfNodes)
                    pfv.Allocate(NoOfItems, NoOfNodes);
            }

            AllocateComponentBuffer(ref m_UxVComponentBuffer, m_VolumeForm_UxV, NoOfItems, NoOfNodes, D, 3, false);
            AllocateComponentBuffer(ref m_UxGradVComponentBuffer, m_VolumeForm_UxGradV, NoOfItems, NoOfNodes, D, 4, false);
            AllocateComponentBuffer(ref m_GradUxVComponentBuffer, m_VolumeForm_GradUxV, NoOfItems, NoOfNodes, D, 4, false);
            AllocateComponentBuffer(ref m_GradUxGradVComponentBuffer, m_VolumeForm_GradUxGradV, NoOfItems, NoOfNodes, D, 5, false);
            //AllocateComponentBuffer(ref m_VComponentBuffer, m_VolumeSource_V, NoOfItems, NoOfNodes, D, 2, true);
            //AllocateComponentBuffer(ref m_GradVComponentBuffer, m_VolumeSource_GradV, NoOfItems, NoOfNodes, D, 3, true);
            
            AllocateSumBuffer(ref m_UxVSumBuffer, m_VolumeForm_UxV, NoOfItems, NoOfNodes, D, 2, false);
            AllocateSumBuffer(ref m_UxGradVSumBuffer, m_VolumeForm_UxGradV, NoOfItems, NoOfNodes, D, 3, false);
            AllocateSumBuffer(ref m_GradUxVSumBuffer, m_VolumeForm_GradUxV, NoOfItems, NoOfNodes, D, 3, false);
            AllocateSumBuffer(ref m_GradUxGradVSumBuffer, m_VolumeForm_GradUxGradV, NoOfItems, NoOfNodes, D, 4, false);
            AllocateSumBuffer(ref m_VSumBuffer, m_VolumeSource_V, NoOfItems, NoOfNodes, D, 2, true);
            AllocateSumBuffer(ref m_GradVSumBuffer, m_VolumeSource_GradV, NoOfItems, NoOfNodes, D, 3, true);

            m_Trf_UxVSumBuffer = m_UxVSumBuffer;  // no transformation required for UxV-components
            AllocateSumBuffer(ref m_Trf_UxGradVSumBuffer, m_VolumeForm_UxGradV, NoOfItems, NoOfNodes, D, 3, false);
            AllocateSumBuffer(ref m_Trf_GradUxVSumBuffer, m_VolumeForm_GradUxV, NoOfItems, NoOfNodes, D, 3, false);
            AllocateSumBuffer(ref m_Trf_GradUxGradVSumBuffer, m_VolumeForm_GradUxGradV, NoOfItems, NoOfNodes, D, 4, false);
            m_Trf_VSumBuffer = m_VSumBuffer;      // no transformation required for V-components
            AllocateSumBuffer(ref m_Trf_GradVSumBuffer, m_VolumeSource_GradV, NoOfItems, NoOfNodes, D, 3, true);
        }

        private void AllocateSumBuffer<Ee>(ref MultidimensionalArray[,] SumBuffer, EquationComponentArgMapping<Ee>[] FormMapping, int NoOfItems, int NoOfNodes, int D, int ArrayDim, bool affine)
            where Ee : IEquationComponent 
        {
            if (SumBuffer == null) {
                SumBuffer = new MultidimensionalArray[GAMMA, affine ? 1 : DELTA];
            }

            int[] LenToAlloc = new int[ArrayDim];
            LenToAlloc[0] = NoOfItems;
            LenToAlloc[1] = NoOfNodes;
            for(int k = 2; k < ArrayDim; k++) {
                LenToAlloc[k] = D;
            }

            bool[] UsedMarker = new bool[DELTA];
            for (int gamma = this.GAMMA - 1; gamma >= 0; gamma--) {
                var ecq = FormMapping[gamma];

                UsedMarker.Clear();
                for (int i = 0; i < ecq.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components
                    var comp = ecq.m_AllComponentsOfMyType[i];
                    int NoOfArgs = comp.ArgumentOrdering.Count;
                    if(affine) {
                        UsedMarker[0] = true;
                    } else {
                        for(int j = 0; j < NoOfArgs; j++) {
                            int targ = ecq.AllToSub[i, j];
                            UsedMarker[targ] = true;
                        }
                    }
                }

                for (int delta = affine ? 0 : this.DELTA - 1; delta >= 0; delta--) {
                    if (UsedMarker[delta]) {
                        if (SumBuffer[gamma, delta] == null) {
                            SumBuffer[gamma, delta] = new MultidimensionalArray(ArrayDim);
                        }

                        var mda = SumBuffer[gamma,delta];

                        if (mda.GetLength(0) != NoOfItems || mda.GetLength(1) != NoOfNodes) {
                            mda.Allocate(LenToAlloc);
                        }
                    }
                }
            }
        }


        private void AllocateComponentBuffer<Ee>(ref MultidimensionalArray[][] CompBuffer, EquationComponentArgMapping<Ee>[] FormMapping, int NoOfItems, int NoOfNodes, int D, int ArrayDim, bool affine)
            where Ee: IEquationComponent
        {
            if (CompBuffer == null) {
                CompBuffer = FormMapping.Select(E => E.m_AllComponentsOfMyType.Select(C => new MultidimensionalArray(ArrayDim)).ToArray()).ToArray();
            }

            for (int gamma = this.GAMMA - 1; gamma >= 0; gamma--) {
                var ecq = FormMapping[gamma];

                for (int i = 0; i < ecq.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components

                    var comp = ecq.m_AllComponentsOfMyType[i];
                    int C = comp.ArgumentOrdering.Count();

                    var mda = CompBuffer[gamma][i];
                    int Dim = mda.Dimension;
                    if (mda.GetLength(0) != NoOfItems || mda.GetLength(1) != NoOfNodes) {
                        int[] LenToAlloc = new int[ArrayDim];
                        LenToAlloc[0] = NoOfItems;
                        LenToAlloc[1] = NoOfNodes;
                        if(affine) {
                            for(int k = 2; k < ArrayDim; k++)
                                LenToAlloc[k] = D;
                        } else {
                            LenToAlloc[2] = C;
                            for(int k = 3; k < ArrayDim; k++)
                                LenToAlloc[k] = D;
                        }

                        mda.Allocate(LenToAlloc);
                    }
                }
            }
        }


        protected void EvaluateEx(int i0, int Length, QuadRule qr, MultidimensionalArray QuadResult) {

            // intitial checks
            // ===============
            NodeSet qrNodes = qr.Nodes;
            int NoOfNodes = qr.NoOfNodes;
            bool Affine;
            int D = this.m_GridDat.SpatialDimension;
            bool bLinearRequired = this.LinearRequired;
            bool bAffineRequired = this.AffineRequired; // 
            int[] geom2log = m_GridDat.iGeomCells.GeomCell2LogicalCell;
            {

                Affine = m_GridDat.iGeomCells.IsCellAffineLinear(i0);
                
#if DEBUG
                int iKref = m_GridDat.iGeomCells.GetRefElementIndex(i0);


                for(int j = 0; j < Length; j++) {
                    int jCell = j + i0;
                    Debug.Assert(m_GridDat.iGeomCells.GetRefElementIndex(jCell) == iKref);
                    Debug.Assert(m_GridDat.iGeomCells.IsCellAffineLinear(jCell) == Affine);

                    for(int gamma = 0; gamma < GAMMA; gamma++)
                        Debug.Assert(this.m_RowL[gamma] == m_Vfunctions[gamma].GetLength(geom2log != null ? geom2log[jCell] : jCell));

                    for(int delta = 0; delta < DELTA; delta++)
                        Debug.Assert(this.m_ColL[delta] == m_Ufunctions[delta].GetLength(geom2log != null ? geom2log[jCell] : jCell));
                }
#endif
            }


            // evaluate parameters
            // ===================
            #region PARAMETERS
            this.ParametersAndNormals.Start();
            MultidimensionalArray globalNodes;
            {
                Debug.Assert(m_ParameterFieldsValues.Length == m_ParameterFields.Length);
                for(int i = 0; i < m_ParameterFields.Length; i++) {
                    if(m_ParameterFields[i] != null)
                        m_ParameterFields[i].Evaluate(i0, Length, qrNodes, m_ParameterFieldsValues[i]);
                    else
                        m_ParameterFieldsValues[i].Clear();
                }
                globalNodes = this.m_GridDat.GlobalNodes.GetValue_Cell(qrNodes, i0, Length);

            }
            this.ParametersAndNormals.Stop();
            #endregion

            // Evaluate test and trial space basis
            // ===================================
            #region BASIS_EVAL
            this.Basis_Eval.Start();
            MultidimensionalArray[] V_XquadWgt = new MultidimensionalArray[GAMMA]; // test function values, multiplied by quadrature weights
            MultidimensionalArray[] GradV_XquadWgt = new MultidimensionalArray[GAMMA]; // test function gradients, multiplied by quadrature weights
            MultidimensionalArray[] U = new MultidimensionalArray[DELTA];       // trial function values 
            MultidimensionalArray[] GradU = new MultidimensionalArray[DELTA];       // trial function gradients
            int NZ3 = -1, NZ4 = -1, MR = -1, NR = -1; // dimensions for intermediate result buffers
            int maxDeg = -1;
            bool JacobiRequired;
            {
                bool[] TestReq = new bool[GAMMA];
                bool[] TestGradientreq = new bool[GAMMA];
                bool[] TrialReq = new bool[DELTA];
                bool[] TrialGradientreq = new bool[DELTA];

                // maximum required degree for basis values/basis gradient values...
                int MaxDeg_Test = -1, MaxDeg_TestGrad = -1;
                int MaxDeg_Trial = -1, MaxDeg_TrialGrad = -1;
                int MaxN_Test = -1, MaxN_TestGrad = -1;
                int MaxN_Trial = -1, MaxN_TrialGrad = -1;

                // detect which test and trial functions are required....
                if(bLinearRequired || bAffineRequired) {
                    TestFunctionRequired(TestReq, m_UxVSumBuffer, ref MaxDeg_Test, ref MaxN_Test);
                    TestFunctionRequired(TestReq, m_GradUxVSumBuffer, ref MaxDeg_Test, ref MaxN_Test);
                    TestFunctionRequired(TestGradientreq, m_UxGradVSumBuffer, ref MaxDeg_TestGrad, ref MaxN_TestGrad);
                    TestFunctionRequired(TestGradientreq, m_GradUxGradVSumBuffer, ref MaxDeg_TestGrad, ref MaxN_TestGrad);
                    TestFunctionRequired(TestReq, m_VSumBuffer, ref MaxDeg_Test, ref MaxN_Test);
                    TestFunctionRequired(TestGradientreq, m_GradVSumBuffer, ref MaxDeg_TestGrad, ref MaxN_TestGrad);
                }

                TrialFunctionRequired(TrialReq, m_UxVSumBuffer, bLinearRequired, ref MaxDeg_Trial, ref MaxN_Trial);
                TrialFunctionRequired(TrialReq, m_UxGradVSumBuffer, bLinearRequired, ref MaxDeg_Trial, ref MaxN_Trial);
                TrialFunctionRequired(TrialGradientreq, m_GradUxVSumBuffer, bLinearRequired, ref MaxDeg_TrialGrad, ref MaxN_TrialGrad);
                TrialFunctionRequired(TrialGradientreq, m_GradUxGradVSumBuffer, bLinearRequired, ref MaxDeg_TrialGrad, ref MaxN_TrialGrad);

                // evaluate
                var MasterBasis = m_GridDat.ChefBasis;
                MultidimensionalArray Values = null, Grad = null;
                int MaxDeg_Values = Math.Max(MaxDeg_Test, MaxDeg_Trial);
                int MaxDeg_Grad = Math.Max(MaxDeg_TestGrad, MaxDeg_TrialGrad);
                if(MaxDeg_Values >= 0)
                    Values = MasterBasis.BasisValues.GetValues(qrNodes, MaxDeg_Values);
                if(MaxDeg_Grad >= 0)
                    Grad = MasterBasis.BasisGradientValues.GetValues(qrNodes, MaxDeg_Grad);

                // multiply test functions with quadrature weights
                MultidimensionalArray ValuesXquadWgt = null, GradsXquadWgt = null;
                if(MaxN_Test >= 0) {
                    ValuesXquadWgt = MultidimensionalArray.Create(NoOfNodes, MaxN_Test);
                    ValuesXquadWgt.Multiply(1.0, qr.Weights, Values.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NoOfNodes - 1, MaxN_Test - 1 }), 0.0, ref mp_kn_k_kn);
                }
                if(MaxN_TestGrad >= 0) {
                    GradsXquadWgt = MultidimensionalArray.Create(NoOfNodes, MaxN_TestGrad, D);
                    GradsXquadWgt.Multiply(1.0, qr.Weights, Grad.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { NoOfNodes - 1, MaxN_TestGrad - 1, D - 1 }), 0.0, ref mp_knd_k_knd);
                }


                for(int gamma = 0; gamma < GAMMA; gamma++) {
                    Debug.Assert(object.ReferenceEquals(m_Vfunctions[gamma].GridDat.ChefBasis, MasterBasis), "We assume that all test functions redict to the same master basis");

                    if(TestReq[gamma])
                        V_XquadWgt[gamma] = ValuesXquadWgt.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NoOfNodes - 1, this.m_RowL[gamma] - 1 });
                    if(TestGradientreq[gamma])
                        GradV_XquadWgt[gamma] = GradsXquadWgt.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { NoOfNodes - 1, this.m_RowL[gamma] - 1, D - 1 });
                }
                for(int delta = 0; delta < DELTA; delta++) {
                    Debug.Assert(object.ReferenceEquals(m_Ufunctions[delta].GridDat.ChefBasis, MasterBasis), "We assume that all trial functions redict to the same master basis");

                    if(TrialReq[delta])
                        U[delta] = Values.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NoOfNodes - 1, this.m_ColL[delta] - 1 });
                    if(TrialGradientreq[delta])
                        GradU[delta] = Grad.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { NoOfNodes - 1, this.m_ColL[delta] - 1, D - 1 });
                }

                NZ3 = Math.Max(MaxN_Test, Math.Max(MaxN_TestGrad, MaxN_TrialGrad));
                NZ4 = MaxN_TrialGrad;
                MR = Math.Max(MaxN_Test, MaxN_TestGrad);
                NR = Math.Max(MaxN_Trial, MaxN_TrialGrad);
                maxDeg = Math.Max(MaxDeg_Test, Math.Max(MaxDeg_TestGrad, Math.Max(MaxDeg_Trial, MaxDeg_TrialGrad)));
                JacobiRequired = (MaxN_TestGrad >= 0) || (MaxDeg_TrialGrad >= 0);
            }
            this.Basis_Eval.Stop();
            #endregion

            // evaluate forms/fluxes
            // =====================
            #region FLUXEVAL
            this.Flux_Eval.Start();
            {
                VolumFormParams vfp = default(VolumFormParams);
                vfp.j0 = i0;
                vfp.Len = Length;
                vfp.GridDat = m_GridDat;
                vfp.Xglobal = globalNodes;
                vfp.time = this.m_time;
                
                if(bLinearRequired) {
                    EvalNSumForm(ref vfp, this.m_VolumeForm_UxV, m_UxVComponentBuffer, m_UxVSumBuffer, (C, M) => C.Form(ref vfp, M), false, this.m_VolumeForm_UxV_Watches);
                    EvalNSumForm(ref vfp, this.m_VolumeForm_UxGradV, m_UxGradVComponentBuffer, m_UxGradVSumBuffer, (C, M) => C.Form(ref vfp, M), false, this.m_VolumeForm_UxGradV_Watches);
                    EvalNSumForm(ref vfp, this.m_VolumeForm_GradUxV, m_GradUxVComponentBuffer, m_GradUxVSumBuffer, (C, M) => C.Form(ref vfp, M), false, this.m_VolumeForm_GradUxV_Watches);
                    EvalNSumForm(ref vfp, this.m_VolumeForm_GradUxGradV, m_GradUxGradVComponentBuffer, m_GradUxGradVSumBuffer, (C, M) => C.Form(ref vfp, M), false, this.m_VolumeForm_GradUxGradV_Watches);
                }
                if(bAffineRequired) {
                    EvalNSumForm(ref vfp, this.m_VolumeSource_V, null, m_VSumBuffer, (C, M) => C.Form(ref vfp, M), true, this.m_VolumeSource_V_Watches);
                    EvalNSumForm(ref vfp, this.m_VolumeSource_GradV, null, m_GradVSumBuffer, (C, M) => C.Form(ref vfp, M), true, this.m_VolumeSource_GradV_Watches);
                }
            }
            this.Flux_Eval.Stop();
            #endregion

            // transform forms/fluxes
            // ======================
            #region FLUX_TRAFO
            this.FluxTrafo.Start();
            {
                MultidimensionalArray invJac; // inverse Jacobi matrix -- transformation of gradients between reference and physical space
                MultidimensionalArray JacDet; // integral transformation metric betewwn reference and physical cell

                if(Affine) {
                    invJac = JacobiRequired ? this.m_GridDat.iGeomCells.InverseTransformation.ExtractSubArrayShallow(new int[] { i0, 0, 0 }, new int[] { i0 + Length - 1, D - 1, D - 1 }) : null;
                    JacDet = null; // for affine-linear cells, this is constant per cell, so it can be applied after quadrature sum (see below)
                } else {
                    invJac = JacobiRequired ? this.m_GridDat.InverseJacobian.GetValue_Cell(qrNodes, i0, Length) : null;
                    JacDet = this.m_GridDat.JacobianDeterminat.GetValue_Cell(qrNodes, i0, Length);
                }
                

                // linear forms
                if(bLinearRequired) {
                    for(int gamma = 0; gamma < GAMMA; gamma++) {
                        for(int delta = 0; delta < DELTA; delta++) {
                            // UXV components 
                            if(m_UxVSumBuffer[gamma, delta] != null) {
                                Debug.Assert(object.ReferenceEquals(m_Trf_UxVSumBuffer[gamma, delta], m_UxVSumBuffer[gamma, delta]));
                                if(!Affine)
                                    m_Trf_UxVSumBuffer[gamma, delta].Multiply(1.0, JacDet, m_Trf_UxVSumBuffer[gamma, delta], 0.0, "jk", "jk", "jk");
                            }

                            // GradUxGradV components
                            if(m_GradUxGradVSumBuffer[gamma, delta] != null) {
                                MultidimensionalArray.MultiplyProgram mp_prog1 = (Affine ? mp_ikae_iad_ikde : mp_ikae_ikad_ikde);
                                m_Trf_GradUxGradVSumBuffer[gamma, delta].Multiply(1.0, invJac, m_GradUxGradVSumBuffer[gamma, delta], 0.0, ref mp_prog1);

                                // swap the buffers for the second transform ...
                                MultidimensionalArray tmp = m_GradUxGradVSumBuffer[gamma, delta];
                                m_GradUxGradVSumBuffer[gamma, delta] = m_Trf_GradUxGradVSumBuffer[gamma, delta];
                                m_Trf_GradUxGradVSumBuffer[gamma, delta] = tmp;

                                MultidimensionalArray.MultiplyProgram mp_prog2 = (Affine ? mp_ikab_ikae_ibe : mp_ikab_ikae_ikbe);
                                m_Trf_GradUxGradVSumBuffer[gamma, delta].Multiply(1.0, m_GradUxGradVSumBuffer[gamma, delta], invJac, 0.0, ref mp_prog2);

                                if(!Affine)
                                    m_Trf_GradUxGradVSumBuffer[gamma, delta].Multiply(1.0, JacDet, m_Trf_GradUxGradVSumBuffer[gamma, delta], 0.0, "jked", "jk", "jked");


                            }

                            // GradUxV components
                            if(m_GradUxVSumBuffer[gamma, delta] != null) {
                                MultidimensionalArray.MultiplyProgram mp_prog = (Affine ? mp_ikb_ikd_ibd : mp_ikb_ikd_ikbd);
                                m_Trf_GradUxVSumBuffer[gamma, delta].Multiply(1.0, m_GradUxVSumBuffer[gamma, delta], invJac, 0.0, ref mp_prog);
                                if(!Affine)
                                    m_Trf_GradUxVSumBuffer[gamma, delta].Multiply(1.0, JacDet, m_Trf_GradUxVSumBuffer[gamma, delta], 0.0, "jkd", "jk", "jkd");
                            }

                            // UxGradV components
                            if(m_UxGradVSumBuffer[gamma, delta] != null) {
                                MultidimensionalArray.MultiplyProgram mp_prog = (Affine ? mp_ika_iad_ikd : mp_ika_ikad_ikd);
                                m_Trf_UxGradVSumBuffer[gamma, delta].Multiply(1.0, invJac, m_UxGradVSumBuffer[gamma, delta], 0.0, ref mp_prog);
                                if(!Affine)
                                    m_Trf_UxGradVSumBuffer[gamma, delta].Multiply(1.0, JacDet, m_Trf_UxGradVSumBuffer[gamma, delta], 0.0, "jkd", "jk", "jkd");
                            }
                        }
                    }
                }

                // affine forms
                if(bAffineRequired) {
                    for(int gamma = 0; gamma < GAMMA; gamma++) {
                        // v components
                        if(m_VSumBuffer[gamma, 0] != null) {
                            Debug.Assert(object.ReferenceEquals(m_Trf_VSumBuffer[gamma, 0], m_VSumBuffer[gamma, 0]));
                            if(!Affine)
                                m_Trf_VSumBuffer[gamma, 0].Multiply(1.0, JacDet, m_Trf_VSumBuffer[gamma, 0], 0.0, "jk", "jk", "jk");
                        }

                        // Gradv components
                        if(m_GradVSumBuffer[gamma, 0] != null) {
                            MultidimensionalArray.MultiplyProgram mp_prog = (Affine ? mp_ika_iad_ikd : mp_ika_ikad_ikd);
                            m_Trf_GradVSumBuffer[gamma, 0].Multiply(1.0, invJac, m_GradVSumBuffer[gamma, 0], 0.0, ref mp_prog);
                                
                            if(!Affine)
                                m_Trf_GradVSumBuffer[gamma, 0].Multiply(1.0, JacDet, m_Trf_GradVSumBuffer[gamma, 0], 0.0, "jkd", "jk", "jkd");
                        }
                    }
                }

            }
            this.FluxTrafo.Stop();
#endregion

            // summation
            // =========
            #region SUMMATION_LOOPS
            this.Loops.Start();
            {
                int iBufZ3 = 0, iBufZ4 = 0, iTrafoBuf = 0, iQlBuf = 0, iQaBuf = 0;
                MultidimensionalArray Z3 = null, Z4 = null, _Z, trafo = null, _R, _Q, Qa = null, Ql = null;
                if(Affine) {
#if DEBUG
                    for(int i = 0; i < Length; i++) {
                        int jCell = i + i0;

                        double detJac_jCell = m_GridDat.iGeomCells.JacobiDet[jCell];
                        double basisScale_jcell = m_GridDat.ChefBasis.Scaling[jCell];

                        // for the Matrix-Part of the operator, scaling of trial and test function and integral transform metric should cancel out.
                        Debug.Assert(Math.Abs(detJac_jCell * basisScale_jcell * basisScale_jcell - 1.0) <= 1.0e-10);
                    }
#endif
                    if(bAffineRequired) {
                        trafo = TempBuffer.GetTempMultidimensionalarray(out iTrafoBuf, Length);

                        MultidimensionalArray scale = m_GridDat.ChefBasis.Scaling;
                        for(int i = 0; i < Length; i++) {
                            int jCell = i + i0;
                            trafo[i] = 1.0 / scale[jCell];
                        }

                    } else {
                        trafo = null; // for linear parts, no transformation is required, since 'detJac_jCell * basisScale_jcell * basisScale_jcell == 1'
                    }

                } else {
                    if(maxDeg >= 0) {
                        trafo = m_GridDat.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(i0, Length, maxDeg);
                    }
                }



                int I0Row = 0;
                for(int gamma = 0; gamma < GAMMA; gamma++) { // loop over codomain variables ...
                    int I0Col = 0;

                    if(bLinearRequired) {
                        for(int delta = 0; delta < DELTA; delta++) { // loop over domain variables ...
                            GetQRbufferLinear(Length, QuadResult, Affine, MR, NR, ref iQlBuf, ref Ql, I0Row, gamma, I0Col, delta, out _R, out _Q);

                            double cF = 0.0;

                            if(m_UxVSumBuffer[gamma, delta] != null) {
                                GetZbuffer(Length, NoOfNodes, -1, NZ3, ref iBufZ3, ref Z3, out _Z, this.m_ColL[delta]);

                                _Z.Multiply(1.0, m_Trf_UxVSumBuffer[gamma, delta], U[delta], 0.0, "ikn", "ik", "kn");
                                _R.Multiply(1.0, V_XquadWgt[gamma], _Z, cF, "imn", "km", "ikn");
                                cF = 1.0;
                            }
                            if(m_GradUxGradVSumBuffer[gamma, delta] != null) {
                                GetZbuffer(Length, NoOfNodes, D, NZ4, ref iBufZ4, ref Z4, out _Z, this.m_ColL[delta]);

                                _Z.Multiply(1.0, m_Trf_GradUxGradVSumBuffer[gamma, delta], GradU[delta], 0.0, "ikna", "ikab", "knb");
                                _R.Multiply(1.0, GradV_XquadWgt[gamma], _Z, cF, "imn", "kma", "ikna");
                                cF = 1.0;
                            }
                            if(m_GradUxVSumBuffer[gamma, delta] != null) {
                                GetZbuffer(Length, NoOfNodes, -1, NZ3, ref iBufZ3, ref Z3, out _Z, this.m_ColL[delta]);

                                _Z.Multiply(1.0, m_Trf_GradUxVSumBuffer[gamma, delta], GradU[delta], 0.0, "ikn", "ikb", "knb");
                                _R.Multiply(1.0, V_XquadWgt[gamma], _Z, cF, "imn", "km", "ikn");
                                cF = 1.0;
                            }
                            if(m_UxGradVSumBuffer[gamma, delta] != null) {
                                GetZbuffer(Length, NoOfNodes, -1, NZ3, ref iBufZ3, ref Z3, out _Z, this.m_RowL[gamma]);

                                _Z.Multiply(1.0, m_Trf_UxGradVSumBuffer[gamma, delta], GradV_XquadWgt[gamma], 0.0, "ikm", "ika", "kma");
                                _R.Multiply(1.0, U[delta], _Z, cF, "imn", "kn", "ikm");
                                cF = 1.0;
                            }

                            I0Col += m_ColL[delta];
                            
                            if(Affine) {
                                // nop: no more trafo required
                                
                            } else {
                                _Q.Multiply(1.0, _R, ExtractTrafo(trafo, this.m_ColL[delta]), 0.0, "jbn", "jba", "jan"); // multiply from right with trafo
                                _R.Multiply(1.0, ExtractTrafo(trafo, this.m_RowL[gamma]), _Q, 0.0, "jmn", "jbm", "jbn"); // multiply from left with trafo^T
                            }
                        }

                    }
                    if(bAffineRequired) {
                        GetRQbufferAffine(Length, QuadResult, Affine, MR, ref iQaBuf, ref Qa, I0Row, gamma, I0Col, out _R, out _Q);
                        double cF = 0.0;

                        if(m_VSumBuffer[gamma, 0] != null) {
                            _R.Multiply(1.0, m_Trf_VSumBuffer[gamma, 0], V_XquadWgt[gamma], cF, "im", "ik", "km");
                            cF = 1.0;
                        }
                        if(m_GradVSumBuffer[gamma, 0] != null) {
                            _R.Multiply(1.0, m_Trf_GradVSumBuffer[gamma, 0], GradV_XquadWgt[gamma], cF, "im", "ika", "kma");
                        }

                        I0Col += 1;

                        if(Affine) {
                            _R.Multiply(1.0, trafo, _R, 0.0, "im", "i", "im");
                        } else {
                            _Q.Multiply(1.0, ExtractTrafo(trafo, this.m_RowL[gamma]), _R, 0.0, "im", "imn", "in"); 
                        }
                    }

                    I0Row += m_RowL[gamma];
                }

                // free the temp buffers
                if(Z3 != null)
                    TempBuffer.FreeTempBuffer(iBufZ3);
                if(Z4 != null)
                    TempBuffer.FreeTempBuffer(iBufZ4);
                if(iTrafoBuf > 0)
                    TempBuffer.FreeTempBuffer(iTrafoBuf);
                if(iQaBuf > 0)
                    TempBuffer.FreeTempBuffer(iQaBuf);
                if(iQlBuf > 0)
                    TempBuffer.FreeTempBuffer(iQlBuf);
            }
            this.Loops.Stop();
            #endregion
        }

        private MultidimensionalArray ExtractTrafo(MultidimensionalArray trafo, int N) {
            if(trafo.GetLength(1) == N)
                return trafo;
            else
                return trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { trafo.GetLength(0) - 1, N - 1, N - 1 });
        }

        private void GetQRbufferLinear(int Length, MultidimensionalArray QuadResult, bool Affine, int MR, int NR, ref int iQlBuf, ref MultidimensionalArray Ql, int I0Row, int gamma, int I0Col, int delta, out MultidimensionalArray _R, out MultidimensionalArray _Q) {
            _R = QuadResult.ExtractSubArrayShallow(new int[] { 0, I0Row, I0Col }, new int[] { Length - 1, I0Row + m_RowL[gamma] - 1, I0Col + m_ColL[delta] - 1 });
            if(Affine) {
                _Q = null;
            } else {
                if(Ql == null)
                    Ql = TempBuffer.GetTempMultidimensionalarray(out iQlBuf, Length, MR, NR);

                if(m_RowL[gamma] != MR || m_ColL[delta] != NR)
                    _Q = Ql.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, m_RowL[gamma] - 1, m_ColL[delta] - 1 });
                else
                    _Q = Ql;

            }
        }

        private void GetRQbufferAffine(int Length, MultidimensionalArray QuadResult, bool Affine, int MR, ref int iQaBuf, ref MultidimensionalArray Qa, int I0Row, int gamma, int I0Col, out MultidimensionalArray _R, out MultidimensionalArray _Q) {
            if(Affine) {
                _R = QuadResult.ExtractSubArrayShallow(new int[] { 0, I0Row, I0Col }, new int[] { Length - 1, I0Row + m_RowL[gamma] - 1, I0Col - 1 });
                _Q = null;
            } else {
                if(Qa == null)
                    Qa = TempBuffer.GetTempMultidimensionalarray(out iQaBuf, Length, MR);

                if(m_RowL[gamma] != MR)
                    _R = Qa.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Length - 1, m_RowL[gamma] - 1 });
                else
                    _R = Qa;

                _Q = QuadResult.ExtractSubArrayShallow(new int[] { 0, I0Row, I0Col }, new int[] { Length - 1, I0Row + m_RowL[gamma] - 1, I0Col - 1 });
            }
        }


        static MultidimensionalArray.MultiplyProgram mp_ikae_iad_ikde = MultidimensionalArray.MultiplyProgram.Compile("ikae", "iad", "ikde");
        static MultidimensionalArray.MultiplyProgram mp_ikae_ikad_ikde = MultidimensionalArray.MultiplyProgram.Compile("ikae", "ikad", "ikde");
        static MultidimensionalArray.MultiplyProgram mp_ikab_ikae_ibe = MultidimensionalArray.MultiplyProgram.Compile("ikab", "ikae", "ibe");
        static MultidimensionalArray.MultiplyProgram mp_ikab_ikae_ikbe = MultidimensionalArray.MultiplyProgram.Compile("ikab", "ikae", "ikbe");
        static MultidimensionalArray.MultiplyProgram mp_knde_k_knde = MultidimensionalArray.MultiplyProgram.Compile("knde", "k", "knde");
        static MultidimensionalArray.MultiplyProgram mp_ika_iad_ikd = MultidimensionalArray.MultiplyProgram.Compile("ika", "iad", "ikd");
        static MultidimensionalArray.MultiplyProgram mp_ika_ikad_ikd = MultidimensionalArray.MultiplyProgram.Compile("ika", "ikad", "ikd");
        static MultidimensionalArray.MultiplyProgram mp_ikb_ikd_ibd = MultidimensionalArray.MultiplyProgram.Compile("ikb", "ikd", "ibd");
        static MultidimensionalArray.MultiplyProgram mp_ikb_ikd_ikbd = MultidimensionalArray.MultiplyProgram.Compile("ikb", "ikd", "ikbd");
        static MultidimensionalArray.MultiplyProgram mp_kn_k_kn = MultidimensionalArray.MultiplyProgram.Compile("kn", "k", "kn");
        static MultidimensionalArray.MultiplyProgram mp_knd_k_knd = MultidimensionalArray.MultiplyProgram.Compile("knd", "k", "knd");
        
        private void GetZbuffer(int Length, int NoOfNodes, int D, int NZ, ref int iBufZ, ref MultidimensionalArray Z, out MultidimensionalArray _Z, int N) {
            Debug.Assert(N >= 0);
            Debug.Assert(NZ >= 0);
            Debug.Assert(N <= NZ);

            if(Z == null) {
                if(D < 0) {
                    Z = TempBuffer.GetTempMultidimensionalarray(out iBufZ, Length, NoOfNodes, NZ);
                } else {
                    Z = TempBuffer.GetTempMultidimensionalarray(out iBufZ, Length, NoOfNodes, NZ, D);
                }
            } 
            
            if(Z.GetLength(2) != N)
                if(D < 0) {
                    Debug.Assert(Z.Dimension == 3);
                    _Z = Z.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, N - 1 });
                } else {
                    Debug.Assert(Z.Dimension == 4);
                    _Z = Z.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, N - 1, D - 1 });
                } 
            else
                _Z = Z;
        }

        private void EvalNSumForm<EE>(ref VolumFormParams vfp, EquationComponentArgMapping<EE>[] Comps, MultidimensionalArray[][] CompBuffer, MultidimensionalArray[,] SumBuffer, Action<EE, MultidimensionalArray> evalForm, bool affine, Stopwatch[][] watches)
                where EE : IEquationComponent //
        {
            int DELTA = affine ? 1 : this.DELTA;
            for(int gamma = GAMMA - 1; gamma >= 0; gamma--) {
                var ecq = Comps[gamma];

                for(int delta = 0; delta < DELTA; delta++) {
                    if(SumBuffer[gamma, delta] != null)
                        SumBuffer[gamma, delta].Clear();
                }

                Stopwatch[] watches_gamma = watches[gamma];

                for (int i = 0; i < ecq.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components

                    var comp = ecq.m_AllComponentsOfMyType[i];
                    int NoOfArgs = ecq.NoOfArguments[i];
                    Debug.Assert(NoOfArgs == comp.ArgumentOrdering.Count);
                    int NoOfParams = ecq.NoOfParameters[i];
                    Debug.Assert(NoOfParams == ((comp.ParameterOrdering != null) ? comp.ParameterOrdering.Count : 0));

                    // map parameters
                    vfp.ParameterVars = new MultidimensionalArray[NoOfParams];
                    for(int c = 0; c < NoOfParams; c++) {
                        int targ = ecq.AllToSub[i, c + NoOfArgs] - this.DELTA;
                        Debug.Assert(targ >= 0);
                        vfp.ParameterVars[c] = m_ParameterFieldsValues[targ];
                    }

                    // evaluate component
                    if(affine) {
                        watches_gamma[i].Start();
                        evalForm(comp, SumBuffer[gamma, 0]);
                        watches_gamma[i].Stop();
                    } else {
                        var CompBuffer_gamma_i = CompBuffer[gamma][i];
                        CompBuffer_gamma_i.Clear();
                        watches_gamma[i].Start();
                        evalForm(comp, CompBuffer_gamma_i);
                        watches_gamma[i].Stop();
                    }

                    // sum up
                    if(affine) {


                    } else {
                        var CompBuffer_gamma_i = CompBuffer[gamma][i];
                        for(int c = 0; c < NoOfArgs; c++) {
                            int targ = ecq.AllToSub[i, c];
                            int[] Sel = new int[CompBuffer_gamma_i.Dimension];
                            Sel.SetAll(-1);
                            Sel[2] = c;
                            SumBuffer[gamma, targ].Acc(1.0, CompBuffer[gamma][i].ExtractSubArrayShallow(Sel));
                        }
                    }
                }
            }
        }

        protected void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
            bool bLinearRequired = LinearRequired;
            bool bAffineRequired = AffineRequired;
            int M = m_RowMap.NoOfCoordinatesPerCell;
            int N = m_ColMap.NoOfCoordinatesPerCell;
            int offset = bLinearRequired ? N : 0;
            int NoUpdate = m_GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int[] geom2log = m_GridDat.iGeomCells.GeomCell2LogicalCell;

            for(int i = 0; i < Length; i++) {
                int jCell;
                if(geom2log != null)
                    jCell = geom2log[i + i0];
                else
                    jCell = i + i0;

                if(bLinearRequired) {
                    var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, 0, 0 }, new int[] { i - 1, M - 1, N - 1 });

                    int m0 = (int)this.m_RowMap.GlobalUniqueCoordinateIndex(0, jCell, 0);
                    int n0 = (int)this.m_ColMap.GlobalUniqueCoordinateIndex(0, jCell, 0);

                    //for (int m = 0; m < M; m++) {
                    //    for (int n = 0; n < N; n++) {
                    //        //m_Matrix[m0 + m, n0 + n] += BlockRes[m, n];
                    //    }
                    //}
                    m_Matrix.AccBlock(m0, n0, 1.0, BlockRes);
                }
                
                if(bAffineRequired) {
                    var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, 0, offset }, new int[] { i - 1, M - 1, offset - 1 });


                    int m0 = (int)this.m_RowMap.LocalUniqueCoordinateIndex(0, jCell, 0);

                    for(int m = 0; m < M; m++)
                        m_Vector[m0 + m] += BlockRes[m];
                }
                 
            }
        }

    }
}
