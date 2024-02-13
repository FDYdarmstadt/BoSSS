﻿/* =======================================================================
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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP.LinSolvers;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Quadrature;
using ilPSP.LinSolvers.PARDISO;
using ilPSP;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using BoSSS.Foundation.Grid.Classic;
using MPI.Wrappers;

namespace BoSSS.Solution.LevelSetTools.Smoothing {
    /// <summary>
    /// Reduces Jumps in a Field by penalizing them as in the SIP Method.
    /// </summary>
    public class JumpPenalization {


        double m_penalty;

        public JumpPenalization(jumpPenalizationTerms terms, double penalty) {
            penaltyTerms = terms;
            m_penalty = penalty;
        }


        jumpPenalizationTerms penaltyTerms;


        public enum jumpPenalizationTerms {

            /// <summary>
            /// none
            /// </summary>
            None,

            /// <summary>
            /// activates only <see cref="JumpForm"/>
            /// </summary>
            Jump,

            /// <summary>
            /// activates only <see cref="GradientJumpForm"/>
            /// </summary>
            GradJump,

            /// <summary>
            /// activates only <see cref="GradientJumpForm2"/>
            /// </summary>
            GradJump2,

            /// <summary>
            /// activates <see cref="JumpForm"/> and <see cref="GradientJumpForm"/>
            /// </summary>
            JumpGradJump,

            /// <summary>
            /// activates <see cref="JumpForm"/> and <see cref="GradientJumpForm2"/>
            /// </summary>
            JumpGradJump2

        }

        class JumpForm : IEdgeForm {

            public double m_pnltyBase;


            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return TermActivationFlags.None;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.UxV;
                }
            }

            MultidimensionalArray h_min;

            public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
                double R = 0;

                if(h_min == null) {
                    h_min = ((GridData)(inp.GridDat)).Cells.h_min;
                }

                int jCell1 = inp.jCellIn;
                int jCell2 = inp.jCellOut;
                double h1 = h_min[jCell1];
                double h2 = h_min[jCell2];
                double penalty = m_pnltyBase / Math.Min(h1, h2);  // 100.0 / Math.Min(h1, h2);

                R -= (_uA[0] - _uB[0]) * (_vA - _vB) * penalty;

                return R;
            }

            public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
                return 0.0;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }


        class GradientJumpForm : IEdgeForm {

            public double m_pnltyBase;

            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return TermActivationFlags.None;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
                }
            }


            public bool BTerm = false;

            public bool ATerm = false;

            MultidimensionalArray h_min;

            public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
                double R = 0;
                int D = inp.D;

                for(int d = 0; d < D; d++) {
                    if(ATerm)
                        R -= (_Grad_uA[0, d] - _Grad_uB[0, d]) * inp.Normal[d] * (_vA - _vB);
                    if(BTerm)
                        R -= (_uA[0] - _uB[0]) * inp.Normal[d] * (_Grad_vA[d] - _Grad_vB[d]);
                }

                if (h_min == null) {
                    h_min = ((GridData)(inp.GridDat)).Cells.h_min;
                }

                int jCell1 = inp.jCellIn;
                int jCell2 = inp.jCellOut;
                double h1 = h_min[jCell1];
                double h2 = h_min[jCell2];
                double penalty = m_pnltyBase / Math.Min(h1, h2);

                R *= penalty;

                return R;
            }

            public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
                return 0.0;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }


        public class GradientJumpForm2 : IEdgeForm {

            public double m_pnltyBase;

            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return TermActivationFlags.None;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV;
                }
            }

            MultidimensionalArray h_min;

            public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
                double R = 0;
                int D = inp.D;

                for(int d = 0; d < D; d++) {
                    R -= (_Grad_uA[0, d] - _Grad_uB[0, d]) * (_Grad_vA[d] - _Grad_vB[d]);
                }

                if (h_min == null) {
                    h_min = ((GridData)(inp.GridDat)).Cells.h_min;
                }

                int jCell1 = inp.jCellIn;
                int jCell2 = inp.jCellOut;
                double h1 = h_min[jCell1];
                double h2 = h_min[jCell2];
                double penalty = m_pnltyBase / Math.Min(h1, h2);

                R *= penalty;

                return R;
            }

            public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
                return 0.0;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }


            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }


        public void Evaluate(EdgeMask em, SinglePhaseField inp_LevSet, SinglePhaseField outp_Result) {



            MsrMatrix Matrix;

            Matrix = PenaltyMatrix(em, inp_LevSet.Basis, outp_Result.Basis);


            Matrix.SpMV(1.0, inp_LevSet.CoordinateVector, 0.0, outp_Result.CoordinateVector);

        }

        private MsrMatrix PenaltyMatrix(EdgeMask em, Basis LevSetBasis, Basis JumpBasis) {

            var OpA = new DifferentialOperator(1, 0, 1, QuadOrderFunc.Linear(), "Phi", "c1");
            switch (penaltyTerms) {
                case jumpPenalizationTerms.Jump: {
                        OpA.EquationComponents["c1"].Add(new JumpForm() { m_pnltyBase = m_penalty});
                        break;
                    }
                case jumpPenalizationTerms.GradJump: {
                        OpA.EquationComponents["c1"].Add(new GradientJumpForm() { m_pnltyBase = m_penalty, ATerm = true, BTerm = true });
                        break;
                    }
                case jumpPenalizationTerms.GradJump2: {
                        OpA.EquationComponents["c1"].Add(new GradientJumpForm2() { m_pnltyBase = m_penalty });
                        break;
                    }
                case jumpPenalizationTerms.JumpGradJump: {
                        OpA.EquationComponents["c1"].Add(new JumpForm() { m_pnltyBase = m_penalty });
                        OpA.EquationComponents["c1"].Add(new GradientJumpForm() { m_pnltyBase = m_penalty, ATerm = true, BTerm = true });
                        break;
                    }
                case jumpPenalizationTerms.JumpGradJump2: {
                        OpA.EquationComponents["c1"].Add(new JumpForm() { m_pnltyBase = m_penalty });
                        OpA.EquationComponents["c1"].Add(new GradientJumpForm2() { m_pnltyBase = m_penalty });
                        break;
                    }
                default:
                    throw new ArgumentException();
            }
            OpA.Commit();
            
            //var OpB = new SpatialOperator(1, 0, 1, "Phi", "c1");
            //Op.EquationComponents["c1"].Add(new JumpForm());
            //OpB.EquationComponents["c1"].Add(new GradientJumpForm() { BTerm = true });
            //Op.EquationComponents["c1"].Add(new GradientJumpForm2());
            //OpB.Commit();

            var inp_LevSet_Mapping = new UnsetteledCoordinateMapping(LevSetBasis);
            var outp_Result_Mapping = new UnsetteledCoordinateMapping(JumpBasis);

            MsrMatrix MatrixA;
            MatrixA = new MsrMatrix(outp_Result_Mapping, inp_LevSet_Mapping);
            double[] AffineA = new double[inp_LevSet_Mapping.LocalLength];
            OpA.ComputeMatrixEx(inp_LevSet_Mapping, null, outp_Result_Mapping,
                    MatrixA, AffineA, OnlyAffine: false,
                    edgeQuadScheme: new EdgeQuadratureScheme(true, em),
                    volQuadScheme: new CellQuadratureScheme(true, CellMask.GetEmptyMask(em.GridData)));
            MatrixA.CheckForNanOrInfM(typeof(JumpPenalization) + ": ");

            //MsrMatrix MatrixB;
            //MatrixB = new MsrMatrix(outp_Result_Mapping, inp_LevSet_Mapping);
            //double[] AffineB = new double[inp_LevSet_Mapping.LocalLength];
            //OpB.ComputeMatrixEx(inp_LevSet_Mapping, null, outp_Result_Mapping,
            //    MatrixB, AffineB, OnlyAffine: false
            //    );//,
            //    //edgeQrCtx: new EdgeQuadratureScheme(true, em),
            //    //volQrCtx: new CellQuadratureScheme(true, CellMask.GetEmptyMask(em.GridData)));

            //Debug.Assert(AffineB.L2Norm() == 0);

            //var Err = MatrixA.Transpose();
            //Err.Acc(-1.0, MatrixB);

            //double errnorm = Err.InfNorm();
            //Debug.Assert(errnorm < 1.0e-10);
            ////Console.WriteLine("Errnorm:" + errnorm);


            return MatrixA;
        }

        public void Evaluate2(SubGrid S, SinglePhaseField inp_LevSet, SinglePhaseField outp_Result) {

            var Op = new DifferentialOperator(1, 0, 1, QuadOrderFunc.Linear(), "Phi", "c1");
            Op.EquationComponents["c1"].Add(new JumpForm());
            //Op.EquationComponents["c1"].Add(new GradientJumpForm() { BTerm = true });
            Op.EquationComponents["c1"].Add(new GradientJumpForm2());
            Op.Commit();
            

            var inp_LevSet_Mapping = inp_LevSet.Mapping;
            var outp_Result_Mapping = outp_Result.Mapping;


            Op.Evaluate(1.0, 1.0,
                inp_LevSet_Mapping, new DGField[0], outp_Result_Mapping);//,
                //qInsEdge: new EdgeQuadratureScheme(true, S.InnerEdgesMask),
                //qInsVol: new CellQuadratureScheme(true, CellMask.GetEmptyMask(S._GridData)),
                //bndMode: SubGridBoundaryModes.InnerEdge,
                //sgrd: S);                
        }


        public void ImplicitEuler(double dt, SubGrid S, SinglePhaseField inout_Levset) {
            var VolMsk = S.VolumeMask;
            var EdgMsk = S.InnerEdgesMask;
            UnsetteledCoordinateMapping Map = inout_Levset.Mapping;

            if (dt <= 0.0)
                throw new ArgumentOutOfRangeException("Timestep size must be greater than 0.");

            MsrMatrix Pmtx = PenaltyMatrix(EdgMsk, inout_Levset.Basis, inout_Levset.Basis);
            Pmtx.Scale(-1.0);

            long[] SubVecIdx = Map.GetSubvectorIndices(S, true, new int[] { 0 });
            int L = SubVecIdx.Length;

            MsrMatrix SubMtx = new MsrMatrix(L, L);
            Pmtx.AccSubMatrixTo(1.0, SubMtx, SubVecIdx, default(long[]), SubVecIdx, default(long[]));

            SubMtx.AccEyeSp(1.0 / dt);

            double[] RHS = new double[L];
            double[] SOL = new double[L];

            int[] SubVecIdxLoc = new int[SubVecIdx.Length];
            for (int idx = 0; idx < SubVecIdxLoc.Length; idx++) {
                SubVecIdxLoc[idx] = Map.TransformIndexToLocal(SubVecIdx[idx]);
            }
            RHS.AccV(1.0 / dt, inout_Levset.CoordinateVector, default(int[]), SubVecIdxLoc);
            
            using (var solver = new PARDISOSolver()) {
                solver.DefineMatrix(SubMtx);
                solver.Solve(SOL, RHS);
            }

            inout_Levset.CoordinateVector.ClearEntries(SubVecIdxLoc);
            inout_Levset.CoordinateVector.AccV(1.0, SOL, SubVecIdxLoc, default(int[]));
        }


    }
}

