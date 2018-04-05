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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {
    
    /// <summary>
    /// Solver for the Reinit-problem on one cell, based on the elliptic approach of Basting&Kuzmin .
    /// </summary>
    class LocalSolver_Elliptic {

        GridData GridDat;
        UnsetteledCoordinateMapping LevelSetMapping;
        Basis LevelSetBasis;
                

        /// <summary>
        /// ctor
        /// </summary>
        public LocalSolver_Elliptic(Basis __LevelSetBasis) {
            this.GridDat = (GridData)(__LevelSetBasis.GridDat);
            this.LevelSetMapping = new UnsetteledCoordinateMapping(__LevelSetBasis);
            this.LevelSetBasis = __LevelSetBasis;

            this.m_LaplaceMatrix = new MsrMatrix(this.LevelSetMapping, this.LevelSetMapping);
            this.m_LaplaceAffine = new double[this.LevelSetMapping.LocalLength];
        }
        

        /// <summary>
        /// Laplace operator for the elliptic reinitialization procedure.
        /// </summary>
        class EllipticReinitForm : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {

            /// <summary>
            /// ctor.
            /// </summary>
            public EllipticReinitForm(BitArray _AcceptedBitmask, int __jCell, double __penaltyBase, MultidimensionalArray __cj) {
                this.AcceptedBitmask = _AcceptedBitmask;
                this.penaltyBase = __penaltyBase;
                this.jCell = __jCell;
                this.cj = __cj;
                if(this.AcceptedBitmask[this.jCell])
                    throw new ArgumentException("Cannot work on accepted cells.");
            }

            double penaltyBase;
            BitArray AcceptedBitmask;
            int jCell;


            /// <summary>
            /// Turns the right-hand-side, i.e. the term $ -\divergence{ \frac{ \nabla \phi_0 }{ | \nabla \phi_0 | }$ on or off.
            /// </summary>
            public double RhsSwitch = 1.0;

            /// <summary>
            /// Turns the left-hand-side, i.e. the term $ \divergence{ \nabla \phi_1 | }$ on or off.
            /// </summary>
            public double LhsSwitch = 1.0;



            public TermActivationFlags VolTerms {
                get {
                    return (TermActivationFlags.GradUxGradV | TermActivationFlags.GradV);
                }
            }

            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return InnerEdgeTerms;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.V | TermActivationFlags.GradV;
                }
            }

            static double d1Clamp(double s) {
                Debug.Assert(s >= 0);
                double Pot;
                if(s < 1.0e-14) {
                    Pot = -1.0e14;
                } else {
                    Pot = 1 - 1 / s;
                }
                Pot = Math.Max(-5.0, Pot);
                return Pot;
            }

            static double d3(double s) {
                Debug.Assert(s >= 0);
                if(s <= 1) {
                    return (2.0 * s * s - 3.0 * s + 1.0);
                } else {
                    return (1.0 - 1.0 / s);
                }
            }

            private double OneOverGradPhi0(double[] Params) {

                double absGradPhi0 = 0;
                int D = Params.Length - 1;
                for(int d = 0; d < D; d++) {
                    absGradPhi0 += Params[1 + d].Pow2();
                }
                absGradPhi0 = Math.Sqrt(absGradPhi0);


                //double R = 1 - d1Clamp(absGradPhi0);
                double R = 1 - d3(absGradPhi0);

                return R; // includes the switch for the right-hand-side
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] Phi, double[,] GradPhi, double V, double[] GradV) {
                int D = cpv.D;
                Debug.Assert(GradPhi.GetLength(0) == 1);
                Debug.Assert(GradPhi.GetLength(1) == D);
                Debug.Assert(AcceptedBitmask[cpv.jCell] == false, "Cannot work on accepted cells");

                double Acc = 0;

                double ooGradPhi0 = OneOverGradPhi0(cpv.Parameters);

                for(int d = 0; d < D; d++) {
                    Acc -= GradPhi[0, d] * GradV[d] * LhsSwitch;
                    Acc -= ooGradPhi0 * cpv.Parameters[1 + d] * GradV[d] * RhsSwitch;
                }

                Debug.Assert(!(double.IsNaN(Acc) || double.IsInfinity(Acc)));
                return Acc;
            }



            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return new string[] { "Phi0", "dPhi0_dx0", "dPhi0_dx1" };
                }
            }

            MultidimensionalArray cj;

            public double InnerEdgeForm(ref CommonParams inp, double[] PhiIn, double[] PhiOt, double[,] Grad_PhiIn, double[,] Grad_PhiOt, double _vIn, double _vOt, double[] _Grad_vIn, double[] _Grad_vOt) {
                Debug.Assert((this.AcceptedBitmask[inp.jCellIn] && this.AcceptedBitmask[inp.jCellOut]) == false, "In and Out - cell accepted: cannot be!");
                int D = inp.D;
                double penalty = cj[this.jCell] * this.penaltyBase;

                double Acc = 0;

                if(this.AcceptedBitmask[inp.jCellIn] == true && this.AcceptedBitmask[inp.jCellOut] == false) {
                    // IN-cell is accepted // OUT-cell should be computed
                    // => the boundary value is given as IN-parameter
                    // => flux penalizes the OUT-Cell
                    Debug.Assert(jCell == inp.jCellOut);

                    Acc = -PenalizedEdge(PhiOt, Grad_PhiOt, _vOt, _Grad_vOt, inp.Normale, inp.Parameters_OUT, inp.Parameters_IN, -penalty);


                } else if(this.AcceptedBitmask[inp.jCellIn] == false && this.AcceptedBitmask[inp.jCellOut] == true) {
                    // ... vice-versa
                    Debug.Assert(jCell == inp.jCellIn);

                    Acc = +PenalizedEdge(PhiIn, Grad_PhiIn, _vIn, _Grad_vIn, inp.Normale, inp.Parameters_IN, inp.Parameters_OUT, penalty);

                } else {
                    // free edge

                    if(jCell == inp.jCellIn) {
                        Acc = +FreeEdge(Grad_PhiIn, _vIn, inp.Normale, inp.Parameters_IN);
                    } else if(jCell == inp.jCellOut) {
                        Acc = -FreeEdge(Grad_PhiOt, _vOt, inp.Normale, inp.Parameters_OUT);
                    } else {
                        Debug.Assert(false);
                    }
                }

                return Acc;
            }



            private double FreeEdge(double[,] Grad_Phi, double _v, double[] Normale, double[] Parameters_IN) {
                double Acc = 0;

                /*
                int D = Normale.Length;
                double ooGradPhi0 = OneOverGradPhi0(Parameters_IN);

                for(int d = 0; d < D; d++) {
                    Acc -= Parameters_IN[1 + d] * Normale[d] * _v * ooGradPhi0;
                }

                for(int d = 0; d < D; d++) {
                    Acc += Grad_Phi[0, d] * Normale[d] * _v * LhsSwitch*0;
                }
                Debug.Assert(!(double.IsNaN(Acc) || double.IsInfinity(Acc)));
                */

                return Acc;
            }

            private double PenalizedEdge(double[] Phi, double[,] Grad_Phi, double _v, double[] _Grad_v,
                double[] Normale, double[] Parameters_IN, double[] Parameters_OUT, double penalty) {
                double Acc = 0;
                int D = Normale.Length;

                double ooGradPhi0 = OneOverGradPhi0(Parameters_IN);

                for(int d = 0; d < D; d++) {
                    Acc += Parameters_IN[1 + d] * Normale[d] * _v * ooGradPhi0 * RhsSwitch;
                }

                for(int d = 0; d < D; d++) {
                    // central difference 
                    Acc += 0.5 * (Grad_Phi[0, d] + Parameters_OUT[1 + d]) * Normale[d] * _v * LhsSwitch;
                    Acc += 0.5 * _Grad_v[d] * Normale[d] * (Phi[0] - Parameters_OUT[0]) * LhsSwitch; // symmetry term

                    // inner values
                    //Acc += Grad_Phi[0, d] * Normale[d] * _v * LhsSwitch;
                    //Acc += _Grad_v[d] * Normale[d] * (Phi[0] - Parameters_OUT[0]) * LhsSwitch; // symmetry term
                }

                Acc -= penalty * (Phi[0] - Parameters_OUT[0]) * _v * LhsSwitch;
                Debug.Assert(!(double.IsNaN(Acc) || double.IsInfinity(Acc)));

                return Acc;
            }


            public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] PhiIn, double[,] Grad_PhiIn, double _vIn, double[] _Grad_vIn) {
                // free edge

                double Acc = FreeEdge(Grad_PhiIn, _vIn, inp.Normale, inp.Parameters_IN);
                return Acc;
            }
        }

        MsrMatrix m_LaplaceMatrix;
        double[] m_LaplaceAffine;

        Stopwatch Stpw_Rhs = new Stopwatch();
        Stopwatch Stpw_Mtx = new Stopwatch();
        Stopwatch Stpw_tot = new Stopwatch();

        public void PrintInstrumentation() {
            double tot = Stpw_tot.Elapsed.TotalSeconds;
            double rhs = Stpw_Rhs.Elapsed.TotalSeconds;
            double mtx = Stpw_Mtx.Elapsed.TotalSeconds;
           

            Console.WriteLine(" Reinit total runtime: " + tot);
            Console.WriteLine("    right-hand-side:   {0:0.##E-00} {1:0.##E-00}%", rhs, 100 * rhs / tot);
            Console.WriteLine("    matrix:            {0:0.##E-00} {1:0.##E-00}%", mtx, 100 * mtx / tot);

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="jCell">
        /// Local index of cell which should be recalculated.
        /// </param>
        /// <param name="AcceptedMask">
        /// Bitmask which marks accepted cells - if any neighbor of
        /// <paramref name="jCell"/> is _accepted_, this defines a Dirichlet
        /// boundary; otherwise, the respective cell face is a free boundary.
        /// </param>
        /// <param name="Phi">
        /// Input and output:
        /// - input for _accepted_ cells, i.e. Dirichlet boundary values
        /// - input and output for <paramref name="jCell"/>: an initial value for the iterative procedure, resp. on exit the result of the iteration.
        /// </param>
        /// <param name="gradPhi">
        /// Auxillary variable to store the gradient of the level-set-field.
        /// </param>
        /// <returns></returns>
        public bool LocalSolve(int jCell, BitArray AcceptedMask, SinglePhaseField Phi, VectorField<SinglePhaseField> gradPhi) {
            Stpw_tot.Start();

            int N = this.LevelSetBasis.GetLength(jCell);
            int i0G = this.LevelSetMapping.GlobalUniqueCoordinateIndex(0, jCell, 0);
            int i0L = this.LevelSetMapping.LocalUniqueCoordinateIndex(0, jCell, 0);
            double penaltyBase = ((double)(this.LevelSetBasis.Degree + 2)).Pow2();
            double CellVolume = this.GridDat.Cells.GetCellVolume(jCell);
            int p = Phi.Basis.Degree;

            // subgrid on which we are working, consisting only of one cell
            SubGrid jCellGrid = new SubGrid(new CellMask(this.GridDat, Chunk.GetSingleElementChunk(jCell)));
            var VolScheme = new CellQuadratureScheme(domain: jCellGrid.VolumeMask);
            var EdgScheme = new EdgeQuadratureScheme(domain: jCellGrid.AllEdgesMask);
            var VolRule = VolScheme.SaveCompile(GridDat, 3 * p);
            var EdgRule = EdgScheme.SaveCompile(GridDat, 3 * p);

            // parameter list for operator
            DGField[] Params = new DGField[] { Phi };
            gradPhi.ForEach(f => f.AddToArray(ref Params));

            // build operator
            var comp = new EllipticReinitForm(AcceptedMask, jCell, penaltyBase, this.GridDat.Cells.cj);
            comp.LhsSwitch = 1.0;  // matrix is constant -- Lhs Matrix only needs to be computed once
            comp.RhsSwitch = -1.0; //                   Rhs must be updated in every iteration
            var op = comp.Operator();


            // iteration loop:
            MultidimensionalArray Mtx = MultidimensionalArray.Create(N, N);
            double[] Rhs = new double[N];
            double ChangeNorm = 0;
            int iIter;
            for(iIter = 0; iIter < 100; iIter++) {

                // update gradient
                gradPhi.Clear(jCellGrid.VolumeMask);
                gradPhi.Gradient(1.0, Phi, jCellGrid.VolumeMask);

                // assemble matrix and rhs
                {
                    // clear
                    for(int n = 0; n < N; n++) {
                        if(iIter == 0)
                            m_LaplaceMatrix.ClearRow(i0G + n);
                        this.m_LaplaceAffine[i0L + n] = 0;
                    }

                                        // compute matrix (only in iteration 0) and rhs
                    if(iIter == 0)
                        Stpw_Mtx.Start();
                    Stpw_Rhs.Start();
                    op.ComputeMatrixEx(this.LevelSetMapping, Params, this.LevelSetMapping,
                        iIter == 0 ? this.m_LaplaceMatrix : null, this.m_LaplaceAffine,
                            OnlyAffine: iIter > 0,
                            //volQrCtx: VolScheme, edgeQrCtx: EdgScheme,
                            volRule: VolRule, edgeRule: EdgRule,
                            ParameterMPIExchange:false);
                    //op.Internal_ComputeMatrixEx(this.GridDat,
                    //    this.LevelSetMapping, Params, this.LevelSetMapping,
                    //    iIter == 0 ? this.m_LaplaceMatrix : default(MsrMatrix), this.m_LaplaceAffine, iIter > 0,
                    //    0.0,
                    //    EdgRule, VolRule,
                    //    null, false);


                    if(iIter == 0)
                        Stpw_Mtx.Stop();
                    Stpw_Rhs.Stop();

                    // extract matrix for 'jCell'
                    for(int n = 0; n < N; n++) {
#if DEBUG
                        int Lr;
                        int[] row_cols = null;
                        double[] row_vals = null;
                        Lr = this.m_LaplaceMatrix.GetRow(i0G + n, ref row_cols, ref row_vals);
                        for (int lr = 0; lr < Lr; lr++) {
                            int ColIndex = row_cols[lr];
                            double Value = row_vals[lr];
                            Debug.Assert((ColIndex >= i0G && ColIndex < i0G + N) || (Value == 0.0), "Matrix is expected to be block-diagonal.");
                        }
#endif
                        if(iIter == 0) {
                            for(int m = 0; m < N; m++) {
                                Mtx[n, m] = this.m_LaplaceMatrix[i0G + n, i0G + m];
                            }
                        } else {
#if DEBUG
                            for(int m = 0; m < N; m++) {
                                Debug.Assert(Mtx[n, m] == this.m_LaplaceMatrix[i0G + n, i0G + m]);
                            }
#endif
                        }
                        Rhs[n] = -this.m_LaplaceAffine[i0L + n];
                    }

                }

                // solve
                double[] sol = new double[N];
                Mtx.Solve(sol, Rhs);
                ChangeNorm = GenericBlas.L2Dist(sol, Phi.Coordinates.GetRow(jCell));
                Phi.Coordinates.SetRow(jCell, sol);
                
                if(ChangeNorm / CellVolume < 1.0e-10)
                    break;
            }

            
            //Console.WriteLine("Final change norm: {0}, \t iter: {1}", ChangeNorm, iIter );

            if(ChangeNorm > 1.0e-6)
                Console.WriteLine("  local solver funky in cell: " + jCell + ", last iteration change norm = " + ChangeNorm);


            Stpw_tot.Stop();
            return true;
        }
        
    }
}
