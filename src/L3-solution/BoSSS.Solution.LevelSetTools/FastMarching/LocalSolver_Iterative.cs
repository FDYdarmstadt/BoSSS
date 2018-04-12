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
using BoSSS.Solution.Timestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {
    class LocalSolver_Iterative {

        GridData GridDat;
        UnsetteledCoordinateMapping LevelSetMapping;
        Basis LevelSetBasis;
        GradientModule gm;

        /// <summary>
        /// ctor
        /// </summary>
        public LocalSolver_Iterative(Basis __LevelSetBasis, GradientModule __gm) {
            this.GridDat = (GridData)(__LevelSetBasis.GridDat);
            this.LevelSetMapping = new UnsetteledCoordinateMapping(__LevelSetBasis);
            this.LevelSetBasis = __LevelSetBasis;
            this.gm = __gm;
        }
        
        /// <summary>
        /// Artificial viscosity/diffusion that is added when the local iterative solver does not converge.
        /// </summary>
        class ArtificialViscosity : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {
            public ArtificialViscosity(BitArray _AcceptedBitmask, double __penaltyBase, MultidimensionalArray __h_min, int __jCell, double __DiffusionCoeff) {
                this.AcceptedBitmask = _AcceptedBitmask;
                this.penaltyBase = __penaltyBase;
                this.jCell = __jCell;
                this.DiffusionCoeff = __DiffusionCoeff;
                this.h_min = __h_min;
            }

            double penaltyBase;
            BitArray AcceptedBitmask;
            int jCell;
            double DiffusionCoeff;

            public TermActivationFlags VolTerms {
                get {
                    return (TermActivationFlags.GradUxGradV);
                }
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] Phi, double[,] GradPhi, double V, double[] GradV) {

                int D = cpv.D;
                Debug.Assert(GradPhi.GetLength(0) == 1);
                Debug.Assert(GradPhi.GetLength(1) == D);
                Debug.Assert(AcceptedBitmask[cpv.jCell] == false, "Cannot work on accepted cells");

                double Diffusion = 0;
                for(int d = 0; d < D; d++) {
                    Diffusion -= GradPhi[0, d] * GradV[d];
                }

                return Diffusion * DiffusionCoeff;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return new string[] { "Phi0", "dPhi0_dx", "dPhi0_dy" };
                }
            }

            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return InnerEdgeTerms;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.GradV | TermActivationFlags.V;
                }
            }

            MultidimensionalArray h_min;

            public double InnerEdgeForm(ref CommonParams inp, double[] PhiIn, double[] PhiOt, double[,] GradPhiIn, double[,] GradPhiOt, double vIn, double vOt, double[] Grad_vIn, double[] Grad_vOt) {
                Debug.Assert((this.AcceptedBitmask[inp.jCellIn] && this.AcceptedBitmask[inp.jCellOut]) == false, "In and Out - cell accepted: cannot be!");
                int D = inp.D;
                double penalty = penaltyBase / h_min[this.jCell];


                if(this.AcceptedBitmask[inp.jCellIn] == true && this.AcceptedBitmask[inp.jCellOut] == false) {
                    // IN-cell is accepted 
                    // OUT-cell should be computed
                    // => the boundary value is given as IN-parameter
                    // => flux penalizes the OUT-Cell

                    Debug.Assert(inp.jCellOut == jCell);

                    double bndValue = inp.Parameters_IN[0];

                    double Diffusion = 0;
                    for(int d = 0; d < D; d++) {
                        Diffusion += GradPhiOt[0, d] * (0 - vOt) * inp.Normale[d];//            consistency
                        Diffusion += Grad_vOt[d] * (bndValue - PhiOt[0]) * inp.Normale[d]; //   symmetry
                    }
                    Diffusion -= (bndValue - PhiOt[0]) * (0 - vOt) * penalty; // penalty

                    return Diffusion * DiffusionCoeff;

                } else if(this.AcceptedBitmask[inp.jCellIn] == false && this.AcceptedBitmask[inp.jCellOut] == true) {
                    // ... vice-versa

                    Debug.Assert(inp.jCellIn == jCell);

                    double bndValue = inp.Parameters_OUT[0];

                    double Diffusion = 0;
                    for(int d = 0; d < D; d++) {
                        Diffusion += GradPhiIn[0, d] * (vIn - 0) * (inp.Normale[d]); //           consistency
                        Diffusion += Grad_vIn[d] * (PhiIn[0] - bndValue) * (inp.Normale[d]); //   symmetry
                    }
                    Diffusion -= (PhiIn[0] - bndValue) * (vIn - 0) * penalty; // penalty

                    return Diffusion * DiffusionCoeff;

                } else {


                    if(this.jCell == inp.jCellIn) {
                        double g_Neu = 0;
                        for(int d = 0; d < D; d++) {
                            g_Neu += inp.Parameters_IN[1 + d] * inp.Normale[d];
                        }

                        Clamp_g_Neu(ref g_Neu);

                        return g_Neu * vIn * DiffusionCoeff; // free edge: no penalty, inhomogeneous Neumann bndy cond (if Diffusion is used)
                    } else if(this.jCell == inp.jCellOut) {
                        double g_Neu = 0;

                        for(int d = 0; d < D; d++) {
                            g_Neu += inp.Parameters_OUT[1 + d] * inp.Normale[d];
                        }

                        Clamp_g_Neu(ref g_Neu);

                        return -g_Neu * vOt * DiffusionCoeff;
                    } else {
                        throw new ApplicationException();
                    }
                }


            }

            static void Clamp_g_Neu(ref double g_Neu) {
                if(g_Neu < -1.0)
                    g_Neu = -1.0;
                if(g_Neu > 1.0)
                    g_Neu = 1.0;
                //g_Neu = 0;
            }


            public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] PhiIn, double[,] Grad_PhiIn, double vIn, double[] Grad_vIn) {
                //return 0.0; // free edge: no penalty, homogeneous Neumann bndy cond (if Diffusion is used)
                int D = inp.D;
                double g_Neu = 0;
                for(int d = 0; d < D; d++) {
                    g_Neu += inp.Parameters_IN[1 + d] * inp.Normale[d];
                }
                g_Neu *= DiffusionCoeff;

                Clamp_g_Neu(ref g_Neu);

                return g_Neu * vIn * DiffusionCoeff;
            }
        }

        class ReinitOperator : BoSSS.Foundation.IVolumeForm {


            public TermActivationFlags VolTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.V;
                }
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                int D = cpv.D;

                double NormGradPhi = 0;
                for(int d = 0; d < D; d++) {
                    NormGradPhi += cpv.Parameters[d].Pow2();
                }
                //NormGradPhi = Math.Sqrt(NormGradPhi);

                return (NormGradPhi - 1) * V;
            }
            
            public IList<string> ArgumentOrdering {
                get {
                    return new string[0];
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return new string[] { "dPhi_dx0", "dPhi_dx1" };
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="AcceptedMask"></param>
        /// <param name="Phi"></param>
        /// <param name="gradPhi"></param>
        /// <param name="__DiffusionCoeff">Output: if artificial diffusion is turned</param>
        /// <param name="MaxAllowedPhi">Input: upper threshold for the values of <paramref name="Phi"/> in cell <see cref="jCell"/>.</param>
        /// <param name="MinAllowedPhi">Input: lower threshold for the values of <paramref name="Phi"/> in cell <see cref="jCell"/>.</param>
        /// <returns></returns>
        public bool LocalSolve_Iterative(int jCell, BitArray AcceptedMask, SinglePhaseField Phi, VectorField<SinglePhaseField> gradPhi, SinglePhaseField __DiffusionCoeff, double MaxAllowedPhi, double MinAllowedPhi) {
            //this.LocalSolve_Geometric(jCell, AcceptedMask, Phi, +1, out MinAllowedPhi, out MaxAllowedPhi);) {
            int N = this.LevelSetBasis.GetLength(jCell);
            int i0G = this.LevelSetMapping.GlobalUniqueCoordinateIndex(0, jCell, 0);
            int i0L = this.LevelSetMapping.LocalUniqueCoordinateIndex(0, jCell, 0);

            SinglePhaseField __AcceptedMask = new SinglePhaseField(new Basis(this.GridDat, 0), "accepted");
            for(int j = 0; j < AcceptedMask.Length; j++) {
                __AcceptedMask.SetMeanValue(j, AcceptedMask[j] ? 1.0 : 0.0);
            }


            // subgrid on which we are working, consisting only of one cell
            SubGrid jCellGrid = new SubGrid(new CellMask(this.GridDat, Chunk.GetSingleElementChunk(jCell)));

            // create spatial operator
            SpatialOperator.Evaluator evo;
            {
                SpatialOperator op = new SpatialOperator(1, 2, 1, QuadOrderFunc.NonLinear(2), "Phi", "dPhi_dx0", "dPhi_dx1", "cod1");
                op.EquationComponents["cod1"].Add(new ReinitOperator());
                op.Commit();

                evo = op.GetEvaluatorEx(Phi.Mapping.Fields, gradPhi.Mapping.Fields, Phi.Mapping,
                    edgeQrCtx: (new EdgeQuadratureScheme(domain: EdgeMask.GetEmptyMask(this.GridDat))),
                    volQrCtx: (new CellQuadratureScheme(domain: jCellGrid.VolumeMask)),
                    subGridBoundaryTreatment: SpatialOperator.SubGridBoundaryModes.InnerEdge);
            }

            // create artificial diffusion operator
            MultidimensionalArray DiffMtx;
            double[] DiffRhs;
            SpatialOperator dop;
            {
                double penaltyBase = this.LevelSetBasis.Degree + 2;
                penaltyBase = penaltyBase.Pow2();

                dop = (new ArtificialViscosity(AcceptedMask, penaltyBase, GridDat.Cells.h_min, jCell, -1.0)).Operator(1);

                MsrMatrix _DiffMtx = new MsrMatrix(this.LevelSetMapping, this.LevelSetMapping);
                double[] _DiffRhs = new double[this.LevelSetMapping.LocalLength];

                dop.ComputeMatrixEx(this.LevelSetMapping, new DGField[] { Phi, null, null }, this.LevelSetMapping,
                    _DiffMtx, _DiffRhs, OnlyAffine: false,
                    edgeQuadScheme: (new EdgeQuadratureScheme(domain: jCellGrid.AllEdgesMask)),
                    volQuadScheme: (new CellQuadratureScheme(domain: jCellGrid.VolumeMask)));

                // extract matrix for 'jCell'
                DiffMtx = MultidimensionalArray.Create(N, N);
                DiffRhs = new double[N];
                for(int n = 0; n < N; n++) {
#if DEBUG
                    int Lr;
                    int[] row_cols = null;
                    double[] row_vals = null;
                    Lr = _DiffMtx.GetRow(i0G + n, ref row_cols, ref row_vals);
                    for (int lr = 0; lr < Lr; lr++) {
                        int ColIndex = row_cols[lr];
                        double Value = row_vals[lr];
                        Debug.Assert((ColIndex >= i0G && ColIndex < i0G + N) || (Value == 0.0), "Matrix is expected to be block-diagonal.");
                    }
#endif
                    for(int m = 0; m < N; m++) {
                        DiffMtx[n, m] = _DiffMtx[i0G + n, i0G + m];
                    }
                    DiffRhs[n] = _DiffRhs[i0L + n];
                }

#if DEBUG
                var Test = DiffMtx.Transpose();
                Test.Acc(-1.0, DiffMtx);
                Debug.Assert(Test.InfNorm() <= 1.0e-8);
#endif

            }

            // find 'good' initial value by geometric solve AND
            // thresholds for the maximum an minimal value of Phi
            double Range = MaxAllowedPhi - MinAllowedPhi;
            MinAllowedPhi -= 0.1 * Range;
            MaxAllowedPhi += 0.1 * Range;


            // timestep for pseudo-timestepping
            double dt = 0.5 * this.GridDat.Cells.h_min[jCell] / (((double)(this.LevelSetBasis.Degree)).Pow2());

            
            DGField[] PlotFields = new DGField[] { Phi, gradPhi[0], gradPhi[1], __DiffusionCoeff, __AcceptedMask };
            //Tecplot.Tecplot.PlotFields(Params, "itt_0", "EllipicReinit", 0, 3);

            double[] PrevVal = new double[N];
            double[] NextVal = new double[N];
            Phi.Coordinates.GetRow(jCell, PrevVal);

            // pseudo-timestepping
            //if(jCell == 80)
            //    Tecplot.Tecplot.PlotFields(PlotFields, "itt_0", "EllipicReinit", 0, 3);
            //Console.Write("  Local solve cell " + jCell + " ... ");

            bool converged = false;
            double DiffusionCoeff = 0;
            int IterGrowCount = 0; // number of iterations in which the residual grew
            double LastResi = double.NaN;
            for(int iIter = 0; iIter < 1000; iIter++) {
                //Console.Write("  Local solve iteration " + iIter + " ... ");
                PerformRKstep(dt, jCell, AcceptedMask, Phi, gradPhi, evo);

                __DiffusionCoeff.SetMeanValue(jCell, DiffusionCoeff);

                if(jCell == 80) {
                    DiffusionCoeff = 0.1;
                }
                if(DiffusionCoeff > 0) {
                    //Console.WriteLine(" Diffusion on.");


                    double[] _DiffRhs = new double[this.LevelSetMapping.LocalLength];


                    dop.ComputeMatrixEx(this.LevelSetMapping, new DGField[] { Phi, gradPhi[0], gradPhi[1]
                    }, this.LevelSetMapping,
                        default(MsrMatrix), _DiffRhs, OnlyAffine: true,
                        edgeQuadScheme: (new EdgeQuadratureScheme(domain: jCellGrid.AllEdgesMask)),
                        volQuadScheme: (new CellQuadratureScheme(domain: CellMask.GetEmptyMask(this.GridDat))));

                    // extract matrix for 'jCell'
                    for(int n = 0; n < N; n++) {
                        DiffRhs[n] = _DiffRhs[i0L + n];
                    }

                    PerformArtificialDiffusion(dt * DiffusionCoeff, jCell, Phi, DiffMtx, DiffRhs);
                }
                Phi.Coordinates.GetRow(jCell, NextVal);
                double resi = Math.Sqrt(GenericBlas.L2DistPow2(NextVal, PrevVal) / GenericBlas.L2NormPow2(PrevVal));
                double[] A = NextVal;
                NextVal = PrevVal;
                PrevVal = A;
                if(iIter > 0 && resi > LastResi)
                    IterGrowCount++;
                else
                    IterGrowCount = 0;
                LastResi = resi;


                if(resi < 1.0e-10) {
                    converged = true;
                    break;
                }

                double maxPhi, minPhi;
                Phi.GetExtremalValuesInCell(out minPhi, out maxPhi, jCell);

                bool MinAlarm = minPhi < MinAllowedPhi;
                bool Maxalarm = maxPhi > MaxAllowedPhi;
                bool GrowAlarm = IterGrowCount > 4;
                bool IterAlarm = iIter >= 50;

                if(MinAlarm || Maxalarm || GrowAlarm) {
                    // Diffusion coefficient should be increased
                    if(DiffusionCoeff == 0) {
                        DiffusionCoeff = 1.0e-2;
                    } else {
                        if(DiffusionCoeff < 1.0e3)
                            DiffusionCoeff *= 2;
                    }
                    //Console.WriteLine("   increasing Diffusion: {0}, Alarms : {1}{2}{3}{4}", DiffusionCoeff, MinAlarm ? 1 : 0, Maxalarm ? 1 : 0, GrowAlarm ? 1 : 0, IterAlarm ? 1 : 0);
                }


                //if(jCell == 80 && iIter < 100)
                //    Tecplot.Tecplot.PlotFields(PlotFields, "itt_" + (iIter + 1), "EllipicReinit", iIter + 1, 3);
            }


            return converged;
        }

        /// <summary>
        /// Runge-Kutta scheme uses in local solver.
        /// </summary>
        RungeKuttaScheme RKsch = RungeKuttaScheme.TVD3;

        double[] temp = null;

        void ComputeChangeRate(double[] k, int jCell, BitArray AcceptedMask, SinglePhaseField Phi, VectorField<SinglePhaseField> gradPhi, SpatialOperator.Evaluator Evaluator) {
            int N = this.LevelSetBasis.GetLength(jCell);
            int i0L = this.LevelSetMapping.LocalUniqueCoordinateIndex(0, jCell, 0);
            Debug.Assert(k.Length == N);

            if(temp == null) {
                temp = new double[this.LevelSetMapping.LocalLength];
            }

            gm.GradientUpdate(jCell, AcceptedMask, Phi, gradPhi);

            for(int n = 0; n < N; n++)
                temp[n + i0L] = 0;
#if DEBUG
            double[] oldtemp = temp.CloneAs();
#endif

            Evaluator.Evaluate<double[]>(1.0, 1.0, temp, 0.0, MPIexchange: false);

#if DEBUG
            for(int i = 0; i < oldtemp.Length; i++) {
                if(i < i0L && i >= (i0L + N))
                    Debug.Assert(oldtemp[i] == temp[i]); // check that the locality of the operator is implemented correctly
            }
#endif

            for(int n = 0; n < N; n++)
                k[n] = temp[n + i0L];
        }

        void PerformRKstep(double dt, int jCell, BitArray AcceptedMask, SinglePhaseField Phi, VectorField<SinglePhaseField> gradPhi, SpatialOperator.Evaluator Evaluator) {
            int N = this.LevelSetBasis.GetLength(jCell);
            int i0G = this.LevelSetMapping.GlobalUniqueCoordinateIndex(0, jCell, 0);
            int i0L = this.LevelSetMapping.LocalUniqueCoordinateIndex(0, jCell, 0);

            double[][] k = new double[RKsch.Stages][];
            for(int i = 0; i < RKsch.Stages; i++)
                k[i] = new double[N];

            double[] y0 = Phi.Coordinates.GetRow(jCell);
            Debug.Assert(y0.Length == N);

            ComputeChangeRate(k[0], jCell, AcceptedMask, Phi, gradPhi, Evaluator);

            for(int s = 1; s < RKsch.Stages; s++) {

                Phi.Coordinates.SetRow(jCell, y0);
                for(int r = 0; r < s; r++) {
                    if(RKsch.a[s, r] != 0.0) {
                        Phi.Coordinates.AccRow(jCell, -dt * RKsch.a[s, r], k[r]);
                    }
                }

                ComputeChangeRate(k[s], jCell, AcceptedMask, Phi, gradPhi, Evaluator);
            }

            // next timestep
            for(int s = 0; s < RKsch.Stages; s++) {
                if(RKsch.b[s] != 0.0) {
                    y0.AccV(-RKsch.b[s] * dt, k[s]);
                }
            }
            Phi.Coordinates.SetRow(jCell, y0);
        }

        void PerformArtificialDiffusion(double dt, int jCell, SinglePhaseField Phi, MultidimensionalArray DiffMtx, double[] B) {
            int N = this.LevelSetBasis.GetLength(jCell);
            int i0L = this.LevelSetMapping.LocalUniqueCoordinateIndex(0, jCell, 0);

            // extract data from jCell
            double[] u0 = Phi.Coordinates.GetRow(jCell);

            // RHS = (1/dt)*u0 - B
            double[] RHS = u0;
            RHS.ScaleV(1 / dt);
            RHS.AccV(-1.0, B);

            // LHS = (1/dt)*I + DiffMtx
            MultidimensionalArray LHS = DiffMtx.CloneAs();
            LHS.AccEye(1 / dt);

            // solve
            double[] u1 = new double[N];
            LHS.Solve(u1, RHS);

            // return
            Phi.Coordinates.SetRow(jCell, u1);
        }


    }
}
