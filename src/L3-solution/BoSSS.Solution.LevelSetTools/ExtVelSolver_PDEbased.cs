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
    /// Solves an extension problem on a cell-by-cell basis
    /// </summary>
    class ExtVelSolver_PDEbased {

        public ExtVelSolver_PDEbased(Basis ExtPropertyBasis) {
            var map = new UnsetteledCoordinateMapping(ExtPropertyBasis);
            m_ExtvelMatrix = new MsrMatrix(map, map);
            m_ExtvelAffine = new double[m_ExtvelMatrix.RowPartitioning.LocalLength];
        }
        
        MsrMatrix m_ExtvelMatrix;
        double[] m_ExtvelAffine;

        /// <summary>
        /// Solution of the extension problem on a single cell in the far-region
        /// </summary>
        /// <param name="Phi">Input;</param>
        /// <param name="GradPhi">Input;</param>
        /// <param name="ExtProperty"></param>
        /// <param name="ExtPropertyMin">Input/Output: lower threshold.</param>
        /// <param name="ExtPropertyMax">Input/Output: upper threshold.</param>
        /// <param name="jCell"></param>
        /// <param name="Accepted"></param>
        /// <param name="signMod"></param>
        public void ExtVelSolve_Far(SinglePhaseField Phi, VectorField<SinglePhaseField> GradPhi, ConventionalDGField ExtProperty, ref double ExtPropertyMin, ref double ExtPropertyMax, int jCell, CellMask Accepted, double signMod) {
            GridData gDat = (GridData)(Phi.GridDat);
            Debug.Assert(signMod.Abs() == 1.0);
            Debug.Assert(ExtPropertyMin <= ExtPropertyMax);
            

            // define cell- and edge-mask for re-compute
            // =========================================

            CellMask cM = new CellMask(gDat, Chunk.GetSingleElementChunk(jCell));

            int[] Edges = gDat.iLogicalCells.Cells2Edges[jCell].CloneAs();
            for(int i = 0; i < Edges.Length; i++) {
                Edges[i] = Math.Abs(Edges[i]) - 1;
            }

            EdgeMask eM = new EdgeMask(gDat, FromIndEnum(Edges, gDat.iLogicalEdges.Count)); // won't scale.

            // solve the linear extension problem for 'jCell', eventually increase
            // diffusion until we are satisfied with the solution
            // ===================================================================
                        
            bool MaximumPrincipleFulfilled = false;
            double DiffusionCoeff = 0; // initially: try without diffusion
            double mini = double.NaN, maxi = double.NaN;
            int count = 0;
            while(MaximumPrincipleFulfilled == false) { // as long as we are satisfied with the solution
                count++;

                // compute operator in 'jCell'
                // ---------------------------

                int N = ExtProperty.Basis.GetLength(jCell);
                int i0G = ExtProperty.Mapping.GlobalUniqueCoordinateIndex(0, jCell, 0);
                int i0L = ExtProperty.Mapping.GlobalUniqueCoordinateIndex(0, jCell, 0);

                for(int n = 0; n < N; n++) {
                    this.m_ExtvelMatrix.ClearRow(i0G + n);
                    this.m_ExtvelAffine[i0L + n] = 0;
                }

                double penaltyBase = ((double)(ExtProperty.Basis.Degree + 1)).Pow2();
                var Flux = new ExtensionVelocityForm(Accepted.GetBitMask(), signMod, penaltyBase, gDat.Cells.h_min, jCell, DiffusionCoeff);
                var op = (Flux).Operator(DegreeOfNonlinearity: 2);

                // increase diffusion coefficient for next cycle
                if(DiffusionCoeff == 0)
                    DiffusionCoeff = 1.0e-3; // should this be minus or plus?
                else
                    DiffusionCoeff *= 10;

                op.ComputeMatrixEx(ExtProperty.Mapping, ArrayTools.Cat<DGField>(new DGField[] { ExtProperty, Phi }, GradPhi), ExtProperty.Mapping,
                    this.m_ExtvelMatrix, this.m_ExtvelAffine,
                    OnlyAffine: false,
                    volQuadScheme: (new CellQuadratureScheme(true, cM)),
                    edgeQuadScheme: (new EdgeQuadratureScheme(true, eM)),
                    ParameterMPIExchange: false);

                // extract operator matrix and RHS
                // -------------------------------

                // the matrix must only have entries in the block-diagonal!

                MultidimensionalArray Mtx = MultidimensionalArray.Create(N, N);
                //MultidimensionalArray rhs = MultidimensionalArray.Create(N);
                double[] rhs = new double[N];

                for(int n = 0; n < N; n++) {
#if DEBUG
                    int Lr;
                    int[] row_cols = null;
                    double[] row_vals = null;
                    Lr = this.m_ExtvelMatrix.GetRow(i0G + n, ref row_cols, ref row_vals);
                    for(int lr = 0; lr < Lr; lr++) {
                        int ColIndex = row_cols[lr];
                        double Value = row_vals[lr];
                        Debug.Assert((ColIndex >= i0G && ColIndex < i0G + N) || (Value == 0.0), "Matrix is expected to be block-diagonal.");
                    }
#endif
                    for(int m = 0; m < N; m++) {
                        Mtx[n, m] = this.m_ExtvelMatrix[i0G + n, i0G + m];
                    }
                    rhs[n] = -this.m_ExtvelAffine[i0L + n];
                }
                
                // Solve the system, i.e. the local extension-velocity equation
                // ------------------------------------------------------------
                
                double[] ep = new double[N];
                Mtx.Solve(ep, rhs);

                for(int n = 0; n < N; n++) {
                    ExtProperty.Coordinates[jCell, n] = ep[n];
                }

                // detect violation of maximum principle
                // -------------------------------------

                ExtProperty.GetExtremalValuesInCell(out mini, out maxi, jCell);

                // define relaxed bounds...
                double compExtPropertyMin = ExtPropertyMin - (1.0e-8 + ExtPropertyMax - ExtPropertyMin) * 0.2;
                double compExtPropertyMax = ExtPropertyMax + (1.0e-8 + ExtPropertyMax - ExtPropertyMin) * 0.2;

                // test if extension velocity solution is within bounds
                MaximumPrincipleFulfilled = (mini >= compExtPropertyMin) && (maxi <= compExtPropertyMax);

                if(count > 5)
                    break;
            }
            if(count > 4)
                Console.WriteLine(" ExtVel, cell #{0}, Diff coeff {1}, extremal p. holds? {2} (min/max soll = {3:e4}/{4:e4}, ist = {5:e4}/{6:e4})", jCell, DiffusionCoeff, MaximumPrincipleFulfilled, ExtPropertyMin, ExtPropertyMax, mini, maxi);

            // record maxima and minima
            // ========================
            ExtPropertyMax = Math.Max(ExtPropertyMax, maxi);
            ExtPropertyMin = Math.Min(ExtPropertyMin, mini);
        }

        /// <summary>
        /// Extension-Velocity form with optional diffusion.
        /// </summary>
        class ExtensionVelocityForm : EllipticExtension.EllipticExtVelForm {
            public ExtensionVelocityForm(BitArray _AcceptedBitmask, double __signMod, double __penaltyBase, MultidimensionalArray __h_min, int __jCell, double __DiffusionCoeff):base(__penaltyBase, __DiffusionCoeff,null) {
                this.AcceptedBitmask = _AcceptedBitmask;
                this.signMod = __signMod;
                this.penaltyBase = __penaltyBase;
                this.jCell = __jCell;
                this.DiffusionCoeff = __DiffusionCoeff;
                this.h_min = __h_min;
            }

            double penaltyBase;
            double signMod;
            BitArray AcceptedBitmask;
            int jCell;
            double DiffusionCoeff;


            public new IList<string> ParameterOrdering {
                get {
                    return new string[] { "s_accepted", "Phi", "dPhi_dx0", "dPhi_dx1" };
                }
            }

            public new TermActivationFlags BoundaryEdgeTerms {
                get {
                    return TermActivationFlags.None;
                }
            }

            public new TermActivationFlags InnerEdgeTerms {
                get {
                    //return TermActivationFlags.GradUxV;
                    return TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradUxV;
                }
            }

            MultidimensionalArray h_min;

            public new double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] sA, double[,] Grad_sA, double _vA, double[] _Grad_vA) {
                return 0.0;
            }

            public override double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIn, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIn, double[] Grad_vOut) {
                Debug.Assert((this.AcceptedBitmask[inp.jCellIn] && this.AcceptedBitmask[inp.jCellOut]) == false, "In and Out - cell accepted: cannot be!");
                int D = inp.D;
                Debug.Assert(signMod.Abs() == 1);

                double PhiM = (inp.Parameters_IN[1] + inp.Parameters_OUT[1]) * 0.5 * signMod;
                double phi_sign = Math.Sign(PhiM);

                double Acc = 0;

                // levelSet*gradPhi*n

                double GradUInGradPhi = 0;
                double GradUOutGradPhi = 0;
                double GradVInGradPhi = 0;
                double GradVOutGradPhi = 0;
                double nGradPhiIn = 0;
                double nGradPhiOut = 0;
                for (int d = 0; d < inp.D; d++) {
                    GradUInGradPhi += Grad_uIn[0, d] * inp.Parameters_IN[d+2];
                    GradVInGradPhi += Grad_vIn[d] * inp.Parameters_IN[d+2];
                    GradVOutGradPhi += Grad_vOut[d] * inp.Parameters_OUT[d+2];
                    GradUOutGradPhi += Grad_uOut[0, d] * inp.Parameters_OUT[d+2];
                    nGradPhiIn += inp.Normale[d] * inp.Parameters_IN[d+2];
                    nGradPhiOut += inp.Normale[d] * inp.Parameters_OUT[d+2];
                }


                //penaltyFactor in normal direction, see Ern,Stephansen,Zunino 2008 - eqns 2.12 - 2.14
                double deltaKA = 0;
                double deltaKB = 0;
                for (int d = 0; d < D; d++) {
                    deltaKA += inp.Parameters_IN[d+2] * inp.Normale[d];
                    deltaKB += inp.Parameters_OUT[d+2] * inp.Normale[d];
                }
                deltaKA *= deltaKA;
                deltaKB *= deltaKB;
                //harmonic mean
                double Viscosity = (deltaKA == 0 && deltaKB == 0) ? 0 : deltaKA * deltaKB / (deltaKA + deltaKB);

                if (this.AcceptedBitmask[inp.jCellIn] == true && this.AcceptedBitmask[inp.jCellOut] == false) {
                    // IN-cell is accepted // OUT-cell should be computed
                    // => the boundary value is given as IN-parameter
                    // => flux penalizes the OUT-Cell

                    Debug.Assert(inp.jCellOut == jCell);
                    double penalty = penaltyBase / h_min[inp.jCellOut];

                    Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (-vOut) * 0.5 * (nGradPhiIn + nGradPhiOut);
                    Acc += (GradVOutGradPhi) * (uIn[0] - uOut[0]) * 0.5 * (nGradPhiIn + nGradPhiOut);

                    Acc -= (uIn[0] - uOut[0]) * (-vOut) * penalty * Viscosity;


                }
                else if (this.AcceptedBitmask[inp.jCellIn] == false && this.AcceptedBitmask[inp.jCellOut] == true) {
                    // ... vice-versa

                    Debug.Assert(inp.jCellIn == jCell);

                    double penalty = penaltyBase / h_min[inp.jCellIn];

                    Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (vIn) * 0.5 * (nGradPhiIn + nGradPhiOut);
                    Acc += (GradVInGradPhi) * (uIn[0] - uOut[0]) * 0.5 * (nGradPhiIn + nGradPhiOut);

                    Acc -= (uIn[0] - uOut[0]) * (vIn) * penalty * Viscosity;

                }
                else {
                    /// Do nothing
                }
                return -Acc;
            }

            public new double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                Debug.Assert(GradU.GetLength(0) == 1);
                Debug.Assert(GradU.GetLength(1) == D);
                double Acc = 0;

                double GradPhiGradU = 0;
                double GradPhiGradV = 0;
                for (int d = 0; d < D; d++) {
                    GradPhiGradU += cpv.Parameters[d+2] * GradU[0, d];
                    GradPhiGradV += cpv.Parameters[d+2] * GradV[d];
                }
                Acc += GradPhiGradU * GradPhiGradV;

                return Acc;
            }
        }
             
        static BitArray FromIndEnum(IEnumerable<int> Indices, int NoOfE) {
            BitArray b = new BitArray(NoOfE);
            foreach(int i in Indices) {
                Debug.Assert(b[i] == false);
                b[i] = true;
            }
            return b;
        }

        public void ExtVelSolve_Cut(SinglePhaseField Phi, VectorField<SinglePhaseField> GradPhi, ConventionalDGField ExtProperty, ref double ExtPropertyMin, ref double ExtPropertyMax, int jCell, CellMask Accepted, double signMod) {


        }

       

    }
}
