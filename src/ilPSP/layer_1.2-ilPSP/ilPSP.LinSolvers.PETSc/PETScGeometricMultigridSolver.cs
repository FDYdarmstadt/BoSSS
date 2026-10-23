/* =======================================================================
Copyright 2026 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// One explicit fine-to-coarse transition in a BoSSS geometric multigrid hierarchy.
    /// </summary>
    /// <remarks>
    /// <see cref="Interpolation"/> maps the coarse vector into the fine space and
    /// <see cref="Restriction"/> maps the fine vector into the coarse space. They
    /// are deliberately independent matrices: BoSSS change-of-basis transforms can
    /// make restriction different from interpolation transposed.
    /// </remarks>
    public sealed class PETScMultigridTransition {

        /// <summary>Physically reassembled operator for the coarse level.</summary>
        public IMutableMatrixEx CoarseOperator { get; set; }

        /// <summary>Explicit coarse-to-fine interpolation operator.</summary>
        public IMutableMatrixEx Interpolation { get; set; }

        /// <summary>Explicit fine-to-coarse restriction operator.</summary>
        public IMutableMatrixEx Restriction { get; set; }
    }

    /// <summary>
    /// A finest-to-coarsest sequence of explicit BoSSS multigrid transitions.
    /// </summary>
    public sealed class PETScGeometricMultigridHierarchy {
        readonly List<PETScMultigridTransition> m_Transitions = new List<PETScMultigridTransition>();

        /// <summary>Transitions ordered from the finest level towards the coarsest one.</summary>
        public IList<PETScMultigridTransition> Transitions => m_Transitions;
    }

    /// <summary>
    /// Direct-ABI PETSc geometric multigrid solver using an explicit <c>PCMG</c> hierarchy.
    /// </summary>
    /// <remarks>
    /// The finest matrix is supplied through <see cref="PETScSolver.DefineMatrix"/>.
    /// Each transition supplies BoSSS's explicitly reassembled coarse operator and
    /// explicit restriction/interpolation matrices; PETSc never constructs a Galerkin
    /// coarse operator. All hierarchy matrices are copied into PETSc, making their
    /// ownership independent from the source BoSSS matrices.
    /// </remarks>
    public sealed class PETScGeometricMultigridSolver : PETScSolver {
        // PCMGGalerkinType.PC_MG_GALERKIN_NONE in PETSc 3.26.
        const int PcMgGalerkinNone = 3;
        readonly List<PETScMPIAIJMatrix> m_HierarchyMatrices = new List<PETScMPIAIJMatrix>();

        /// <summary>PETSc command-line-style outer-KSP and PCMG options.</summary>
        public string Options { get; set; }

        /// <summary>Explicit BoSSS geometric multigrid hierarchy.</summary>
        public PETScGeometricMultigridHierarchy Hierarchy { get; set; }

        /// <inheritdoc/>
        protected override PETScMatrixInput SelectMatrixInput(int mpiSize) {
            // Fine-level transfer remains configurable. All P/R/coarse matrices are
            // copied because their source matrices are temporary hierarchy products.
            return base.SelectMatrixInput(mpiSize);
        }

        /// <inheritdoc/>
        protected override void ConfigureKsp() {
            if(Hierarchy == null || Hierarchy.Transitions.Count == 0)
                throw new InvalidOperationException("PETSc geometric multigrid requires at least one fine-to-coarse transition.");

            Check(m_Native.GetKspPreconditioner(m_Ksp, out IntPtr pc), "PETSc KSPGetPC");
            Check(m_Native.SetPcType(pc, "mg"), "PETSc PCSetType(mg)");
            ConfigureHierarchy(pc);
            ApplyKspOptions(Options);
            Check(m_Native.ObjectTypeCompare(pc, "mg", out int isMultigrid), "PETSc PetscObjectTypeCompare(PCMG)");
            if(isMultigrid == 0)
                throw new ArgumentException("PETSc geometric multigrid requires '-pc_type mg' when PETSc options specify a PC type.", nameof(Options));
        }

        void ConfigureHierarchy(IntPtr pc) {
            int numberOfLevels = checked(Hierarchy.Transitions.Count + 1);
            // Match the historical native wrapper exactly. PCMG uses this
            // communicator only while constructing internal level state; NULL
            // is PETSc's supported "use the PC communicator" value and avoids
            // passing BoSSS's foreign MPI_Comm representation through this API.
            Check(m_Native.SetPcMgLevels(pc, numberOfLevels, IntPtr.Zero), "PETSc PCMGSetLevels");
            // BoSSS provides all coarse operators explicitly below. Do not let
            // PETSc replace them with a Galerkin R*A*P construction.
            Check(m_Native.SetPcMgGalerkin(pc, PcMgGalerkinNone), "PETSc PCMGSetGalerkin(none)");

            for(int transitionIndex = 0; transitionIndex < Hierarchy.Transitions.Count; transitionIndex++) {
                PETScMultigridTransition transition = Hierarchy.Transitions[transitionIndex]
                    ?? throw new ArgumentException($"PETSc multigrid transition {transitionIndex} is null.", nameof(Hierarchy));
                if(transition.CoarseOperator == null || transition.Interpolation == null || transition.Restriction == null)
                    throw new ArgumentException($"PETSc multigrid transition {transitionIndex} must supply coarse operator, interpolation, and restriction matrices.", nameof(Hierarchy));

                ValidateTransition(transition, transitionIndex);
                int finePetscLevel = numberOfLevels - 1 - transitionIndex;
                PETScMPIAIJMatrix interpolation = CreateCopiedMatrix(transition.Interpolation);
                PETScMPIAIJMatrix restriction = CreateCopiedMatrix(transition.Restriction);
                PETScMPIAIJMatrix coarseOperator = CreateCopiedMatrix(transition.CoarseOperator);
                try {
                    Check(m_Native.SetPcMgInterpolation(pc, finePetscLevel, interpolation.Handle), "PETSc PCMGSetInterpolation");
                    Check(m_Native.SetPcMgRestriction(pc, finePetscLevel, restriction.Handle), "PETSc PCMGSetRestriction");
                    Check(m_Native.GetPcMgSmoother(pc, finePetscLevel - 1, out IntPtr coarseKsp), "PETSc PCMGGetSmoother");
                    Check(m_Native.SetKspOperators(coarseKsp, coarseOperator.Handle, coarseOperator.Handle), "PETSc KSPSetOperators(coarse level)");
                    m_HierarchyMatrices.Add(interpolation);
                    m_HierarchyMatrices.Add(restriction);
                    m_HierarchyMatrices.Add(coarseOperator);
                } catch {
                    interpolation.Dispose();
                    restriction.Dispose();
                    coarseOperator.Dispose();
                    throw;
                }
            }
        }

        PETScMPIAIJMatrix CreateCopiedMatrix(IMutableMatrixEx matrix) {
            return PETScMPIAIJMatrix.Create(m_Native, m_Communicator, PETScMatrixInput.CopiedCsr, matrix);
        }

        static void ValidateTransition(PETScMultigridTransition transition, int transitionIndex) {
            long fineRows = transition.Interpolation.RowPartitioning.TotalLength;
            long coarseRows = transition.CoarseOperator.RowPartitioning.TotalLength;
            if(transition.CoarseOperator.RowPartitioning.TotalLength != transition.CoarseOperator.ColPartitioning.TotalLength)
                throw new ArgumentException($"PETSc multigrid coarse operator {transitionIndex} is not square.");
            if(transition.Interpolation.RowPartitioning.TotalLength != transition.Restriction.ColPartitioning.TotalLength
                || transition.Interpolation.ColPartitioning.TotalLength != transition.Restriction.RowPartitioning.TotalLength
                || transition.Interpolation.ColPartitioning.TotalLength != coarseRows
                || transition.Restriction.RowPartitioning.TotalLength != coarseRows
                || fineRows != transition.Restriction.ColPartitioning.TotalLength)
                throw new ArgumentException($"PETSc multigrid transition {transitionIndex} has incompatible interpolation/restriction dimensions.");
        }

        /// <inheritdoc/>
        protected override void DestroyNativeObjects() {
            try {
                base.DestroyNativeObjects();
            } finally {
                for(int index = m_HierarchyMatrices.Count - 1; index >= 0; index--)
                    m_HierarchyMatrices[index].Dispose();
                m_HierarchyMatrices.Clear();
            }
        }
    }
}
