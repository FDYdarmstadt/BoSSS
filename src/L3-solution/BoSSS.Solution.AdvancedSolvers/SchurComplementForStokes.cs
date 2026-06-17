using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using NUnit.Framework.Constraints;

using ilPSP;
using ilPSP.LinSolvers;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using System.Diagnostics;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// A global Schur complement solver for Stokes problems.
    /// </summary>
    public class SchurComplementForStokes : ISolverSmootherTemplate, IProgrammableTermination {

        [Serializable]
        public class Config : IterativeSolverConfig {
            override public string Name => "Schur Complement";

            override public string Shortname => "SchurComp";

            [DataMember]
            public ISolverFactory ConvectionDiffusionSolver = new AdvancedSolvers.OrthoMGSchwarzConfig() {
#if DEBUG
                CoarseKickIn = 100 // for debugging purposes, enfoce a later switch to direct solver, so that we have a non-trivial multigrid depth on coarse meshes.
#endif
            };

            override public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new SchurComplementForStokes();
                instance.Init(level);
                return instance;
            }

            
        }

        Config m_config = new Config();


        public int IterationsInNested {
            get;
            private set;
        }


        public int ThisLevelIterations {
            get;
            private set;
        }

        public bool Converged {
            get;
            private set;
        }
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }

        public object Clone() {
            throw new NotImplementedException();
        }

        public void Dispose() {

        }

        public void Init(MultigridOperator op) {
            int D = op.GridData.SpatialDimension;
            int NoOfVars = op.Mapping.NoOfVariables;
            if(NoOfVars != D + 1) {
                throw new NotSupportedException($"expecting {D + 1} variables, got {NoOfVars}");
            }


            // create the sub-multigrid solver
            int[] VelVarIndices = D.ForLoop(d => d);
            var cdMultigridOp = op.GetSubOperator(VelVarIndices);

            // we init the convection-diffusion solver on the top level;
            //  -> this should create a sequence of (recursive) multigrid solvers
            //  -> we pass this sequence along the sequence for the Schur-Complement solver
            var __ConvectionDiffusionSolver = this.m_config.ConvectionDiffusionSolver.CreateInstance(cdMultigridOp);

            InitCoarseRecursive(op, __ConvectionDiffusionSolver, cdMultigridOp);
        }

        public void InitCoarseRecursive(MultigridOperator op, ISolverSmootherTemplate opSolver, MultigridOperator cdMultigridOp) {
            this.ConvectionDiffusionSolver = opSolver;
            int D = op.GridData.SpatialDimension;

            int[] VelVarIndices = D.ForLoop(d => d);
            int[] PrsVarIndices = new int[] { D };
            VelocityDofIndices = op.Mapping.GetSubvectorIndices(VelVarIndices);
            PressureDofIndices = op.Mapping.GetSubvectorIndices(PrsVarIndices);

            // extract the sub-matrices
            PressureGradient = op.OperatorMatrix.GetSubMatrix(VelocityDofIndices, PressureDofIndices);
            VelocityDivergence  = op.OperatorMatrix.GetSubMatrix(PressureDofIndices, VelocityDofIndices);
            M22Matrix = op.OperatorMatrix.GetSubMatrix(PressureDofIndices, PressureDofIndices);



            // recursion
            var coarseSolver = GetCoarseConvectionDiffusionSolver();
            if(coarseSolver != null) {
                InitCoarseRecursive(op.CoarserLevel, coarseSolver, cdMultigridOp.CoarserLevel);
            }
        }


        ISolverSmootherTemplate GetCoarseConvectionDiffusionSolver() {
            var current = this.ConvectionDiffusionSolver;
            if(current is OrthonormalizationMultigrid omg) {
                return omg.CoarserLevelSolver;
            } else {

                return null;
            }
        }

        ISolverSmootherTemplate ConvectionDiffusionSolver;


        long[] VelocityDofIndices;
        long[] PressureDofIndices;


        BlockMsrMatrix PressureGradient;

        BlockMsrMatrix VelocityDivergence;

        /// <summary>
        /// lower right block of the system matrix;
        /// typically, this is zero; otherwise it might contain stabilization terms.
        /// </summary>
        BlockMsrMatrix M22Matrix;




        public void ResetStat() {
            IterationsInNested = 0;
            ThisLevelIterations = 0;
            Converged = false;
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            throw new NotImplementedException();
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}
