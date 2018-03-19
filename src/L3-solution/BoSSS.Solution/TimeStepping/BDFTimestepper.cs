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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Timestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;


namespace BoSSS.Solution.TimeStepping {
    
    /// <summary>
    /// Implicit Timestepping for General Spatial Operators
    /// also supports implicit and Explicit Euler and Crank-Nicholson
    /// </summary>
    public class BDFTimestepper : ITimeStepper {
        /// <summary>
        /// The spatial operator this time stepper is based on
        /// </summary>
        public SpatialOperator Operator {
            get;
            private set;
        }

        /// <summary>
        /// mapping of DG fields that are affected by this timestepper
        /// </summary>
        public CoordinateMapping Mapping {
            get;
            private set;
        }

        /// <summary>
        /// mapping of parameter fields
        /// </summary>
        public CoordinateMapping ParameterMapping {
            get;
            private set;
        }

        /// <summary>
        /// The DG coordinates of the <see cref="BoSSS.Foundation.DGField"/>'s
        /// in <see cref="Mapping"/> (see
        /// <see cref="BoSSS.Foundation.CoordinateMapping.Fields"/>).
        /// </summary>
        public CoordinateVector CurrentState {
            get;
            private set;
        }

       
        /// <summary>
        /// The affine offset b from F(x) = Mx + b
        /// </summary>
        /// <remarks>
        /// If boundary conditions are time-dependent, this vector may change over 
        /// time; Re-calculation can be implemented e.g. in <see cref="BeforeTimeStep"/>.
        /// At a given time <em>t</em>, for given initial values
        /// it is 
        /// the convention is that this vector represents the inhomogeouos b.c.
        /// at time <em>t</em>+ <em>dt</em>, i.e. at the next timestep.
        /// </remarks>
        protected double[] m_AffineOffset1;

        /// <summary>
        /// The Matrix for the SpatialOperator
        /// </summary>
        BlockMsrMatrix SpatialMatrix;

        /// <summary>
        /// An Array of BDSchemes from implicit euler(entry[i]) up to the desired order (entry[0])
        /// </summary>
        private BDFSchemeCoeffs[] TSCchain;

        private int PopulatedStackDepth;
        
        /// <summary>
        /// Coefficients of the currently used scheme
        /// </summary>
        BDFSchemeCoeffs Tsc;

        /// <summary>
        /// The sparse solver used to solve the equation system in
        /// <see cref="PerformTimeStep"/>
        /// </summary>
        /// <remarks>
        /// The Matrix of the system is stored within the sparse solver;
        /// </remarks>
        protected Func<ISparseSolver> SolverFactory;
        private BlockMsrMatrix[] Stack_OpMatrix;
        private double[][] Stack_OpAffine;

        private CoordinateVector[] Stack_u;

        SubGrid subGrid;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialOp">
        /// The Spatial Discretization,
        /// also supports nonconstant operators,
        /// which are linear in the unknown variable
        /// </param>
        /// <param name="UnknownFields"></param>
        /// <param name="ParameterFields"></param>
        /// <param name="BDForder">
        /// The Temporal Discretization Order
        /// >=1 : BDF-Scheme
        /// 0: Explicit Euler
        /// -1: Crank-Nicholson
        /// </param>
        /// <param name="SolverFactory">
        /// The Linear solver for the resulting systems
        /// </param>
        /// <param name="DelayInit">TODO: Delayed Initialization for restarts etc.</param>
        /// <param name="subGrid">TODO: Perform Time-Marching only on Substep</param>
        public BDFTimestepper(SpatialOperator spatialOp, IEnumerable<DGField> UnknownFields, IEnumerable<DGField> ParameterFields, int BDForder, Func<ISparseSolver> SolverFactory, bool DelayInit, SubGrid subGrid = null) {
            using (new ilPSP.Tracing.FuncTrace()) {

                if (spatialOp.ContainsNonlinear) { throw new NotImplementedException("No Inversion of Nonlinear Operators implemented, yet."); };

                if (DelayInit) throw new NotImplementedException();
                if (subGrid != null) throw new NotImplementedException();

                this.Mapping = new CoordinateMapping(UnknownFields.ToArray());
                this.ParameterMapping = new CoordinateMapping(ParameterFields.ToArray());

                // verify input
                // ============
                TimeStepperCommon.VerifyInput(spatialOp, Mapping, ParameterMapping);

                // Initialize
                CurrentState = new CoordinateVector(UnknownFields.ToArray());
                Operator = spatialOp;
                this.SolverFactory = SolverFactory;

                this.subGrid = subGrid;

                // Initialize Matrix
                SpatialMatrix = new BlockMsrMatrix(Mapping);

                TSCchain = BDFCommon.GetChain(BDForder);
                int S = TSCchain[0].S;
                Debug.Assert(S == TSCchain.Length);

                // Set coarsest BDF-Scheme to CrankNicholson, since this is more accurate than Implicit Euler
                Console.WriteLine("Warning! - 1st Initialization Step set to Crank-Nicholson!");
                TSCchain[S-1] = BDFSchemeCoeffs.CrankNicolson();

                InitStacks(Mapping, S);

                // other stuff
                // -----------

                if (!DelayInit)
                    InitTimestepping(true);

            }
        }

        /// <summary>
        /// Initilializes the History-stack for the Unknown u and the Operator Matrix and Affine Part
        /// </summary>
        /// <param name="Mapping"></param>
        /// <param name="StackLength"></param>
        private void InitStacks(CoordinateMapping Mapping, int StackLength) {
            // operator matrix - stack
            // -----------------------

            Stack_OpMatrix = new BlockMsrMatrix[2]; // only required for Crank.Nic. and Expl. Euler,
            Stack_OpAffine = new double[2][]; //      in this case 'theta0' is unequal 0.0.

            Stack_OpMatrix[1] = new BlockMsrMatrix(Mapping);
            Stack_OpAffine[1] = new double[Mapping.LocalLength];


            // Stack for the unknown field
            //----------------------------
            Stack_u = new CoordinateVector[StackLength];
            for (int i = 0; i < StackLength; i++) {
                Stack_u[i] = new CoordinateVector(Mapping.Fields.Select(dgf => dgf.CloneAs()).ToArray());
            }
            Stack_u[0].Clear();
            Stack_u[0].AccV(1.0, CurrentState);
            
            // Initialize the Operator Matrix for TimeStep 0, i.e. Calculate it and Push it to its Position
            UpdateOperatorMatrix();
            PushStacks();
        }


        /// <summary>
        /// final initialization for the BDF timestepper scheme, all necessary timesteps have to be initialized
        /// </summary>
        /// <param name="OpInit"></param>
        private void InitTimestepping(bool OpInit) {


            // operator matrix update
            // ----------------------

            if (OpInit && TSCchain[0].theta0 != 0.0) {             
                if (Stack_OpAffine[0] == null) {
                    Stack_OpAffine[0] = new double[Stack_u[0].Mapping.LocalLength];
                }

                Debug.Assert(Stack_OpMatrix[0].InfNorm() == 0);
                Debug.Assert(Stack_OpAffine[0].L2Norm() == 0);

                Tsc = TSCchain.Last();

                


            }
        }



        /// <summary>
        /// Initialization for a multi-step method, e.g. BDF4 requires 4 timesteps.
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="TimestepNo"></param>
        /// <param name="dt"></param>
        /// <param name="SetTimestep"></param>
        public void MultiInit(double physTime, int TimestepNo, double dt, Action<int, double, DGField[]> SetTimestep) {
            //the code below has not been tested or debugged, yet
            throw new NotImplementedException("Not yet tested");
            using (new FuncTrace()) {
                if (dt <= 0)
                    throw new ArgumentOutOfRangeException();

                int S = TSCchain[0].S;
                Debug.Assert(S == TSCchain.Length);

                for (int iStage = 0; iStage < S; iStage++) {
                    int TimeIndex = -S + iStage + 1;
                    double time = physTime + TimeIndex * dt;
                    SetTimestep(TimeIndex + TimestepNo, time, Stack_u[0].Mapping.Fields.ToArray()); // the push-operation will init the other steps.

                    // all the other stuff (cut-cell-metrics, ...)
                    InitTimestepping(iStage == (S - 1));

                    if (iStage < (S - 1))
                        PushStacks();
                }
            }

        }


        void UpdateOperatorMatrix() {
            Stack_OpMatrix[1].Clear();
            Stack_OpAffine[1].Clear();

            // Assemble matrix and affine offset
            Operator.ComputeMatrixEx(
                Mapping, ParameterMapping, Mapping,
                Stack_OpMatrix[1], Stack_OpAffine[1],
                false,
                 0.0,
                 new EdgeQuadratureScheme(true, subGrid?.AllEdgesMask),
                 new CellQuadratureScheme(true, subGrid?.VolumeMask)
                );

            //Debug.Assert(Stack_OpMatrix[1].InfNorm() > 0);
        }

        public double Time {
            get;
            private set;
        }

        /// <summary>
        /// Performs a single BDF Timestep
        /// This does NOT include any nonlinear iterations
        /// </summary>
        /// <param name="dt"></param>
        /// <returns></returns>
        public double Perform(double dt) {
            if (dt <= 0)
                throw new ArgumentOutOfRangeException();

            double[] RHS;
            BlockMsrMatrix myMatrix;
            AssembleMatrix(dt, out myMatrix, out RHS);
            using (var slv = SolverFactory()  ) {
                slv.DefineMatrix(myMatrix);
                slv.Solve(CurrentState, RHS);
                slv.Dispose();
            }
            return dt;
        }

        private void AssembleMatrix(double dt, out BlockMsrMatrix SystemMatrix, out double[] SystemAffine) {
            // Init Matrix and Affine Part 
            SystemMatrix = new BlockMsrMatrix(Mapping);
            SystemAffine = new double[Mapping.LocalLength];

            // choose TimeStepping-Scheme, based on what has been pushed to the stack,yet
            int Smax = TSCchain[0].S;
            Debug.Assert(Smax == TSCchain.Length);
            Tsc = TSCchain[Smax - PopulatedStackDepth];

            UpdateOperatorMatrix();





            //Implicit Part of RHS
            SystemMatrix.Acc(Tsc.theta1, Stack_OpMatrix[1]);
            SystemAffine.AccV(Tsc.theta1, Stack_OpAffine[1]);



            //Implicit Part of LHS
            SystemMatrix.AccEyeSp(1 / dt);

            // Explicit part of RHS
            Stack_OpMatrix[0].SpMV(-Tsc.theta0, CurrentState, 1.0, SystemAffine);
            SystemAffine.AccV(-Tsc.theta0, Stack_OpAffine[0]);

            //Explicit parts of LHS
            for (int i = 0;  i< Tsc.beta.Length; i++) {
                SystemAffine.AccV(Tsc.beta[i]*1/dt, Stack_u[i]);
            }

            Debug.Assert(SystemMatrix.InfNorm() > 0);
            Debug.Assert(SystemAffine.L2Norm() > 0);

            if (subGrid != null) {
                int[] SubVecIdx = Mapping.GetSubvectorIndices(subGrid, true, new int[] { 0 });
                int L = SubVecIdx.Length;

                for (int i=0; i < L; i++) {
                    SystemMatrix.ClearRow(SubVecIdx[i]);
                    SystemMatrix[SubVecIdx[i], SubVecIdx[i]] = 1;
                    SystemAffine[SubVecIdx[i]] = 0;
                } 
            }

        }


        public void FinishTimeStep() {
            PushStacks();
        }

        private void PushStacks() {
            //increase the Stack length up to the maximum Value
            
            if (PopulatedStackDepth < TSCchain[0].S) {
                PopulatedStackDepth++;
            }
            
            // Push Operator-Part
            Stack_OpMatrix[0] = Stack_OpMatrix[1].CloneAs();
            Stack_OpAffine[0] = Stack_OpAffine[1].CloneAs();
            
            // Push Unknowns
            for (int i = Stack_u.Length-1; i == 1; i--) {
                Stack_u[i].Clear();
                Stack_u[i].Acc(1.0, Stack_u[i - 1]);
            }
            Stack_u[0].Clear();
            Stack_u[0].Acc(1.0, CurrentState);
        }

        public void ResetTime(double NewTime, int timestepNumber) {
            Time = NewTime;
        }
    }
}
