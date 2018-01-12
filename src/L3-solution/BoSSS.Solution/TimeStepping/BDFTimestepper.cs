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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Timestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;


namespace BoSSS.Solution.TimeStepping {
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
        protected Func<ISparseSolver> solver;
        private BlockMsrMatrix[] Stack_OpMatrix;
        private double[][] Stack_OpAffine;

        private CoordinateVector[] Stack_u;

        public BDFTimestepper(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int BDForder, Func<ISparseSolver> solver, bool DelayInit) {
            using (new ilPSP.Tracing.FuncTrace()) {

                // verify input
                // ============
                TimeStepperCommon.VerifyInput(spatialOp, Fieldsmap, Parameters);

                this.Mapping = Fieldsmap;
                this.ParameterMapping = Parameters;
                CurrentState = new CoordinateVector(Mapping);

                Operator = spatialOp;

                this.solver = solver;

                // Initialize Matrix
                SpatialMatrix = new BlockMsrMatrix(Mapping);
                
                TSCchain = BDFCommon.GetChain(BDForder);
                int S = TSCchain[0].S;
                Debug.Assert(S == TSCchain.Length);

                InitStacks(Mapping, S);

                // other stuff
                // -----------

                if (!DelayInit)
                    InitTimestepping(true);

            }
        }

        private void InitStacks(CoordinateMapping Mapping, int StackLength) {
            // operator matrix - stack
            // -----------------------

            Stack_OpMatrix = new BlockMsrMatrix[2]; // only required for Crank.Nic. and Expl. Euler,
            Stack_OpAffine = new double[2][]; //      in this case 'theta0' is unequal 0.0.

            Stack_OpMatrix[0] = new BlockMsrMatrix(Mapping);
            Stack_OpAffine[0] = new double[Mapping.LocalLength];


            // Stack for the unknown field
            //----------------------------
            Stack_u = new CoordinateVector[StackLength];
            for (int i = 0; i < StackLength; i++) {
                Stack_u[i] = new CoordinateVector(Mapping.ToArray());
            }
            Stack_u[0].Clear();
            Stack_u[0].AccV(1.0, CurrentState);
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
                ComputeOperatorMatrix();

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


        void ComputeOperatorMatrix() {
            

            // Assemble matrix and affine offset
            Operator.ComputeMatrixEx(
                Mapping, ParameterMapping, Mapping,
                Stack_OpMatrix[0], Stack_OpAffine[0],
                false
                );

            Debug.Assert(Stack_OpMatrix[0].InfNorm() > 0);
            //Debug.Assert(Stack_OpAffine[0].L2Norm() > 0);

        }

        public double Time {
            get;
            private set;
        }

        public double Perform(double dt) {
            if (dt <= 0)
                throw new ArgumentOutOfRangeException();

            double[] RHS;
            BlockMsrMatrix myMatrix;
            AssembleMatrix(dt, out myMatrix, out RHS);
            using (var slv = solver()) {
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

            ComputeOperatorMatrix();





            //Implicit Part of RHS
            SystemMatrix.Acc(Tsc.theta1, Stack_OpMatrix[0]);
            SystemAffine.AccV(Tsc.theta1, Stack_OpAffine[0]);

            //Implicit Part of LHS
            SystemMatrix.AccEyeSp(1 / dt);

            // Explicit part of RHS
            Stack_OpMatrix[1].SpMV(Tsc.theta0, Stack_u[1], 1.0, SystemAffine);
            SystemAffine.AccV(Tsc.theta0, Stack_OpAffine[1]);

            //Explicit parts of LHS
            for (int i = 0;  i< Tsc.beta.Length; i++) {
                SystemAffine.AccV(Tsc.beta[i], Stack_u[i]);
            }
        }

        public void FinishTimeStep() {
            PushStacks();
        }

        private void PushStacks() {
            //increase the Stack length up to the maximum Value
            PopulatedStackDepth++;
            if (PopulatedStackDepth > TSCchain[0].S) {
                PopulatedStackDepth = TSCchain[0].S;
            }
            
            // Push Operator-Part
            Stack_OpMatrix[1] = Stack_OpMatrix[0].CloneAs();
            Stack_OpAffine[1] = Stack_OpAffine[0].CloneAs();
            
            // Push Unknowns
            for (int i = 1; i< Stack_u.Length; i++) {
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
