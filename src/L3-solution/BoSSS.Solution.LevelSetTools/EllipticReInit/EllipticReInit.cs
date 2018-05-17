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
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.EllipticReInit {
    /// <summary>
    /// Potential Functions - Used for Control Files
    /// </summary>
    public enum ReInitPotential {
        /// <summary>
        /// Standard Potential: (1-s)^2
        /// </summary>
        BastingSingleWell,
        
        /// <summary>
        /// Allows gradient to become zero: s^2*(1-s)^2
        /// leads to most stable results
        /// </summary>
        BastingDoubleWell,

        /// <summary>
        /// Trial Option for Ideas
        /// </summary>
        P4DoubleWell, // Diverges without Underrelaxation

        /// <summary>
        /// (1-s)^2, but limited to Near Region,
        /// Elswhere heat equation is solved
        /// </summary>
        SingleWellNear,

        /// <summary>
        /// (1-s)^2, but limited to CutCells,
        /// Elswhere allows gradient to become zero: s^2*(1-s)^2
        /// leads to most stable results
        /// </summary>
        SingleWellOnCutDoubleWellElse
    }



    /// <summary>
    /// Elliptic ReInit
    /// Based on the Ideas by 
    /// C. Basting and D. Kuzmin, 
    /// “A minimization-based finite element formulation for interface-preserving level set reinitialization”,
    /// Computing, vol. 95, no. 1, pp. 13–25, Dec. 2012.
    /// </summary>
    public class EllipticReInit : IReinitialisationAlgorithm {
        int D;

        EllipticReInitAlgoControl Control;

        private LevelSetTracker LevelSetTracker;
        SpatialOperator Operator_bulk;
        XSpatialOperator Operator_interface;
        SpatialOperator Operator_RHS;
        SinglePhaseField Phi;

        SinglePhaseField Residual;
        SinglePhaseField OldPhi;
        SinglePhaseField NewPhi;
        public VectorField<SinglePhaseField> LevelSetGradient;
        public VectorField<SinglePhaseField> MeanLevelSetGradient;

        double[] OldDirection;
        double[] NewDirection;
        bool UpdateDirection = true;

        SinglePhaseField[] parameterFields;

        double ConvergenceCriterion;
        int MaxIteration;
        double underrelaxation;


        public MsrMatrix LHSMatrix
        {
            get { return OpMatrix; }
        }

        MsrMatrix OpMatrix;
        double[] OpAffine;
        public SinglePhaseField RHSField;


        MsrMatrix OpMatrix_interface;
        double[] OpAffine_interface;

        MsrMatrix OpMatrix_bulk;
        double[] OpAffine_bulk;

        MsrMatrix SubMatrix;
        double[] SubRHS;
        double[] SubSolution;
        int[] SubVecIdx;

        ISparseSolver slv;

        /// <summary>
        /// Perform ReInit
        /// </summary>
        /// <returns>
        /// The Number of Iterations and 
        /// The Change-Rate in the ReInit in the least Iteration
        /// </returns>
        public Tuple<int,double> ReInitializeWithOutput(LevelSet LS = null, SubGrid Restriction = null) {
            if (LS != null) throw new NotImplementedException("TODO");

            if (Control.FastMarchingPrecond) {
                //Build FastMarching Solver
                Reinit.FastMarch.FastMarchReinit ReInit = new Reinit.FastMarch.FastMarchReinit(Phi.Basis);

                //Calculate LevelSet
                ReInit.FirstOrderReinit(LevelSetTracker, "A");
            }



            UpdateOperators(Restriction);
            double Residual = double.MaxValue;
            int IterationCounter = 0;
            double OldResidual = double.MaxValue;
            bool EndIteration = false;
            if (Control.PrintIterations) Console.WriteLine("Performing ReInit");
            while (!EndIteration) {
                ReInitSingleStep(out Residual, Restriction);
                IterationCounter++;
                if (Control.PrintIterations) Console.WriteLine("EllipticReInit:Step {0} \t ChangeRate-L2 {1}", IterationCounter, Residual);
                if (Residual < ConvergenceCriterion || IterationCounter >= MaxIteration) EndIteration = true;
                if (Residual > 5*OldResidual) {
                    Console.WriteLine("!!!!!!! \t Iteration Diverged, ABORTING !!!! ");
                    EndIteration = true;
                }
                OldResidual = Residual;
            }
            slv.Dispose();
            return new Tuple<int, double> (IterationCounter,Residual);
        }

        public void ReInitialize(LevelSet LS = null, SubGrid Restriction = null){
            ReInitializeWithOutput(LS, Restriction);
        }

            // Local Variables for Iteration

            // <summary>
            // Counter for Iteration Steps
            // </summary>


            //double OldResidual = double.MaxValue;
            //int divergencecounter = 0;
            ///// <summary>
            ///// Checks for Reaching Max. Number of Iterations and Divergence of Algorithm
            ///// </summary>
            ///// <param name="Residual">Change Rate of the Algorithm</param>
            ///// <returns>Reaching Max Iterations, Aborts when diverged</returns>
            //public bool CheckAbortCriteria(double Residual, int IterationCounter) {
            //    if (Residual <= ConvergenceCriterion) {
            //        Console.WriteLine("EllipticReInit converged after {0} Iterations ", IterationCounter);
            //        return true;
            //    }
            //    if (Residual >= OldResidual) divergencecounter++;
            //    else divergencecounter = 0;
            //    if (IterationCounter >= MaxIteration) {
            //        Console.WriteLine("Elliptic Reinit Max Iterations Reached");
            //        return true;
            //    };
            //    if (divergencecounter > MaxIteration / 2) {
            //        Console.WriteLine("Elliptic Reinit diverged - Aborting");
            //        throw new ApplicationException();
            //    }

            //    OldResidual = Residual;
            //    IterationCounter++;
            //    return false;
            //}


        //bool PreviouslyOnSubgrid = false;

        /// <summary>
        /// Updates the Operator Matrix after level-set motion
        /// </summary>
        /// <param name="Restriction">
        /// The subgrid, on which the ReInit is performed
        /// </param>
        /// <param name="IncludingInterface">
        /// !! Not yet functional !!
        /// True, if the subgrid contains the interface, this causes all external edges of the subgrid to be treated as boundaries
        /// False, for the rest of the domain, thus the flux to the adjacent cells wil be evaluated
        /// </param>
        public void UpdateOperators(SubGrid Restriction = null, bool IncludingInterface = true) {
            if (!IncludingInterface) { throw new NotImplementedException("Untested, not yet functional!"); }
            using (new FuncTrace()) {
                //using (var slv = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
                //using (var slv = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
                //using (var slv = new ilPSP.LinSolvers.HYPRE.GMRES()) {

                if (Control.Upwinding) {
                    OldPhi.Clear();
                    OldPhi.Acc(1.0, Phi);
                    //Calculate
                    LevelSetGradient.Clear();
                    LevelSetGradient.Gradient(1.0, Phi, Restriction?.VolumeMask);
                    //LevelSetGradient.Gradient(1.0, Phi);
                    
                    //LevelSetGradient.GradientByFlux(1.0, Phi);
                    MeanLevelSetGradient.Clear();
                    MeanLevelSetGradient.AccLaidBack(1.0, LevelSetGradient, Restriction?.VolumeMask);
                    //MeanLevelSetGradient.AccLaidBack(1.0, LevelSetGradient);
                }

                if (slv != null) slv.Dispose();

                slv = Control.solverFactory();

                OpMatrix_interface.Clear();
                OpAffine_interface.Clear();


                /// Build the Quadrature-Scheme for the interface operator
                /// Note: The HMF-Quadrature over a surface is formally a volume quadrature, since it uses the volume quadrature nodes.
                XSpatialOperatorExtensions.ComputeMatrixEx(Operator_interface,
                //Operator_interface.ComputeMatrixEx(
                    LevelSetTracker,
                    Phi.Mapping,
                    null,
                    Phi.Mapping,
                    OpMatrix_interface,
                    OpAffine_interface,
                    false,
                    0,
                    false,
                    subGrid:Restriction,
                    whichSpc: LevelSetTracker.GetSpeciesId("A")
                    );

                // Regenerate OpMatrix for subgrid -> adjacent cells must be trated as boundary
                if (Restriction != null){

                    OpMatrix_bulk.Clear();
                    OpAffine_bulk.Clear();

                    //Operator_bulk.ComputeMatrix(
                    //    Phi.Mapping,
                    //    parameterFields,
                    //    Phi.Mapping,
                    //    OpMatrix_bulk, OpAffine_bulk,
                    //    OnlyAffine: false, sgrd: Restriction);
                    EdgeQuadratureScheme edgescheme;
                    //if (Control.Upwinding) {
                    //    edgescheme = new EdgeQuadratureScheme(true, IncludingInterface ? Restriction.AllEdgesMask : null);
                    //}
                    //else {
                        edgescheme = new EdgeQuadratureScheme(true, IncludingInterface ? Restriction.InnerEdgesMask : null);
                    //}
                    Operator_bulk.ComputeMatrixEx(Phi.Mapping,
                            parameterFields,
                            Phi.Mapping, OpMatrix_bulk, OpAffine_bulk, false, 0,
                            edgeQuadScheme: edgescheme,
                            volQuadScheme: new CellQuadratureScheme(true, IncludingInterface ? Restriction.VolumeMask : null)
                            );
                    //PreviouslyOnSubgrid = true;
                }
                // recalculate full Matrix
                //else if (PreviouslyOnSubgrid) {
                else { 
                    OpMatrix_bulk.Clear();
                    OpAffine_bulk.Clear();


                    Operator_bulk.ComputeMatrixEx(Phi.Mapping,
                        parameterFields,
                        Phi.Mapping, OpMatrix_bulk, OpAffine_bulk, false, 0
                        );
                    //PreviouslyOnSubgrid = false;
                }


                /// Compose the Matrix
                /// This is symmetric due to the symmetry of the SIP and the penalty term
                OpMatrix.Clear();
                OpMatrix.Acc(1.0, OpMatrix_bulk);
                OpMatrix.Acc(1.0, OpMatrix_interface);
                OpMatrix.AssumeSymmetric = !Control.Upwinding;
                //OpMatrix.AssumeSymmetric = false;

                /// Compose the RHS of the above operators. (-> Boundary Conditions)
                /// This does NOT include the Nonlinear RHS, which will be added later
                OpAffine.Clear();
                OpAffine.AccV(1.0, OpAffine_bulk);
                OpAffine.AccV(1.0, OpAffine_interface);


#if Debug
                ilPSP.Connectors.Matlab.BatchmodeConnector matlabConnector;
                matlabConnector = new BatchmodeConnector();
#endif

                if (Restriction != null) {
                    SubVecIdx = Phi.Mapping.GetSubvectorIndices(Restriction, true, new int[] { 0 });
                    int L = SubVecIdx.Length;
                    SubMatrix = new MsrMatrix(L);
                    SubRHS = new double[L];
                    SubSolution = new double[L];

                    OpMatrix.AccSubMatrixTo(1.0, SubMatrix, SubVecIdx, default(int[]), SubVecIdx, default(int[]));

                    slv.DefineMatrix(SubMatrix);
#if Debug
                    Console.WriteLine("ConditionNumber of ReInit-Matrix is " + SubMatrix.condest().ToString("E"));
#endif
                }
                else {
                    slv.DefineMatrix(OpMatrix);
#if Debug
                    Console.WriteLine("ConditionNumber of ReInit-Matrix is " + OpMatrix.condest().ToString("E"));
#endif
                }



                
            }
        }




        
        
        
        /// <summary>
        /// 
        /// </summary>
        public void AnalyzeOperators(out MultidimensionalArray ret) {
            ret = MultidimensionalArray.Create(1, 7);
            // MultidimensionalArray ret = MultidimensionalArray.Create(1, 4);
            Console.WriteLine("Calling MATLAB/Octave...");
            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(OpMatrix, "OpMatrix");
                bmc.Cmd("condNoOpMatrix = condest(OpMatrix)");
                bmc.Cmd("OpRank = rank(full(OpMatrix))");
                bmc.Cmd("OpSize = size(OpMatrix)");
                bmc.Cmd("eigiMaxi = eigs(OpMatrix,1,'lm')");
                bmc.Cmd("eigiMini = eigs(OpMatrix,1,'sm')");
                bmc.Cmd("lasterr");
                bmc.Cmd("[V,r]=chol(OpMatrix);");
                bmc.Cmd("ret = [condNoOpMatrix, eigiMaxi, eigiMini, r, OpRank, OpSize]");
                bmc.GetMatrix(ret, "ret");

                bmc.Execute(true);
            }

            double condNoOpMatrix = ret[0, 0];
            double eigiMaxi = ret[0, 1];
            double eigiMini = ret[0, 2];
            double posDef = ret[0, 3] == 0 ? 1 : 0;

            Console.WriteLine("Condition number operator: {0:0.####E-00}", condNoOpMatrix);

            if (posDef == 0.0)
                Console.WriteLine("WARNING: Operator matrix is not positive definite.");
            else
                Console.WriteLine("Good news: Operator matrix seems to be positive definite.");

        }

        /// <summary>
        /// One Iteration of the ReInitialization
        /// Operators must be built first
        /// </summary>
        /// <param name="ChangeRate">
        /// L2-Norm of the Change-Rate in the level set in this reinit step
        /// </param>
        /// <param name="Restriction">
        /// The subgrid, on which the ReInit is performed
        /// </param>
        /// <param name="IncludingInterface">
        /// !! Not yet functional !!
        /// True, if the subgrid contains the interface, this causes all external edges of the subgrid to be treated as boundaries
        /// False, for the rest of the domain, thus the flux to the adjacent cells wil be evaluated
        /// </param>
        /// <returns></returns>
        public void ReInitSingleStep(out double ChangeRate,SubGrid Restriction = null, bool IncludingInterface = true) {
            if (!IncludingInterface) { throw new NotImplementedException("Untested, not yet functional!"); }

            using (new FuncTrace()) {
                /// Init Residuals
                Residual.Clear();
                Residual.Acc(1.0, Phi);
                OldPhi.Clear();
                OldPhi.Acc(1.0, Phi);
                NewPhi.Clear();
                NewPhi.Acc(1.0, Phi);

                CellMask RestrictionMask = Restriction == null ? null : Restriction.VolumeMask;


                //if (Control.Upwinding && UpdateDirection && IterationCounter % 10 == 0) {

                if (false && Control.Upwinding && UpdateDirection) {
                //if (Control.Upwinding && UpdateDirection) {
                    UpdateBulkMatrix(Restriction);
                }
                
                UpdateDirection = false;
                // RHS part
                RHSField.CoordinateVector.Clear();
                //Operator_RHS.Evaluate(NewPhi.Mapping, RHSField.Mapping);
                Operator_RHS.Evaluate(double.NaN, IncludingInterface ? Restriction : null, IncludingInterface ? SpatialOperator.SubGridBoundaryModes.BoundaryEdge : SpatialOperator.SubGridBoundaryModes.InnerEdge , ArrayTools.Cat(new DGField[] { Phi }, parameterFields, new DGField[] { RHSField }));
#if DEBUG
            RHSField.CheckForNanOrInf();
#endif
                // solve
                // =====
                double[] RHS = OpAffine.CloneAs();
                RHS.ScaleV(-1.0);
                RHS.AccV(1.0, RHSField.CoordinateVector);

                
                SolverResult Result;
                if (Restriction != null) {
                    SubRHS.Clear();
                    SubSolution.Clear();

                    SubRHS.AccV(1.0, RHS, default(int[]), SubVecIdx);
                    SubSolution.AccV(1.0, NewPhi.CoordinateVector, default(int[]), SubVecIdx);

                    Result = slv.Solve(SubSolution, SubRHS);

                    NewPhi.Clear(RestrictionMask);
                    NewPhi.CoordinateVector.AccV(1.0, SubSolution, SubVecIdx, default(int[]) );

                }
                else {
                    Result = slv.Solve(NewPhi.CoordinateVector, RHS);
                }
#if Debug
            OpMatrix.SpMV(-1.0, NewPhi.CoordinateVector, 1.0, RHS);
            Console.WriteLine("LinearSolver: {0} Iterations, Converged={1}, Residual = {2} ", Result.NoOfIterations, Result.Converged, RHS.L2Norm());
#endif





                /// Apply underrelaxation

                Phi.Clear(RestrictionMask);
                Phi.Acc(1 - underrelaxation, OldPhi, RestrictionMask);
                Phi.Acc(underrelaxation, NewPhi, RestrictionMask);
                Residual.Acc(-1.0, Phi, RestrictionMask);
                ChangeRate = Residual.L2Norm(RestrictionMask);


                //Calculate
                LevelSetGradient.Clear();
                LevelSetGradient.Gradient(1.0, Phi, RestrictionMask);
                //LevelSetGradient.GradientByFlux(1.0, Phi);
                MeanLevelSetGradient.Clear();
                MeanLevelSetGradient.AccLaidBack(1.0, LevelSetGradient, RestrictionMask);

                if (Control.Upwinding) {    
                    //RestrictionMask.GetBitMask();
                    for (int i = 0; i < MeanLevelSetGradient.CoordinateVector.Length; i++) {
                        NewDirection[i] = Math.Sign(MeanLevelSetGradient.CoordinateVector[i]);
                        //NewDirection[i] = MeanLevelSetGradient.CoordinateVector[i];
                        OldDirection[i] -= NewDirection[i];
                    }

                    double MaxDiff = OldDirection.L2Norm();
                    //if (MaxDiff > 1E-20 && IterationCounter % 10 == 0 ) {
                    //if (MaxDiff > 1E-20) {
                    //    Console.WriteLine("Direction Values differ by {0}", MaxDiff);
                    if (MaxDiff > 0.2 ) {
                        //UpdateDirection = true;
                        //Console.WriteLine("Direction Values differ by {0} => Updating ReInit-Matrix", MaxDiff);
                    };
                    //}

                    //Console.WriteLine("HACK!!! Updating Upwind Matrix everytime!");
                    //UpdateDirection = true;

                    // Reset Value
                    OldDirection.Clear();
                    OldDirection.AccV(1.0, NewDirection);
                }

            }
        }


        private void UpdateBulkMatrix(SubGrid Restriction = null) {
            //if (Restriction != null) throw new NotImplementedException();

            if (Control.Upwinding) {
                LevelSetGradient.Clear();
                LevelSetGradient.Gradient(1.0, Phi, Restriction?.VolumeMask);
                //LevelSetGradient.GradientByFlux(1.0, Phi);
                MeanLevelSetGradient.Clear();
                MeanLevelSetGradient.AccLaidBack(1.0, LevelSetGradient, Restriction?.VolumeMask);
            }
            /// Matrix and RHS for the Bulk component
            Console.WriteLine("DirectionMatrix has changed. Recomputing LHS");
                OpMatrix_bulk.Clear();
                OpAffine_bulk.Clear();
                // bulk part of the matrix
                Operator_bulk.ComputeMatrix(
                        Phi.Mapping,
                        parameterFields,
                        Phi.Mapping,
                        OpMatrix_bulk, OpAffine_bulk,
                        OnlyAffine: false, sgrd: Restriction);

                /// Compose the Matrix
                /// This is symmetric due to the symmetry of the SIP and the penalty term
                OpMatrix.Clear();
                OpMatrix.Acc(1.0, OpMatrix_bulk);
                OpMatrix.Acc(1.0, OpMatrix_interface);
                OpMatrix.AssumeSymmetric = false;

                /// Compose the RHS of the above operators. (-> Boundary Conditions)
                /// This does NOT include the Nonlinear RHS, which will be added later
                OpAffine.Clear();
                OpAffine.AccV(1.0, OpAffine_bulk);
                OpAffine.AccV(1.0, OpAffine_interface);
        }


        /// <summary>
        /// Based on the Ideas by 
        /// C. Basting and D. Kuzmin, 
        /// “A minimization-based finite element formulation for interface-preserving level set reinitialization”,
        /// Computing, vol. 95, no. 1, pp. 13–25, Dec. 2012.
        /// Create Spatial Operators and build the corresponding Matrices
        /// For the Left-Hand Side of the ReInitProblem
        /// RHS is computed on the fly in <see cref="ReInitSingleStep"/>
        /// The Bulk component is constant unless the grid changes, thus it is computed in <see cref="BuildOperators(CellQuadratureScheme)"/>.
        /// The Interface component changes with its motion.
        /// This component is calculated in <see cref="UpdateOperators(CellQuadratureScheme)"/>.
        /// </summary>
        /// <param name="LSTrck"></param>
        /// <param name="Control">various parameters <paramref name="EllipticReinitControl"/></param>
        /// <param name="HMFOrder">order of tghe interface quadrature</param>
        public EllipticReInit(LevelSetTracker LSTrck, EllipticReInitAlgoControl Control, SinglePhaseField LevelSetForReInit = null) {
            this.Control = Control;
            this.LevelSetTracker = LSTrck;
            if (LevelSetForReInit == null) {
                Phi =  LevelSetTracker.LevelSets[0] as SinglePhaseField;
            }
            else {
                Phi = LevelSetForReInit;
            }
            this.underrelaxation = Control.underrelaxation;

            Residual = new SinglePhaseField(Phi.Basis);
            OldPhi = new SinglePhaseField(Phi.Basis);
            NewPhi = new SinglePhaseField(Phi.Basis);
            foreach (SinglePhaseField f in new List<SinglePhaseField> { Residual, OldPhi, NewPhi }) {
                f.Clear();
                f.Acc(1.0, Phi);
            }


            this.D = LevelSetTracker.GridDat.SpatialDimension;

            this.ConvergenceCriterion = Control.ConvergenceCriterion;
            this.MaxIteration = Control.MaxIt;

            double PenaltyBase = ((double)((Phi.Basis.Degree + 1) * (Phi.Basis.Degree + D))) / ((double)D);


            /// Choose Forms according to Upwinding or Central Fluxes
            string[] paramNames;
            int noOfParamFields;
            
            IEquationComponent BulkForm;
            RHSForm myRHSForm;
            LevelSetGradient = new VectorField<SinglePhaseField>(D, Phi.Basis, "LevelSetGradient", SinglePhaseField.Factory);
            MeanLevelSetGradient = new VectorField<SinglePhaseField>(D, new Basis(Phi.GridDat, 0), "MeanLevelSetGradient", SinglePhaseField.Factory);

            if (Control.Upwinding) {
                paramNames = new string[] {"OldLevelSet", "MeanLevelSetGradient[0]", "MeanLevelSetGradient[1]" };
                noOfParamFields = D;
                LevelSetGradient.Clear();
                LevelSetGradient.Gradient(1.0, Phi);
                //LevelSetGradient.GradientByFlux(1.0, Phi);
                MeanLevelSetGradient.Clear();
                MeanLevelSetGradient.AccLaidBack(1.0, LevelSetGradient);

                parameterFields = ArrayTools.Cat(new SinglePhaseField[] { OldPhi }, MeanLevelSetGradient.ToArray());
                //throw new NotImplementedException("ToDO");
                BulkForm = new EllipticReInitUpwindForm_Laplace(Control.PenaltyMultiplierFlux*PenaltyBase, LSTrck);
                myRHSForm = new EllipticReInitUpwindForm_RHS(Control.PenaltyMultiplierFlux*PenaltyBase, LSTrck);

                OldDirection = new double[MeanLevelSetGradient.CoordinateVector.ToArray().Length];
                for (int i = 0; i < MeanLevelSetGradient.CoordinateVector.Length; i++) {
                    OldDirection[i] = Math.Sign(MeanLevelSetGradient.CoordinateVector[i]);
                }
                NewDirection = OldDirection.CloneAs();

            }
            else
            {
                paramNames = new string[] { };
                noOfParamFields = 0;
                parameterFields = new SinglePhaseField[] { };
                BulkForm = new CentralDifferencesLHSForm(Control.PenaltyMultiplierFlux*PenaltyBase, LSTrck.GridDat.Cells.cj);
                myRHSForm = new CentralDifferencesRHSForm(Control.PenaltyMultiplierFlux*PenaltyBase, LSTrck);
            }


            // SIP for the bulk Phase
            //this.Operator_bulk = new SpatialOperator(1, noOfParamFields, 1, QuadOrderFunc.SumOfMaxDegrees(1, RoundUp: false), variableNames);
            this.Operator_bulk = BulkForm.Operator();
            


            /// Zero at the Interface
            /// 
            /// Calculate Quadrature Order
            Func<int[], int[], int[], int> InterfaceQuadOrder;
            InterfaceQuadOrder = QuadOrderFunc.FixedOrder(Phi.Basis.Degree * 2 + 2);

            /// Generate Interface Operator
            this.Operator_interface = (new EllipticReInitInterfaceForm(Control.PenaltyMultiplierInterface*PenaltyBase, LSTrck)).XOperator(InterfaceQuadOrder);
            
            // Nonlinear Part on the RHS
            /// switch for the potential functions
            switch (Control.Potential) {
                case ReInitPotential.BastingDoubleWell: {
                        myRHSForm.DiffusionRate = ((d, b) => DiffusionRates.DoubleWell(d,b));
                        break;
                    };
                case ReInitPotential.BastingSingleWell: {
                        myRHSForm.DiffusionRate = ((d, b) => DiffusionRates.SingleWell(d,b));
                        break;
                    };
                case ReInitPotential.SingleWellNear: {
                        myRHSForm.DiffusionRate = ((d, b) => DiffusionRates.SingleWellNear(d, b));
                        break;
                    };
                case ReInitPotential.P4DoubleWell: {
                        Console.WriteLine("Warning - This Option for Elliptic ReInit does not work well");
                        myRHSForm.DiffusionRate = ((d, b) => DiffusionRates.DoubleWellAlternative(d, b));
                        break;
                    };
                case ReInitPotential.SingleWellOnCutDoubleWellElse: {
                        myRHSForm.DiffusionRate = ((d, b) => DiffusionRates.SingleWellOnCutDoubleWellElse(d, b));
                        break;
                    }
            }
            Operator_RHS = myRHSForm.Operator(QuadOrderFunc.SumOfMaxDegrees(2, RoundUp: true));


            /// The result of the nonlinear part on the rhs is projected on a single-phase field
            RHSField = new SinglePhaseField(Phi.Basis,"RHS");

            OpMatrix = new MsrMatrix(this.Phi.Mapping, this.Phi.Mapping);
            OpAffine = new double[OpMatrix.RowPartitioning.LocalLength];

            /// Matrix and RHS for the Bulk component
            OpMatrix_bulk = new MsrMatrix(this.Phi.Mapping, this.Phi.Mapping);
            OpAffine_bulk = new double[OpMatrix.RowPartitioning.LocalLength];

            /// Matrix and RHS for the Interface Penalty
            OpMatrix_interface = new MsrMatrix(this.Phi.Mapping, this.Phi.Mapping);
            OpAffine_interface = new double[OpMatrix.RowPartitioning.LocalLength];

            // Init Parameter Fields
            OldPhi.Clear();
            OldPhi.Acc(1.0, Phi);

            // Compute Matrices
            UpdateBulkMatrix();            
        }


    }
}
