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
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.LinSolvers.MUMPS;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.EllipticExtension {


    /// <summary>
    /// Elliptic Extension Velocity
    /// Using a similar idea as used for the Elliptic ReInit for the Extension velocity.
    /// 
    /// See
    /// 
    /// Utz,Kummer,Oberlack - A high-order DG method for the extension problem (2017)
    /// 
    /// 
    /// Based on Ideas from 
    /// 
    /// C. Basting and D. Kuzmin, 
    /// “A minimization-based finite element formulation for interface-preserving level set reinitialization”,
    /// Computing, vol. 95, no. 1, pp. 13–25, Dec. 2012.
    /// And
    /// Utz, Kummer, Oberlack - Interface-preserving level-set reinitialization for DG-FEM (2016)
    /// 
    /// 
    /// </summary>
    public class Extender {
        int D;
        private LevelSetTracker LevelSetTracker;
        SpatialOperator Operator_bulk;
        XSpatialOperatorMk2 Operator_interface;
        LevelSet Phi;
        public SinglePhaseField MeanLevelSet;

        SinglePhaseField Extension;
        //SinglePhaseField InterfaceValue;
        VectorField<SinglePhaseField> LevelSetGradient;
        public VectorField<SinglePhaseField> MeanLevelSetGradient;

        MsrMatrix OpMatrix;
        double[] OpAffine;

        MsrMatrix OpMatrix_bulk;
        double[] OpAffine_bulk;
        IList<DGField> BulkParams;

        MsrMatrix OpMatrix_interface;
        double[] OpAffine_interface;
        List<DGField> InterfaceParams;


        EllipticExtVelAlgoControl Control;





        /// <summary>
        /// Calculate Extension
        /// </summary>
        public void ConstructExtension(IList<DGField> InterfaceValue=null, bool nearfield = false) {
            using (new FuncTrace()) {
                Extension.Clear(nearfield ? LevelSetTracker.Regions.GetNearFieldMask(1) : null);
                ComputeMatrices(InterfaceValue ?? InterfaceParams,nearfield);

                //Solve System
                double[] RHS = OpAffine.CloneAs();
                RHS.ScaleV(-1.0);
                
                ISparseSolver slv = Control.solverFactory();

                if (nearfield) {
                    SubGrid subGrid = LevelSetTracker.Regions.GetNearFieldSubgrid(1);
                    int[] SubVecIdx = Extension.Mapping.GetSubvectorIndices(subGrid, true, new int[] { 0 });
                    int L = SubVecIdx.Length;
                    MsrMatrix SubMatrix = new MsrMatrix(L);
                    double[] SubRHS = new double[L];
                    double[] SubSolution = new double[L];

                    OpMatrix.AccSubMatrixTo(1.0, SubMatrix, SubVecIdx, default(int[]), SubVecIdx, default(int[]));


                    SubRHS.Clear();
                    SubSolution.Clear();

                    SubRHS.AccV(1.0, RHS, default(int[]), SubVecIdx);
                    SubSolution.AccV(1.0, Extension.CoordinateVector, default(int[]), SubVecIdx);

                    slv.DefineMatrix(SubMatrix);
                    slv.Solve(SubSolution, SubRHS);

                    Extension.Clear(subGrid.VolumeMask);
                    Extension.CoordinateVector.AccV(1.0, SubSolution, SubVecIdx, default(int[]));

                }
                else {
                    slv.DefineMatrix(OpMatrix);
                    slv.Solve(Extension.CoordinateVector, RHS);
                }
                slv.Dispose();
            }
        }


        /// <summary>
        /// Calculate Extension
        /// </summary>
        public void ConstructExtension(VectorField<SinglePhaseField> NewExtension, VectorField<SinglePhaseField> InterfaceValue,ref SubGrid subGrid)
        {
            using (new FuncTrace()) {
                NewExtension.Clear();
                ComputeMatrices(new List<DGField> { InterfaceValue[0] }, Control.subGridRestriction);
                ISparseSolver slv = Control.solverFactory();

                int[] SubVecIdx = null;
                int L;
                MsrMatrix SubMatrix;
                double[] SubRHS = null;
                double[] SubSolution = null;

                if (subGrid != null) {
                    SubVecIdx = Extension.Mapping.GetSubvectorIndices(subGrid, true, new int[] { 0 });
                    L = SubVecIdx.Length;
                    SubMatrix = new MsrMatrix(L);
                    SubRHS = new double[L];
                    SubSolution = new double[L];

                    OpMatrix.AccSubMatrixTo(1.0, SubMatrix, SubVecIdx, default(int[]), SubVecIdx, default(int[]));
                    slv.DefineMatrix(SubMatrix);
                }
                else {
                    slv.DefineMatrix(OpMatrix);
                }

                for (int d = 0; d < D; d++) {
                    UpdateRHS(NewExtension[d], InterfaceValue[d], Control.subGridRestriction);
                    double[] RHS = OpAffine.CloneAs();
                    RHS.ScaleV(-1.0);
                    if (subGrid != null) {
                        SubRHS.Clear();
                        SubSolution.Clear();
                        SubRHS.AccV(1.0, RHS, default(int[]), SubVecIdx);
                        SubSolution.AccV(1.0, NewExtension[d].CoordinateVector, default(int[]), SubVecIdx);
                        slv.Solve(SubSolution, SubRHS);

                        NewExtension[d].Clear(subGrid.VolumeMask);
                        NewExtension[d].CoordinateVector.AccV(1.0, SubSolution, SubVecIdx, default(int[]));
                    }
                    else {
                        slv.Solve(NewExtension[d].CoordinateVector, RHS);
                    }
                }
                slv.Dispose();
            }
                
                
        }
        


        /// <summary>
        /// Create Spatial Operators and build the corresponding Matrices
        /// </summary>
        public void ComputeMatrices(IList<DGField> InterfaceParams,bool nearfield) {
            OpMatrix = new MsrMatrix(this.Extension.Mapping, this.Extension.Mapping);
            OpAffine = new double[OpMatrix.RowPartitioning.LocalLength];

            OpMatrix_bulk = new MsrMatrix(this.Extension.Mapping, this.Extension.Mapping);
            OpAffine_bulk = new double[OpMatrix.RowPartitioning.LocalLength];

            OpMatrix_interface = new MsrMatrix(this.Extension.Mapping, this.Extension.Mapping);
            OpAffine_interface = new double[OpMatrix.RowPartitioning.LocalLength];


            //LevelSetTracker.GetLevelSetGradients(0,);

            // bulk part of the matrix
            //Operator_bulk.ComputeMatrix(
            //    Extension.Mapping,
            //    LevelSetGradient.ToArray(),
            //    Extension.Mapping,
            //    OpMatrix_bulk, OpAffine_bulk,
            //    OnlyAffine: false, sgrd: null);

            switch (Control.FluxVariant) {
                case FluxVariant.GradientBased:
                    // Flux Direction based on Mean Level Set Gradient
                    BulkParams = new List<DGField> { }; // Hack, to make ArrayTools.Cat produce a List of DGFields
                    // second Hack: Does only work, when InterfaceParams is according to a single component flux,
                    // else, we will have to change the boundary edge flux
                    BulkParams = ArrayTools.Cat(BulkParams,LevelSetGradient.ToArray(), Phi, MeanLevelSetGradient.ToArray(),InterfaceParams.ToArray());
                    MeanLevelSetGradient.Clear();
                    MeanLevelSetGradient.AccLaidBack(1.0, LevelSetGradient);
                    break;
                case FluxVariant.ValueBased:
                    // Flux Direction Based on Cell-Averaged Level-Set Value
                    BulkParams = ArrayTools.Cat(LevelSetGradient.ToArray(), Phi, MeanLevelSet);
                    MeanLevelSet.Clear();
                    MeanLevelSet.AccLaidBack(1.0, Phi);
                    break;
                case FluxVariant.SWIP:
                    BulkParams = LevelSetGradient.ToArray();
                    break;
                default:
                    throw new Exception();
            }

            // Build Operator

            Operator_bulk.ComputeMatrixEx(Extension.Mapping,
                BulkParams,
                Extension.Mapping,
                OpMatrix_bulk, OpAffine_bulk,
                OnlyAffine: false,
                time: 0.0,
                edgeQuadScheme: new EdgeQuadratureScheme(true, nearfield ? LevelSetTracker.Regions.GetNearFieldSubgrid(1).InnerEdgesMask : null),
                volQuadScheme: new CellQuadratureScheme(true, nearfield ? LevelSetTracker.Regions.GetNearFieldSubgrid(1).VolumeMask : null)
                );

            

            //Operator_interface.ComputeMatrixEx(
            //    LevelSetTracker,
            //    Extension.Mapping,
            //    InterfaceParams,
            //    Extension.Mapping,
            //    OpMatrix_interface,
            //    OpAffine_interface,
            //    OnlyAffine: false,
            //    time: 0,
            //    MPIParameterExchange: false,
            //    whichSpc: LevelSetTracker.GetSpeciesId("A")
            //    );
            XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = Operator_interface.GetMatrixBuilder(LevelSetTracker,
                Extension.Mapping, InterfaceParams, Extension.Mapping, LevelSetTracker.GetSpeciesId("A"));

            MultiphaseCellAgglomerator dummy = LevelSetTracker.GetAgglomerator(LevelSetTracker.SpeciesIdS.ToArray(), 2 * Extension.Basis.Degree + 2, 0.0);
            mtxBuilder.SpeciesOperatorCoefficients[LevelSetTracker.GetSpeciesId("A")].CellLengthScales = dummy.CellLengthScales[LevelSetTracker.GetSpeciesId("A")];

            mtxBuilder.time = 0;
            mtxBuilder.MPITtransceive = false;
            mtxBuilder.ComputeMatrix(OpMatrix_interface, OpAffine_interface);

#if DEBUG
            OpMatrix_bulk.CheckForNanOrInfM();
            OpAffine_bulk.CheckForNanOrInfV();

            OpMatrix_interface.CheckForNanOrInfM();
            OpAffine_interface.CheckForNanOrInfV();
#endif

            //Only for Debugging purposes

            Debug.Assert(OpMatrix_interface.GetDiagVector().L2Norm() > 0, "L2-Norm of Diagonal of InterfaceOperator is 0");
            Debug.Assert(OpMatrix_bulk.GetDiagVector().L2Norm() > 0, "L2-Norm of Diagonal of BulkOperator is 0");
#if DEBUG
            //Console.WriteLine( "L2-Norm of Diagonal of InterfaceOperator is {0}", OpMatrix_interface.GetDiagVector().L2Norm() );
#endif
            OpMatrix.Clear();
            OpMatrix.Acc(1.0, OpMatrix_bulk);
            OpMatrix.Acc(1.0, OpMatrix_interface);
            //Console.WriteLine("Op-Matrix Symmetry-Deviation: {0}", OpMatrix.SymmetryDeviation());
            OpMatrix.AssumeSymmetric = false;

            OpAffine.Clear();
            OpAffine.AccV(1.0, OpAffine_bulk);
            OpAffine.AccV(1.0, OpAffine_interface);
#if DEBUG
            //Console.WriteLine("Condition Number of Extension Operator {0}", OpMatrix.condest());
#endif
        }


        /// <summary>
        /// Create Spatial Operators and build the corresponding Matrices
        /// </summary>
        public void UpdateRHS(SinglePhaseField Extension, SinglePhaseField InterfaceValue,bool nearfield) {

            this.Extension = Extension;

            OpAffine.Clear();

            //XSpatialOperatorExtensions.ComputeMatrixEx(Operator_interface,
            ////Operator_interface.ComputeMatrixEx(
            //   LevelSetTracker,
            //    Extension.Mapping,
            //    new List<DGField> { InterfaceValue },
            //    Extension.Mapping,
            //    OpMatrix_interface,
            //    OpAffine_interface,
            //    OnlyAffine: true,
            //    time: 0,
            //    MPIParameterExchange: false,
            //    whichSpc: LevelSetTracker.GetSpeciesId("A"),
            //    subGrid: nearfield ? LevelSetTracker.Regions.GetNearFieldSubgrid(1) : null
            //    );
            XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = Operator_interface.GetMatrixBuilder(LevelSetTracker, 
                Extension.Mapping, new List<DGField> { InterfaceValue }, Extension.Mapping, LevelSetTracker.GetSpeciesId("A"));
            mtxBuilder.time = 0;
            mtxBuilder.MPITtransceive = false;
            mtxBuilder.ComputeAffine(OpAffine_interface);

            if (OpAffine.L2Norm() == 0) Console.WriteLine("RHS of Bulk equation is empty as expected.");

            OpAffine.Clear();
            OpAffine.AccV(1.0, OpAffine_bulk);
            OpAffine.AccV(1.0, OpAffine_interface);

        }


        /// <summary>
        /// 
        /// </summary>
        public void AnalyzeOperators(out MultidimensionalArray ret) {
            int length = OpMatrix.NoOfRows / Extension.Basis.Length;
            //int[,] PatternMatrix = new int[length, length]; // lässt sich besser anschauen?
            MsrMatrix PatternMatrix = new MsrMatrix(length);
            for (int i=0; i<length; i++) {
                for (int j = 0; j < length; j++) {
                    if (OpMatrix[i*Extension.Basis.Length, j*Extension.Basis.Length] != 0) {
                        PatternMatrix[i, j] = 1;
                    }
                }

                }



            ret = MultidimensionalArray.Create(1, 4);
            Console.WriteLine("Calling MATLAB/Octave...");
            MultidimensionalArray L = MultidimensionalArray.Create(PatternMatrix.NoOfRows, PatternMatrix.NoOfRows);
            MultidimensionalArray U = MultidimensionalArray.Create(PatternMatrix.NoOfRows, PatternMatrix.NoOfRows);
            MultidimensionalArray P = MultidimensionalArray.Create(PatternMatrix.NoOfRows, PatternMatrix.NoOfRows);
            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(PatternMatrix, "PatternMatrix");
                bmc.Cmd("[L, U, P] = lu(PatternMatrix);");


                bmc.PutSparseMatrix(OpMatrix, "OpMatrix");
                bmc.Cmd("condNoOpMatrix = condest(OpMatrix);");
                bmc.Cmd("eigiMaxi = eigs(OpMatrix,1,'lm')");
                bmc.Cmd("eigiMini = eigs(OpMatrix,1,'sm')");
                bmc.Cmd("lasterr");
                bmc.Cmd("[V,r]=chol(OpMatrix);");
                bmc.Cmd("ret = [condNoOpMatrix, eigiMaxi, eigiMini, r]");

                bmc.GetMatrix(L, "full(L)");
                bmc.GetMatrix(U, "full(U)");
                bmc.GetMatrix(P, "full(P)");
                bmc.GetMatrix(ret, "ret");

                bmc.Execute(true);
            }
            //PatternMatrix.SaveToTextFile("Patternmatrix.txt");
            //PatternMatrix.SaveToTextFileSparse("PatternmatrixSparse.txt");

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
        /// Based on the Ideas by 
        /// C. Basting and D. Kuzmin, 
        /// “A minimization-based finite element formulation for interface-preserving level set reinitialization”,
        /// Computing, vol. 95, no. 1, pp. 13–25, Dec. 2012.
        /// </summary>
        /// <param name="LSTrck"></param>
        /// <param name="bcmap">Boundary Conditions for the LevelSet Equations</param>
        /// <param name="Control">various parameters <paramref name="EllipticReinitControl"/></param>
        public Extender(SinglePhaseField Extension, LevelSetTracker LSTrck, ILevelSetForm InterfaceFlux, List<DGField> InterfaceParams, VectorField<SinglePhaseField> LevelSetGradient, EllipticExtVelAlgoControl Control) {

            if (InterfaceFlux.ParameterOrdering.Count != InterfaceParams.Count) throw new ArgumentException("Missmatch in Number of Parameters and expected amount in the given flux.");
            this.InterfaceParams = InterfaceParams;

            this.Control = Control;
            int D = LSTrck.GridDat.SpatialDimension;

            //this.InterfaceValue = InterfaceValue;
            //InterfaceValue.Identification = "InterfaceValue";
            this.Extension = Extension;

            this.LevelSetTracker = LSTrck;
            Phi = LevelSetTracker.LevelSets[0] as LevelSet;
            this.LevelSetGradient = LevelSetGradient;

            switch (Control.FluxVariant) {
                case FluxVariant.GradientBased:
                    // Flux Direction based on Mean Level Set Gradient
                    MeanLevelSetGradient = new VectorField<SinglePhaseField>(
                    D.ForLoop(
                            d => new SinglePhaseField(
                                  new Basis(LSTrck.GridDat, 0), VariableNames.MeanLevelSetGradientComponent(d)
                                )
                             )
                    );
                    BulkParams = new List<DGField> { };
                    BulkParams = ArrayTools.Cat(BulkParams,LevelSetGradient.ToArray(), Phi, MeanLevelSetGradient.ToArray(), InterfaceParams.ToArray());
                    break;
                case FluxVariant.ValueBased:
                    // Flux Direction Based on Cell-Averaged Level-Set Value
                    MeanLevelSet = new SinglePhaseField(new Basis(LSTrck.GridDat, 0), "MeanLevelSet");
                    BulkParams = ArrayTools.Cat(LevelSetGradient.ToArray(), Phi, MeanLevelSet);
                    break;
                case FluxVariant.SWIP:
                    BulkParams = LevelSetGradient.ToArray();
                    break;
                default:
                    throw new NotImplementedException();
            }


            this.D = LevelSetTracker.GridDat.SpatialDimension;

            double PenaltyBase = Control.PenaltyMultiplierFlux * ((double)((Phi.Basis.Degree + 1) * (Phi.Basis.Degree + D))) / ((double)D);

            DefineBulkOperator(LSTrck, InterfaceFlux, D, PenaltyBase);

            Operator_interface = InterfaceFlux.XOperator(QuadOrderFunc.FixedOrder(2 * Extension.Basis.Degree + 2) );
        }

        private void DefineBulkOperator(LevelSetTracker LSTrck, ILevelSetForm InterfaceFlux, int D, double PenaltyBase) {

            foreach (Foundation.Grid.RefElements.RefElement r in LevelSetTracker.GridDat.Grid.RefElements) {
                if (r is Foundation.Grid.RefElements.Triangle && this.Control.FluxVariant == FluxVariant.ValueBased) { throw new NotSupportedException("Level-Set Based flux direction does not work on triangular grids "); }
            }
            switch (Control.FluxVariant) {
                case FluxVariant.GradientBased:
                    /// Flux Based on Mean LevelSetGradient
                    this.Operator_bulk = (new EllipticExtVelFormDirected(PenaltyBase, Control.IsotropicViscosity, InterfaceFlux,  LSTrck)).Operator(2);
                    break;
                case FluxVariant.ValueBased:
                    /// Flux Based on Cell-Averaged Level-Set Value
                    this.Operator_bulk = (new EllipticExtVelFormLevelSetBased(PenaltyBase, Control.IsotropicViscosity, LSTrck)).Operator(2);
                    break;
                case FluxVariant.SWIP:
                    this.Operator_bulk = (new EllipticExtVelFormSWIP(PenaltyBase, Control.IsotropicViscosity, LSTrck)).Operator(0);
                    break;
                default:
                    throw new NotImplementedException();
            }
        }


    }




}
