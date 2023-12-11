using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using BoSSS.Solution.Utils;
using SAIDT.Fluxes;
using ApplicationWithIDT;
using BoSSS.Solution.AdvancedSolvers;

namespace SAIDT {
    /// <summary>
    /// Implements XDG space-time Scalar Advection in (1D in space) which is solved by the routines defined in <see cref="ApplicationWithIDT"/>
    /// Naming: S(calar) - A(dvection) - I(mplict) - D(iscontinuity) - T(racking)
    /// Concrete configurations of solver (Initial Guess, optimization parameters,...) are set in a <see cref="SAIDTControl.cs"/> object, e.g. boundary conditions are set by  by the property <see cref="IDTControl.DirichletBoundaryMap"/>
    /// Fluxes are implemented in <see cref="SAIDT.Fluxes"/>, so far only upwind flux is supported
    /// 
    /// Author: Jakob Vandergrift 
    /// Date of Creation/Maintenance: 08-2022 until at least 08-2024
    /// </summary>
    public class SAIDTMain : ApplicationWithIDT<SAIDTControl> {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="args">string pointing to a control file, i.e. `cs:SAIDT.SAIDTHardCodedControl.CurvedShock_Eccomas22()` </param>
        static void Main(string[] args) {
            //SAIDT.Tests.SAIDTTestProgram.StraightShock_p0_SInglePhaseFieldLS();
            //SAIDT.Tests.SAIDTTestProgram.CurvedShock_Eccomas22();
            SAIDTMain._Main(args, false, () => new SAIDTMain());
        }

        /// <summary>
        /// The state c
        /// </summary>
        public XDGField Concentration { get; private set; }

        /// <summary>
        /// Creates the XDG Fields for the conserved quantities, here the Concentration c
        /// </summary>
        /// <param name="in_LsTrk"></param>
        /// <param name="in_DgDegree"></param>
        public override void CreateConservativeFields(LevelSetTracker in_LsTrk, int in_DgDegree) {
            #region Create mandatory conservative fields
            int D = in_LsTrk.GridDat.SpatialDimension;
            this.Concentration = new XDGField(new XDGBasis(in_LsTrk, in_DgDegree), "c");
            Concentration.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            ConservativeFields = new XDGField[1];
            ConservativeFields[0] = Concentration;
            #endregion
        }
        /// <summary>
        /// Loads the Grid
        /// </summary>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        protected override IGrid CreateOrLoadGrid() {
            IGrid grid = null;
            grid = base.CreateOrLoadGrid();
            CompressibleEnvironment.Initialize(grid.SpatialDimension);
            return grid;
        }
        /// <summary>
        /// Initializes the optimization problem and computes and initial residual, also if chosen initial value is computed by P0 projection
        /// 1. Here the objects for the discretized equation are constructed by constructing a Spatial Operator and adding the fluxes defined in <see cref="SAIDT.Fluxes"/>
        /// 2. The functional needed for the objective function is constructed
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            #region Init operator
            

            GridData gridData = (GridData)this.GridData;
            this.XSpatialOperator = new XDifferentialOperatorMk2(new string[] { "c" }, null, new string[] { "c" }, Control.quadOrderFunc, this.SpeciesToEvaluate);
            this.Op_obj = new XDifferentialOperatorMk2(new string[] { "c" }, null, new string[] { "c" }, Control.quadOrderFunc, this.SpeciesToEvaluate);
            #endregion

            #region add EquationComponents
            this.XSpatialOperator.EquationComponents["c"].Add(new SAIDT.Fluxes.ScalarAdvectionUpwindFlux("L", Control.ShockPos, Control.LeftValue, Control.RightValue, Control.FlowFunc, Control.DirichletBoundaryMap));
            this.XSpatialOperator.EquationComponents["c"].Add(new SAIDT.Fluxes.ScalarAdvectionUpwindFlux("R", Control.ShockPos, Control.LeftValue, Control.RightValue, Control.FlowFunc, Control.DirichletBoundaryMap));
            XSpatialOperator.EquationComponents["c"].Add(new ScalarAdvectionUpwindFlux_Interface(LsTrk, Control.FlowFunc));
            XSpatialOperator.LinearizationHint = LinearizationHint.FDJacobi;
            #endregion

            #region construct Objective function Operator
            switch(Control.optProblemType) {
                case OptProblemType.RankineHugoniotFull:
                //Operator for Objective
                this.Op_obj.EquationComponents["c"].Add(new SAIDT.Fluxes.ScalarAdvectionRHFlux("L", Control.ShockPos, Control.LeftValue, Control.RightValue, Control.FlowFunc, Control.DirichletBoundaryMap));
                this.Op_obj.EquationComponents["c"].Add(new SAIDT.Fluxes.ScalarAdvectionRHFlux("R", Control.ShockPos, Control.LeftValue, Control.RightValue, Control.FlowFunc, Control.DirichletBoundaryMap));
                this.Op_obj.EquationComponents["c"].Add(new ScalarAdvectionRHFlux_Interface(LsTrk, Control.FlowFunc));
                break;
                case OptProblemType.RankineHugoniotOnlyInterface:
                //Operator for Objective
                this.Op_obj.EquationComponents["c"].Add(new ScalarAdvectionRHFlux_Interface(LsTrk, Control.FlowFunc));
                break;
                default:
                Op_obj = XSpatialOperator.CloneAs();
                break;
            }
            switch(Control.Linearization) {
                case Linearization.FD:
                XSpatialOperator.LinearizationHint = LinearizationHint.FDJacobi;
                Op_obj.LinearizationHint = LinearizationHint.FDJacobi;
                break;
                case Linearization.JacobiOperator:
                XSpatialOperator.LinearizationHint = LinearizationHint.GetJacobiOperator;
                Op_obj.LinearizationHint = LinearizationHint.GetJacobiOperator;
                break;
                case Linearization.Adhoc:
                XSpatialOperator.LinearizationHint = LinearizationHint.AdHoc;
                Op_obj.LinearizationHint = LinearizationHint.AdHoc;
                break;
            }
            //init OptProblem
            ChooseOptProblem();
            #endregion

            #region commit Operators
            Op_obj.Commit();
            XSpatialOperator.Commit();
            #endregion 

            #region necessary initializations for Operator evaluation
            LsTBO = LevelSet;
            LevelSetOpti.AssembleTransMat(LsTBO);
            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            LsTrk.UpdateTracker(CurrentStepNo);
            LsTrk.PushStacks();
            //note that the operator is assembled we can compute the p0 solution
            if(Control.GetInitialValue != GetInitialValue.FromFunctionPerSpecies) {
                ComputeP0Solution();
            }
            //Initialize empty vectors and matrices
            InitializeMatricesAndVectors();
            //// Cell agglomeration 
            UpdateAgglomerator();
            #endregion

            #region  Compute Residual and Derived Quantities
            ComputeResiduals();
            InitResNorm = res_l2;
            Init_obj_f = obj_f;
            ResNorms.Add(res_l2);
            obj_f_vals.Add(obj_f);
            UpdateDerivedVariables();
            #endregion

        }
        public Type GetSolverType() {
            return typeof(SAIDTMain);
        }
        /// <summary>
        /// Solves the linear P=0 problem, by assembling the operator matrix of the residual and solving the linear system. 
        /// </summary>
        public void ComputeP0Solution() {
            this.Concentration.Clear();
            XDGBasis basis_p0 = new XDGBasis(LsTrk, 0);
            XDGField c_p0 = new XDGField(basis_p0, "c_p0_initial");
            XDGField[] ConservativeVarsP0 = new XDGField[1];
            ConservativeVarsP0[0] = c_p0;
            var mappingP0 = new CoordinateMapping(ConservativeVarsP0);
            var Matbuilder = XSpatialOperator.GetMatrixBuilder(mappingP0, null, mappingP0);
            var Mat = new BlockMsrMatrix(mappingP0);
            double[] offset = new double[mappingP0.TotalLength];
            Matbuilder.ComputeMatrix(Mat, offset);
            double[] nullVec = new double[mappingP0.TotalLength];
            double[] solVec = new double[mappingP0.TotalLength];
            offset.ScaleV(-1.0);
            SimpleSolversInterface.Solve_Direct(Mat, solVec, offset);

            c_p0.CoordinateVector.SetV(solVec, 1.0);
            //c_p0.ProjectField((x,t) => ExactSol(x,t,FlowFunc));
            ConservativeFields[0].AccLaidBack(1.0, c_p0);
        }
        /// <summary>
        /// Sets the Initial Value for State and Level Set
        /// </summary>
        /// <param name="t"></param>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        protected override void SetInitial(double t) {
            switch(Control.OptiLevelSetType) {
                case OptiLevelSetType.SinglePhaseField:
                case OptiLevelSetType.SpecFemField:
                case OptiLevelSetType.SplineLevelSet:
                LevelSetOpti.ProjectFromFunction(Control.LevelSetTwoInitialValue);
                break;
                case OptiLevelSetType.GlobalLevelSet:
                break;

                default: throw new ArgumentOutOfRangeException(nameof(Control.OptiLevelSetType));
            }
            //We project the LevelSetOpti object onto the DG LsTBO
            LevelSetOpti.AssembleTransMat(LsTBO);
            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            LsTrk.UpdateTracker(CurrentStepNo);
            LsTrk.PushStacks();
            switch (Control.GetInitialValue) {
                case GetInitialValue.FromFunctionPerSpecies:
                foreach(string spec in this.SpeciesToEvaluate) {
                    Concentration.GetSpeciesShadowField(spec).ProjectField(Control.GetInitialValueFunc(spec));
                }

                break;
            }
            

            //The initial values are set elsewhere as we use the p0 projection and thus we need to assemble the operator first

        }
        /// <summary>
        /// Updates derived quantities, but does nothing here as there are no derived variables defined so far for scalar advection
        /// </summary>
        public override void UpdateDerivedVariables() {
        }
        /// <summary>
        /// Some BoSSS structure needed for agglomeration
        /// </summary>
        public override void InitializeMultiGridOpConfig() {
            MultiGridOperatorConfig = new MultigridOperator.ChangeOfBasisConfig[][]{
                    new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] { 0 }, mode = MultigridOperator.Mode.Eye, DegreeS = new int[] { Concentration.Basis.Degree}
                        }
                    }
            };
        }
    }
}

