using ApplicationWithIDT;
using ApplicationWithIDT.OptiLevelSets;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Utils;
using BUIDT.Fluxes;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;


namespace BUIDT {
    /// <summary>
    /// Implements XDG space-timeBurgers equation in (1D in space) which is solved by the routines defined in <see cref="ApplicationWithIDT"/>
    /// Naming: BU(rgers) - I(mplict) - D(iscontinuity) - T(racking)
    /// 
    /// Concrete configurations of solver (initial guess, optimization parameters,...) are set in a <see cref="BUIDTControl.cs"/> object, e.g. boundary conditions are set by the property <see cref="IDTControl.DirichletBoundaryMap"/>
    /// Fluxes are implemented in <see cref="BUIDT.Fluxes"/>, so far only upwind flux is supported
    /// 
    /// Author: Jakob Vandergrift 
    /// Date of Creation/Maintenance: 08-2022 until at least 08-2024
    /// </summary>
    public class BUIDTMain : ApplicationWithIDT<BUIDTControl> {
        
        /// <summary>
        /// </summary>
        /// <param name="args">string pointing to a control file, i.e. `cs:BUIDT.BUIDTHardCodedControl.StraightShockCurvedStart_Eccomas22()` </param>
        static void Main(string[] args) {
            //Tests.BUIDTTestProgram.StraightShockCurvedStart_Eccomas22();
            //Tests.BUIDTTestProgram.AcceleratingShock();
            BUIDTMain._Main(args, false, () => new BUIDTMain());
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
        public override void CreateConservativeFields(LevelSetTracker in_LsTrk, int dgDegree) {
            #region Create mandatory conservative fields
            int D = in_LsTrk.GridDat.SpatialDimension;
            this.Concentration = new XDGField(new XDGBasis(in_LsTrk, dgDegree), "c");
            Concentration.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            ConservativeFields = new XDGField[1];
            ConservativeFields[0] = Concentration;
            #endregion

        }
        /// <summary>
        /// creates the grid
        /// </summary>
        /// <returns></returns>
        protected override IGrid CreateOrLoadGrid() {
            IGrid grid = null;
            grid = base.CreateOrLoadGrid();
            CompressibleEnvironment.Initialize(grid.SpatialDimension);

            return grid;
        }
        /// <summary>
        /// Initializes the Optimization Problem and computes and initial residual, also if chosen initial value is computed by P0 projection
        /// 1. Here the objects for the discretized equation are constructed by constructing a Spatial Operator and adding the fluxes defined in <see cref="BUIDT.Fluxes"/>
        /// 2. The Functional needed for the objective function is constructed
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            #region Init operator
            

            GridData gridData = (GridData)this.GridData;
            this.XSpatialOperator = new XDifferentialOperatorMk2(new string[] { "c" }, null, new string[] { "c" }, Control.quadOrderFunc, this.SpeciesToEvaluate);
            this.Op_obj = new XDifferentialOperatorMk2(new string[] { "c" }, null, new string[] { "c" }, Control.quadOrderFunc, this.SpeciesToEvaluate);
            #endregion
            this.XSpatialOperator.EquationComponents["c"].Add(new BUIDT.Fluxes.STBurgersUpwindFlux("L", Control.is_nf_smth, Control.s_alpha, Control.DirichletBoundaryMap));
            this.XSpatialOperator.EquationComponents["c"].Add(new BUIDT.Fluxes.STBurgersUpwindFlux("R", Control.is_nf_smth, Control.s_alpha, Control.DirichletBoundaryMap));
            this.XSpatialOperator.EquationComponents["c"].Add(new BurgersUpwindFlux_Interface(Control.is_nf_smth, Control.s_alpha));

            switch(Control.optProblemType) {
                case OptProblemType.RankineHugoniotFull:
                //Operator for Objective
                this.Op_obj.EquationComponents["c"].Add(new BUIDT.Fluxes.STBurgersRHFlux("L", Control.is_nf_smth, Control.s_alpha, Control.DirichletBoundaryMap));
                this.Op_obj.EquationComponents["c"].Add(new BUIDT.Fluxes.STBurgersRHFlux("R", Control.is_nf_smth, Control.s_alpha, Control.DirichletBoundaryMap));
                this.Op_obj.EquationComponents["c"].Add(new BurgersRHFlux_Interface(Control.is_nf_smth, Control.s_alpha));
                break;
                case OptProblemType.RankineHugoniotOnlyInterface:
                this.Op_obj.EquationComponents["c"].Add(new BurgersRHFlux_Interface(Control.is_nf_smth, Control.s_alpha));
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

            Op_obj.Commit();
            XSpatialOperator.Commit();

            r_JacobiOperator = XSpatialOperator.GetJacobiOperator(SpatialDimension: 2);
            R_JacobiOperator = XSpatialOperator.GetJacobiOperator(SpatialDimension: 2);

            switch (Control.GetInitialValue)
            {
                case GetInitialValue.FromP0Timestepping:
                    ComputeP0Solution();
                    break;
                case GetInitialValue.OneFullNewtonStep:
                    ComputeP0SolutionOneNewtonStep();
                    break;
                default:
                    break;
            }

            //init OptProblem
            ChooseOptProblem();
            //Initialize empty vectors and matrices
            InitializeMatricesAndVectors();
            //// Cell agglomerator (cell length scales are needed for diffusive AV fluxes)
            UpdateAgglomerator();
            (res_l2,obj_f,res_L2)=ComputeResiduals();
            InitResNorm = res_l2;
            Init_obj_f = obj_f;
            ResNorms.Add(res_l2);
            obj_f_vals.Add(obj_f);
            UpdateDerivedVariables();

        }
        public Type GetSolverType() {
            return typeof(BUIDTMain);
        }
        ///// <summary>
        ///// Approximates the linear P=0 problem, doing one Newton step. 
        ///// </summary>
        public void ComputeP0SolutionOneNewtonStep()
        {

            // Clear fields (just to be sure)
            this.Concentration.Clear();

            // Initial conditions are set by computing the p0 solution
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

            //Func<double[], double> exact_sol_func = x => exact_sol_linear(x[0], x[1]);
            //c_p0.ProjectField(exact_sol_func);
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
                switch(Control.GetLevelSet) {
                    case GetLevelSet.FromReconstructionFromPoints:
                    var Points = IMatrixExtensions.LoadFromTextFile(Control.SLSPointPath);
                    if(LevelSetOpti is SplineOptiLevelSet splineLS) {
                        splineLS.GetSplineOverDetermined(Points);
                        SplineOptiLevelSet.EmbeddInLevelSet(splineLS.Spline, splineLS);
                    } else {
                        throw new ArgumentException("something went wrong, as we expect LeveSetOpti to be of type SplineLS");
                    }

                    break;
                    default:
                    LevelSetOpti.ProjectFromFunction(Control.LevelSetTwoInitialValue);
                    break;
                }

                break;
                case OptiLevelSetType.GlobalLevelSet:
                break;

                default: throw new ArgumentOutOfRangeException(nameof(Control.OptiLevelSetType));
            }
            switch(Control.GetInitialValue) {
                case GetInitialValue.FromFunctionPerSpecies:
                XDGField ConcentrationP0 = new XDGField(new XDGBasis(LsTrk, 0));

                foreach(string spec in this.SpeciesToEvaluate) {
                    ConcentrationP0.GetSpeciesShadowField(spec).ProjectField(Control.GetInitialValueFunc(spec));
                    //centration.GetSpeciesShadowField(spec).ProjectField(Control.GetInitialValueFunc(spec));
                }
                Concentration.AccLaidBack(1.0, ConcentrationP0);
                break;
            }
            //We project the LevelSetOpti object onto the DG LsTBO
            LevelSetOpti.AssembleTransMat(LsTBO);
            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            LsTrk.UpdateTracker(CurrentStepNo);
            LsTrk.PushStacks();

        }
        /// <summary>
        /// Updates Derived Quantities, but does nothing here as there are no derived variables defined so far for Scalar Advection
        /// </summary>
        public override void UpdateDerivedVariables() {
        }

        /// <summary>
        /// Some BoSSS structure needed for Agglomeration
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


