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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Platform;
using ilPSP;
using System.Linq;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Solution.Queries;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.Timestepping;
using System.Runtime.CompilerServices;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet {
    /// <summary>
    /// Base class for a level set based on the Cahn-Hilliard Phasefield equation
    /// </summary>
    public partial class Phasefield : DgApplicationWithSolver<PhasefieldControl> {
#pragma warning disable 649
        /// <summary>
        /// concentration aka the level set
        /// </summary>
        protected SinglePhaseField phi;

        /// <summary>
        /// concentration (linearization point)
        /// </summary>
        protected SinglePhaseField phi0;

        /// <summary>
        /// Transport velocity
        /// </summary>
        protected VectorField<SinglePhaseField> gradPhi0;

        /// <summary>
        /// curvature
        /// </summary>
        protected SinglePhaseField Curvature;

        /// <summary>
        /// curvature, directly evaluated
        /// </summary>
        protected SinglePhaseField DCurvature;

        /// <summary>
        /// chemical potential
        /// </summary>
        protected SinglePhaseField mu;

        /// <summary>
        /// residual of 'c(phi)'-equation
        /// </summary>
        protected SinglePhaseField phi_Resi;

        /// <summary>
        /// residual of 'mu'-equation
        /// </summary>
        protected SinglePhaseField mu_Resi;

        /// <summary>
        /// residual of 'mu'-equation
        /// </summary>
        protected SinglePhaseField curvature_Resi;

        /// <summary>
        /// Transport velocity
        /// </summary>
        protected SinglePhaseField[] Velocity;
#pragma warning restore 649

        /// <summary>
        /// LevelSet
        /// </summary>
        LevelSet LevSet;

        /// <summary>
        /// DG LevelSet
        /// </summary>
        protected SinglePhaseField DGLevSet;

        ///// <summary>
        ///// LevelSetTracker
        ///// </summary>
        //LevelSetTracker LsTrk;

        ///// <summary>
        ///// DummyLevelSet, as the Phasefield is in DG not XDG
        ///// </summary>
        //LevelSet DummyLevSet;

        ///// <summary>
        ///// DummyLevelSetTracker
        ///// </summary>
        //LevelSetTracker DummyLsTrk;

        LevelSet CorrectionLevSet;

        /// <summary>
        /// LevelSetTracker, used for algebraic correction operations
        /// </summary>
        LevelSetTracker CorrectionLsTrk;

        int m_HMForder;

        ///// <summary>
        ///// Settings for solving the Phasefield equations
        ///// </summary>
        //PhasefieldControl Control;

        /// <summary>
        /// Cahn-Hilliard spatial Operator
        /// </summary>
        DifferentialOperator CHOp;

        /// <summary>
        /// Operator for Jacobian
        /// </summary>
        DifferentialOperator JacobiOp;

        /// <summary>
        /// Boundary Condition map for Cahn Hilliard
        /// </summary>
        BoundaryCondMap<BoundaryType> m_bcMap;

        /// <summary>
        /// Level Set Timestepper
        /// </summary>
        Solution.XdgTimestepping.XdgBDFTimestepping m_Timestepper;

        /// <summary>
        /// GridData
        /// </summary>
        IGridData GridData;

        /// <summary>
        /// MultigridSequence
        /// </summary>
        AggregationGridData[] mgSeq;

        /// <summary>
        /// Control of base solver that calls this level set method
        /// </summary>
        AppControl ParentControl;

        /// <summary>
        /// Phasefield instantiation, constructor
        /// </summary>
        public Phasefield(PhasefieldControl _Control, LevelSet _LevSet, SinglePhaseField _DGLevSet, LevelSetTracker _LsTrk, SinglePhaseField[] _Velocity, IGridData _GridData, AppControl _control, AggregationGridData[] _mgSeq) {
            if(_Control != null)
                this.Init(_Control);
            else
                this.Init(new PhasefieldControl());

            LevSet = _LevSet;
            //LsTrk = _LsTrk;
            DGLevSet = _DGLevSet;
            Velocity = _Velocity;
            GridData = _GridData;
            ParentControl = _control;
            mgSeq = _mgSeq;

        }

        /// <summary>
        /// Updating Phasefield after changes in Grid and/or Basis
        /// </summary>
        public void UpdateFields(LevelSet _LevSet, SinglePhaseField _DGLevSet, LevelSetTracker _LsTrk, SinglePhaseField[] _Velocity, IGridData _GridData, AppControl _control, AggregationGridData[] _mgSeq) {
            // check if signature of external and local Grid changed
            if(this.GridData != _GridData) {
                Console.WriteLine("Grid changed, adapting Phasefield");
                LevSet = _LevSet;
                //LsTrk = _LsTrk;               
                DGLevSet = _DGLevSet;
                Velocity = _Velocity;
                GridData = _GridData;
                ParentControl = _control;
                mgSeq = _mgSeq;

                double Cahn_Old = this.Control.cahn;

                CreateFields();

                if(Math.Abs(Cahn_Old - this.Control.cahn) > 1e-3 * Cahn_Old)
                    //RelaxationStep();
                    ReInit(Cahn_Old, this.Control.cahn);

                m_SOperator = GetOperatorInstance(GridData.SpatialDimension);
                CreateEquationsAndSolvers(null);

                // remember last cahn number for potential reinit
                Cahn_Reinit = this.Control.cahn;
            }
        }

   

        ///// <summary>
        ///// DG coordinates of <see cref="CurrentState"/> in a single vector
        ///// </summary>
        //public override CoordinateVector CurrentStateVector {
        //    get;
        //    protected set;
        //}

        ///// <summary>
        ///// DG coordinates of <see cref="CurrentResidual"/> in a single vector
        ///// </summary>
        //public override CoordinateVector CurrentResidualVector {
        //    get;
        //    protected set;
        //}
    

        /// <summary>
        /// Create Fields with same basis as DG Level Set
        /// </summary>
        protected override void CreateFields() {

            // create fields
            phi0 = new SinglePhaseField(DGLevSet.Basis, "phi0");
            gradPhi0 = new VectorField<SinglePhaseField>(DGLevSet.GridDat.SpatialDimension.ForLoop(d => new SinglePhaseField(DGLevSet.Basis, "dPhiDG_dx[" + d + "]")));
            Curvature = new SinglePhaseField(new Basis(phi0.GridDat, 0), VariableNames.Curvature);
            DCurvature = new SinglePhaseField(new Basis(phi0.GridDat, 0), "D" + VariableNames.Curvature);
            phi = new SinglePhaseField(DGLevSet.Basis, "phi");
            phi.Acc(1.0, DGLevSet);
            mu = new SinglePhaseField(DGLevSet.Basis, "mu");
            phi_Resi = new SinglePhaseField(DGLevSet.Basis, "phi_Resi");
            mu_Resi = new SinglePhaseField(DGLevSet.Basis, "mu_Resi");
            curvature_Resi = new SinglePhaseField(Curvature.Basis, "curvature_Resi");

            // residuals:
            var solFields = InstantiateSolutionFields();
            CurrentStateVector = new CoordinateVector(solFields);

            // residuals:
            var resFields = InstantiateResidualFields();
            CurrentResidualVector = new CoordinateVector(resFields);

            //// Dummy Level Set
            //DummyLevSet = new LevelSet(new Basis(this.GridData, 1), "Levset");
            //DummyLevSet.AccConstant(-1);
            //this.DummyLsTrk = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B" }, DummyLevSet);
            //this.DummyLsTrk.UpdateTracker(0.0);

            // Actual Level Set used for correction operations
            CorrectionLevSet = new LevelSet(phi.Basis, "Levset");
            this.CorrectionLsTrk = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 2, new string[] { "A", "B" }, CorrectionLevSet);
            CorrectionLevSet.Clear();
            CorrectionLevSet.Acc(1.0, phi);
            this.CorrectionLsTrk.UpdateTracker(0.0);

            // set coefficients
            SetCHCoefficents();
        }

        DifferentialOperator m_SOperator;

        /// <summary>
        /// Cache for <see cref="GetOperatorInstance"/>
        /// </summary>
        public override DifferentialOperator SOperator {
            get {
                if(m_SOperator == null) {
                    m_SOperator = GetOperatorInstance(this.DGLevSet.GridDat.SpatialDimension);
                }
                return m_SOperator;
            }
        }

        protected override DifferentialOperator GetOperatorInstance(int D) {
            // create operator
            // ===============
            {

                //BoundaryCondMap<BoundaryType> PoissonBcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "T");

                m_bcMap = new BoundaryCondMap<BoundaryType>(this.GridData, BoundaryTranslator(ParentControl.BoundaryValues), "phi");

                #region variables

                //create Parameter and Variablelists
                string[] paramVar = VariableNames.VelocityVector(D).Cat("phi0");
                string[] domainVar = new string[] { "phi" };
                string[] codomainVar = new string[] { "Res_phi" };

                switch(this.Control.ModTyp) {
                    case PhasefieldControl.ModelType.modelA:
                    break;
                    case PhasefieldControl.ModelType.modelB:
                    domainVar = domainVar.Cat("mu");
                    codomainVar = codomainVar.Cat("Res_mu");
                    break;
                    case PhasefieldControl.ModelType.modelC:
                    default:
                    throw new NotImplementedException();
                    //break;
                }

                //switch(this.Control.CurvatureCorrectionType) {
                //    case PhasefieldControl.CurvatureCorrection.FullyCoupled:
                //    domainVar = domainVar.Cat(VariableNames.Curvature);
                //    codomainVar = codomainVar.Cat("Res_" + VariableNames.Curvature);
                //    break;
                //    case PhasefieldControl.CurvatureCorrection.DirectCoupledIterative:
                //    case PhasefieldControl.CurvatureCorrection.DirectCoupledOnce:
                //    domainVar = domainVar.Cat(VariableNames.Curvature);
                //    codomainVar = codomainVar.Cat("Res_" + VariableNames.Curvature);
                //    paramVar = paramVar.Cat("D" + VariableNames.Curvature);
                //    break;
                //    case PhasefieldControl.CurvatureCorrection.None:
                //    default:
                //    break;
                //}

                #endregion

                CHOp = new DifferentialOperator(
                            domainVar,
                            paramVar,
                            codomainVar,
                            (DomainVariableDegrees, ParameterDegrees, CodomainVariableDegrees) => 4 * DomainVariableDegrees[0]//QuadOrderFunc.NonLinear(3)
                            );

                CHOp.ParameterUpdates.Add(CompleteParameterUpdate);
                CHOp.ParameterFactories.Add(AllocateParameters);

                #region equation components

                // convection term
                CHOp.EquationComponents["Res_phi"].Add(
                new phi_Flux(D, () => Velocity, m_bcMap)
                );

                switch(this.Control.ModTyp) {
                    case PhasefieldControl.ModelType.modelA:

                    CHOp.EquationComponents["Res_phi"].Add(
                    new phi_Source(this.Control.diff)
                    );


                    CHOp.EquationComponents["Res_phi"].Add(
                        new mu_Diffusion(D, this.Control.penalty_poisson, this.Control.cahn * this.Control.diff.Sqrt(), m_bcMap)
                        );

                    //switch(this.Control.CurvatureCorrectionType) {
                    //    case PhasefieldControl.CurvatureCorrection.FullyCoupled:
                    //    CHOp.EquationComponents["Res_phi"].Add(
                    //        new phi_CurvatureCorrection(D, this.Control.cahn * this.Control.diff.Sqrt())
                    //        );

                    //    CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                    //        new curvature_Source(D)
                    //        );

                    //    CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                    //        new curvature_Divergence(D, this.Control.penalty_poisson, 0.001 / this.Control.cahn)
                    //        );
                    //    break;
                    //    case PhasefieldControl.CurvatureCorrection.DirectCoupledIterative:
                    //    case PhasefieldControl.CurvatureCorrection.DirectCoupledOnce:
                    //    CHOp.EquationComponents["Res_phi"].Add(
                    //        new phi_CurvatureCorrection(D, this.Control.cahn * this.Control.diff.Sqrt())
                    //        );

                    //    CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                    //        new curvature_Direct(D)
                    //        );
                    //    break;
                    //    case PhasefieldControl.CurvatureCorrection.None:
                    //    default:
                    //    break;
                    //}
                    break;
                    case PhasefieldControl.ModelType.modelB:

                    CHOp.EquationComponents["Res_phi"].Add(
                        new phi_Diffusion(D, this.Control.penalty_poisson, this.Control.diff, this.Control.lambda, m_bcMap)
                        );

                    CHOp.EquationComponents["Res_mu"].Add(
                        new mu_Diffusion(D, this.Control.penalty_poisson, this.Control.cahn, m_bcMap)
                        );

                    CHOp.EquationComponents["Res_mu"].Add(
                        new mu_Source()
                        );

                    //switch(this.Control.CurvatureCorrectionType) {
                    //    case PhasefieldControl.CurvatureCorrection.FullyCoupled:
                    //    CHOp.EquationComponents["Res_mu"].Add(
                    //        new phi_CurvatureCorrection(D, this.Control.cahn)
                    //        );

                    //    CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                    //        new curvature_Source(D)
                    //        );

                    //    CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                    //        new curvature_Divergence(D, this.Control.penalty_poisson, 0.001 / this.Control.cahn)
                    //        );
                    //    break;
                    //    case PhasefieldControl.CurvatureCorrection.DirectCoupledIterative:
                    //    case PhasefieldControl.CurvatureCorrection.DirectCoupledOnce:
                    //    CHOp.EquationComponents["Res_mu"].Add(
                    //        new phi_CurvatureCorrection(D, this.Control.cahn)
                    //        );

                    //    CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                    //        new curvature_Direct(D)
                    //        );
                    //    break;
                    //    case PhasefieldControl.CurvatureCorrection.None:
                    //    default:
                    //    break;
                    //}
                    break;
                    
                    case PhasefieldControl.ModelType.modelC:
                    throw new NotImplementedException();
                    //break;
                    
                    default:
                    throw new ArgumentOutOfRangeException();
                    //break;
                }

                #endregion

                // temporal derivative
                double[] MassScales = new double[domainVar.Length];
                MassScales[0] = 1.0;
                CHOp.TemporalOperator = new ConstantTemporalOperator(CHOp, MassScales);

                CHOp.LinearizationHint = LinearizationHint.GetJacobiOperator;

                CHOp.Commit();

                return CHOp;
            }
        }

        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            DGField[] SolutionFields = new DGField[] { phi };

            switch(this.Control.ModTyp) {
                case PhasefieldControl.ModelType.modelB:
                SolutionFields = SolutionFields.Cat(mu);
                break;
                case PhasefieldControl.ModelType.modelA:
                case PhasefieldControl.ModelType.modelC:
                default:
                break;
            }

            //if(this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None) {
            //    SolutionFields.Cat(Curvature);
            //}

            return SolutionFields;
        }


        //protected override IEnumerable<DGField> InstantiateParameterFields() {
        //    DGField[] ParameterFields = Velocity.ToArray().Cat(phi0);
        //    return ParameterFields;
        //}


        public override IEnumerable<DGField> InstantiateResidualFields() {
            DGField[] ResidualFields = new DGField[] { phi_Resi };

            switch(this.Control.ModTyp) {
                case PhasefieldControl.ModelType.modelB:
                ResidualFields = ResidualFields.Cat(mu_Resi);
                break;
                case PhasefieldControl.ModelType.modelA:
                case PhasefieldControl.ModelType.modelC:
                default:
                break;
            }

            //if(this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None) {
            //    ResidualFields.Cat(curvature_Resi);
            //}

            return ResidualFields;
        }

        CoordinateVector m_CurrentSolution = null;

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution {
            get {
                m_CurrentSolution = new CoordinateVector(this.phi);

                switch(Control.ModTyp) {
                    case PhasefieldControl.ModelType.modelA:
                    break;
                    case PhasefieldControl.ModelType.modelB:
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.mu));
                    break;
                    case PhasefieldControl.ModelType.modelC:
                    default:
                    break;
                }

                //if(this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None) {
                //    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.Curvature));
                //}

                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual {
            get {
                m_CurrentResidual = new CoordinateVector(this.phi_Resi);

                switch(Control.ModTyp) {
                    case PhasefieldControl.ModelType.modelA:
                    break;
                    case PhasefieldControl.ModelType.modelB:
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.mu_Resi));
                    break;
                    case PhasefieldControl.ModelType.modelC:
                    default:
                    break;
                }

                //if(this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None) {
                //    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.curvature_Resi));
                //}

                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// default plotting
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            string caseStr = "";

            DGField[] Fields = new DGField[0];
            SinglePhaseField phiDist = ComputeDistanceField();
            Fields = Fields.Cat(this.phi, this.mu, this.Velocity, this.gradPhi0, this.Curvature, this.DCurvature, phiDist);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "Phasefield-" + timestepNo + caseStr, phystime, superSampling);
        }

        private SinglePhaseField ComputeDistanceField() {
            SinglePhaseField phiDist = new SinglePhaseField(phi.Basis);
            GridData GridDat = (GridData)(phi.GridDat);

            // calculate distance field phiDist = 0.5 * log(Max(1+phi, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn
            // ===================
            phiDist.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Phi
                    // -----------------------------
                    phi.Evaluate(j0, Len, NS, result);

                    // compute the pointwise values of the new level set
                    // -----------------------------

                    result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * Math.Sqrt(2) * this.Control.cahn);
                }
            );

            return phiDist;
        }

        void CompleteParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            UpdateParameter(t);
        }

        private void UpdateParameter(double time) {
            int D = this.GridData.SpatialDimension;
            this.phi0.Clear();
            this.phi0.Acc(1.0, this.phi);

            this.gradPhi0.Clear();
            this.gradPhi0.Gradient(1.0, this.phi0);

            //if(this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledIterative) {
            //    VectorField<SinglePhaseField> filtgrad;
            //    CurvatureAlgorithmsForLevelSet.CurvatureDriver(
            //                 CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected,
            //                 CurvatureAlgorithmsForLevelSet.FilterConfiguration.Phasefield,
            //                 this.DCurvature, out filtgrad, CorrectionLsTrk,
            //                 this.DCurvature.Basis.Degree * 2,
            //                 this.phi);
            //}
        }

        (string ParameterName, DGField ParamField)[] AllocateParameters(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            DGField[] prms = Velocity.ToArray().Cat(phi0);

            //if(this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledIterative || this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledOnce)
            //    prms.Cat(DCurvature);


            if(prms.Length != this.CHOp.ParameterVar.Count)
                throw new ApplicationException("mismatch between params in the operator and allocated fields.");

            var ret = new ValueTuple<string, DGField>[prms.Length];
            for(int iParam = 0; iParam < prms.Length; iParam++) {
                ret[iParam] = (this.CHOp.ParameterVar[iParam], prms[iParam] as DGField);
            }

            return ret;
        }
    }
}

