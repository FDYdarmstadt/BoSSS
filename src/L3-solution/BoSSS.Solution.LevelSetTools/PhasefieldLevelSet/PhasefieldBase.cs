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

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Base class for a level set based on the Cahn-Hilliard Phasefield equation
    /// Insert Latex style equations here
    /// </summary>
    public partial class Phasefield : Application
    {
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
        protected VectorField<SinglePhaseField> Velocity;
#pragma warning restore 649

        /// <summary>
        /// LevelSet
        /// </summary>
        LevelSet LevSet;

        /// <summary>
        /// DG LevelSet
        /// </summary>
        protected SinglePhaseField DGLevSet;

        /// <summary>
        /// LevelSetTracker
        /// </summary>
        LevelSetTracker LsTrk;

        /// <summary>
        /// DummyLevelSet, as the Phasefield is in DG not XDG
        /// </summary>
        LevelSet DummyLevSet;

        /// <summary>
        /// DummyLevelSetTracker
        /// </summary>
        LevelSetTracker DummyLsTrk;

        LevelSet CorrectionLevSet;

        /// <summary>
        /// LevelSetTracker, used for algebraic correction operations
        /// </summary>
        LevelSetTracker CorrectionLsTrk;

        int m_HMForder;

        /// <summary>
        /// Settings for solving the Phasefield equations
        /// </summary>
        PhasefieldControl Control;

        /// <summary>
        /// Cahn-Hilliard spatial Operator
        /// </summary>
        SpatialOperator CHOp;

        /// <summary>
        /// Operator for Jacobian
        /// </summary>
        SpatialOperator JacobiOp;

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
        public Phasefield(PhasefieldControl _Control, LevelSet _LevSet, SinglePhaseField _DGLevSet,  LevelSetTracker _LsTrk, VectorField<SinglePhaseField> _Velocity, IGridData _GridData, AppControl _control, AggregationGridData[] _mgSeq)
        {
            if (_Control != null)
                Control = _Control;
            else
                Control = new PhasefieldControl();

            LevSet = _LevSet;
            LsTrk = _LsTrk;
            DGLevSet = _DGLevSet;
            Velocity = _Velocity;
            GridData = _GridData;
            ParentControl = _control;
            mgSeq = _mgSeq;

        }

        /// <summary>
        /// Updating Phasefield after changes in Grid and/or Basis
        /// </summary>
        public void UpdateFields(LevelSet _LevSet, SinglePhaseField _DGLevSet, LevelSetTracker _LsTrk, VectorField<SinglePhaseField> _Velocity, IGridData _GridData, AppControl _control, AggregationGridData[] _mgSeq)
        {
            // check if signature of external and local Grid changed
            if (this.GridData != _GridData)
            {
                Console.WriteLine("Grid changed, adapting Phasefield");
                LevSet = _LevSet;
                LsTrk = _LsTrk;
                DGLevSet = _DGLevSet;
                Velocity = _Velocity;
                GridData = _GridData;
                ParentControl = _control;
                mgSeq = _mgSeq;

                double Cahn_Old = this.Control.cahn;

                CreateFields();

                if (Math.Abs(Cahn_Old - this.Control.cahn) > 1e-3 * Cahn_Old)
                    //RelaxationStep();
                    ReInit(Cahn_Old, this.Control.cahn);

                CreateEquationsAndSolvers(null);                

                // remember last cahn number for potential reinit
                Cahn_Reinit = this.Control.cahn;
            }
        }

        /// <summary>
        /// Create Fields with same basis as DG Level Set
        /// </summary>
        protected override void CreateFields()
        {
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

            // Dummy Level Set
            DummyLevSet = new LevelSet(new Basis(this.GridData, 1), "Levset");
            DummyLevSet.AccConstant(-1);
            this.DummyLsTrk = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B" }, DummyLevSet);
            this.DummyLsTrk.UpdateTracker();

            // Actual Level Set used for correction operations
            CorrectionLevSet = new LevelSet(phi.Basis, "Levset");
            this.CorrectionLsTrk = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 2, new string[] { "A", "B" }, CorrectionLevSet);
            CorrectionLevSet.Clear();
            CorrectionLevSet.Acc(1.0, phi);
            this.CorrectionLsTrk.UpdateTracker();

            // set coefficients
            SetCHCoefficents();
        }

        /// <summary>
        /// Includes assembly of the matrix.
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        {
            using (FuncTrace tr = new FuncTrace())
            {
                // create operator
                // =============
                this.CreateOperator();

                // create solver
                // =============
                this.CreateTimestepper();               
            }            
        }

        private void CreateOperator()
        {
            // create operator
            // ===============
            {

                int D = this.GridData.SpatialDimension;
                m_HMForder = phi.Basis.Degree * 2;
                double _D = D;
                double penalty_base = (phi.Basis.Degree + 1) * (phi.Basis.Degree + _D) / _D;

                // Get this from where?
                double penalty_factor = this.Control.penalty_poisson * penalty_base;

                //BoundaryCondMap<BoundaryType> PoissonBcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "T");

                MultidimensionalArray LengthScales;
                if (this.GridData is GridData)
                {
                    LengthScales = ((GridData)this.GridData).Cells.cj;
                }
                else if (this.GridData is AggregationGridData)
                {
                    LengthScales = ((AggregationGridData)this.GridData).AncestorGrid.Cells.cj;
                }
                else
                {
                    throw new NotImplementedException();
                }

                m_bcMap = new BoundaryCondMap<BoundaryType>(this.GridData, BoundaryTranslator(ParentControl.BoundaryValues), "phi");

                #region variables

                //create Parameter and Variablelists
                string[] paramVar = VariableNames.VelocityVector(D).Cat("phi0");
                string[] domainVar = new string[] { "phi" };
                string[] codomainVar = new string[] { "Res_phi" };

                switch (this.Control.ModTyp)
                {
                    case PhasefieldControl.ModelType.modelA:
                        break;
                    case PhasefieldControl.ModelType.modelB:
                        domainVar = domainVar.Cat("mu");
                        codomainVar = codomainVar.Cat("Res_mu");
                        break;
                    case PhasefieldControl.ModelType.modelC:
                    default:
                        throw new NotImplementedException();
                        break;
                }

                switch (this.Control.CurvatureCorrectionType)
                {
                    case PhasefieldControl.CurvatureCorrection.FullyCoupled:
                        domainVar = domainVar.Cat(VariableNames.Curvature);
                        codomainVar = codomainVar.Cat("Res_" + VariableNames.Curvature);
                        break;
                    case PhasefieldControl.CurvatureCorrection.DirectCoupledIterative:
                    case PhasefieldControl.CurvatureCorrection.DirectCoupledOnce:
                        domainVar = domainVar.Cat(VariableNames.Curvature);
                        codomainVar = codomainVar.Cat("Res_" + VariableNames.Curvature);
                        paramVar = paramVar.Cat("D" + VariableNames.Curvature);
                        break;
                    case PhasefieldControl.CurvatureCorrection.None:
                    default:
                        break;                    
                }

                #endregion

                CHOp = new SpatialOperator(
                            domainVar,
                            paramVar,
                            codomainVar,
                            QuadOrderFunc.NonLinear(3)
                            );

                CHOp.ParameterUpdate = new PFParameterUpdate(this);

                #region equation components

                // convection term
                CHOp.EquationComponents["Res_phi"].Add(
                new phi_Flux(D, m_bcMap)
                );

                switch (this.Control.ModTyp)
                {
                    case PhasefieldControl.ModelType.modelA:

                        CHOp.EquationComponents["Res_phi"].Add(
                        new phi_Source(this.Control.diff)
                        );


                        CHOp.EquationComponents["Res_phi"].Add(
                            new mu_Diffusion(D, penalty_factor, LengthScales, this.Control.cahn * this.Control.diff.Sqrt(), m_bcMap)
                            );

                        switch (this.Control.CurvatureCorrectionType)
                        {
                            case PhasefieldControl.CurvatureCorrection.FullyCoupled:
                                CHOp.EquationComponents["Res_phi"].Add(
                                    new phi_CurvatureCorrection(D, this.Control.cahn * this.Control.diff.Sqrt())
                                    );

                                CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                    new curvature_Source(D)
                                    );

                                CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                    new curvature_Divergence(D, penalty_factor, 0.001 / this.Control.cahn, LengthScales)
                                    );
                                break;
                            case PhasefieldControl.CurvatureCorrection.DirectCoupledIterative:
                            case PhasefieldControl.CurvatureCorrection.DirectCoupledOnce:
                                CHOp.EquationComponents["Res_phi"].Add(
                                    new phi_CurvatureCorrection(D, this.Control.cahn * this.Control.diff.Sqrt())
                                    );

                                CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                    new curvature_Direct(D)
                                    );
                                break;
                            case PhasefieldControl.CurvatureCorrection.None:
                            default:
                                break;
                        }
                        break;
                    case PhasefieldControl.ModelType.modelB:

                        CHOp.EquationComponents["Res_phi"].Add(
                            new phi_Diffusion(D, penalty_factor, LengthScales, this.Control.diff, this.Control.lambda, m_bcMap)
                            );

                        CHOp.EquationComponents["Res_mu"].Add(
                            new mu_Diffusion(D, penalty_factor, LengthScales, this.Control.cahn, m_bcMap)
                            );

                        CHOp.EquationComponents["Res_mu"].Add(
                            new mu_Source()
                            );

                        switch (this.Control.CurvatureCorrectionType)
                        {
                            case PhasefieldControl.CurvatureCorrection.FullyCoupled:
                                CHOp.EquationComponents["Res_mu"].Add(
                                    new phi_CurvatureCorrection(D, this.Control.cahn)
                                    );

                                CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                    new curvature_Source(D)
                                    );

                                CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                    new curvature_Divergence(D, penalty_factor, 0.001 / this.Control.cahn, LengthScales)
                                    );
                                break;
                            case PhasefieldControl.CurvatureCorrection.DirectCoupledIterative:
                            case PhasefieldControl.CurvatureCorrection.DirectCoupledOnce:
                                CHOp.EquationComponents["Res_mu"].Add(
                                    new phi_CurvatureCorrection(D, this.Control.cahn)
                                    );

                                CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                    new curvature_Direct(D)
                                    );
                                break;
                            case PhasefieldControl.CurvatureCorrection.None:
                            default:
                                break;
                        }
                        break;
                    case PhasefieldControl.ModelType.modelC:
                        throw new NotImplementedException();
                        break;
                    default:
                        throw new ArgumentOutOfRangeException();
                        break;
                }

                #endregion
                
                CHOp.ParameterUpdate = new PFParameterUpdate(this);

                CHOp.Commit();
                JacobiOp = CHOp._GetJacobiOperator(D);
            }
        }

        RungeKuttaScheme rksch = null;
        int bdfOrder = -1000;
        private void CreateTimestepper()
        {
            switch (this.Control.TimeSteppingScheme)
            {
                case TimeSteppingScheme.RK_ImplicitEuler:
                    {
                        rksch = RungeKuttaScheme.ImplicitEuler;
                        break;
                    }
                case TimeSteppingScheme.RK_CrankNic:
                    {
                        rksch = RungeKuttaScheme.CrankNicolson;
                        break;
                    }
                case TimeSteppingScheme.CrankNicolson:
                    {
                        //do not instantiate rksch, use bdf instead
                        bdfOrder = -1;
                        break;
                    }
                case TimeSteppingScheme.ImplicitEuler:
                    {
                        //do not instantiate rksch, use bdf instead
                        bdfOrder = 1;
                        break;
                    }
                default:
                    {
                        if (this.Control.TimeSteppingScheme.ToString().StartsWith("BDF"))
                        {
                            //do not instantiate rksch, use bdf instead
                            bdfOrder = Convert.ToInt32(this.Control.TimeSteppingScheme.ToString().Substring(3));
                            break;
                        }
                        else
                            throw new NotImplementedException();
                    }

            }


            if (rksch == null)
            {
                m_Timestepper = new XdgBDFTimestepping(
                    this.CurrentSolution.Mapping.Fields,
                    this.CurrentResidual.Mapping.Fields,
                    this.DummyLsTrk,
                    false,
                    DelComputeOperatorMatrix, null, null,
                    (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? bdfOrder : 1,
                    LevelSetHandling.None, MassMatrixShapeandDependence.IsTimeDependent, SpatialOperatorType.Nonlinear,
                    this.MassScale,
                    this.MgConfig, mgSeq,
                    new[] { this.DummyLsTrk.GetSpeciesId("A") }, phi.Basis.Degree * 2,
                    0.0,
                    false,
                    this.Control.NonLinearSolver,
                    this.Control.LinearSolver
                    );
                m_Timestepper.m_ResLogger = new ResidualLogger(this.MPIRank, null, new Guid());
                m_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                m_Timestepper.m_ResLogger.WriteResidualsToTextFile = false;
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        CoordinateVector m_CurrentSolution = null;

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution
        {
            get
            {
                m_CurrentSolution = new CoordinateVector(this.phi);

                switch (Control.ModTyp)
                {
                    case PhasefieldControl.ModelType.modelA:
                        break;
                    case PhasefieldControl.ModelType.modelB:
                        m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.mu));
                        break;
                    case PhasefieldControl.ModelType.modelC:
                    default:
                        break;
                }

                if (this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None)
                {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.Curvature));
                }
                
                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual
        {
            get
            {                
                m_CurrentResidual = new CoordinateVector(this.phi_Resi);

                switch (Control.ModTyp)
                {
                    case PhasefieldControl.ModelType.modelA:
                        break;
                    case PhasefieldControl.ModelType.modelB:
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.mu_Resi));
                        break;
                    case PhasefieldControl.ModelType.modelC:
                    default:
                        break;
                }

                if (this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None)
                {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.curvature_Resi));
                }
                
                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// Computation of operator matrix
        /// </summary>
        void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time)
        {

            //create mappings
            var codMap = Mapping;
            var domMap = Mapping;

            DGField[] prms;
            prms = ArrayTools.Cat(this.Velocity, phi0);

            if (this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledIterative ||
                this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledOnce)
            {
                prms = prms.Cat(this.DCurvature);
            }

            if (OpMtx != null)
            {
                // ++++++++++++++++++++++++++++++++
                // create matrix and affine vector:
                // ++++++++++++++++++++++++++++++++

                IEvaluatorLinear mtxBuilder;

                // for completeness, FDJacobian is slower, can be switched on if needed
                bool UseFDJacobian = false;

                if (UseFDJacobian)
                {
                    mtxBuilder = CHOp.GetFDJacobianBuilder(CurrentState, prms, codMap);
                    mtxBuilder.time = time;
                    mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                }
                else
                {
                    var JacParams = JacobiOp.ParameterUpdate;
                    var TmpParams = JacParams.AllocateParameters(CurrentState, prms);
                    var map = Mapping;// new CoordinateMapping(CurrentState);

                    var JacBuilder = JacobiOp.GetMatrixBuilder(map, TmpParams, map);
                    JacobiOp.ParameterUpdate.PerformUpdate(CurrentState, TmpParams);
                    JacParams.PerformUpdate(CurrentState, TmpParams);
                    JacBuilder.ComputeMatrix(OpMtx, OpAffine);
                }

            }
            else
            {
                // ++++++++++++++++++++++++++++++++
                // evaluate the operator
                // ++++++++++++++++++++++++++++++++
                var eval = CHOp.GetEvaluatorEx(CurrentState, prms, codMap);
                eval.time = time;
                eval.Evaluate(1.0, 1.0, OpAffine);
            }

            try
            {
                if (OpMtx != null)
                    OpMtx.CheckForNanOrInfM();
                OpAffine.CheckForNanOrInfV();
            }
            catch (ArithmeticException ae)
            {
                Console.WriteLine("Found NAN");

                foreach (DGField f in CurrentState)
                {
                    f.GetExtremalValues(out double min, out double max);
                    Console.WriteLine($"  Field {f.Identification} extremal values: {min}  ---  {max}");
                }

                throw ae;
            }

            #region deprecated, only Picard-Support
            /*
            SinglePhaseField Current_phi = (SinglePhaseField)(CurrentState[0]);

            phi0.Clear();
            phi0.Acc(1.0, Current_phi);
            gradPhi0.Clear();
            gradPhi0.Gradient(1.0, phi0);
            
            OpMtx.Clear();
            OpAffine.ClearEntries();
            var mb = CHOp.GetMatrixBuilder(Mapping, this.Velocity.ToArray().Cat(phi0).Cat(gradPhi0.ToArray()).Cat(Curvature), Mapping);
            mb.ComputeMatrix(OpMtx, OpAffine);
            */
            #endregion
        }

        MultigridOperator.ChangeOfBasisConfig[][] MgConfig
        {
            get
            {
                int p = this.phi.Basis.Degree;
                int NoOfLevels = mgSeq.Length;
                var config = new MultigridOperator.ChangeOfBasisConfig[NoOfLevels][];
                int m = 0;

                for (int iLevel = 0; iLevel < NoOfLevels; iLevel++)
                {

                    switch (Control.ModTyp)
                    {
                        case PhasefieldControl.ModelType.modelA:
                            m = 1;
                            config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m];
                            break;
                        case PhasefieldControl.ModelType.modelB:
                            m = 2;
                            config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m];
                            break;
                        case PhasefieldControl.ModelType.modelC:
                            throw new NotImplementedException();
                            break;
                        default:
                            throw new ArgumentOutOfRangeException();
                            break;
                    }

                    if (this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None)
                    {
                        config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m + 1];
                    }

                    config[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig()
                    {
                        VarIndex = new int[] { 0 },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                        DegreeS = new int[] { Math.Max(1, p - iLevel) }
                    };

                    switch (Control.ModTyp)
                    {
                        case PhasefieldControl.ModelType.modelA:
                            break;
                        case PhasefieldControl.ModelType.modelB:
                            config[iLevel][1] = new MultigridOperator.ChangeOfBasisConfig()
                            {
                                VarIndex = new int[] { 1 },
                                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                                DegreeS = new int[] { Math.Max(1, p - iLevel) }
                            };

                            break;
                        case PhasefieldControl.ModelType.modelC:
                            throw new NotImplementedException();
                            break;
                        default:
                            throw new ArgumentOutOfRangeException();
                            break;
                    }

                    if (this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None)
                    {
                        config[iLevel][m] = new MultigridOperator.ChangeOfBasisConfig()
                        {
                            VarIndex = new int[] { m },
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            DegreeS = new int[] { Math.Max(0, this.Curvature.Basis.Degree - iLevel) }
                        };
                    }
                }

                return config;
            }

        }

        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale
        {
            get
            {

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                double[] scale = new double[1];

                switch (Control.ModTyp)
                {
                    case PhasefieldControl.ModelType.modelA:
                        break;
                    case PhasefieldControl.ModelType.modelB:
                        scale = new double[2];
                        break;
                    case PhasefieldControl.ModelType.modelC:
                        throw new NotImplementedException();
                        break;
                    default:
                        throw new ArgumentOutOfRangeException();
                        break;
                }

                if (this.Control.CurvatureCorrectionType != PhasefieldControl.CurvatureCorrection.None)
                {
                    scale = new double[scale.Length + 1];
                }

                scale[0] = 1.0;

                R.Add(this.DummyLsTrk.GetSpeciesId("A"), scale);

                return R;
            }
        }        

        /// <summary>
        /// default plotting
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0)
        {
            string caseStr = "";
            
            DGField[] Fields = new DGField[0];
            SinglePhaseField phiDist = ComputeDistanceField();
            Fields = Fields.Cat(this.phi, this.mu, this.Velocity, this.gradPhi0, this.Curvature, this.DCurvature, phiDist);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "Phasefield-" + timestepNo + caseStr, phystime, superSampling);
        }

        private SinglePhaseField ComputeDistanceField()
        {
            SinglePhaseField phiDist = new SinglePhaseField(phi.Basis);
            GridData GridDat = (GridData)(phi.GridDat);

            // calculate distance field phiDist = 0.5 * log(Max(1+phi, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn
            // ===================
            phiDist.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result)
                { // ScalarFunction2
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

        private class PFParameterUpdate : IParameterUpdate
        {
            private Phasefield phasefield;

            public PFParameterUpdate(Phasefield _phasefield)
            {
                this.phasefield = _phasefield;
            }

            public DGField[] AllocateParameters(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar)
            {
                throw new NotImplementedException();
            }

            public void PerformUpdate(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar)
            {
                this.phasefield.UpdateParameter();
            }
        }

        private void UpdateParameter()
        {
            int D = this.GridData.SpatialDimension;
            this.phi0.Clear();
            this.phi0.Acc(1.0, this.phi);

            this.gradPhi0.Clear();
            this.gradPhi0.Gradient(1.0, this.phi0);

            if (this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledIterative)
            {
                VectorField<SinglePhaseField> filtgrad;
                CurvatureAlgorithmsForLevelSet.CurvatureDriver(
                             CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                             CurvatureAlgorithmsForLevelSet.FilterConfiguration.Phasefield,
                             this.DCurvature, out filtgrad, LsTrk,
                             this.DCurvature.Basis.Degree * 2,
                             this.phi);
            }
        }
    }

    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'phi'
    /// </summary>
    class phi_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace
    {

        public phi_Diffusion(int D, double penalty_const, MultidimensionalArray cj, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, cj, "mu") // note: in the equation for 'phi', we have the Laplacian of 'mu'
        {
            m_D = D;
            m_diff = __diff;
            m_boundaryCondMap = __boundaryCondMap;
            m_lambda = __lambda;
        }

        double m_diff;
        double m_lambda;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => null;

        protected override double g_Diri(ref CommonParamsBnd inp)
        {
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += (inp.Parameters_IN[d + 1]) * inp.Normal[d];
            }

            double v;
            if (UxN >= 0)
            {
                v = 1.0;
            }
            else
            {
                v = 0.0;
            }

            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp)
        {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp)
        {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType)
            {
                case BoundaryType.Wall:
                case BoundaryType.Slip:
                case BoundaryType.SlipSymmetry:
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                case BoundaryType.Inlet:
                case BoundaryType.Outlet:
                case BoundaryType.Outflow:
                case BoundaryType.Pressure_Dirichlet:
                default:
                    return true;
            }
        }

        public override double Nu(double[] x, double[] p, int jCell)
        {

            //double n = 0.0;
            //double D = 0.0;

            //for (int d = 0; d < m_D; d++)
            //{
            //    n += p[1 + m_D + d].Pow2();
            //}

            //D = 0.01 * Math.Exp(-8.0 * n * m_diff.Pow2());
            //D = 1e-7;

            if (m_lambda == 0.0)
            {
                return -m_diff; // -D;
            }
            else if (m_lambda > 0.0 && m_lambda <= 1.0)
            {
                return -m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0), 0.0);//-D * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0),0.0);
            }
            else
            {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp))
            {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                // inflow / outflow
                if (g_D == 0)
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                }
                else
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        //Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency 
                        //Acc += (nuA * 0.0) * (_vA) * nd;        // consistency, homogenous Neumann
                    }
                    Acc *= this.m_alpha;

                }
            }
            else
            {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }

    /// <summary>
    /// Transport flux for Cahn-Hilliard
    /// </summary>
    class phi_Flux : IVolumeForm, IEdgeForm, ISupportsJacobianComponent
    {
        public phi_Flux(int D, BoundaryCondMap<BoundaryType> __boundaryCondMap)
        {
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap.bndFunction["phi"];
        }

        int m_D;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ArgumentOrdering => new string[] { "phi" };

        public IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);

        public TermActivationFlags BoundaryEdgeTerms => InnerEdgeTerms | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN)
        {
            // expand for treatment of input functions, for now hardcode to -1.0
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double phi;
            if (UxN >= 0)
            {
                phi = _uIN[0];
            }
            else
            {
                phi = -1.0;//m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }

            return phi * UxN * _vIN;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
        {
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];
            }

            double phi;
            if (UxN >= 0)
            {
                phi = _uIN[0];
            }
            else
            {
                phi = _uOUT[0];
            }

            return phi * UxN * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double acc = 0;
            double phi = U[0];
            for (int d = 0; d < m_D; d++)
            {
                acc += phi * cpv.Parameters[d] * GradV[d];
            }

            return -acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new[] { this };
        }
    }

    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'mu'
    /// </summary>
    class mu_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace
    {

        public mu_Diffusion(int D, double penalty_const, MultidimensionalArray cj, double __cahn, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, cj, "phi") // note: in the equation for 'mu', we have the Laplacian of 'phi'
        {
            m_D = D;
            m_cahn = __cahn * __cahn;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap.bndFunction["phi"];
        }

        double m_cahn;
        int m_D;
        public override IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        protected override double g_Diri(ref CommonParamsBnd inp)
        {
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double v;
            if (UxN >= 0)
            {
                v = 0.0;
            }
            else
            {
                // for now set inflow value to -1.0
                v = -1.0; // m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp)
        {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp)
        {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType)
            {
                case BoundaryType.Wall:
                case BoundaryType.Slip:
                case BoundaryType.SlipSymmetry:
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                case BoundaryType.Inlet:
                case BoundaryType.Outlet:
                case BoundaryType.Outflow:
                case BoundaryType.Pressure_Dirichlet:
                default:
                    return true;
            }
        }

        public override double Nu(double[] x, double[] p, int jCell)
        {
            return -m_cahn;
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp))
            {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);                

                // inflow / outflow
                if (g_D != 0)
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                }
                else
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        //Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency 
                        //Acc += (nuA * 0.0) * (_vA) * nd;        // consistency Neumann
                    }

                    //double D0 = 1.0;
                    //Acc -= (nuA * _vA) * D0 * (_uA[0]-inp.Parameters_IN[inp.D])/ 0.01;        // Outflow boundary condition see S. Dong

                    Acc *= this.m_alpha;

                }
            }
            else
            {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }

    /// <summary>
    /// nonlinear source term in the 'phi'-equation
    /// </summary>
    class mu_Source : IVolumeForm, ISupportsJacobianComponent
    {

        public mu_Source()
        {

        }

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "mu", "phi" };

        public IList<string> ParameterOrdering => null; // new[] { "phi0" };

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {

            double mu = U[0];
            double phi = U[1];
            //double phi0 = cpv.Parameters[0];

            double Acc = 0;

            Acc += -mu;

            //Acc += (3 * phi0 * phi0 - 1.0) * phi - 2 * Math.Pow(phi0, 3.0); // linearized around phi0 (Taylor expansion)
            Acc += phi.Pow(3.0) - phi; // when using Newton with Jacobian linearization is not needed

            return Acc * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new IEquationComponent[] { new jacobi_mu_Source() };
        }

        private class jacobi_mu_Source : IVolumeForm
        {

            public jacobi_mu_Source()
            {
            }
            public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

            public IList<string> ArgumentOrdering => new[] { "mu", "phi" };

            public IList<string> ParameterOrdering => new[] { "phi0" };

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
            {
                double mu = U[0];
                double phi = U[1];
                double phi0 = cpv.Parameters[0];

                double Acc = 0;

                Acc -= mu;
                Acc += 3 * phi0.Pow2() * phi - phi; // linearized around c0 (Taylor expansion)

                return Acc * V;
            }
        }
    }

    /// <summary>
    /// source term of phi in Model A
    /// </summary>
    class phi_Source : IVolumeForm, ISupportsJacobianComponent
    {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public phi_Source(double _diff = 0.0)
        {
            m_diff = _diff;
        }



        //double m_lambda;
        //double m_epsilon;
        double m_diff;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "phi" };
        public IList<string> ParameterOrdering => null; //new[] { "phi0" };

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {

            double phi = U[0];
            //double phi0 = cpv.Parameters[0];
            // values seem shifted without this offset hack

            double Acc = 0;

            //Acc += (3.0 * phi0 * phi0 - 1.0) * phi - 2 * Math.Pow(phi0, 3.0); // linearized around (Taylor expansion)
            Acc += phi.Pow(3.0) - phi; // when using Newton with Jacobian linearization is not needed

            return m_diff * Acc * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new IEquationComponent[] { new jacobi_phi_Source(m_diff) };
        }

        private class jacobi_phi_Source : IVolumeForm
        {
            private double m_diff;

            public jacobi_phi_Source(double m_diff)
            {
                this.m_diff = m_diff;
            }

            public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

            public IList<string> ArgumentOrdering => new[] { "phi" };

            public IList<string> ParameterOrdering => new[] { "phi0" };

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
            {
                double phi = U[0];
                double phi0 = cpv.Parameters[0];

                double Acc = 0;

                Acc += 3 * phi0.Pow2() * phi - phi; // linearized around c0 (Taylor expansion)

                return Acc * V;
            }
        }
    }

    /// <summary>
    /// Correction term to counter along the interface diffusion
    /// </summary>
    class phi_CurvatureCorrection : IVolumeForm, ISupportsJacobianComponent
    {

        public phi_CurvatureCorrection(int D, double _cahn = 0.0)
        {
            m_cahn = _cahn.Pow2();
            m_D = D;
        }

        int m_D;
        double m_cahn;

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "phi", VariableNames.Curvature };
        public IList<string> ParameterOrdering => null;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double Acc = 0.0;
            double[] grad = new double[m_D];
            for (int d = 0; d < m_D; d++)
            {
                grad[d] = GradU[0, d];
            }
            Acc += grad.L2Norm();
            Acc *= m_cahn * U[1];

            // sign minus should be correct, plus produces more sensual results sign of curvature?
            return Acc * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            var VolDiff = new VolumeFormDifferentiator(this, m_D);            
            return new IEquationComponent[] { VolDiff };
        }
    }

    class curvature_Direct : IEquationComponent, IVolumeForm, ISupportsJacobianComponent
    {

        public curvature_Direct(int _D)
        {
            m_D = _D;
        }
        int m_D;

        public IList<string> ArgumentOrdering => new[] { VariableNames.Curvature };

        public IList<string> ParameterOrdering => new[] { "D" + VariableNames.Curvature };

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double Acc = 0.0;

            Acc += 1e-8 * (cpv.Parameters[0] - U[0]) * V;

            return Acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new[] { this };
        }

    }

    class curvature_Source : IEquationComponent, IVolumeForm, ISupportsJacobianComponent
    {

        public curvature_Source(int _D)
        {
            m_D = _D;
        }
        int m_D;

        public IList<string> ArgumentOrdering => new[] { VariableNames.Curvature };

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double Acc = 0.0;

            Acc += -U[0] * V;

            return Acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new[] { this };
        }
    }

    class curvature_Divergence : IVolumeForm, IEdgeForm, IEquationComponent, ISupportsJacobianComponent
    {
        public curvature_Divergence(int D, double penalty, double limit, MultidimensionalArray InverseLengthScales)
        {
            m_D = D;
            m_limit = limit;
            m_penalty = penalty;
            this.InverseLengthScales = InverseLengthScales;
        }

        int m_D;
        double m_limit;
        double m_penalty;

        public CellMask m_cells;

        public IList<string> ArgumentOrdering => new[] { "phi", VariableNames.Curvature };

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms => TermActivationFlags.GradUxGradV;

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;

        /// <summary>
        /// a little switch...
        /// </summary>
        protected double m_alpha = 1.0;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double acc = 0;

            double[] grad = new double[m_D];
            for (int d = 0; d < cpv.D; d++)
            {
                grad[d] = GradU[0, d];
                acc -= GradU[0, d] * GradV[d] * this.m_alpha;
            }
            double norm = Math.Max(grad.L2Norm(), m_limit);
            acc *= 1 / norm;

            return acc;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
        {
            double Acc = 0.0;

            double pnlty = this.GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);

            Acc += Flux(_uIN, _uOUT, _Grad_uIN, _Grad_uOUT, inp.Normal) * (_vIN - _vOUT);
            Acc *= this.m_alpha;

            return Acc;
        }

        private enum FluxType
        {
            Central,

            LaxFriedrich,

            Upwind,

            Godunov
        }

        private double Flux(double[] uIN, double[] uOUT, double[,] grad_uIN, double[,] grad_uOUT, Vector normal)
        {
            var FType = FluxType.Central;

            double Acc = 0.0;

            double[] MeanGrad = new double[m_D];
            double norm;

            switch (FType)
            {

                case FluxType.Central:

                    for (int d = 0; d < m_D; d++)
                        MeanGrad[d] = 0.5 * (grad_uIN[0, d] + grad_uOUT[0, d]);

                    norm = Math.Max(MeanGrad.L2Norm(), m_limit);

                    for (int d = 0; d < m_D; d++)
                        Acc += MeanGrad[d] / norm * normal[d];

                    break;
                case FluxType.Upwind:

                    // How to calculate upwind direction
                    double P = 0.0;

                    if (P > 0)
                    {
                        for (int d = 0; d < m_D; d++)
                            MeanGrad[d] = grad_uIN[0, d];
                    }
                    else
                    {
                        for (int d = 0; d < m_D; d++)
                            MeanGrad[d] = grad_uOUT[0, d];
                    }

                    norm = Math.Max(MeanGrad.L2Norm(), m_limit);

                    for (int d = 0; d < m_D; d++)
                        Acc += MeanGrad[d] / norm * normal[d];

                    break;
                case FluxType.LaxFriedrich:

                    double[] GradIN = new double[m_D];
                    double[] GradOUT = new double[m_D];
                    double gamma = 0.0;

                    for (int d = 0; d < m_D; d++)
                    {
                        GradIN[d] = grad_uIN[0, d];
                        GradOUT[d] = grad_uOUT[0, d];
                        MeanGrad[d] = GradIN[d] - GradOUT[d];
                    }

                    double normIN = Math.Max(GradIN.L2Norm(), m_limit);
                    double normOUT = Math.Max(GradOUT.L2Norm(), m_limit);
                    gamma = 0.1 / Math.Min(normIN, normOUT);

                    for (int d = 0; d < m_D; d++)
                        Acc += 0.5 * (GradIN[d] / normIN + GradOUT[d] / normOUT) * normal[d];

                    Acc -= gamma * MeanGrad.L2Norm();

                    break;
                case FluxType.Godunov:
                    break;
                default:
                    throw new NotSupportedException();
            }

            return Acc;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);

            double[] grad = new double[m_D];
            for (int d = 0; d < inp.D; d++)
            {
                grad[d] = _Grad_uA[0, d];
            }
            double norm = Math.Max(grad.L2Norm(), m_limit);

            bool IsDirichlet = false;

            if (IsDirichlet)
            {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = 0.0; //this.g_Diri(ref inp);

                for (int d = 0; d < inp.D; d++)
                {
                    double nd = inp.Normal[d];
                    Acc += (_Grad_uA[0, d] / norm) * (_vA) * nd;        // consistency
                    Acc += (_Grad_vA[d] / norm) * (_uA[0] - g_D) * nd;  // symmetry
                }
                Acc *= this.m_alpha;

                Acc -= (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty

            }
            else
            {

                double g_N = grad.InnerProd(inp.Normal) / norm;

                Acc += g_N * _vA * this.m_alpha;
            }
            return Acc;
        }

        /// <summary>
        /// Length scales used in <see cref="GetPenalty"/>
        /// </summary>
        protected MultidimensionalArray InverseLengthScales;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected virtual double GetPenalty(int jCellIn, int jCellOut)
        {
            double cj_in = InverseLengthScales[jCellIn];
            double mu = m_penalty * cj_in;
            if (jCellOut >= 0)
            {
                double cj_out = InverseLengthScales[jCellOut];
                mu = Math.Max(mu, m_penalty * cj_out);
            }

            return mu;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            var EdgeDiff = new EdgeFormDifferentiator(this, m_D);
            var VolDiff = new VolumeFormDifferentiator(this, m_D);
            return new IEquationComponent[] { EdgeDiff, VolDiff };
        }
    }


}

