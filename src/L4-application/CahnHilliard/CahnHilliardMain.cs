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
using NUnit.Framework;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.CahnHilliard {

    /// <summary>
    /// Benchmark application, solves a Poisson problem using the symmetric interior penalty (SIP) method.
    /// </summary>
    public class CahnHilliardMain : BoSSS.Solution.XdgTimestepping.DgApplicationWithSolver<CahnHilliardControl> {

#pragma warning disable 649
        /// <summary>
        /// concentration
        /// </summary>
        [InstantiateFromControlFile("c", "c", IOListOption.Always)]
        protected SinglePhaseField c;

        /// <summary>
        /// concentration (linearization point)
        /// </summary>
        [InstantiateFromControlFile("c0", "c", IOListOption.Always)]
        protected SinglePhaseField c0;

        /// <summary>
        /// Transport velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.LevelSetGradient0, VariableNames.LevelSetGradient1, VariableNames.LevelSetGradient2 },
            new[] { "c", "c", "c" },
            true, true,
            IOListOption.Always)]
        protected VectorField<SinglePhaseField> gradc0;

        /// <summary>
        /// hessian of the phasefield
        /// </summary>
        protected SinglePhaseField[,] hessc0;

        /// <summary>
        /// curvature
        /// </summary>
        [InstantiateFromControlFile(VariableNames.Curvature, VariableNames.Curvature, IOListOption.Always)]
        protected SinglePhaseField Curvature;

        /// <summary>
        /// curvature direct evaluation
        /// </summary>
        [InstantiateFromControlFile("D" + VariableNames.Curvature, VariableNames.Curvature, IOListOption.Always)]
        protected SinglePhaseField DCurvature;

        /// <summary>
        /// potential
        /// </summary>
        [InstantiateFromControlFile("phi", "c", IOListOption.Always)]
        protected SinglePhaseField phi;

        /// <summary>
        /// residual of 'c'-equation
        /// </summary>
        [InstantiateFromControlFile("c_Resi", "c", IOListOption.Always)]
        protected SinglePhaseField c_Resi;

        /// <summary>
        /// residual of 'phi'-equation
        /// </summary>
        [InstantiateFromControlFile("phi_Resi", "c", IOListOption.Always)]
        protected SinglePhaseField phi_Resi;

        /// <summary>
        /// residual of 'curvature'-equation
        /// </summary>
        [InstantiateFromControlFile("curvature_Resi", VariableNames.Curvature, IOListOption.Always)]
        protected SinglePhaseField curvature_Resi;

        /// <summary>
        /// Transport velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            null,
            true, true,
            IOListOption.Always)]
        protected VectorField<SinglePhaseField> Velocity;

        /// <summary>
        /// Transport velocity gradients
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY },
            null,
            true, true,
            IOListOption.Always)]
        protected VectorField<SinglePhaseField> VelocityGradX;

        /// <summary>
        /// Transport velocity gradients
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY },
            null,
            true, true,
            IOListOption.Always)]
        protected VectorField<SinglePhaseField> VelocityGradY;

        /// <summary>
        /// exact solution, to determine L2-Error, see also <see cref="CahnHilliardControl.ExactSolution_provided"/>.
        /// </summary>
        [InstantiateFromControlFile("cex", "cex", IOListOption.Always)]
        protected SinglePhaseField cex;
#pragma warning restore 649

        ///// <summary>
        ///// Not used presently, initialized to -1;
        ///// </summary>
        //LevelSet DummyLevset;

        /// <summary>
        /// Actual LevelSet, used for calculating Benchmark Quantities
        /// </summary>
        LevelSet RealLevSet;

        LevelSetTracker RealTracker;

        int m_HMForder;

        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            DGField[] SolutionFields = new DGField[] { c };

            switch(this.Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelB:
                SolutionFields = SolutionFields.Cat(phi);
                break;
                case CahnHilliardControl.ModelType.modelA:
                case CahnHilliardControl.ModelType.modelC:
                default:
                break;
            }

            if(this.Control.CurvatureCorrection) {
                SolutionFields.Cat(Curvature);
            }

            return SolutionFields;
        }

        /*
        protected override IEnumerable<DGField> InstantiateParameterFields() {
            return Velocity.Cat<DGField>(c0);
        }
        */

        public override IEnumerable<DGField> InstantiateResidualFields() {
            DGField[] ResidualFields = new DGField[] { c_Resi };

            switch(this.Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelB:
                ResidualFields = ResidualFields.Cat(phi_Resi);
                break;
                case CahnHilliardControl.ModelType.modelA:
                case CahnHilliardControl.ModelType.modelC:
                default:
                break;
            }

            if(this.Control.CurvatureCorrection) {
                ResidualFields.Cat(curvature_Resi);
            }

            return ResidualFields;
        }


        /// <summary>
        /// DG field instantiation
        /// </summary>
        protected override void CreateFields() {

            base.CreateFields();

            /*
            DummyLevset = new LevelSet(c.Basis, "Levset");
            this.LsTrk = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B" }, DummyLevset);
            DummyLevset.Clear();
            DummyLevset.AccConstant(-1.0);
            this.LsTrk.UpdateTracker(0.0);
            */

            RealLevSet = new LevelSet(c.Basis, "Levset");
            this.RealTracker = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 2, new string[] { "A", "B" }, RealLevSet);
            RealLevSet.Clear();
            RealLevSet.Acc(1.0, c);
            this.RealTracker.UpdateTracker(0.0);
        }


        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            //InitMPI(args);
            //BoSSS.Application.CahnHilliard.Tests.TestProgram.TestCartesian();
            //Assert.True(false);
            _Main(args, false, () => new CahnHilliardMain());
        }

        /// <summary>
        /// Sets the multigrid coloring
        /// </summary>
        protected override void SetInitial(double t) {


            base.SetInitial(t);

            RealLevSet.Clear();
            RealLevSet.Acc(1.0, c);
            this.RealTracker.UpdateTracker(t);

            InitLogFile(this.CurrentSessionInfo.ID);
            WriteLogLine(0, t);

            //VectorField<SinglePhaseField> filtgrad;
            //CurvatureAlgorithmsForLevelSet.CurvatureDriver(
            //                        CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected,
            //                        CurvatureAlgorithmsForLevelSet.FilterConfiguration.Phasefield,
            //                        this.DCurvature, out filtgrad, RealTracker,
            //                        this.Curvature.Basis.Degree * 2,
            //                        c);
        }

        SpatialOperator CHOp;

        BoundaryCondMap<BoundaryType> m_bcMap;

        //Solution.XdgTimestepping.XdgBDFTimestepping m_Timestepper;

        ///// <summary>
        ///// Includes assembly of the matrix.
        ///// 
        ///// </summary>
        ///// <param name="L"></param>
        //protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        //{
        //    using (FuncTrace tr = new FuncTrace())
        //    {

        //        this.CreateTimestepper();
        //    }
        //}

        protected override SpatialOperator GetOperatorInstance(int D) {
            double _D = D;
            //double penalty_base = (c.Basis.Degree + 1) * (c.Basis.Degree + _D) / _D;
            double penalty_factor = base.Control.penalty_poisson;// * penalty_base;

            //BoundaryCondMap<BoundaryType> PoissonBcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "T");

            MultidimensionalArray LengthScales;
            if(this.GridData is GridData) {
                LengthScales = ((GridData)GridData).Cells.cj;
            } else if(this.GridData is AggregationGridData) {
                LengthScales = ((AggregationGridData)GridData).AncestorGrid.Cells.cj;
            } else {
                throw new NotImplementedException();
            }


            m_bcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "c");

            #region variables

            //create Parameter and Variablelists
            string[] paramVar = VariableNames.VelocityVector(D).Cat("c0").Cat(VariableNames.LevelSetGradient(D));
            string[] domainVar = new string[] { "c" };
            string[] codomainVar = new string[] { "Res_c" };

            switch(Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelA:
                break;
                case CahnHilliardControl.ModelType.modelB:
                domainVar = domainVar.Cat("phi");
                codomainVar = codomainVar.Cat("Res_phi");
                break;
                case CahnHilliardControl.ModelType.modelC:
                default:
                break;
            }

            if(this.Control.CurvatureCorrection) {
                domainVar = domainVar.Cat(VariableNames.Curvature);
                codomainVar = codomainVar.Cat("Res_" + VariableNames.Curvature);
            }

            if(this.Control.UseDirectCurvature) {
                paramVar = paramVar.Cat("D" + VariableNames.Curvature);
            }

            #endregion

            CHOp = new SpatialOperator(
                        domainVar,
                        paramVar,
                        codomainVar,
                        (DomainVariableDegrees, ParameterDegrees, CodomainVariableDegrees) => 3 * DomainVariableDegrees[0]//QuadOrderFunc.NonLinear(3)
                        );

            CHOp.ParameterUpdates.Add(CompleteParameterUpdate);

            #region equation components

            // convection term
            if(this.Control.includeConvection == true) {
                CHOp.EquationComponents["Res_c"].Add(
                new c_Flux(D, () => this.Velocity.ToArray(), m_bcMap)
                );
            }

            switch(Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelA:

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                    new c_Source(Control.diff)
                    );
                }

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                        new phi_Diffusion(D, penalty_factor, Control.cahn * Control.diff.Sqrt(), m_bcMap)
                        );
                }

                if(Control.CurvatureCorrection == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                        new phi_CurvatureCorrection(D, Control.cahn * Control.diff.Sqrt(), this.Control.UseDirectCurvature)
                        );

                    if(!this.Control.UseDirectCurvature) {
                        CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                            new curvature_Source(D)
                            );

                        CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                            new curvature_Divergence(D, penalty_factor, 0.001 / this.Control.cahn, LengthScales)
                            );
                    } else {
                        CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                            new curvature_Direct(D)
                            );
                    }
                }

                break;
                case CahnHilliardControl.ModelType.modelB:

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                    new c_Diffusion(D, penalty_factor, Control.diff, Control.lambda, m_bcMap)
                    );
                }

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_phi"].Add(
                        new phi_Diffusion(D, penalty_factor, Control.cahn, m_bcMap)
                        );
                }

                CHOp.EquationComponents["Res_phi"].Add(
                    //new phi_Source(Control.kappa, Control.lambda)
                    new phi_Source(this.Control.includeDiffusion, Control.cahn)
                    );

                if(Control.CurvatureCorrection == true) {
                    CHOp.EquationComponents["Res_phi"].Add(
                        new phi_CurvatureCorrection(D, Control.cahn)
                        );

                        if (!this.Control.UseDirectCurvature) {
                            CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                new curvature_Source(D)
                                );

                            CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                new curvature_Divergence(D, penalty_factor, 0.001 / this.Control.cahn, LengthScales)
                                );
                        } else {
                            CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                                new curvature_Direct(D)
                                );
                        }
                }

                break;

                case CahnHilliardControl.ModelType.modelC:
                throw new NotImplementedException();
                //break;
                
                default:
                throw new ArgumentOutOfRangeException();
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

        //RungeKuttaScheme rksch = null;
        //int bdfOrder = -1000;
        //private void CreateTimestepper()
        //{
        //    switch (this.Control.TimeSteppingScheme)
        //    {
        //        case TimeSteppingScheme.RK_ImplicitEuler:
        //            {
        //                rksch = RungeKuttaScheme.ImplicitEuler;
        //                break;
        //            }
        //        case TimeSteppingScheme.RK_CrankNic:
        //            {
        //                rksch = RungeKuttaScheme.CrankNicolson;
        //                break;
        //            }
        //        case TimeSteppingScheme.CrankNicolson:
        //            {
        //                //do not instantiate rksch, use bdf instead
        //                bdfOrder = -1;
        //                break;
        //            }
        //        case TimeSteppingScheme.ImplicitEuler:
        //            {
        //                //do not instantiate rksch, use bdf instead
        //                bdfOrder = 1;
        //                break;
        //            }
        //        default:
        //            {
        //                if (this.Control.TimeSteppingScheme.ToString().StartsWith("BDF"))
        //                {
        //                    //do not instantiate rksch, use bdf instead
        //                    bdfOrder = Convert.ToInt32(this.Control.TimeSteppingScheme.ToString().Substring(3));
        //                    break;
        //                }
        //                else
        //                    throw new NotImplementedException();
        //            }

        //    }


        //    if (rksch == null)
        //    {
        //        m_Timestepper = new XdgBDFTimestepping(
        //            this.CurrentSolution.Mapping.Fields,
        //            this.CurrentResidual.Mapping.Fields,
        //            base.LsTrk,
        //            true,
        //            DelComputeOperatorMatrix, null, null,
        //            (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? bdfOrder : 1,
        //            LevelSetHandling.None, MassMatrixShapeandDependence.IsTimeDependent, SpatialOperatorType.Nonlinear,
        //            this.MassScale,
        //            this.MgConfig, base.MultigridSequence,
        //            new[] { base.LsTrk.GetSpeciesId("A") }, c.Basis.Degree * 2,
        //            0.0,
        //            false,
        //            this.Control.NonLinearSolver,
        //            this.Control.LinearSolver
        //            );
        //        m_Timestepper.m_ResLogger = base.ResLogger;
        //        m_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
        //        m_Timestepper.m_ResLogger.WriteResidualsToTextFile = false;
        //    }
        //    else
        //    {
        //        throw new NotSupportedException();                
        //    }
        //}

        ///// <summary>
        ///// Computation of operator matrix
        ///// </summary>
        //void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time) {

        //    #region FDJacobian 0
        //    //BlockMsrMatrix OutMtx;
        //    //double[] OutAffine;

        //    //AssembleMatrix(out OutMtx, out OutAffine, CurrentState, OpMtx != null);
        //    //if (OpMtx != null)
        //    //{
        //    //    OpMtx.Clear();
        //    //    OpMtx.Acc(1.0, OutMtx);
        //    //}

        //    //OpAffine.Clear();
        //    //OpAffine.AccV(1.0, OutAffine);
        //    #endregion

        //    #region FDJacobian 1

        //    //create mappings
        //    var codMap = Mapping;
        //    var domMap = Mapping;

        //    DGField[] prms;
        //    prms = ArrayTools.Cat(this.Velocity, c0, gradc0);

        //    if (this.Control.CurvatureCorrection && this.Control.UseDirectCurvature)
        //    {
        //        prms = prms.Cat(this.DCurvature);
        //    }
        //    //prms = ArrayTools.Cat(this.Velocity, c0, gradc0);

        //    if (OpMtx != null)
        //    {
        //        // ++++++++++++++++++++++++++++++++
        //        // create matrix and affine vector:
        //        // ++++++++++++++++++++++++++++++++

        //        IEvaluatorLinear mtxBuilder;

        //        if (this.Control.UseFDJacobian)
        //        {
        //            mtxBuilder = CHOp.GetFDJacobianBuilder(CurrentState, prms, codMap);
        //            mtxBuilder.time = time;
        //            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
        //        }
        //        else
        //        {
        //            var JacParams = JacobiOp.ParameterUpdate;
        //            var TmpParams = JacParams.AllocateParameters(CurrentState, prms);
        //            var map = Mapping;// new CoordinateMapping(CurrentState);

        //            var JacBuilder = JacobiOp.GetMatrixBuilder(map, TmpParams, map);
        //            JacobiOp.ParameterUpdate.PerformUpdate(CurrentState, TmpParams);
        //            JacParams.PerformUpdate(CurrentState, TmpParams);
        //            JacBuilder.ComputeMatrix(OpMtx, OpAffine);
        //        }

        //    }
        //    else
        //    {
        //        // ++++++++++++++++++++++++++++++++
        //        // evaluate the operator
        //        // ++++++++++++++++++++++++++++++++
        //        var eval = CHOp.GetEvaluatorEx(CurrentState, prms, codMap);
        //        eval.time = time;
        //        eval.Evaluate(1.0, 1.0, OpAffine);
        //    }

        //    try
        //    {
        //        if (OpMtx != null)
        //            OpMtx.CheckForNanOrInfM();
        //        OpAffine.CheckForNanOrInfV();
        //    }
        //    catch (ArithmeticException ae)
        //    {
        //        Console.WriteLine("Found NAN");

        //        foreach (DGField f in CurrentState)
        //        {
        //            f.GetExtremalValues(out double min, out double max);
        //            Console.WriteLine($"  Field {f.Identification} extremal values: {min}  ---  {max}");
        //        }

        //        throw ae;
        //    }

        //    #endregion

        //    #region old
        //    //SinglePhaseField Current_c = (SinglePhaseField)(CurrentState[0]);

        //    //c0.Clear();
        //    //c0.Acc(1.0, Current_c);
        //    //gradc0.Clear();
        //    //gradc0.Gradient(1.0, c0);

        //    //if (Control.CurvatureCorrection)
        //    //{
        //    //    VectorField<SinglePhaseField> filtgrad;
        //    //    CurvatureAlgorithms.CurvatureDriver(
        //    //                 SurfaceStressTensor_IsotropicMode.Curvature_Projected,
        //    //                 CurvatureAlgorithms.FilterConfiguration.Phasefield,
        //    //                 this.Curvature, out filtgrad, LsTrk,
        //    //                 this.Curvature.Basis.Degree * 2,
        //    //                 c0);
        //    //}

        //    //var domainvar = ArrayTools.Cat<DGField>(CurrentState);
        //    //var parameter = ArrayTools.Cat<DGField>(this.Velocity, c0, gradc0, Curvature);
        //    //var mb = CHOp.GetMatrixBuilder(Mapping, this.Velocity.ToArray().Cat(c0).Cat(gradc0.ToArray()), Mapping);

        //    //mb.ComputeMatrix(OpMtx, OpAffine);
        //    #endregion

        //}

        //private void AssembleMatrix(out BlockMsrMatrix outMtx, out double[] outAffine, DGField[] CurrentState, bool Linearization)
        //{
        //    // parameters
        //    //============================================================
        //    DGField[] Params;
        //    Params = ArrayTools.Cat(this.Velocity, c0, gradc0, Curvature);

        //    // create mappings
        //    //==========================================================
        //    var codMap = this.CurrentResidual.Mapping;
        //    var domMap = this.CurrentSolution.Mapping;

        //    // provide a linearization of the operator
        //    //===========================================================
        //    if (Linearization)
        //    {
        //        // create matrix and affine vector:
        //        outMtx = new BlockMsrMatrix(codMap, domMap);
        //        outAffine = new double[codMap.LocalLength];

        //        var FDbuilder = CHOp.GetFDJacobianBuilder(domMap, Params, codMap);

        //        FDbuilder.ComputeMatrix(outMtx, outAffine);

        //        // FDJacobian has (Mx +b) as RHS, for unsteady calc. we must subtract Mx for real affine Vector!
        //        outMtx.SpMV(-1.0, new CoordinateVector(CurrentState), 1.0, outAffine);


        //        outMtx.CheckForNanOrInfM();
        //        outAffine.CheckForNanOrInfV();
        //    }
        //    else
        //    {

        //        // explicit evaluation of the operator
        //        //========================================================
        //        outMtx = null;
        //        outAffine = new double[codMap.LocalLength];
        //        var eval = CHOp.GetEvaluatorEx(CurrentState, Params, codMap);
        //        this.UpdateParameter();

        //        eval.Evaluate(1.0, 1.0, outAffine);

        //    }

        //}

        int reinit = 0;

        [InstantiateFromControlFile("cDist", "c", IOListOption.Always)]
        SinglePhaseField cDist;

        private void ComputeDistanceField() {
            GridData GridDat = (GridData)(c.GridDat);

            // compute and project 
            // step one calculate distance field phiDist = 0.5 * log(Max(1+phi, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn_old
            // step two project the new phasefield phiNew = tanh(phiDist/(sqrt(2) * Cahn_new))
            // here done in one step, with default quadscheme
            // ===================
            cDist.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Phi
                    // -----------------------------
                    c.Evaluate(j0, Len, NS, result);

                    // compute the pointwise values of the new level set
                    // -----------------------------

                    result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * this.Control.cahn);
                }
            );

            reinit++;
            //PlotCurrentState(0.0, reinit);
        }

        [InstantiateFromControlFile("Correction", "c", IOListOption.Always)]
        SinglePhaseField Correction;

        private void ComputeCorrection() {
            GridData GridDat = (GridData)(c.GridDat);

            SinglePhaseField[,] cHess;
            SinglePhaseField[] cGrad = this.gradc0.CloneAs().ToArray();

            // Compute Hessian
            CurvatureAlgorithmsForLevelSet.ComputeHessian(this.c, out cHess);

            // buffers:
            MultidimensionalArray C = new MultidimensionalArray(2);
            MultidimensionalArray GradC = new MultidimensionalArray(3);
            MultidimensionalArray HessC = new MultidimensionalArray(4);

            MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
            MultidimensionalArray NormGrad = new MultidimensionalArray(3);
            MultidimensionalArray Q = new MultidimensionalArray(3);

            // compute and project 
            // correction = n * cHess * n
            // ===================
            Correction.Clear();
            Correction.ProjectField(1.0,
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes
                    int D = c.GridDat.SpatialDimension;

                    // alloc buffers
                    // -------------

                    if(C.GetLength(0) != Len || C.GetLength(1) != K) {
                        C.Allocate(Len, K);
                        GradC.Allocate(Len, K, D);
                        HessC.Allocate(Len, K, D, D);
                        ooNormGrad.Allocate(Len, K);
                        NormGrad.Allocate(Len, K, D);
                        Q.Allocate(Len, K, D);
                    } else {
                        C.Clear();
                        GradC.Clear();
                        HessC.Clear();
                        ooNormGrad.Clear();
                        NormGrad.Clear();
                        Q.Clear();
                    }

                    // evaluate Gradient and Hessian
                    // -----------------------------

                    {
                        for(int d = 0; d < D; d++)
                            cGrad[d].Evaluate(j0, Len, NS, GradC.ExtractSubArrayShallow(-1, -1, d));
                    }
                    {
                        for(int d1 = 0; d1 < D; d1++)
                            for(int d2 = 0; d2 < D; d2++)
                                cHess[d1, d2].Evaluate(j0, Len, NS, HessC.ExtractSubArrayShallow(-1, -1, d1, d2));
                    }

                    // compute the monstrous formula
                    // -----------------------------

                    // norm of Gradient:
                    for(int d = 0; d < D; d++) {
                        var GradC_d = GradC.ExtractSubArrayShallow(-1, -1, d);
                        ooNormGrad.Multiply(1.0, GradC_d, GradC_d, 1.0, "ik", "ik", "ik");
                    }
                    ooNormGrad.ApplyAll(x => x > 1e-6 ? 1.0 / Math.Sqrt(x) : 0.0); // prohibit the norm to be very small

                    NormGrad.Multiply(1.0, GradC, ooNormGrad, 0.0, "ikj", "ikj", "ik");

                    // normal part of cHess: Q = n * cHess
                    Q.Multiply(1.0, NormGrad, HessC, 0.0, "ikj", "ikl", "iklj");

                    // result = Q * n
                    result.Multiply(1.0, Q, NormGrad, 0.0, "ik", "ikj", "ikj");

                    // result = laplace/(|GradC|)^2
                    //result.Multiply(1.0, Laplace, ooNormGrad, 0.0, "ik", "ik", "ik");

                }, (new CellQuadratureScheme()).SaveCompile(GridDat, Correction.Basis.Degree * 2 + 2));
            //PlotCurrentState(0.0, reinit);
        }

        void CompleteParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            UpdateParameter();
        }

        void UpdateParameter() {
            //c0.Clear();
            //c0.Acc(1.0, DomainVar.First());
            //gradc0.Clear();
            //gradc0.Gradient(1.0, c0);
            int D = this.GridData.SpatialDimension;
            this.c0.Clear();
            this.c0.Acc(1.0, this.c);

            ComputeDistanceField();

            gradc0.Clear();
            gradc0.Gradient(1.0, this.c0);

            CurvatureAlgorithmsForLevelSet.ComputeHessian(this.c, out hessc0);
            if(Control.CurvatureCorrection) {
                VectorField<SinglePhaseField> filtgrad;
                CurvatureAlgorithmsForLevelSet.CurvatureDriver(
                                CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                                CurvatureAlgorithmsForLevelSet.FilterConfiguration.Phasefield,
                                this.DCurvature, out filtgrad, RealTracker,
                                this.c.Basis.Degree * 2,
                                this.c/*this.c0*/);

                ComputeCorrection();
            }

        }

        CoordinateVector m_CurrentSolution = null;

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution {
            get {
                if(m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(this.c);

                    switch(Control.ModTyp) {
                        case CahnHilliardControl.ModelType.modelA:
                        break;
                        case CahnHilliardControl.ModelType.modelB:
                        m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.phi));
                        break;
                        case CahnHilliardControl.ModelType.modelC:
                        default:
                        break;
                    }

                    if(this.Control.CurvatureCorrection) {
                        m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.Curvature));
                    }
                }
                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual {
            get {
                if(m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(this.c_Resi);

                    switch(Control.ModTyp) {
                        case CahnHilliardControl.ModelType.modelA:
                        break;
                        case CahnHilliardControl.ModelType.modelB:
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.phi_Resi));
                        break;
                        case CahnHilliardControl.ModelType.modelC:
                        default:
                        break;
                    }

                    if(this.Control.CurvatureCorrection) {
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.curvature_Resi));
                    }
                }
                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                double[] scale = new double[1];

                switch(Control.ModTyp) {
                    case CahnHilliardControl.ModelType.modelA:
                    break;
                    case CahnHilliardControl.ModelType.modelB:
                    scale = new double[2];
                    break;
                    case CahnHilliardControl.ModelType.modelC:
                    throw new NotImplementedException();
                    break;
                    default:
                    throw new ArgumentOutOfRangeException();
                    break;
                }

                if(this.Control.CurvatureCorrection) {
                    scale = new double[scale.Length + 1];
                }

                scale[0] = 1.0;

                R.Add(this.LsTrk.GetSpeciesId("A"), scale);

                return R;
            }
        }


        /// <summary>
        /// control of mesh adaptation
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            if(this.Control.AdaptiveMeshRefinement && TimestepNo > 1) {


                long oldJ = this.GridData.CellPartitioning.TotalLength;

                //double LocNormPow2 = this.ResiualKP1.CoordinateVector.L2NormPow2(); // norm of residual on this processor
                //double TotNormPow2 = LocNormPow2.MPISum(); //                          norm of residual over all processors
                //double MeanNormPow2PerCell = TotNormPow2 / oldJ; //                    mean norm per cell


                int MyLevelIndicator(int j, int CurrentLevel) {
                    /*
                    double CellNorm = this.ResiualKP1.Coordinates.GetRow(j).L2NormPow2();


                    if (j == 0)
                        CurrentLevel = CurrentLevel + 1;

                    if (CellNorm > MeanNormPow2PerCell * 1.1)
                        return CurrentLevel + 1;
                    else
                        return CurrentLevel;
                    */

                    return 0;
                }


                GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData, null);
                bool AnyChange = gridRefinementController.ComputeGridChange(MyLevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if(AnyChange) {
                    int[] glb = (new int[] {
                        CellsToRefineList.Count,
                        Coarsening.Sum(L => L.Length),
                        //0, 0
                    }).MPISum();

                    NoOfCellsToRefine = glb[0];
                    NoOfCellsToCoarsen = glb[1];
                }
                //*/


                // Update Grid
                // ===========

                if(AnyChange) {


                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                    newGrid = ((GridData)(this.GridData)).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }
            } else {

                newGrid = null;
                old2NewGrid = null;
            }
        }


        /// <summary>
        /// Single run of the solver
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using(new FuncTrace()) {

                dt = this.Control.dtFixed;
                Console.WriteLine("Starting with timestep #{0}, t = {1}, dt = {2} ...", TimestepNo, phystime, dt);

                // Perform timestep
                // ================

                //this.m_Timestepper.m_ResLogger = base.ResLogger;
                //this.m_Timestepper.m_ResLogger.WriteResidualsToTextFile = false;
                //switch (this.Control.ModTyp)
                //{
                //    case CahnHilliardControl.ModelType.modelA:
                //        this.m_Timestepper.m_ResidualNames = new string[] { "Res_c", "Res_" + VariableNames.Curvature };
                //        break;
                //    case CahnHilliardControl.ModelType.modelB:
                //        this.m_Timestepper.m_ResidualNames = new string[] { "Res_c", "Res_phi", "Res_" + VariableNames.Curvature };
                //        break;
                //    case CahnHilliardControl.ModelType.modelC:
                //    default:
                //        throw new NotImplementedException();
                //}


                //c0.Clear();
                //c0.Acc(1.0, c);
                //gradc0.Clear();
                //gradc0.Gradient(1.0, c0);
                //VelocityGradX.Clear();
                //VelocityGradY.Clear();
                //VelocityGradX.Gradient(1.0, this.Velocity[0]);
                //VelocityGradY.Gradient(1.0, this.Velocity[1]);

                //if (Control.CurvatureCorrection)
                //{
                //    VectorField<SinglePhaseField> filtgrad;
                //    CurvatureAlgorithmsForLevelSet.CurvatureDriver(
                //                    CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //                    CurvatureAlgorithmsForLevelSet.FilterConfiguration.Phasefield,
                //                    this.DCurvature, out filtgrad, RealTracker,
                //                    this.c.Basis.Degree * 2 + 2,
                //                    this.c);
                //}

                var Qnts_old = ComputeBenchmarkQuantities();

                base.Timestepping.Solve(phystime, dt);

                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                this.RealTracker.UpdateTracker(0.0);

                // algebraic correction
                switch(this.Control.CorrectionType) {
                    case CahnHilliardControl.Correction.Concentration:
                    ConservativityCorrection(Qnts_old);
                    break;
                    case CahnHilliardControl.Correction.Mass:
                    MassCorrection(Qnts_old);
                    break;
                    case CahnHilliardControl.Correction.None:
                    default:
                    break;
                }

                WriteLogLine(TimestepNo, phystime + dt);

                // return
                // ======

                Console.WriteLine("done with timestep #{0}.", TimestepNo);
                return dt;
            }
        }

        /// shift the concentration field to account for mass loss or gain in the tracked phase
        private void MassCorrection(double[] Qnts_old) {
            double[] Qnts = ComputeBenchmarkQuantities();
            double mass = Qnts[0];
            double surface = Qnts[1];
            double shape = Qnts[2];
            double massDiff = Qnts_old[0] - mass;

            // we assume the current phasefield is close to the equilibrium tangenshyperbolicus form
            SinglePhaseField cNew = new SinglePhaseField(phi.Basis);
            GridData GridDat = (GridData)(phi.GridDat);
            double mass_uc = mass;

            int i = 0;
            while(massDiff.Abs() > 1e-6 && i < 10) {
                // calculated for a cone, one could include the shape e.g. by using the circularity
                // correction guess
                double correction = surface / (4 * mass * Math.PI) * massDiff;

                // take the correction guess and calculate a forward difference to approximate the derivative
                cNew.ProjectField(
                    (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                        Debug.Assert(result.Dimension == 2);
                        Debug.Assert(Len == result.GetLength(0));
                        int K = result.GetLength(1); // number of nodes

                        // evaluate Phi
                        // -----------------------------
                        c.Evaluate(j0, Len, NS, result);

                        // compute the pointwise values of the new level set
                        // -----------------------------

                        result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * Math.Sqrt(2) * this.Control.cahn);
                        result.ApplyAll(x => Math.Tanh((x + correction) / (Math.Sqrt(2) * this.Control.cahn)));
                    }
                );

                // update LsTracker
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, cNew);
                this.RealTracker.UpdateTracker(0.0);

                Qnts = ComputeBenchmarkQuantities();
                mass = Qnts[0];

                correction = -(massDiff) / ((Qnts_old[0] - mass - massDiff) / (correction));

                // compute and project 
                // step one calculate distance field phiDist = 0.5 * log(Max(1+c, eps)/Max(1-c, eps)) * sqrt(2) * Cahn
                // step two project the new phasefield phiNew = tanh((cDist + correction)/(sqrt(2) * Cahn))
                // ===================
                cNew.ProjectField(
                    (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                        Debug.Assert(result.Dimension == 2);
                        Debug.Assert(Len == result.GetLength(0));
                        int K = result.GetLength(1); // number of nodes

                        // evaluate Phi
                        // -----------------------------
                        c.Evaluate(j0, Len, NS, result);

                        // compute the pointwise values of the new level set
                        // -----------------------------

                        result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * Math.Sqrt(2) * this.Control.cahn);
                        result.ApplyAll(x => Math.Tanh((x + correction) / (Math.Sqrt(2) * this.Control.cahn)));
                    }
                );

                // update field
                c.Clear();
                c.Acc(1.0, cNew);

                // update LsTracker
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                this.RealTracker.UpdateTracker(0.0);

                Qnts = ComputeBenchmarkQuantities();
                mass = Qnts[0];
                surface = Qnts[1];
                shape = Qnts[2];

                massDiff = Qnts_old[0] - mass;
                i++;
            }

            Qnts = ComputeBenchmarkQuantities();

            Console.WriteLine($"Performed Mass Correction in {i} iteratins: \n" +
                $"\told mass:           {Qnts_old[0]:N4}\n" +
                $"\tuncorrected mass:   {mass_uc:N4}\n" +
                $"\tcorrected mass:     {Qnts[0]:N4}");
        }

        private void ConservativityCorrection(double[] Qnts_old) {
            throw new NotImplementedException();
        }

        MultigridOperator.ChangeOfBasisConfig[][] MgConfig {
            get {
                int p = this.c.Basis.Degree;
                int NoOfLevels = this.MultigridSequence.Length;
                var config = new MultigridOperator.ChangeOfBasisConfig[NoOfLevels][];
                int m = 0;

                for(int iLevel = 0; iLevel < NoOfLevels; iLevel++) {

                    switch(Control.ModTyp) {
                        case CahnHilliardControl.ModelType.modelA:
                        m = 1;
                        config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m];
                        break;
                        case CahnHilliardControl.ModelType.modelB:
                        m = 2;
                        config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m];
                        break;
                        case CahnHilliardControl.ModelType.modelC:
                        throw new NotImplementedException();
                        break;
                        default:
                        throw new ArgumentOutOfRangeException();
                        break;
                    }

                    if(this.Control.CurvatureCorrection) {
                        config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m + 1];
                    }

                    config[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        VarIndex = new int[] { 0 },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                        DegreeS = new int[] { Math.Max(1, p - iLevel) }
                    };

                    switch(Control.ModTyp) {
                        case CahnHilliardControl.ModelType.modelA:
                        break;
                        case CahnHilliardControl.ModelType.modelB:
                        config[iLevel][1] = new MultigridOperator.ChangeOfBasisConfig() {
                            VarIndex = new int[] { 1 },
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            DegreeS = new int[] { Math.Max(1, p - iLevel) }
                        };

                        break;
                        case CahnHilliardControl.ModelType.modelC:
                        throw new NotImplementedException();
                        break;
                        default:
                        throw new ArgumentOutOfRangeException();
                        break;
                    }

                    if(this.Control.CurvatureCorrection) {
                        config[iLevel][m] = new MultigridOperator.ChangeOfBasisConfig() {
                            VarIndex = new int[] { m },
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            DegreeS = new int[] { Math.Max(0, this.Curvature.Basis.Degree - iLevel) }
                        };
                    }
                }

                return config;
            }

        }

        public double[] ComputeBenchmarkQuantities() {

            int order = 0;
            if(RealTracker.GetCachedOrders().Count > 0) {
                order = RealTracker.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = RealTracker.GetXDGSpaceMetrics(RealTracker.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // area of bubble
            double area = 0.0;
            SpeciesId spcId = RealTracker.SpeciesIdS[1];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, RealTracker.GridDat,
                vqs.Compile(RealTracker.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            area = area.MPISum();

            // surface
            double surface = 0.0;
            //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId, 0);
            CellQuadrature.GetQuadrature(new int[] { 1 }, RealTracker.GridDat,
                surfElemVol.Compile(RealTracker.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            surface = surface.MPISum();

            // circularity
            double diamtr_c = Math.Sqrt(4 * area / Math.PI);
            double perimtr_b = surface;

            double circ = Math.PI * diamtr_c / perimtr_b;

            // total concentration
            double concentration = 0.0;
            var tqs = new CellQuadratureScheme();
            CellQuadrature.GetQuadrature(new int[] { 1 }, c.GridDat,
                tqs.Compile(c.GridDat, c.Basis.Degree * 2 + 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    c.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        concentration += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            concentration = concentration.MPISum();

            // total mixing energy
            double energy = 0.0;
            var eqs = new CellQuadratureScheme();

            int D = c.GridDat.SpatialDimension;
            SinglePhaseField[] cGrad = new SinglePhaseField[D];
            for(int d = 0; d < D; d++) {
                cGrad[d] = new SinglePhaseField(c.Basis, string.Format("G_{0}", d));
                cGrad[d].Derivative(1.0, c, d);
            }

            MultidimensionalArray Phi = new MultidimensionalArray(2);
            MultidimensionalArray GradPhi = new MultidimensionalArray(3);
            MultidimensionalArray NormGrad = new MultidimensionalArray(2);

            CellQuadrature.GetQuadrature(new int[] { 1 }, c.GridDat,
                eqs.Compile(c.GridDat, c.Basis.Degree * 2 + 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    int K = EvalResult.GetLength(1);
                    // alloc buffers
                    // -------------

                    if(Phi.GetLength(0) != Length || Phi.GetLength(1) != K) {
                        Phi.Allocate(Length, K);
                        GradPhi.Allocate(Length, K, D);
                        NormGrad.Allocate(Length, K);
                    } else {
                        Phi.Clear();
                        GradPhi.Clear();
                        NormGrad.Clear();
                    }

                    // chemical potential
                    c.Evaluate(i0, Length, QR.Nodes, Phi.ExtractSubArrayShallow(-1, -1));
                    Phi.ApplyAll(x => 0.25 / (this.Control.cahn.Pow2()) * (x.Pow2() - 1.0).Pow2());

                    for(int d = 0; d < D; d++) {
                        cGrad[d].Evaluate(i0, Length, QR.Nodes, GradPhi.ExtractSubArrayShallow(-1, -1, d));
                    }

                    // free surface energy
                    for(int d = 0; d < D; d++) {
                        var GradPhi_d = GradPhi.ExtractSubArrayShallow(-1, -1, d);
                        NormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                    }
                    NormGrad.ApplyAll(x => 0.5 * x);

                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, Phi);
                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, NormGrad);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        energy += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            energy = energy.MPISum();

            return new double[] { area, surface, circ, concentration, energy };
        }

        /// <summary>
        /// testcase specific LogFile
        /// </summary>
        TextWriter Log;

        /// <summary>
        /// initializes the format of the Log File
        /// </summary>
        /// <param name="sessionID"></param>
        public void InitLogFile(Guid sessionID) {
            if(this.Control.savetodb) {
                Log = base.DatabaseDriver.FsDriver.GetNewLog("Phasefield_Quantities.txt", sessionID);
            } else {
                Log = new StreamWriter("Phasefield_Quantities_MPI" + this.MPIRank + ".txt");
            }

            return;
        }

        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {
            double[] BQnts = ComputeBenchmarkQuantities();

            string line = String.Format($"{TimestepNo}\t{phystime}\t{BQnts[0]}\t{BQnts[1]}\t{BQnts[2]}\t{BQnts[3]}\t{BQnts[4]}");
            Log.WriteLine(line);
            Log.Flush();

            return;
        }

        /// <summary>
        /// Shutdown function
        /// </summary>
        protected override void Bye() {
            object SolL2err;
            if(this.QueryHandler.QueryResults.TryGetValue("SolL2err", out SolL2err)) {
                Console.WriteLine("Value of Query 'SolL2err' " + SolL2err.ToString());
            } else {
                Console.WriteLine("query 'SolL2err' not found.");
            }
        }

        /// <summary>
        /// default plotting
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            string caseStr = "";
            if(base.Control.Paramstudy_CaseIdentification != null) {
                var pstudy_case = base.Control.Paramstudy_CaseIdentification.FirstOrDefault(tt => tt.Item1 == "pstudy_case");
                if(pstudy_case != null) {
                    caseStr = "." + pstudy_case.Item2;
                }
            }

            DGField[] Fields = new DGField[0];
            Fields = Fields.Cat(this.cex, this.c, this.phi, this.Velocity, this.gradc0, this.Curvature, this.DCurvature, this.c_Resi, this.phi_Resi, this.curvature_Resi, this.cDist, this.Correction);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "CahnHilliard-" + timestepNo + caseStr, phystime, superSampling);
        }

    }



    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'c'
    /// </summary>
    public class c_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace, IVolumeForm, ISupportsJacobianComponent {

        public c_Diffusion(int D, double penalty_const, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, "phi") // note: in the equation for 'c', we have the Laplacian of 'phi'
        {
            m_D = D;
            m_diff = __diff;
            m_lambda = __lambda;
            m_boundaryCondMap = __boundaryCondMap;
        }

        double m_diff;
        double min = double.MaxValue;
        double max = 0.0;
        double m_lambda;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => new[] { "c0" }.Cat(VariableNames.VelocityVector(m_D)).Cat(VariableNames.LevelSetGradient(m_D));



        protected override double g_Diri(ref CommonParamsBnd inp) {
            double UxN = (new Vector(inp.Parameters_IN, 1, inp.D))*inp.Normal;

            double v;
            if(UxN >= 0) {
                // outflow
                v = 1.0;
            } else {
                // inflow
                v = 0.0;
            }

            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgeType) {
                case BoundaryType.Wall:
                // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                return false;

                case BoundaryType.Flow:
                // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                return true;

                default:
                throw new NotImplementedException();
            }
        }

        public override double Nu(double[] x, double[] p, int jCell) {
            double n = 0.0;

            for(int d = 0; d < m_D; d++) {
                // n += p[1 + m_D + d].Pow2();
                n += p[1 + m_D + d]*p[1 + m_D + d];
            }

            //double D = 0.0;
            //if (n.Sqrt() < 1.0 / (Math.Sqrt(2) * m_diff))
            //{
            //    D = 0.027 / (1 + n * m_diff.Pow2());
            //}
            //else
            //{
            // D = Math.Exp(-Math.Pow(2.0, 3.0) * n * Math.Pow(m_diff, 2.0));
            //}

            //if (D < min || D > max) Console.WriteLine("min/max: " + D);
            //max = Math.Max(D, max);
            //min = Math.Min(D, min);

            if(m_lambda == 0.0) {
                return -m_diff;
                //return -D;//-m_diff * U;//-gradUxN.L2Norm() * m_diff;
            } else if(m_lambda > 0.0 && m_lambda <= 1.0) {
                // return -m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0), 0.0);
                double ret = -m_diff * 1 - m_lambda * p[0]*p[0];
                if (ret > 0){
                    return -m_diff * 1 - m_lambda * p[0]*p[0];
                } else {
                    return 0.0;
                }
                //return -(Math.Abs(p[0]) - Math.Min(p[0].Pow2(),1.0));//-1.0 * Math.Abs(p[0]);//-m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0), 0.0);
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if(this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                if(g_D == 0) {
                    for(int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                } else {
                    for(int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency                        
                    }

                    Acc *= 0.0;            //switch, 0.0 seems stable, 1.0 explodes
                    Acc *= this.m_alpha;

                }
            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }

        // linear term return self
        public override IEquationComponent[] GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { this };
        }
    }

    /// <summary>
    /// Transport flux for Cahn-Hilliard
    /// </summary>
    public class c_Flux : IVolumeForm, IEdgeForm, ISupportsJacobianComponent, IParameterHandling {
        public c_Flux(int D, Func<DGField[]> GetVelVector, BoundaryCondMap<BoundaryType> __boundaryCondMap) {
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap?.bndFunction["c"];
            m_GetVelVector = GetVelVector;
        }

        Func<DGField[]> m_GetVelVector;

        int m_D;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ArgumentOrdering => new string[] { "c" };

        public IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            if(Arguments.Length != 1)
                throw new ArgumentException();
            if(Parameters.Length != m_D)
                throw new ArgumentException();

            // Velocity Vector is provided externally -> no update required
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            if(Arguments.Length != 1)
                throw new ArgumentException();

            var Vel = m_GetVelVector();
            if (Vel == null || Vel.Length != m_D)
                throw new ArgumentException();

            return Vel;
        }


        public TermActivationFlags BoundaryEdgeTerms => InnerEdgeTerms | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double UxN = 0;
            for(int d = 0; d < m_D; d++) {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double c;
            if(UxN >= 0) {
                c = _uIN[0];
            } else {
                //c =m_bndFunc[inp.EdgeTag](inp.X, inp.time);
                c = -1.0;
            }

            return c * UxN * _vIN;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double UxN = 0;
            for(int d = 0; d < m_D; d++) {
                UxN += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];
            }

            double c;
            if(UxN >= 0) {
                c = _uIN[0];
            } else {
                c = _uOUT[0];
            }

            return c * UxN * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            double c = U[0];
            for(int d = 0; d < m_D; d++) {
                acc += c * cpv.Parameters[d] * GradV[d];
            }

            return -acc;
        }

        // linear term return self
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { this };
        }
    }



    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'c'
    /// </summary>
    public class phi_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace, ISupportsJacobianComponent {

        public phi_Diffusion(int D, double penalty_const, double __cahn, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, "c") // note: in the equation for 'phi', we have the Laplacian of 'c'
        {
            // m_cahn = __cahn * __cahn;
            m_cahn = __cahn;
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap?.bndFunction["c"];
        }

        double m_cahn;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);
        Func<double[], double, double>[] m_bndFunc;


        protected override double g_Diri(ref CommonParamsBnd inp) {
            double UxN = (new Vector(inp.Parameters_IN, 0, inp.D))*inp.Normal;

            double v;
            if(UxN >= 0) {
                // outflow
                v = 0.0;
            } else {
                // inflow
                v = m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgeType) {
                case BoundaryType.Wall:
                return false;

                case BoundaryType.Flow:
                return true;

                default:
                throw new NotImplementedException();
            }
        }

        public override double Nu(double[] x, double[] p, int jCell) {
            return -m_cahn;
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if(this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                if(g_D != 0) {
                    for(int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                } else {
                    for(int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency    
                    }

                    Acc *= 0.0;                 //switch, 0.0 seems stable, 1.0 explodes
                    Acc *= this.m_alpha;

                }
            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }

        /*
        // linear term return self
        public override IEquationComponent[] GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { this };
        }
        */
    }


    /// <summary>
    /// nonlinear source term in the 'phi'-equation
    /// </summary>
    public class phi_Source : IVolumeForm, ISupportsJacobianComponent, IParameterHandling {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public phi_Source(bool __inc, double _cahn = 0.0) {
            m_inc = __inc;
            m_scale = _cahn;
        }


        //double m_lambda;
        //double m_epsilon;
        bool m_inc;
        double m_scale;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "phi", "c" };

        public IList<string> ParameterOrdering => new[] { "c0" };

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var c = Arguments[1];
            var c0 = Parameters[0];

            if(object.ReferenceEquals(c0, c))
                return;

            c0.Clear();
            c0.Acc(1.0, c);
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[1];
            return new[] { c };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double phi = U[0];
            double c = U[1];
            double c0 = cpv.Parameters[0];

            double Acc = 0;
            if(m_inc == false) {
                Acc += -phi;
            } else {
                Acc += -phi;

                // Acc += (3.0 * c0 * c0 - 1.0) * c - 2 * Math.Pow(c0, 3.0); // linearized around c0 (Taylor expansion)
                // Acc += (3.0 * c0 * c0 - 1.0) * c - 2 * c0*c0*c0; // linearized around c0 (Taylor expansion)
                Acc += c.Pow(3.0) - c; // for newton with jacobian no linearization is needed
            }


            return Acc * V;
        }

        // already linearized term return self
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { new jacobi_phi_Source(m_inc, m_scale) };
        }
    }

    /// <summary>
    /// nonlinear source term in the 'phi'-equation
    /// </summary>
    class jacobi_phi_Source : IVolumeForm, IParameterHandling {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public jacobi_phi_Source(bool __inc, double _cahn = 0.0) {
            m_inc = __inc;
            m_scale = _cahn;
        }


        //double m_lambda;
        //double m_epsilon;
        bool m_inc;
        double m_scale;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "phi", "c" };

        public IList<string> ParameterOrdering => new[] { "c0" };

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var c = Arguments[1];
            var c0 = Parameters[0];

            if(object.ReferenceEquals(c0, c))
                return;

            c0.Clear();
            c0.Acc(1.0, c);
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[1];
            return new[] { c };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double phi = U[0];
            double c = U[1];
            double c0 = cpv.Parameters[0];

            double Acc = 0;
            if(m_inc == false) {
                Acc += -phi;
            } else {
                Acc += -phi;

                //Acc += (m_lambda / m_epsilon) * (c0 * c0 - 1) * c; // linearized around c0
                //Acc += ((2*c0 - 1) * (2*c0 - 1) - 1) * (c - 0.5); // linearized around c0, 0<c<1
                // Acc += (c0 * c0 - 1) * c; // linearized around c0
                // Acc += 3 * c0.Pow2() * c - c; // linearized around c0 (Taylor expansion)
                Acc += 3 * c0*c0 * c - c; // linearized around c0 (Taylor expansion)
                //Acc += c.Pow(3) - c; // for newton with jacobian no linearization is needed
            }


            return Acc * V;
        }
    }

    /// <summary>
    /// linear "source" term of phi in c equation in Model A
    /// </summary>
    public class c_Source : IVolumeForm, ISupportsJacobianComponent, IParameterHandling {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public c_Source(double _diff = 0.0) {
            m_diff = _diff;
        }



        //double m_lambda;
        //double m_epsilon;
        double m_diff;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "c" };
        public IList<string> ParameterOrdering => new[] { "c0" };

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var c = Arguments[0];
            var c0 = Parameters[0];

            if(object.ReferenceEquals(c0, c))
                return;

            c0.Clear();
            c0.Acc(1.0, c);
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[0];
            return new[] { c };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double c = U[0];
            double c0 = cpv.Parameters[0];
            // values seem shifted without this offset hack

            double Acc = 0;

            //Acc += (3.0 * c0 * c0 - 1.0) * c - 2 * Math.Pow(c0, 3.0); // linearized around (Taylor expansion)
            Acc += c.Pow(3.0) - c; // without linearization

            return m_diff * Acc * V;
        }

        // already linearized term return self
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { new jacobi_c_Source(m_diff) };
        }
    }

    /// <summary>
    /// linear "source" term of phi in c equation in Model A
    /// </summary>
    class jacobi_c_Source : IVolumeForm, IParameterHandling {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public jacobi_c_Source(double _diff = 0.0) {
            m_diff = _diff;
        }



        //double m_lambda;
        //double m_epsilon;
        double m_diff;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "c" };
        public IList<string> ParameterOrdering => new[] { "c0" };

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var c = Arguments[0];
            var c0 = Parameters[0];

            if(object.ReferenceEquals(c0, c))
                return;

            c0.Clear();
            c0.Acc(1.0, c);
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[0];
            return new[] { c };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double c = U[0];
            double c0 = cpv.Parameters[0];
            // values seem shifted without this offset hack

            double Acc = 0;

            Acc += 3 * c0.Pow2() * c - c; // linearized around c0 (Taylor expansion)

            return m_diff * Acc * V;
        }
    }

    /// <summary>
    /// Correction term to counter along the interface diffusion for Model A
    /// </summary>
    class phi_CurvatureCorrection : IVolumeForm, ISupportsJacobianComponent {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public phi_CurvatureCorrection(int D, double _cahn = 0.0, bool direct = false) {
            // m_cahn = _cahn.Pow2();
            m_cahn = _cahn;
            m_D = D;
            m_direct = direct;
        }



        //double m_lambda;
        //double m_epsilon;
        int m_D;
        double m_cahn;
        bool m_direct;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "c", VariableNames.Curvature };
        //public IList<string> ArgumentOrdering => new[] { "c"};
        //public IList<string> ParameterOrdering => new[] { "D" + VariableNames.Curvature }.Cat(VariableNames.LevelSetGradient(m_D));
        public IList<string> ParameterOrdering => null;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0.0;
            double[] grad = new double[m_D];
            for(int d = 0; d < m_D; d++) {
                grad[d] = GradU[0, d];
                //Acc += cpv.Parameters[1+d].Pow2();
            }
            //Acc = Acc.Sqrt();
            Acc += grad.L2Norm();

            Acc *= m_cahn * U[1];


            // sign minus should be correct, plus produces more sensual results, (sign of curvature?)
            return Acc * V;
        }

        // Use VolumeFormDifferentiator (FD-like)
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            var JacVol = new VolumeFormDifferentiator(this, m_D);
            return new IEquationComponent[] { JacVol };
            //return new IEquationComponent[] { new jacobi_phi_CurvatureCorrection(m_D,m_cahn,m_direct) };
        }
    }

    /// <summary>
    /// Correction term to counter diffusion along the interface for Model A
    /// not correct, do not use
    /// </summary>
    class jacobi_phi_CurvatureCorrection : IVolumeForm {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public jacobi_phi_CurvatureCorrection(int D, double _cahn = 0.0, bool direct = false) {
            m_cahn = _cahn;
            m_D = D;
            m_direct = direct;
        }



        //double m_lambda;
        //double m_epsilon;
        int m_D;
        double m_cahn;
        bool m_direct;

        public TermActivationFlags VolTerms => TermActivationFlags.GradUxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "c", VariableNames.Curvature };
        public IList<string> ParameterOrdering => new[] { "D" + VariableNames.Curvature }.Cat(VariableNames.LevelSetGradient(m_D));

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0.0;
            double[] grad = new double[m_D];
            for(int d = 0; d < m_D; d++) {
                grad[d] = GradU[0, d];
                //Acc += cpv.Parameters[1+d].Pow2();
            }
            //Acc = Acc.Sqrt();
            Acc = grad.L2Norm();
            if(m_direct)
                Acc *= m_cahn * cpv.Parameters[0];
            else
                Acc *= m_cahn * U[1];


            // sign minus should be correct, plus produces more sensual results, (sign of curvature?)
            return Acc * V;
        }
    }

    class curvature_Direct : IEquationComponent, IVolumeForm, ISupportsJacobianComponent {

        public curvature_Direct(int _D) {
            m_D = _D;
        }
        int m_D;

        public IList<string> ArgumentOrdering => new[] { VariableNames.Curvature };

        public IList<string> ParameterOrdering => new[] { "D" + VariableNames.Curvature };

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[1];
            return new[] { c };
        }

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0.0;

            Acc += 1e-8 * (cpv.Parameters[0] - U[0]) * V;

            return Acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }

    }

    class curvature_Source : IEquationComponent, IVolumeForm, ISupportsJacobianComponent {

        public curvature_Source(int _D) {
            m_D = _D;
        }
        int m_D;

        public IList<string> ArgumentOrdering => new[] { VariableNames.Curvature };

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        //public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV;

        //public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0.0;

            Acc += -U[0] * V;

            return Acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }

    }

    class curvature_Divergence : IVolumeForm, IEdgeForm, IEquationComponent, ISupportsJacobianComponent {
        public curvature_Divergence(int D, double penalty, double limit, MultidimensionalArray InverseLengthScales) {
            m_D = D;
            m_limit = limit;
            m_penalty = penalty;
            this.InverseLengthScales = InverseLengthScales;
        }

        int m_D;
        double m_limit;
        double m_penalty;

        public CellMask m_cells;

        public IList<string> ArgumentOrdering => new[] { "c", VariableNames.Curvature };

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms => TermActivationFlags.GradUxGradV;

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;

        /// <summary>
        /// a little switch...
        /// </summary>
        protected double m_alpha = 1.0;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;

            double[] grad = new double[m_D];
            for(int d = 0; d < cpv.D; d++) {
                grad[d] = GradU[0, d];
                acc -= GradU[0, d] * GradV[d] * this.m_alpha;
            }
            double norm = Math.Max(grad.L2Norm(), m_limit);
            acc *= 1 / norm;

            return acc;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double Acc = 0.0;

            double pnlty = this.GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);

            //double[] gradIN = new double[m_D];
            //double[] gradOUT = new double[m_D];
            //for (int d = 0; d < inp.D; d++)
            //{
            //    gradIN[d] = _Grad_uIN[0, d];
            //    gradOUT[d] = _Grad_uOUT[0, d];
            //} 

            //double normIN = Math.Max(gradIN.L2Norm(), m_limit);
            //double normOUT = Math.Max(gradOUT.L2Norm(), m_limit);

            //normOUT = gradIN.InnerProd(gradOUT) < 0.0 ? -normOUT : normOUT;

            //for (int d = 0; d < inp.D; d++)
            //{
            //    //Acc += 0.5 * (_Grad_uIN[0, d] / normIN + _Grad_uOUT[0, d] / normOUT) * (_vIN - _vOUT) * inp.Normal[d];  // consistency term
            //    //Acc += 0.5 * (_Grad_vIN[d]  + _Grad_vOUT[d]) * (_uIN[0] - _uOUT[0]) * inp.Normal[d];  // symmetry term         
            //}

            Acc += Flux(_uIN, _uOUT, _Grad_uIN, _Grad_uOUT, inp.Normal) * (_vIN - _vOUT);
            Acc *= this.m_alpha;

            //Acc -= (_uIN[1] - _uOUT[1]) * (_vIN - _vOUT) * pnlty; // penalty term
            //Acc -= (normIN - normOUT) * (_vIN - _vOUT) * pnlty; // penalty term

            return Acc;
        }

        private enum FluxType {
            Central,

            LaxFriedrich,

            Upwind,

            Godunov
        }

        private double Flux(double[] uIN, double[] uOUT, double[,] grad_uIN, double[,] grad_uOUT, Vector normal) {
            var FType = FluxType.Central;

            double Acc = 0.0;

            double[] MeanGrad = new double[m_D];
            double norm;

            switch(FType) {

                case FluxType.Central:

                for(int d = 0; d < m_D; d++)
                    MeanGrad[d] = 0.5 * (grad_uIN[0, d] + grad_uOUT[0, d]);

                norm = Math.Max(MeanGrad.L2Norm(), m_limit);

                for(int d = 0; d < m_D; d++)
                    Acc += MeanGrad[d] / norm * normal[d];

                break;
                case FluxType.Upwind:

                // How to calculate upwind direction
                double P = 0.0;

                if(P > 0) {
                    for(int d = 0; d < m_D; d++)
                        MeanGrad[d] = grad_uIN[0, d];
                } else {
                    for(int d = 0; d < m_D; d++)
                        MeanGrad[d] = grad_uOUT[0, d];
                }

                norm = Math.Max(MeanGrad.L2Norm(), m_limit);

                for(int d = 0; d < m_D; d++)
                    Acc += MeanGrad[d] / norm * normal[d];

                break;
                case FluxType.LaxFriedrich:

                double[] GradIN = new double[m_D];
                double[] GradOUT = new double[m_D];
                double gamma = 0.0;

                for(int d = 0; d < m_D; d++) {
                    GradIN[d] = grad_uIN[0, d];
                    GradOUT[d] = grad_uOUT[0, d];
                    MeanGrad[d] = GradIN[d] - GradOUT[d];
                }

                double normIN = Math.Max(GradIN.L2Norm(), m_limit);
                double normOUT = Math.Max(GradOUT.L2Norm(), m_limit);
                gamma = 0.1 / Math.Min(normIN, normOUT);

                for(int d = 0; d < m_D; d++)
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

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);

            double[] grad = new double[m_D];
            for(int d = 0; d < inp.D; d++) {
                grad[d] = _Grad_uA[0, d];
            }
            double norm = Math.Max(grad.L2Norm(), m_limit);

            bool IsDirichlet = false;

            if(IsDirichlet) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = 0.0; //this.g_Diri(ref inp);

                for(int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
                    Acc += (_Grad_uA[0, d] / norm) * (_vA) * nd;        // consistency
                    Acc += (_Grad_vA[d] / norm) * (_uA[0] - g_D) * nd;  // symmetry
                }
                Acc *= this.m_alpha;

                Acc -= (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty

            } else {

                double g_N = grad.InnerProd(inp.Normal) / norm;
                //Acc -= (_uA[1] - 0.0) * (_vA - 0) * pnlty; // penalty
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
        protected virtual double GetPenalty(int jCellIn, int jCellOut) {
            double cj_in = InverseLengthScales[jCellIn];
            double mu = m_penalty * cj_in;
            if(jCellOut >= 0) {
                double cj_out = InverseLengthScales[jCellOut];
                mu = Math.Max(mu, m_penalty * cj_out);
            }

            return mu;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var EdgeDiff = new EdgeFormDifferentiator(this, m_D);
            var VolDiff = new VolumeFormDifferentiator(this, m_D);
            return new IEquationComponent[] { EdgeDiff, VolDiff };
        }
    }

}
