using BoSSS.Application.XNSE_Solver.LoadBalancing;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using CommandLine;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Reflection;
using ilPSP.LinSolvers;
using BoSSS.Solution.Gnuplot;
using static System.Reflection.Metadata.BlobBuilder;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// Multiphase-XDG-solver, with features:
    /// - incompressible two-phase flows.
    /// - solid immersed boundaries (planned).
    /// - three phase contact lines at the domain boundary
    /// - three phase contact lines at the intersection of the immersed solid boundary 
    /// - the generic control parameter <typeparamref name="T"/> allows derivations of this solver
    /// </summary>
    /// <remarks>
    /// Development history:
    /// - Current (jan2021) Maintainers: Beck, Rieckmann, Kummer
    /// - successor of the old XNSE solver <see cref="XNSE_SolverMain"/>, which was mainly used for SFB 1194 and PhD thesis of M. Smuda.
    /// - Quadrature order: Saye algorithm can be regarded as a nonlinear transformation to the [-1,1] reference Element. 
    ///   We transform $` \int f dx $` to the reference Element, $` \int f dx = \int f(T) |det D(T)| d\hat{x} $`
    ///   Suppose f has degree n and suppose the transformation T has degree p, then the integrand in reference space
    ///   has approximately degree <= n * p + (p - 1)
    ///   This is problematic, because we need to find sqrt(n * p + (p - 1)) roots of the level set function, if we want to integrate f exactly.
    ///   This goes unnoticed when verifying the quadrature method via volume/surface integrals with constant f = 1.
    ///   When evaluating a constant function, n = 0, the degree of the integrand immensely simplifies to (p - 1).
    /// - see also: Extended discontinuous Galerkin methods for two-phase flows: the spatial discretization, F. Kummer, IJNME 109 (2), 2017. 
    /// 
    /// Quick start instructions:
    /// - The Fluid-Fluid-Interphase is level-set No. 0 and is named `Phi` 
    ///   in the initial values <see cref="AppControl.InitialValues"/>, see <see cref="DGField.Identification"/>)
    /// - The immersed boundary interface has the name `Phi2`
    ///   `Phi` is always activated and cannot be turned of; `Phi2` can be activated/deactivated via <see cref="XNSE_Control.UseImmersedBoundary"/>.
    /// - One should not need to set an initial value for `Phi` for simulations where no fluid-fluid interface is required:
    ///   If `Phi` it is identical to 0.0 in the entire domain, it will be set to -1, i.e. the entire domain will be assigned to species "A".
    /// - The vectorial velocity of `Phi2` is set via the fields `VelocityX@Phi2`, `VelocityY@Phi2` and `VelocityZ@Phi2` (in 3D). 
    ///   By specifying a vectorial velocity, one can also specify tangential velocities of the surface. 
    ///   The normal component should match the normal velocity which can be computed from `Phi2`,
    ///   yielding following compatibility condition: 
    /// ```math
    ///     \vec{v} \cdot \frac{\nabla \varphi_2}{| \nabla \varphi_2 |} =  \frac{- \partial_t \varphi_2}{| \nabla \varphi_2 |} 
    /// ```
    /// </remarks>
    public class XNSE : XNSE<XNSE_Control> {

        // ===========
        //  Main file
        // ===========
        static void Main(string[] args) {

            //ilPSP.Environment.NumThreads = 8;
            //InitMPI();
            //BoSSS.Application.XNSE_Solver.Tests.RestartTest.Run_RestartTests(false, LevelSetHandling.Coupled_Once, TimeSteppingScheme.BDF2, false, 3);
            ////BoSSS.Application.XNSE_Solver.Tests.LevelSetUnitTests.LevelSetAdvectionTest2D_reverse(2, 0, LevelSetEvolution.FastMarching, LevelSetHandling.LieSplitting);
            ////DeleteOldPlotFiles();
            ////BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ChannelTest(1, 0.0d, ViscosityMode.FullySymmetric, 0.0d, true, XQuadFactoryHelper.MomentFittingVariants.Saye, NonLinearSolverCode.Picard);
            ////BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ViscosityJumpTest(2, 1, 0.1, ViscosityMode.Standard, XQuadFactoryHelper.MomentFittingVariants.Saye, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);
            //NUnit.Framework.Assert.IsTrue(false, "remove me");

            /*
            var plots = new List<Plot2Ddata>();
            for(int i = 0; i < 4; i++) {
                var p = new Plot2Ddata();
                var x = GenericBlas.Linspace(1, 20, 100);
                var y = x.Select(x_i => Math.Sin(x_i + i)).ToArray();
                p.AddDataGroup("pli" + i, x, y);

                plots.Add(p);
            }

            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            if (rank == 0) {
                Plot2Ddata[,] multi = new Plot2Ddata[1, plots.Count];
                int iCol = 0;
                foreach (var kv in plots) {
                    multi[0, iCol] = kv;
                    iCol++;
                    //var CL = kv.Value.ToGnuplot().PlotCairolatex(xSize: 14, ySize: 12);
                    //CL.WriteMinimalCompileableExample(Path.Combine(OutputDir, "plot_" + kv.Key + ".tex"), kv.Key + ".tex");
                    //kv.Value.SavePgfplotsFile_WA(Path.Combine(OutputDir, kv.Key + ".tex");
                }

                multi.SaveToGIF("waterfall." + DateTime.Now.ToString("yyyyMMMdd_HHmmss") + ".png", 600*4, 600);

            }

            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            throw new Exception("exit2");

            //*/

            /*
            var mtxa = IMatrixExtensions.LoadFromTextFile(@"..\..\..\bin\release\net5.0\weirdo\indef.txt");

            for (int NN = 10; NN <= mtxa.NoOfRows; NN++) {
                Console.WriteLine("NN = " + NN);
                var mtx = mtxa.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NN - 1, NN - 1 }).CloneAs();

                var UpTri = MultidimensionalArray.Create(NN, NN);
                var UpTri2 = MultidimensionalArray.Create(NN, NN);
                mtx.SymmetricLDLInversion(UpTri, null);
                mtx.GramSchmidt(UpTri2, null);


                Console.WriteLine("UpTri vs. LDL " + UpTri2.Storage.L2Dist(UpTri.Storage));

                var LoTri = UpTri.TransposeTo();
                var Eye = LoTri.GEMM(mtx, UpTri);
                var LoTri2 = UpTri2.TransposeTo();
                var Eye2 = LoTri2.GEMM(mtx, UpTri2);
                //L.SaveToTextFile(@"..\..\..\bin\release\net5.0\weirdo\GS.txt");
                //Eye.SaveToTextFile(@"..\..\..\bin\release\net5.0\weirdo\ID.txt");
                Eye.AccEye(-1.0);
                Eye2.AccEye(-1.0);
                Console.WriteLine("Ortho Error: " + Eye.InfNorm() + "    " + Eye2.InfNorm());

            }
            //DeleteOldPlotFiles();
            //BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ChannelTest(3, 0.0d, ViscosityMode.Standard, 1.0471975511965976d, XQuadFactoryHelper.MomentFittingVariants.Saye, NonLinearSolverCode.Newton);
            //BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.TaylorCouetteConvergenceTest_2Phase_Curvature_Proj_Soff_p2(NonLinearSolverCode.Picard);
            NUnit.Framework.Assert.IsTrue(false, "remove me"); 
            */

            
            {
                XNSE._Main(args, false, delegate () {
                    var p = new XNSE();
                    return p;
                });

                //Tmeas.Memtrace.Flush();
                //Tmeas.Memtrace.Close();
            }
            //*/
            ilPSP.LinSolvers.BlockMsrMatrix.PrintPerfStat();

            
        }
    }

    /// <summary>
    /// Generic versions which should be used for derivatives 
    /// </summary>
    public class XNSE<T> : SolverWithLevelSetUpdater<T> where T : XNSE_Control, new() {

        public override void Init(AppControl control) {
                

            base.Init(control);
            var ctrl = (control as XNSE_Control);


            if (ctrl.Rigidbody.IsInitialized())
                ctrl.Rigidbody.ArrangeAll(ctrl);
        }

        /// <summary>
        /// - 3x the velocity degree if convection is included (quadratic term in convection times test function yields triple order)
        /// - 2x the velocity degree in the Stokes case
        /// </summary>
        /// <remarks>
        /// Note: 
        /// Sayes algorithm can be regarded as a nonlinear transformation to the [-1,1] reference Element. 
        /// We transform $`\int f dx $` to the reference Element, $`\int f dx = \int f(T) |det D(T)| d\hat{x} $`
        /// Suppose $`f$` has degree $`n$` and suppose the transformation $`T$` has degree $`p$`, then the integrand in reference space
        /// has approximately degree $`\leq n * p + (p - 1) $`
        /// This is problematic, because we need to find $`\sqrt(n * p + (p - 1))$` roots of the level set function, if we want to integrate $`f$` exactly.
        /// This goes unnoticed when verifying the quadrature method via volume/surface integrals with constant $`f = 1$`.
        /// When evaluating a constant function, $`n = 0$`, the degree of the integrand immensely simplifies to $`(p - 1)$`.        
        /// </remarks>
        override public int QuadOrder() {
            if(Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye
               && Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes) {
                throw new ArgumentException($"The XNSE solver is only verified for cut-cell quadrature rules " +
                    $"{XQuadFactoryHelper.MomentFittingVariants.Saye} and {XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes}; " +
                    $"you have set {Control.CutCellQuadratureType}, so you are notified that you reach into unknown territory; " +
                    $"If you do not know how to remove this exception, you should better return now!");
            }

            //QuadOrder
            int degU = VelocityDegree();
            int quadOrder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if(this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                //See remarks
                quadOrder *= 2;
                quadOrder += 1;
            } 

            return quadOrder;
        }


        /// <summary>
        /// Current velocity
        /// </summary>
        public VectorField<XDGField> Velocity {
            get {
                int D = this.GridData.SpatialDimension;
                return new VectorField<XDGField>(this.CurrentState.Fields.Take(D).Select(f => ((XDGField)f)).ToArray());
            }
        }

        /// <summary>
        /// Current Pressure
        /// </summary>
        public XDGField Pressure {
            get {
                int D = this.GridData.SpatialDimension;
                return (XDGField)this.CurrentState.Fields[D];
            }
        }



        /// <summary>
        /// Usually, the term "DG order of the calculation" means the velocity degree.
        /// </summary>
        protected int VelocityDegree() {
            int pVel;
            if(this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v)) {
                pVel = v.Degree;
            } else if(this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1)) {
                pVel = v1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        /// <summary>
        /// 
        /// </summary>
        protected int PressureDegree() {
            int pPres;
            if(Pressure != null)
                return Pressure.Basis.Degree;
            if(this.Control.FieldOptions.TryGetValue(VariableNames.Pressure, out FieldOpts v)) {
                pPres = v.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Pressure not found");
            }
            return pPres;
        }

        

        private IncompressibleMultiphaseBoundaryCondMap m_boundaryMap;

        /// <summary>
        /// Relation between 
        /// - edge tags (<see cref="Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>, passed to equation components via <see cref="BoSSS.Foundation.CommonParams.EdgeTag"/>)
        /// - boundary conditions specified in the control object (<see cref="AppControl.BoundaryValues"/>)
        /// </summary>
        protected IncompressibleMultiphaseBoundaryCondMap boundaryMap {
            get {
                if(m_boundaryMap == null)
                    m_boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" });
                return m_boundaryMap;
            }
        }

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            return boundaryMap;
        }


        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;

            if(iLevSet == 0) {
                // +++++++++++++++++++
                // the fluid interface 
                // +++++++++++++++++++

                // averaging at interface:
                ILevelSetParameter levelSetVelocity = new LevelSetVelocity(VariableNames.LevelSetCG, D, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                return levelSetVelocity;
            } else if(iLevSet == 1) {
                // +++++++++++++++++++++
                // the immersed boundary
                // +++++++++++++++++++++

                string[] VelocityNames = VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.VelocityVector(D)).ToArray();
                ScalarFunctionTimeDep[] VelFuncs = new ScalarFunctionTimeDep[D];
                for(int d = 0; d < D; d++) {
                    Control.InitialValues_EvaluatorsVec.TryGetValue(VelocityNames[d], out VelFuncs[d]);
                }

                ILevelSetParameter levelSetVelocity = new ExplicitLevelSetVelocity(VariableNames.LevelSetCGidx(1), VelFuncs);
                return levelSetVelocity;
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        protected override int NoOfLevelSets {
            get {
                if(Control.UseImmersedBoundary)
                    return 2;
                else
                    return 1;
            }
        }

        /// <summary>
        /// Either fluids A and B; or A, B and solid C.
        /// </summary>
        protected override Array SpeciesTable {
            get {
                if(Control.UseImmersedBoundary) {
                    var r = new string[2, 2];
                    r[0, 0] = "A";
                    r[0, 1] = "C"; // solid
                    r[1, 0] = "B";
                    r[1, 1] = "C"; // also solid
                    return r;
                } else {
                    return new[] { "A", "B" };
                }
            }
        }

        /*
        /// <summary>
        /// Cell-performance classes:
        /// cell performance class equals number of species present in that cell
        /// </summary>
        protected override void GetCellPerformanceClasses(out int NoOfClasses, out int[] CellPerfomanceClasses, int TimeStepNo, double physTime) {
            (NoOfClasses, CellPerfomanceClasses) = CellClassifier.ClassifyCells(this, this.Control.DynamicLoadbalancing_ClassifierType);
        }
        */

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            using (var tr = new FuncTrace()) {
                int pVel = VelocityDegree();
                int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
                int D = this.GridData.SpatialDimension;

                
                if (this.Control.UseSchurBlockPrec) {
                    tr.Info($"pre-precond, level {iLevel}: using {MultigridOperator.Mode.SchurComplement}");


                    // using a Schur complement for velocity & pressure
                    var confMomConti = new MultigridOperator.ChangeOfBasisConfig();
                    for (int d = 0; d < D; d++) {
                        d.AddToArray(ref confMomConti.VarIndex);
                        //Math.Max(1, pVel - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                        pVel.AddToArray(ref confMomConti.DegreeS);
                    }
                    D.AddToArray(ref confMomConti.VarIndex);
                    //Math.Max(0, pPrs - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                    pPrs.AddToArray(ref confMomConti.DegreeS);

                    confMomConti.mode = MultigridOperator.Mode.SchurComplement;

                    configsLevel.Add(confMomConti);
                } else {
                    tr.Info($"pre-precond, level {iLevel}: using {MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite} and {MultigridOperator.Mode.IdMass_DropIndefinite}");


                    // configurations for velocity
                    for (int d = 0; d < D; d++) {
                        var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { pVel },
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                            VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                        };
                        configsLevel.Add(configVel_d);
                    }
                    // configuration for pressure
                    var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pPrs },
                        //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) }, // p-multigrid reduction
                        mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
                    };
                    configsLevel.Add(configPres);
                }
            }
        }
        /// <summary>
        /// Operator/equation assembly
        /// </summary>
        protected override XDifferentialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {

            OperatorFactory opFactory = new OperatorFactory();

            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XDifferentialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            FinalOperatorSettings(XOP, D);
            XOP.FluxesAreNOTMultithreadSafe = false; // enable multi-threaded evaluation of fluxes
            XOP.Commit();

            return XOP;
        }

        /// <summary>
        /// Misc adjustments to the spatial operator before calling <see cref="IDifferentialOperator.Commit"/>
        /// </summary>
        protected virtual void FinalOperatorSettings(XDifferentialOperatorMk2 XOP, int D) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                XOP.FreeMeanValue[VariableNames.Pressure] = !GetBcMap().DirichletPressureBoundary;
                XOP.IsLinear = !(this.Control.PhysicalParameters.IncludeConvection || Control.NonlinearCouplingSolidFluid);
                XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;
                tr.Info("Going with agglomeration threshold: " + XOP.AgglomerationThreshold);


                if (XOP.IsLinear == true) {
                    XOP.LinearizationHint = LinearizationHint.AdHoc;
                } else {
                    if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                        if (UseAdHocLinearization)
                            XOP.LinearizationHint = LinearizationHint.AdHoc;
                        else
                            XOP.LinearizationHint = LinearizationHint.GetJacobiOperator;
                        //XOP.LinearizationHint = LinearizationHint.FDJacobi;
                    } else {
                        XOP.LinearizationHint = LinearizationHint.AdHoc;
                    }

                    //(UseAdHocLinearization ? LinearizationHint.AdHoc : LinearizationHint.GetJacobiOperator);
                }
                tr.Info("Linearization hint: " + XOP.LinearizationHint);

                // elementary checks on operator
                if (XOP.CodomainVar.IndexOf(EquationNames.ContinuityEquation) != D)
                    throw new ApplicationException("Operator configuration messed up.");
                if (XOP.DomainVar.IndexOf(VariableNames.Pressure) != D)
                    throw new ApplicationException("Operator configuration messed up.");
                for (int d = 0; d < D; d++) {
                    if (XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationComponent(d)) != d)
                        throw new ApplicationException("Operator configuration messed up.");
                    if (XOP.DomainVar.IndexOf(VariableNames.Velocity_d(d)) != d)
                        throw new ApplicationException("Operator configuration messed up.");
                }

                PrintConfiguration();
            }
        }

        bool UseAdHocLinearization {
            get {
                //return true;

                
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard)
                    return true;
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton)
                    return false;
                throw new ArgumentException("unkonwn Nonlinear solver: " + this.Control.NonLinearSolver.SolverCode);
                //*/
            }
        }


        /// <summary>
        /// Setup of the incompressible two-phase Navier-Stokes equation
        /// </summary>
        protected virtual void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            int quadOrder = QuadOrder();
            GetBcMap();

            XNSE_OperatorConfiguration config = new XNSE_OperatorConfiguration(this.Control);

            // === momentum equations === //
            for (int d = 0; d < D; ++d) {
                DefineMomentumEquation(opFactory, config, d, D);

                // Add Gravitation
                if(config.isGravity) {
                    var GravA = Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A, Control.GetGravity("A", d));
                    opFactory.AddParameter(GravA);
                    var GravB = Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B, Control.GetGravity("B", d));
                    opFactory.AddParameter(GravB);
                }

                // Add additional volume forces
                if(config.isVolForce) {
                    var VolForceA = VolumeForce.CreateFrom("A", d, D, Control, Control.GetVolumeForce("A", d));
                    opFactory.AddParameter(VolForceA);
                    var VolForceB = VolumeForce.CreateFrom("B", d, D, Control, Control.GetVolumeForce("B", d));
                    opFactory.AddParameter(VolForceB);
                }
            }

            // === continuity equation === //
            if (config.isContinuity) {
                DefineContinuityEquation(opFactory, config, D);
            }

            // === additional parameters === //
            opFactory.AddCoefficient(new SlipLengths(config, VelocityDegree()));
            Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            
            if (config.physParams.IncludeConvection && config.isTransport) {
                opFactory.AddParameter(v0Mean);
                if(this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard)
                    opFactory.AddParameter(new Velocity0(D));
            }

            // === level set related parameters === //
            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);

            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, v0Mean);
            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, normalsParameter);
            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                    MaxSigma maxSigmaParameter = new MaxSigma(Control.PhysicalParameters, Control.AdvancedDiscretizationOptions, QuadOrder(), Control.dtFixed);
                    opFactory.AddParameter(maxSigmaParameter);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, maxSigmaParameter);
                    BeltramiGradient lsBGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsBGradient);
                    break;

                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                    BeltramiGradient lsGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradient);
                    break;

                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                    BeltramiGradientAndCurvature lsGradientAndCurvature =
                        FromControl.BeltramiGradientAndCurvature(Control, "Phi", quadOrder, D);
                    opFactory.AddParameter(lsGradientAndCurvature);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradientAndCurvature);
                    break;

                case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                    FourierLevelSet ls = (FourierLevelSet)lsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet;
                    var fourier = new FourierEvolver(
                        VariableNames.LevelSetCG,
                        ls,
                        Control.FourierLevSetControl,
                        Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Curvature].Degree);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, fourier);
                    //lsUpdater.AddEvolver(VariableNames.LevelSetCG, fourier);
                    opFactory.AddParameter(fourier);
                    break;

                default:
                    throw new NotImplementedException($"option {Control.AdvancedDiscretizationOptions.SST_isotropicMode} is not handled.");

            }

            if (Control.UseImmersedBoundary)
                DefineSystemImmersedBoundary(D, opFactory, lsUpdater);
        }

        /// <summary>
        /// Override this method to customize the assembly of the momentum equation
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="d">Momentum component index</param>
        /// <param name="D">Spatial dimension (2 or 3)</param>
        virtual protected void DefineMomentumEquation(OperatorFactory opFactory, XNSE_OperatorConfiguration config, int d, int D) {

            // === linearized or parameter free variants, difference only in convective term === //
            if (UseAdHocLinearization) {
                opFactory.AddEquation(new NavierStokes("A", d, D, boundaryMap, config));
                opFactory.AddEquation(new NavierStokes("B", d, D, boundaryMap, config));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, config, config.isMovingMesh));
            } else {
                if (this.Control.NonLinearSolver.SolverCode != NonLinearSolverCode.Newton)
                    throw new ApplicationException("illegal configuration");
                opFactory.AddEquation(new NavierStokes_Newton("A", d, D, boundaryMap, config));
                opFactory.AddEquation(new NavierStokes_Newton("B", d, D, boundaryMap, config));
                opFactory.AddEquation(new NSEInterface_Newton("A", "B", d, D, boundaryMap, config, config.isMovingMesh));
            } 

            opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, config));
        }

        /// <summary>
        /// Override this method to customize the assembly of the continuity equation
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="D">Spatial dimension (2 or 3)</param>
        virtual protected void DefineContinuityEquation(OperatorFactory opFactory, XNSE_OperatorConfiguration config, int D) {
            opFactory.AddEquation(new Continuity("A", config, D, boundaryMap));
            opFactory.AddEquation(new Continuity("B", config, D, boundaryMap));
            opFactory.AddEquation(new InterfaceContinuity("A", "B", config, D, config.isMatInt));
        }

        /// <summary>
        /// Definition of the boundary condition on the immersed boundary (fluid-solid boundary, level-set 1), 
        /// <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// Override to customize.
        /// </summary>
        protected virtual void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            XNSE_OperatorConfiguration config = new XNSE_OperatorConfiguration(this.Control);

           
            if (this.Control.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature) 
                BoSSS.Foundation.XDG.Quadrature.BruteForceSettingsOverride.doubleCutCellOverride = true;

            for (int d = 0; d < D; ++d) {
                // so far only no slip!
                if (UseAdHocLinearization) {
                    opFactory.AddEquation(new NSEimmersedBoundary("A", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                    opFactory.AddEquation(new NSEimmersedBoundary("B", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                } else {
                    opFactory.AddEquation(new NSEimmersedBoundary_Newton("A", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                    opFactory.AddEquation(new NSEimmersedBoundary_Newton("B", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                }

                // surface tension on IBM
                if (config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine && config.physParams.Sigma != 0.0) {
                    opFactory.AddEquation(new NSEimmersedBoundary_SurfaceTension("A", "B", d, D, 1));
                }

                // GNBC
                if (config.dntParams.IBM_BoundaryType != IBM_BoundaryType.NoSlip) {
                    opFactory.AddEquation(new NSEimmersedBoundary_GNBC("A", "B", d, D, config.getPhysParams, 1));
                }
            }

            opFactory.AddEquation(new ImmersedBoundaryContinuity("A", "C", 1, config, D));
            opFactory.AddEquation(new ImmersedBoundaryContinuity("B", "C", 1, config, D));

            //throw new NotImplementedException("todo");
            opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));

            if (config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine || config.dntParams.IBM_BoundaryType != IBM_BoundaryType.NoSlip) {
                var normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[1]).Basis.Degree, VariableNames.LevelSetCGidx(1));
                opFactory.AddParameter(normalsParameter);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCGidx(1), normalsParameter);
            }
        }

        public List<(double, double)> ImbalanceTrack;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (var f = new FuncTrace()) {
                if ((int)this.Control.TimeSteppingScheme >= 100 && this.Control.TimesteppingMode != AppControl._TimesteppingMode.Steady) {
                    
                    // this is a RK scheme, set here the maximum 
                    dt = this.GetTimestep();
                    if (phystime + dt > this.Control.Endtime) {
                        Console.WriteLine("restricting time-step to reach end-time");
                        dt = this.Control.Endtime - phystime;
                    }
                    if (TimestepNo == 1) dt = Math.Min(dt, 1e-3 * this.Control.dtFixed); // small start timestep, to get the levelset rolling
                } else {
                    // this is a BDF or non-adaptive scheme, use the base implementation, i.e. the fixed timestep
                    dt = base.GetTimestep();
                }

                int NoOfCutCells = this.LsTrk.Regions.GetNearFieldMask(1).NoOfItemsLocally;
                int NoOfCells = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int NoOfCutCells_Min = NoOfCutCells.MPIMin();
                int NoOfCutCells_Max = NoOfCutCells.MPIMax();
                int NoOfCells_Min = NoOfCells.MPIMin();
                int NoOfCells_Max = NoOfCells.MPIMax();
                int NoOfCutCellsTot = NoOfCutCells.MPISum();
                long NoOfCellsTot = this.GridData.CellPartitioning.TotalLength;
                double AvgCutCells = NoOfCutCellsTot / (double)this.MPISize;
                double AvgCells = NoOfCellsTot / (double)this.MPISize;
                double RelCellsInbalance = (double)(NoOfCells_Max - NoOfCells_Min) / NoOfCells_Max;
                double RelCutCellsInbalance = NoOfCutCells_Max == 0 ? 0 : (double)(NoOfCutCells_Max - NoOfCutCells_Min) / NoOfCutCells_Max;

                Console.WriteLine($"All Cells: min={NoOfCells_Min} max={NoOfCells_Max} avg={AvgCells:G5} inb={RelCellsInbalance:G4} tot={NoOfCellsTot}");
                Console.WriteLine($"Cut Cells: min={NoOfCutCells_Min} max={NoOfCutCells_Max} avg={AvgCutCells:G5} inb={RelCutCellsInbalance:G4}, tot={NoOfCutCellsTot}");

                if (ImbalanceTrack != null) ImbalanceTrack.Add((RelCellsInbalance, RelCutCellsInbalance));

                Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt} ...");
                bool success = Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

                Console.WriteLine($"Done with time step {TimestepNo}; solver success: {success}");
                GC.Collect();
                if (Control.FailOnSolverFail && !success) {
                    PlotCurrentState(phystime, TimestepNo, this.Control.SuperSampling);
                    SaveToDatabase(TimestepNo, phystime);
                    throw new ArithmeticException("Solver did not converge.");
                }

                return dt;
            }
        }
        protected virtual List<DGField> GetInterfaceVelocity(){
            MPICollectiveWatchDog.Watch();
            int D = this.GridData.SpatialDimension;
            var cm = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
            var PhysParam = this.Control.PhysicalParameters;
            var gdat = this.GridData;
            var Tracker = this.LsTrk;

            VectorField<XDGField> Velocity = new VectorField<XDGField>(D.ForLoop(d => (XDGField)m_RegisteredFields.SingleOrDefault(s => s.Identification == VariableNames.Velocity_d(d))));
            var LevelSet = this.m_RegisteredFields.Where(s => s.Identification == VariableNames.LevelSetCG).SingleOrDefault();
            int p = LevelSet.Basis.Degree;

            Basis b = new Basis(gdat, 2 * p);
            XDGBasis xb = new XDGBasis(Tracker, 2 * p);

            VectorField<DGField> Normal = new VectorField<DGField>(D.ForLoop(d => m_RegisteredFields.SingleOrDefault(s => s.Identification == VariableNames.NormalVectorComponent(d)))).CloneAs();
            // this normal is the direct gradient of phi, we normalize first, but only in cut cells!!
            {
                Normal.Clear();
                Normal.Gradient(1.0, LevelSet);
                SinglePhaseField Normalizer = new SinglePhaseField(b);
                for (int d = 0; d < D; d++) {
                    Normalizer.ProjectPow(1.0, Normal[d], 2.0, cm);
                }
                var NormalizerTemp = Normalizer.CloneAs();
                Normalizer.Clear();
                Normalizer.ProjectPow(1.0, NormalizerTemp, -0.5, cm);

                Normal.ScalePointwise(1.0, Normalizer);
            }

            List<XDGField> xInterfaceVelocity = new List<XDGField>(); // = this.m_RegisteredFields.Where(s => s.Identification.Contains("IntefaceVelocityXDG")).ToList();
            List<DGField> InterfaceVelocity = new List<DGField>(); // = this.m_RegisteredFields.Where(s => s.Identification.Contains("IntefaceVelocityDG")).ToList();
            for (int d = 0; d < D; d++) {
                var xField = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == $"IntefaceVelocityXDG_{d}").SingleOrDefault();
                if (xField == null) {
                    xField = new XDGField(xb, $"IntefaceVelocityXDG_{d}");
                    this.RegisterField(xField);
                }
                xInterfaceVelocity.Add(xField);

                var Field = (DGField)this.m_RegisteredFields.Where(s => s.Identification == $"IntefaceVelocityDG_{d}").SingleOrDefault();
                if (Field == null) {
                    Field = new SinglePhaseField(b, $"IntefaceVelocityDG_{d}");
                    this.RegisterField(Field);
                }
                InterfaceVelocity.Add(Field);
            }

            for (int d = 0; d < D; d++) {
                var xField = xInterfaceVelocity[d];
                xField.Clear();
                xField.AccLaidBack(1.0, Velocity[d], cm);

                var Field = InterfaceVelocity[d];
                Field.Clear();
                Field.AccLaidBack(PhysParam.rho_A, xField.GetSpeciesShadowField("A"), cm);
                Field.AccLaidBack(PhysParam.rho_B, xField.GetSpeciesShadowField("B"), cm);
                Field.Scale(1.0 / (PhysParam.rho_A + PhysParam.rho_B), cm);
            }

            //// Project normal
            //var NormalInterfaceVelocity = InterfaceVelocity.Select(s => s.CloneAs()).ToArray();
            //NormalInterfaceVelocity.ForEach(s => s.Clear());
            //for (int i = 0; i < D; i++) {
            //    var temp = NormalInterfaceVelocity[i].CloneAs();
            //    temp.Clear();
            //    for (int j = 0; j < D; j++) {
            //        temp.ProjectProduct(1.0, Normal[j], InterfaceVelocity[i], cm, true);
            //    }
            //    NormalInterfaceVelocity[i].ProjectProduct(1.0, temp, Normal[i], cm, false);
            //}
            //for (int i = 0; i < D; i++) {
            //    InterfaceVelocity[i].Clear();
            //    InterfaceVelocity[i].Acc(1.0, NormalInterfaceVelocity[i]);
            //}
            return InterfaceVelocity;
        }

        public override double GetTimestep() {
            double dt = this.Control.dtFixed;
            double s = 0.5; // safety factor
            int D = this.GridData.SpatialDimension;

            List<DGField> LevelSetVelocity;
            IList<string> LevelSetVelocityNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(VariableNames.LevelSetCG, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            if (this.RegisteredFields.Any(s => LevelSetVelocityNames.Any(x => x == s.Identification)))
                LevelSetVelocity = GetInterfaceVelocity();
            else
                LevelSetVelocity = null;


            var CC = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
            if (CC.NoOfItemsLocally > 0) {
                // search minimum grid width in cutcells
                var hmin = double.MaxValue;
                foreach (var chunk in CC) {
                    for (int i = chunk.i0; i < chunk.i0 + chunk.Len; i++) {
                        hmin = Math.Min(this.GridData.iGeomCells.h_min[i], hmin);
                    }
                }

                // get level set degree
                int p = 0;
                try {
                    p = ((DGField)LsTrk.LevelSets[0]).Basis.Degree;
                } catch {
                    Console.WriteLine("Cannot determine DG order of level set");
                }

                // capillary timestep
                {
                    var physParams = this.Control.PhysicalParameters;
                    if (physParams.Sigma != 0) {
                        double dt_cap = s * Math.Sqrt((physParams.rho_A + physParams.rho_B) * Math.Pow(hmin * (p + 1), 3.0) / (2 * Math.PI * Math.Abs(physParams.Sigma)));
                        if (dt_cap < dt) {
                            dt = Math.Min(dt_cap, dt);
                            Console.WriteLine("Restricting time step size to: {0}, due to capillary timestep restriction", dt_cap);
                        }
                    }
                }

                // level set cfl
                {   
                    if (LevelSetVelocity != null) {
                        // At this point the level set velocity is not updated to the correct value
                        //VectorField<DGField> LevelSetVelocity = new VectorField<DGField>(D.ForLoop(d => RegisteredFields.SingleOrDefault(s => s.Identification == LevelSetVelocityNames[d])));

                        double dt_cfl = this.GridData.ComputeCFLTime(LevelSetVelocity.ToArray(), 10000, CC);
                        dt_cfl *= s / Math.Pow(p, 2);
                        if (dt_cfl < dt) {
                            dt = Math.Min(dt_cfl, dt);
                            Console.WriteLine("Restricting time step size to: {0}, due to level set cfl", dt_cfl);
                        }
                    }
                }
            }

            // determine Minimum timestep over all processess
            dt = MPIExtensions.MPIMin(dt);

            return dt;
        }

        bool PrintOnlyOnce = true;
        private void PrintConfiguration() {
            if (PrintOnlyOnce) {
                PrintOnlyOnce = false;
                XNSE_OperatorConfiguration config = new XNSE_OperatorConfiguration(this.Control);

                Console.WriteLine("=============== {0} ===============", "Operator Configuration");
                PropertyInfo[] properties = typeof(XNSE_OperatorConfiguration).GetProperties();
                foreach (PropertyInfo property in properties) {
                    if (property.PropertyType == typeof(bool)) {
                        bool s = (bool)property.GetValue(config);
                        Console.WriteLine("     {0,-30}:{1,3}", property.Name, "[" + (s == true ? "x" : " ") + "]");
                    }
                }

                Console.WriteLine("=============== {0} ===============", "Linear Solver Configuration");
                Console.WriteLine("     {0,-30}:{1}", "Solvercode", this.Control.LinearSolver.Name);
                
                Console.WriteLine("=============== {0} ===============", "Nonlinear Solver Configuration");
                Console.WriteLine("     {0,-30}:{1}", "Solvercode", this.Control.NonLinearSolver.SolverCode);
                Console.WriteLine("     {0,-30}:{1}", "Convergence Criterion", this.Control.NonLinearSolver.ConvergenceCriterion);
                Console.WriteLine("     {0,-30}:{1}", "Globalization", this.Control.NonLinearSolver.Globalization);
                Console.WriteLine("     {0,-30}:{1}", "Minsolver Iterations", this.Control.NonLinearSolver.MinSolverIterations);
                Console.WriteLine("     {0,-30}:{1}", "Maxsolver Iterations", this.Control.NonLinearSolver.MaxSolverIterations);
                Console.WriteLine("==============={0}===============", "========================");


            }
        }


        /// <summary>
        /// delegate for the initialization of previous timesteps from an analytic solution
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="Time"></param>
        /// <param name="St"></param>
        protected override void BDFDelayedInitSetIntial(int TimestepIndex, double Time, DGField[] St) {
            using (new FuncTrace()) {
                Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, Time);

                // level-set
                // ---------
                this.LsUpdater.LevelSets["Phi"].DGLevelSet.ProjectField(X => this.Control.Phi(X, Time));
                //this.LsUpdater.LevelSets[0].CGLevelSet..ProjectField(X => this.Control.Phi(X, Time));

                //this.LsTrk.UpdateTracker(Time, incremental: true);

                // solution
                // --------
                int D = this.LsTrk.GridDat.SpatialDimension;

                for (int d = 0; d < D; d++) {
                    XDGField _u = (XDGField)St[d];
                    _u.Clear();
                    _u.GetSpeciesShadowField("A").ProjectField(X => this.Control.ExactSolutionVelocity["A"][d](X, Time));
                    _u.GetSpeciesShadowField("B").ProjectField((X => this.Control.ExactSolutionVelocity["B"][d](X, Time)));
                }
                XDGField _p = (XDGField)St[D];
                _p.Clear();
                _p.GetSpeciesShadowField("A").ProjectField(X => this.Control.ExactSolutionPressure["A"](X, Time));
                _p.GetSpeciesShadowField("B").ProjectField((X => this.Control.ExactSolutionPressure["B"](X, Time)));
            }
        }



    }
}