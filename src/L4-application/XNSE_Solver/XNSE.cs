using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Multiphase-XDG-solver, with features:
    /// - incompressible two-phase flows.
    /// - solid immersed boundaries (planned).
    /// - three phase contact lines at the domain boundary
    /// - three phase contact lines at the intersection of the immersed solid boundary 
    /// </summary>
    /// <remarks>
    /// Development history:
    /// - Current (jan2021) Maintainers: Beck, Rieckmann, Kummer
    /// - successor of the old XNSE solver <see cref="XNSE_SolverMain"/>, which was mainly used for SFB 1194 and PhD thesis of M. Smuda.
    /// - Quadrature order: saye algorithm can be regarded as a nonlinear transformation to the [-1,1] reference Element. 
    ///   We transform \int f dx to the reference Element, \int f dx = \int f(T) |det D(T)| d\hat{x}
    ///   Suppose f has degree n and suppose the transformation T has degree p, then the integrand in reference space
    ///   has approximately degree <= n * p + (p - 1)
    ///   This is problematic, because we need to find sqrt(n * p + (p - 1)) roots of the level set function, if we want to integrate f exactly.
    ///   This goes unnoticed when verifying the quadrature method via volume/surface integrals with constant f = 1.
    ///   When evaluating a constant function, n = 0, the degree of the integrand immensely simplifies to (p - 1).
    /// </remarks>
    public class XNSE : SolverWithLevelSetUpdater<XNSE_Control> {

        //===========
        // Main file
        //===========
        static void Main(string[] args) {

            //InitMPI();
            //DeleteOldPlotFiles();
            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.ScalingStaticDropletTest_p3_Standard_OneStepGaussAndStokes();
            //BoSSS.Application.XNSE_Solver.Tests.LevelSetUnitTest.LevelSetAdvectiontTest(2, 2, LevelSetEvolution.FastMarching, LevelSetHandling.LieSplitting);
            //BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.MovingDropletTest_rel_p2_Saye_Standard(0.01d, true, SurfaceStressTensor_IsotropicMode.Curvature_Projected, 0.69711d, true, false);
            //BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.BasicThreePhaseTest();
            //Tests.ASUnitTest.HeatDecayTest(r: 0.8598,
            //                                q: -50,
            //                                deg: 3,
            //                                AgglomerationTreshold: 0,
            //                                SolverMode_performsolve: true,
            //                                CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes,
            //                                stm: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);
            //throw new Exception("Remove me");

            void KatastrophenPlot(DGField[] dGFields) {
                Tecplot.PlotFields(dGFields, "AgglomerationKatastrophe", 0.0, 3);
            }

            MultiphaseCellAgglomerator.Katastrophenplot = KatastrophenPlot;
            _Main(args, false, delegate () {
                var p = new XNSE();
                return p;
            });
        }

        /// <summary>
        /// - 3x the velocity degree if convection is included (quadratic term in convection times test function yields tripple order)
        /// - 2x the velocity degree in the Stokes case
        /// </summary>
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
        IncompressibleMultiphaseBoundaryCondMap boundaryMap;


        

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            if(boundaryMap == null)
                boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" });
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


        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            int pVel = VelocityDegree();
            int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
            int D = this.GridData.SpatialDimension;

            if (this.Control.UseSchurBlockPrec) {
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
                // configurations for velocity
                for (int d = 0; d < D; d++) {
                    var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pVel },
                        //DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                    };
                    configsLevel.Add(configVel_d);
                }
                // configuration for pressure
                var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pPrs },
                    //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
                };
                configsLevel.Add(configPres);
            }
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {
            OperatorFactory opFactory = new OperatorFactory();
            
            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            XOP.FreeMeanValue[VariableNames.Pressure] = !GetBcMap().DirichletPressureBoundary;
            XOP.LinearizationHint = LinearizationHint.AdHoc;
            XOP.IsLinear = !(this.Control.PhysicalParameters.IncludeConvection);
            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;
            XOP.Commit();

            return XOP;
        }

        /// <summary>
        /// Setup of the incompressible two-phase Navier-Stokes equation
        /// </summary>
        protected virtual void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            int quadOrder = QuadOrder();
            GetBcMap();

            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            for (int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NavierStokes("A", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A));
                opFactory.AddEquation(new NavierStokes("B", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, LsTrk, config));
            }
            opFactory.AddCoefficient(new SlipLengths(config, VelocityDegree()));
            Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            if(config.physParams.IncludeConvection && config.isTransport) {
                opFactory.AddParameter(new Velocity0(D));
                opFactory.AddParameter(v0Mean);
            }

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);

            if (config.isContinuity) {
                opFactory.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
                opFactory.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
                opFactory.AddEquation(new InterfaceContinuity(config, D, LsTrk, config.isMatInt));
            }

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
                lsUpdater.AddEvolver(VariableNames.LevelSetCG, fourier);
                opFactory.AddParameter(fourier);
                break;

                default:
                throw new NotImplementedException($"option {Control.AdvancedDiscretizationOptions.SST_isotropicMode} is not handled.");
                
            }

            if(Control.UseImmersedBoundary)
                DefineSystemImmersedBoundary(D, opFactory, lsUpdater);
        }

        /// <summary>
        /// Definition of the boundary condition on the immersed boundary, <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// Override to customize.
        /// </summary>
        protected virtual void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            for (int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NSEimmersedBoundary("A", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                opFactory.AddEquation(new NSEimmersedBoundary("B", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
            }

            opFactory.AddEquation(new ImersedBoundaryContinuity("A", "C", 1, config, D, LsTrk));
            opFactory.AddEquation(new ImersedBoundaryContinuity("B", "C", 1, config, D, LsTrk));

            //throw new NotImplementedException("todo");
            opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {

            DGField[] plotFields = this.m_RegisteredFields.ToArray();
            void AddPltField(DGField f) {
                bool add = true;
                foreach(var ff in plotFields) {
                    if(object.ReferenceEquals(f, ff) || (f.Identification == ff.Identification)) {
                        add = false;
                        break;
                    }
                }
                if(add) {
                    f.AddToArray(ref plotFields);
                }
            }
            void AddPltFields(IEnumerable<DGField> fs) {
                foreach(var f in fs)
                    AddPltField(f);
            }

            if (Timestepping?.Parameters != null) {
                AddPltFields(Timestepping.Parameters);
            }
            if (LsUpdater?.Parameters != null) {
                AddPltFields(LsUpdater.Parameters.Values);
                AddPltField(LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet);
                AddPltFields(LsUpdater.InternalFields.Values);
            }

            Tecplot.PlotFields(plotFields, "XNSE_Solver-" + timestepNo, physTime, superSampling);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetFixedTimestep();
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt}");
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"done with time step {TimestepNo}");
            return dt;
        }


    }
}