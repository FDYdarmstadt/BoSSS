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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
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
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Statistic.QuadRules;


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
        public SinglePhaseField c;

        /// <summary>
        /// potential
        /// </summary>
        [InstantiateFromControlFile("mu", "c", IOListOption.Always)]
        public SinglePhaseField mu;

        /// <summary>
        /// residual of 'c'-equation
        /// </summary>
        [InstantiateFromControlFile("c_Resi", "c", IOListOption.Always)]
        protected SinglePhaseField c_Resi;

        /// <summary>
        /// residual of 'mu'-equation
        /// </summary>
        [InstantiateFromControlFile("mu_Resi", "c", IOListOption.Always)]
        protected SinglePhaseField mu_Resi;

        ///// <summary>
        ///// residual of 'curvature'-equation
        ///// </summary>
        //[InstantiateFromControlFile("curvature_Resi", VariableNames.Curvature, IOListOption.Always)]
        //protected SinglePhaseField curvature_Resi;

        /// <summary>
        /// Transport velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            null,
            true, true,
            IOListOption.Always)]
        protected VectorField<SinglePhaseField> Velocity;

        ///// <summary>
        ///// Transport velocity gradients
        ///// </summary>
        //[InstantiateFromControlFile(new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY },
        //    null,
        //    true, true,
        //    IOListOption.Always)]
        //protected VectorField<SinglePhaseField> VelocityGradX;

        ///// <summary>
        ///// Transport velocity gradients
        ///// </summary>
        //[InstantiateFromControlFile(new string[] { VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY },
        //    null,
        //    true, true,
        //    IOListOption.Always)]
        //protected VectorField<SinglePhaseField> VelocityGradY;

        /// <summary>
        /// exact solution, to determine L2-Error, see also <see cref="CahnHilliardControl.ExactSolution_provided"/>.
        /// </summary>
        [InstantiateFromControlFile("cex", "cex", IOListOption.Always)]
        protected SinglePhaseField cex;

        [InstantiateFromControlFile("cDist", "c", IOListOption.Always)]
        SinglePhaseField cDist;
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



        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            DGField[] SolutionFields = new DGField[] { c };

            switch(this.Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelB:
                SolutionFields = SolutionFields.Cat(mu);
                break;
                case CahnHilliardControl.ModelType.modelA:
                case CahnHilliardControl.ModelType.modelC:
                default:
                break;
            }

            //if(this.Control.CurvatureCorrection) {
            //    SolutionFields.Cat(Curvature);
            //}

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
                ResidualFields = ResidualFields.Cat(mu_Resi);
                break;
                case CahnHilliardControl.ModelType.modelA:
                case CahnHilliardControl.ModelType.modelC:
                default:
                break;
            }

            //if(this.Control.CurvatureCorrection) {
            //    ResidualFields.Cat(curvature_Resi);
            //}

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
            base.LsTrk = RealTracker;
        }

        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            //InitMPI(args);
            //DeleteOldPlotFiles();

            // var ctrl = Examples.EllipticDroplet(xRes: res, yRes: res, pDG: pDeg);
            //BoSSS.Application.CahnHilliard.Tests.TestProgram.TestEllipticDroplet();
            //Assert.True(false);
            //using (var solver = new CahnHilliardMain())
            //{
            //    // var C = Examples.EllipticDroplet(20,20,2);
            //    var C = Examples.EllipticDropletPseudo3D(20,20,2);
            //    C.ImmediatePlotPeriod = 1;
            //    C.SuperSampling = 2;
            //    C.NoOfTimesteps = 10;
            //    solver.Init(C);
            //    solver.RunSolverMode();
            //}

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
            string[] paramVar = VariableNames.VelocityVector(D).Cat("c0");//.Cat(VariableNames.LevelSetGradient(D));
            string[] domainVar = new string[] { "c" };
            string[] codomainVar = new string[] { "Res_c" };

            switch(Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelA:
                break;
                case CahnHilliardControl.ModelType.modelB:
                domainVar = domainVar.Cat("mu");
                codomainVar = codomainVar.Cat("Res_mu");
                break;
                case CahnHilliardControl.ModelType.modelC:
                default:
                break;
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
            if (this.Control.includeConvection == true) {
                CHOp.EquationComponents["Res_c"].Add(
                    new __c_Flux(D, () => this.Velocity.ToArray(), m_bcMap)
                );
            }

            switch(Control.ModTyp) {
                case CahnHilliardControl.ModelType.modelA:

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                        new __c_Source(Control.diff)
                    );
                }

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                        new mu_Diffusion(D, penalty_factor, Control.cahn * Control.diff.Sqrt(), m_bcMap, "c")
                        );
                }

                //if(Control.CurvatureCorrection == true) {
                //    CHOp.EquationComponents["Res_c"].Add(
                //        new mu_CurvatureCorrection(D, Control.cahn * Control.diff.Sqrt(), this.Control.UseDirectCurvature)
                //        );

                //    if(!this.Control.UseDirectCurvature) {
                //        CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                //            new curvature_Source(D)
                //            );

                //        CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                //            new curvature_Divergence(D, penalty_factor, 0.001 / this.Control.cahn, LengthScales)
                //            );
                //    } else {
                //        CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                //            new curvature_Direct(D)
                //            );
                //    }
                //}

                break;
                case CahnHilliardControl.ModelType.modelB:

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_c"].Add(
                    new __c_Diffusion(D, penalty_factor, Control.diff, Control.lambda, m_bcMap) // new (L3)
                    //new c_Diffusion(D, penalty_factor, Control.diff, Control.lambda, m_bcMap) // old (L4)
                    );
                }

                if(this.Control.includeDiffusion == true) {
                    CHOp.EquationComponents["Res_mu"].Add(
                         new mu_Diffusion(D, penalty_factor, Control.cahn, m_bcMap, "c") // new (L3)
                         //new phi_Diffusion(D, penalty_factor, Control.cahn, m_bcMap) // old (L4)
                        );
                }

                CHOp.EquationComponents["Res_mu"].Add(
                    new mu_Source("c") // new (L3)
                    //new phi_Source(this.Control.includeDiffusion, Control.cahn) // old (L4)
                    );

                //if(Control.CurvatureCorrection == true) {
                //    CHOp.EquationComponents["Res_mu"].Add(
                //        new mu_CurvatureCorrection(D, Control.cahn)
                //        );

                //        if (!this.Control.UseDirectCurvature) {
                //            CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                //                new curvature_Source(D)
                //                );

                //            CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                //                new curvature_Divergence(D, penalty_factor, 0.001 / this.Control.cahn, LengthScales)
                //                );
                //        } else {
                //            CHOp.EquationComponents["Res_" + VariableNames.Curvature].Add(
                //                new curvature_Direct(D)
                //                );
                //        }
                //}

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

        int reinit = 0;



        private void ComputeDistanceField() {
            GridData GridDat = (GridData)(c.GridDat);

            // compute and project 
            // step one calculate distance field muDist = 0.5 * log(Max(1+mu, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn_old
            // step two project the new phasefield muNew = tanh(muDist/(sqrt(2) * Cahn_new))
            // here done in one step, with default quadscheme
            // ===================
            cDist.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Mu
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

        void CompleteParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            UpdateParameter();
        }

        void UpdateParameter() {
            int D = this.GridData.SpatialDimension;

            ComputeDistanceField();

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
                        m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.mu));
                        break;
                        case CahnHilliardControl.ModelType.modelC:
                        default:
                        break;
                    }

                    //if(this.Control.CurvatureCorrection) {
                    //    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.Curvature));
                    //}
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
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.mu_Resi));
                        break;
                        case CahnHilliardControl.ModelType.modelC:
                        default:
                        break;
                    }

                    //if(this.Control.CurvatureCorrection) {
                    //    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.curvature_Resi));
                    //}
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
                    //break;
                    default:
                    throw new ArgumentOutOfRangeException();
                    //break;
                }

                //if(this.Control.CurvatureCorrection) {
                //    scale = new double[scale.Length + 1];
                //}

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
        private void MassCorrection(BenchmarkQuantities Qnts_old) {
            var Qnts = ComputeBenchmarkQuantities();
            double mass = Qnts.area;
            double surface = Qnts.surface;
            double shape = Qnts.circ;
            double massDiff = Qnts_old.area - mass;

            // we assume the current phasefield is close to the equilibrium tangenshyperbolicus form
            SinglePhaseField cNew = new SinglePhaseField(mu.Basis);
            GridData GridDat = (GridData)(mu.GridDat);
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

                        // evaluate Mu
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
                mass = Qnts.area;

                correction = -(massDiff) / ((Qnts_old.area - mass - massDiff) / (correction));

                // compute and project 
                // step one calculate distance field muDist = 0.5 * log(Max(1+c, eps)/Max(1-c, eps)) * sqrt(2) * Cahn
                // step two project the new phasefield muNew = tanh((cDist + correction)/(sqrt(2) * Cahn))
                // ===================
                cNew.ProjectField(
                    (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                        Debug.Assert(result.Dimension == 2);
                        Debug.Assert(Len == result.GetLength(0));
                        int K = result.GetLength(1); // number of nodes

                        // evaluate Mu
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
                mass = Qnts.area;
                surface = Qnts.surface;
                shape = Qnts.circ;

                massDiff = Qnts_old.area - mass;
                i++;
            }

            Qnts = ComputeBenchmarkQuantities();

            Console.WriteLine($"Performed Mass Correction in {i} iteratins: \n" +
                $"\told mass:           {Qnts_old.area:N4}\n" +
                $"\tuncorrected mass:   {mass_uc:N4}\n" +
                $"\tcorrected mass:     {Qnts.area:N4}");
        }

        private void ConservativityCorrection(BenchmarkQuantities Qnts_old) {
            throw new NotImplementedException();
        }

        /*
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
                        //break;
                        default:
                        throw new ArgumentOutOfRangeException();
                        //break;
                    }

                    //if(this.Control.CurvatureCorrection) {
                    //    config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[m + 1];
                    //}

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
                        //break;
                        default:
                        throw new ArgumentOutOfRangeException();
                        //break;
                    }

                    //if(this.Control.CurvatureCorrection) {
                    //    config[iLevel][m] = new MultigridOperator.ChangeOfBasisConfig() {
                    //        VarIndex = new int[] { m },
                    //        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                    //        DegreeS = new int[] { Math.Max(0, this.Curvature.Basis.Degree - iLevel) }
                    //    };
                    //}
                }

                return config;
            }

        }
        */

        public int m_HMForder {
            get {
                return Math.Max(this.RealLevSet.Basis.Degree, this.c.Basis.Degree)*2;
            }
        }

        public struct BenchmarkQuantities {
            // public BenchmarkQuantities() {
            //     area = 0;
            //     surface = 0;
            //     circ = 0;
            //     concentration = 0;
            //     energy = 0;
            // }
            public BenchmarkQuantities(double __area=0, double __surface=0, double __circ=0, double __concentration=0, double __energy=0)  {
                area = __area;
                surface = __surface;
                circ = __circ;
                concentration = __concentration;
                energy = __energy;
            }


            public double area;
            public double surface;
            public double circ;
            public double concentration;
            public double energy;
        }

        public BenchmarkQuantities ComputeBenchmarkQuantities() {

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

            MultidimensionalArray Mu = new MultidimensionalArray(2);
            MultidimensionalArray GradMu = new MultidimensionalArray(3);
            MultidimensionalArray NormGrad = new MultidimensionalArray(2);

            CellQuadrature.GetQuadrature(new int[] { 1 }, c.GridDat,
                eqs.Compile(c.GridDat, c.Basis.Degree * 2 + 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    int K = EvalResult.GetLength(1);
                    // alloc buffers
                    // -------------

                    if(Mu.GetLength(0) != Length || Mu.GetLength(1) != K) {
                        Mu.Allocate(Length, K);
                        GradMu.Allocate(Length, K, D);
                        NormGrad.Allocate(Length, K);
                    } else {
                        Mu.Clear();
                        GradMu.Clear();
                        NormGrad.Clear();
                    }

                    // chemical potential
                    c.Evaluate(i0, Length, QR.Nodes, Mu.ExtractSubArrayShallow(-1, -1));
                    Mu.ApplyAll(x => 0.25 / (this.Control.cahn.Pow2()) * (x.Pow2() - 1.0).Pow2());

                    for(int d = 0; d < D; d++) {
                        cGrad[d].Evaluate(i0, Length, QR.Nodes, GradMu.ExtractSubArrayShallow(-1, -1, d));
                    }

                    // free surface energy
                    for(int d = 0; d < D; d++) {
                        var GradMu_d = GradMu.ExtractSubArrayShallow(-1, -1, d);
                        NormGrad.Multiply(1.0, GradMu_d, GradMu_d, 1.0, "ik", "ik", "ik");
                    }
                    NormGrad.ApplyAll(x => 0.5 * x);

                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, Mu);
                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, NormGrad);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        energy += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            energy = energy.MPISum();

            return new BenchmarkQuantities( area, surface, circ, concentration, energy );
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
            if (Control.StoreBenchmarkQinFile) {
                if (this.Control.savetodb) {
                    Log = base.DatabaseDriver.FsDriver.GetNewLog("Phasefield_Quantities.txt", sessionID);
                } else {
                    Log = new StreamWriter("Phasefield_Quantities_MPI" + this.MPIRank + ".txt");
                }
            }
        }

        internal List<(double time, BenchmarkQuantities bq)> LogHistory = new List<(double time, BenchmarkQuantities bq)>();


        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {
            var BQnts = ComputeBenchmarkQuantities();

            LogHistory.Add((phystime, BQnts));

            if (Log != null) {
                string line = String.Format($"{TimestepNo}\t{phystime}\t{BQnts.area}\t{BQnts.surface}\t{BQnts.circ}\t{BQnts.concentration}\t{BQnts.energy}");
                Log.WriteLine(line);
                Log.Flush();
            }

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

            if (Log != null) {
                Log.Flush();
                Log.Close();
                Log.Dispose();
                Log = null;
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
            Fields = Fields.Cat(this.cex, this.c, this.mu, this.Velocity,
                //this.gradc0, this.Curvature, this.DCurvature, this.curvature_Resi, this.Correction,
                this.c_Resi, this.mu_Resi, this.cDist);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "CahnHilliard-" + timestepNo + caseStr, phystime, superSampling);
        }

        // the purpose of these fluxes is only consistency in the naming of variables: here, we prefer c, but in the Level-Set context of L3, phi is more suitable
        public class c_Flux : phi_Flux {
            public c_Flux(int D, Func<DGField[]> VelocityGetter, BoundaryCondMap<BoundaryType> __boundaryCondMap) : base(D, VelocityGetter, __boundaryCondMap, "c") {}
        }

        // the purpose of these fluxes is only consistency in the naming of variables: here, we prefer c, but in the Level-Set context of L3, phi is more suitable
        public class c_Source : phi_Source {

            public c_Source(double _diff = 0.0) : base(_diff, "c") {}
        }

        // the purpose of these fluxes is only consistency in the naming of variables: here, we prefer c, but in the Level-Set context of L3, phi is more suitable
        public class c_Diffusion : phi_Diffusion {

            public c_Diffusion(int D, double penalty_const, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
                : base(D, penalty_const, __diff, __lambda, __boundaryCondMap) {}
        }
    }
}
