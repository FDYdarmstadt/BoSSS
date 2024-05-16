using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using ZwoLevelSetSolver.Boundary;
using ZwoLevelSetSolver.SolidPhase;
using ZwoLevelSetSolver.ContactLine;
using NSEVariableNames = BoSSS.Solution.NSECommon.VariableNames;
using NSEEquationNames = BoSSS.Solution.NSECommon.EquationNames;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using BoSSS.Solution.Utils;
using NUnit.Framework;
using BoSSS.Solution;
using System.Diagnostics;
using System.ComponentModel;
using BoSSS.Solution.Tecplot;

namespace ZwoLevelSetSolver {

    /// <summary>
    /// ZLS: Two-Level-Set-Solver;
    /// - Main purpose seems to be Euler/Euler fluid-structure interaction
    /// - Original author: Lauritz Beck
    /// </summary>
    public class ZLS : XNSE<ZLS_Control> {

        /// <summary>
        /// DG polynomial degree of the displacement variable.
        /// 
        /// Note: Usually, the term "DG order of the calculation" means the velocity degree.
        /// despite dis, displacement fields or the pressure might have other DG polynomial degrees.
        /// </summary>
        protected int DisplacementDegree() {
            int pDspl;
            if(this.Control.FieldOptions.TryGetValue("Displacement*", out FieldOpts v)) {
                pDspl = v.Degree;
            } else if(this.Control.FieldOptions.TryGetValue(VariableNames.DisplacementX, out FieldOpts v1)) {
                pDspl = v1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pDspl;
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            int D = this.GridData.SpatialDimension;
            int pVel = VelocityDegree();
            int pPrs = PressureDegree();
            int pDispl = DisplacementDegree();

            // configurations for velocity
            for (int d = 0; d < D; d++) {
                var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pVel },
                    mode = MultigridOperator.Mode.IdMass,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(NSEVariableNames.VelocityVector(D)[d]) }
                };
                configsLevel.Add(configVel_d);
            }
            // configuration for pressure
            var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pPrs },
                mode = MultigridOperator.Mode.IdMass,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(NSEVariableNames.Pressure) }
            };
            configsLevel.Add(configPres);
            // configuration for displacements
            for (int d = 0; d < D; d++) {
                var configDisplacement = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pDispl},
                    mode = MultigridOperator.Mode.IdMass,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.DisplacementVector(D)[d]) }
                };
                configsLevel.Add(configDisplacement);
            }
        }

        protected override void SetInitial(double t) {
            base.SetInitial(t);
            //
        }

        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystem(D, opFactory, lsUpdater);
            if(this.Control.NonLinearSolver.SolverCode != NonLinearSolverCode.Newton) {
                throw new NotSupportedException();
            }
            DefineSolidPhase(D, opFactory, lsUpdater);
        }

        void DefineSolidPhase(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            IncompressibleMultiphaseBoundaryCondMap boundaryMap = this.boundaryMap;

            for(int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NavierCauchy("C", Control.Material, d, D, boundaryMap));
                opFactory.AddEquation(new DisplacementEvolution("C", Control.Material, d, D, boundaryMap));
                opFactory.AddEquation(new Dummy("A", VariableNames.DisplacementVector(D)[d], EquationNames.DisplacementEvolutionComponent(d)));
                opFactory.AddEquation(new Dummy("B", VariableNames.DisplacementVector(D)[d], EquationNames.DisplacementEvolutionComponent(d)));
                opFactory.AddParameter(Gravity.CreateFrom("C", d, D, Control, Control.Material.Density, Control.GetGravity("C", d)));
            }
            var continuityEquation = new SolidPhase.Continuity("C", D, Control.Material);
            
            opFactory.AddEquation( continuityEquation);

        }

        protected override void FinalOperatorSettings(XDifferentialOperatorMk2 XOP, int D) {
            base.FinalOperatorSettings(XOP, D);
            //XOP.FreeMeanValue[NSEVariableNames.Pressure] = false;
            XOP.IsLinear = false;
        }

        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            XNSE_OperatorConfiguration config = new XNSE_OperatorConfiguration(this.Control);

            for(int d = 0; d < D; ++d) {
                opFactory.AddEquation(new Boundary.NavierCauchyBoundary("A", "C", d, D, Control.Material, config.physParams.rho_A, config.physParams.mu_A));
                opFactory.AddEquation(new Boundary.NavierCauchyBoundary("B", "C", d, D, Control.Material, config.physParams.rho_B, config.physParams.mu_B));
                //opFactory.AddEquation(new DisplacementBoundary(LsTrk, "A", "C", d, D, Control.ArtificialViscosity, config.physParams.mu_A, Control.Material));
                //opFactory.AddEquation(new DisplacementBoundary(LsTrk, "B", "C", d, D, Control.ArtificialViscosity, config.physParams.mu_B, Control.Material));
            }


            opFactory.AddEquation(new ContinuityBoundary("A", "C", D));
            opFactory.AddEquation(new ContinuityBoundary("B", "C", D));


            //*
            if(config.dntParams.SST_isotropicMode == BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                for(int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new EquilibriumContactLine(d, D, config.physParams.betaL, config.physParams.theta_e));
                    //opFactory.AddEquation(new EquilibriumContactLine1(d, D, config.physParams.betaL, config.physParams.theta_e));
                }
            }
            //*/
            var normalsParameter = new BoSSS.Solution.XNSECommon.Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[1]).Basis.Degree, VariableNames.SolidLevelSetCG);
            opFactory.AddParameter(normalsParameter);
            lsUpdater.AddLevelSetParameter(VariableNames.SolidLevelSetCG, normalsParameter);
        }

        internal bool LastSolverSuccess;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetTimestep();
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt}");
            LastSolverSuccess = Timestepping.Solve(phystime, dt, this.Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"done with time step {TimestepNo}, Solver success? {LastSolverSuccess}");
            Assert.IsTrue(LastSolverSuccess, "Solver did not converge");

            //OperatorAnalysis();
            //base.TerminationKey = true;
            //Debugger.Launch();
            {
                var Phi0 = this.LsUpdater.LevelSets.ElementAt(0).Value.C0LevelSet;
                var Phi1 = this.LsUpdater.LevelSets.ElementAt(1).Value.C0LevelSet;

                Console.WriteLine("  Len of Phi0: " + Phi0.CoordinateVector.Mapping.GlobalCount);
                Console.WriteLine("  Len of Phi1: " + Phi1.CoordinateVector.Mapping.GlobalCount);

                var LsChecker = new TestingIO(this.GridData, $"LevelSets-{TimestepNo}.csv", TestingIO.DataCorrelation.GeometricalCode, 1);
                LsChecker.AddDGField(Phi0);
                LsChecker.AddDGField(Phi1);

                var Species = new[] { "A", "B", "C" };
                int order = this.LsTrk.GetCachedOrders().Max();
                var metrics = this.LsTrk.GetXDGSpaceMetrics(Species, order).CutCellMetrics;

                int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;

                var cutLineViz = new List<ConventionalDGField>();

                foreach(string spc in Species) {
                    var spcId = LsTrk.GetSpeciesId(spc);
                    var vol = metrics.CutCellVolumes[spcId].To1DArray().GetSubVector(0, J);
                    var LsArea = metrics.InterfaceArea[spcId].To1DArray().GetSubVector(0, J);
                    var BndyArea = metrics.CellSurface[spcId].To1DArray().GetSubVector(0, J);

                    var cutLine = metrics.CutLineLength[spcId].To1DArray().GetSubVector(0, J);
                    var intersectLine = metrics.IntersectionLength[spcId].To1DArray().GetSubVector(0, J);

                    SinglePhaseField cutLineDG = new SinglePhaseField(new Basis(this.Grid, 0), "cutline-Spc#" + spc);
                    cutLineDG.CoordinateVector.SetV(cutLine);
                    cutLineViz.Add(cutLineDG);

                    LsChecker.AddVector("Volume-" + spc, vol);
                    LsChecker.AddVector("Interface-" + spc, LsArea);
                    LsChecker.AddVector("CellBndy-" + spc, BndyArea);
                    LsChecker.AddVector("cutLine-" + spc, cutLine);
                    LsChecker.AddVector("intersectLine-" + spc, intersectLine);
                    LsChecker.AddDGField(cutLineDG);
                }
                LsChecker.AddDGField(this.Velocity[0]);
                LsChecker.AddDGField(this.Velocity[1]);

                LsChecker.DoIOnow();

                foreach(string spc in Species) {
                    var name = "cutline-Spc#" + spc;
                    var err = LsChecker.LocalError(cutLineViz.Single(f => f.Identification == name));
                    cutLineViz.Add(err);

                    var refValue = cutLineViz.Single(f => f.Identification == name).CloneAs();
                    LsChecker.OverwriteDGField(refValue);
                    refValue.Identification = "Reference-" + refValue.Identification;
                    cutLineViz.Add(refValue);
                }


                Assert.Less(LsChecker.AbsError(Phi0), 1.0e-8, "Mismatch in level-set 0 between single-core and parallel run.");
                Assert.Less(LsChecker.AbsError(Phi1), 1.0e-8, "Mismatch in level-set 1 between single-core and parallel run.");

                foreach(string spc in Species) {
                    Console.WriteLine("    Volume comparison error    : " + LsChecker.AbsError("Volume-" + spc));
                    Console.WriteLine("    Surface comparison error   : " + LsChecker.AbsError("Interface-" + spc));
                    Console.WriteLine("    Cell surf comparison error : " + LsChecker.AbsError("CellBndy-" + spc));
                    Console.WriteLine("    Cut line comparison error  : " + LsChecker.AbsError("cutLine-" + spc));
                    Console.WriteLine("    Intersection error         : " + LsChecker.AbsError("intersectLine-" + spc));
                }

                Console.WriteLine("    VelocityX comparison error     : " + LsChecker.LocalError(this.Velocity[0]).L2NormAllSpecies());
                Console.WriteLine("    VelocityX comparison error     : " + LsChecker.LocalError(this.Velocity[1]).L2NormAllSpecies());


                Tecplot.PlotFields(cutLineViz, "cutLineMeasures", 0, 0);
            }


            return dt;
        }

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;

            if(iLevSet == 0) {
                string[] lsVelocityName = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(
                    NSEVariableNames.LevelSetCG, NSEVariableNames.VelocityVector(D)).ToArray();
                ILevelSetParameter levelSetVelocity = new SinglePhaseFieldVariableCopy("A", NSEVariableNames.VelocityVector(D), lsVelocityName);
                return levelSetVelocity;
            } else if(iLevSet == 1) {
                string[] lsVelocityName = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(
                    VariableNames.SolidLevelSetCG, NSEVariableNames.VelocityVector(D)).ToArray();
                ILevelSetParameter levelSetVelocity = new SinglePhaseFieldVariableCopy("C" ,NSEVariableNames.VelocityVector(D), lsVelocityName);
                return levelSetVelocity;
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// automatized analysis of condition number 
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis(OperatorAnalysisConfig config) {
            int D = this.Grid.SpatialDimension;

            int[] varGroup_all = (2 * D + 1).ForLoop(i => i);

            int[][] groups = new[] {
                varGroup_all
            };
            

            var res = this.Timestepping.OperatorAnalysis(config, groups);

            // filter only those results that we want;
            // this is a DG app, but it uses the LevelSetTracker; therefore, we want to filter analysis results for cut cells and only return uncut cells resutls
            var ret = new Dictionary<string, double>();
            foreach(var kv in res) {
                if(kv.Key.ToLowerInvariant().Contains("innercut") || kv.Key.ToLowerInvariant().Contains("bndycut")) {
                    // ignore
                } else {
                    ret.Add(kv.Key, kv.Value);
                }
            }

            return ret;
        }
    }
}

