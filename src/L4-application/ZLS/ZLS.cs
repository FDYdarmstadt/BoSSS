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

namespace ZwoLevelSetSolver {

    public class ZLS : XNSE<ZLS_Control> {
        
        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            base.AddMultigridConfigLevel(configsLevel, iLevel);
            
            int pVel = VelocityDegree();
            int D = this.GridData.SpatialDimension;
            for(int d = 0; d < D; d++) {
                var configDisplacement = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pVel },
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.DisplacementVector(D)[d]) }
                };
                configsLevel.Add(configDisplacement);
            }
        }

        protected override void SetInitial(double t) {
            base.SetInitial(t);
            BoSSS.Solution.Application.DeleteOldPlotFiles();
        }

        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystem(D, opFactory, lsUpdater);
            DefineSolidPhase(D, opFactory, lsUpdater);
        }

        void DefineSolidPhase(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            
            for(int d = 0; d < D; ++d) {
                if(this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new LinearNavierCauchy("C", Control.Material, d, D));
                    opFactory.AddEquation(new LinearDisplacementEvolution("C", d, D));
                }else if(this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    opFactory.AddEquation(new NavierCauchy("C", Control.Material, d, D));
                    opFactory.AddEquation(new DisplacementEvolution("C", d, D));
                } else {
                    throw new NotSupportedException();
                }

                opFactory.AddEquation(new Dummy("A", VariableNames.DisplacementVector(D)[d], EquationNames.DisplacementEvolutionComponent(d)));
                opFactory.AddEquation(new Dummy("B", VariableNames.DisplacementVector(D)[d], EquationNames.DisplacementEvolutionComponent(d)));
            }
            opFactory.AddEquation(new SolidPhase.Continuity("C", D));
        }

        protected override void FinalOperatorSettings(XSpatialOperatorMk2 XOP) {
            base.FinalOperatorSettings(XOP);
            XOP.IsLinear = false;
        }

        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            //*
            for(int d = 0; d < D; ++d) {
                if(this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new LinearDisplacementBoundary(LsTrk, "A", "C", d, D));
                    opFactory.AddEquation(new LinearDisplacementBoundary(LsTrk, "B", "C", d, D));
                    opFactory.AddEquation(new LinearNavierCauchyBoundary("A", "C", d, D, Control.Material, config.physParams.rho_A, config.physParams.mu_A));
                    opFactory.AddEquation(new LinearNavierCauchyBoundary("B", "C", d, D, Control.Material, config.physParams.rho_A, config.physParams.mu_B));
                } else if(this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    opFactory.AddEquation(new DisplacementBoundary(LsTrk, "A", "C", d, D));
                    opFactory.AddEquation(new DisplacementBoundary(LsTrk, "B", "C", d, D));
                    opFactory.AddEquation(new NavierCauchyBoundary("A", "C", d, D, Control.Material, config.physParams.rho_A, config.physParams.mu_A));
                    opFactory.AddEquation(new NavierCauchyBoundary("B", "C", d, D, Control.Material, config.physParams.rho_A, config.physParams.mu_B));
                } else {
                    throw new NotSupportedException();
                }
            }

            if(this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                int quadOrder = QuadOrder();
                Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
                Velocity0 v0 = new Velocity0(D);
                opFactory.AddParameter(v0);
                opFactory.AddParameter(v0Mean);
                lsUpdater.AddLevelSetParameter(VariableNames.SolidLevelSetCG, v0Mean);
            }

            //*
            opFactory.AddEquation(new FluidSolidContinuity("A", "C", D));
            opFactory.AddEquation(new FluidSolidContinuity("B", "C", D));
            //*/

            /*
            for(int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NSEFreeSlipBoundary("A", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                opFactory.AddEquation(new NSEFreeSlipBoundary("B", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
            }

            opFactory.AddEquation(new ImmersedBoundaryContinuity("A", "C", 1, config, D));
            opFactory.AddEquation(new ImmersedBoundaryContinuity("B", "C", 1, config, D));

            //throw new NotImplementedException("todo");
            opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));
            */
            //*
            if(config.dntParams.SST_isotropicMode == BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                for(int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new SlipContactLine(d, D, config.physParams.betaL, config.physParams.theta_e));
                }
                //ContactLine
                //=====================
                var normalsParameter = new BoSSS.Solution.XNSECommon.Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[1]).Basis.Degree, VariableNames.SolidLevelSetCG);
                opFactory.AddParameter(normalsParameter);
                lsUpdater.AddLevelSetParameter(VariableNames.SolidLevelSetCG, normalsParameter);
            }
            //*/
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetFixedTimestep();
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt}");
            Timestepping.Solve(phystime, dt, false);
            Console.WriteLine($"done with time step {TimestepNo}");
            return dt;
        }

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;

            if(iLevSet == 0) {
                //ILevelSetParameter levelSetVelocity = new LevelSetVelocity(NSEVariableNames.LevelSetCG, D, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                ILevelSetParameter levelSetVelocity = new SinglePhaseFieldVariableCopy(NSEVariableNames.LevelSetCG, "A", NSEVariableNames.VelocityVector(D));
                return levelSetVelocity;
            } else if(iLevSet == 1) {
                ILevelSetParameter levelSetVelocity = new SinglePhaseFieldVariableCopy(VariableNames.SolidLevelSetCG, "C" ,NSEVariableNames.VelocityVector(D));
                return levelSetVelocity;
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }
    }
}
