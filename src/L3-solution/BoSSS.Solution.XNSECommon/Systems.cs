using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XNSECommon {
    
    public interface ISystem {
        void DefineSystem( OperatorFactory operatorFactory, LevelSetUpdater levelSetUpdater);
    }

    public class XNSESystem : ISystem {

        int velocityDegree;

        IncompressibleMultiphaseBoundaryCondMap boundaryMap;

        IXNSE_Configuration config;

        AppControl control;

        int QuadOrder;
        
        int D;

        public XNSESystem(int D, int quadOrder, IXNSE_Configuration config, AppControl Control, IncompressibleMultiphaseBoundaryCondMap boundaryMap) {
            this.QuadOrder = quadOrder;
            this.velocityDegree = VelocityDegree(Control);
            this.boundaryMap = boundaryMap;
            this.D = D;
            this.control = Control;
            this.config = config;
        }

        static int VelocityDegree(AppControl Control) {
            int pVel;
            if (Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v)) {
                pVel = v.Degree;
            } else if (Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1)) {
                pVel = v1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        public void DefineSystem(OperatorFactory OpFactory, LevelSetUpdater LsUpdater) {
            LevelSetTracker LsTrk = LsUpdater.Tracker;
            for (int d = 0; d < D; ++d) {
                OpFactory.AddEquation(new NavierStokes("A", d, LsTrk, D, boundaryMap, config));
                OpFactory.AddParameter(Gravity.CreateFrom("A", d, D, control, config.getPhysParams.rho_A));
                OpFactory.AddEquation(new NavierStokes("B", d, LsTrk, D, boundaryMap, config));
                OpFactory.AddParameter(Gravity.CreateFrom("B", d, D, control, config.getPhysParams.rho_B));
                OpFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                OpFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, LsTrk, config));
            }
            OpFactory.AddCoefficient(new SlipLengths(config, velocityDegree));
            Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, QuadOrder);
            if (config.getPhysParams.IncludeConvection && config.isTransport) {
                OpFactory.AddParameter(new Velocity0(D));
                OpFactory.AddParameter(v0Mean);
            }
            int lsDegree = LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet.Basis.Degree;
            Normals normalsParameter = new Normals(D, lsDegree);
            OpFactory.AddParameter(normalsParameter);

            if (config.isContinuity) {
                OpFactory.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
                OpFactory.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
                OpFactory.AddEquation(new InterfaceContinuity(config, D, LsTrk, config.isMatInt));
            }

            LsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, v0Mean);
            LsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, normalsParameter);

            switch (config.getDntParams.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                MaxSigma maxSigmaParameter = new MaxSigma(config.getPhysParams, config.getDntParams, QuadOrder, control.dtFixed);
                OpFactory.AddParameter(maxSigmaParameter);
                LsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, maxSigmaParameter);
                BeltramiGradient lsBGradient = new BeltramiGradient( D, config.getDntParams, lsDegree);
                LsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsBGradient);
                break;

                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                BeltramiGradient lsGradient = new BeltramiGradient(D, config.getDntParams, lsDegree);
                LsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradient);
                break;

                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                int curvatureDegree;
                if (control.FieldOptions.TryGetValue(VariableNames.Curvature, out FieldOpts opts)) {
                    curvatureDegree = opts.Degree;
                } else {
                    throw new Exception("Curvature options not found in FieldOptions");
                }
                BeltramiGradientAndCurvature lsGradientAndCurvature = new BeltramiGradientAndCurvature(curvatureDegree, lsDegree, QuadOrder, config.getDntParams, D);
                OpFactory.AddParameter(lsGradientAndCurvature);
                LsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradientAndCurvature);
                break;

                default:
                throw new NotImplementedException($"option {config.getDntParams.SST_isotropicMode} is not implemented.");
            }
        }
    }
}
