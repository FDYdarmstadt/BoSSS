using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using System;
using System.Linq;
namespace BoSSS.Application.XNSERO_Solver {
    class Equations {

        /// <summary>
        /// Incompressible, Newtonian momentum equation, (fluid/solid) immersed boundary
        /// </summary>
        public class NSEROimmersedBoundary : SurfaceEquation {
            string m_codomainName;
            string m_fluidPhase;
            string m_solidPhase;
            int m_iLevSet;

            //Methode aus der XNSF_OperatorFactory
            public NSEROimmersedBoundary(
                string fluidPhase,
                string solidPhase,
                int iLevSet,
                int d,
                int D,
                IncompressibleMultiphaseBoundaryCondMap boundaryMap,
                LevelSetTracker LsTrk,
                INSE_Configuration config,
                bool isMovingMesh) : base() //
            {
                m_fluidPhase = fluidPhase;
                m_solidPhase = solidPhase;
                m_iLevSet = iLevSet;
                m_codomainName = EquationNames.MomentumEquationComponent(d);
                AddInterfaceNSE(D, d, boundaryMap, LsTrk, config, isMovingMesh);
                AddVariableNames(Solution.NSECommon.VariableNames.VelocityVector(D).Cat(Solution.NSECommon.VariableNames.Pressure));
                AddParameter(Solution.NSECommon.VariableNames.AsLevelSetVariable(Solution.NSECommon.VariableNames.LevelSetCGidx(m_iLevSet), Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray());
                AddParameter(Solution.NSECommon.VariableNames.AsLevelSetVariable(Solution.NSECommon.VariableNames.LevelSetCGidx(m_iLevSet), Solution.NSECommon.VariableNames.SurfaceForceVector(D)).ToArray());
            }


            void AddInterfaceNSE(
                int D,
                int d,
                IncompressibleMultiphaseBoundaryCondMap boundaryMap,
                LevelSetTracker LsTrk,
                INSE_Configuration config,
                bool isMovingMesh) {
                PhysicalParameters physParams = config.getPhysParams;
                DoNotTouchParameters dntParams = config.getDntParams;

                // set species arguments
                double rho, LFF, mu;
                switch(this.m_fluidPhase) {
                    case "A":
                    rho = physParams.rho_A;
                    LFF = dntParams.LFFA;
                    mu = physParams.mu_A;
                    break;

                    case "B":
                    rho = physParams.rho_B;
                    LFF = dntParams.LFFB;
                    mu = physParams.mu_B;
                    break;

                    default: throw new NotSupportedException($"Unknown fluid species: {this.m_fluidPhase}");
                }

                // convective operator
                // ===================
                if(physParams.IncludeConvection && config.isTransport) {
                    var ConvIB = new ConvectionAtIB( d, D, LsTrk, LFF, boundaryMap, rho, isMovingMesh, m_iLevSet, m_fluidPhase, m_solidPhase, true);
                    AddComponent(ConvIB);
                }
                if(isMovingMesh && (physParams.IncludeConvection && config.isTransport == false)) {
                    // if Moving mesh, we need the interface transport term somehow
                    throw new NotImplementedException("Something missing here.");
                }

                // pressure gradient
                // =================
                if(config.isPressureGradient) {
                    var presLs = new BoSSS.Solution.NSECommon.Operator.Pressure.PressureFormAtIB(d, D, LsTrk, m_iLevSet, m_fluidPhase, m_solidPhase);
                    AddComponent(presLs);
                }

                // viscous operator
                // ================
                if(config.isViscous && (mu != 0.0)) {
                    double penalty = dntParams.PenaltySafety;
                    switch(dntParams.ViscosityMode) {
                        case ViscosityMode.Standard:
                        case ViscosityMode.TransposeTermMissing:
                        AddComponent(new ViscosityAtIB(d, D, LsTrk, penalty, mu, m_iLevSet, m_fluidPhase, m_solidPhase, true));
                        break;
                        case ViscosityMode.FullySymmetric:
                        case ViscosityMode.Viscoelastic:
                        throw new NotImplementedException("todo");
                        default:
                        throw new NotImplementedException();
                    }
                }
            }

            public override string FirstSpeciesName {
                get { return m_fluidPhase; }
            }

            public override string SecondSpeciesName {
                get { return m_solidPhase; }
            }

            public override string CodomainName {
                get {
                    return m_codomainName;
                }
            }
        }
    }
}
