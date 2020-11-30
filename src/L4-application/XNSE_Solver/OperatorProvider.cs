using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {    

    static class XHeatOperatorProvider{

        public static void SetOperatorEquations(int D, OperatorFactory opFactory, int quadOrder, ThermalMultiphaseBoundaryCondMap boundaryMap, LevelSetUpdater lsUpdater, XNSE_Control control, XNSFE_OperatorConfiguration config) {            

            // add Heat equation components
            // ============================
            opFactory.AddEquation(new Heat("A", lsUpdater.Tracker, D, boundaryMap, config));
            opFactory.AddEquation(new Heat("B", lsUpdater.Tracker, D, boundaryMap, config));

            if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                for (int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new HeatFlux("A", d, lsUpdater.Tracker, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFlux("B", d, lsUpdater.Tracker, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, boundaryMap, lsUpdater.Tracker, config));
                }
            }
            opFactory.AddEquation(new HeatInterface("A", "B", D, boundaryMap, lsUpdater.Tracker, config));      

        }

        public static void SetOperatorParameter(int D, OperatorFactory opFactory, int quadOrder, ThermalMultiphaseBoundaryCondMap boundaryMap, LevelSetUpdater lsUpdater, XNSE_Control Control, XNSFE_OperatorConfiguration config) {

            // Add Velocity parameters as prescribed variables
            for (int d = 0; d < D; d++)
                opFactory.AddParameter(Velocity0Prescribed.CreateFrom(lsUpdater.Tracker, d, D, Control));

            Velocity0MeanPrescribed v0Mean = new Velocity0MeanPrescribed(D, lsUpdater.Tracker, quadOrder);
            opFactory.AddParameter(v0Mean);
            lsUpdater.AddLevelSetParameter("Phi", v0Mean);

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);
            lsUpdater.AddLevelSetParameter("Phi", normalsParameter);
            
        }

    }

    static class XNSEOperatorProvider {

        public static void SetOperatorEquations(int D, OperatorFactory opFactory, int quadOrder, IncompressibleMultiphaseBoundaryCondMap boundaryMap, LevelSetUpdater lsUpdater, XNSE_Control Control, XNSFE_OperatorConfiguration config) {   

            ///Build Equations
            for (int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NavierStokes("A", d, lsUpdater.Tracker, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A));
                opFactory.AddEquation(new NavierStokes("B", d, lsUpdater.Tracker, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, lsUpdater.Tracker, config));
                opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, lsUpdater.Tracker, config));
            }          

            if (config.isContinuity) {
                opFactory.AddEquation(new Continuity(config, D, "A", lsUpdater.Tracker.GetSpeciesId("A"), boundaryMap));
                opFactory.AddEquation(new Continuity(config, D, "B", lsUpdater.Tracker.GetSpeciesId("B"), boundaryMap));
                opFactory.AddEquation(new InterfaceContinuity(config, D, lsUpdater.Tracker));
            }

        }

        public static void SetOperatorParameter(int D, OperatorFactory opFactory, int quadOrder, ThermalMultiphaseBoundaryCondMap boundaryMap, LevelSetUpdater lsUpdater, XNSE_Control Control, XNSFE_OperatorConfiguration config) {

            ///Build Equations
            for (int d = 0; d < D; ++d) {
                opFactory.AddParameter(Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A));
                opFactory.AddParameter(Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B));
            }

            opFactory.AddParameter(new Velocity0(D));
            Velocity0Mean v0Mean = new Velocity0Mean(D, lsUpdater.Tracker, quadOrder);
            opFactory.AddParameter(v0Mean);
            lsUpdater.AddLevelSetParameter("Phi", v0Mean);

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);
            lsUpdater.AddLevelSetParameter("Phi", normalsParameter);
        
        }

    }
}
