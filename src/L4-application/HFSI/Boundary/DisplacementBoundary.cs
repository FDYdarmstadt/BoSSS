using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using HFSISolver.ContactLine;
using HFSISolver.SolidPhase;

namespace HFSISolver.Boundary {
    
    /// <summary>
    /// Solid Fluid Interface for Displacement Transport Equation
    /// </summary>
    class DisplacementBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;


        public DisplacementBoundary(LevelSetTracker LsTrkr, string fluidSpecies, string solidSpecies, int d, int D, double artificialViscosity, double fluidViscosity, Solid material) {
            codomainName = HFSISolver.EquationNames.DisplacementEvolutionComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(HFSISolver.VariableNames.DisplacementVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            
            AddComponent(new InnerNonLinearSolidConvectionForm(HFSISolver.VariableNames.DisplacementVector(D), 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), 
                1.0, d, 1, 
                fluidSpecies, 
                solidSpecies));
            //*/
            /*
            AddComponent(new ParameterTransportBoundaryForm( HFSISolver.VariableNames.Displacement0Vector(D)[d],
                BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D),
                1.0, D,
                fluidSpecies,
                solidSpecies,
                1
                ));
            AddParameter(HFSISolver.VariableNames.Displacement0Vector(D)[d]);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            
            //AddComponent(new ConvectionDivergenceBoundaryForm(d, D, 1, fluidSpecies, solidSpecies, 1000));

            if (artificialViscosity != 0.0)
            {
                //AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, d, D, 1, artificialViscosity, fluidViscosity, material.Viscosity));
                AddComponent(new ParameterSIPBoundaryForm(fluidSpecies, solidSpecies, d, D, 1, HFSISolver.VariableNames.ArtificialViscosity));
            }
            //AddSurfaceComponent(new SurfacePenaltyForm(HFSISolver.VariableNames.DisplacementComponent(d), -100));
            */
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}
