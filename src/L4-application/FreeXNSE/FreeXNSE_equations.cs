using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace FreeXNSE {
    internal class BulkContinuity : BulkEquation {
        string speciesName;

        string codomainName;

        double massScale;

        public BulkContinuity(
            int D,
            string[] spcName,
            FreeXNSE_BoundaryCondMap BcMap,
            FreeXNSE_Control Control) {

            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            speciesName = spcName[0];
            massScale = 0;

            // Velocity Divergence
            switch(Control.ActiveTerms.VelocityDivergence) {
                case VelocityDivergence.Off: {
                    break;
                }
                default:
                case VelocityDivergence.Central: {
                    AddComponent(new BulkComponent_VelocityDivergence_Central(D, spcName, BcMap));
                    break;
                }
            }
            if(Control.EqualOrder) {
                if(DimensionlessNumbers.IsPhysical(Control.DimensionlessNumbers.Re)) {
                    AddComponent(new BulkComponent_PressurePenalty(D, spcName, BcMap, Control.DimensionlessNumbers.Re));
                }
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => massScale;

        public override string CodomainName => codomainName;
    }

    internal class InterfaceContinuity : SurfaceEquation {
        string codomainName;
        string phaseA, phaseB;

        //Methode aus der XNSF_OperatorFactory
        public InterfaceContinuity(
            int D,
            string[] spcName,
            FreeXNSE_BoundaryCondMap BcMap,
            FreeXNSE_Control Control) {

            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));

            this.phaseA = spcName[0];
            this.phaseB = spcName[1];

            /// Velocity Divergence
            switch(Control.ActiveTerms.VelocityDivergence) {
                case VelocityDivergence.Off: {
                    break;
                }
                default:
                case VelocityDivergence.Central: {
                    AddComponent(new InterfaceComponent_VelocityDivergence_Central(D, spcName, BcMap));
                    break;
                }
            }
            if(Control.EqualOrder) {
                // nothing to do
            }
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    internal class BulkMomentum : BulkEquation {
        string speciesName;

        string codomainName;

        double massScale;

        public BulkMomentum(
            int d,
            int D,
            string[] spcName,
            FreeXNSE_BoundaryCondMap BcMap,
            FreeXNSE_Control Control) {
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            codomainName = EquationNames.MomentumEquationComponent(d);
            speciesName = spcName[0];

            // Temporal
            switch(Control.ActiveTerms.Temporal) {
                case Temporal.Off: {
                    massScale = 0;
                    break;
                }
                default:
                case Temporal.On: {
                    massScale = 1;
                    break;
                }
            }
            // Convective
            switch(Control.ActiveTerms.Convective) {
                case Convective.Off: {
                    break;
                }
                default:
                case Convective.LaxFriedrich: {
                    AddComponent(new BulkComponent_Convective_LaxFriedrich(d, D, spcName, BcMap));
                    AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
                    break;
                }
                case Convective.Temam: {
                    AddComponent(new BulkComponent_Convective_Temam(d, D, spcName, BcMap));
                    break;
                }
                case Convective.ConservativeTemam: {
                    Console.WriteLine("Not sure if ConservativeTemam is implemented correctly!");
                    AddComponent(new BulkComponent_Convective_ConservativeTemam(d, D, spcName, BcMap));
                    break;
                }
            }
            // Pressure Gradient
            switch(Control.ActiveTerms.PressureGradient) {
                case PressureGradient.Off: {
                    break;
                }
                default:
                case PressureGradient.Central: {
                    AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
                    AddComponent(new BulkComponent_PressureGradient_Central(d, D, spcName, BcMap));                    
                    break;
                }
            }
            // Viscous
            if(DimensionlessNumbers.IsPhysical(Control.DimensionlessNumbers.Re)) {
                switch(Control.ActiveTerms.Viscous) {
                    case Viscous.Off: {
                        break;
                    }
                    default:
                    case Viscous.SIP: {
                        AddComponent(new BulkComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, 1, 4));
                        AddCoefficient(Coefficientnames.Bulkfriction);
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                    case Viscous.NIP: {
                        AddComponent(new BulkComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, -1, 4));
                        AddCoefficient(Coefficientnames.Bulkfriction);
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                    case Viscous.IIP: {
                        AddComponent(new BulkComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, 0, 4));
                        AddCoefficient(Coefficientnames.Bulkfriction);
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                    case Viscous.OBB: {
                        AddComponent(new BulkComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, -1, 0));
                        AddCoefficient(Coefficientnames.Bulkfriction);
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                }                
            }

            if(Control.VolumeForce != null) {
                if(DimensionlessNumbers.IsPhysical(Control.DimensionlessNumbers.Fr))
                    AddComponent(new BulkComponent_VolumeForce(d, Control.DimensionlessNumbers.Fr, Control.VolumeForce, D, spcName, BcMap));
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => massScale;

        public override string CodomainName => codomainName;
    }

    internal class InterfaceMomentum : SurfaceEquation {
        string codomainName;
        string phaseA, phaseB;

        //Methode aus der XNSF_OperatorFactory
        public InterfaceMomentum(
            int d,
            int D,
            string[] spcName,
            FreeXNSE_BoundaryCondMap BcMap,
            FreeXNSE_Control Control) {

            codomainName = EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));

            this.phaseA = spcName[0];
            this.phaseB = spcName[1];

            // Convective
            switch(Control.ActiveTerms.Convective) {
                case Convective.Off: {
                    break;
                }
                default:
                case Convective.LaxFriedrich: {
                    AddComponent(new InterfaceComponent_Convective_LaxFriedrich(d, D, spcName, BcMap));
                    break;
                }
                case Convective.Temam: {
                    AddComponent(new InterfaceComponent_Convective_Temam(d, D, spcName, BcMap));
                    break;
                }
                case Convective.ConservativeTemam: {
                    AddComponent(new InterfaceComponent_Convective_ConservativeTemam(d, D, spcName, BcMap));
                    break;
                }
            }
            // Pressure Gradient
            switch(Control.ActiveTerms.PressureGradient) {
                case PressureGradient.Off: {
                    break;
                }
                default:
                case PressureGradient.Central: {
                    AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
                    AddComponent(new InterfaceComponent_PressureGradient_Central(d, D, spcName, BcMap));
                    break;
                }
            }
            // Surface Tension
            if(DimensionlessNumbers.IsPhysical(Control.DimensionlessNumbers.We)) {
                switch(Control.ActiveTerms.SurfaceTension) {
                    case SurfaceTension.Off: {
                        break;
                    }
                    default:
                    case SurfaceTension.LaplaceBeltrami: {
                        AddSurfaceComponent(new InterfaceComponent_SurfaceTension_LaplaceBeltrami(d, D, spcName, BcMap, Control.DimensionlessNumbers.We));
                        AddParameter(BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[d]);
                        AddCoefficient(Coefficientnames.Contactangle);
                        AddCoefficient(Coefficientnames.Contactangle + "Adv");
                        AddCoefficient(Coefficientnames.Contactangle + "Rec");
                        AddCoefficient(Coefficientnames.Contactlinefriction);
                        AddCoefficient(Coefficientnames.Surfacetensionfield);
                        break;
                    }
                    case SurfaceTension.LaplaceBeltrami_BoussinesqScriven: {
                        AddSurfaceComponent(new InterfaceComponent_SurfaceTension_LaplaceBeltrami_BoussinesqScriven(d, D, spcName, BcMap, Control.DimensionlessNumbers.We));
                        AddParameter(BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[d]);
                        AddCoefficient(Coefficientnames.Contactangle);
                        AddCoefficient(Coefficientnames.Contactlinefriction);
                        AddCoefficient(Coefficientnames.Surfacetensionfield);
                        AddCoefficient(Coefficientnames.Surfaceshearviscosityfield);
                        AddCoefficient(Coefficientnames.Surfacedilatationalviscosityfield);
                        break;
                    }
                }
            }

            // Viscous
            if(DimensionlessNumbers.IsPhysical(Control.DimensionlessNumbers.Re)) { 
                switch(Control.ActiveTerms.Viscous) {
                    case Viscous.Off: {
                        break;
                    }
                    default:
                    case Viscous.SIP: {
                        AddComponent(new InterfaceComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, 1, 4));
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                    case Viscous.NIP: {
                        AddComponent(new InterfaceComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, -1, 4));
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                    case Viscous.IIP: {
                        AddComponent(new InterfaceComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, 0, 4));
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                    case Viscous.OBB: {
                        AddComponent(new InterfaceComponent_Viscous_IP(d, D, spcName, BcMap, Control.DimensionlessNumbers.Re, -1, 0));
                        AddCoefficient(Coefficientnames.Bulkviscosityfield);
                        break;
                    }
                }
            }
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    internal class BulkDummy : BulkEquation {
        string spc;
        public override string SpeciesName => spc;

        public override double MassScale => 0.0;
        string codom;
        public override string CodomainName => codom;

        internal BulkDummy(string spc, string codom, string dom) { 
            this.spc = spc;
            this.codom = codom;
            AddVariableNames(dom);
            AddComponent(new BulkComponent_Dummy(SpeciesName, dom));
        }
    }

}
