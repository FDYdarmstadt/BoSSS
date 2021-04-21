using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XNSECommon {


    /// <summary>
    /// Used for setting the field everywhere equal to one in the system of equations.
    /// </summary>
    public class identityTerm : Idsource, ISpeciesFilter {
        public identityTerm(string spcsName, string varname) : base(varname) {
            ValidSpecies = spcsName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    public class Idsource : LinearSource {
        public Idsource(string _var) {
            m_var = _var;
        }

        string m_var;

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { m_var };
            }
        }

        protected override double Source(double[] x, double[] parameters, double[] U) {
            return U[0];
        }
    }



    /// <summary>
    /// 
    /// </summary>
    public class LowMach_Gravity : BuoyancyJacobi, ISpeciesFilter {
        public LowMach_Gravity(string spcsName, Vector GravityDirection, int SpatialComponent, double Froude, PhysicsMode physicsMode, MaterialLaw EoS, int noOfChemComponents) : base(GravityDirection, SpatialComponent, Froude, physicsMode, EoS, noOfChemComponents) {
            ValidSpecies = spcsName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    /// <summary>
    /// Heat release term on energy equation
    /// </summary>
    public class LowMach_HeatSource : ReactionHeatSourceJacobi, ISpeciesFilter {
        public LowMach_HeatSource(string spcName, double HeatReleaseFactor, double[] ReactionRateConstants, double[] molarmasses, MaterialLaw EoS, double TRef, double cpRef, bool VariableOneStepParameters) : base(HeatReleaseFactor, ReactionRateConstants, molarmasses, EoS, TRef, cpRef, VariableOneStepParameters) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }



    /// <summary>
    /// Reaction term on mass fraction equation
    /// </summary>
    public class LowMach_MassFractionSource : ReactionSpeciesSourceJacobi, ISpeciesFilter {
        public LowMach_MassFractionSource(string spcName, double[] ReactionRateConstants, double[] StoichiometricCoefficients, double[] MolarMasses, MaterialLaw EoS, int NumberOfReactants, int SpeciesIndex, double TRef, double cpRef, bool VariableOneStepParameters) : base(ReactionRateConstants, StoichiometricCoefficients, MolarMasses, EoS, NumberOfReactants, SpeciesIndex, TRef, cpRef, VariableOneStepParameters) {

            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    /// <summary>
    /// Manufactured solution 
    /// </summary>
    public class LowMach_ContiManSolution : RHSManuSourceDivKonti, ISpeciesFilter {
        public LowMach_ContiManSolution(string spcName, double Reynolds, double[] MolarMasses, PhysicsMode physicsMode, bool rhoOne, Func<double[], double, double> _sourceFunc = null) : base(Reynolds, MolarMasses, physicsMode, rhoOne, _sourceFunc) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }
    /// <summary>
    /// Manufactured solution 
    /// </summary>
    public class LowMach_MomentumManSolution : RHSManuSourceNS, ISpeciesFilter {
        public LowMach_MomentumManSolution(string spcName, double Reynolds, double Froude, double[] MolarMasses, string direction, PhysicsMode physMode, bool rhoOne, Func<double[], double, double> _SourceTerm) : base(Reynolds, Froude, MolarMasses, direction, physMode, rhoOne, _SourceTerm) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }

    }
    /// <summary>
    /// Manufactured solution 
    /// </summary>
    public class LowMach_ScalarManSolution : RHSManuSourceTransportEq, ISpeciesFilter {
        public LowMach_ScalarManSolution(string spcName, double HeatRelease, double Reynolds, double Prandtl, double Schmidt, double[] StoichiometricCoefficients, double[] ReactionRateConstants, double[] MolarMasses, MaterialLaw EoS, string EqType, PhysicsMode physicsMode, int SpeciesIndex = -1, bool chemReactionOK = true, bool rhoOne = false) : base(HeatRelease, Reynolds, Prandtl, Schmidt, StoichiometricCoefficients, ReactionRateConstants, MolarMasses, EoS, EqType, physicsMode, SpeciesIndex, chemReactionOK, rhoOne) {

            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }

    /// <summary>
    /// Time derivative contribution in the continuity equation.
    /// </summary>
    public class LowMach_TimeDerivativeConti_term1 : ContiTimeDerivativeActual, ISpeciesFilter {
        public LowMach_TimeDerivativeConti_term1(string spcName, MaterialLaw EoS, double dt, int NumberOfChemicalComponents) : base(EoS, dt, NumberOfChemicalComponents) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    /// <summary>
    /// implementation of the term rho^n_(t)/dt
    /// </summary>
    public class ContiTimeDerivativeActual : IVolumeForm, ISupportsJacobianComponent, IEquationComponentCoefficient {

        double m_dt;
        string[] m_ArgumentOrdering;
        MaterialLaw m_EoS;

        public ContiTimeDerivativeActual(MaterialLaw EoS, double dt, int NumberOfChemicalComponents) {
            m_EoS = EoS;
            m_dt = dt;
            m_ArgumentOrdering = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(NumberOfChemicalComponents)); // Variables for the density evaluation
        }
        /// <summary>
        /// CHECK
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                //return (TermActivationFlags.V | TermActivationFlags.UxV);
                return TermActivationFlags.AllOn;
            }
        }

        /// <summary>
        ///  
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }
        /// <summary>
        ///  
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            // For now no update. In future dt could be dynamicaly changed here
            return;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] {
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double rho = m_EoS.GetDensity(U);
            return rho / m_dt * V;

        }
    }


    /// <summary>
    /// Time derivative contribution in the continuity equation.
    /// </summary>
    public class LowMach_TimeDerivativeConti_term2 : ContiTimeDerivativePrevious, ISpeciesFilter {
        public LowMach_TimeDerivativeConti_term2(string spcName, MaterialLaw EoS, double dt, int NumberOfChemicalComponents) : base(EoS, dt, NumberOfChemicalComponents) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    /// <summary>
    /// implementation of the term rho^n_(t)/dt
    /// </summary>
    public class ContiTimeDerivativePrevious : IVolumeForm, ISupportsJacobianComponent, IEquationComponentCoefficient {

        double m_dt;
        string[] m_ParameterOrdering;
        MaterialLaw m_EoS;

        public ContiTimeDerivativePrevious(MaterialLaw EoS, double dt, int NumberOfChemicalComponents) {
            m_EoS = EoS;
            m_dt = dt;

            m_ParameterOrdering = ArrayTools.Cat(new string[] { VariableNames.Temperature + "_t0" }, VariableNames.MassFractions_t0(NumberOfChemicalComponents)); // Variables for the density evaluation
        }
        /// <summary>
        /// CHECK
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.AllOn;
            }
        }

        /// <summary>
        ///  
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }
        /// <summary>
        ///  
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            // For now no update. In future dt could be dynamicaly changed here
            return;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] {
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double[] densityParameters = cpv.Parameters;
            double rho = m_EoS.GetDensity(densityParameters);
            return -1 * rho / m_dt * V;

        }
    }
}
