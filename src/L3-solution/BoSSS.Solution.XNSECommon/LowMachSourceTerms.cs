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
            return (U[0] - 1.0);
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
        public LowMach_HeatSource(string spcName, double HeatReleaseFactor, double[] ReactionRateConstants, double[] molarmasses, MaterialLaw_MultipleSpecies EoS, double TRef, double cpRef, bool VariableOneStepParameters, int NoOfChemSpecies) : base(HeatReleaseFactor, ReactionRateConstants, molarmasses, EoS, TRef, cpRef, VariableOneStepParameters, NoOfChemSpecies) {
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
        public LowMach_MassFractionSource(string spcName, double[] ReactionRateConstants, double[] StoichiometricCoefficients, double[] MolarMasses, MaterialLaw_MultipleSpecies EoS, int NumberOfReactants, int SpeciesIndex, double TRef, double cpRef, bool VariableOneStepParameters) : base(ReactionRateConstants, StoichiometricCoefficients, MolarMasses, EoS, NumberOfReactants, SpeciesIndex, TRef, cpRef, VariableOneStepParameters) {

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
    public class LowMach_ManSolution : RHSManuSource, ISpeciesFilter {
        public LowMach_ManSolution(string spcName, Func<double[], double, double> _sourceFunc ) : base(_sourceFunc) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    ///// <summary>
    ///// Manufactured solution 
    ///// </summary>
    //public class LowMach_ContiManSolution : RHSManuSourceDivKonti, ISpeciesFilter {
    //    public LowMach_ContiManSolution(string spcName, double Reynolds, double[] MolarMasses, PhysicsMode physicsMode, bool rhoOne, Func<double[], double, double> _sourceFunc = null) : base(Reynolds, MolarMasses, physicsMode, rhoOne, _sourceFunc) {
    //        ValidSpecies = spcName;
    //    }

    //    public string ValidSpecies {
    //        get;
    //        private set;
    //    }
    //}
    ///// <summary>
    ///// Manufactured solution 
    ///// </summary>
    //public class LowMach_MomentumManSolution : RHSManuSourceNS, ISpeciesFilter {
    //    public LowMach_MomentumManSolution(string spcName, double Reynolds, double Froude, double[] MolarMasses, string direction, PhysicsMode physMode, bool rhoOne, Func<double[], double, double> _SourceTerm) : base(Reynolds, Froude, MolarMasses, direction, physMode, rhoOne, _SourceTerm) {
    //        ValidSpecies = spcName;
    //    }

    //    public string ValidSpecies {
    //        get;
    //        private set;
    //    }

    //}
    ///// <summary>
    ///// Manufactured solution 
    ///// </summary>
    //public class LowMach_ScalarManSolution : RHSManuSourceTransportEq, ISpeciesFilter {
    //    public LowMach_ScalarManSolution(string spcName, double HeatRelease, double Reynolds, double Prandtl, double Schmidt, double[] StoichiometricCoefficients, double[] ReactionRateConstants, double[] MolarMasses, MaterialLaw EoS, string EqType, PhysicsMode physicsMode, int SpeciesIndex = -1, bool chemReactionOK = true, bool rhoOne = false) : base(HeatRelease, Reynolds, Prandtl, Schmidt, StoichiometricCoefficients, ReactionRateConstants, MolarMasses, EoS, EqType, physicsMode, SpeciesIndex, chemReactionOK, rhoOne) {

    //        ValidSpecies = spcName;
    //    }

    //    public string ValidSpecies {
    //        get;
    //        private set;
    //    }
    //}




    /// <summary>
    /// Time derivative contribution in the continuity equation.
    /// </summary>
    public class LowMach_TimeDerivativeConti : IVolumeForm, ISupportsJacobianComponent, IEquationComponentCoefficient, ISpeciesFilter {
        public LowMach_TimeDerivativeConti(string spcName, MaterialLaw EoS, double dt, int NumberOfChemicalComponents){
            ValidSpecies = spcName;
            m_EoS = EoS;
            m_dt = dt;

            if (m_EoS is MaterialLawMixtureFractionNew) {
                m_ArgumentOrdering = new string[] { VariableNames.MixtureFraction };

            } else {
                m_ArgumentOrdering = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(NumberOfChemicalComponents)); // Variables for the density evaluation

            }


            m_ParameterOrdering = new string[] { "Density_t00", "Density_t0", "Density" }; // Density is only added as a dummy  in order to be able to recognize it in the operator
        }

        public string ValidSpecies {
            get;
            private set;
        }


        double m_dt;
        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;
        MaterialLaw m_EoS;

 
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


            int timestepNo = (int) ( cpv.time / m_dt);

         

            double rho_t;
            if (m_EoS is MaterialLawMixtureFractionNew) {
                 rho_t = m_EoS.getDensityFromZ(U[0]);
            } else {
                rho_t = m_EoS.GetDensity(U);
            }
            double rho_t_00 = cpv.Parameters[0];
            double rho_t_0 = cpv.Parameters[1];

            double res = 0;
            if (timestepNo == 0) {
                throw new Exception("something went wrong");
            } else if (timestepNo == 1) {
              res = (rho_t - rho_t_0) / m_dt;  // implicit euler => drho/dt = (rho^t - rho^t-1) /2

            } else if (timestepNo == 2) {
                res = (3*rho_t - 4*rho_t_0+rho_t_00) /(2* m_dt);   //  2° order => drho/dt = (rho^t - rho^t-1) /2

            }


            return res* V;

        }
    }

 
     









    /// <summary>
    /// Thermodynamic pressure time derivative contribution in the energy equation.
    ///  -p0^{n}_{t}   / delta t
    /// </summary>
    public class LowMach_TimeDerivativep0 : IVolumeForm, ISupportsJacobianComponent, IEquationComponentCoefficient, ISpeciesFilter {
        public LowMach_TimeDerivativep0(string spcName, double dt)  {
            ValidSpecies = spcName;
            m_dt = dt;
            //m_ParameterOrdering = new string[] { VariableNames.ThermodynamicPressure, VariableNames.ThermodynamicPressure + "_t0" };
            m_ParameterOrdering = new string[] { "dp0dt"};
            m_ArgumentOrdering = new string[] { };// no arguments
        }

        double m_dt;
        string[] m_ArgumentOrdering;

        string[] m_ParameterOrdering;

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
                return m_ArgumentOrdering;
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
            // For now no update. In future dt could be dynamically changed here
            return;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] {
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            //double p0_Actual = cpv.Parameters[0];
            //double p0_Old = cpv.Parameters[1];
            //double dp0_dt = (p0_Actual - p0_Old) / m_dt;


            double dp0_dt = cpv.Parameters[0];

            return -dp0_dt * V;
        }
        public string ValidSpecies {
            get;
            private set;
        }
    }

 










}
