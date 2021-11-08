using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.NSECommon {
    //[DataContract]
    //[Serializable]


    public class ReactionParameters : ICloneable {





        /// <summary>
        /// For now only considering a one-step reaction
        /// </summary>
        /// <param name="ChemicalComponents"></param>
        /// <param name=""></param>
        public ReactionParameters(string[] ChemicalComponents, double[] StoichiometricCoefficients ) {
            if (ChemicalComponents.Length != StoichiometricCoefficients.Length)
                throw new Exception("lengths of the component arrays and stoichiomnetric coefficients dont coincide");

            m_ChemicalComponents = ChemicalComponents;
            m_MolarMasses = MolarMassVector(ChemicalComponents);
            m_StoichiometricCoefficients = StoichiometricCoefficients;



            // which reaction model?
            //reactionModel = OneStepModel;
        }

        /// <summary>
        /// 
        /// </summary>
        Dictionary<string, double> MolarMassDictionary = new Dictionary<string, double>  {
            { "CH4", 16 },
            { "O2", 32},
            { "CO2", 44 },
            { "H2O", 18},
            { "N2", 28 }
        };


        private double[] MolarMassVector(string[] ChemicalComponents) {
            double[] MMVec = new double[ChemicalComponents.Length];
            for (int i = 0; i < ChemicalComponents.Length; i++)
                MMVec[i] = MolarMassDictionary[ChemicalComponents[i]];
            return MMVec;
        }


        string[] m_ChemicalComponents;
        public int totalNumChemComponents {
            get {
                if (m_ChemicalComponents == null || m_ChemicalComponents.Length == 0)
                    throw new Exception("Chemical composition is wrong");

                return m_ChemicalComponents.Length;
            }
        }

        /// <summary>
        /// Array of length <see cref="totalNumChemComponents"/> containing the Stoichiometric coefficients
        /// used for each chemical reaction
        /// </summary>
        double[] m_StoichiometricCoefficients;


        /// <summary>
        /// Array of length <see cref="totalNumChemComponents"/> containing all molar masses   
        /// </summary>
        double[] m_MolarMasses;


        // <summary>
        // One-step kinetic model constants for an arrhenius expression A*exp(-Ta/T)*c_f^a*c_o^b
        // 1) A
        // 2) Ta
        // 3) a
        // 4)b
        // </summary>



        ///// <summary>
        ///// Base class for implementing different material laws.
        ///// </summary>
        //public abstract class ReactionModel {


        //}
        //class MethaneCombustion_OneStepModel : ReactionModel {




        //    /// <summary>
        //    /// Arrhenius preExponential factor of the one-Step combustion of hydrocarbon.
        //    /// </summary>
        //    [DataMember]
        //    public double PreExponentialFactor = 6.9e11; // Pre-exponential factor,  6.9e14 cm3/(mol s) =>  6.9e11m3/(kmol s)
        //    /// <summary>
        //    /// Activation temperature of the one-Step combustion of hydrocarbon.
        //    /// </summary>
        //    [DataMember]
        //    public double Ta = 15900; // Activation TEMPERATURE, K


        //    /// <summary>
        //    /// Heat release of combustion for a one-step chemistry with (phi < 1), per unit mass
        //    /// </summary>
        //    [DataMember]
        //    public double HeatReleaseMass = 50100; //  KJ/(Kg fuel)

        //    /// <summary>
        //    /// Heat release of combustion for a one-step chemistry with (phi < 1), per molar unit
        //    /// </summary>
        //    [DataMember]
        //    public double HeatReleaseMolar = 802400; //  KJ/(kmol fuel)


        //    double a { get; set; } = 1;
        //    double b { get; set; } = 1;
        //    /// <summary>
        //    /// Reaction of methane combustion:
        //    /// CH4 + 2O2 -> CO2 + 2H2O 
        //    /// the fifth component corresponds to an inert (Nitrogen)
        //    /// </summary>
        //    double[] stoichiometricArray { get; set; } = new double[] { -1, -2, 1, 2, 0 };



        //}

        /// <summary>
        /// Arrhenius chemical reaction parameters
        /// </summary>
        interface ChemicalReaction {
            double ActivationTemperature { get; set; }
            double PreExponentialFactor { get; set; }
            double a { get; set; }
            double b { get; set; }
        }

        public object Clone() {
            throw new NotImplementedException();
        }
    }

    public class ChemicalReaction {
        /// <summary>
        /// base constructor
        /// </summary>
        public ChemicalReaction() {
            ArrheniusParams = new double[] { PreExponentialFactor, ActivationTemperature, a, b };
        }
        /// <summary>     
        /// </summary>
        public double[] m_StoichiometricCoefficients { get; set; }

        double[] ArrheniusParams { get; set; }

        public double ActivationTemperature { get; set; }
        public double PreExponentialFactor { get; set; }
        public double a { get; set; }
        public double b { get; set; }

        public string[] ComponentNames { get; set; }


    }

    public class OneStepMethaneCombustion : ChemicalReaction {
        public OneStepMethaneCombustion() :base()  {
            base.ActivationTemperature = 15900;
            base.PreExponentialFactor = 100000;
            base.a = 1;
            base.b = 2;
            base.m_StoichiometricCoefficients = new double[] { -1, -2, 1, 2, 0 }; ;
        }


    }







}




