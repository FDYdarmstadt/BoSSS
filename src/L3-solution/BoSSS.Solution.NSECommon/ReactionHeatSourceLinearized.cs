/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
//using System.Diagnostics.Eventing.Reader;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP.Utils;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Reaction heat source in temperature equation.
    /// </summary>
    public class ReactionHeatSourceLinearized : BoSSS.Solution.Utils.LinearSource, IEquationComponentCoefficient {
        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;
        double ReactionRate;
        double HeatReleaseFactor;
        double[] ReactionRateConstants;
        double OneOverMolarMass0MolarMass1;

        MaterialLaw EoS;
        double rho;
        double m_Da;
        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="HeatReleaseFactor">Heat release computed from the sum of the product of the stoichiometric coefficient, partial heat capacity and molar mass of species alpha for all species. I.e.: sum(alpha = 1.. ns)[v_\alpha cp_alpha M_alpha]. Must be computed locally for non-constant partial heat capacities in later iterations of the code.</param>   
        /// <param name="ReactionRateConstants">0. PreExpFactor/Damköhler number, 1. ActivationTemperature, 2. MassFraction0Exponent, 3. MassFraction1Exponent</param>  
        /// <param name="OneOverMolarMass0MolarMass1"> 1/(M_infty^(a + b -1) * MolarMassFuel^a * MolarMassOxidizer^b). M_infty is the reference for the molar mass steming from non-dimensionalisation of the governing equations.</param>  
        /// <param name="EoS">MaterialLawCombustion</param>  
        public ReactionHeatSourceLinearized(double HeatReleaseFactor, double[] ReactionRateConstants, double OneOverMolarMass0MolarMass1, MaterialLaw EoS) {
            m_ArgumentOrdering = new string[] { VariableNames.Temperature };
            m_ParameterOrdering = new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 };
            this.HeatReleaseFactor = HeatReleaseFactor;
            this.ReactionRateConstants = ReactionRateConstants;
            this.OneOverMolarMass0MolarMass1 = OneOverMolarMass0MolarMass1;
            this.EoS = EoS;
            m_Da = ReactionRateConstants[0]; // Damköhler number 


        }


        /// <summary>
        /// Temperature
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }

        /// <summary>
        /// Temperature, MassFraction0, MassFraction1, MassFraction2, ..., MassFraction_ns at linearization point.
        /// </summary>
        public override IList<string> ParameterOrdering {
            get { return m_ParameterOrdering; }
        }


        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Damkoehler"))
                m_Da = (double)cs.UserDefinedValues["Damkoehler"];
        }

        /// <summary>
        /// 
        /// </summary>
        protected override double Source(double[] x, double[] parameters, double[] U) {
            rho = EoS.GetDensity(parameters);
            Debug.Assert(!double.IsNaN(rho));
            Debug.Assert(!double.IsInfinity(rho));

            ReactionRate = m_Da * Math.Exp(-ReactionRateConstants[1] / parameters[0]) * OneOverMolarMass0MolarMass1 * Math.Pow(rho * parameters[1], ReactionRateConstants[2]) * Math.Pow(rho * parameters[2], ReactionRateConstants[3]);

            return HeatReleaseFactor * U[0] * ReactionRate;

        }
    }



    /// <summary>
    /// Reaction heat source in temperature equation.
    /// </summary>
    public class ReactionHeatSourceJacobi : IVolumeForm, IEquationComponentCoefficient, ISupportsJacobianComponent {
        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;
        double ReactionRate;
        double HeatReleaseFactor;
        double[] ReactionRateConstants;
        double[] molarMasses;

        MaterialLaw_MultipleSpecies EoS;
        double rho;
        double m_Da;
        double TRef;
        double cpRef;
        bool VariableOneStepParameters;
        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="HeatReleaseFactor">Heat release computed from the sum of the product of the stoichiometric coefficient, partial heat capacity and molar mass of species alpha for all species. I.e.: sum(alpha = 1.. ns)[v_\alpha cp_alpha M_alpha]. Must be computed locally for non-constant partial heat capacities in later iterations of the code.</param>   
        /// <param name="ReactionRateConstants">0. PreExpFactor/Damköhler number, 1. ActivationTemperature, 2. MassFraction0Exponent, 3. MassFraction1Exponent</param>  
        /// <param name="EoS">MaterialLawCombustion</param>  
        public ReactionHeatSourceJacobi(double HeatReleaseFactor, double[] ReactionRateConstants, double[] molarmasses, MaterialLaw_MultipleSpecies EoS, double TRef, double cpRef, bool VariableOneStepParameters, int NumberOfReactants) {
            m_ArgumentOrdering = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(NumberOfReactants));// Y4 is not a variable!!!!;
            m_ParameterOrdering = new string[] { "kReact" };
            this.HeatReleaseFactor = HeatReleaseFactor;
            this.ReactionRateConstants = ReactionRateConstants;
            this.molarMasses = molarmasses;
            this.EoS = EoS;
            m_Da = ReactionRateConstants[0]; // Damköhler number 
            this.TRef = TRef;
            this.cpRef = cpRef;
            this.VariableOneStepParameters = VariableOneStepParameters;
        }


        /// <summary>
        /// Temperature
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }

        /// <summary>
        /// Temperature, MassFraction0, MassFraction1, MassFraction2, ..., MassFraction_ns at linearization point.
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get { return m_ParameterOrdering; }
        }

        public virtual TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Damkoehler"))
                m_Da = (double)cs.UserDefinedValues["Damkoehler"];
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var SourceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { SourceDerivVol };
        }


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this.Source(cpv.Xglobal, cpv.Parameters, U) * V;
        }

        /// <summary>
        /// 
        /// </summary>
        protected double Source(double[] x, double[] parameters, double[] U) {
       
            Debug.Assert(!double.IsNaN(rho));
            Debug.Assert(!double.IsInfinity(rho));

            double Ta = ReactionRateConstants[1];
 
            double Temperature = U[0];
            double YF = U[1];
            double YO = U[2];


            ////===================================================
            //// Limit value of variables using known bounds
            ////====================================================

            Temperature = Temperature > 10 ? 10 : Temperature;
            Temperature = Temperature < 0.7 ? 0.7 : Temperature;

            YF = YF > 1.0 ? 1.0 : YF;
            YF = YF < 0.0 ? 0.0 : YF;

            YO = YO > 1.0 ? 1.0 : YO;
            YO = YO < 0.0 ? 0.0 : YO;


            double cp = 1.0;// EoS.GetMixtureHeatCapacity(U);

            if (YF * YO > 1e-8 && VariableOneStepParameters) {//  calculate one-Step model parameters
                Ta = EoS.m_ChemModel.getTa(YF, YO) / TRef;
                HeatReleaseFactor = EoS.m_ChemModel.getHeatRelease(YF, YO) / (cpRef * TRef);
            }


            rho = EoS.GetDensity(U);
            double PM_CH4 = molarMasses[0];
            double PM_O2 = molarMasses[1];
            ReactionRate = m_Da * Math.Exp(-Ta / Temperature) * (rho * YF / PM_CH4) * (rho * YO / PM_O2);

            if (double.IsInfinity(ReactionRate)) {
                Console.WriteLine("Infinite found");
                Console.WriteLine("Temperature: {0}", Temperature);
                Console.WriteLine("rho:{0}", rho);
                Console.WriteLine("ExponentialTerm:{0}", Math.Exp(-Ta / Temperature));

            }

            if (double.IsNaN(ReactionRate)) {
                Console.WriteLine("Nan found");
                Console.WriteLine("Temperature:{0}", Temperature);
                Console.WriteLine("rho:{0}", rho);
                Console.WriteLine("ExponentialTerm:{0}", Math.Exp(-Ta / Temperature));
            }


        

            return -HeatReleaseFactor * ReactionRate * PM_CH4 / cp;

        }
    }
}
