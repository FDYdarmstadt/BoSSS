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
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.NSECommon {


    /// <summary>
    /// Auxillary class to compute a source term from a manufactured solution for the scalar transport equations. (I.e. temperature ans species transport equations.)
    /// Current implementation only supports 2D flows.
    /// Current manufactured solutions used is T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
    /// Constant physical properties are asumed ( viscosity (mu), thermal conductivity(lambda),mass diffusivity(D), heat capacity(Cp)) 
    /// See also ControlManuSol() control function.
    /// </summary>
    //public class RHSManuSourceTransportEq : BoSSS.Solution.Utils.LinearSource {
    public class RHSManuSourceTransportEq : IVolumeForm, ISupportsJacobianComponent ,IEquationComponentCoefficient  {

        double HeatRelease;
        double ReynoldsNumber;
        double PrandtlNumber;
        double SchmidtNumber;
        MaterialLaw Eos;

        String EqType;
        int SpeciesIndex = -1;
        double StoichiometricCoeff;


        double Da;
        double Ta;
        double a;
        double b;
        double time;
        double[] MolarMasses;
        double OneOverMolarMass0MolarMass1;
        PhysicsMode physicsMode;
        bool chemReactionOK;
        bool rhoOne;
        Func<double[], double, double> sourceFunc = null;
        /// <summary>
        /// Da number used within the homotopie algorithm
        /// </summary>
        /// <param name="cs"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Damkoehler"))
                Da = (double)cs.UserDefinedValues["Damkoehler"];
            if (cs.UserDefinedValues.Keys.Contains("time"))
                time = (double)cs.UserDefinedValues["time"];
        }

        /// <summary>
        /// Ctor.
        /// <param name="HeatRelease">Heat release computed from the sum of the product of the stoichiometric coefficient, partial heat capacity and molar mass of species alpha for all species. I.e.: sum(alpha = 1.. ns)[v_\alpha cp_alpha M_alpha]. Must be computed locally for non-constant partial heat capacities in later iterations of the code.</param>
        /// <param name="Reynolds">The Reynolds number</param>  
        /// <param name="Prandtl">The Prandtl number</param>
        /// <param name="Schmidt">The Schmidt number</param>
        /// <param name="StoichiometricCoefficients">0. of fuel, 1. of oxidizer, 2. of CO2, 3. of H2O</param>  
        /// <param name="ReactionRateConstants">0. PreExpFactor/Damköhler number, 1. ActivationTemperature, 2. MassFraction0Exponent, 3. MassFraction1Exponent</param>  
        /// <param name="MolarMasses">Array of molar masses of fuel, oxidizer, CO2 and H2O</param>   
        /// <param name="EoS">Material law</param>
        /// <param name="EqType">Temperature" for temperature equation. "MassFraction" for a MassFraction balance</param>
        /// <param name="physicsMode"></param>
        /// <param name="phystime"></param>
        /// <param name="unsteady"></param>
        /// <param name="SpeciesIndex">(optional). Necessary for "MassFraction" EqType. Species index: 0 for fuel, 1 for oxidizer, 2 for CO2 and 3 for H2O.</param>
        /// <param name="chemReactionOK"></param>
        /// <param name="rhoOne"></param>
        /// </summary>
        public RHSManuSourceTransportEq(double HeatRelease, double Reynolds, double Prandtl, double Schmidt, double[] StoichiometricCoefficients, double[] ReactionRateConstants, double[] MolarMasses, MaterialLaw EoS, String EqType, PhysicsMode physicsMode, int SpeciesIndex = -1, bool chemReactionOK = true, bool rhoOne = false, Func<double[], double, double> _sourceFunc = null) {
            this.HeatRelease = HeatRelease;
            this.ReynoldsNumber = Reynolds;
            this.PrandtlNumber = Prandtl;
            this.SchmidtNumber = Schmidt;
            this.Eos = EoS;
            this.EqType = EqType;
            this.physicsMode = physicsMode;
            this.chemReactionOK = chemReactionOK;
            this.rhoOne = rhoOne;
            if(EqType == "MassFraction") {
                this.SpeciesIndex = SpeciesIndex;
                this.StoichiometricCoeff = StoichiometricCoefficients[SpeciesIndex];
            }

            this.Da = ReactionRateConstants[0]; // Damköhler number 
            this.Ta = ReactionRateConstants[1];
            this.a =  ReactionRateConstants[2];
            this.b =  ReactionRateConstants[3];

            this.MolarMasses = MolarMasses;

            this.sourceFunc = _sourceFunc;
        }


        /// <summary>
        /// None
        /// </summary>
        public IList<string> ArgumentOrdering {
            get { return new string[0]; }
        }

        /// <summary>
        /// None
        /// </summary>
        public IList<string> ParameterOrdering {
            get { return null; }
        }
        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }


        /// <summary>
        /// None
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.AllOn;
            }
        }


        public double getAlpha(int SpeciesIndex) {
            double alpha1 = 0.3;
            double alpha2 = 0.4;
            double alpha3 = 0.1;
            double alpha4 = 0.2;
            double alpha;
            switch (SpeciesIndex) {
                case 0:
                    alpha = alpha1;
                    break;
                case 1:
                    alpha = alpha2;
                    break;
                case 2:
                    alpha = alpha3;
                    break;
                case 3:
                    alpha = alpha4;
                    break;
                default:
                    throw new NotImplementedException("wrong index");
            }
            return alpha;

        }
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            ////Manufactured solution for T = cos(x*y), Y0 = 0.2 cos(x*y), Y1 = 0.3 cos(x*y), Y2 = 0.1 cos(x*y), Y3 = 0.1 cos(x*y) => Y4 = 0.3 cos(x*y)
            /// u = -cos(x), v = -cos(y), p = sin(x*y).
            double[] x = cpv.Xglobal;
            double x_ = x[0];
            double y_ = x[1];
            double p0 = 1.0;

            double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
            double alpha1 = 0.3;
            double alpha2 = 0.4;
            double alpha3 = 0.1;
            double alpha4 = 0.2;
            double[] Coefficients = new double[] { alpha1, alpha2, alpha3, alpha4 };


            double ConvectionTerm = 0;
            double ReactionRate = 0;
            double SourceTerm = 0;
            double DiffussionTerm = 0;

            double ConvectionTermSwitch = 1.0;
            double DiffussionTermSwitch = 1.0;
            double res = 0.0;
            double alpha;
            switch (EqType) {
                case "Temperature":
                    switch (physicsMode) {
                        case PhysicsMode.LowMach:
                            ConvectionTerm = p0 * Math.Sin(x_) + p0 * Math.Sin(y_); // This is the MS if the discretized convective term of the energy equation is div(rho*u*T)
                            ReactionRate = 0;
                            DiffussionTerm = 1.0 / (ReynoldsNumber * PrandtlNumber) * Math.Cos(x_ * y_) * (Math.Pow(x_, 2) + Math.Pow(y_, 2)); // OK for lowmach AND combustion
                            SourceTerm = -HeatRelease * ReactionRate * M1;
                            res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                            break;
                        case PhysicsMode.Combustion:
                            ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                            ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);
                            if (rhoOne) {
                                ConvectionTerm = 0.10e1 * Math.Sin(x_) * Math.Cos(x_ * y_) + 0.10e1 * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + 0.10e1 * Math.Sin(y_) * Math.Cos(x_ * y_) + 0.10e1 * Math.Cos(y_) * x_ * Math.Sin(x_ * y_);
                                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(0.10e1 * alpha1 * Math.Cos(x_ * y_), a) * Math.Pow(0.10e1 * alpha2 * Math.Cos(x_ * y_), b);
                            }

                            if (!chemReactionOK)
                                ReactionRate = 0.0;
                            DiffussionTerm = 1.0 / (ReynoldsNumber * PrandtlNumber) * Math.Cos(x_ * y_) * (Math.Pow(x_, 2) + Math.Pow(y_, 2)); // OK for lowmach AND combustion
                            SourceTerm = -HeatRelease * ReactionRate /** M1*/;
                            res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                            break;
                        case PhysicsMode.MixtureFraction:
                            res = sourceFunc(x, time);
                            break;
                        default:
                            throw new NotImplementedException("wrong switch");
                    }
                    break;
                case "MassFraction":
                    switch (physicsMode) {
                        case PhysicsMode.LowMach:
                            ConvectionTerm = p0 * Math.Sin(x_) + p0 * Math.Sin(y_); // This is the MS if the discretized convective term of the energy equation is div(rho*u*T)
                            ReactionRate = 0;
                            if (SpeciesIndex == -1)
                                throw new ArgumentException("Species index needs to be specified");
                            ConvectionTerm *= Coefficients[SpeciesIndex];
                            alpha = getAlpha(SpeciesIndex);
                               
                            // The diffusive term corresponds to div(Diffusivity*rho*gradY0)! not to div(grad(rho*Y0))
                            DiffussionTerm = 1.0 / (ReynoldsNumber * SchmidtNumber) * alpha * y_ * y_ * Math.Cos(x_ * y_) + alpha * x_ * x_ * Math.Cos(x_ * y_);

                            SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                            res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                            break;
                        case PhysicsMode.Combustion:
                            ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                            ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);
                            if (rhoOne) {
                                ConvectionTerm = 0.10e1 * Math.Sin(x_) * Math.Cos(x_ * y_) + 0.10e1 * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + 0.10e1 * Math.Sin(y_) * Math.Cos(x_ * y_) + 0.10e1 * Math.Cos(y_) * x_ * Math.Sin(x_ * y_);
                                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(0.10e1 * alpha1 * Math.Cos(x_ * y_), a) * Math.Pow(0.10e1 * alpha2 * Math.Cos(x_ * y_), b);
                            }

                            if (!chemReactionOK)
                                ReactionRate = 0.0;
                            if (SpeciesIndex == -1)
                                throw new ArgumentException("Species index needs to be specified");
                            ConvectionTerm *= Coefficients[SpeciesIndex];
                            alpha = getAlpha(SpeciesIndex);
                            // The diffusive term corresponds to div(Diffusivity*rho*gradY0)! not to div(grad(rho*Y0))
                            DiffussionTerm = 1.0 / (ReynoldsNumber * SchmidtNumber) * alpha * y_ * y_ * Math.Cos(x_ * y_) + alpha * x_ * x_ * Math.Cos(x_ * y_);

                            SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                            res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                            break;

                        case PhysicsMode.MixtureFraction:                 
                            res = sourceFunc(x, time);
                            break;

                        default:
                            throw new NotImplementedException("wrong switch");
                    }

                    break;
                case "MixtureFraction":
                    res = sourceFunc(x, time);
                    break;
            }


            Debug.Assert(!double.IsNaN(SourceTerm) && !double.IsInfinity(SourceTerm));

           
            return -(res) * V;
        }
    }
}
 
 