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
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.NSECommon {


    /// <summary>
    /// Auxillary class to compute a source term from a manufactured solution for the scalar transport equations. (I.e. temperature ans species transport equations.)
    /// Current implementation only supports 2D flows.
    /// Current manufactured solutions used is T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
    /// See also ControlManuSol() control function.
    /// </summary>
    public class RHSManuSourceTransportEq : BoSSS.Solution.Utils.LinearSource {

        double HeatReleaseFactor;
        double ReynoldsNumber;
        double PrandtlNumber;
        double SchmidtNumber;
        MaterialLaw Eos;

        String EqType;
        int SpeciesIndex = -1;
        double StoichiometricCoeff;


        double PreExpFactor;
        double ActivationTemperature;
        double MassFraction0Exponent;
        double MassFraction1Exponent;

        double[] MolarMasses;
        double OneOverMolarMass0MolarMass1;
        PhysicsMode physicsMode;
        double phystime;
        bool unsteady;
        SinglePhaseField ThermodynamicPressure;
        /// <summary>
        /// Ctor.
        /// <param name="HeatReleaseFactor">Heat release computed from the sum of the product of the stoichiometric coefficient, partial heat capacity and molar mass of species alpha for all species. I.e.: sum(alpha = 1.. ns)[v_\alpha cp_alpha M_alpha]. Must be computed locally for non-constant partial heat capacities in later iterations of the code.</param>
        /// <param name="Reynolds">The Reynolds number</param>  
        /// <param name="Prandtl">The Prandtl number</param>
        /// <param name="Schmidt">The Schmidt number</param>
        /// <param name="StoichiometricCoefficients">0. of fuel, 1. of oxidizer, 2. of CO2, 3. of H2O</param>  
        /// <param name="ReactionRateConstants">0. PreExpFactor/Damköhler number, 1. ActivationTemperature, 2. MassFraction0Exponent, 3. MassFraction1Exponent</param>  
        /// <param name="MolarMasses">Array of molar masses of fuel, oxidizer, CO2 and H2O</param>  
        /// <param name="OneOverMolarMass0MolarMass1"> 1/(M_infty^(a + b -1) * MolarMassFuel^a * MolarMassOxidizer^b). M_infty is the reference for the molar mass steming from non-dimensionalisation of the governing equations.</param>
        /// <param name="EoS">Material law</param>
        /// <param name="EqType">Temperature" for temperature equation. "MassFraction" for a MassFraction balance</param>
        /// <param name="physicsMode"></param>
        /// <param name="phystime"></param>
        /// <param name="unsteady"></param>
        /// <param name="SpeciesIndex">(optional). Necessary for "MassFraction" EqType. Species index: 0 for fuel, 1 for oxidizer, 2 for CO2 and 3 for H2O.</param>
        /// </summary>
        public RHSManuSourceTransportEq(double HeatReleaseFactor, double Reynolds, double Prandtl, double Schmidt, double[] StoichiometricCoefficients, double[] ReactionRateConstants, double[] MolarMasses, double OneOverMolarMass0MolarMass1, MaterialLaw EoS, String EqType, PhysicsMode physicsMode,double phystime, bool unsteady, SinglePhaseField ThermodynamicPressure, int SpeciesIndex = -1) {
            this.HeatReleaseFactor = HeatReleaseFactor;
            this.ReynoldsNumber = Reynolds;
            this.PrandtlNumber = Prandtl;
            this.SchmidtNumber = Schmidt;
            this.Eos = EoS;
            this.EqType = EqType;
            this.physicsMode = physicsMode;
            this.phystime = phystime;
            this.unsteady = unsteady;

            if (EqType == "MassFraction") {
                this.SpeciesIndex = SpeciesIndex;
                this.StoichiometricCoeff = StoichiometricCoefficients[SpeciesIndex];
            }

            this.PreExpFactor = ReactionRateConstants[0];
            this.ActivationTemperature = ReactionRateConstants[1];
            this.MassFraction0Exponent = ReactionRateConstants[2];
            this.MassFraction1Exponent = ReactionRateConstants[3];

            this.MolarMasses = MolarMasses;
            this.OneOverMolarMass0MolarMass1 = OneOverMolarMass0MolarMass1;
            this.ThermodynamicPressure = ThermodynamicPressure;

        }


        /// <summary>
        /// None
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[0]; }
        }

        /// <summary>
        /// Temperature, MassFraction0, MassFraction1, MassFraction2, MassFraction3 at linearization point.
        /// </summary>
        public override IList<string> ParameterOrdering {
            get { return new string[] { /*VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0*/ }; }
        }


        //Manufactured solution for T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
        protected override double Source(double[] x, double[] parameters, double[] U) {

            double x_ = x[0];
            double y_ = x[1];
            double t_ = 0.0;
            double p0 = ThermodynamicPressure.GetMeanValue(3);

            double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3];
            double alpha1 = 0.3;
            double alpha2 = 0.6;
            double alpha3 = 0.1;
            double[] Coefficients = new double[] { alpha1, alpha2, alpha3 };

            bool unsteady = false;
            double ConvectionTerm;
            double ReactionRate;
            double SourceTerm;
            double DiffussionTerm;
            double unsteadyTerm = 0.0;


            if (unsteady) {

                switch (physicsMode) {
                    case PhysicsMode.LowMach:
                        ConvectionTerm = p0 * t_ * Math.Sin(x_ * t_) + p0 * t_ * Math.Sin(y_ * t_);
                        ReactionRate = 0;
                        unsteadyTerm = 0.0; // 0.0 is the correct MS... ((p0/T)*T)' = (p0)' = 0.0
                        break;
                    case PhysicsMode.Combustion:
                        ConvectionTerm = 0.0;//TODO

                        ReactionRate = 0.0;//TODO

                        break;
                    default:
                        throw new NotImplementedException("wrong switch");
                }
                switch (EqType) {
                    case "Temperature":
                        DiffussionTerm = 1.0 / (ReynoldsNumber * PrandtlNumber) * (y_ * y_ * t_ * t_ * Math.Cos(x_ * y_ * t_) + x_ * x_ * t_ * t_ * Math.Cos(x_ * y_ * t_));

                        SourceTerm = 0.0;//TODO

                        return -(unsteadyTerm + ConvectionTerm + DiffussionTerm + SourceTerm);

                    case "MassFraction":
                        if (SpeciesIndex == -1)
                            throw new ArgumentException("Species index needs to be specified");
                        ConvectionTerm *= Coefficients[SpeciesIndex];

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
                            default:
                                throw new NotImplementedException("wrong index");

                        }
                        switch (physicsMode) {
                            case PhysicsMode.LowMach:
                                DiffussionTerm = 0.0;//TODO
                                break;
                            case PhysicsMode.Combustion:
                                DiffussionTerm = 0.0;//TODO
                                break;
                            default:
                                throw new NotImplementedException("wrong switch");
                        }


                        SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                        return -(ConvectionTerm + -DiffussionTerm + SourceTerm); // TODO CHECK SIGNS

                    default:
                        throw new NotImplementedException();
                }

            }
            else {
                switch (physicsMode) {
                    case PhysicsMode.LowMach:
                        ConvectionTerm = p0 * Math.Sin(x_) + p0 * Math.Sin(y_);
                        ReactionRate = 0;
                        break;
                    case PhysicsMode.Combustion:
                        ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_)) / M4) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_)) / M4) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(y_);
                        ReactionRate = PreExpFactor * OneOverMolarMass0MolarMass1 * Math.Exp(-0.1e1 / Math.Cos(x_ * y_)) * Math.Pow(p0, 0.3e1) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.3e1) * alpha1 * alpha2 * alpha2;
                        break;
                    default:
                        throw new NotImplementedException("wrong switch");
                }
                switch (EqType) {
                    case "Temperature":
                        DiffussionTerm = 1.0 / (ReynoldsNumber * PrandtlNumber) * Math.Cos(x_ * y_) * (Math.Pow(x_, 2) + Math.Pow(y_, 2)); // OK
                        SourceTerm = HeatReleaseFactor * Math.Cos(x_ * y_) * ReactionRate;
                        return -(ConvectionTerm + DiffussionTerm + SourceTerm);

                    case "MassFraction":
                        if (SpeciesIndex == -1)
                            throw new ArgumentException("Species index needs to be specified");
                        ConvectionTerm *= Coefficients[SpeciesIndex];

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
                            default:
                                throw new NotImplementedException("wrong index");

                        }
                        switch (physicsMode) {
                            case PhysicsMode.LowMach:
                                DiffussionTerm = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * alpha * y_ * y_ * Math.Pow(Math.Sin(x_ * y_), 0.2e1) - p0 * alpha * y_ * y_ - p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * alpha * x_ * x_ * Math.Pow(Math.Sin(x_ * y_), 0.2e1) - p0 * alpha * x_ * x_;
                                break;
                            case PhysicsMode.Combustion:
                                DiffussionTerm = 1.0 / (ReynoldsNumber * SchmidtNumber) * (-p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * alpha * y_ * y_ * Math.Pow(Math.Sin(x_ * y_), 0.2e1) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * alpha * y_ * Math.Sin(x_ * y_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_)) / M4) - p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * alpha * y_ * y_ - p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * alpha * x_ * x_ * Math.Pow(Math.Sin(x_ * y_), 0.2e1) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * alpha * x_ * Math.Sin(x_ * y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_)) / M4) - p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * alpha * x_ * x_);
                                break;
                            default:
                                throw new NotImplementedException("wrong switch");
                        }


                        SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                        return -(ConvectionTerm + -DiffussionTerm + SourceTerm); // TODO CHECK SIGNS

                    default:
                        throw new NotImplementedException();
                }

            }
        }
    }
}
