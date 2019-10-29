/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {

    class ParticleUnderrelaxation {
        /// <summary>
        /// Constructor for the underrelaxation. 
        /// </summary>
        /// <param name="parameter">
        /// <see cref="ParticleUnderrelaxationParam"/>
        /// </param>
        /// <param name="averageForce">
        /// A mean value for the forces used to reduce the impact of computation errors.
        /// </param>
        internal ParticleUnderrelaxation(ParticleUnderrelaxationParam parameter, double averageForce) {
            m_ConvergenceLimit = parameter.m_ConvergenceLimit;
            m_RelaxationFactor = parameter.m_UnderrelaxationFactor;
            m_UseAdaptiveUnderrelaxation = parameter.m_UseAdaptiveUnderrelaxation;
            m_AverageForce = averageForce;
        }
        [DataMember]
        private readonly double m_ConvergenceLimit;
        [DataMember]
        private readonly double m_RelaxationFactor;
        [DataMember]
        private readonly double m_AverageForce;
        [DataMember]
        private readonly bool m_UseAdaptiveUnderrelaxation;

        /// <summary>
        /// This method underrelaxates the hydrodynamic forces. The Underrelaxation
        /// factor is either predefined or is calculated dynamically.
        /// </summary>
        /// <param name="forces">
        /// The hydrodynamic forces.
        /// </param>
        /// <param name="forcesPreviousIteration">
        /// The hydrodynamic forces at the previous iteration.
        /// </param>
        internal void Forces(ref double[] forces, List<double[]> forcesPreviousIteration, ref double[] oldOmega, List<double[]> rawForcesPrev) {
            for (int d = 0; d < forces.Length; d++) {
                List<double> tempForcesPrev = new List<double>();
                List<double> tempRawForcesPrev = new List<double>();
                for (int i = 0; i < forcesPreviousIteration.Count; i++) {
                    tempForcesPrev.Add(forcesPreviousIteration[i][d]);
                    tempRawForcesPrev.Add(rawForcesPrev[i][d]);
                }
                double[] Omega = new double[] { 0, oldOmega[d]};
                AitkenUnderrelaxation(ref forces[d], tempForcesPrev, Omega, rawForcesPrev[1][d]);
                oldOmega[d] = Omega[0];
            }
        }

        internal void Forces(ref double[] forces, List<double[]> forcesPreviousIteration) {
            for (int d = 0; d < forces.Length; d++) {
                List<double> tempForcesPrev = new List<double>();
                for (int i = 0; i < forcesPreviousIteration.Count; i++) {
                    tempForcesPrev.Add(forcesPreviousIteration[i][d]);
                }
                double underrelaxationCoeff = m_UseAdaptiveUnderrelaxation == true
                    ? CalculateAdaptiveUnderrelaxation(forces[d], tempForcesPrev, m_AverageForce, m_ConvergenceLimit, m_RelaxationFactor)
                    : m_RelaxationFactor;
                forces[d] = underrelaxationCoeff * forces[d] + (1 - underrelaxationCoeff) * forcesPreviousIteration[0][d];
            }
        }

        /// <summary>
        /// This method underrelaxates the hydrodynamic torque. The Underrelaxation
        /// factor is either predefined or is calculated dynamically.
        /// </summary>
        /// <param name="torque">
        /// The hydrodynamic torque.
        /// </param>
        /// <param name="torquePreviousIteration">
        /// The hydrodynamic torque at the previous iteration.
        /// </param>
        internal void Torque(ref double torque, List<double> torquePreviousIteration, ref double[] oldOmega, List<double> rawTorquePrev) {
            double[] Omega = new double[] { 0, oldOmega[2] };
            AitkenUnderrelaxation(ref torque, torquePreviousIteration, Omega, rawTorquePrev[1]);
            oldOmega[2] = Omega[0];
        }
        internal void Torque(ref double torque, List<double> torquePreviousIteration) {
            double underrelaxationCoeff = m_UseAdaptiveUnderrelaxation == true
                ? CalculateAdaptiveUnderrelaxation(torque, torquePreviousIteration, m_AverageForce, m_ConvergenceLimit, m_RelaxationFactor)
                : m_RelaxationFactor;
            torque = underrelaxationCoeff * torque + (1 - underrelaxationCoeff) * torquePreviousIteration[0];

        }

        internal void ForcesAndTorque(ref double[] forces, List<double[]> forcesPreviousIteration, ref double oldOmega, List<double[]> rawForcesPrev) {
            double[] Omega = new double[] { 0, oldOmega };
            AitkenUnderrelaxation(ref forces, forcesPreviousIteration, Omega, rawForcesPrev[1]);
            oldOmega = Omega[0];
        }

        /// <summary>
        /// This method calculates the underrelaxation factor dynamically.
        /// </summary>
        /// <param name="variable">
        /// The variable to be relaxated.
        /// </param>
        /// <param name="variableAtPrevIteration">
        /// The variable to be relaxated at the previous iteration.
        /// </param>
        /// <param name="convergenceLimit">
        /// The predefined convergence limit for the fully coupled system.
        /// </param>
        /// <param name="predefinedFactor">
        /// The predefined relaxation factor used a base to calculate a dynamic factor.
        /// </param>
        /// <param name="averageValueOfVar">
        /// In case that there are multiple vars to be underrelaxated, this is the average value.
        /// </param>
        /// <param name="iterationCounter">
        /// No. of iterations.
        /// </param>
        private double CalculateAdaptiveUnderrelaxation(double variable, List<double> variableAtPrevIteration, double averageValueOfVar, double convergenceLimit, double predefinedFactor) {
            double UnderrelaxationCoeff;
            double denominator = variableAtPrevIteration.Count;
            for (int i = 1; i < variableAtPrevIteration.Count; i++) {
                if (Math.Sign(variable - variableAtPrevIteration[0]) != Math.Sign(variableAtPrevIteration[i-1] - variableAtPrevIteration[i])) {
                    predefinedFactor /= denominator;
                    break;
                }
                else
                    denominator -= 1;
            }
            UnderrelaxationCoeff = predefinedFactor * 1 / (Math.Abs((variable - variableAtPrevIteration[0]) / variableAtPrevIteration[0]));
            if (UnderrelaxationCoeff >= 10 * predefinedFactor)
                UnderrelaxationCoeff = 10 * predefinedFactor;

            if (UnderrelaxationCoeff >= 2)
                UnderrelaxationCoeff = 2;

            if (Math.Abs(averageValueOfVar) > 1e4 * Math.Abs(variable) || Math.Abs(variable) < 1e-10) {
                UnderrelaxationCoeff = 1e-20;
            }
            else if (UnderrelaxationCoeff < predefinedFactor * 1e-2)
                UnderrelaxationCoeff = predefinedFactor * 1e-2;
            double GlobalStateBuffer = UnderrelaxationCoeff.MPIMin();
            UnderrelaxationCoeff = GlobalStateBuffer;
            return UnderrelaxationCoeff;
        }

        private void AitkenUnderrelaxation(ref double variable, List<double> variableAtPrevIteration, double[] Omega, double rawForcesPrev) {
            double[] residual = new double[] { (variable - variableAtPrevIteration[0]), (rawForcesPrev - variableAtPrevIteration[1]) };
            Omega[0] = -Omega[1] * residual[1] / (residual[0] - residual[1]);
            variable = Omega[0] * (variable - variableAtPrevIteration[0]) + variableAtPrevIteration[0];
        }

        private void AitkenUnderrelaxation(ref double[] variable, List<double[]> variableAtPrevIteration, double[] Omega, double[] rawForcesPrev) {
            double[][] residual = new double[variable.Length][];
            double[] residualDiff = new double[variable.Length];
            double residualScalar = 0;
            double sumVariable = 0;
            for (int i = 0; i < variable.Length; i++) {
                residual[i] = new double[] { (variable[i] - variableAtPrevIteration[0][i]), (rawForcesPrev[i] - variableAtPrevIteration[1][i]) };
                residualDiff[i] = residual[i][0] - residual[i][1];
                residualScalar += residual[i][1] * residualDiff[i];
                sumVariable += variable[i];
            }
            Omega[0] = -Omega[1] * residualScalar / residualDiff.L2Norm().Pow2();
            Console.WriteLine("OmegaOld: " + Omega[1] + " OmegaNew: " + Omega[0] + " forcesOld 0: " + variable[0] + " forcesOld 1: " + variable[1] + " torqueOld: " + variable[2]);
            for (int i = 0; i < variable.Length; i++) {
                variable[i] = Omega[0] * (variable[i] - variableAtPrevIteration[0][i]) + variableAtPrevIteration[0][i];
            }
            Console.WriteLine("forcesNew 0: " + variable[0] + " forcesNew 1: " + variable[1] + " torqueNew: " + variable[2]);
        }

        private void MinimalPolynomialExtrapolation(ref double variable, List<double> variableAtPrevIteration) {
            double[] polynomialCoefficiants = new double[variableAtPrevIteration.Count];
            polynomialCoefficiants[0] = NewPolynomialCoefficiant(variable, variableAtPrevIteration);
            for (int p = 1; p < polynomialCoefficiants.Length; p++) {
                polynomialCoefficiants[p] = 2 / (double)p;
            }
            double[] normalizedCoefficients = polynomialCoefficiants.CloneAs();
            for (int p = 0; p < normalizedCoefficients.Length; p++) {
                normalizedCoefficients[p] /= polynomialCoefficiants.Sum();
            }
            double test = normalizedCoefficients.Sum();
            variable = normalizedCoefficients[0] * variable;
            for (int k = 1; k < variableAtPrevIteration.Count; k++) {
                variable += normalizedCoefficients[k] * variableAtPrevIteration[k];
            }
        }

        private double NewPolynomialCoefficiant(double variable, List<double> variableAtPrevIteration) {
            double[] residual = new double[variableAtPrevIteration.Count];
            residual[0] = variable - variableAtPrevIteration[1];
            double residualSum = residual[0];
            for (int r = 2; r < variableAtPrevIteration.Count; r++) {
                residual[r] = variableAtPrevIteration[r - 1] - variableAtPrevIteration[r];
                residualSum += Math.Pow(-1, r) * residual[r] / r;
            }
            return -residualSum / residual[0] - 1;
        }
    }
}
