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

using MPI.Wrappers;
using System;
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
            m_ConvergenceLimit = parameter.ConvergenceLimit;
            m_RelaxationFactor = parameter.UnderrelaxationFactor;
            m_UseAdaptiveUnderrelaxation = parameter.UsaAddaptiveUnderrelaxation;
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
        internal void Forces(ref double[] forces, double[] forcesPreviousIteration) {
            for (int d = 0; d < forces.Length; d++) {
                double underrelaxationCoeff = m_UseAdaptiveUnderrelaxation == true
                    ? CalculateAdaptiveUnderrelaxation(forces[d], forcesPreviousIteration[d], m_AverageForce, m_ConvergenceLimit, m_RelaxationFactor)
                    : m_RelaxationFactor;
                forces[d] = underrelaxationCoeff * forces[d] + (1 - underrelaxationCoeff) * forcesPreviousIteration[d];
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
        internal void Torque(ref double torque, double torquePreviousIteration) {
            double underrelaxationCoeff = m_UseAdaptiveUnderrelaxation == true
                ? CalculateAdaptiveUnderrelaxation(torque, torquePreviousIteration, m_AverageForce, m_ConvergenceLimit, m_RelaxationFactor)
                : m_RelaxationFactor;
            torque = underrelaxationCoeff * torque + (1 - underrelaxationCoeff) * torquePreviousIteration;

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
        private double CalculateAdaptiveUnderrelaxation(double variable, double variableAtPrevIteration, double averageValueOfVar, double convergenceLimit, double predefinedFactor) {
            double ConvergenceHelperFactor = 1;
            double UnderrelaxationCoeff = predefinedFactor * 1e-1;
            double UnderrelaxationExponent = 0;

            while (Math.Abs(UnderrelaxationCoeff * variable) > 0.75 * Math.Abs(variableAtPrevIteration) && UnderrelaxationCoeff > 1e-20) {
                UnderrelaxationExponent -= 1;
                UnderrelaxationCoeff = predefinedFactor * Math.Pow(10, UnderrelaxationExponent);
            }

            if (Math.Abs(UnderrelaxationCoeff) < convergenceLimit * 1 && 1000 * Math.Abs(variable) > Math.Abs(averageValueOfVar))
                UnderrelaxationCoeff = predefinedFactor * convergenceLimit * 1;

            if (UnderrelaxationCoeff >= predefinedFactor * 1e-1)
                UnderrelaxationCoeff = predefinedFactor * 1e-1;

            double GlobalStateBuffer = UnderrelaxationCoeff.MPIMin();
            UnderrelaxationCoeff = GlobalStateBuffer;

            return UnderrelaxationCoeff * ConvergenceHelperFactor;
        }
    }
}
