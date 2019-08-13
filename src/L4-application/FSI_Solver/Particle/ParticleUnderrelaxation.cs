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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    class ParticleUnderrelaxation
    {
        /// <summary>
        /// This method underrelaxates the hydrodynamic forces and torque. The Underrelaxation
        /// factor is either predefined or is calculated dynamically.
        /// </summary>
        /// <param name="forces">
        /// The hydrodynamic forces.
        /// </param>
        /// <param name="torque">
        /// The hydrodynamic torque.
        /// </param>
        /// <param name="forcesAtPrevIteration">
        /// The hydrodynamic forces at the previous iteration.
        /// </param>
        /// <param name="torqueAtPrevIteration">
        /// The hydrodynamic torque at the previous iteration.
        /// </param>
        /// <param name="convergenceLimit">
        /// The predefined convergence limit for the fully coupled system.
        /// </param>
        /// <param name="relaxationFactor">
        /// The predefined relaxation factor, either used directly or as a base to calculate
        /// a dynamic factor.
        /// </param>
        /// <param name="clearSmallValues">
        /// Bool to check whether small components should be set to zero (tool to stabilize the
        /// simulation).
        /// </param>
        /// <param name="UseAdaptiveUnderrelaxation">
        /// Bool to check whether the relaxation factor should be calculated dynamically or not.
        /// </param>
        /// <param name="averageDistance">
        /// The average Lengthscale of the particle.
        /// </param>
        /// <param name="iterationCounter">
        /// No. of iterations.
        /// </param>
        internal void Forces(ref double[] forces, double[] forcesAtPrevIteration, double convergenceLimit, double relaxationFactor, bool UseAdaptiveUnderrelaxation, double averageForce)
        {
            int spatialDim = forces.Length;
            double[] underrelaxationCoeff = new double[spatialDim];

            for (int d = 0; d < spatialDim; d++)
            {
                underrelaxationCoeff[d] = UseAdaptiveUnderrelaxation == true
                    ? CalculateAdaptiveUnderrelaxation(forces[d], forcesAtPrevIteration[d], averageForce, convergenceLimit, relaxationFactor)
                    : relaxationFactor;
            }
            for (int d = 0; d < spatialDim; d++)
            {
                forces[d] = underrelaxationCoeff[d] * forces[d] + (1 - underrelaxationCoeff[d]) * forcesAtPrevIteration[d];
            }
        }

        internal void Torque(ref double torque, double torqueAtPrevIteration, double convergenceLimit, double relaxationFactor, bool UseAdaptiveUnderrelaxation, double averageForce)
        {
            double underrelaxationCoeff = UseAdaptiveUnderrelaxation == true
                ? CalculateAdaptiveUnderrelaxation(torque, torqueAtPrevIteration, averageForce, convergenceLimit, relaxationFactor)
                : relaxationFactor;
            torque = underrelaxationCoeff * torque + (1 - underrelaxationCoeff) * torqueAtPrevIteration;

        }

        /// <summary>
        /// Does what it says.
        /// </summary>
        /// <param name="forces">
        /// The hydrodynamic forces.
        /// </param>
        /// <param name="torque">
        /// The hydrodynamic torque.
        /// </param>
        /// <param name="averageDistance">
        /// The average Lengthscale of the particle.
        /// </param>
        public void CalculateAverageForces(double[] forces, double torque, double averageDistance, out double averageForces)
        {
            averageForces = Math.Abs(torque) / averageDistance;
            for (int d = 0; d < forces.Length; d++)
            {
                averageForces += forces[d];
            }
            averageForces /= 3;
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
        private double CalculateAdaptiveUnderrelaxation(double variable, double variableAtPrevIteration, double averageValueOfVar, double convergenceLimit, double predefinedFactor)
        {
            double ConvergenceHelperFactor = 1;
            double UnderrelaxationCoeff = predefinedFactor * 1e-1;
            double UnderrelaxationExponent = 0;

            while (Math.Abs(UnderrelaxationCoeff * variable) > 0.75 * Math.Abs(variableAtPrevIteration) && UnderrelaxationCoeff > 1e-20)
            {
                UnderrelaxationExponent -= 1;
                UnderrelaxationCoeff = predefinedFactor * Math.Pow(10, UnderrelaxationExponent);
            }

            if (Math.Abs(UnderrelaxationCoeff * variable) < convergenceLimit * 100 && 1000 * Math.Abs(variable) > Math.Abs(averageValueOfVar))
                UnderrelaxationCoeff = predefinedFactor * convergenceLimit * 10;

            if (UnderrelaxationCoeff >= predefinedFactor * 1e-1)
                UnderrelaxationCoeff = predefinedFactor * 1e-1;

            double GlobalStateBuffer = UnderrelaxationCoeff.MPIMin();
            UnderrelaxationCoeff = GlobalStateBuffer;

            return UnderrelaxationCoeff * ConvergenceHelperFactor;
        }
    }
}
