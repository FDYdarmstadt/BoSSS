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
        /// <param name="RelaxationFactor">
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
        internal double[] RelaxatedForcesAndTorque(double[] forces, double torque, double[] forcesAtPrevIteration, double torqueAtPrevIteration, double convergenceLimit, double RelaxationFactor, bool clearSmallValues, bool UseAdaptiveUnderrelaxation, double averageDistance, int iterationCounter)
        {
            int spatialDim = forces.Length;
            double[] ForcesAndTorque = new double[spatialDim + 1];
            double[] underrelaxationCoeff = new double[spatialDim + 1];
            double averageForce = CalculateAverageForces(forces, torque, averageDistance);

            for (int d = 0; d < spatialDim; d++)
            {
                if (UseAdaptiveUnderrelaxation == true)
                {
                    underrelaxationCoeff[d] = CalculateAdaptiveUnderrelaxation(forces[d], forcesAtPrevIteration[d], averageForce, convergenceLimit, iterationCounter, RelaxationFactor);
                }
                else
                {
                    underrelaxationCoeff[d] = RelaxationFactor;
                }
            }
            if (UseAdaptiveUnderrelaxation == true)
            {
                underrelaxationCoeff[spatialDim] = CalculateAdaptiveUnderrelaxation(torque, torqueAtPrevIteration, averageForce, convergenceLimit, iterationCounter, RelaxationFactor);
            }
            else
            {
                underrelaxationCoeff[spatialDim] = RelaxationFactor;
            }
            Console.WriteLine("ForcesUnderrelaxation[0]  " + underrelaxationCoeff[0] + ", ForcesUnderrelaxation[1]: " + underrelaxationCoeff[1] + ", TorqueUnderrelaxation " + underrelaxationCoeff[spatialDim]);
            Console.WriteLine("tempfForces[0]  " + forces[0] + ", temp_Forces[1]: " + forces[1] + ", tempTorque " + torque);
            for (int d = 0; d < spatialDim; d++)
            {
                ForcesAndTorque[d] = underrelaxationCoeff[d] * forces[d] + (1 - underrelaxationCoeff[d]) * forcesAtPrevIteration[d];
                if (Math.Abs(ForcesAndTorque[d]) < convergenceLimit * 1e-2 && clearSmallValues == true)
                {
                    ForcesAndTorque[d] = 0;
                }
            }

            ForcesAndTorque[spatialDim] = underrelaxationCoeff[spatialDim] * torque + (1 - underrelaxationCoeff[spatialDim]) * torqueAtPrevIteration;
            if (Math.Abs(ForcesAndTorque[spatialDim]) < convergenceLimit * 1e-2 && clearSmallValues == true)
            {
                ForcesAndTorque[spatialDim] = 0;
            }

            return ForcesAndTorque;
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
        private double CalculateAverageForces(double[] forces, double torque, double averageDistance = 1)
        {
            double averageForces = Math.Abs(torque) / averageDistance;
            for (int d = 0; d < forces.Length; d++)
            {
                averageForces += forces[d];
            }
            return averageForces / 3;
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
        private double CalculateAdaptiveUnderrelaxation(double variable, double variableAtPrevIteration, double averageValueOfVar, double convergenceLimit, int iterationCounter, double predefinedFactor)
        {
            int IterationHelper = 1;
            double ConvergenceHelperFactor = 1;
            double UnderrelaxationCoeff = predefinedFactor * 1e-1;
            double UnderrelaxationExponent = 0;

            if (iterationCounter / IterationHelper == 30000 * IterationHelper)
            {
                IterationHelper += 1;
                ConvergenceHelperFactor *= 1e-1;
                if (ConvergenceHelperFactor < 1e-5)
                    throw new ArithmeticException("I can not reach convergence even with very small underrelaxation factors");
                
            }

            while (Math.Abs(UnderrelaxationCoeff * variable) > 0.75 * Math.Abs(variableAtPrevIteration) && UnderrelaxationCoeff > 1e-20)
            {
                UnderrelaxationExponent -= 1;
                UnderrelaxationCoeff = predefinedFactor * Math.Pow(10, UnderrelaxationExponent);
            }

            //if (Math.Abs(UnderrelaxationCoeff * variable) < convergenceLimit * 100 && 10000 * Math.Abs(variable) > Math.Abs(averageValueOfVar))
            //    UnderrelaxationCoeff = convergenceLimit * 10;

            if (UnderrelaxationCoeff >= predefinedFactor * 1e-1)
                UnderrelaxationCoeff = predefinedFactor * 1e-1;

            double GlobalStateBuffer = UnderrelaxationCoeff.MPIMin();
            UnderrelaxationCoeff = GlobalStateBuffer;

            return UnderrelaxationCoeff * ConvergenceHelperFactor;
        }
    }
}
