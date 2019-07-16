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

using BoSSS.Application.FSI_Solver;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    internal class ParticleForceIntegration
    {
        /// <summary>
        /// Main method the integral over the level set to obtain the hydrodynamic forces.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="SpatialDim">
        /// The spatial dimensions.
        /// simulation).
        /// </param>
        /// <param name="currentDimension">
        /// The current dimension to be calculated.
        /// </param>
        public double CalculateStressTensor(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int SpatialDim, int currentDimension)
        {
            double temp;
            switch (SpatialDim)
            {
                case 2:
                    temp = CalculateStressTensor2D(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j, currentDimension);
                    break;
                case 3:
                    temp = CalculateStressTensor3D(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j, currentDimension);
                    break;
                default:
                    throw new NotSupportedException("Unknown particle dimension: SpatialDim = " + SpatialDim);
            }
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to calculate the particle stress tensor");
            return temp;
        }

        /// <summary>
        /// This method calculates the stress tensor in case of a 2D-probem
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="currentDimension">
        /// The current dimension to be calculated.
        /// </param>
        private double CalculateStressTensor2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int currentDimension)
        {
            double acc;
            switch (currentDimension)
            {
                case 0:
                    acc = CalculateStressTensorX(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
                    break;

                case 1:
                    acc = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        /// <summary>
        /// This method performs the integration in x-direction
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        private double CalculateStressTensorX(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j)
        {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
            return SummationWithNeumaier(SummandsVelGradient, SummandsPressure, FluidViscosity);
        }

        /// <summary>
        /// This method performs the integration in y-direction
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        private double CalculateStressTensorY(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j)
        {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
            return SummationWithNeumaier(SummandsVelGradient, SummandsPressure, FluidViscosity);
        }

        /// <summary>
        /// This method calculates the stress tensor in case of a 3D-probem
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="currentDimension">
        /// The current dimension to be calculated.
        /// </param>
        private double CalculateStressTensor3D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int currentDimension)
        {
            double acc = 0.0;
            double[] SummandsVelGradient = new double[5];
            double SummandsPressure;
            switch (currentDimension)
            {
                case 0:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 2];
                    acc += SummationWithNeumaier(SummandsVelGradient, SummandsPressure, FluidViscosity);
                    break;
                case 1:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 2];
                    acc += SummationWithNeumaier(SummandsVelGradient, SummandsPressure, FluidViscosity);
                    break;
                case 2:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 2];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 2, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 1];
                    acc += SummationWithNeumaier(SummandsVelGradient, SummandsPressure, FluidViscosity);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        /// <summary>
        /// Main method the integral over the level set to obtain the hydrodynamic torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="NodeSetClone">
        /// The node set.
        /// simulation).
        /// </param>
        /// <param name="currentPosition">
        /// The current position of the particle.
        /// </param>
        public double CalculateTorqueFromStressTensor2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, MultidimensionalArray NodeSetClone, double FluidViscosity, int k, int j, double[] currentPosition)
        {
            double temp1;
            double temp2;
            temp1 = CalculateStressTensorX(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
            temp1 *= -NormalVector[j, k, 1] * (currentPosition[1] - NodeSetClone[k, 1]).Abs();
            if (double.IsNaN(temp1) || double.IsInfinity(temp1))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            temp2 = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
            temp2 *= NormalVector[j, k, 0] * (currentPosition[0] - NodeSetClone[k, 0]).Abs();
            if (double.IsNaN(temp2) || double.IsInfinity(temp2))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            return temp1 + temp2;
        }

        /// <summary>
        /// This method performs the Neumaier algorithm form the sum of the entries of an array.
        /// It is specifically designed to sum up the velocity gradient and the pressure to 
        /// calculate the hydrodynamic forces.
        /// </summary>
        /// <param name="SummandsVelGradient">
        /// The array of the velocity gradient.
        /// </param>
        /// <param name="SummandsPressure">
        /// The pressure.
        /// </param>
        /// <param name="muA">
        /// The fluid viscosity.
        /// </param>
        static internal double SummationWithNeumaier(double[] SummandsVelGradient, double SummandsPressure, double muA)
        {
            double sum = SummandsVelGradient[0];
            double naiveSum;
            double c = 0;
            for (int i = 1; i < SummandsVelGradient.Length; i++)
            {
                naiveSum = sum + SummandsVelGradient[i];
                if (Math.Abs(sum) >= SummandsVelGradient[i])
                {
                    c += (sum - naiveSum) + SummandsVelGradient[i];
                }
                else
                {
                    c += (SummandsVelGradient[i] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            sum *= muA;
            c *= muA;
            naiveSum = sum + SummandsPressure;
            if (Math.Abs(sum) >= SummandsPressure)
            {
                c += (sum - naiveSum) + SummandsPressure;
            }
            else
            {
                c += (SummandsPressure - naiveSum) + sum;
            }
            sum = naiveSum;
            return sum + c;
        }
    }
}
