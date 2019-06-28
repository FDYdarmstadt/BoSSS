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

namespace BoSSS.Application.FSI_Solver
{
    class ParticleAcceleration
    {
        internal double[,] CalculateCoefficients(double[,] AddedDampingTensor, double ParticleMass, double InertialMoment, double Timestep, double AddedDampingCoefficient = 1)
        {
            double[,] temp = new double[3, 3];
            double[,] MassMatrix = new double[3, 3];
            MassMatrix[0, 0] = MassMatrix[1, 1] = ParticleMass;
            MassMatrix[2, 2] = InertialMoment;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (i == j)
                        temp[i, j] = MassMatrix[i, j] + Timestep * AddedDampingCoefficient * AddedDampingTensor[i, j];
                    else
                        temp[i, j] = 0;
                }
            }
            return temp;
        }

        internal double CalculateDenominator(double[,] CoefficientMatrix)
        {
            return (CoefficientMatrix[0, 0] * CoefficientMatrix[1, 1] * CoefficientMatrix[2, 2]) - (CoefficientMatrix[0, 0] * CoefficientMatrix[1, 2] * CoefficientMatrix[2, 1]) - (CoefficientMatrix[0, 1] * CoefficientMatrix[1, 0] * CoefficientMatrix[2, 2]) + (CoefficientMatrix[0, 1] * CoefficientMatrix[1, 2] * CoefficientMatrix[2, 0]) + (CoefficientMatrix[0, 2] * CoefficientMatrix[1, 0] * CoefficientMatrix[2, 1]) - (CoefficientMatrix[0, 2] * CoefficientMatrix[1, 1] * CoefficientMatrix[2, 0]);
        }

        internal double[] Translational(double[,] CoefficientMatrix, double Denominator, double[] Forces, double Torque)
        {
            double[] temp = new double[2];
            temp[0] = Forces[0] * (CoefficientMatrix[1, 1] * CoefficientMatrix[2, 2] - CoefficientMatrix[1, 2] * CoefficientMatrix[2, 1]);
            temp[0] += Forces[1] * (-CoefficientMatrix[0, 1] * CoefficientMatrix[2, 2] + CoefficientMatrix[0, 2] * CoefficientMatrix[2, 1]);
            temp[0] += Torque * (CoefficientMatrix[0, 1] * CoefficientMatrix[1, 2] - CoefficientMatrix[0, 2] * CoefficientMatrix[1, 1]);
            temp[0] = temp[0] / Denominator;

            temp[1] = Forces[0] * (-CoefficientMatrix[1, 0] * CoefficientMatrix[2, 2] + CoefficientMatrix[1, 2] * CoefficientMatrix[2, 0]);
            temp[1] += Forces[1] * (CoefficientMatrix[0, 0] * CoefficientMatrix[2, 2] - CoefficientMatrix[0, 2] * CoefficientMatrix[2, 0]);
            temp[1] += Torque * (-CoefficientMatrix[0, 0] * CoefficientMatrix[1, 2] + CoefficientMatrix[0, 2] * CoefficientMatrix[1, 0]);
            temp[1] = temp[1] / Denominator;
            return temp;
        }

        internal double Rotational(double[,] CoefficientMatrix, double Denominator, double[] Forces, double Torque)
        {
            double temp = 0;
            temp += CoefficientMatrix[0, 0] * CoefficientMatrix[1, 1] * Torque - CoefficientMatrix[0, 1] * CoefficientMatrix[1, 0] * Torque;
            temp += CoefficientMatrix[0, 1] * CoefficientMatrix[2, 0] * Forces[1] - CoefficientMatrix[0, 0] * CoefficientMatrix[2, 1] * Forces[1];
            temp += CoefficientMatrix[1, 0] * CoefficientMatrix[2, 1] * Forces[0] - CoefficientMatrix[1, 1] * CoefficientMatrix[2, 0] * Forces[0];
            return temp / Denominator;
        }
    }
}
