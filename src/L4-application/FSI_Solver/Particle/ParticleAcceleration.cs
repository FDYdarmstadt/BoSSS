using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    class ParticleAcceleration
    {
        internal double[,] CalculateCoefficients(double[,] AddedDampingTensor, double ParticleMass, double ParticleIntertalMoment, double Timestep, double AddedDampingCoefficient = 1)
        {
            double[,] temp = new double[3, 3];
            double[,] MassMatrix = new double[3, 3];
            MassMatrix[0, 0] = MassMatrix[1, 1] = ParticleMass;
            MassMatrix[2, 2] = ParticleIntertalMoment;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    temp[i, j] = MassMatrix[i,j] + Timestep * AddedDampingCoefficient * AddedDampingTensor[i, j];
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
            temp[0] = (CoefficientMatrix[0, 1] * CoefficientMatrix[1, 2] * Torque - CoefficientMatrix[0, 1] * CoefficientMatrix[2, 2] * Forces[1] - CoefficientMatrix[0, 2] * CoefficientMatrix[1, 1] * Torque + CoefficientMatrix[0, 2] * CoefficientMatrix[2, 1] * Forces[1] + CoefficientMatrix[1, 1] * CoefficientMatrix[2, 2] * Forces[0] - CoefficientMatrix[1, 2] * CoefficientMatrix[2, 1] * Forces[0]) / Denominator;
            temp[1] = -(CoefficientMatrix[0, 0] * CoefficientMatrix[1, 2] * Torque - CoefficientMatrix[0, 0] * CoefficientMatrix[2, 2] * Forces[1] - CoefficientMatrix[0, 2] * CoefficientMatrix[1, 0] * Torque + CoefficientMatrix[0, 2] * CoefficientMatrix[2, 0] * Forces[1] + CoefficientMatrix[1, 0] * CoefficientMatrix[2, 2] * Forces[0] - CoefficientMatrix[1, 2] * CoefficientMatrix[2, 0] * Forces[0]) / Denominator;
            return temp;
        }

        internal double Rotational(double[,] CoefficientMatrix, double Denominator, double[] Forces, double Torque)
        {
            return (CoefficientMatrix[0, 0] * CoefficientMatrix[1, 1] * Torque - CoefficientMatrix[0, 0] * CoefficientMatrix[2, 1] * Forces[1] - CoefficientMatrix[0, 1] * CoefficientMatrix[1, 0] * Torque + CoefficientMatrix[0, 1] * CoefficientMatrix[2, 0] * Forces[1] + CoefficientMatrix[1, 0] * CoefficientMatrix[2, 1] * Forces[0] - CoefficientMatrix[1, 1] * CoefficientMatrix[2, 0] * Forces[0]) / Denominator;
        }
    }
}
