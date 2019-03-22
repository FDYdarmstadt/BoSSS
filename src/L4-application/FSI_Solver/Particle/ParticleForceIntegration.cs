using BoSSS.Application.FSI_Solver;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    class ParticleForceIntegration
    {
        ParticleAuxillary Aux = new ParticleAuxillary();
        public double CalculateStressTensorX(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j)
        {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
            return ParticleAuxillary.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
        }

        public double CalculateStressTensorY(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j)
        {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
            return ParticleAuxillary.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
        }

        public double CalculateStressTensor2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j, int currentDimension)
        {
            double acc;
            switch (currentDimension)
            {
                case 0:
                    acc = CalculateStressTensorX(Grad_UARes, pARes, NormalVector, muA, k, j);
                    break;

                case 1:
                    acc = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, muA, k, j);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        public double CalculateTorqueFromStressTensor2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, MultidimensionalArray NodeSetClone, double muA, int k, int j, double[] currentPosition)
        {
            double temp1;
            double temp2;
            temp1 = CalculateStressTensorX(Grad_UARes, pARes, NormalVector, muA, k, j);
            temp1 *= -NormalVector[j, k, 1] * (currentPosition[1] - NodeSetClone[k, 1]).Abs();
            if (double.IsNaN(temp1) || double.IsInfinity(temp1))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            temp2 = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, muA, k, j);
            temp2 *= NormalVector[j, k, 0] * (currentPosition[0] - NodeSetClone[k, 0]).Abs();
            if (double.IsNaN(temp2) || double.IsInfinity(temp2))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            return temp1 + temp2;
        }

        public double CalculateStressTensor3D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j, int currentDimension)
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
                    acc += ParticleAuxillary.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                    break;
                case 1:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 2];
                    acc += ParticleAuxillary.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                    break;
                case 2:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 2];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 2, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 1];
                    acc += ParticleAuxillary.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        public double CalculateStressTensor(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j, int Dimensionality, int currentDimension)
        {
            double temp;
            switch (Dimensionality)
            {
                case 2:
                    temp = CalculateStressTensor2D(Grad_UARes, pARes, NormalVector, muA, k, j, currentDimension);
                    break;
                case 3:
                    temp = CalculateStressTensor3D(Grad_UARes, pARes, NormalVector, muA, k, j, currentDimension);
                    break;
                default:
                    throw new NotSupportedException("Unknown particle dimension: m_Dim = " + Dimensionality);
            }
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to calculate the particle stress tensor");
            return temp;
        }

        
    }
}
