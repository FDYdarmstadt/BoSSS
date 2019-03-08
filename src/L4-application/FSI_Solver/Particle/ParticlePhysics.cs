using BoSSS.Application.FSI_Solver;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    class ParticlePhysics
    {
        ParticleAuxillary aux = new ParticleAuxillary();

        public double CalculateStressTensorX(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j)
        {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
            return aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
        }

        public double CalculateStressTensorY(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j)
        {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
            return aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
        }

        public double CalculateStressTensor2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j)
        {
            int spatialDim = 2;
            double acc;
            switch (spatialDim)
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
            temp2 = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, muA, k, j);
            temp2 *= NormalVector[j, k, 0] * (currentPosition[0] - NodeSetClone[k, 0]).Abs();
            return temp1 + temp2;
        }

        public double CalculateStressTensor3D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double muA, int k, int j)
        {
            int spatialDim = 3;
            double acc = 0.0;
            double[] SummandsVelGradient = new double[5];
            double SummandsPressure;
            switch (spatialDim)
            {
                case 0:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 2];
                    acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                    break;
                case 1:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 2];
                    acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                    break;
                case 2:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 2];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 2, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 1];
                    acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }
    }
}
