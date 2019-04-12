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
        private double CalculateAverageForces(double[] forces, double torque, double averageDistance = 1)
        {
            double averageForces = Math.Abs(torque) / averageDistance;
            for (int d = 0; d < forces.Length; d++)
            {
                averageForces += forces[d];
            }
            return averageForces / 3;
        }

        private int IterationHelper = 1;
        private double ConvergenceHelperFactor = 1;

        private double CalculateAdaptiveUnderrelaxation(double variable, double variableAtPrevIteration, double averageForce, double convergenceLimit, int iterationCounter, double predefinedFactor = 1)
        {
            double UnderrelaxationCoeff = predefinedFactor * 1e-1;
            double UnderrelaxationExponent = 0;
            if (iterationCounter < 10)
            {
                IterationHelper = 1;
                ConvergenceHelperFactor = 1;
            }
            if (iterationCounter / IterationHelper == 60 * IterationHelper)
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
            if (Math.Abs(UnderrelaxationCoeff * variable) < convergenceLimit * 100 && 10000 * Math.Abs(variable) > Math.Abs(averageForce))
            {
                UnderrelaxationCoeff = convergenceLimit * 100;
            }
            if (UnderrelaxationCoeff >= predefinedFactor * 1e-1)
            {
                UnderrelaxationCoeff = predefinedFactor * 1e-1;
            }
            double GlobalStateBuffer = UnderrelaxationCoeff.MPIMin();
            UnderrelaxationCoeff = GlobalStateBuffer;
            return UnderrelaxationCoeff * ConvergenceHelperFactor;
        }

        internal double[] RelaxatedForcesAndTorque(double[] forces, double torque , double[] forcesAtPrevIteration, double torqueAtPrevIteration, double convergenceLimit, double RelaxationFactor, bool clearSmallValues, bool UseAdaptiveUnderrelaxation, double averageDistance, int iterationCounter)
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
    }
}
