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

        private double CalculateAdaptiveUnderrelaxation(double torque, double torqueAtPrevIteration, double averageForce, double convergenceLimit, double predefinedFactor = 1)
        {
            double UnderrelaxationCoeff = predefinedFactor * 1e-1;
            double UnderrelaxationExponent = 0;
            while (Math.Abs(UnderrelaxationCoeff * torque) > 0.75 * Math.Abs(torqueAtPrevIteration))
            {
                UnderrelaxationExponent -= 1;
                UnderrelaxationCoeff = predefinedFactor * Math.Pow(10, UnderrelaxationExponent);
            }
            if (Math.Abs(UnderrelaxationCoeff * torque) < convergenceLimit * 1000 && 100 * Math.Abs(torque) > averageForce)
            {
                UnderrelaxationCoeff = convergenceLimit * 1000;
            }
            if (UnderrelaxationCoeff >= predefinedFactor * 1e-1)
            {
                UnderrelaxationCoeff = predefinedFactor * 1e-1;
            }
            return UnderrelaxationCoeff;
        }

        internal double[] RelaxatedForcesAndTorque(double[] forces, double torque , double[] forcesAtPrevIteration, double torqueAtPrevIteration, double convergenceLimit, double RelaxationFactor, bool clearSmallValues, bool UseAdaptiveUnderrelaxation)
        {
            int spatialDim = forces.Length;
            double[] ForcesAndTorque = new double[spatialDim + 1];
            double[] underrelaxationCoeff = new double[spatialDim + 1];
            double averageForce = CalculateAverageForces(forces, torque);

            for (int d = 0; d < spatialDim; d++)
            {
                if (UseAdaptiveUnderrelaxation == true)
                {
                    underrelaxationCoeff[d] = CalculateAdaptiveUnderrelaxation(forces[d], forcesAtPrevIteration[d], averageForce, convergenceLimit, RelaxationFactor);
                }
                else
                {
                    underrelaxationCoeff[d] = RelaxationFactor;
                }
            }
            if (UseAdaptiveUnderrelaxation == true)
            {
                underrelaxationCoeff[spatialDim] = CalculateAdaptiveUnderrelaxation(torque, torqueAtPrevIteration, averageForce, convergenceLimit, RelaxationFactor);
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
