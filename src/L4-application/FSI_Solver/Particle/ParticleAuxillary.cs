using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    class ParticleAuxillary
    {
        internal void SaveValueToList(List<double> variable, double value, int listPosition = 0)
        {
            variable.Insert(listPosition, value);
            variable.Remove(variable.Last());
        }

        internal void SaveValueOfLastTimestep(List<double> variable)
        {
            SaveValueToList(variable, variable[0], 1);
        }

        internal void SaveMultidimValueToList(List<double[]> variable, double[] value, int listPosition = 0)
        {
            variable.Insert(listPosition, value);
            variable.Remove(variable.Last());
        }

        internal void SaveMultidimValueOfLastTimestep(List<double[]> variable)
        {
            SaveMultidimValueToList(variable, variable[0], 1);
        }

        internal double ApproxTorqueForActiveParticles()
        {
            return 0;
        }

        internal double[] ApproxForcesFromActiveStress(double Circumference_P, double active_stress_P, int spatialDim, double particleAngle, double muA)
        {
            double[] forces = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++)
            {
                forces[d] = 0.0125 * Circumference_P * active_stress_P.Pow2() * Math.Cos(particleAngle) / muA;
            }
            return forces;
        }

        internal double[] ApproxForcesFromPreviousTimestep(double[] hydrodynForcesAtTimestep)
        {
            double[] forces = new double[hydrodynForcesAtTimestep.Length];
            for (int d = 0; d < hydrodynForcesAtTimestep.Length; d++)
            {
                forces[d] = 0.5 * hydrodynForcesAtTimestep[d];
            }
            return forces;
        }

        internal double CalculateAverageForces(double[] forces, double torque, double averageDistance)
        {
            double averageForces = Math.Abs(torque) / averageDistance;
            for (int d = 0; d < forces.Length; d++)
            {
                averageForces += forces[d];
            }
            return averageForces / 3;
        }

        internal double[] CalculateAdaptiveForceUnderrelaxation(double[] forces, double[] forcesAtPrevIteration, double averageForce, double convergenceLimit, int predefinedFactor = 1)
        {
            double[] UnderrelaxationCoeff = new double[forces.Length];
            double UnderrelaxationExponent = 0;
            int spatialDim = forces.Length;
            for (int d = 0; d < spatialDim; d++)
            {
                while (Math.Abs(UnderrelaxationCoeff[d] * forces[d]) > 0.75 * Math.Abs(forcesAtPrevIteration[d]))
                {
                    UnderrelaxationExponent -= 1;
                    UnderrelaxationCoeff[d] = predefinedFactor * Math.Pow(10, UnderrelaxationExponent);
                }
                if (Math.Abs(UnderrelaxationCoeff[d] * forces[d]) < convergenceLimit * 1000 && 100 * Math.Abs(forces[d]) > averageForce)
                {
                    UnderrelaxationCoeff[d] = convergenceLimit * 1000;
                }
                if (UnderrelaxationCoeff[d] >= predefinedFactor * 1e-1)
                {
                    UnderrelaxationCoeff[d] = predefinedFactor * 1e-1;
                }
            }
            return UnderrelaxationCoeff;
        }

        internal double CalculateAdaptiveTorqueUnderrelaxation(double torque, double torqueAtPrevIteration, double averageForce, double convergenceLimit, int predefinedFactor = 1)
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

        internal double RelaxatedForce(double[] underrelaxationCoeff, double[] forces, double[] forcesAtPrevIteration, bool clearSmallValues = false)
        {
            int spatialDim = forces.Length;
            for (int d = 0; d < spatialDim; d++)
            {
                forces_underR[d] = ForcesUnderrelaxation[d] * forces[d] + (1 - ForcesUnderrelaxation[d]) * hydrodynForcesAtIteration[0][d];
                // kill all forces smaller than a certain value (increases stability)
                if (Math.Abs(forces_underR[d]) < forceAndTorque_convergence * 1e-2 && deleteSmallValues == true)
                {
                    forces_underR[d] = 0;
                }
            }
            return;
        }

    }
}
