using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    class ParticleAuxillary {
        internal void SaveValueToList(List<double> variable, double value, int listPosition = 0) {
            variable.Insert(listPosition, value);
            variable.RemoveAt(variable.Count - 1);
        }

        internal void SaveValueOfLastTimestep(List<double> variable) {
            SaveValueToList(variable, variable[0], 1);
        }

        internal void SaveMultidimValueToList(List<double[]> variable, double[] value, int listPosition = 0) {
            variable.Insert(listPosition, value);
            variable.RemoveAt(variable.Count - 1);
        }

        internal void SaveMultidimValueOfLastTimestep(List<double[]> variable) {
            SaveMultidimValueToList(variable, variable[0], 1);
        }

        internal double ApproxTorqueForActiveParticles() {
            return 0;
        }

        internal double[] ApproxForcesFromActiveStress(double Circumference_P, double active_stress_P, int spatialDim, double particleAngle, double muA) {
            double[] forces = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                forces[d] = 0.0125 * Circumference_P * active_stress_P.Pow2() * Math.Cos(particleAngle) / muA;
            }
            return forces;
        }

        internal double[] ApproxForcesFromPreviousTimestep(double[] hydrodynForcesAtTimestep) {
            double[] forces = new double[hydrodynForcesAtTimestep.Length];
            for (int d = 0; d < hydrodynForcesAtTimestep.Length; d++) {
                forces[d] = 0.5 * hydrodynForcesAtTimestep[d];
            }
            return forces;
        }

        internal double ForceTorqueSummationWithNeumaierArray(double ForcesTorque, MultidimensionalArray Summands, double Length)
        {
            double sum = ForcesTorque;
            double naiveSum;
            double c = 0.0;
            for (int i = 0; i < Length; i++)
            {
                naiveSum = sum + Summands[i, 0];
                if (Math.Abs(sum) >= Math.Abs(Summands[i, 0]))
                {
                    c += (sum - naiveSum) + Summands[i, 0];
                }
                else
                {
                    c += (Summands[i, 0] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            return sum + c;
        }
        internal double SummationWithNeumaier(double[] SummandsVelGradient, double SummandsPressure, double muA)
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
            // Neumaier pressure term
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
