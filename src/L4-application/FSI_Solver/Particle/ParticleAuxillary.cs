using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    class ParticleAuxillary {

        internal void SaveValueOfLastTimestep(List<double> variable) {
            variable.Insert(0, new double());
            variable[0] = 0;
            variable.RemoveAt(variable.Count - 1);
        }

        internal void SaveMultidimValueOfLastTimestep(List<double[]> variable) {
            int Dim = variable[0].Length;
            variable.Insert(0, new double[Dim]);
            for (int d = 0; d < Dim; d++)
            {
                variable[0][d] = 0;
            }
            variable.RemoveAt(variable.Count - 1);
        }

        static internal double ForceTorqueSummationWithNeumaierArray(double ForcesTorque, MultidimensionalArray Summands, double Length)
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
