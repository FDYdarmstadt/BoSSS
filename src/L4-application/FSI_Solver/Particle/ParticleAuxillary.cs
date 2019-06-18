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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    class ParticleAuxillary {

        /// ====================================================================================
        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for onedimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        /// ====================================================================================
        internal void SaveValueOfLastTimestep(List<double> variable) {
            variable.Insert(0, new double());
            variable[0] = 0;
            variable.RemoveAt(variable.Count - 1);
        }

        /// ====================================================================================
        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for multidimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        /// ====================================================================================
        internal void SaveMultidimValueOfLastTimestep(List<double[]> variable) {
            int Dim = variable[0].Length;
            variable.Insert(0, new double[Dim]);
            for (int d = 0; d < Dim; d++)
            {
                variable[0][d] = 0;
            }
            variable.RemoveAt(variable.Count - 1);
        }

        /// ====================================================================================
        /// <summary>
        /// This method performs the Neumaier algorithm form the sum of the entries of an array.
        /// </summary>
        /// <param name="ResultVariable">
        /// The variable where  the sum will be saved.
        /// </param>
        /// <param name="Summands">
        /// The array of summands
        /// </param>
        /// <param name="Length">
        /// The number of summands.
        /// </param>
        /// ====================================================================================
        static internal double ForceTorqueSummationWithNeumaierArray(double ResultVariable, MultidimensionalArray Summands, double Length)
        {
            double sum = ResultVariable;
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
        /// ====================================================================================
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
        /// ====================================================================================
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
