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
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    class ParticleUnderrelaxation {

        /// <summary>
        /// Constructor for underrelaxation. 
        /// </summary>
        /// <param name="variable">
        /// The variable to be relaxated
        /// </param>  
        /// <param name="variablePreviousIteration">
        /// The variable at the previous iteration.
        /// </param>  
        internal ParticleUnderrelaxation(double[] variable, List<double[]> variablePreviousIteration) {
            m_Variable = variable;
            m_VariablePreviousIteration = variablePreviousIteration;
        }

        [DataMember]
        private readonly double[] m_Variable;
        [DataMember]
        private readonly List<double[]> m_VariablePreviousIteration;

        /// <summary>
        /// This method underrelaxates the hydrodynamic forces and torque. The underrelaxation coefficient is static.
        /// Used for init of the Jacobian underrelaxation.
        /// </summary>
        internal double[] StaticUnderrelaxation() {
            double[] returnVariable = m_Variable.CloneAs();
            for (int d = 0; d < m_Variable.Length; d++) {
                returnVariable[d] = 0.1 * m_Variable[d] + (1 - 0.1) * m_VariablePreviousIteration[1][d];
            }
            return returnVariable;
        }

        internal double[] AitkenUnderrelaxation(ref double Omega, List<double[]> variableWithoutRelaxationPreviousIteration) {
            double[][] residual = new double[m_Variable.Length][];
            double[] residualDiff = new double[m_Variable.Length];
            double residualScalar = 0;
            double sumVariable = 0;
            for (int i = 0; i < m_Variable.Length; i++) {
                residual[i] = new double[] { (m_Variable[i] - m_VariablePreviousIteration[0][i]), (variableWithoutRelaxationPreviousIteration[1][i] - m_VariablePreviousIteration[1][i]) };
                residualDiff[i] = residual[i][0] - residual[i][1];
                residualScalar += residual[i][1] * residualDiff[i];
                sumVariable += m_Variable[i];
            }
            Omega = -Omega * residualScalar / residualDiff.L2Norm().Pow2();
            double[] outVar = m_Variable.CloneAs();
            for (int i = 0; i < m_Variable.Length; i++) {
                outVar[i] = Omega * (m_Variable[i] - m_VariablePreviousIteration[0][i]) + m_VariablePreviousIteration[0][i];
            }
            return outVar;
        }

        /// <summary>
        /// This method underrelaxates the hydrodynamic forces and torque. The Underrelaxation
        /// factor is calculated dynamically by employing a gradient procedure.
        /// </summary>
        /// <param name="variableWithoutRelaxationPreviousIteration">
        /// The hydrodynamics of the previous iterations without underrelaxation.
        /// </param>  
        internal double[] JacobianUnderrelaxation(List<double[]> variableWithoutRelaxationPreviousIteration) {
            // ursprüngliche idee: Fixed-point fluid–structure interaction solvers with dynamic relaxation Ulrich Küttler, Wolfgang Wall
            double[] returnVariable = m_Variable.CloneAs();
            double[,] jacobian = ApproximateRelaxJacobian(variableWithoutRelaxationPreviousIteration);
            double[] residual = new double[m_Variable.Length];
            double residualScalar = 0;
            double[] residualJacFirst = new double[m_Variable.Length];
            double residualJacSecond = 0;
            for (int i = 0; i < m_Variable.Length; i++) {
                residual[i] = m_Variable[i] - m_VariablePreviousIteration[0][i];
                residualScalar += residual[i].Pow2();
            }
            for (int i = 0; i < m_Variable.Length; i++) {
                for (int j = 0; j < m_Variable.Length; j++) {
                    residualJacFirst[i] += residual[j] * jacobian[j, i];
                }
            }
            for (int i = 0; i < m_Variable.Length; i++) {
                residualJacSecond += residualJacFirst[i] * residual[i];
            }
            double underrelaxationCoefficient = -residualScalar / residualJacSecond;
            for (int i = 0; i < m_Variable.Length; i++) {
                returnVariable[i] = underrelaxationCoefficient * (m_Variable[i] - m_VariablePreviousIteration[0][i]) + m_VariablePreviousIteration[0][i];
            }
            return returnVariable;
        }

        private double[,] ApproximateRelaxJacobian(List<double[]> variableWithoutRelaxationPreviousIteration) {
            double[][] residual = new double[m_Variable.Length][];
            double[] residualDiff = new double[m_Variable.Length];
            double[] forceDiff = new double[m_Variable.Length];
            for (int i = 0; i < m_Variable.Length; i++) {
                residual[i] = new double[] { (m_Variable[i] - m_VariablePreviousIteration[0][i]), (variableWithoutRelaxationPreviousIteration[1][i] - m_VariablePreviousIteration[1][i]), (variableWithoutRelaxationPreviousIteration[2][i] - m_VariablePreviousIteration[2][i]) };//n+1, n, n-1
                residualDiff[i] = residual[i][1] - residual[i][2];
                forceDiff[i] = m_VariablePreviousIteration[1][i] - m_VariablePreviousIteration[2][i];
            }
            double[,] jacobian = new double[m_Variable.Length, m_Variable.Length];
            for (int i = 0; i < m_Variable.Length; i++) {
                for (int j = 0; j < m_Variable.Length; j++) {
                    if (forceDiff[j] == 0)
                        jacobian[i, j] = 0; // residualDiff should be 0 in this case anyway
                    else
                        jacobian[i,j] = residualDiff[i] / forceDiff[j];
                }
            }
            return jacobian;
        }
    }
}
