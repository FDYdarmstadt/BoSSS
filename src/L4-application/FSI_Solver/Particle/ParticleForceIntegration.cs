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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    public class ParticleHydrodynamicsIntegration {

        public ParticleHydrodynamicsIntegration(int spatialDim, VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker levelSetTracker, double fluidViscosity) {
            m_SpatialDim = spatialDim;
            m_RequiredOrder = U[0].Basis.Degree * 3 + 2;
            m_U = U.ToArray();
            m_P = P;
            m_LevelSetTracker = levelSetTracker;
            m_GridData = m_LevelSetTracker.GridDat;
            m_FluidViscosity = fluidViscosity;
        }

        [DataMember]
        private readonly int m_SpatialDim;
        [DataMember]
        readonly int m_RequiredOrder;
        [DataMember]
        readonly SinglePhaseField[] m_U;
        [DataMember]
        private readonly ConventionalDGField m_P;
        [DataMember]
        private readonly LevelSetTracker m_LevelSetTracker;
        [DataMember]
        private readonly GridData m_GridData;
        [DataMember]
        private readonly double m_FluidViscosity;

        /// <summary>
        /// Calculate Forces acting from fluid onto the particle
        /// </summary>
        internal double[] Forces(out List<double[]>[] stressToPrintOut, CellMask cutCells) {
            double[] tempForces = new double[m_SpatialDim];
            double[] IntegrationForces = tempForces.CloneAs();
            List<double[]>[] stressToPrint = new List<double[]>[m_SpatialDim];
            stressToPrint[0] = new List<double[]>();
            stressToPrint[1] = new List<double[]>();
            for (int d = 0; d < m_SpatialDim; d++) {
                void ErrFunc(int CurrentCellID, int Length, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1);
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Length, K, m_SpatialDim, m_SpatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Length, K);
                    MultidimensionalArray Normals = m_LevelSetTracker.DataHistories[0].Current.GetLevelSetNormals(Ns, CurrentCellID, Length);
                    for (int i = 0; i < m_SpatialDim; i++) {
                        m_U[i].EvaluateGradient(CurrentCellID, Length, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    m_P.Evaluate(CurrentCellID, Length, Ns, pARes);
                    for (int j = 0; j < Length; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = StressTensor(Grad_UARes, pARes, Normals, m_FluidViscosity, k, j, m_SpatialDim, d);
                            double t = Math.PI * (1 - Math.Sign(Normals[j, k, 1])) / 2 + Math.Acos(Normals[j, k, 0]);
                            stressToPrint[d].Add(new double[] { t, result[j, k] });
                        }
                    }
                }

                int[] noOfIntegrals = new int[] { 1 };
                XQuadSchemeHelper SchemeHelper = m_LevelSetTracker.GetXDGSpaceMetrics(new[] { m_LevelSetTracker.GetSpeciesId("A") }, m_RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, cutCells);
                ICompositeQuadRule<QuadRule> surfaceRule = cqs.Compile(m_LevelSetTracker.GridDat, m_RequiredOrder);

                CellQuadrature.GetQuadrature(noOfIntegrals, m_GridData, surfaceRule,
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0)); },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { IntegrationForces[d] = ForceTorqueSummationWithNeumaierArray(IntegrationForces[d], ResultsOfIntegration, Length); }
                ).Execute();
            }
            stressToPrintOut = stressToPrint.CloneAs();
            return tempForces = IntegrationForces.CloneAs();
        }

        /// <summary>
        /// Calculate Forces acting from fluid onto the particle
        /// </summary>
        internal double Torque(double[] position, CellMask cutCells) {
            double tempTorque = new double();
            void ErrFunc2(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1);
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, m_SpatialDim, m_SpatialDim); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray Normals = m_LevelSetTracker.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                for (int i = 0; i < m_SpatialDim; i++) {
                    m_U[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }
                m_P.Evaluate(j0, Len, Ns, pARes);
                for (int j = 0; j < Len; j++) {
                    MultidimensionalArray Ns_Global = Ns.CloneAs();
                    m_LevelSetTracker.GridDat.TransformLocal2Global(Ns, Ns_Global, j0 + j);
                    for (int k = 0; k < K; k++) {
                        result[j, k] = TorqueStressTensor(Grad_UARes, pARes, Normals, Ns_Global, m_FluidViscosity, k, j, position);
                    }
                }
            }
            var SchemeHelper2 = m_LevelSetTracker.GetXDGSpaceMetrics(new[] { m_LevelSetTracker.GetSpeciesId("A") }, m_RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, cutCells);
            CellQuadrature.GetQuadrature(new int[] { 1 }, m_LevelSetTracker.GridDat, cqs2.Compile(m_LevelSetTracker.GridDat, m_RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    tempTorque = ForceTorqueSummationWithNeumaierArray(tempTorque, ResultsOfIntegration, Length);
                }
            ).Execute();
            return tempTorque;
        }

        /// <summary>
        /// Main method the integral over the level set to obtain the hydrodynamic forces.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="SpatialDim">
        /// The spatial dimensions.
        /// simulation).
        /// </param>
        /// <param name="currentDimension">
        /// The current dimension to be calculated.
        /// </param>
        private double StressTensor(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int SpatialDim, int currentDimension) {
            double temp;
            switch (SpatialDim) {
                case 2:
                    temp = CalculateStressTensor2D(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j, currentDimension);
                    break;
                case 3:
                    temp = CalculateStressTensor3D(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j, currentDimension);
                    break;
                default:
                    throw new NotSupportedException("Unknown particle dimension: SpatialDim = " + SpatialDim);
            }
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to calculate the particle stress tensor");
            return temp;
        }

        /// <summary>
        /// Main method the integral over the level set to obtain the hydrodynamic torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="NodeSetClone">
        /// The node set.
        /// </param>
        /// <param name="currentPosition">
        /// The current position of the particle.
        /// </param>
        private double TorqueStressTensor(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, MultidimensionalArray NodeSetClone, double FluidViscosity, int k, int j, double[] currentPosition) {
            double temp1;
            double temp2;
            temp1 = CalculateStressTensorX(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
            //temp1 *= -NormalVector[j, k, 1] * (currentPosition[1] - NodeSetClone[k, 1]).Abs();
            temp1 *= -(NodeSetClone[k, 1] - currentPosition[1]);
            if (double.IsNaN(temp1) || double.IsInfinity(temp1))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            temp2 = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
            //temp2 *= NormalVector[j, k, 0] * (currentPosition[0] - NodeSetClone[k, 0]).Abs();
            temp2 *= (NodeSetClone[k, 0] - currentPosition[0]);
            if (double.IsNaN(temp2) || double.IsInfinity(temp2))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            return temp1 + temp2;
        }

        /// <summary>
        /// This method calculates the stress tensor in case of a 2D-probem
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="currentDimension">
        /// The current dimension to be calculated.
        /// </param>
        private double CalculateStressTensor2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int currentDimension) {
            double acc;
            switch (currentDimension) {
                case 0:
                    acc = CalculateStressTensorX(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
                    break;
                case 1:
                    acc = CalculateStressTensorY(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        /// <summary>
        /// This method performs the integration in x-direction
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        private double CalculateStressTensorX(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j) {
            double[] SummandsVelGradient = new double[3];
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
            double SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
            return NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
        }

        /// <summary>
        /// This method performs the integration in y-direction
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        private double CalculateStressTensorY(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j) {
            double[] SummandsVelGradient = new double[3];
            double SummandsPressure;
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
            SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
            return NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
        }

        /// <summary>
        /// This method calculates the stress tensor in case of a 3D-probem
        /// torque.
        /// </summary>
        /// <param name="Grad_UARes">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="pARes">
        /// The pressure.
        /// </param>
        /// <param name="NormalVector">
        /// The normal vector at the current node in the current cell.
        /// </param>
        /// <param name="FluidViscosity">
        /// The viscosity of the fluid.
        /// </param>
        /// <param name="k">
        /// The current node ID.
        /// </param>
        /// <param name="j">
        /// The current cell ID
        /// </param>
        /// <param name="currentDimension">
        /// The current dimension to be calculated.
        /// </param>
        private double CalculateStressTensor3D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int currentDimension) {
            double acc = 0.0;
            double[] SummandsVelGradient = new double[5];
            double SummandsPressure;
            switch (currentDimension) {
                case 0:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 2];
                    acc += NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
                    break;
                case 1:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 2];
                    acc += NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
                    break;
                case 2:
                    SummandsPressure = pARes[j, k] * NormalVector[j, k, 2];
                    SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 2, 2] * NormalVector[j, k, 2];
                    SummandsVelGradient[1] = -Grad_UARes[j, k, 2, 0] * NormalVector[j, k, 0];
                    SummandsVelGradient[2] = -Grad_UARes[j, k, 2, 1] * NormalVector[j, k, 1];
                    SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 2] * NormalVector[j, k, 0];
                    SummandsVelGradient[4] = -Grad_UARes[j, k, 1, 2] * NormalVector[j, k, 1];
                    acc += NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        /// <summary>
        /// This method performs the Neumaier algorithm form the sum of the entries of an array.
        /// </summary>
        /// <param name="resultVariable">
        /// The variable where  the sum will be saved.
        /// </param>
        /// <param name="summands">
        /// The array of summands
        /// </param>
        /// <param name="noOfSummands">
        /// The number of summands.
        /// </param>
        private double ForceTorqueSummationWithNeumaierArray(double resultVariable, MultidimensionalArray summands, double noOfSummands) {
            double sum = resultVariable;
            double naiveSum;
            double c = 0.0;
            for (int i = 0; i < noOfSummands; i++) {
                naiveSum = sum + summands[i, 0];
                if (Math.Abs(sum) >= Math.Abs(summands[i, 0])) {
                    c += (sum - naiveSum) + summands[i, 0];
                }
                else {
                    c += (summands[i, 0] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            return sum + c;
        }

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
        private double NeumaierSummation(double[] SummandsVelGradient, double SummandsPressure, double muA) {
            double sum = SummandsVelGradient[0];
            double naiveSum;
            double c = 0;
            for (int i = 1; i < SummandsVelGradient.Length; i++) {
                naiveSum = sum + SummandsVelGradient[i];
                if (Math.Abs(sum) >= SummandsVelGradient[i]) {
                    c += (sum - naiveSum) + SummandsVelGradient[i];
                }
                else {
                    c += (SummandsVelGradient[i] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            sum *= muA;
            c *= muA;
            naiveSum = sum + SummandsPressure;
            if (Math.Abs(sum) >= SummandsPressure) {
                c += (sum - naiveSum) + SummandsPressure;
            }
            else {
                c += (SummandsPressure - naiveSum) + sum;
            }
            sum = naiveSum;
            return sum + c;
        }
    }
}
