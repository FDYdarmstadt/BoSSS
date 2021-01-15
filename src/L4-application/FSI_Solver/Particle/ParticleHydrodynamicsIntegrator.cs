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
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    public class ParticleHydrodynamicsIntegration {

        /// <summary>
        /// Constructor for the integrator. Calculates the forces and torque at the particle boundary.
        /// </summary>
        /// <param name="SpatialDim">Spatial dimension of the problem</param>
        /// <param name="U">Velocity field</param>
        /// <param name="P">Pressure field</param>
        /// <param name="LevelSetTracker">The level set tracker</param>
        /// <param name="FluidViscosity">FLuid viscosity</param>
        public ParticleHydrodynamicsIntegration(int SpatialDim, VectorField<XDGField> U, XDGField P, LevelSetTracker LevelSetTracker, double FluidViscosity) {
            this.SpatialDim = SpatialDim;
            RequiredOrder = U[0].Basis.Degree * 3 + 2;
            this.U = U.ToArray();
            this.P = P;
            this.LevelSetTracker = LevelSetTracker;
            GridData = this.LevelSetTracker.GridDat;
            this.FluidViscosity = FluidViscosity;
        }

        [DataMember]
        private readonly int SpatialDim;
        [DataMember]
        private readonly int RequiredOrder;
        [DataMember]
        private readonly XDGField[] U;
        [DataMember]
        private readonly XDGField P;
        [DataMember]
        private readonly LevelSetTracker LevelSetTracker;
        [DataMember]
        private readonly GridData GridData;
        [DataMember]
        private readonly double FluidViscosity;

        /// <summary>
        /// Calculate Forces acting from fluid onto the particle
        /// </summary>
        /// <param name="CutCells">
        /// Cut cells of a single particle.
        /// </param>
        internal double[] Forces(CellMask CutCells) {
            double[] tempForces = new double[SpatialDim];
            double[] IntegrationForces = tempForces.CloneAs();
            for (int d = 0; d < SpatialDim; d++) {
                void ErrFunc(int CurrentCellID, int Length, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1);
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Length, K, SpatialDim, SpatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Length, K);
                    MultidimensionalArray Normals = LevelSetTracker.DataHistories[0].Current.GetLevelSetNormals(Ns, CurrentCellID, Length);
                    for (int i = 0; i < SpatialDim; i++) {
                        U[i].EvaluateGradient(CurrentCellID, Length, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    P.Evaluate(CurrentCellID, Length, Ns, pARes);
                    for (int j = 0; j < Length; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = CalculateStressVector(Grad_UARes, pARes, Normals, FluidViscosity, k, j, SpatialDim, d);
                        }
                    }
                }

                int[] noOfIntegrals = new int[] { 1 };
                XQuadSchemeHelper SchemeHelper = LevelSetTracker.GetXDGSpaceMetrics(new[] { LevelSetTracker.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, CutCells);
                ICompositeQuadRule<QuadRule> surfaceRule = cqs.Compile(GridData, RequiredOrder);

                CellQuadrature.GetQuadrature(noOfIntegrals, GridData, surfaceRule,
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0)); },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { IntegrationForces[d] = ForceTorqueSummationWithNeumaierArray(IntegrationForces[d], ResultsOfIntegration, Length); }
                ).Execute();
            }
            return tempForces = IntegrationForces.CloneAs();
        }

        /// <summary>
        /// Calculate Forces acting from fluid onto the particle
        /// </summary>
        /// <param name="Position">Position of the centre of mass of the current particle</param>
        /// <param name="CutCells">Cut cells of a single particle.</param>
        internal double Torque(double[] Position, CellMask CutCells) {
            double tempTorque = new double();
            void ErrFunc2(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1);
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, SpatialDim, SpatialDim); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray Normals = LevelSetTracker.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                for (int i = 0; i < SpatialDim; i++) {
                    U[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }
                P.Evaluate(j0, Len, Ns, pARes);
                for (int j = 0; j < Len; j++) {
                    MultidimensionalArray Ns_Global = Ns.CloneAs();
                    LevelSetTracker.GridDat.TransformLocal2Global(Ns, Ns_Global, j0 + j);
                    for (int k = 0; k < K; k++) {
                        result[j, k] = CalculateTorque(Grad_UARes, pARes, Normals, Ns_Global, FluidViscosity, k, j, Position);
                    }
                }
            }
            var SchemeHelper2 = LevelSetTracker.GetXDGSpaceMetrics(new[] { LevelSetTracker.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, CutCells);
            CellQuadrature.GetQuadrature(new int[] { 1 }, GridData, cqs2.Compile(GridData, RequiredOrder),
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
        private double CalculateStressVector(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int SpatialDim, int currentDimension) {
            double stressVector;
            switch (SpatialDim) {
                case 2:
                    stressVector = CalculateStressVector2D(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j, currentDimension);
                    break;
                case 3:
                    stressVector = CalculateStressTensor3D(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j, currentDimension);
                    break;
                default:
                    throw new NotSupportedException("Unknown particle dimension: SpatialDim = " + SpatialDim);
            }
            if (double.IsNaN(stressVector) || double.IsInfinity(stressVector))
                throw new ArithmeticException("Error trying to calculate the particle stress tensor");
            return stressVector;
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
        private double CalculateTorque(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, MultidimensionalArray NodeSetClone, double FluidViscosity, int k, int j, double[] currentPosition) {
            double torqueFromXStress;
            double torqueFromYStress;
            torqueFromXStress = CalculateStressVectorX(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
            torqueFromXStress *= -(NodeSetClone[k, 1] - currentPosition[1]);
            if (double.IsNaN(torqueFromXStress) || double.IsInfinity(torqueFromXStress))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            torqueFromYStress = CalculateStressVectorY(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
            torqueFromYStress *= (NodeSetClone[k, 0] - currentPosition[0]);
            if (double.IsNaN(torqueFromYStress) || double.IsInfinity(torqueFromYStress))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            return torqueFromXStress + torqueFromYStress;
        }

        /// <summary>
        /// This method calculates the stress tensor in case of a 2D-problem
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
        private double CalculateStressVector2D(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j, int currentDimension) {
            double acc;
            switch (currentDimension) {
                case 0:
                    acc = CalculateStressVectorX(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
                    break;
                case 1:
                    acc = CalculateStressVectorY(Grad_UARes, pARes, NormalVector, FluidViscosity, k, j);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return acc;
        }

        /// <summary>
        /// This method performs the integration in x-direction
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
        private double CalculateStressVectorX(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j) {
            double[] SummandsVelGradient = new double[3];
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 1];
            double SummandsPressure = pARes[j, k] * NormalVector[j, k, 0];
            return NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
        }

        /// <summary>
        /// This method performs the integration in y-direction
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
        private double CalculateStressVectorY(MultidimensionalArray Grad_UARes, MultidimensionalArray pARes, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j) {
            double[] SummandsVelGradient = new double[3];
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * NormalVector[j, k, 1];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 0] * NormalVector[j, k, 0];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * NormalVector[j, k, 0];
            double SummandsPressure = pARes[j, k] * NormalVector[j, k, 1];
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
        /// <param name="FluidViscosity">
        /// The fluid viscosity.
        /// </param>
        private double NeumaierSummation(double[] SummandsVelGradient, double SummandsPressure, double FluidViscosity) {
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
            sum *= FluidViscosity;
            c *= FluidViscosity;
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
