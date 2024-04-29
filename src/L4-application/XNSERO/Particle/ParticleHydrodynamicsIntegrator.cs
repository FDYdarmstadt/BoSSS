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

namespace BoSSS.Application.XNSERO_Solver {
    public class ParticleHydrodynamicsIntegration {

        /// <summary>
        /// Constructor for the integrator. Calculates the forces and torque at the particle boundary.
        /// </summary>
        /// <param name="SpatialDim">Spatial dimension of the problem</param>
        /// <param name="U">Velocity field</param>
        /// <param name="P">Pressure field</param>
        /// <param name="LevelSetTracker">The level set tracker</param>
        /// <param name="FluidViscosity">FLuid viscosity</param>
        public ParticleHydrodynamicsIntegration(int SpatialDim, XDGField[] U, XDGField P, LevelSetTracker LevelSetTracker, double[] FluidViscosity) {
            this.SpatialDim = SpatialDim;
            RequiredOrder = U[0].Basis.Degree * 3 + 2;
            this.U = U;
            this.P = P;
            this.LevelSetTracker = LevelSetTracker;
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
        public readonly LevelSetTracker LevelSetTracker;
        [DataMember]
        private readonly double[] FluidViscosity;

        /// <summary>
        /// Calculate Forces acting from fluid onto the particle
        /// </summary>
        /// <param name="CutCells">
        /// Cut cells of a single particle.
        /// </param>
        internal double[] Main(double[] Position, CellMask CutCells, string FluidSpecies) {
            int fluidSpeciesID = GetSpeciesID(FluidSpecies);
            double[] tempForces = new double[SpatialDim];
            double[] IntegrationForces = tempForces.CloneAs();
            double[] forcesAndTorque = new double[SpatialDim + 1];
            if (!CutCells.IsSubMaskOf(CellMask.GetFullMask(LevelSetTracker.GridDat))){
                throw new Exception("Cut cell mask not a sub mask of full mask.");
            }

            {
                XQuadSchemeHelper SchemeHelper = LevelSetTracker.GetXDGSpaceMetrics(new[] { LevelSetTracker.GetSpeciesId(FluidSpecies) }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(1, CutCells);
                CellQuadrature.GetQuadrature(new int[] { SpatialDim + 1 }, LevelSetTracker.GridDat,
                    cqs.Compile(LevelSetTracker.GridDat, RequiredOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray result) {
                        NodeSet Ns = QR.Nodes;
                        int K = result.GetLength(1);
                        MultidimensionalArray GradU = MultidimensionalArray.Create(Length, K, SpatialDim, SpatialDim);
                        MultidimensionalArray pressure = MultidimensionalArray.Create(Length, K);
                        MultidimensionalArray Normals = LevelSetTracker.DataHistories[1].Current.GetLevelSetNormals(Ns, i0, Length);
                        for(int d = 0; d < SpatialDim; d++) {
                            U[d].EvaluateGradient(i0, Length, Ns, GradU.ExtractSubArrayShallow(-1, -1, d, -1), 0, 1);
                        }
                        P.Evaluate(i0, Length, Ns, pressure);
                        MultidimensionalArray Ns_Global = Ns.CloneAs();
                        for(int j = 0; j < Length; j++) {
                            LevelSetTracker.GridDat.TransformLocal2Global(Ns, Ns_Global, j + i0);
                            for(int k = 0; k < K; k++) {
                                result[j, k, 0] = CalculateStressVectorX(GradU, pressure, Normals, FluidViscosity[fluidSpeciesID], k, j);
                                result[j, k, 1] = CalculateStressVectorY(GradU, pressure, Normals, FluidViscosity[fluidSpeciesID], k, j);
                                result[j, k, 2] = CalculateTorqueFromStress(new double[] { result[j, k, 0], result[j, k, 1] }, Ns_Global, k, j, Position);
                            }
                        }
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++)
                            for(int d = 0; d < SpatialDim + 1; d++) {
                                // I am quite sure that this way of using Neumair's summation is not working as intended
                                forcesAndTorque[d] = ForceTorqueSummationWithNeumaierArray(forcesAndTorque[d], ResultsOfIntegration, Length, d);
                            }
                    }
                ).Execute();
            }
            return tempForces = forcesAndTorque.CloneAs();

        }

        private static int GetSpeciesID(string FluidSpecies) {
            int fluidSpeciesID;
            if(FluidSpecies.Length == 1)// save some ms by skipping the string comparison
                return 0;
            switch(FluidSpecies) {
                case "A":
                fluidSpeciesID = 0;
                break;
                case "B":
                fluidSpeciesID = 1;
                break;
                default:
                throw new NotSupportedException("Unknown fluid species");
            }
            return fluidSpeciesID;
        }

        private double CalculateTorqueFromStress(double[] stressVector, MultidimensionalArray NodeSetClone, int k, int j, double[] currentPosition) {
            double torqueFromXStress = stressVector[0];
            double torqueFromYStress = stressVector[1];
            torqueFromXStress *= -(NodeSetClone[k, 1] - currentPosition[1]);
            if(double.IsNaN(torqueFromXStress) || double.IsInfinity(torqueFromXStress))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            torqueFromYStress *= (NodeSetClone[k, 0] - currentPosition[0]);
            if(double.IsNaN(torqueFromYStress) || double.IsInfinity(torqueFromYStress))
                throw new ArithmeticException("Error trying to calculate the particle torque");
            return torqueFromXStress + torqueFromYStress;
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
            Vector normal = new Vector(NormalVector[j, k, 0], NormalVector[j, k, 1]);
            normal = normal / normal.Abs();
            SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * normal[0];
            SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * normal[1];
            SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * normal[1];
            double SummandsPressure = pARes[j, k] * normal[0];
            return NeumaierSummation(SummandsVelGradient, SummandsPressure, FluidViscosity);
        }

        /// <summary>
        /// This method performs the integration in y-direction
        /// </summary>
        /// <param name="GradU">
        /// The gradient of the velocity.
        /// </param>
        /// <param name="Pressure">
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
        private double CalculateStressVectorY(MultidimensionalArray GradU, MultidimensionalArray Pressure, MultidimensionalArray NormalVector, double FluidViscosity, int k, int j) {
            double[] SummandsVelGradient = new double[3];
            Vector normal = new Vector(NormalVector[j, k, 0], NormalVector[j, k, 1]);
            normal = normal / normal.Abs();
            SummandsVelGradient[0] = -2 * GradU[j, k, 1, 1] * normal[1];
            SummandsVelGradient[1] = -GradU[j, k, 1, 0] * normal[0];
            SummandsVelGradient[2] = -GradU[j, k, 0, 1] * normal[0];
            double SummandsPressure = Pressure[j, k] * normal[1];
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
        private double ForceTorqueSummationWithNeumaierArray(double resultVariable, MultidimensionalArray summands, double noOfSummands, int d) {
            double sum = resultVariable;
            double naiveSum;
            double c = 0.0;
            for (int i = 0; i < noOfSummands; i++) {
                naiveSum = sum + summands[i, d];
                if (Math.Abs(sum) >= Math.Abs(summands[i, d])) {
                    c += (sum - naiveSum) + summands[i, d];
                } else {
                    c += (summands[i, d] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            return sum + c;
        }
    }
}
