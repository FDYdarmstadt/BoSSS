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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;

namespace BoSSS.Application.XNSERO_Solver {
    class ParticleAddedDamping {
        /// <summary>
        /// Calculates the added damping tensor by integrating over the level set of the particle.
        /// </summary>
        /// <param name="particle">
        /// The current particle.
        /// </param>
        /// <param name="levelSetTracker">
        /// The level set tracker.
        /// </param>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        /// <param name="dt"></param>
        /// <param name="currentPosition"></param>
        /// <returns></returns>
        internal double[,] IntegrationOverLevelSet(Particle particle, LevelSetTracker levelSetTracker, double fluidViscosity, double fluidDensity, double dt, double[] currentPosition) {
            double[,] addedDampingTensor = new double[6, 6];
            double alpha = 0.5;
            int RequiredOrder = 4;
            for (int DampingTensorID = 0; DampingTensorID < 4; DampingTensorID++) {
                for (int d1 = 0; d1 < 3; d1++) {
                    for (int d2 = 0; d2 < 3; d2++) {
                        void evalfD(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                            int K = result.GetLength(1);
                            MultidimensionalArray Normals = levelSetTracker.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                            MultidimensionalArray NodeSetGlobal = Ns.CloneAs();

                            if (levelSetTracker.GridDat.SpatialDimension == 2) {
                                for (int j = 0; j < Len; j++) {
                                    for (int k = 0; k < K; k++) {
                                        levelSetTracker.GridDat.TransformLocal2Global(Ns, NodeSetGlobal, j0 + j);
                                        double dh = CalculateNormalMeshSpacing(levelSetTracker, Ns, Normals, j, k);
                                        double delta = dh * Math.Sqrt(fluidDensity) / (Math.Sqrt(alpha * fluidViscosity * dt));
                                        double dn = dh / (1 - Math.Exp(-delta));
                                        double[] R = new double[3];
                                        R[0] = NodeSetGlobal[k, 0] - currentPosition[0];
                                        R[1] = NodeSetGlobal[k, 1] - currentPosition[1];
                                        R[2] = 0;
                                        double[] NormalComponent = new double[3];
                                        double test = NodeSetGlobal[k, 0];
                                        double test2 = NodeSetGlobal[k, 1];
                                        NormalComponent[0] = Normals[j, k, 0];
                                        NormalComponent[1] = Normals[j, k, 1];
                                        NormalComponent[2] = 0;
                                        switch (DampingTensorID) {
                                            case 0://D^{vv}
                                                result[j, k] = d1 == d2 ? (1 - NormalComponent[d1] * NormalComponent[d2]) * fluidViscosity / dn : -NormalComponent[d1] * NormalComponent[d2] * fluidViscosity / dn;
                                                break;
                                            case 1://D^{vw}
                                                if (d1 == 2 && d2 != 2) {
                                                    result[j, k] = R[1 - d2] * Math.Pow(-1, d2) * fluidViscosity / dn;
                                                }
                                                else if (d1 != 2 && d2 == 2) {
                                                    result[j, k] = ((1 - NormalComponent[d1] * NormalComponent[d1]) * (-R[1 - d1]) - NormalComponent[d1] * NormalComponent[1 - d1] * R[d1]) * Math.Pow(-1, d1) * fluidViscosity / dn;
                                                }
                                                else result[j, k] = 0;
                                                break;
                                            case 2://D^{wv}
                                                if (d2 == 2 && d1 != 2) {
                                                    result[j, k] = R[1 - d1] * Math.Pow(-1, d1) * fluidViscosity / dn;
                                                }
                                                else if (d2 != 2 && d1 == 2) {
                                                    result[j, k] = ((1 - NormalComponent[d2] * NormalComponent[d2]) * (-R[1 - d2]) - NormalComponent[d2] * NormalComponent[1 - d2] * R[d2]) * Math.Pow(-1, d2) * fluidViscosity / dn;
                                                }
                                                else result[j, k] = 0;
                                                break;
                                            case 3://D^{ww}
                                                if (d1 == d2 && d1 != 2) {
                                                    result[j, k] = R[1 - d1].Pow2() * fluidViscosity / dn;
                                                }
                                                else if (d1 != d2 && d1 != 2 && d2 != 2) {
                                                    result[j, k] = -R[0] * R[1] * fluidViscosity / dn;
                                                }
                                                else if (d1 == 2 && d2 == 2) {
                                                    result[j, k] = (((1 - NormalComponent[0] * NormalComponent[0]) * R[1] + NormalComponent[0] * NormalComponent[1] * R[0]) * R[1] + ((1 - NormalComponent[1] * NormalComponent[1]) * R[0] + NormalComponent[0] * NormalComponent[1] * R[1]) * R[0]) * fluidViscosity / dn;
                                                }
                                                else {
                                                    result[j, k] = 0;
                                                }
                                                break;
                                        }
                                    }
                                }
                            }
                            else {
                                throw new NotImplementedException("Currently the calculation of the Damping tensors is only available for 2D");
                            }
                        }

                        CellMask allCutCells = particle.LsTrk.Regions.GetCutCellMask();
                        var SchemeHelper = levelSetTracker.GetXDGSpaceMetrics(new[] { levelSetTracker.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, particle.ParticleCutCells(levelSetTracker, allCutCells));
                        CellQuadrature.GetQuadrature(new int[] { 1 }, levelSetTracker.GridDat,
                            cqs.Compile(levelSetTracker.GridDat, RequiredOrder),
                            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                                evalfD(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                            },
                            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                for (int l = 0; l < Length; l++) {
                                    switch (DampingTensorID) {
                                        case 0:
                                            addedDampingTensor[d1, d2] += ResultsOfIntegration[l, 0];
                                            break;
                                        case 1:
                                            addedDampingTensor[d1, d2 + 3] += ResultsOfIntegration[l, 0];
                                            break;
                                        case 2:
                                            addedDampingTensor[d1 + 3, d2] += ResultsOfIntegration[l, 0];
                                            break;
                                        case 3:
                                            addedDampingTensor[d1 + 3, d2 + 3] += ResultsOfIntegration[l, 0];
                                            break;
                                    }
                                }
                            }
                        ).Execute();
                    }
                }
            }
            if (levelSetTracker.GridDat.SpatialDimension == 2) {
                return ModifyDampingTensor2D(addedDampingTensor);
            }
            else {
                throw new NotImplementedException("Currently the calculation of the Damping tensors is only available for 2D");
            }
        }

        /// <summary>
        /// Rotates the added damping tensor to obtain an updated tensor.
        /// </summary>
        /// <param name="currentAngle"></param>
        /// <param name="startingAngle"></param>
        /// <param name="addedDampingTensor"></param>
        /// <returns></returns>
        internal double[,] RotateTensor(double currentAngle, double startingAngle, double[,] addedDampingTensor) {
            // form rotation matrix R=EpEp^T, where Ep is the matrix of the principle axis of inertia
            // symmetry axis are always axis of inertia:
            double[,,] principleAxisOfInertiaMatrix = new double[2, 3, 3];
            double[] angle = new double[2];
            angle[0] = startingAngle;
            angle[1] = currentAngle;
            for (int i = 0; i < 2; i++) {
                principleAxisOfInertiaMatrix[i, 0, 0] = Math.Cos(angle[i]);
                principleAxisOfInertiaMatrix[i, 1, 0] = Math.Sin(angle[i]);
                principleAxisOfInertiaMatrix[i, 0, 1] = -Math.Sin(angle[i]);
                principleAxisOfInertiaMatrix[i, 1, 1] = Math.Cos(angle[i]);
                principleAxisOfInertiaMatrix[i, 2, 2] = 1;
            }
            double[,] rotationMatrix = new double[3, 3];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        rotationMatrix[i, j] += principleAxisOfInertiaMatrix[1, i, k] * principleAxisOfInertiaMatrix[0, j, k];
                    }
                }
            }
            double[,] temp = new double[3, 3];
            for (int i = 0; i < 3; i++) {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        temp[i, k] += rotationMatrix[i, m] * addedDampingTensor[m, k];
                    }
                }
            }

            double[,] rotatedTensor = new double[3, 3];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        rotatedTensor[i, j] += temp[i, k] * rotationMatrix[j, k];
                    }
                }
            }
            return rotatedTensor;
        }

        private double[,] ModifyDampingTensor2D(double[,] AddedDampingTensor) {
            double[,] Tensor2D = new double[3, 3];
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    Tensor2D[i, j] = AddedDampingTensor[i, j];
                }
            }
            Tensor2D[0, 2] = AddedDampingTensor[0, 5];
            Tensor2D[1, 2] = AddedDampingTensor[1, 5];
            Tensor2D[2, 0] = AddedDampingTensor[5, 0];
            Tensor2D[2, 1] = AddedDampingTensor[5, 1];
            Tensor2D[2, 2] = AddedDampingTensor[5, 5];
            return Tensor2D;
        }

        private double CalculateNormalMeshSpacing(LevelSetTracker LsTrk, NodeSet Ns, MultidimensionalArray Normals, int CellID, int NodeID) {
            double CellLength = Math.Sqrt(LsTrk.GridDat.iGeomCells.GetCellVolume(CellID));
            return Math.Abs((1 - Ns[NodeID, 0]) * Normals[CellID, NodeID, 0] * CellLength);
        }
    }
}