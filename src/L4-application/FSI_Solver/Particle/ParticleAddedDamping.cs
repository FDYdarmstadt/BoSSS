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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver
{
    class ParticleAddedDamping
    {
        public double[,] IntegrationOverLevelSet(LevelSetTracker LsTrk, double muA, double rhoA, double dt, double[] currentPosition)
        {
            double[,] addedDampingTensor = new double[6, 6];
            double alpha = 0.5;
            int RequiredOrder = 2;
            for (int DampingTensorID = 0; DampingTensorID < 4; DampingTensorID++)
            {
                for (int d1 = 0; d1 < 3; d1++)
                {
                    for (int d2 = 0; d2 < 3; d2++)
                    {
                        void evalfD(int j0, int Len, NodeSet Ns, MultidimensionalArray result)
                        {
                            int K = result.GetLength(1);
                            // Normal vector
                            var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                            MultidimensionalArray NodeSetGlobal = Ns.CloneAs();

                            if (LsTrk.GridDat.SpatialDimension == 2)
                            {
                                for (int j = 0; j < Len; j++)
                                {
                                    for (int k = 0; k < K; k++)
                                    {
                                        LsTrk.GridDat.TransformLocal2Global(Ns, NodeSetGlobal, j0 + j);
                                        double dh = CalculateNormalMeshSpacing(LsTrk, Ns, Normals, j, k);
                                        double delta = dh * Math.Sqrt(rhoA) / (Math.Sqrt(alpha * muA * dt));
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
                                        switch (DampingTensorID)
                                        {
                                            case 0:
                                                result[j, k] = d1 == d2 ? (1 - NormalComponent[d1] * NormalComponent[d2]) * muA / dn : -NormalComponent[d1] * NormalComponent[d2] * muA / dn;
                                                break;
                                            case 1:
                                                if (d1 == 2 && d2 != 2)
                                                {
                                                    result[j, k] = R[1 - d2] * Math.Pow(-1, d2) * muA / dn;
                                                }
                                                else if (d1 != 2 && d2 == 2)
                                                {
                                                    result[j, k] = ((1 - NormalComponent[d1] * NormalComponent[d1]) * (-R[1 - d1]) - NormalComponent[d1] * NormalComponent[1 - d1] * R[d1]) * Math.Pow(-1, d1) * muA / dn;
                                                }
                                                else result[j, k] = 0;
                                                break;
                                            case 2:
                                                if (d2 == 2 && d1 != 2)
                                                {
                                                    result[j, k] = R[1 - d1] * Math.Pow(-1, d1) * muA / dn;
                                                }
                                                else if (d2 != 2 && d1 == 2)
                                                {
                                                    result[j, k] = ((1 - NormalComponent[d2] * NormalComponent[d2]) * (-R[1 - d2]) - NormalComponent[d2] * NormalComponent[1 - d2] * R[d2]) * Math.Pow(-1, d2) * muA / dn;
                                                }
                                                else result[j, k] = 0;
                                                break;
                                            case 3:
                                                if (d1 == d2 && d1 != 2)
                                                {
                                                    result[j, k] = R[1 - d1].Pow2() * muA / dn;
                                                }
                                                else if (d1 != d2 && d1 != 2 && d2 != 2)
                                                {
                                                    result[j, k] = -R[0] * R[1] * muA / dn;
                                                }
                                                else if (d1 == 2 && d2 == 2)
                                                {
                                                    result[j, k] = (((1 - NormalComponent[0] * NormalComponent[0]) * R[1] + NormalComponent[0] * NormalComponent[1] * R[0]) * R[1] + ((1 - NormalComponent[1] * NormalComponent[1]) * R[0] + NormalComponent[0] * NormalComponent[1] * R[1]) * R[0]) * muA / dn;
                                                }
                                                else
                                                {
                                                    result[j, k] = 0;
                                                }
                                                break;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                throw new NotImplementedException("Currently the calculation of the Damping tensors is only available for 2D");
                            }
                        }

                        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                        //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );

                        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
                        //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, ParticleCutCells);

                        CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                            cqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult)
                            {
                                evalfD(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                            },
                            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                            {
                                for (int l = 0; l < Length; l++)
                                {
                                    switch (DampingTensorID)
                                    {
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
            if (LsTrk.GridDat.SpatialDimension == 2)
            {
                return ModifyDampingTensor2D(addedDampingTensor);
            }
            else
            {
                throw new NotImplementedException("Currently the calculation of the Damping tensors is only available for 2D");
            }
        }

        private double[,] ModifyDampingTensor2D(double[,] AddedDampingTensor)
        {
            double[,] temp = new double[3, 3];
            for (int i = 0; i < 2; i++)
            {
                for (int j= 0; j < 2; j++)
                {
                    temp[i, j] = AddedDampingTensor[i, j];
                }
            }
            temp[0, 2] = AddedDampingTensor[0, 5];
            temp[1, 2] = AddedDampingTensor[1, 5];
            temp[2, 0] = AddedDampingTensor[5, 0];
            temp[2, 1] = AddedDampingTensor[5, 1];
            temp[2, 2] = AddedDampingTensor[5, 5];
            return temp;
        }

        internal double[,] RotateTensor(double CurrentAngle, double StartingAngle, double[,] AddedDampingTensor)
        {
            // form rotation matrix R=EpEp^T, where Ep is the matrix of the principle axis of inertia
            // symmetry axis are always axis of inertia:
            double[,,] Ep = new double[2, 3, 3];
            double[] angle = new double[2];
            angle[0] = StartingAngle;
            angle[1] = CurrentAngle;
            for (int i = 0; i < 2; i++)
            {
                Ep[i, 0, 0] = Math.Cos(angle[i]);
                Ep[i, 1, 0] = Math.Sin(angle[i]);
                Ep[i, 0, 1] = -Math.Sin(angle[i]);
                Ep[i, 1, 1] = Math.Cos(angle[i]);
                Ep[i, 2, 2] = 1;
            }
            double[,] R = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        R[i, j] += Ep[1, i, k] * Ep[0, j, k];
                    }
                }
            }
            double[,] temp = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        temp[i, k] += R[i, m] * AddedDampingTensor[m, k];
                    }
                }
            }
            double[,] temp2 = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        temp2[i, j] += temp[i, k] * R[j, k];
                    }
                }
            }
            return temp2;
        }

        private double CalculateNormalMeshSpacing(LevelSetTracker LsTrk, NodeSet Ns, MultidimensionalArray Normals, int CellID, int NodeID)
        {
            double CellLength = Math.Sqrt(LsTrk.GridDat.iGeomCells.GetCellVolume(CellID));
            return Math.Abs((1 - Ns[NodeID, 0]) * Normals[CellID, NodeID, 0] * CellLength);
        }
    }
}