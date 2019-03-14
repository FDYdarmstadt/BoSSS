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
        public double[,] IntegrationOverLevelSet(int DampingTensorID, LevelSetTracker LsTrk, double muA, double rhoA, double dt, double[] currentPosition, CellMask ParticleCutCells, bool neglectAddedDamping = true)
        {
            
            if (neglectAddedDamping == true)
            {
                return null;
            }
            int D = LsTrk.GridDat.SpatialDimension;
            double[,] addedDampingTensor = new double[D, D];
            double alpha = 0.5;
            int RequiredOrder = 2;
            for (int d1 = 0; d1 < D; d1++)
            {
                for (int d2 = 0; d2 < D; d2++)
                {
                    ScalarFunctionEx evalfD = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result)
                    {
                        int K = result.GetLength(1);
                        // Normal vector
                        var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                        MultidimensionalArray NodeSetClone = Ns.CloneAs();

                        if (LsTrk.GridDat.SpatialDimension == 2)
                        {
                            for (int j = 0; j < Len; j++)
                            {
                                for (int k = 0; k < K; k++)
                                {
                                    double dh = CalculateNormalMeshSpacing(LsTrk, Ns, Normals, j, k);
                                    double delta = dh * Math.Sqrt(rhoA) / (Math.Sqrt(alpha * muA * dt));
                                    double dn = dh / (1 - Math.Exp(-delta));
                                    double[] R = new double[D];
                                    R[0] = Math.Abs(NodeSetClone[k, 0] - currentPosition[0]);
                                    R[1] = Math.Abs(NodeSetClone[k, 1] - currentPosition[1]);
                                    switch (DampingTensorID)
                                    {
                                        case 0:
                                            result[j, k] = d1 == d2 ? (1 - Normals[j, k, d1] * Normals[j, k, d2]) * muA / dn : (0 - Normals[j, k, d1] * Normals[j, k, d2]) * muA / dn;
                                            break;
                                        case 1:
                                            result[j, k] = 0;
                                            break;
                                        case 2:
                                            result[j, k] = 0;
                                            break;
                                        case 3:
                                            result[j, k] = d1 == d2 ? -(R[1 - d1] * R[1 - d2]) * muA / dn : R[1 - d1] * R[1 - d2] * muA / dn;
                                            break;
                                    }
                                }
                            }
                        }
                    };

                    var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                    //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );

                    //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
                    CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, ParticleCutCells);

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
                                addedDampingTensor[d1, d2] += ResultsOfIntegration[l, 0];
                            }
                        }
                    ).Execute();
                }
            }
            return addedDampingTensor;
        }

        internal double[,] RotateTensor(double currentAngle, double[,] addedDampingTensor)
        {
            // form rotation matrix R=EpEp^T, where Ep is the matrix of the principle axis of inertia
            // symmetry axis are always axis of inertia:
            double[,] Ep = new double[2, 2];
            Ep[0, 0] = Math.Cos(currentAngle);
            Ep[1, 0] = Math.Sin(currentAngle);
            Ep[0, 1] = -Math.Sin(currentAngle);
            Ep[1, 1] = Math.Cos(currentAngle);

            double[,] R = new double[2, 2];
            R[0, 0] = Ep[0, 0].Pow2() + Ep[0, 1] * Ep[1, 0];
            R[1, 0] = Ep[0, 0] * Ep[0, 1] + Ep[0, 1] * Ep[1, 1];
            R[0, 1] = Ep[1, 0] * Ep[0, 0] + Ep[1, 1] * Ep[1, 0];
            R[1, 1] = Ep[1, 0] * Ep[0, 1] + Ep[1, 1].Pow2();

            addedDampingTensor[0, 0] = R[0, 0].Pow2() * addedDampingTensor[0, 0] + R[1, 0] * R[0, 0] * addedDampingTensor[0, 1] + R[0, 0] * R[1, 0] * addedDampingTensor[1, 0] + R[1, 0].Pow2() * addedDampingTensor[1, 1];
            addedDampingTensor[1, 0] = R[0, 0] * R[0, 1] * addedDampingTensor[0, 0] + R[1, 0] * R[0, 1] * addedDampingTensor[0, 1] + R[0, 0] * R[1, 1] * addedDampingTensor[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensor[1, 1];
            addedDampingTensor[0, 1] = R[0, 0] * R[0, 1] * addedDampingTensor[0, 0] + R[0, 0] * R[1, 1] * addedDampingTensor[0, 1] + R[0, 0] * R[0, 1] * addedDampingTensor[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensor[1, 1];
            addedDampingTensor[1, 1] = R[0, 1].Pow2() * addedDampingTensor[0, 0] + R[1, 1] * R[0, 1] * addedDampingTensor[0, 1] + R[0, 1] * R[1, 1] * addedDampingTensor[1, 0] + R[1, 1].Pow2() * addedDampingTensor[1, 1];

            return addedDampingTensor;
        }

        private double CalculateNormalMeshSpacing(LevelSetTracker LsTrk, NodeSet Ns, MultidimensionalArray Normals, int CellID, int NodeID)
        {
            double CellLength = Math.Sqrt(LsTrk.GridDat.iGeomCells.GetCellVolume(CellID));
            return Math.Abs((1 - Ns[NodeID, 0]) * Normals[CellID, NodeID, 0] * CellLength);
        }
    }
}
