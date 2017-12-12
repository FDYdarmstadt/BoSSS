using System;

using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;

namespace BoSSS.Application.CDG_ProjectionTest {


    class CDGprojectionMain : BoSSS.Solution.Application{

        static void Main(string[] args) {
            var AUT = new BoSSS.Application.CDG_ProjectionTest.AllUpTest();
            AUT.TestFixtureSetUp();
            AUT.AllUp();
            AUT.TestFixtureTearDown();
        }


        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if (MPISize <= 1)
                return null;
            //int i0 = this.Grid.CellPartitioning.i0;
            //int iE = this.Grid.CellPartitioning.iE;
            //int[] CellDist = new int[iE - i0];
            //for (int i = i0; i < iE; i++) {
            //    switch (MPIRank) {
            //        case 0: {
            //                if (i == i0) {
            //                    CellDist[i - i0] = 1;
            //                } else {
            //                    CellDist[i - i0] = 0;
            //                }
            //                break;
            //            }
            //        case 1: {
            //                CellDist[i - i0] = 1;
            //                break;
            //            }
            //    }
            //}
            //return CellDist;

            return base.ComputeNewCellDistribution(TimeStepNo, physTime);
        }


        protected override GridCommons CreateOrLoadGrid() {

            double[] Xnodes = GenericBlas.Linspace(0, 4, 3);
            double[] Ynodes = GenericBlas.Linspace(0, 4, 3);
            var grid = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

            //double[] Xnodes = GenericBlas.Linspace(0, 3, 5);
            //double[] Ynodes = GenericBlas.Linspace(0, 3, 5);
            //double[] Znodes = GenericBlas.Linspace(0, 3, 2);
            //var grid = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

            //int MeshPara = 32;
            //double[] nodesX = GenericBlas.Linspace(-2, 2, MeshPara + 1);
            //double[] nodesY = GenericBlas.Linspace(-3, 3, MeshPara / 2 + 1);
            //var grid = Grid2D.Cartesian2DGrid(nodesX, nodesY);

            //int kelem = 9;
            //double[] Xnodes = GenericBlas.Linspace(0, 1, kelem + 1);
            //double[] Ynodes = GenericBlas.Linspace(-2, 2, kelem + 1);
            //var grid = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX:true, periodicY:true);

            //int kelem = 9;
            //double[] Xnodes = Grid1D.TanhSpacing(0, 6, 9, 1.5, false);
            //double[] Ynodes = Grid1D.TanhSpacing(0, 3, 9, 1.5, false);
            //var grid = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

            //var box_outer_p1 = new double[2] { 0, 0 };
            //var box_outer_p2 = new double[2] { 3, 3 };
            //var box_outer = new GridCommons.GridBox(box_outer_p1, box_outer_p2, 3, 3);

            //var box_inner_p1 = new double[2] { 1, 1 };
            //var box_inner_p2 = new double[2] { 2, 2 };
            //var box_inner = new GridCommons.GridBox(box_inner_p1, box_inner_p2, 2, 2);

            //var grid = Grid2D.HangingNodes2D(box_outer, box_inner);

            //this.m_GridPartitioningType = GridPartType.ParMETIS;

            return grid;

        }

        SinglePhaseField origin;
        SinglePhaseField result;

        ContinuousDGField cdgField;

        SpecFemField specField;
        SinglePhaseField specFieldDG;

        protected override void CreateFields() {

            int degree = 1;

            //SpecFemBasis specBasis = new SpecFemBasis(this.GridData, degree);

            //Basis dgBasis = specBasis.ContainingDGBasis;
            //Basis cdgBasis = dgBasis;

            //specField = new SpecFemField(specBasis);
            //specFieldDG = new SinglePhaseField(dgBasis, "specFEM");

            Basis dgBasis = new Basis(this.GridData, degree);
            Basis cdgBasis = new Basis(this.GridData, degree);

            origin = new SinglePhaseField(dgBasis, "origin");
            cdgField = new ContinuousDGField(cdgBasis);
            result = new SinglePhaseField(cdgBasis, "result");

           

        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
           
        }

        public bool passed = false;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            origin.ProjectField((x, y) => (Math.Sin(x) + Math.Cos(x) + x - (Math.Cos(y) + 1)));   // 2D
            //origin.ProjectField((x, y, z) => (Math.Sin(x) + Math.Cos(x) + x - (y + 1) + Math.Sin(z))); // 3D
            //origin.ProjectField((x, y, z) => Math.Sqrt(x.Pow2() + y.Pow2() + z.Pow2()) - 1);
            //origin.ProjectField((x, y) => x * x + y * y * y - x * y);
            //origin.ProjectField((x,y) => x + y);
            //origin.ProjectField((x, y) => Math.Sin(2 * Math.PI * (x / 3.0)));


            CellMask msk2D = CellMask.GetCellMask(this.GridData, X => (X[0] > 0.0 && X[0] < 4.0 && X[1] > 0.0 && X[1] < 1.0));
            //|| (X[0] > 1.0 && X[0] < 3.0 && X[1] > 1.0 && X[1] < 2.0)
            //|| (X[0] > 2.0 && X[0] < 4.0 && X[1] > 2.0 && X[1] < 3.0)
            //|| (X[0] > 3.0 && X[0] < 4.0 && X[1] > 3.0 && X[1] < 4.0));
            //CellMask msk3D = CellMask.GetCellMask(this.GridData, X => (X[0] > 0.0 && X[0] < 3.0 && X[1] > 0.0 && X[1] < 3.0 && X[2] > 0.0 && X[2] < 1.5)
            //|| (X[0] > 0.0 && X[0] < 1.5 && X[1] > 0.0 && X[1] < 3.0 && X[2] > 0.0 && X[2] < 3.0)
            //|| (X[0] > 0.0 && X[0] < 3.0 && X[1] > 0.0 && X[1] < 1.5 && X[2] > 0.0 && X[2] < 3.0));
            CellMask test = null;

            cdgField.ProjectDGField(1.0, origin, test);
            cdgField.AccToDGField(1.0, result);

            //specField.ProjectDGField(1.0, origin, test);
            //specField.AccToDGField(1.0, specFieldDG);

            //MultidimensionalArray n4 = MultidimensionalArray.Create(3, 2);
            //n4[0, 0] = 0.3;
            //n4[0, 1] = -1;
            //n4[1, 0] = 0.6;
            //n4[1, 1] = -1;
            //n4[2, 0] = 0.9;
            //n4[2, 1] = -1;

            //MultidimensionalArray n11 = MultidimensionalArray.Create(3, 2);
            //n11[0, 0] = -0.4;
            //n11[0, 1] = 1;
            //n11[1, 0] = 0.2;
            //n11[1, 1] = 1;
            //n11[2, 0] = 0.8;
            //n11[2, 1] = 1;

            //NodeSet ns4 = new NodeSet(Grid.GetRefElement(0), n4);
            //NodeSet ns11 = new NodeSet(Grid.GetRefElement(0), n11);

            //MultidimensionalArray res4 = MultidimensionalArray.Create(1, 3);
            //MultidimensionalArray res11 = MultidimensionalArray.Create(1, 3);

            //result.Evaluate(4,1,ns4, res4);
            //result.Evaluate(11, 1, ns11, res11);



            PlotCurrentState(0.0, 0, 4);

            //var err = origin.CloneAs();
            //err.Acc(-1.0, result);

            double L2err = 0.0; //err.L2Norm();
            double L2jump = JumpNorm(result, test);

            Console.WriteLine("");          
            //double L2jump_specFEM = JumpNorm(specFieldDG, test);             
            //Console.WriteLine("L2 error =  " + L2err);
            Console.WriteLine("L2 Norm of [[u]] = " + L2jump);
            //Console.WriteLine("L2 Norm of [[u]] = " + L2jump_specFEM);
            if (L2err < 1.0e-10 && L2jump < 1.0e-12) {
                Console.WriteLine("Test PASSED");
                passed = true;
            } else {
                Console.WriteLine("Test FAILED");
                passed = false;
            }

            base.TerminationKey = true;
            return 0.0;

        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { this.origin, this.result }, "CDGproj_" + timestepNo, physTime, superSampling);
        }


        static double JumpNorm(DGField f, CellMask mask = null)
        {
            GridData grd = (GridData)f.GridDat;
            int D = grd.SpatialDimension;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;

            if (mask == null) {
                mask = CellMask.GetFullMask(grd);
            }
            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            f.MPIExchange();

            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { D + 1 }, grd,
                (new EdgeQuadratureScheme(true, innerEM)).Compile(grd, f.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++)
                    {
                        int iEdge = j + i0;
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if (jCell_OT >= 0)
                        {

                            int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT);

                            f.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            f.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);

                            //Console.WriteLine("Diff at edge {0} between cell {1} and cell {2}: {3}", iEdge, jCell_IN, jCell_OT, uDiff.L2Norm());
                        }
                        else
                        {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum();

            return Unorm.Sqrt();
        }


    }


}
