/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using System;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System.Linq;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;

namespace BoSSS.Application.ZwoLsTest {

    /// <summary>
    /// This guy tests the basic functionality of the XDG framework
    /// if more than one level-set is involved.
    /// It is mainly relevant for surface equations (Christina's work).
    /// </summary>
    class ZwoLsTestMain : BoSSS.Solution.Application {

        private static readonly XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant
            = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            //AllUpTest.SetUp();
            //BoSSS.Application.ZwoLsTest.AllUpTest.AllUp(0.3d, 0, true);
            //AllUpTest.Teardown();
            //Assert.IsTrue(false, "Remove me");
            //return;


            BoSSS.Solution.Application._Main(
                args,
                true,
                () => new ZwoLsTestMain() { DEGREE = 1, THRESHOLD = 0.0 });
        }

        protected override GridCommons CreateOrLoadGrid() {
            var t = Triangle.Instance;


            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-3, 3, 13), GenericBlas.Linspace(-3, 3, 13));
            Array.Sort(grd.Cells, (C1, C2) => (int)(C1.GlobalID - C2.GlobalID));

            int J = grd.Cells.Length;
            int j0 = grd.CellPartitioning.i0;
            for (int jCell = 0; jCell < J; jCell++) {
                Debug.Assert(jCell + j0 == grd.Cells[jCell].GlobalID);
            }

            return grd;
        }

        public override void Init(BoSSS.Solution.Control.AppControl control) {
            control.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            base.Init(control);
        }


        LevelSet Phi0;

        LevelSet Phi1;

        SinglePhaseField u;

        SinglePhaseField du_dx;

        SinglePhaseField du_dx_Exact;

        SinglePhaseField ERR;

        XDGField XERR;

        SinglePhaseField Amarker;
        SinglePhaseField Bmarker;
        SinglePhaseField Xmarker;


        protected override void CreateFields() {
            Phi0 = new LevelSet(new Basis(this.GridData, 2), "Phi_0");
            Phi1 = new LevelSet(new Basis(this.GridData, 2), "Phi_1");


            {
                string[,] speciesTable = new string[2, 2];
                speciesTable[0, 0] = "A"; // rechter Rand von A
                speciesTable[0, 1] = "B"; // Species zwischen den LevelSets
                speciesTable[1, 0] = "X"; // 'verbotene' Species: sollte in der geg. LevelSet-Konstellation nicht vorkommen!
                speciesTable[1, 1] = "A"; // linker Rand von A
               
                base.LsTrk = new LevelSetTracker(this.GridData, MomentFittingVariant, 1, speciesTable, Phi0, Phi1);
            }

            u = new SinglePhaseField(new Basis(this.GridData, DEGREE), "U");
            du_dx = new SinglePhaseField(new Basis(this.GridData, DEGREE), "du_dx");
            du_dx_Exact = new SinglePhaseField(new Basis(this.GridData, DEGREE), "du_dx_exact");
            ERR = new SinglePhaseField(new Basis(this.GridData, DEGREE), "ERROR");
            XERR = new XDGField(new XDGBasis(LsTrk, DEGREE), "ERROR_xdg");

            Amarker = new SinglePhaseField(new Basis(this.GridData, 0), "Amarker");
            Bmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Bmarker");
            Xmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Xmarker");

            base.m_RegisteredFields.Add(Phi0);
            base.m_RegisteredFields.Add(Phi1);
            base.m_RegisteredFields.Add(u);
            base.m_RegisteredFields.Add(du_dx);
            base.m_RegisteredFields.Add(du_dx_Exact);
            base.m_RegisteredFields.Add(ERR);
            base.m_RegisteredFields.Add(XERR);
            base.m_RegisteredFields.Add(Amarker);
            base.m_RegisteredFields.Add(Bmarker);
            base.m_RegisteredFields.Add(Xmarker);

        }


        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 1;

        /// <summary>
        /// Cell Agglomeration threshold
        /// </summary>
        internal double THRESHOLD = 0.1;

        /// <summary>
        /// Turn dynamic load balancing on/off
        /// </summary>
        internal bool DYNAMIC_BALANCE = true;


        void LsUpdate(double t) {
            double offset = t;

            Console.WriteLine("LSUpdate t = " + t);


            Phi0.ProjectField((x, y) => -(x - offset).Pow2() - y.Pow2() + (0.85).Pow2());
            Phi1.ProjectField((x, y) => -(x - offset).Pow2() - y.Pow2() + (2.4).Pow2());
            LsTrk.UpdateTracker();

            if (LsTrk.Regions.GetSpeciesSubGrid("X").GlobalNoOfCells > 0)
                throw new ApplicationException("there should be no X-species");
        }


        protected override void SetInitial() {
            this.LsUpdate(0.0);

            u.ProjectField((x, y) => x);
            du_dx_Exact.ProjectField((x, y) => 1.0);

        }

        XSpatialOperator Op;


        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            Op = new XSpatialOperator(1, 0, 1, QuadOrderFunc.SumOfMaxDegrees(RoundUp: true), "u", "c1");

            Op.EquationComponents["c1"].Add(new DxFlux()); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetFlx_phi0(this.LsTrk)); // flux am lev-set 0
            Op.EquationComponents["c1"].Add(new LevSetFlx_phi1(this.LsTrk)); // flux am lev-set 1

            Op.Commit();
        }

        
        void TestAgglomeration_Extraploation(MultiphaseCellAgglomerator Agg) {

            _2D SomePolynomial;

            if (this.u.Basis.Degree == 0) {
                SomePolynomial = (x, y) => 3.2;
            } else if (this.u.Basis.Degree == 1) {
                SomePolynomial = (x, y) => 3.2 + x - 2.4*y;
            } else {
                SomePolynomial = (x, y) => x * x + y * y;
            }

            // ========================================================
            // extrapolate polynomial function for species B of a XDG-field
            // ========================================================


            var Bmask = LsTrk.Regions.GetSpeciesMask("B");

            XDGField xt = new XDGField(new XDGBasis(this.LsTrk, this.u.Basis.Degree), "XDG-test");
            xt.GetSpeciesShadowField("B").ProjectField(1.0,
            SomePolynomial.Vectorize(), 
            new CellQuadratureScheme(true, Bmask));
                      
            
            double[] bkupVec = xt.CoordinateVector.ToArray();

            Agg.ClearAgglomerated(xt.Mapping);
            Agg.Extrapolate(xt.Mapping);

            double Err = GenericBlas.L2DistPow2(bkupVec, xt.CoordinateVector).MPISum().Sqrt();
            Console.WriteLine("Agglom: extrapolation error: " + Err);
            Assert.LessOrEqual(Err, 1.0e-10);

            // ========================================================
            // extrapolate polynomial function for a single-phase-field
            // ========================================================

            SinglePhaseField xt2 = new SinglePhaseField(this.u.Basis, "DG-test");
            xt2.ProjectField(SomePolynomial);

            double[] bkupVec2 = xt2.CoordinateVector.ToArray();

            Agg.ClearAgglomerated(xt2.Mapping);
            Agg.Extrapolate(xt2.Mapping);

            double Err2 = GenericBlas.L2DistPow2(bkupVec, xt.CoordinateVector).MPISum().Sqrt();
            Console.WriteLine("Agglom: extrapolation error: " + Err2);
            Assert.LessOrEqual(Err2, 1.0e-10);
        }


        void TestAgglomeration_Projection(int quadOrder, MultiphaseCellAgglomerator Agg) {

            var Bmask = LsTrk.Regions.GetSpeciesMask("B");
            var BbitMask = Bmask.GetBitMask();
            var XDGmetrics = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, quadOrder, 1);

            int degree = Math.Max(0, this.u.Basis.Degree - 1);

            _2D SomePolynomial;
            if (degree == 0) {
                SomePolynomial = (x, y) => 3.2;
            } else if (degree == 1) {
                SomePolynomial = (x, y) => 3.2 + x - 2.4 * y;
            } else {
                SomePolynomial = (x, y) => x * x + y * y;
            }
            
            // ------------------------------------------------
            // project the polynomial onto a single-phase field
            // ------------------------------------------------
            SinglePhaseField NonAgglom = new SinglePhaseField(new Basis(this.GridData, degree), "NonAgglom");
            NonAgglom.ProjectField(1.0,
                SomePolynomial.Vectorize(),
                new CellQuadratureScheme(true, Bmask));

            // ------------------------------------------------
            // project the polynomial onto the aggomerated XDG-field
            // ------------------------------------------------

            var qsh = XDGmetrics.XQuadSchemeHelper;

            // Compute the inner product of 'SomePolynomial' with the cut-cell-basis, ...
            SinglePhaseField xt = new SinglePhaseField(NonAgglom.Basis, "test");
            xt.ProjectField(1.0,
                SomePolynomial.Vectorize(),
                qsh.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("B")).Compile(this.GridData, Agg.CutCellQuadratureOrder));
                       
            CoordinateMapping map = xt.Mapping;

            // ... change to the cut-cell-basis, ...
            double[] xVec = xt.CoordinateVector.ToArray();
            Agg.ManipulateMatrixAndRHS(default(MsrMatrix), xVec, map, null);

            {
                double[] xVec2 = xVec.CloneAs();
                Agg.ClearAgglomerated(xVec2, map); // clearing should have no effect now.

                double Err3 = GenericBlas.L2DistPow2(xVec, xVec2).MPISum().Sqrt();
                Assert.LessOrEqual(Err3, 0.0);
            }

            // ... multiply with the inverse mass-matrix, ...
            var Mfact = XDGmetrics.MassMatrixFactory;
            BlockMsrMatrix MassMtx = Mfact.GetMassMatrix(map, false);
            Agg.ManipulateMatrixAndRHS(MassMtx, default(double[]), map, map);
            var InvMassMtx = MassMtx.InvertBlocks(OnlyDiagonal: true, Subblocks: true, ignoreEmptyBlocks: true, SymmetricalInversion: true);
            double[] xAggVec = new double[xVec.Length];

            InvMassMtx.SpMV(1.0, xVec, 0.0, xAggVec);
            
            // ... and extrapolate. 
            // Since we projected a polynomial, the projection is exact and after extrapolation, must be equal
            // to the polynomial.
            xt.CoordinateVector.SetV(xAggVec);
            Agg.Extrapolate(xt.Mapping);

            // ------------------------------------------------
            // check the error
            // ------------------------------------------------


            //double Err2 = xt.L2Error(SomePolynomial.Vectorize(), new CellQuadratureScheme(domain: Bmask));
            double Err2 = NonAgglom.L2Error(xt);
            Console.WriteLine("Agglom: projection and extrapolation error: " + Err2);
            Assert.LessOrEqual(Err2, 1.0e-10);
        }


        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime);

            //phystime = 1.8;
            LsUpdate(phystime);


            // operator-matrix assemblieren
            MsrMatrix OperatorMatrix = new MsrMatrix(u.Mapping, u.Mapping);
            double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];
            MultiphaseCellAgglomerator Agg;
            MassMatrixFactory Mfact;

            // Agglomerator setup
            int quadOrder = Op.QuadOrderFunction(new int[] { u.Basis.Degree }, new int[0], new int[] { u.Basis.Degree });
            //Agg = new MultiphaseCellAgglomerator(new CutCellMetrics(MomentFittingVariant, quadOrder, LsTrk, ), this.THRESHOLD, false);
            Agg = LsTrk.GetAgglomerator(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, quadOrder, this.THRESHOLD);

            Console.WriteLine("Inter-Process agglomeration? " + Agg.GetAgglomerator(LsTrk.GetSpeciesId("B")).AggInfo.InterProcessAgglomeration);
            if (this.THRESHOLD > 0.01) {
                TestAgglomeration_Extraploation(Agg);
                TestAgglomeration_Projection(quadOrder, Agg);
            }

            // operator matrix assembly
            Op.ComputeMatrixEx(LsTrk,
                u.Mapping, null, u.Mapping,
                OperatorMatrix, Affine, false, 0.0, true,
                Agg.CellLengthScales, null,
                LsTrk.GetSpeciesId("B"));
            Agg.ManipulateMatrixAndRHS(OperatorMatrix, Affine, u.Mapping, u.Mapping);

            // mass matrix factory
            Mfact = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, quadOrder, 1).MassMatrixFactory;// new MassMatrixFactory(u.Basis, Agg);
                        
            // Mass matrix/Inverse Mass matrix
            //var MassInv = Mfact.GetMassMatrix(u.Mapping, new double[] { 1.0 }, true, LsTrk.GetSpeciesId("B"));
            var Mass = Mfact.GetMassMatrix(u.Mapping, new double[] { 1.0 }, false, LsTrk.GetSpeciesId("B"));
            Agg.ManipulateMatrixAndRHS(Mass, default(double[]), u.Mapping, u.Mapping);
            var MassInv = Mass.InvertBlocks(OnlyDiagonal: true, Subblocks: true, ignoreEmptyBlocks: true, SymmetricalInversion: false);
            

            // test that operator depends only on B-species values
            double DepTest = LsTrk.Regions.GetSpeciesSubGrid("B").TestMatrixDependency(OperatorMatrix, u.Mapping, u.Mapping);
            Console.WriteLine("Matrix dependency test: " + DepTest);
            Assert.LessOrEqual(DepTest, 0.0);

            // diagnostic output
            Console.WriteLine("Number of Agglomerations (all species): " + Agg.TotalNumberOfAgglomerations);
            Console.WriteLine("Number of Agglomerations (species 'B'): " + Agg.GetAgglomerator(LsTrk.GetSpeciesId("B")).AggInfo.SourceCells.NoOfItemsLocally.MPISum());

            // operator auswerten:
            double[] x = new double[Affine.Length];
            BLAS.daxpy(x.Length, 1.0, Affine, 1, x, 1);
            OperatorMatrix.SpMVpara(1.0, u.CoordinateVector, 1.0, x);
            MassInv.SpMV(1.0, x, 0.0, du_dx.CoordinateVector);
            Agg.GetAgglomerator(LsTrk.GetSpeciesId("B")).Extrapolate(du_dx.Mapping);


            // markieren, wo ueberhaupt A und B sind
            Bmarker.AccConstant(1.0, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            Amarker.AccConstant(+1.0, LsTrk.Regions.GetSpeciesSubGrid("A").VolumeMask);
            Xmarker.AccConstant(+1.0, LsTrk.Regions.GetSpeciesSubGrid("X").VolumeMask);

            // compute error
            ERR.Clear();
            ERR.Acc(1.0, du_dx_Exact, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            ERR.Acc(-1.0, du_dx, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            double L2Err = ERR.L2Norm(LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            Console.WriteLine("L2 Error: " + L2Err);

            XERR.Clear();
            XERR.GetSpeciesShadowField("B").Acc(1.0, ERR, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            double xL2Err = XERR.L2Norm();
            Console.WriteLine("L2 Error (in XDG space): " + xL2Err);



            // check error
            if (TimestepNo > 1) {
                if (this.THRESHOLD > 0.01) {
                    // without agglomeration, the error in very tiny cut-cells may be large over the whole cell
                    // However, the error in the XDG-space should be small under all circumstances
                    Assert.LessOrEqual(L2Err, 1.0e-6);
                }
                Assert.LessOrEqual(xL2Err, 1.0e-6);
            }

            bool IsPassed = ((L2Err <= 1.0e-6 || this.THRESHOLD <= 0.01) && xL2Err <= 1.0e-7);
            if (IsPassed) {
                Console.WriteLine("Test PASSED");
            } else {
                Console.WriteLine("Test FAILED: check errors.");
            }

            // return/Ende
            base.NoOfTimesteps = 17;
            //base.NoOfTimesteps = 2;
            dt = 0.3;
            return dt;
        }
        
        /// <summary>
        /// A dummy routine in order to test cell dynamic load balancing 
        /// (**not** a good balancing, but triggers redistributon).
        /// </summary>
        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if(DYNAMIC_BALANCE && MPISize == 4) {
                int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int[] Part = new int[J];

                for(int j = 0; j < J; j++) {
                    double[] X = this.GridData.Cells.CellCenter.GetRow(j);
                    double x = X[0];
                    double y = X[1];

                    int px = x < Math.Min(physTime * 0.7, 2) ? 0 : 1;
                    int py = y < 0 ? 0 : 1;

                    Part[j] = px*2 + py;
                }
                
                return Part;
            } else {
                return null;
            }
        }


        /// <summary>
        /// some helper for manual debugging
        /// </summary>
        private void MatrixTests(MsrMatrix OpMatrix) {
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;
            int E = this.GridData.Edges.Count;
            int[,] e2c = this.GridData.Edges.CellIndices;

            MsrMatrix ConMatrix = new MsrMatrix(new Partitioning(J));
            MsrMatrix ConMatrix2 = new MsrMatrix(new Partitioning(J));
            var map = new UnsetteledCoordinateMapping(new Basis(this.GridData, 0));
            var FConMatrix = new XSpatialOperator.SpeciesFrameMatrix<MsrMatrix>(ConMatrix2, this.LsTrk.Regions, this.LsTrk.GetSpeciesId("B"), map, map);

            int jCell0 = this.GridData.CellPartitioning.i0;



            //for (int e = 0; e < E; e++) {
            foreach (int e in this.LsTrk.Regions.GetSpeciesSubGrid("B").InnerEdgesMask.ItemEnum) {
                int j0 = e2c[e, 0];
                int j1 = e2c[e, 1];

                int j0G;
                if (j0 >= J) {
                    j0G = (int)(this.GridData.Parallel.GlobalIndicesExternalCells[j0 - J]);
                } else {
                    j0G = j0 + jCell0;
                }

                int j1G;
                if (j1 >= 0) {
                    if (j1 >= J) {
                        j1G = (int)(this.GridData.Parallel.GlobalIndicesExternalCells[j1 - J]);
                    } else {
                        j1G = j1 + jCell0;
                    }
                } else {
                    j1G = -1;
                }

                ConMatrix[j0G, j0G] += 1.0;
                FConMatrix[j0G, j0G] += 1.0;

                if (j1G >= 0) {
                    ConMatrix[j0G, j1G] += 1.0;
                    FConMatrix[j0G, j1G] += 1.0;

                    if (this.GridData.CellPartitioning.IsInLocalRange(j1G)) {
                        ConMatrix[j1G, j1G] += 1.0;
                        FConMatrix[j1G, j1G] += 1.0;
                        ConMatrix[j1G, j0G] += 1.0;
                        FConMatrix[j1G, j0G] += 1.0;
                    }
                }
            }


            for (int jCell = 0; jCell < J; jCell++) {
                Debug.Assert(jCell + jCell0 == this.GridData.Cells.GetCell(jCell).GlobalID);
            }


            if (GridData.MpiSize == 1) {
                OpMatrix.SaveToFile("matrix.bin");
                ConMatrix.SaveToFile("conMtx.bin");
                ConMatrix2.SaveToFile("conMtx2.bin");

                //OpMatrix.SaveToTextFileSparse("C:\\tmp\\matrix1.txt");
            } else {
                var compOpMatrix = MsrMatrix.LoadFromFile("matrix.bin", csMPI.Raw._COMM.WORLD, OpMatrix.RowPartitioning, OpMatrix.ColPartition);
                var compConMatrix = MsrMatrix.LoadFromFile("conMtx.bin", csMPI.Raw._COMM.WORLD, ConMatrix.RowPartitioning, ConMatrix.ColPartition);
                var compConMatrix2 = MsrMatrix.LoadFromFile("conMtx2.bin", csMPI.Raw._COMM.WORLD, ConMatrix2.RowPartitioning, ConMatrix2.ColPartition);

                //OpMatrix.SaveToTextFileSparse("C:\\tmp\\matrix2.txt");

                compOpMatrix.Acc(-1.0, OpMatrix);
                compConMatrix.Acc(-1.0, ConMatrix);
                compConMatrix2.Acc(-1.0, ConMatrix2);

                double OpErrNrm = compOpMatrix.InfNorm();
                double ConErrNrm = compConMatrix.InfNorm();
                double Con2ErrNrm = compConMatrix.InfNorm();


                Console.WriteLine("Matrix Comparison: {0}, {1}, {2}", ConErrNrm, Con2ErrNrm, OpErrNrm);

            }
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "ZwoLsTest." + timestepNo;
            Tecplot.PlotFields(new DGField[] { u, du_dx, Phi0, Phi1, Amarker, Bmarker, Xmarker, du_dx_Exact, ERR, XERR },  filename, 0, superSampling);
        }
    }
}
