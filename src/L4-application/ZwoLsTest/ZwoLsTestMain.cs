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
using System.Collections.Generic;

namespace BoSSS.Application.ZwoLsTest {

    /// <summary>
    /// This guy tests the basic functionality of the XDG framework
    /// if more than one level-set is involved.
    /// </summary>
    class ZwoLsTestMain : BoSSS.Solution.Application {

        internal XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            InitMPI();
            BoSSS.Application.ZwoLsTest.AllUpTest.AllUp(0.0d, 1, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, true);
            //BoSSS.Application.ZwoLsTest.AllUpTest.AllUp(0.0d, 1, XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, true);
            Assert.IsTrue(false, "Remove me");
                        
            BoSSS.Solution.Application._Main(
                args,
                true,
                () => new ZwoLsTestMain() { DEGREE = 3, THRESHOLD = 0.3, MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye, DYNAMIC_BALANCE = true });
            
        }

        protected override IGrid CreateOrLoadGrid() {
            var t = Triangle.Instance;


            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-3, 3, 13), GenericBlas.Linspace(-3, 3, 13));
            Array.Sort(grd.Cells, (C1, C2) => (int)(C1.GlobalID - C2.GlobalID));

            int J = grd.Cells.Length;
            long j0 = grd.CellPartitioning.i0;
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

        bool usePhi0 = true;
        bool usePhi1 = true;

        protected override void CreateFields() {
            Phi0 = new LevelSet(new Basis(this.GridData, 2), "Phi_0");
            Phi1 = new LevelSet(new Basis(this.GridData, 2), "Phi_1");


            if(usePhi0 && usePhi1) {
                string[,] speciesTable = new string[2, 2];
                speciesTable[0, 0] = "A"; // rechter Rand von A
                speciesTable[0, 1] = "B"; // Species zwischen den LevelSets
                speciesTable[1, 0] = "X"; // 'verbotene' Species: sollte in der geg. LevelSet-Konstellation nicht vorkommen!
                speciesTable[1, 1] = "A"; // linker Rand von A

                base.LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData), MomentFittingVariant, 1, speciesTable, Phi0, Phi1);
            } else if(!usePhi0 && usePhi1) {
                string[] speciesTable = new string[2];
                speciesTable[0] = "A"; 
                speciesTable[1] = "B"; 

                base.LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData), MomentFittingVariant, 1, speciesTable, Phi1);
            } else if(usePhi0 && !usePhi1) {
                string[] speciesTable = new string[2];
                speciesTable[0] = "B"; 
                speciesTable[1] = "A"; 

                base.LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData), MomentFittingVariant, 1, speciesTable, Phi0);
            } else {
                throw new NotImplementedException();
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
        
        
        /// <summary>
        /// Triggers a comparison between single-core and MPI-parallel run, <see cref="TestLengthScales"/>
        /// </summary>
        internal bool SER_PAR_COMPARISON = false;


        void LsUpdate(double t) {
            double offset = t;

            Console.WriteLine("LSUpdate t = " + t);


            Phi0.ProjectField((x, y) => -(x - offset).Pow2() - y.Pow2() + (0.85).Pow2());
            Phi1.ProjectField((x, y) => -(x - offset).Pow2() - y.Pow2() + (2.4).Pow2());
            LsTrk.UpdateTracker(t);
            LsTrk.PushStacks();

            if(usePhi0 && usePhi1)
                if (LsTrk.Regions.GetSpeciesSubGrid("X").GlobalNoOfCells > 0)
                    throw new ApplicationException("there should be no X-species");
        }


        protected override void SetInitial(double t) {
            this.LsUpdate(t);

            u.ProjectField((x, y) => x);
            du_dx_Exact.ProjectField((x, y) => 1.0);

        }

        XDifferentialOperatorMk2 Op;


        int QuadOrder {
            get {
                return this.DEGREE * 2 + 2;
            }
        }

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            Op = new XDifferentialOperatorMk2(1, 0, 1,
                QuadOrderFunc: (int[] DomDegs, int[] ParamDegs, int[] CoDomDegs) => QuadOrder,
                __Species: new [] { "B" },
                __varnames: new[] { "u", "c1" });
            Op.AgglomerationThreshold = this.THRESHOLD;

            Op.EquationComponents["c1"].Add(new DxFlux()); // Flux in Bulk Phase;
            if(usePhi0)
                Op.EquationComponents["c1"].Add(new LevSetFlx_phi0()); // flux am lev-set 0
            if(usePhi1)
                Op.EquationComponents["c1"].Add(new LevSetFlx_phi1()); // flux am lev-set 1

            //Op.EquationComponents["c1"].Add(new DxBroken());

            Op.FluxesAreNOTMultithreadSafe = true;
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


        /// <summary>
        /// Mainly, comparison of single-core vs MPI-parallel run
        /// </summary>
        void TestLengthScales(int quadOrder, int TimestepNo) {
            string name_disc = $"t{TimestepNo}-alpha{this.THRESHOLD}-p{this.DEGREE}-q{MomentFittingVariant}";
            var spcA = LsTrk.GetSpeciesId("A");
            var spcB = LsTrk.GetSpeciesId("B");

            var species = new[] { spcA, spcB };
            //MultiphaseCellAgglomerator.CheckFile = $"InsideMpagg-{name_disc}.csv";
            MultiphaseCellAgglomerator Agg = LsTrk.GetAgglomerator(species, quadOrder, this.THRESHOLD);

            int RefMPIsize = 1;


            // check level-set coordinates
            // ===========================
            {
                var LsChecker = new TestingIO(this.GridData, $"LevelSets-{name_disc}.csv", true, RefMPIsize);
                LsChecker.AddDGField(this.Phi0);
                LsChecker.AddDGField(this.Phi1);
                LsChecker.DoIOnow();

                Assert.Less(LsChecker.AbsError(this.Phi0), 1.0e-8, "Mismatch in level-set 0 between single-core and parallel run.");
                Assert.Less(LsChecker.AbsError(this.Phi1), 1.0e-8, "Mismatch in level-set 1 between single-core and parallel run.");

            }

           
            // check equality of agglomeration
            // ===============================

            // Note: agglomeration in parallel and serial mode is not necessarily equal - but very often it is.
            //       Therefore, we need to check whether the agglomeration is equal or not,
            //       in order to know whether agglomerated length scales should be compared or not.

            bool[] equalAggAsinReferenceRun;
            {
                var aggoCheck = new TestingIO(this.GridData, $"Agglom-{name_disc}.csv", true, RefMPIsize);
                long[] GiDs = GridData.CurrentGlobalIdPermutation.Values;
                int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                long[] extGiDs = GridData.iParallel.GlobalIndicesExternalCells;

                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var spc = species[iSpc];
                    var spcN = LsTrk.GetSpeciesName(spc);
                    var ai = Agg.GetAgglomerator(spc).AggInfo;
                    var srcMask = ai.SourceCells.GetBitMask();
                    var aggPairs = ai.AgglomerationPairs;

                    aggoCheck.AddColumn($"SourceCells{spcN}", delegate (double[] X, int j, int jG) {
                        if(srcMask[j]) {
                            var pair = aggPairs.Single(cap => cap.jCellSource == j);
                            if(pair.jCellTarget < J)
                                return (double)GiDs[pair.jCellTarget];
                            else
                                return (double)extGiDs[pair.jCellTarget - J];
                        }
                        return 0.0;
                    });

                }

                aggoCheck.DoIOnow();

                equalAggAsinReferenceRun = new bool[species.Length];
                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var spc = species[iSpc];
                    var spcN = LsTrk.GetSpeciesName(spc);

                    equalAggAsinReferenceRun[iSpc] = aggoCheck.AbsError($"SourceCells{spcN}") < 0.1;
                    if(equalAggAsinReferenceRun[iSpc] == false)
                        Console.WriteLine("Different agglomeration between single-core and parallel run for species " + spcN + ".");
                }
            
                //Agg.PlotAgglomerationPairs($"-aggPairs-{name_disc}-MPI{MPIRank + 1}of{MPISize}");
            }


            // compare length scales
            // =====================

            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);

            foreach(ICutCellMetrics ccm in new ICutCellMetrics[] { Agg.NonAgglomeratedMetrics, Agg }) { // loop over non-agglom and agglomerated metrics
                string name;
                bool Agglom;
                if(object.ReferenceEquals(ccm,Agg)) {
                    name = "Agglomerated";
                    Agglom = true;
                } else {
                    name = "Nonagglom";
                    Agglom = false;
                }

                MultidimensionalArray CellSurfaceA = ccm.CellSurface[spcA];
                MultidimensionalArray CellVolumeA = ccm.CutCellVolumes[spcA];
                MultidimensionalArray CellSurfaceB = ccm.CellSurface[spcB];
                MultidimensionalArray CellVolumeB = ccm.CutCellVolumes[spcB];


                string FileName = $"{name}LengthScales-{name_disc}.csv";
                var Checker = new TestingIO(this.GridData, FileName, true, RefMPIsize);
                Checker.AddColumn("CellSurfA", (double[] X, int j, int jG) => CellSurfaceA[j]);
                Checker.AddColumn("CellVolA", (double[] X, int j, int jG) => CellVolumeA[j]);
                Checker.AddColumn("CellSurfB", (double[] X, int j, int jG) => CellSurfaceB[j]);
                Checker.AddColumn("CellVolB", (double[] X, int j, int jG) => CellVolumeB[j]);
                Checker.DoIOnow();

                if(this.MPISize == 1) {
                    var Checker2 = new TestingIO(this.GridData, FileName, true, int.MaxValue);
                    Checker2.AddColumn("CellSurfA", (double[] X, int j, int jG) => CellSurfaceA[j]);
                    Checker2.AddColumn("CellVolA", (double[] X, int j, int jG) => CellVolumeA[j]);
                    Checker2.AddColumn("CellSurfB", (double[] X, int j, int jG) => CellSurfaceB[j]);
                    Checker2.AddColumn("CellVolB", (double[] X, int j, int jG) => CellVolumeB[j]);
                    Checker2.DoIOnow();

                    foreach(string s in Checker2.ColumnNames) {
                        Assert.Less(Checker2.RelError(s), 1.0e-10, "'TestingUtils.Compare' itself is fucked up.");
                    }

                }

                double srfA = Checker.RelError("CellSurfA") * ((!Agglom || equalAggAsinReferenceRun[0]) ? 1.0 : 0.0);
                double volA = Checker.RelError("CellVolA") * ((!Agglom || equalAggAsinReferenceRun[0]) ? 1.0 : 0.0);
                double srfB = Checker.RelError("CellSurfB") * ((!Agglom || equalAggAsinReferenceRun[1]) ? 1.0 : 0.0);
                double volB = Checker.RelError("CellVolB") * ((!Agglom || equalAggAsinReferenceRun[1]) ? 1.0 : 0.0);

                if(srfA + volA + srfB + volB > 0.001) {
                  
                    Console.WriteLine($"Mismatch in {name} cell surface for species A between single-core and parallel run: {srfA}.");
                    Console.WriteLine($"Mismatch in {name} cell volume  for species A between single-core and parallel run: {volA}.");
                    Console.WriteLine($"Mismatch in {name} cell surface for species B between single-core and parallel run: {srfB}.");
                    Console.WriteLine($"Mismatch in {name} cell volume  for species B between single-core and parallel run: {volB}.");

                }


                Assert.Less(srfA, BLAS.MachineEps.Sqrt(), $"Mismatch in {name} cell surface for species A between single-core and parallel run.");
                Assert.Less(volA, BLAS.MachineEps.Sqrt(), $"Mismatch in {name} cell volume  for species A between single-core and parallel run.");
                Assert.Less(srfB, BLAS.MachineEps.Sqrt(), $"Mismatch in {name} cell surface for species B between single-core and parallel run.");
                Assert.Less(volB, BLAS.MachineEps.Sqrt(), $"Mismatch in {name} cell volume  for species B between single-core and parallel run.");
            }
        }


        /// <summary>
        /// Checks MPI exchange
        /// </summary>
        void CheckExchange(bool ExchangeShadow) {

            // its important to test more than one field!
            var xf = new XDGField(new XDGBasis(this.LsTrk, 0), "xf");
            var xf2 = new XDGField(new XDGBasis(this.LsTrk, 0), "xf2");
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int JE = this.GridData.iLogicalCells.Count;
            //long[] GiD = this.GridData.CurrentGlobalIdPermutation.Values;

            long GetGid(int j) { 
                return this.GridData.iLogicalCells.GetGlobalID(j);
            }

            {
                var st = new SinglePhaseField(new Basis(this.GridData, 1), "muh");
                for(int j = 0; j < J; j++) {
                    st.Coordinates[j, 0] = GetGid(j);
                }

                var trx = new BoSSS.Foundation.Comm.Transceiver(st);
                trx.TransceiveStartImReturn();
                trx.TransceiveFinish();


                for(int j = 0; j < JE; j++) {
                    double should = GetGid(j);
                    Assert.IsTrue(st.Coordinates[j, 0] == should, "everything is fucked");
                }

            }

            // ---------------------
            foreach(string s in this.LsTrk.SpeciesNames) {
                double sv = s[0] * 0.01;
                var xfs = xf.GetSpeciesShadowField(s);
                var xfs2 = xf2.GetSpeciesShadowField(s);
                for(int j = 0; j < J; j++) {
                    xfs.Coordinates[j, 0] = GetGid(j) + sv;
                    xfs2.Coordinates[j, 0] = 1.0e5 + GetGid(j) + sv;
                }
            }

            double[] xfB4 = xf.CoordinateVector.ToArray();
            double[] xf2B4 = xf2.CoordinateVector.ToArray();

            // ---------------------
            if(ExchangeShadow) {
                // exchange shadow fields separately
                foreach(string s in this.LsTrk.SpeciesNames) {
                    var xfs = xf.GetSpeciesShadowField(s);
                    var xfs2 = xf2.GetSpeciesShadowField(s);

                    //var xfsCopy = new SinglePhaseField(xfs.Basis, "muh");
                    //xfsCopy.Acc(1.0, xfs);
                    //dd.Add(s, xfsCopy);

                    var trx = new BoSSS.Foundation.Comm.Transceiver(xfs, xfs2);
                    trx.TransceiveStartImReturn();
                    trx.TransceiveFinish();
                }
            } else {
                // exchange XDG fields
                var trx2 = new BoSSS.Foundation.Comm.Transceiver(xf, xf2);
                trx2.TransceiveStartImReturn();
                trx2.TransceiveFinish();
            }

            foreach(string s in this.LsTrk.SpeciesNames) {
                double sv = s[0] * 0.01;
                var xfs = xf.GetSpeciesShadowField(s);
                var xfs2 = xf2.GetSpeciesShadowField(s);
                
                for(int j = 0; j < JE; j++) {
                    double should = GetGid(j) + sv;
                    double should2 = 1.0e5 + GetGid(j) + sv;

                    string nerv = ExchangeShadow ? "shadow" : "XDG";
                    Assert.IsTrue(xfs.Coordinates[j, 0] == 0.0 || xfs.Coordinates[j, 0] == should, $"error in MPI exchange of {nerv} fields");
                    Assert.IsTrue(xfs2.Coordinates[j, 0] == 0.0 || xfs2.Coordinates[j, 0] == should2, $"error in MPI exchange of {nerv} fields (2nd field)");
                }
            }
        }



        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime);

            //phystime = 1.8;
            LsUpdate(phystime);


            // operator-matrix assemblieren
            BlockMsrMatrix OperatorMatrixPar = new BlockMsrMatrix(u.Mapping, u.Mapping);
            BlockMsrMatrix OperatorMatrix = new BlockMsrMatrix(u.Mapping, u.Mapping);
            double[] AffinePar = new double[OperatorMatrixPar.RowPartitioning.LocalLength];
            double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];

            // Agglomerator setup
            MultiphaseCellAgglomerator Agg = LsTrk.GetAgglomerator(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, QuadOrder, this.THRESHOLD);

            // plausibility of cell length scales 
            if(SER_PAR_COMPARISON)
                TestLengthScales(QuadOrder, TimestepNo);

            Console.WriteLine("Inter-Process agglomeration? " + Agg.GetAgglomerator(LsTrk.GetSpeciesId("B")).AggInfo.InterProcessAgglomeration);
            if (this.THRESHOLD > 0.01) {
                TestAgglomeration_Extraploation(Agg);
                TestAgglomeration_Projection(QuadOrder, Agg);
            }
            CheckExchange(true);
            CheckExchange(false);

            // operator matrix assembly
            XDifferentialOperatorMk2.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(base.LsTrk, u.Mapping, null, u.Mapping);
            mtxBuilder.time = 0.0;
            //ilPSP.Environment.NumThreads = 1;
            mtxBuilder.ComputeMatrix(OperatorMatrixPar, AffinePar);
            int oldNumThreads = ilPSP.Environment.NumThreads;
            ilPSP.Environment.NumThreads = 1;
            mtxBuilder.ComputeMatrix(OperatorMatrix, Affine);
            ilPSP.Environment.NumThreads = oldNumThreads;

            var diff = OperatorMatrixPar.CloneAs();
            diff.Acc(-1.0, OperatorMatrix);
            double MatrixDiffNorm = diff.InfNorm();
            double AffineDiffNorm = AffinePar.MPI_L2Dist(Affine);
            Console.WriteLine($"MatrixDiffNorm: {MatrixDiffNorm:-0.####e-00}, AffineDiffNorm: {AffineDiffNorm:-0.####e-00}");

            Assert.Less(MatrixDiffNorm, 1.0e-10, "Difference in multi-thread and single-thread matrix computation.");
            Assert.Less(AffineDiffNorm, 1.0e-10, "Difference in multi-thread and single-thread affine computation.");


            Agg.ManipulateMatrixAndRHS(OperatorMatrix, Affine, u.Mapping, u.Mapping);

            // mass matrix factory
            var Mfact = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, QuadOrder, 1).MassMatrixFactory;// new MassMatrixFactory(u.Basis, Agg);
                        
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
            OperatorMatrix.SpMV(1.0, u.CoordinateVector, 1.0, x);
            MassInv.SpMV(1.0, x, 0.0, du_dx.CoordinateVector);
            Agg.GetAgglomerator(LsTrk.GetSpeciesId("B")).Extrapolate(du_dx.Mapping);

            // markieren, wo ueberhaupt A und B sind
            Bmarker.AccConstant(1.0, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            Amarker.AccConstant(+1.0, LsTrk.Regions.GetSpeciesSubGrid("A").VolumeMask);
            if(usePhi0 && usePhi1)
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
            double ErrorThreshold = 1.0e-1;
            if(this.MomentFittingVariant == XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)
                ErrorThreshold = 1.0e-6; // HMF is designed for such integrands and should perform close to machine accuracy; on general integrands, the precision is different.


            bool IsPassed = ((L2Err <= ErrorThreshold || this.THRESHOLD <= ErrorThreshold) && xL2Err <= ErrorThreshold);
            if (IsPassed) {
                Console.WriteLine("Test PASSED");
            } else {
                Console.WriteLine("Test FAILED: check errors.");
                //PlotCurrentState(phystime, TimestepNo, 3);
            }

            if (TimestepNo > 1) {
                if (this.THRESHOLD > ErrorThreshold) {
                    // without agglomeration, the error in very tiny cut-cells may be large over the whole cell
                    // However, the error in the XDG-space should be small under all circumstances
                    Assert.LessOrEqual(L2Err, ErrorThreshold, "DG L2 error of computing du_dx");
                }
                Assert.LessOrEqual(xL2Err, ErrorThreshold, "XDG L2 error of computing du_dx");
            }

            

            // return/Ende
            base.NoOfTimesteps = 17;
            //base.NoOfTimesteps = 2;
            dt = 0.3;
            LsTrk.UpdateTracker(phystime + dt); // required to ensure the internal time-levels are correct

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
                    double[] X = this.GridData.iLogicalCells.GetCenter(j);
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
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int E = this.GridData.iLogicalEdges.Count;
            int[,] e2c = this.GridData.iLogicalEdges.CellIndices;

            MsrMatrix ConMatrix = new MsrMatrix(new Partitioning(J));
            MsrMatrix ConMatrix2 = new MsrMatrix(new Partitioning(J));
            var map = new UnsetteledCoordinateMapping(new Basis(this.GridData, 0));
            var FConMatrix = new XDifferentialOperatorMk2.SpeciesFrameMatrix<MsrMatrix>(ConMatrix2, this.LsTrk.Regions, this.LsTrk.GetSpeciesId("B"), map, map);

            long jCell0 = this.GridData.CellPartitioning.i0;



            //for (int e = 0; e < E; e++) {
            foreach (int e in this.LsTrk.Regions.GetSpeciesSubGrid("B").InnerEdgesMask.ItemEnum) {
                int j0 = e2c[e, 0];
                int j1 = e2c[e, 1];

                long j0G;
                if (j0 >= J) {
                    j0G = (this.GridData.iParallel.GlobalIndicesExternalCells[j0 - J]);
                } else {
                    j0G = j0 + jCell0;
                }

                long j1G;
                if (j1 >= 0) {
                    if (j1 >= J) {
                        j1G = (this.GridData.iParallel.GlobalIndicesExternalCells[j1 - J]);
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
                Debug.Assert(jCell + jCell0 == this.GridData.iLogicalCells.GetGlobalID(jCell));
            }


            if (GridData.MpiSize == 1) {
                OpMatrix.SaveToFile("matrix.bin");
                ConMatrix.SaveToFile("conMtx.bin");
                ConMatrix2.SaveToFile("conMtx2.bin");

            } else {
                var compOpMatrix = MsrMatrix.LoadFromFile("matrix.bin", csMPI.Raw._COMM.WORLD, OpMatrix.RowPartitioning, OpMatrix.ColPartition);
                var compConMatrix = MsrMatrix.LoadFromFile("conMtx.bin", csMPI.Raw._COMM.WORLD, ConMatrix.RowPartitioning, ConMatrix.ColPartition);
                var compConMatrix2 = MsrMatrix.LoadFromFile("conMtx2.bin", csMPI.Raw._COMM.WORLD, ConMatrix2.RowPartitioning, ConMatrix2.ColPartition);


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
