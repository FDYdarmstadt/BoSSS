using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.XNSECommon;
using System.Diagnostics;
using MPI.Wrappers;

namespace AdvancedSolverTests {
    /// <summary>
    /// Test control parameter
    /// </summary>
    public enum XDGusage
    {
        none,
        //mixed types do not work with ChangeOfBasis,
        //Exception: Basis types have to be the same
        mixed1,
        mixed2,
        all
    }

    public enum MatrixShape
    {
        laplace, // matrix of a laplacian operator
        full, // a matrix without coupling of species and variables
        full_var, // a matrix without coupling of species
        full_spec, // a matrix without coupling of variables
        full_var_spec, // a matrix
        diagonal, // a matrix without coupling of cells, species and variables
        diagonal_var, // a matrix without coupling of cells and species
        diagonal_var_spec, // a matrix without coupling of cells
        diagonal_spec, // a matrix without coupling of cells and variables
    }

    public class SubBlockTestSolver2Var : Application {

        internal XDGusage m_UseXdg;
        internal MatrixShape m_Mshape;
        internal int m_DGorder = 1;

        protected override IGrid CreateOrLoadGrid() {
            //base.Control.ImmediatePlotPeriod = 1;
            DeleteOldPlotFiles();

            var comm = csMPI.Raw._COMM.WORLD;
            int size;
            csMPI.Raw.Comm_Size(comm,out size);

            
            double xMax = 1.0, yMax = 1.0;
            double xMin = -1.0, yMin = -1.0;
            Func<double[], int> MakeMyPartioning = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                int separation = 2;
                double xspan = (xMax - xMin) / separation;
                double yspan = (yMax - yMin) / separation;
                int rank = int.MaxValue;
                int icore = 0;
                for (int i = 0; i < separation; i++) {
                    for (int j = 0; j < separation; j++) {
                        bool xtrue = x <= xspan * (i + 1) + xMin;
                        bool ytrue = y <= yspan * (j + 1) + yMin;
                        if (xtrue && ytrue) {
                            rank = icore;
                            return rank;
                        }
                        icore++;
                    }
                }

                return rank;
            };

            m_grid = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, 1, m_Res + 1), GenericBlas.Linspace(-1, 1, m_Res + 1));
            if (size == 4) {
                //base.Control.GridPartType = BoSSS.Foundation.Grid.GridPartType.Hilbert;
                base.Control.GridPartType = BoSSS.Foundation.Grid.GridPartType.Predefined;
                base.Control.GridPartOptions = "hallo";
                ((Grid2D)m_grid).AddPredefinedPartitioning("hallo", MakeMyPartioning);
            } else {
                if (size != 1)
                    throw new NotSupportedException("supported MPIsize is 1 or 4");
                base.Control.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            }
            return m_grid;
        }



        LevelSet Phi;

        DGField u1;
        DGField u2;
        DGField MPIrank;

        IGrid m_grid;

        public int m_Res;

        public AggregationGridData[] MgSeq;

        MultigridMapping MG_Mapping;

        public MultigridMapping GetMapping {
            get { return MG_Mapping; }
        }

        SinglePhaseField Amarker;
        SinglePhaseField Bmarker;

        protected override void CreateFields() {
            Phi = new LevelSet(new Basis(this.GridData, 2), "Phi");

            LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData)this.GridData, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 1, new string[] { "A", "B" }, Phi);

            if (m_DGorder < 1)
                throw new ArgumentException();

            if (m_UseXdg == XDGusage.all) {
                u1 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder), "u1");
                u2 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder - 1), "u2");
            } else if (m_UseXdg == XDGusage.none) {
                u1 = new SinglePhaseField(new Basis(this.GridData, m_DGorder), "u1");
                u2 = new SinglePhaseField(new Basis(this.GridData, m_DGorder - 1), "u2");
            } else if (m_UseXdg == XDGusage.mixed1) {
                u1 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder), "u1");
                u2 = new SinglePhaseField(new Basis(this.GridData, m_DGorder - 1), "u2");
            } else if (m_UseXdg == XDGusage.mixed2) {
                u1 = new SinglePhaseField(new Basis(this.GridData, m_DGorder), "u1");
                u2 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder - 1), "u2");
            } else {
                throw new NotImplementedException();
            }

            Amarker = new SinglePhaseField(new Basis(this.GridData, 0), "Amarker");
            Bmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Bmarker");

            MPIrank = new SinglePhaseField(new Basis(this.GridData, 0), "MPIRank");
            MPIrank.AccConstant(this.MPIRank);

        }

        /// <summary>
        /// Cell Agglomeration threshold
        /// </summary>
        internal double THRESHOLD = 0.01;
        //internal double THRESHOLD = 0;

        void LsUpdate(double t) {
            double offset = t;
            Phi.ProjectField((x, y) => -(x - offset).Pow2() - y.Pow2() + (0.3).Pow2());
            LsTrk.UpdateTracker(t);
        }


        protected override void SetInitial(double t) {
            this.LsUpdate(t);

            u1.ProjectField((x, y) => x);
            u2.ProjectField((x, y) => x);

        }

        XDifferentialOperatorMk2 Op;
        int m_quadOrder;

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            m_quadOrder = u1.Basis.Degree * 2;

            
            Setup();

            //XLaplaceBCs xLaplaceBCs = new XLaplaceBCs();
            //xLaplaceBCs.g_Diri = ((CommonParamsBnd inp) => 0.0);
            //xLaplaceBCs.IsDirichlet = (inp => true);
            //double penalty_base = (m_DGorder + 1) * (m_DGorder + 2) / 2;
            //var lengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)GridData).Cells.PenaltyLengthScales;
            //var lap = new XLaplace_Bulk(this.LsTrk, 2.0 * penalty_base, "u1", xLaplaceBCs, 1.0, 1, 1000, lengthScales, XLaplace_Interface.Mode.SIP);

            Op = new XDifferentialOperatorMk2(2, 0, 2, (A, B, c) => m_quadOrder, LsTrk.SpeciesNames, "u1", "u2", "c1", "c2");
            //Op = new XSpatialOperatorMk2(1, 0, 1, (A, B, c) => m_quadOrder, LsTrk.SpeciesIdS.ToArray(), "u1","c1");

            switch (m_Mshape) {
                case MatrixShape.laplace:
                    int p = u1.Basis.Degree;
                    int D=this.GridData.SpatialDimension;
                    double penalty_base = (p + 1) * (p + D) / D;
                    double MU_A = 1;
                    double MU_B = 10;
                    
                    Op.EquationComponents["c1"].Add(new XLaplace_Bulk(MU_A, MU_B, penalty_base * 2, "u1"));      // Bulk form
                    Op.EquationComponents["c1"].Add(new XLaplace_Interface( MU_A, MU_B, penalty_base * 2, "u1"));   // coupling form
                    Op.EquationComponents["c1"].Add(new XLaplace_Bulk(MU_A, MU_B, penalty_base * 2, "u2"));      // Bulk form
                    Op.EquationComponents["c1"].Add(new XLaplace_Interface( MU_A, MU_B, penalty_base * 2, "u2"));   // coupling form
                    
                    Op.EquationComponents["c2"].Add(new XLaplace_Bulk(MU_A, MU_B, penalty_base * 2, "u1"));      // Bulk form
                    Op.EquationComponents["c2"].Add(new XLaplace_Interface( MU_A, MU_B, penalty_base * 2, "u1"));   // coupling form
                    Op.EquationComponents["c2"].Add(new XLaplace_Bulk(MU_A, MU_B, penalty_base * 2, "u2"));      // Bulk form
                    Op.EquationComponents["c2"].Add(new XLaplace_Interface( MU_A, MU_B, penalty_base * 2, "u2"));   // coupling form
                    
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 11)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new SourceTest("u2", 11)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u1", 11)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", 11)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.full:
                    //tested: shape valid for testing
                    Op.EquationComponents["c1"].Add(new DxFlux("u1", 3)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new DxFlux("u2", 4)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.full_var:
                    //tested: shape valid for testing
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 1)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new SourceTest("u2", 2)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new DxFlux("u1", 3)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new DxFlux("u2", 4)); // Flux in Bulk Phase;

                    Op.EquationComponents["c2"].Add(new SourceTest("u1", -5)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", -6)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new DxFlux("u1", -7)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new DxFlux("u2", -8)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.full_spec:
                    //tested: no spec coupling in the secondary diagonals (obviously)
                    Op.EquationComponents["c1"].Add(new LevSetFlx( "u1", -1));
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 1)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new DxFlux("u1", 10)); // Flux in Bulk Phase;

                    Op.EquationComponents["c2"].Add(new LevSetFlx( "u2", -4));
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", -2)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new DxFlux("u2", -11)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.full_var_spec:
                    //tested: no spec coupling in the secondary diagonals (obviously)
                    Op.EquationComponents["c1"].Add(new LevSetFlx( "u1", -1));
                    Op.EquationComponents["c1"].Add(new LevSetFlx( "u2", -1));
                    Op.EquationComponents["c1"].Add(new DxFlux("u1", 3)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new DxFlux("u2", 4)); // Flux in Bulk Phase;

                    Op.EquationComponents["c2"].Add(new LevSetFlx( "u1", -1));
                    Op.EquationComponents["c2"].Add(new LevSetFlx( "u2", -1));
                    Op.EquationComponents["c2"].Add(new DxFlux("u1", 3)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new DxFlux("u2", 4)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.diagonal_var_spec:
                    // block diagonal matrix (ignore cell coupling) 
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 1)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u1", -2)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new SourceTest("u2", 3)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", -4)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new LevSetFlx( "u1", -2));
                    Op.EquationComponents["c2"].Add(new LevSetFlx( "u2", -4));
                    break;
                case MatrixShape.diagonal_spec:
                    // block diagonal matrix (ignore cell and variable coupling) 
                    Op.EquationComponents["c1"].Add(new LevSetFlx( "u1", -2));
                    Op.EquationComponents["c2"].Add(new LevSetFlx( "u2", -4));
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 1)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", -4)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.diagonal:
                    // sparse matrix (ignore cell and variable coupling)
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 1)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", -4)); // Flux in Bulk Phase;
                    break;
                case MatrixShape.diagonal_var:
                    // sparse matrix (ignore cell and variable coupling)
                    Op.EquationComponents["c1"].Add(new SourceTest("u1", 1)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u1", -2)); // Flux in Bulk Phase;
                    Op.EquationComponents["c1"].Add(new SourceTest("u2", 3)); // Flux in Bulk Phase;
                    Op.EquationComponents["c2"].Add(new SourceTest("u2", -4)); // Flux in Bulk Phase;
                    break;
            }




            //Op.EquationComponents["c1"].Add(new LevSetFlx(this.LsTrk, "u1", -1));
            //Op.EquationComponents["c1"].Add(new DxFlux("u2", -1)); // Flux in Bulk Phase;
            //Op.EquationComponents["c2"].Add(new DxFlux("u1", 2)); // Flux in Bulk Phase;
            //Op.EquationComponents["c2"].Add(new DxFlux("u2", -2)); // Flux in Bulk Phase;

            //Op.EquationComponents["c1"].Add(new DxFlux("u1", -3.0)); // Flux in Bulk Phase;
            //Op.EquationComponents["c1"].Add(new LevSetFlx(this.LsTrk, "u1", -3.0));

            //Op.EquationComponents["c2"].Add(new DxFlux("u1", +3.0)); // Flux in Bulk Phase;
            //Op.EquationComponents["c2"].Add(new LevSetFlx(this.LsTrk, "u1", +3.0));
            //Op.EquationComponents["c2"].Add(new DxFlux("u2", 77.7)); // Flux in Bulk Phase;
            //Op.EquationComponents["c2"].Add(new LevSetFlx(this.LsTrk, "u2", 77.7));

            

            Op.Commit();

            Basis maxB = map.BasisS.ElementAtMax(bss => bss.Degree);
            //massFact = new MassMatrixFactory(maxB, agg);
            massFact = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), m_quadOrder).MassMatrixFactory;


            
        }

        MassMatrixFactory massFact;
        AggregationGridBasis[][] XAggB;
        UnsetteledCoordinateMapping map;

        /// <summary>
        /// Operator matrix, i.e. test data for this test.
        /// </summary>
        public BlockMsrMatrix OperatorMatrix;

        public MultigridOperator MGOp;

        /// <summary>
        /// Operator matrix, i.e. test data for this test.
        /// </summary>
        public MsrMatrix AltOperatorMatrix;

        MultigridOperator.ChangeOfBasisConfig[][] OpConfig {
            get {
                return new MultigridOperator.ChangeOfBasisConfig[][] {
                    new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] { 0,1 }, mode = MultigridOperator.Mode.Eye, DegreeS = new int[] { u1.Basis.Degree , u2.Basis.Degree } }
                    }
                };
            }
        }

        private void Setup() {
            MgSeq = CoarseningAlgorithms.CreateSequence(m_grid.iGridData);

            int p = m_DGorder;

            var uMapping = new UnsetteledCoordinateMapping(u1.Basis, u2.Basis);
            //var uMapping = new UnsetteledCoordinateMapping(u1.Basis);
            //var uMapping = new UnsetteledCoordinateMapping(new Basis(m_grid.iGridData, p));

            XAggB = AggregationGridBasis.CreateSequence(MgSeq, uMapping.BasisS);
            var bla = LsTrk.SpeciesIdS.ToArray();
            var agg = LsTrk.GetAgglomerator(bla, m_quadOrder, THRESHOLD, AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);
            XAggB.UpdateXdgAggregationBasis(agg);
            var VarDegrees=uMapping.BasisS.Count.ForLoop(i=> uMapping.BasisS[i].Degree);
            MG_Mapping = new MultigridMapping(uMapping, XAggB[0], VarDegrees);
            map = uMapping;
        }

        private double[] GetRHS(double[] OpAffine, BlockMsrMatrix M) {
            List<int> Rows2Keep = new List<int>();
            for (int iRow=0; iRow < M.RowPartitioning.LocalLength; iRow++) {
                long Row = iRow + M.RowPartitioning.i0;
                Debug.Assert(M.RowPartitioning.IsInLocalRange(Row));
                if (M.GetNoOfNonZerosPerRow(Row) != 0) {
                    Rows2Keep.Add(iRow);
                }
            }

            var rArr=Rows2Keep.ToArray();

            List<double> RHS = new List<double>();

            for(int i = 0; i < rArr.Length; i++) {
                RHS.Add(OpAffine[rArr[i]]);
            }
            return RHS.ToArray();
        }

        public double[] someVec;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            LsUpdate(phystime);

            // operator-matrix assemblieren
            OperatorMatrix = new BlockMsrMatrix(MG_Mapping.ProblemMapping);
            AltOperatorMatrix = new MsrMatrix(MG_Mapping.ProblemMapping);
            double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];
            MultiphaseCellAgglomerator Agg;

            Agg = LsTrk.GetAgglomerator(this.LsTrk.SpeciesIdS.ToArray(), m_quadOrder, __AgglomerationTreshold: this.THRESHOLD);

            XDifferentialOperatorMk2.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(base.LsTrk, MG_Mapping.ProblemMapping, null, MG_Mapping.ProblemMapping);
            mtxBuilder.time = 0.0;
            mtxBuilder.CellLengthScales.AddRange(Agg.CellLengthScales);
            mtxBuilder.ComputeMatrix(OperatorMatrix, Affine);
            Agg.ManipulateMatrixAndRHS(OperatorMatrix, Affine, MG_Mapping.ProblemMapping, MG_Mapping.ProblemMapping);

            Agg.PrintInfo(Console.Out);

            var MamaAgg = this.massFact.GetMassMatrix(map, false);
            Agg.ManipulateMatrixAndRHS(MamaAgg, default(double[]), map, map);

            MGOp = new MultigridOperator(XAggB, map,
                    OperatorMatrix,
                    MamaAgg,
                    OpConfig, this.Op);
            Debug.Assert(MGOp.OperatorMatrix != null);
            Debug.Assert(MGOp.Mapping != null);

            someVec = GetRHS(Affine, OperatorMatrix);

            mtxBuilder.ComputeMatrix(AltOperatorMatrix, Affine);
            Agg.ManipulateMatrixAndRHS(AltOperatorMatrix, Affine, MG_Mapping.ProblemMapping, MG_Mapping.ProblemMapping);


            //LsTrk.GetSpeciesName(((XdgAggregationBasis)MGOp.Mapping.AggBasis[0]).UsedSpecies[1]);
            //LsTrk.GetSpeciesName(((XdgAggregationBasis)MGOp.Mapping.AggBasis[0]).UsedSpecies[0]);

            long nnz = this.OperatorMatrix.GetTotalNoOfNonZeros();
            Console.WriteLine("Number of non-zeros in matrix: " + nnz);

            long nnz2 = this.AltOperatorMatrix.GetTotalNoOfNonZeros();
            Assert.IsTrue(nnz == nnz2, "Number of non-zeros in matrix different for " + OperatorMatrix.GetType() + " and " + AltOperatorMatrix.GetType());
            Console.WriteLine("Number of non-zeros in matrix (reference): " + nnz2);

            MsrMatrix Comp = AltOperatorMatrix.CloneAs();
            Comp.Acc(-1.0, OperatorMatrix);
            double CompErr = Comp.InfNorm();
            double Denom = Math.Max(AltOperatorMatrix.InfNorm(), OperatorMatrix.InfNorm());
            double CompErrRel = Denom > Math.Sqrt(double.Epsilon) ? CompErr / Denom : CompErr;
            Console.WriteLine("Comparison: " + CompErrRel);

            Assert.LessOrEqual(CompErrRel, 1.0e-7, "Huge difference between MsrMatrix and BlockMsrMatrix.");

            base.TerminationKey = true;
            return 0.0;
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "MatrixTest." + timestepNo;
            //Tecplot.PlotFields(new DGField[] { u1, u2, Phi, Amarker, Bmarker, MPIrank }, filename, 0, superSampling);
            Tecplot.PlotFields(new DGField[] { u1, u2, Phi, Amarker, Bmarker, MPIrank }, filename, 0, superSampling);
        }
    }
}
