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

namespace BoSSS.Application.Matrix_MPItest {

    /// <summary>
    /// Test control parameter
    /// </summary>
    public enum XDGusage {
        none,
        mixed1,
        mixed2,
        all
    }


    /// <summary>
    /// We are testing several functions of the <see cref="BlockMsrMatrix"/>
    /// in a MPI-parallel setup. In order to have a sufficiently complicated matrix (variable block-lengths, etc.), 
    /// we use XDG to create the test setup.
    /// </summary>
    class Matrix_MPItestMain : BoSSS.Solution.Application {


        internal XDGusage m_UseXdg;
        internal int m_DGorder = 2;

        protected override IGrid CreateOrLoadGrid() {
            base.Control.GridPartType = BoSSS.Foundation.Grid.GridPartType.METIS;
            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-3, 3, 13), GenericBlas.Linspace(-3, 3, 13));
            //var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-3, 3, 8), GenericBlas.Linspace(-3, 3, 2));
            //Console.WriteLine("testcode");

            return grd;
        }

        LevelSet Phi;
                
        DGField u1;
        DGField u2;
        DGField MPIrank;

        CoordinateMapping m_ProblemMapping;

        /// <summary>
        /// Some Coordinate mapping of a system.
        /// </summary>
        public CoordinateMapping ProblemMapping {
            get {
                if (m_ProblemMapping == null)
                    m_ProblemMapping = new CoordinateMapping(u1, u2);
                return m_ProblemMapping;
            }
        }

        SinglePhaseField Amarker;
        SinglePhaseField Bmarker;


        protected override void CreateFields() {
            Phi = new LevelSet(new Basis(this.gridData, 2), "Phi");
   
            LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData)this.gridData, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 1, new string[] { "A", "B" }, Phi);

            if (m_DGorder < 1)
                throw new ArgumentException();

            if (m_UseXdg == XDGusage.all) {
                u1 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder), "u1");
                u2 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder -1), "u2");
            } else if (m_UseXdg == XDGusage.none) {
                u1 = new SinglePhaseField(new Basis(this.gridData, m_DGorder), "u1");
                u2 = new SinglePhaseField(new Basis(this.gridData, m_DGorder - 1), "u2");
            } else if(m_UseXdg == XDGusage.mixed1) {
                u1 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder), "u1");
                u2 = new SinglePhaseField(new Basis(this.gridData, m_DGorder - 1), "u2");
            } else if(m_UseXdg == XDGusage.mixed2) {
                u1 = new SinglePhaseField(new Basis(this.gridData, m_DGorder), "u1");
                u2 = new XDGField(new XDGBasis(this.LsTrk, m_DGorder - 1), "u2");
            } else {
                throw new NotImplementedException();
            }

            Amarker = new SinglePhaseField(new Basis(this.gridData, 0), "Amarker");
            Bmarker = new SinglePhaseField(new Basis(this.gridData, 0), "Bmarker");

            MPIrank = new SinglePhaseField(new Basis(this.gridData, 0), "MPIRank");
            MPIrank.AccConstant(this.MPIRank);

        }
        
        /// <summary>
        /// Cell Agglomeration threshold
        /// </summary>
        internal double THRESHOLD = 0.3;


        void LsUpdate(double t) {
            double offset = t;
            Phi.ProjectField((x, y) => -(x - offset).Pow2() - y.Pow2() + (2.4).Pow2());
            //Phi.ProjectField((x, y) => x);
            //Console.WriteLine("more testcode");
            LsTrk.UpdateTracker();
        }


        protected override void SetInitial() {
            this.LsUpdate(0.0);

            u1.ProjectField((x, y) => x);
            u2.ProjectField((x, y) => x);

        }

        XSpatialOperator Op;
        int m_quadOrder;

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            m_quadOrder = u1.Basis.Degree * 2;

            Op = new XSpatialOperator(2, 0, 2, (A, B, c) => m_quadOrder, "u1", "u2", "c1", "c2");
            
            Op.EquationComponents["c1"].Add(new DxFlux("u1", -3.0)); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetFlx(this.LsTrk, "u1", -3.0));

            Op.EquationComponents["c2"].Add(new DxFlux("u1", +3.0)); // Flux in Bulk Phase;
            Op.EquationComponents["c2"].Add(new LevSetFlx(this.LsTrk, "u1", +3.0));
            Op.EquationComponents["c2"].Add(new DxFlux("u2", 77.7)); // Flux in Bulk Phase;
            Op.EquationComponents["c2"].Add(new LevSetFlx(this.LsTrk, "u2", 77.7));

            Op.Commit();
        }

        /// <summary>
        /// Operator matrix, i.e. test data for this test.
        /// </summary>
        public BlockMsrMatrix OperatorMatrix;

        /// <summary>
        /// Operator matrix, i.e. test data for this test.
        /// </summary>
        public MsrMatrix AltOperatorMatrix;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            
            
            //phystime = 1.8;
            LsUpdate(phystime);
            
            // operator-matrix assemblieren
            OperatorMatrix = new BlockMsrMatrix(ProblemMapping);
            AltOperatorMatrix = new MsrMatrix(ProblemMapping);
            double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];
            MultiphaseCellAgglomerator Agg;

            // Agglomerator setup
            //Agg = new MultiphaseCellAgglomerator(new CutCellMetrics(MomentFittingVariant, m_quadOrder, LsTrk, LsTrk.GetSpeciesId("B")), this.THRESHOLD, false);
            Agg = LsTrk.GetAgglomerator(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, m_quadOrder, __AgglomerationTreshold: this.THRESHOLD);
            Console.WriteLine("Inter-Process agglomeration? " + Agg.GetAgglomerator(LsTrk.GetSpeciesId("B")).AggInfo.InterProcessAgglomeration);
            
            // operator matrix assembly
            Op.ComputeMatrixEx(LsTrk,
                ProblemMapping, null, ProblemMapping,
                OperatorMatrix, Affine, false, 0.0, true,
                Agg.CellLengthScales, null, null,
                LsTrk.SpeciesIdS.ToArray());
            Agg.ManipulateMatrixAndRHS(OperatorMatrix, Affine, this.ProblemMapping, this.ProblemMapping);

            Op.ComputeMatrixEx(LsTrk,
                ProblemMapping, null, ProblemMapping,
                AltOperatorMatrix, Affine, false, 0.0, true,
                Agg.CellLengthScales, null, null,
                LsTrk.SpeciesIdS.ToArray());
            Agg.ManipulateMatrixAndRHS(AltOperatorMatrix, Affine, this.ProblemMapping, this.ProblemMapping);


            int nnz = this.OperatorMatrix.GetTotalNoOfNonZeros();
            Console.WriteLine("Number of non-zeros in matrix: " + nnz);
           
            int nnz2 = this.AltOperatorMatrix.GetTotalNoOfNonZeros();
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
            Tecplot.PlotFields(new DGField[] { u1, u2, Phi, Amarker, Bmarker, MPIrank },  filename, 0, superSampling);
        }
    }
}
