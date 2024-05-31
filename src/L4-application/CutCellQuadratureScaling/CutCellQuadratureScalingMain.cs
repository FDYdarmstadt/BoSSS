using System;

using ilPSP.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LoadBalancing;
using ilPSP;
using BoSSS.Solution.Tecplot;
using System.Diagnostics;

namespace CutCellQuadratureScaling {


    class CutCellQuadratureScalingMain {

        public static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI(args);

            using(var Ref = new TestSetupSingleLevset2D(1.0, 8)) {
                Ref.Init();
                Ref.RunSolverMode();

                using (var Test = new TestSetupSingleLevset2D(0.5, 8)) {
                    Test.Init();
                    Test.RunSolverMode();

                    Test.CompareTo(Ref);
                }

            }


            BoSSS.Solution.Application.FinalizeMPI();
        }


    }



    class TestSetupSingleLevset2D : Application {

        public TestSetupSingleLevset2D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2) {
            this.MeshScaling = meshScaling;
            this.CutCellQuadratureOrder = cutCellQuadratureOrder;
        }


        public double MeshScaling = 1.0;
        public int CutCellQuadratureOrder = 2;

        protected override IGrid CreateOrLoadGrid() {
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            xNodes.ScaleV(MeshScaling);
            yNodes.ScaleV(MeshScaling);

            GridCommons grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
            return grd;
        }

        LevelSet Phi1;
        //LevelSet Phi2;

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            Phi1 = new LevelSet(new Basis(this.GridData, 2), "Phi1");
            //Phi2 = new LevelSet(new Basis(this.GridData, 2), "Phi2");


            double R = 4.0*this.MeshScaling;
            double Phi1_ana(double x, double y) {
                return 1.0 - (x / R).Pow2() - (y / R).Pow2();
            }
            Debug.Assert(Phi1_ana(R, 0).Abs() <= 1.0e-8);

            Phi1.ProjectField(Phi1_ana);
            //Phi2.ProjectField(((x, y) => y));

            LsTrk = new LevelSetTracker(this.GridData as GridData, XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new [] { "A", "B" }, Phi1);
            LsTrk.UpdateTracker(0.0);
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { Phi1 }, "CutCellQuadratureScaling", timestepNo.MajorNumber, superSampling);
        }

        protected CutCellMetrics latestCCM;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            latestCCM = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS, CutCellQuadratureOrder).CutCellMetrics;
            base.TerminationKey = true;
            return -1;
        }


        public void CompareTo(TestSetupSingleLevset2D othr) {
            foreach (string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var Area_this = this.latestCCM.InterfaceArea[SpcId_this];
                var Area_othr = othr.latestCCM.InterfaceArea[SpcId_othr];
                var Area_err = Area_othr*this.MeshScaling - Area_this*othr.MeshScaling;

                double absErr = Area_err.L2Norm();
                double relErr = absErr / (Area_this.L2Norm() + Area_othr.L2Norm());

                Console.WriteLine($"Level Set Surface, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"Level Set Surface, species {Species} relative error : {relErr:g7}");
            }
        }
    }


    class TestSetupSingleLevset3D : TestSetupSingleLevset2D {


        protected override IGrid CreateOrLoadGrid() {
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] zNodes = GenericBlas.Linspace(-7, +7, 8);
            xNodes.ScaleV(MeshScaling);
            xNodes.ScaleV(MeshScaling);
            xNodes.ScaleV(MeshScaling);

            GridCommons grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
            return grd;
        }

    }

}
