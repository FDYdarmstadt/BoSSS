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
using NUnit.Framework;


namespace BoSSS.Application.CutCellQuadratureScaling {

    /// <summary>
    /// Stumb class for ad-hoc execution of tests
    /// </summary>
    static class CutCellQuadratureScalingMain {

        static CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Algoim;

        public static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI(args);

            using(var Ref = new TestSetupSingleLevset2D(1.0, 8, quadratureType)) {
                Ref.Init();
                Ref.RunSolverMode();

                using (var Test = new TestSetupSingleLevset2D(0.5, 8, quadratureType)) {
                    Test.Init();
                    Test.RunSolverMode();

                    Test.CompareSurfaceTo(Ref);
                    Test.CompareVolumeTo(Ref);
                    Test.CompareEdgeAreaTo(Ref);
                    Test.CompareCutLineTo(Ref);
                }

            }


            BoSSS.Solution.Application.FinalizeMPI();
        }
    }


    /// <summary>
    /// Tests if the 
    /// <see cref="BoSSS.Foundation.Quadrature.ICompositeQuadRule{TQuadRule}.IntegrationMetric"/>
    /// works correctly.
    /// 
    /// Therefore, a reference-integral-evaluation has to be created with a scaling (<see cref="MeshScaling"/>) of 1.
    /// Then, it needs to be compared to a different integral-evaluation has to be created with a non-unit scaling.
    /// This should be done using the methods
    /// - <see cref="TestSetupBase.CompareVolumeTo"/>
    /// - <see cref="TestSetupBase.CompareEdgeAreaTo"/>
    /// - <see cref="TestSetupBase.CompareSurfaceTo"/>
    /// - <see cref="TestSetupBase.CompareCutLineTo"/>
    /// </summary>
    abstract class TestSetupBase : BoSSS.Solution.Application {

        public TestSetupBase(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye) {
            this.MeshScaling = meshScaling;
            this.CutCellQuadratureOrder = cutCellQuadratureOrder;
            this.QuadratureType = quadratureType;
        }


        public readonly double MeshScaling = 1.0;
        public readonly int CutCellQuadratureOrder = 2;
        public readonly CutCellQuadratureMethod QuadratureType = CutCellQuadratureMethod.Saye;

        
        protected CutCellMetrics latestCCM;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            latestCCM = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS, CutCellQuadratureOrder).CutCellMetrics;
            base.TerminationKey = true;
            return -1;
        }


        public void CompareSurfaceTo(TestSetupBase othr) {
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

                Assert.Less(relErr, 1.0e-10, $"relative surface error above threshold for species {Species}" );
            }
        }

        /// <summary>
        /// verifies quadratic/kubic (<see cref="MeshScaling"/> to the power of the spatial dimension) scaling of volume integrals in 2D/3D;
        /// </summary>
        public void CompareVolumeTo(TestSetupBase othr) {
            double D = this.Grid.SpatialDimension;
            
            foreach (string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);



                var Vol_this = this.latestCCM.CutCellVolumes[SpcId_this];
                var Vol_othr = othr.latestCCM.CutCellVolumes[SpcId_othr];
                // Areas should scale with MeshScaling^SpatialDimension
                var Area_err = Vol_othr*(this.MeshScaling.Pow(D)) - Vol_this*(othr.MeshScaling.Pow(D));

                double absErr = Area_err.L2Norm();
                double relErr = absErr / (Vol_this.L2Norm() + Vol_othr.L2Norm());

                Console.WriteLine($"Cut Cell Volume, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"Cut Cell Volume, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, 1.0e-10, $"relative volume error above threshold for species {Species}");
            }
        }


        /// <summary>
        /// verifies quadratic/linear (<see cref="MeshScaling"/> to the power of the spatial dimension - 1) scaling of cell edge area integrals in 2D/3D;
        /// </summary>
        public void CompareEdgeAreaTo(TestSetupBase othr) {
            double D = this.Grid.SpatialDimension;
            
            foreach (string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var Area_this = this.latestCCM.CutEdgeAreas[SpcId_this];
                var Area_othr = othr.latestCCM.CutEdgeAreas[SpcId_othr];
                var Area_err = Area_othr*(this.MeshScaling.Pow(D - 1)) - Area_this*(othr.MeshScaling.Pow(D - 1));

                double absErr = Area_err.L2Norm();
                double relErr = absErr / (Area_this.L2Norm() + Area_othr.L2Norm());

                Console.WriteLine($"Cut Edge area, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"Cut Edge area, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, 1.0e-10, $"relative edge area error above threshold for species {Species}");
            }
        }


        public void CompareCutLineTo(TestSetupBase othr) {
            double D = this.Grid.SpatialDimension;

            foreach(string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var Area_this = this.latestCCM.CutLineLength[SpcId_this];
                var Area_othr = othr.latestCCM.CutLineLength[SpcId_othr];
                //
                var Area_err = Area_othr * (this.MeshScaling.Pow(D - 2)) - Area_this * (othr.MeshScaling.Pow(D - 2));

                double absErr = Area_err.L2Norm();
                double relErr = absErr / (Area_this.L2Norm() + Area_othr.L2Norm());

                Console.WriteLine($"Cut line measure, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"Cut line measure, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, 1.0e-10, $"relative edge area error above threshold for species {Species}");
            }
        }
    }


    abstract class TestSetupSingleLevSetBase : TestSetupBase {

        public TestSetupSingleLevSetBase(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType) { }



        LevelSet Phi1;
        //LevelSet Phi2;

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            Phi1 = new LevelSet(new Basis(this.GridData, 2), "Phi1");
            //Phi2 = new LevelSet(new Basis(this.GridData, 2), "Phi2");


            double R = 4.0 * this.MeshScaling;
            double Phi1_ana(double[] X) {
                double x = X[0], y = X[1];
                return 1.0 - (x / R).Pow2() - (y / R).Pow2();
            }
            var X0 = new double[this.Grid.SpatialDimension]; X0[0] = R;
            Debug.Assert(Phi1_ana(X0).Abs() <= 1.0e-8);

            Phi1.ProjectField(Phi1_ana);
            //Phi2.ProjectField(((x, y) => y));

            LsTrk = new LevelSetTracker(this.GridData as GridData, this.QuadratureType, 1, new[] { "A", "B" }, Phi1);
            LsTrk.UpdateTracker(0.0);
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { Phi1 }, "CutCellQuadratureScaling", timestepNo.MajorNumber, superSampling);
        }


    }

    class TestSetupSingleLevset2D : TestSetupSingleLevSetBase {

        public TestSetupSingleLevset2D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType) { }



        protected override IGrid CreateOrLoadGrid() {
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            xNodes.ScaleV(MeshScaling);
            yNodes.ScaleV(MeshScaling);

            GridCommons grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
            return grd;
        }


    }

    class TestSetupSingleLevset3D : TestSetupSingleLevSetBase {

        public TestSetupSingleLevset3D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType) { }


        protected override IGrid CreateOrLoadGrid() {
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] zNodes = GenericBlas.Linspace(-3, +3, 4);
            xNodes.ScaleV(MeshScaling);
            xNodes.ScaleV(MeshScaling);
            xNodes.ScaleV(MeshScaling);

            GridCommons grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
            return grd;
        }

    }

}
