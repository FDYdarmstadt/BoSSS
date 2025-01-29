using System;
using System.Collections.Generic;
using System.Linq;

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
using BoSSS.Foundation.Grid.RefElements;


namespace BoSSS.Application.CutCellQuadratureScaling {

    /// <summary>
    /// Stumb class for ad-hoc execution of tests
    /// </summary>
    static class CutCellQuadratureScalingMain {

        static CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Classic;
        static int order = 6;

        public static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI(args);

            //BoSSS.Application.CutCellQuadratureScaling.AllTests.OneLevelSet_2D(order, quadratureType);
            BoSSS.Application.CutCellQuadratureScaling.AllTests.OneLevelSet_3D(order, quadratureType);



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
        
        internal CutCellMetrics latestCCM;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            latestCCM = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS, CutCellQuadratureOrder).CutCellMetrics;
            base.TerminationKey = true;
            return -1;
        }

        /// <summary>
        /// threshold for all scaling comparisons (scald vs. un-scaled mesh)
        /// </summary>
        virtual protected double threshold_scaling {
            get {
                return 1.0e-10;
            }
        }


        /// <summary>
        /// Total volume for each species, for <see cref="MeshScaling"/> equal to 1
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        protected Dictionary<string, double> TotalIntegral_Volume4Species  = new Dictionary<string, double>();

        virtual protected double threshold_totVolume {
            get {
                return 1.0e-10;
            }
        }


        /// <summary>
        /// Total surface for each species, for <see cref="MeshScaling"/> equal to 1
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        protected Dictionary<string, double> TotalIntegral_SurfaceSpecies = new Dictionary<string, double>();

        /// <summary>
        /// Total cut line for each species, for <see cref="MeshScaling"/> equal to 1
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        protected Dictionary<string, double> TotalIntegral_CutLine = new Dictionary<string, double>();


        virtual protected double threshold_totSurface {
            get {
                return 1.0e-10;
            }
        }

        virtual protected double threshold_totCutline {
            get {
                return 1.0e-10;
            }
        }

        /// <summary>
        /// verifies that the total volume (cut-cell-volume summed over all cells) is correct
        /// </summary>
        public void CompareTotalVolume() {
            double D = this.Grid.SpatialDimension;
            CompareTotal("Cut Cell volume", this.TotalIntegral_Volume4Species, spcId => this.latestCCM.CutCellVolumes[spcId].Sum(), threshold_totVolume, this.MeshScaling.Pow(D));
        }

        /// <summary>
        /// verifies that the total volume (cut-cell-volume summed over all cells) is correct
        /// </summary>
        public void CompareTotalSurface() {
            double D = this.Grid.SpatialDimension;
            CompareTotal("Cut Cell surface", this.TotalIntegral_SurfaceSpecies, spcId => this.latestCCM.InterfaceArea[spcId].Sum(), threshold_totSurface, this.MeshScaling.Pow(D - 1));
        }


        /// <summary>
        /// verifies that the total volume (cut-cell-volume summed over all cells) is correct
        /// </summary>
        public void CompareTotalCutLine() {
            double D = this.Grid.SpatialDimension;
            CompareTotal("Cut Cell boundary line", this.TotalIntegral_CutLine, spcId => this.latestCCM.CutLineLengthEdge[spcId].Sum(), threshold_totCutline, this.MeshScaling.Pow(D - 2));
        }

        private void CompareTotal(string name, IDictionary<string,double> refDict, Func<SpeciesId, double> getVal, double __threshold, double scaling_D) {
            foreach (var Species in refDict.Keys) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);


                double Vol_this = getVal(SpcId_this);

                double Vol_err = (refDict[Species] - Vol_this/scaling_D).Abs();


                Console.WriteLine($"{name}, species {Species} absolute error : {Vol_err:g7}");
                Console.WriteLine($"      reference value: {Vol_this/scaling_D:g7}   actual result {refDict[Species]:g7}");

                Assert.Less(Vol_err, __threshold, $"{name} error above threshold for species {Species}");
            }
        }

        /// <summary>
        /// verifies linear/quadratic (<see cref="MeshScaling"/> to the power of the spatial dimension) scaling of level-set surface integrals in 2D/3D;
        /// </summary>
        public void CompareSurfaceTo(TestSetupBase othr) {
            double D = this.Grid.SpatialDimension;

            foreach(string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var Area_this = this.latestCCM.InterfaceArea[SpcId_this];
                var Area_othr = othr.latestCCM.InterfaceArea[SpcId_othr];
                var Area_err = Area_othr*(this.MeshScaling.Pow(D - 1)) - Area_this* (othr.MeshScaling.Pow(D - 1));

                double absErr = Area_err.L2Norm();
                double relErr = absErr / (Area_this.L2Norm() + Area_othr.L2Norm());

                Console.WriteLine($"Level Set Surface, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"Level Set Surface, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, threshold_scaling, $"relative surface error above threshold for species {Species}" );
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
                var Vol_err = Vol_othr*(this.MeshScaling.Pow(D)) - Vol_this*(othr.MeshScaling.Pow(D));

                double absErr = Vol_err.L2Norm();
                double relErr = absErr / (Vol_this.L2Norm() + Vol_othr.L2Norm());

                Console.WriteLine($"Cut Cell Volume, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"Cut Cell Volume, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, threshold_scaling, $"relative volume error above threshold for species {Species}");
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

                Assert.Less(relErr, threshold_scaling, $"relative edge area error above threshold for species {Species}");
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

                Assert.Less(relErr, threshold_scaling, $"relative edge area error above threshold for species {Species}");
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

            PlotCurrentState(0.0, new TimestepNumber(0), 4);
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

            double surf2D = 4.0 * 2 * Math.PI; // r*2*pi
            double vol2D_B = 4.0 * 4.0 * Math.PI; // r^2*pi

            base.TotalIntegral_SurfaceSpecies.Add("A", surf2D);
            base.TotalIntegral_SurfaceSpecies.Add("B", surf2D);

            base.TotalIntegral_Volume4Species.Add("A", 14*14 - vol2D_B);
            base.TotalIntegral_Volume4Species.Add("B", vol2D_B);

            base.TotalIntegral_CutLine.Add("A", 16);
            base.TotalIntegral_CutLine.Add("B", 16);


            GridCommons grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
            return grd;
        }
        protected override double threshold_totSurface {
            get {

                if (QuadratureType == CutCellQuadratureMethod.Classic) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 5e-8;
                    }
                }

                if (QuadratureType == CutCellQuadratureMethod.Saye) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 4e-7;
                    }
                }

                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 4e-7;
                    }
                }


                if (QuadratureType == CutCellQuadratureMethod.OneStepGauss) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else if (this.CutCellQuadratureOrder <= 6) {
                        return 4e-5;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 2e-6;
                    }
                }

                if (QuadratureType == CutCellQuadratureMethod.OneStepGaussAndStokes) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 1e-7;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-9;
                    } else {
                        return 4e-11;
                    }
                }


                return base.threshold_totSurface;
            }


        }

        protected override double threshold_totVolume {
            get {

                if (QuadratureType == CutCellQuadratureMethod.Classic) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 3e-8;
                    }
                }

                if(QuadratureType == CutCellQuadratureMethod.Saye) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if(this.CutCellQuadratureOrder <= 9) {
                        return 4e-6;
                    } else {
                        return 1e-8;
                    }
                }

                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if(this.CutCellQuadratureOrder <= 9) {
                        return 4e-6;
                    } else {
                        return 1e-8;
                    }
                }

                if(QuadratureType == CutCellQuadratureMethod.OneStepGauss) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 1e-6;
                    }
                }

                if (QuadratureType == CutCellQuadratureMethod.OneStepGaussAndStokes) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 1e-7;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-9;
                    } else {
                        return 4e-11;
                    }
                }


                return base.threshold_totVolume;
            }


        }
    }


    /// <summary>
    /// The 3D testcase is an extrusion (sic) in z-direction of the 2D testcase, i.e., the circle from the 2D testcase is extruded into a cylinder.
    /// In this fashion, it can be compared to 2D results, by just multiplying the 2D-results with the domain width in z-direction.
    /// Therefore it 
    /// </summary>
    class TestSetupSingleLevset3D : TestSetupSingleLevSetBase {

        public TestSetupSingleLevset3D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType) {
            
        }


        const double zWidht = 6;

        protected override IGrid CreateOrLoadGrid() {
            const int noOfZnodes = 4;
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] zNodes = GenericBlas.Linspace(-3, -3 + zWidht, noOfZnodes);
            xNodes.ScaleV(MeshScaling);
            yNodes.ScaleV(MeshScaling);
            zNodes.ScaleV(MeshScaling);
            GridCommons grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);


            // volume and area of level-set
            // ============================

            double surf3D = 4.0 * 2 * Math.PI*zWidht; // mantle of a cylinder
            double vol3D_B = 4.0 * 4.0 * Math.PI*zWidht; // volume of a cylinder
            base.TotalIntegral_SurfaceSpecies.Add("A", surf3D);
            base.TotalIntegral_SurfaceSpecies.Add("B", surf3D);

            base.TotalIntegral_Volume4Species.Add("A", 14 * 14 * zWidht - vol3D_B); // total mesh volume - cylinder volume
            base.TotalIntegral_Volume4Species.Add("B", vol3D_B);

            double line3D = (4.0 * 2 * Math.PI) * noOfZnodes; // 4 circles in xy-planes
            using(var _2D = new TestSetupSingleLevset2D()) {
                _2D.Init();
                _2D.RunSolverMode();

                line3D += _2D.latestCCM.CutLineLengthEdge[_2D.LsTrk.GetSpeciesId("A")].Sum() * zWidht;
            }
            base.TotalIntegral_CutLine.Add("A", line3D);
            base.TotalIntegral_CutLine.Add("B", line3D);


            return grd;
        }
        protected override double threshold_totSurface {
            get {

                if (QuadratureType == CutCellQuadratureMethod.Classic) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else if (this.CutCellQuadratureOrder <= 8) {
                        return 1e-7;
                    } else {
                        return 5e-8;
                    }
                }

                if (QuadratureType == CutCellQuadratureMethod.Saye) {
                    if (this.CutCellQuadratureOrder <= 3) {
                        return 1e-1;
                    } else if (this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 8e-5;
                    } else if (this.CutCellQuadratureOrder <= 9) {
                        return 3e-6;
                    } else {
                        return 4e-7;
                    }
                }


                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 3) {
                        return 1e-1;
                    } else if(this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 8e-5;
                    } else if(this.CutCellQuadratureOrder <= 9) {
                        return 3e-6;
                    } else {
                        return 4e-7;
                    }
                }


                return base.threshold_totSurface;
            }


        }

        protected override double threshold_totVolume {
            get {

                if (QuadratureType == CutCellQuadratureMethod.Classic) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else  if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 3e-8;
                    }
                }

                if (QuadratureType == CutCellQuadratureMethod.Saye) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 4e-2;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 3e-5;
                    } else if (this.CutCellQuadratureOrder <= 9) {
                        return 4e-6;
                    } else {
                        return 3e-8;
                    }
                }


                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 4e-2;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 3e-5;
                    } else if(this.CutCellQuadratureOrder <= 9) {
                        return 4e-6;
                    } else {
                        return 3e-8;
                    }
                }



                return base.threshold_totVolume;
            }


        }

        protected override double threshold_totCutline {
            get {
                if(QuadratureType == CutCellQuadratureMethod.Classic) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 5e-2;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 5e-5;
                    } else if(this.CutCellQuadratureOrder <= 8) {
                        return 2e-5;
                    } else {
                        return 5e-6;
                    }
                }

                if(QuadratureType == CutCellQuadratureMethod.Saye) {
                    if(this.CutCellQuadratureOrder <= 3) {
                        return 6e-2;
                    } else if(this.CutCellQuadratureOrder <= 4) {
                        return 5e-3;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 5e-5;
                    } else if(this.CutCellQuadratureOrder <= 8) {
                        return 5e-6;
                    } else {
                        return 5e-6;
                    }
                }


                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 5e-3;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 5e-5;
                    } else if(this.CutCellQuadratureOrder <= 8) {
                        return 5e-6;
                    } else {
                        return 5e-6;
                    }
                }



                return base.threshold_totCutline;
            }
        }



        double threshold_2dvs3d {
            get {
                if (QuadratureType == CutCellQuadratureMethod.Classic) {
                    if (this.CutCellQuadratureOrder <= 4)
                        return 1.0e-4;

                    return 1.0e-7;
                }
                

                return 1.0e-10;
            }
        }

        protected override double threshold_scaling {
            get {
                if (QuadratureType == CutCellQuadratureMethod.Saye) {
                    if (this.CutCellQuadratureOrder <= 4)
                        return 1.0e-4;
                    if (this.CutCellQuadratureOrder <= 7)
                        return 1.0e-7;

                    return 1.0e-8;
                }

                return base.threshold_scaling;
            }
        }

        public void CompareSurfaceTo2D(TestSetupSingleLevset2D othr) {
            CompareTo2D("Level Set Surface", othr, test => test.latestCCM.InterfaceArea, -1);
        }
       
        public void CompareVolumeTo2D(TestSetupSingleLevset2D othr) {
            CompareTo2D("Level Set Volume", othr, test => test.latestCCM.CutCellVolumes, 0);
        }

        private void CompareTo2D(string name, TestSetupSingleLevset2D othr, Func<TestSetupBase, IDictionary<SpeciesId, MultidimensionalArray>> propertySelector, int ScalingExponent) {
            double D = othr.Grid.SpatialDimension;
            if(D != 2)
                throw new ApplicationException();

            foreach(string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var totArea_this = propertySelector(this)[SpcId_this].Sum() / (zWidht * MeshScaling);
                var totArea_othr = propertySelector(othr)[SpcId_othr].Sum();

                double absErr = (totArea_othr * (this.MeshScaling.Pow(D - ScalingExponent)) - totArea_this * (othr.MeshScaling.Pow(D - ScalingExponent))).Abs();
                double relErr = absErr / (totArea_this.Abs() + totArea_othr.Abs());

                Console.WriteLine($"{name}, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"{name}, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, threshold_2dvs3d, $"relative {name} error above threshold for species {Species}");
            }
        }
    }
}
