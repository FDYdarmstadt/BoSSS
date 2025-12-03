using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.LoadBalancing;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using static ilPSP.Utils.UnsafeAlgoim;



namespace BoSSS.Application.CutCellQuadratureScaling {

    /// <summary>
    /// Stumb class for ad-hoc execution of tests
    /// </summary>
    static class CutCellQuadratureScalingMain {

        //static CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Classic;
        //static int order = 6;

        public static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI(args);
            //ilPSP.Utils.Algoim.TwoLsIsFucked1();
            //ilPSP.Utils.Algoim.TwoLsIsSuperFucked1();

            BoSSS.Application.CutCellQuadratureScaling.AllTests.OneLevelSet_2D(2, 3, CutCellQuadratureMethod.OneStepGaussAndStokes);
            //BoSSS.Application.CutCellQuadratureScaling.AllTests.TwoLevelSets_2D(2, 8, CutCellQuadratureMethod.Algoim);

            //BoSSS.Application.CutCellQuadratureScaling.AllTests.TwoLevelSets_2D(0, 3, CutCellQuadratureMethod.Algoim);
            //BoSSS.Application.CutCellQuadratureScaling.AllTests.OneLevelSet_3D(0, 3, CutCellQuadratureMethod.Algoim);
            //BoSSS.Application.CutCellQuadratureScaling.AllTests.TwoLevelSet_2Dvs3D(3, CutCellQuadratureMethod.Saye);
            //            BoSSS.Application.CutCellQuadratureScaling.AllTests.TwoLevelSets_3D(3, CutCellQuadratureMethod.Algoim);
            //BoSSS.Application.CutCellQuadratureScaling.AllTests.TwoLevelSets_3D(1, 3, CutCellQuadratureMethod.Algoim);

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
    /// - <see cref="TestSetupBase.CompareTotalVolumeTo"/>
    /// - <see cref="TestSetupBase.CompareTotalEdgeAreaTo"/>
    /// - <see cref="TestSetupBase.CompareTotalSurfaceTo"/>
    /// - <see cref="TestSetupBase.CompareTotalCutLineTo"/>
    /// </summary>
    abstract class TestSetupBase : BoSSS.Solution.Application {

        public TestSetupBase(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0) {
            this.MeshScaling = meshScaling;
            this.CutCellQuadratureOrder = cutCellQuadratureOrder;
            this.QuadratureType = quadratureType;
            this.MeshVariation = meshVariation;
        }

        protected readonly int MeshVariation = 0;
        public readonly double MeshScaling = 1.0;
        public readonly int CutCellQuadratureOrder = 2;

        public readonly CutCellQuadratureMethod QuadratureType = CutCellQuadratureMethod.Saye;
        
        internal CutCellMetrics latestCCM;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //int upgrade = QuadratureType == CutCellQuadratureMethod.Algoim ? 5 : 0;
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
        /// Total volume for each species
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        /// <remarks>
        /// - scales with <see cref="MeshScaling"/> to the power of D
        /// - compares to <see cref="CutCellMetrics.CutCellVolumes"/>
        /// </remarks>
        protected Dictionary<string, double> TotalIntegral_Volume4Species  = new Dictionary<string, double>();

        /// <summary>
        /// Total surface for each species
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        /// <remarks>
        /// - scales with <see cref="MeshScaling"/> to the power of D-1
        /// - compares to <see cref="CutCellMetrics.InterfaceArea"/>
        /// </remarks>
        protected Dictionary<string, double> TotalIntegral_SurfaceSpecies = new Dictionary<string, double>();

        /// <summary>
        /// Total cut line for each species
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        /// <remarks>
        /// - scales with <see cref="MeshScaling"/> to the power of D-2
        /// - compares to <see cref="CutCellMetrics.CutLineLengthEdge"/>
        /// </remarks>
        protected Dictionary<string, double> TotalIntegral_CutLine = new Dictionary<string, double>();


        /// <summary>
        /// Total intersection line for each species
        /// - key: species index
        /// - value: value of the integral over all cells
        /// </summary>
        /// <remarks>
        /// - scales with <see cref="MeshScaling"/> to the power of D-2
        /// - compares to <see cref="CutCellMetrics.IntersectionLength"/>
        /// </remarks>
        protected Dictionary<string, double> TotalIntegral_IntersectionLine = new Dictionary<string, double>();


        virtual protected double threshold_totVolume {
            get {
                return 1.0e-10;
            }
        }

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
        virtual protected double threshold_totIntersectionLine {
            get {
                if(this.LsTrk.NoOfLevelSets <= 0)
                    return 0.0;
                else
                    return 1.0e-10;
            }
        }

        virtual protected double threshold_2dvs3d {
            get {
                throw new NotImplementedException("must be overridden");
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
        /// verifies that the total level-set surface (level-set surface summed over all cells) is correct
        /// </summary>
        public void CompareTotalSurface() {
            double D = this.Grid.SpatialDimension;
            CompareTotal("Cut Cell surface", this.TotalIntegral_SurfaceSpecies, spcId => this.latestCCM.InterfaceArea[spcId].Sum(), threshold_totSurface, this.MeshScaling.Pow(D - 1));
        }

        /// <summary>
        /// verifies that the cut-line-length (cut-line-length summed over all edges) is correct
        /// </summary>
        public void CompareTotalCutLine() {
            double D = this.Grid.SpatialDimension;
            //foreach(var spcName in LsTrk.SpeciesNames) {
            //    var spc = LsTrk.GetSpeciesId(spcName);
            //    (EdgeMask.GetFullMask(GridData, MaskType.Logical)).SaveToTextFile($"CutLine{spcName}.csv", false, (double[] X, int jL, int jG) => this.latestCCM.CutLineLengthEdge[spc][jL]);
            //}
            CompareTotal("Cut Cell boundary line", this.TotalIntegral_CutLine, spcId => this.latestCCM.CutLineLengthEdge[spcId].Sum(), threshold_totCutline, this.MeshScaling.Pow(D - 2));
        }

        /// <summary>
        /// verifies that the intersection-line-length (intersection-line-length summed over all cells) is correct
        /// </summary>
        public void CompareIntersectionLine() {
            double D = this.Grid.SpatialDimension;
            CompareTotal("Intersection line", this.TotalIntegral_IntersectionLine, spcId => this.latestCCM.IntersectionLength[spcId].Sum(), threshold_totIntersectionLine, this.MeshScaling.Pow(D - 2));
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

        #region comparison_to_other
        /// <summary>
        /// verifies linear/quadratic scaling 
        /// (<see cref="MeshScaling"/> to the power of the spatial dimension)
        /// of level-set surface integrals in 2D/3D;
        /// </summary>
        public void CompareTotalSurfaceTo(TestSetupBase othr) {
            CompareTotalTo("Interface Area", othr, (t, spid) => t.latestCCM.InterfaceArea[spid], this.threshold_scaling, -1);
        }

        /// <summary>
        /// verifies quadratic/kubic 
        /// (<see cref="MeshScaling"/> to the power of the spatial dimension)
        /// scaling of volume integrals in 2D/3D;
        /// </summary>
        public void CompareTotalVolumeTo(TestSetupBase othr) {
            CompareTotalTo("Cut Volume", othr, (t, spid) => t.latestCCM.CutCellVolumes[spid], this.threshold_scaling, 0);
        }


        /// <summary>
        /// verifies quadratic/linear 
        /// (<see cref="MeshScaling"/> to the power of the spatial dimension - 1) 
        /// scaling of cell edge area integrals in 2D/3D;
        /// </summary>
        public void CompareTotalEdgeAreaTo(TestSetupBase othr) {
             CompareTotalTo("Cut Edge", othr, (t, spid) => t.latestCCM.CutEdgeAreas[spid], this.threshold_scaling, -1);
        }


        /// <summary>
        /// verifies linear/constant scaling
        ///  (<see cref="MeshScaling"/> to the power of the spatial dimension - 2) 
        /// of surface element edges
        /// </summary>
        public void CompareTotalCutLineTo(TestSetupBase othr) {
            CompareTotalTo("Cut Line", othr, (t, spid) => t.latestCCM.CutLineLength[spid], this.threshold_scaling, -2);
        }

        /// <summary>
        /// verifies linear/constant scaling
        ///  (<see cref="MeshScaling"/> to the power of the spatial dimension - 2) 
        /// of intersection lines
        /// </summary>
        public void CompareTotalIntersectionLineTo(TestSetupBase othr) {
            CompareTotalTo("Intersection Line", othr, (t, spid) => t.latestCCM.IntersectionLength[spid], this.threshold_scaling, -2);
        }

        private void CompareTotalTo(string name, TestSetupBase othr, Func<TestSetupBase, SpeciesId, MultidimensionalArray> getRes, double __threshold, int scaling_D) {
            double D = this.Grid.SpatialDimension;

            foreach(string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var Area_this = getRes(this, SpcId_this);
                var Area_othr = getRes(othr, SpcId_othr);
                //
                var Area_err = Area_othr * (this.MeshScaling.Pow(D + scaling_D)) - Area_this * (othr.MeshScaling.Pow(D + scaling_D));

                double absErr = Area_err.L2Norm();
                double relErr = absErr / (Area_this.L2Norm() + Area_othr.L2Norm());

                Console.WriteLine($"{name} measure, species {Species} absolute error : {absErr:g7}");
                Console.WriteLine($"{name} measure, species {Species} relative error : {relErr:g7}");

                Assert.Less(relErr, threshold_scaling, $"relative {name} error above threshold for species {Species}");
            }
        }
        #endregion


        #region comparison_to_other_perElement
        /// <summary>
        /// verifies linear/quadratic scaling 
        /// (<see cref="MeshScaling"/> to the power of the spatial dimension)
        /// of level-set surface integrals in 2D/3D;
        /// </summary>
        public void CompareElementSurfaceTo(TestSetupBase othr) {
            //this.latestCCM.WriteSurfaceRulesToVtp();
            //othr.latestCCM.WriteSurfaceRulesToVtp();
            CompareElementTo("Interface Area", othr, (t, spid) => t.latestCCM.InterfaceArea[spid], this.threshold_scaling, -1);
        }

        /// <summary>
        /// verifies quadratic/kubic 
        /// (<see cref="MeshScaling"/> to the power of the spatial dimension)
        /// scaling of volume integrals in 2D/3D;
        /// </summary>
        public void CompareElementVolumeTo(TestSetupBase othr) {
            //this.latestCCM.WriteVolumeRulesToVtp();
            //othr.latestCCM.WriteVolumeRulesToVtp();
            CompareElementTo("Cut Volume", othr, (t, spid) => t.latestCCM.CutCellVolumes[spid], this.threshold_scaling, 0);
        }

        /// <summary>
        /// verifies quadratic/linear 
        /// (<see cref="MeshScaling"/> to the power of the spatial dimension - 1) 
        /// scaling of cell edge area integrals in 2D/3D;
        /// </summary>
        public void CompareElementEdgeAreaTo(TestSetupBase othr) {
            CompareElementTo("Cut Edge", othr, (t, spid) => t.latestCCM.CutEdgeAreas[spid], this.threshold_scaling, -1);
        }

        /// <summary>
        /// verifies linear/constant scaling
        ///  (<see cref="MeshScaling"/> to the power of the spatial dimension - 2) 
        /// of surface element edges
        /// </summary>
        public void CompareElementCutLineTo(TestSetupBase othr) {
            CompareElementTo("Cut Line", othr, (t, spid) => t.latestCCM.CutLineLength[spid], this.threshold_scaling, -2);
        }

        /// <summary>
        /// verifies linear/constant scaling
        ///  (<see cref="MeshScaling"/> to the power of the spatial dimension - 2) 
        /// of intersection lines
        /// </summary>
        public void CompareElementIntersectionLineTo(TestSetupBase othr) {
            CompareElementTo("Intersection Line", othr, (t, spid) => t.latestCCM.IntersectionLength[spid], this.threshold_scaling, -2);
        }

        private void CompareElementTo(string name, TestSetupBase othr, Func<TestSetupBase, SpeciesId, MultidimensionalArray> getRes, double __threshold, int scaling_D) {
            double D = this.Grid.SpatialDimension;
            string _name = name.Replace(" ", "_");
            foreach(string Species in this.LsTrk.SpeciesNames) {
                var SpcId_this = this.LsTrk.GetSpeciesId(Species);
                var SpcId_othr = othr.LsTrk.GetSpeciesId(Species);

                var Res_this = getRes(this, SpcId_this).To1DArray();
                var Res_othr = getRes(othr, SpcId_othr).To1DArray();
                //

                double[] ElementError = Res_this.CloneAs();
                ElementError.ScaleV(this.MeshScaling.Pow(D + scaling_D));
                ElementError.AccV(-othr.MeshScaling.Pow(D + scaling_D), Res_othr);

                var l2_err = ElementError.L2Norm();
                Console.WriteLine($"Element-wise {name} error, species {Species}, l2 norm is: " + l2_err);

                int L = ElementError.Length;
                if(L == this.GridData.iLogicalCells.Count) {

                    CellMask.GetFullMask(this.GridData, MaskType.Logical).SaveToTextFile("error_" + _name + "_" + Species + ".csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => ElementError[LogicalItemIndex]);
                } else {
                    EdgeMask.GetFullMask(this.GridData, MaskType.Logical).SaveToTextFile("error_" + _name + "_" + Species + ".csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => ElementError[LogicalItemIndex]);
                }


                //Console.WriteLine($"{name} measure, species {Species} absolute error : {absErr:g7}");
                //Console.WriteLine($"{name} measure, species {Species} relative error : {relErr:g7}");

                //Assert.Less(relErr, threshold_scaling, $"relative {name} error above threshold for species {Species}");
            }
        }        
        #endregion

        #region comparison_2Dvs3D
        public void CompareSurfaceTo2D(TestSetupBase othr) {
            CompareTo2D("Level Set Surface", othr, test => test.latestCCM.InterfaceArea, -1);
        }

        public void CompareVolumeTo2D(TestSetupBase othr) {
            CompareTo2D("Level Set Volume", othr, test => test.latestCCM.CutCellVolumes, 0);
        }

        public void CompareIntersectionLineTo2D(TestSetupBase othr) {
            CompareTo2D("Intersection Line Length", othr, test => test.latestCCM.IntersectionLength, -2);
        }

        private void CompareTo2D(string name, TestSetupBase othr, Func<TestSetupBase, IDictionary<SpeciesId, MultidimensionalArray>> propertySelector, int ScalingExponent) {
            double D = othr.Grid.SpatialDimension;
            if(D != 2)
                throw new ArgumentException("other testcase is supposed to be 2D", nameof(othr));
            if(this.Grid.SpatialDimension != 3)
                throw new NotSupportedException("this testcase is supposed to be 3D");

            var BB = (this.GridData as BoSSS.Foundation.Grid.Classic.GridData).GlobalBoundingBox;
            double zWidht = BB.Max[2] - BB.Min[2];

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
        #endregion


        protected IGrid CreateOrLoadGrid_2D() {
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            switch(MeshVariation) {
                case 0: break;
                case 1: {
                    double dx = xNodes[1] - xNodes[0];
                    xNodes[1] -= dx * 0.3;
                    xNodes[6] -= dx * 0.3;
                    yNodes[3] -= dx * 0.1;
                    break;
                }
                case 2: {
                    ArrayTools.AddToArray(0, ref yNodes);
                    Array.Sort(yNodes);
                    break;
                }
                default:
                throw new NotImplementedException("unknown mesh variation");
            }

            xNodes.ScaleV(MeshScaling);
            yNodes.ScaleV(MeshScaling);


            ReferenceResults_2D();

            GridCommons grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
            return grd;
        }

        protected IGrid CreateOrLoadGrid_3D(double zWidht) {
            const int noOfZnodes = 4;
            double[] xNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] yNodes = GenericBlas.Linspace(-7, +7, 8);
            double[] zNodes = GenericBlas.Linspace(-3, -3 + zWidht, noOfZnodes);
            switch(MeshVariation) {
                case 0: break;
                case 1: {
                    double dx = xNodes[1] - xNodes[0];
                    xNodes[1] -= dx * 0.3;
                    xNodes[6] -= dx * 0.3;
                    yNodes[3] -= dx * 0.1;
                    break;
                }
                case 2: {
                    ArrayTools.AddToArray(0, ref yNodes);
                    Array.Sort(yNodes);
                    break;
                }
                default:
                throw new NotImplementedException("unknown mesh variation");
            }

            xNodes.ScaleV(MeshScaling);
            yNodes.ScaleV(MeshScaling);
            zNodes.ScaleV(MeshScaling);
            GridCommons grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);


            ReferenceResults_3D(zWidht, noOfZnodes);


            return grd;
        }

        protected abstract void ReferenceResults_3D(double zWidht, int noOfZnodes);
        protected abstract void ReferenceResults_2D();


        Dictionary<string, SinglePhaseField> m_SpeciesMarkers;

        void GenerateSpeciesMarkers() {
            m_SpeciesMarkers = new Dictionary<string, SinglePhaseField>();

            // Degree 0 basis for marker fields
            var markerBasis = new Basis(this.GridData, 0);

            foreach(var spc in this.LsTrk.SpeciesNames) {
                var marker = new SinglePhaseField(markerBasis, "SpcMarker-" + spc);
                var spcSubGrid = this.LsTrk.Regions.GetSpeciesMask(spc);
                    marker.AccConstant(1.0, spcSubGrid);

                SpeciesMarkers.Add(spc, marker);
                base.RegisterField(marker, IOListOption.Always);
                
            }
        }

        protected Dictionary<string, SinglePhaseField> SpeciesMarkers {
            get {
                if(m_SpeciesMarkers == null)
                    GenerateSpeciesMarkers();
                return m_SpeciesMarkers;
            }
        }

        internal void DoPlot() {
            PlotCurrentState(0, new TimestepNumber(), 3);
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(this.LsTrk.LevelSets.Select(f => f as DGField).ToArray().Cat(SpeciesMarkers.Values), "CutCellQuadratureScaling", timestepNo.MajorNumber, superSampling);
        }
    }


    abstract class TestSetupSingleLevSetBase : TestSetupBase {

        public TestSetupSingleLevSetBase(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType, meshVariation) { }



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

        protected override void ReferenceResults_2D() {
            double surf2D = 4.0 * 2 * Math.PI; // r*2*pi
            double vol2D_B = 4.0 * 4.0 * Math.PI; // r^2*pi

            TotalIntegral_SurfaceSpecies.Add("A", surf2D);
            TotalIntegral_SurfaceSpecies.Add("B", surf2D);

            TotalIntegral_Volume4Species.Add("A", 14 * 14 - vol2D_B);
            TotalIntegral_Volume4Species.Add("B", vol2D_B);

            if(base.MeshVariation == 2) {
                // one additional mesh node in y-direction at y = 0
                TotalIntegral_CutLine.Add("A", 18);
                TotalIntegral_CutLine.Add("B", 18);
            } else {
                TotalIntegral_CutLine.Add("A", 16);
                TotalIntegral_CutLine.Add("B", 16);
            }

            TotalIntegral_IntersectionLine.Add("A", 0);
            TotalIntegral_IntersectionLine.Add("B", 0);
        }

        protected override void ReferenceResults_3D(double zWidht, int noOfZnodes) {
            // volume and area of level-set
            // ============================

            double surf3D = 4.0 * 2 * Math.PI * zWidht; // mantle of a cylinder
            double vol3D_B = 4.0 * 4.0 * Math.PI * zWidht; // volume of a cylinder
            TotalIntegral_SurfaceSpecies.Add("A", surf3D);
            TotalIntegral_SurfaceSpecies.Add("B", surf3D);

            TotalIntegral_Volume4Species.Add("A", 14 * 14 * zWidht - vol3D_B); // total mesh volume - cylinder volume
            TotalIntegral_Volume4Species.Add("B", vol3D_B);

            double line3D = (4.0 * 2 * Math.PI) * noOfZnodes; // 4 circles in xy-planes
            using(var _2D = new TestSetupSingleLevset2D(meshScaling:1, cutCellQuadratureOrder:this.CutCellQuadratureOrder, quadratureType:this.QuadratureType, meshVariation: MeshVariation)) {
                _2D.Init();
                _2D.RunSolverMode();

                line3D += _2D.latestCCM.CutLineLengthEdge[_2D.LsTrk.GetSpeciesId("A")].Sum() * zWidht;
            }
            TotalIntegral_CutLine.Add("A", line3D);
            TotalIntegral_CutLine.Add("B", line3D);

            TotalIntegral_IntersectionLine.Add("A", 0);
            TotalIntegral_IntersectionLine.Add("B", 0);
        }



    }

    class TestSetupSingleLevset2D : TestSetupSingleLevSetBase {

        public TestSetupSingleLevset2D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType, meshVariation) { }



        protected override IGrid CreateOrLoadGrid() {
            return base.CreateOrLoadGrid_2D();
        }
        protected override double threshold_totSurface {
            get {

                if (QuadratureType == CutCellQuadratureMethod.Classic) {
                    if (this.CutCellQuadratureOrder <= 4) {
                        return 3e-3;
                    } else if (this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else if (this.CutCellQuadratureOrder <= 9) {
                        return 5e-7;
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
                    } else if (this.CutCellQuadratureOrder <= 8) {
                        return 5e-5;
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
                    } else if(this.CutCellQuadratureOrder <= 6) {
                        return 5e-5;
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
    /// </summary>
    class TestSetupSingleLevset3D : TestSetupSingleLevSetBase {

        public TestSetupSingleLevset3D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType, meshVariation) {
           
        }


        const double zWidht = 6;

        protected override IGrid CreateOrLoadGrid() {
            return base.CreateOrLoadGrid_3D(zWidht);
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
                    } else if(this.CutCellQuadratureOrder <= 6) {
                        return 1e-4;
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


        protected override double threshold_2dvs3d {
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


    }




    abstract class TestSetupTwoLevSetsBase : TestSetupBase {

        public TestSetupTwoLevSetsBase(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType, meshVariation) { }



        LevelSet Phi1;
        LevelSet Phi2;

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            Phi1 = new LevelSet(new Basis(this.GridData, 2), "Phi1");
            Phi2 = new LevelSet(new Basis(this.GridData, 2), "Phi2");


            double R = 4.0 * this.MeshScaling;
            double Phi1_ana(double[] X) {
                double x = X[0], y = X[1];
                return 1.0 - (x / R).Pow2() - (y / R).Pow2();
            }
            var X0 = new double[this.Grid.SpatialDimension]; X0[0] = R;
            Debug.Assert(Phi1_ana(X0).Abs() <= 1.0e-8);
            double Phi2_ana(double[] X) {
                double y = X[1] / this.MeshScaling;
                double phi = - y;
                //double phi = y.Pow2() * 0.05 - y;
                return phi;
            }

            Phi1.ProjectField(Phi1_ana);
            Phi2.ProjectField(Phi2_ana);

            string[,] speciesTable = new string[2, 2];
            speciesTable[0, 0] = "A";
            speciesTable[1, 0] = "B";
            speciesTable[0, 1] = "C";
            speciesTable[1, 1] = "C";

            LsTrk = new LevelSetTracker(this.GridData as GridData, this.QuadratureType, 1, speciesTable, Phi1, Phi2);
            LsTrk.UpdateTracker(0.0);

            PlotCurrentState(0.0, new TimestepNumber(0), 4);
        }


        protected override void ReferenceResults_2D() {
            double surf2D = 4.0 * 2 * Math.PI*0.5; // (half-circle)
            double vol2D_B = 4.0 * 4.0 * Math.PI*0.5; // r^2*pi/2

            TotalIntegral_SurfaceSpecies.Add("A", surf2D + 6);
            TotalIntegral_SurfaceSpecies.Add("B", surf2D + 8);
            TotalIntegral_SurfaceSpecies.Add("C", 14);

            TotalIntegral_Volume4Species.Add("A", 14 * 7 - vol2D_B);
            TotalIntegral_Volume4Species.Add("B", vol2D_B);
            TotalIntegral_Volume4Species.Add("C", 14 * 7);

            TotalIntegral_CutLine.Add("A", 12);
            TotalIntegral_CutLine.Add("B", 12);
            TotalIntegral_CutLine.Add("C", 8);

            TotalIntegral_IntersectionLine.Add("A", 2);
            TotalIntegral_IntersectionLine.Add("B", 2);
            TotalIntegral_IntersectionLine.Add("C", 2);

        }

        protected override void ReferenceResults_3D(double zWidht, int noOfZnodes) {
            // volume and area of level-set
            // ============================

            double surf3D = 4.0 * 2 * Math.PI * 0.5 * zWidht; // mantle of a half-cylinder
            double vol3D_B = 4.0 * 4.0 * Math.PI * 0.5 * zWidht; // volume of a half-cylinder
            TotalIntegral_SurfaceSpecies.Add("A", surf3D + 6*zWidht);
            TotalIntegral_SurfaceSpecies.Add("B", surf3D + 8*zWidht);
            TotalIntegral_SurfaceSpecies.Add("C", 14*zWidht);


            TotalIntegral_Volume4Species.Add("A", 14 * 7 * zWidht - vol3D_B); // half mesh volume - cylinder volume
            TotalIntegral_Volume4Species.Add("B", vol3D_B);
            TotalIntegral_Volume4Species.Add("C", 14 * 7 * zWidht); // lower half

            double line3D_A = (4.0 * 2 * Math.PI * 0.5 + 6) * noOfZnodes; // `noOfZnodes` half-circles in xy-planes + 2*3
            double line3D_B = (4.0 * 2 * Math.PI * 0.5 + 8) * noOfZnodes; // `noOfZnodes` half-circles in xy-planes + 8
            double line3D_C = 14 * noOfZnodes;  // `noOfZnodes` lines from left to right
            using(var _2D = new TestSetupTwoLevSets2D(cutCellQuadratureOrder:this.CutCellQuadratureOrder, quadratureType: this.QuadratureType, meshVariation:this.MeshVariation)) {
                _2D.Init();
                _2D.RunSolverMode();

                line3D_A += _2D.latestCCM.CutLineLengthEdge[_2D.LsTrk.GetSpeciesId("A")].Sum() * zWidht;
                line3D_B += _2D.latestCCM.CutLineLengthEdge[_2D.LsTrk.GetSpeciesId("B")].Sum() * zWidht;
                line3D_C += _2D.latestCCM.CutLineLengthEdge[_2D.LsTrk.GetSpeciesId("C")].Sum() * zWidht;
            }
            TotalIntegral_CutLine.Add("A", line3D_A);
            TotalIntegral_CutLine.Add("B", line3D_B);
            TotalIntegral_CutLine.Add("C", line3D_C);


            TotalIntegral_IntersectionLine.Add("A", 2 * zWidht);
            TotalIntegral_IntersectionLine.Add("B", 2 * zWidht);
            TotalIntegral_IntersectionLine.Add("C", 2 * zWidht);
        }
    }


    class TestSetupTwoLevSets2D : TestSetupTwoLevSetsBase {

        public TestSetupTwoLevSets2D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType, meshVariation) { }



        protected override IGrid CreateOrLoadGrid() {
            return base.CreateOrLoadGrid_2D();
        }
        protected override double threshold_totSurface {
            get {

                //if(QuadratureType == CutCellQuadratureMethod.Classic) {
                //    if(this.CutCellQuadratureOrder <= 4) {
                //        return 3e-3;
                //    } else if(this.CutCellQuadratureOrder <= 7) {
                //        return 1e-5;
                //    } else {
                //        return 5e-8;
                //    }
                //}

                if(QuadratureType == CutCellQuadratureMethod.Saye) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 1e-5;
                    } else {
                        return 4e-7;
                    }
                }

                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 4) {
                        return 1e-2;
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 3e-5;
                    } else if(this.CutCellQuadratureOrder <= 8) {
                        return 9e-6;
                    } else if(this.CutCellQuadratureOrder <= 9) {
                        return 4e-6;
                    } else if(this.CutCellQuadratureOrder <= 10) {
                        return 1e-6;
                    } else {
                        return 4e-7;
                    }
                }


                //if(QuadratureType == CutCellQuadratureMethod.OneStepGauss) {
                //    if(this.CutCellQuadratureOrder <= 4) {
                //        return 3e-3;
                //    } else if(this.CutCellQuadratureOrder <= 6) {
                //        return 4e-5;
                //    } else if(this.CutCellQuadratureOrder <= 7) {
                //        return 1e-5;
                //    } else {
                //        return 2e-6;
                //    }
                //}

                //if(QuadratureType == CutCellQuadratureMethod.OneStepGaussAndStokes) {
                //    if(this.CutCellQuadratureOrder <= 4) {
                //        return 1e-7;
                //    } else if(this.CutCellQuadratureOrder <= 7) {
                //        return 1e-9;
                //    } else {
                //        return 4e-11;
                //    }
                //}


                return base.threshold_totSurface;
            }


        }

        protected override double threshold_totVolume {
            get {

                //if(QuadratureType == CutCellQuadratureMethod.Classic) {
                //    if(this.CutCellQuadratureOrder <= 4) {
                //        return 3e-3;
                //    } else if(this.CutCellQuadratureOrder <= 7) {
                //        return 1e-5;
                //    } else {
                //        return 3e-8;
                //    }
                //}

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
                    } else if(this.CutCellQuadratureOrder <= 7) {
                        return 5e-5;
                    } else if(this.CutCellQuadratureOrder <= 9) {
                        return 3e-5;
                    } else if(this.CutCellQuadratureOrder <= 10) {
                        return 3e-6;
                    } else {
                        return 1e-8;
                    }
                }

                //if(QuadratureType == CutCellQuadratureMethod.OneStepGauss) {
                //    if(this.CutCellQuadratureOrder <= 4) {
                //        return 3e-3;
                //    } else if(this.CutCellQuadratureOrder <= 6) {
                //        return 5e-5;
                //    } else if(this.CutCellQuadratureOrder <= 7) {
                //        return 1e-5;
                //    } else {
                //        return 1e-6;
                //    }
                //}

                //if(QuadratureType == CutCellQuadratureMethod.OneStepGaussAndStokes) {
                //    if(this.CutCellQuadratureOrder <= 4) {
                //        return 1e-7;
                //    } else if(this.CutCellQuadratureOrder <= 7) {
                //        return 1e-9;
                //    } else {
                //        return 4e-11;
                //    }
                //}


                return base.threshold_totVolume;
            }


        }

    }


    /// <summary>
    /// The 3D testcase is an extrusion (sic) in z-direction of the 2D testcase, i.e., the circle from the 2D testcase is extruded into a cylinder.
    /// In this fashion, it can be compared to 2D results, by just multiplying the 2D-results with the domain width in z-direction.
    /// </summary>
    class TestSetupTwoLevSets3D : TestSetupTwoLevSetsBase {

        public TestSetupTwoLevSets3D(double meshScaling = 1.0, int cutCellQuadratureOrder = 2, CutCellQuadratureMethod quadratureType = CutCellQuadratureMethod.Saye, int meshVariation = 0)
            : base(meshScaling, cutCellQuadratureOrder, quadratureType, meshVariation) {

        }


        const double zWidht = 6;

        protected override IGrid CreateOrLoadGrid() {
            return base.CreateOrLoadGrid_3D(zWidht);
        }
        protected override double threshold_totSurface {
            get {

                if(QuadratureType == CutCellQuadratureMethod.Saye) {
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

                if(QuadratureType == CutCellQuadratureMethod.Saye) {
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


        protected override double threshold_2dvs3d {
            get {
                //if(QuadratureType == CutCellQuadratureMethod.Classic) {
                //    if(this.CutCellQuadratureOrder <= 4)
                //        return 1.0e-4;
                //    return 1.0e-7;
                //}

                return 1.0e-10;
            }
        }

        protected override double threshold_scaling {
            get {
                if(QuadratureType == CutCellQuadratureMethod.Saye) {
                    if(this.CutCellQuadratureOrder <= 4)
                        return 1.0e-4;
                    if(this.CutCellQuadratureOrder <= 7)
                        return 1.0e-7;

                    return 1.0e-8;
                }

                if(QuadratureType == CutCellQuadratureMethod.Algoim) {
                    if(this.CutCellQuadratureOrder <= 4)
                        return 1.0e-4;
                    if(this.CutCellQuadratureOrder <= 7)
                        return 1.0e-7;

                    return 1.0e-8;
                }

                return base.threshold_scaling;
            }
        }

        
    }




}
