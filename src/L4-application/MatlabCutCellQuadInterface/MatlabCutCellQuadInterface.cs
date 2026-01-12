using System;
using System.Collections.Generic;
using ilPSP.LinSolvers;
using System.Runtime.InteropServices;
using System.Xml;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP;
using ilPSP.Utils;
using ilPSP.Connectors;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Tecplot;
using System.Linq;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Application.ExternalBinding.MatlabCutCellQuadInterface {
    /// <summary>
    /// MatlabCutCellQuadInterface
    /// </summary>
    public class MatlabCutCellQuadInterface {

        /// <summary>
        /// Constructor for export.
        /// </summary>
        [CodeGenExport]
        public MatlabCutCellQuadInterface() {

        }

        // Fields and properties
        GridCommons grd;

        List<_2D> levelsets2D;
        List<_3D> levelsets3D;

        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesA;
        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesB;
        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesInterface;

        LevelSetTracker lsTrk;

        static bool mustFinalizeMPI;

        // Methods

        /// <summary>
        /// the list of 2d level sets
        /// </summary>
        public List<_2D> Levelsets2D { get => levelsets2D ?? (levelsets2D = new List<_2D>()); }

        /// <summary>
        /// the list of 3d level sets
        /// </summary>
        public List<_3D> Levelsets3D { get => levelsets3D ?? (levelsets3D = new List<_3D>()); }

        static void Main(string[] args) {
            
            Console.WriteLine("External binder for Matlab");
			Console.WriteLine("Running an example 2d ellipse test case");
			MatlabCutCellQuadInterfaceTests.ellipse2D();
		}


        /// <summary>
        /// Load/lookup of native libraries
        /// </summary>
        public void BoSSSInitialize() {
            try {
                Console.WriteLine("BoSSS has been initialized");
                mustFinalizeMPI |= BoSSS.Solution.Application.InitMPI();
            } catch (Exception ex) {
                Console.WriteLine(ex);
                throw;
            }
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        public void BoSSSFinalize() {
            if (mustFinalizeMPI)
                MPI.Wrappers.csMPI.Raw.mpiFinalize();

            Console.WriteLine("BoSSS has been finalized");

        }

        /// <summary>
        /// Set the domain in 2D
        /// </summary>
        /// <param name="Dim"></param>
        /// <param name="xNodes"></param>
        /// <param name="yNodes"></param>
        /// <exception cref="ArgumentException"></exception>
        public void SetDomain(int Dim, double[] xNodes, double[] yNodes) {
            if (Dim != 2) {
                throw new ArgumentException("Dimension must be 2 for the given input parameters.");
            }
            grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
        }


        public void SetDomain(int Dim, double[] xNodes, double[] yNodes, double[] zNodes) {
            if (Dim != 3) {
                throw new ArgumentException("Dimension must be 3 for the given input parameters.");
            }
            if (zNodes == null || zNodes.Length == 0) {
                throw new ArgumentException("zNodes must not be empty for a 3D grid.");
            }
            grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
        }

        /// <summary>
        /// When multiple level sets are supplied, this method returns a delegate that gives the maximum value from the list of level sets.
        /// </summary>
        /// <param name="delegates"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private _2D ReturnMaxDelegate(IList<_2D> delegates, string CalculateMethod = "Max") {
            if (delegates.Count == 0)
                throw new Exception("No level sets are submitted, call Submit2DLevelSet() method first!");

            if (delegates.Count == 1)
                return delegates.First();

            switch (CalculateMethod) {
                case "Max":
                    return (x0, x1) => delegates.Max(del => del(x0, x1));
                case "Min":
                    return (x0, x1) => delegates.Min(del => del(x0, x1));
                default:
                    throw new ArgumentException("Invalid CalculateMethod. Use 'Max' or 'Min'.");
			}

        }

        /// <summary>
        /// When multiple level sets are supplied, this method returns a delegate that gives the maximum value from the list of level sets.
        /// </summary>
        /// <param name="delegates"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private _3D ReturnMaxDelegate(IList<_3D> delegates, string CalculateMethod = "Max") {
            if (delegates.Count == 0)
                throw new Exception("No level sets are submitted, call Submit3DLevelSet() method first!");

            if (delegates.Count == 1)
                return delegates.First();

            switch (CalculateMethod) {
                case "Max":
                    return (x0, x1, x2) => delegates.Max(del => del(x0, x1, x2));
                case "Min":
                    return (x0, x1, x2) => delegates.Min(del => del(x0, x1, x2));
                default:
                    throw new ArgumentException("Invalid CalculateMethod. Use 'Max' or 'Min'.");
			}
        }

        /// <summary>
        /// Submit a level set in 2D
        /// </summary>
        /// <param name="inLevelSet">A delegate with (double, double) => double </param>
        public void Submit2DLevelSet(_2D inLevelSet) {
            if(grd == null)
                throw new InvalidOperationException("Grid must be set");

            if (grd.SpatialDimension != 2)
                throw new Exception($"Mismatch in the spatial dimension of the grid ({grd.SpatialDimension}D) with the level set (2D).");

            Levelsets2D.Add(inLevelSet);
            Console.WriteLine("Submitted a 2D level set function.");
        }

        /// <summary>
        /// Submit a level set in 3D
        /// </summary>
        /// <param name="inLevelSet">A delegate with (double, double, double) => double </param>
        public void Submit3DLevelSet(_3D inLevelSet) {
            if(grd == null)
                throw new InvalidOperationException("Grid must be set");

            if (grd.SpatialDimension != 3)
                throw new Exception($"Mismatch in the spatial dimension of the grid ({grd.SpatialDimension}D) with the level set (3D).");

            Levelsets3D.Add(inLevelSet);
            Console.WriteLine("Submitted a 3D level set function.");
        }

        /// <summary>
        /// Project the submitted level sets with the classic HMF quadrature rule
        /// </summary>
        /// <param name="degree">degree of the level set</param>
        public void ProjectLevelSetWithClassic(int degree) {
            Console.WriteLine("Using classic HMF quadrature rule");
            ProjectLevelSet(degree, XQuadFactoryHelperBase.MomentFittingVariants.Classic);
        }

        /// <summary>
        /// Project the submitted level sets with the classic HMF quadrature rule
        /// </summary>
        /// <param name="degree">degree of the level set</param>
        public void ProjectLevelSetWithClassic(int degree) {
            Console.WriteLine("Using classic HMF quadrature rule");
            ProjectLevelSet(degree, CutCellQuadratureMethod.Classic);
        }

        /// <summary>
        /// Project the submitted level sets with the GaussAndStokes HMF quadrature rule
        /// </summary>
        /// <param name="degree">degree of the level set</param>
        public void ProjectLevelSetWithGaussAndStokes(int degree) {
            Console.WriteLine("Using GaussAndStokes HMF quadrature rule");
            ProjectLevelSet(degree, CutCellQuadratureMethod.OneStepGaussAndStokes);
        }

        /// <summary>
        /// Project the submitted level sets with the TwoStepStokesAndGauss HMF quadrature rule
        /// </summary>
        /// <param name="degree">degree of the level set</param>
        public void ProjectLevelWithTwoStepStokesAndGauss(int degree) {
            Console.WriteLine("Using GaussAndStokes HMF quadrature rule");
            ProjectLevelSet(degree, CutCellQuadratureMethod.TwoStepStokesAndGauss);
        }

        /// <summary>
        /// Project the submitted level sets with classic (for compatibility reaons keep the same signature)
        /// </summary>
        /// <param name="degree">degree of the level set</param>
        public void ProjectLevelSet(int degree) {
            Console.WriteLine("Using classic HMF quadrature rule");
            ProjectLevelSet(degree, CutCellQuadratureMethod.Classic);
        }

        /// <summary>
        /// Project the submitted level sets
        /// </summary>
        /// <param name="degree">degree of the level set</param>
        public void ProjectLevelSet(int degree, CutCellQuadratureMethod cellQuadratureMethod) {
            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");

            // Projection
            if (grd.SpatialDimension == 2) {
                _2D inLevelSet = ReturnMaxDelegate(Levelsets2D, CalculateMethod);
                levSet0.ProjectField(inLevelSet);

            } else if (grd.SpatialDimension == 3) {
                _3D inLevelSet = ReturnMaxDelegate(Levelsets3D, CalculateMethod);
                levSet0.ProjectField(inLevelSet);

            } else {
                throw new Exception("Only 2D and 3D meshes are supported.");
            }

            lsTrk = new LevelSetTracker(grd.GridData, cellQuadratureMethod, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);
            Console.WriteLine($"Successful projection of level set with {cellQuadratureMethod.ToString()}");
        }

        /// <summary>
        /// Plot the current state as a .plt file
        /// </summary>
        /// <param name="superSampling">the super sampling level</param>
        public void PlotCurrentState(int superSampling=0) {
            Tecplot tecplot = new Tecplot(grd.GridData, (uint)superSampling);
            string path = Path.Combine(Path.GetFullPath("."), "plot_LS");
            double t = 0.0;
            DGField[] LevelSets = lsTrk.LevelSets.Select(s => (DGField)s).ToArray();

            tecplot.PlotFields(path, t, LevelSets);
        }

        /// <summary>
        /// Checks if it is a cut cell, if so returns true. 
        /// </summary>
        /// <param name="jCell">local cell index</param>
        /// <param name="levelSetIndex">level set index (default: 0)</param>
        /// <returns></returns>
        public bool IsItACutCell(int jCell, int levelSetIndex = 0) {
            var mask = lsTrk.Regions.GetCutCellMask4LevSet(levelSetIndex);                    
            return mask.Contains(jCell);
        }

        /// <summary>
        /// Compile quadrature rules for the given degree and spieces
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="spec">Integer value for the phase:
        /// 1 - external points; -1 - inner points; 0 - boundary points </param>
        public void CompileQuadRules(int deg, int SpeciesId = -1) {
            // If interface, do not bother to calculate others
            if (SpeciesId == 0) {
                CompileLevelsetQuadRules(deg);
                return;
            }

            var spec = SpeciesId == -1 ? lsTrk.GetSpeciesId("A") : lsTrk.GetSpeciesId("B");

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spec }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetVolumeQuadScheme(spec);
            Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rules = scheme.Compile(grd.GridData, deg);

            if (SpeciesId == -1)
                rulesA = rules;
            else
                rulesB = rules;
        }

        /// <summary>
        /// Compile quadrature rules for the given degree and level set index
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="levelSetIndex">Integer value for the level set</param>
        public void CompileLevelsetQuadRules(int deg, int levelSetIndex = 0) {
            var metrics = lsTrk.GetXDGSpaceMetrics();
            var scheme = metrics.XQuadSchemeHelper.GetLevelSetQuadScheme(levelSetIndex, lsTrk.Regions.GetCutCellMask4LevSet(levelSetIndex));
            var rules = scheme.Compile(grd.GridData, deg);

            rulesInterface = rules;
        }


        /// <summary>
        /// Get the quadrature rules for the cell in a multidimensional array
        /// </summary>
        /// <param name="cellNo">local cell index</param>
        /// <param name="spec">Integer value for the phase:
        /// 1 - external points; -1 - inner points; 0 - boundary points </param>
        /// <returns>Returns the quadrature rules with a multidimensional array where the first index is node index,
        /// the second index indicates dimension or weight (when D+1)</returns>
        /// <exception cref="NotSupportedException"></exception>
        public MultidimensionalArray GetQuadRules(int cellNo, int spec = -1) {
            ICompositeQuadRule<QuadRule> rules = spec == 1 ? rulesB : (spec == 0 ? rulesInterface : rulesA);

            MultidimensionalArray ret = null;

            var JacobiDet = grd.GridData.iGeomCells.JacobiDet;

            foreach (var pair in rules) {
                var qr = pair.Rule;

                for (int jCell = pair.Chunk.i0; jCell < pair.Chunk.JE; jCell++) {
                    if (jCell != cellNo)
                        continue;

                    if (!grd.GridData.Cells.IsCellAffineLinear(jCell)) {
                        throw new NotSupportedException("curved cells not supported!");
                    }

                    //qr.OutputQuadratureRuleAsVtpXML("NodesJ" + jCell + ".vtp");
                    //var globTr = qr.CloneAs();
                    var globTr = qr.Nodes.TransformLocal2Global(grd.GridData, jCell);

#if DEBUG
                    globTr.OutputQuadratureRuleAsVtpXML(qr.Weights, "NodestransformedJ" + jCell + ".vtp");
#endif                    

                    double metric_jCell = JacobiDet[jCell];
                    var WeightsGlobal_jCell = qr.Weights.CloneAs();
                    WeightsGlobal_jCell.Scale(metric_jCell);

                    int SpatialDim = grd.SpatialDimension;
                    ret = MultidimensionalArray.Create(qr.NoOfNodes, SpatialDim + 1);

                    for (int n = 0; n < qr.NoOfNodes; n++) {
                        for (int d=0; d < qr.SpatialDim; d++) {
                            ret[n, d] = globTr[n, d];
                        }
                        ret[n, SpatialDim] = WeightsGlobal_jCell[n];
                    }
                    return ret;
                }
            }

            //if (ret == null)
            //    throw new ArgumentOutOfRangeException($"jCell{cellNo} could not be found");

            return ret;

        }

        /// <summary>
        /// Debugging for volume quad rules
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="spec"></param>
        /// <exception cref="NotSupportedException"></exception>
        public void WriteVolQuadRules(int deg, int spec = -1) {
            var spcA = spec == -1 ? lsTrk.GetSpeciesId("A") : lsTrk.GetSpeciesId("B");

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spcA }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetVolumeQuadScheme(spcA);
            var rule = scheme.Compile(grd.GridData, deg);

            var JacobiDet = grd.GridData.iGeomCells.JacobiDet;

            foreach (var pair in rule) {
                var qr = pair.Rule;
                var NodesGlobal = grd.GridData.GlobalNodes.GetValue_Cell(qr.Nodes, pair.Chunk.i0, pair.Chunk.Len);

                for (int jCell = pair.Chunk.i0; jCell < pair.Chunk.JE; jCell++) {
                    int j = jCell - pair.Chunk.i0;
                    if (!grd.GridData.Cells.IsCellAffineLinear(jCell)) {
                        throw new NotSupportedException("curved cells not supported!");
                    }

                    qr.OutputQuadratureRuleAsVtpXML("NodesJ" + jCell + ".vtp");
                    var globTr = qr.Nodes.TransformLocal2Global(grd.GridData, jCell);
                    globTr.OutputQuadratureRuleAsVtpXML(qr.Weights, "NodestransformedJ" + jCell + ".vtp");


                    //var NodesGlobal_jCell = NodesGlobal.ExtractSubArrayShallow(j, -1, -1);
                    //double metric_jCell = JacobiDet[jCell];
                    //var WeightsGlobal_jCell = qr.Weights.CloneAs();
                    //WeightsGlobal_jCell.Scale(metric_jCell);
                }
            }
            Console.WriteLine("Calculated the volume quadrature rule");
        }

        /// <summary>
        /// Debugging for surface quad rules
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="spec"></param>
        /// <exception cref="NotSupportedException"></exception>
        public void WriteSrfQuadRules(int deg) {

            int iLevSet = 0;

            var metrics = lsTrk.GetXDGSpaceMetrics();
            var scheme = metrics.XQuadSchemeHelper.GetLevelSetQuadScheme(iLevSet, lsTrk.Regions.GetCutCellMask4LevSet(iLevSet));
            var rule = scheme.Compile(grd.GridData, deg);

            var JacobiDet = grd.GridData.iGeomCells.JacobiDet;

            foreach (var pair in rule) {
                var qr = pair.Rule;
                var NodesGlobal = grd.GridData.GlobalNodes.GetValue_Cell(qr.Nodes, pair.Chunk.i0, pair.Chunk.Len);

                for (int jCell = pair.Chunk.i0; jCell < pair.Chunk.JE; jCell++) {
                    int j = jCell - pair.Chunk.i0;
                    if (!grd.GridData.Cells.IsCellAffineLinear(jCell)) {
                        throw new NotSupportedException("curved cells not supported!");
                    }

                    var NodesGlobal_jCell = NodesGlobal.ExtractSubArrayShallow(j, -1, -1);


                    double metric_jCell = JacobiDet[jCell];


                    var WeightsGlobal_jCell = qr.Weights.CloneAs();
                    WeightsGlobal_jCell.Scale(metric_jCell);
                }
            }


        }


    }

}

