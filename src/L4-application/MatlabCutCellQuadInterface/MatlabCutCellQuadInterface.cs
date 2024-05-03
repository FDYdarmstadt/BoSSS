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
using ilPSP.Connectors;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Tecplot;
using System.Linq;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Application.ExternalBinding.MatlabCutCellQuadInterface {
    /// <summary>
    /// Initialization stuff
    /// </summary>
    public class MatlabCutCellQuadInterface {

        /// <summary>
        /// Constructor for export.
        /// </summary>
        [CodeGenExport]
        public MatlabCutCellQuadInterface() {

        }


        static void Main(string[] args) {
            
            Console.WriteLine(args.Length);
            Console.WriteLine("External binder for Matlab");
            MatlabCutCellQuadInterfaceTests.circle2D();
        }

        static bool mustFinalizeMPI;

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
        }

        GridCommons grd;


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

        LevelSetTracker lsTrk;

        public void SetLevelSet(int degree, _2D inLevelSet) {

            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");
            levSet0.ProjectField(inLevelSet);

            lsTrk = new LevelSetTracker(grd.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);
        }

        public void SetLevelSet(int degree, _3D inLevelSet) {

            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");
            levSet0.ProjectField(inLevelSet);

            lsTrk = new LevelSetTracker(grd.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);
        }

        List<_2D> levelsets2D;
        List<_3D> levelsets3D;

        /// <summary>
        /// When multiple level sets are supplied, this method returns a delegate that gives the maximum value from the list for a given pint.
        /// </summary>
        /// <param name="delegates"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private _2D ReturnMaxDelegate(IEnumerable<_2D> delegates) {
            if (delegates.Count() < 1)
                throw new Exception("No level sets are submitted, call SubmitLevelSet() method first!");

            if (delegates.Count() == 1)
                return delegates.First();

            return (x0, x1) => {
                double maxResult = double.MinValue;
                foreach (_2D del in delegates) {
                    double result = del(x0, x1);
                    if (result > maxResult) {
                        maxResult = result;
                    }
                }
                return maxResult;
            };
        }

        /// <summary>
        /// When multiple level sets are supplied, this method returns a delegate that gives the maximum value from the list for a given pint.
        /// </summary>
        /// <param name="delegates"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private _3D ReturnMaxDelegate(IEnumerable<_3D> delegates) {
            if (delegates.Count() < 1)
                throw new Exception("No level sets are submitted, call SubmitLevelSet() method first!");

            if (delegates.Count() == 1)
                return delegates.First();

            return (x0, x1, x2) =>
            {
                double maxResult = double.MinValue;
                foreach (_3D del in delegates) {
                    double result = del(x0, x1, x2);
                    if (result > maxResult) {
                        maxResult = result;
                    }
                }
                return maxResult;
            };
        }

        public void SubmitLevelSet(_2D inLevelSet) {
            Levelsets2D.Add(inLevelSet);
        }

        public void SubmitLevelSet(_3D inLevelSet) {
            Levelsets3D.Add(inLevelSet);
        }


        public void ProjectLevelSet(int degree) {
            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");

            // Projection
            if (grd.SpatialDimension == 2) {
                _2D inLevelSet = ReturnMaxDelegate(Levelsets2D);
                levSet0.ProjectField(inLevelSet);

            } else if (grd.SpatialDimension == 3) {
                _3D inLevelSet = ReturnMaxDelegate(Levelsets3D);
                levSet0.ProjectField(inLevelSet);

            } else {
                throw new Exception("Only 2D and 3D meshes are supported.");
            }


            lsTrk = new LevelSetTracker(grd.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);
            Console.WriteLine("Successful creation of level set");
        }


        public void SetLevelSets(int degree, _3D[] inLevelSets) {
            _3D inLevelSet = ReturnMaxDelegate(inLevelSets);
            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");
            levSet0.ProjectField(inLevelSet);

            lsTrk = new LevelSetTracker(grd.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);
        }


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
            // If interface does not bother to calculate others
            if (SpeciesId == 0) {
                CompileLevelsetQuadRules(deg);
                return;
            }

            var spcA = SpeciesId == -1 ? lsTrk.GetSpeciesId("A") : lsTrk.GetSpeciesId("B");

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spcA }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetVolumeQuadScheme(spcA);
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
            var spcA = lsTrk.GetSpeciesId("A");

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spcA }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetLevelSetquadScheme(levelSetIndex, spcA, lsTrk.Regions.GetCutCellMask4LevSet(levelSetIndex));
            var rules = scheme.Compile(grd.GridData, deg);

            rulesInterface = rules;
        }

        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesA;
        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesB;
        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesInterface;

        public List<_2D> Levelsets2D { get => levelsets2D ?? (levelsets2D = new List<_2D>()); }

        public List<_3D> Levelsets3D { get => levelsets3D ?? (levelsets3D = new List<_3D>()); }


        public MultidimensionalArray GetQuadRules(int cellNo, int spec = -1) {
            ICompositeQuadRule<QuadRule> rules = spec == 1 ? rulesB : (spec == 0 ? rulesInterface : rulesA);

            MultidimensionalArray ret = null;

            var JacobiDet = grd.GridData.iGeomCells.JacobiDet;

            foreach (var pair in rules) {
                var qr = pair.Rule;

                for (int jCell = pair.Chunk.i0; jCell < pair.Chunk.JE; jCell++) {
                    if (jCell != cellNo)
                        continue;

                    int j = jCell - pair.Chunk.i0;
                    if (!grd.GridData.Cells.IsCellAffineLinear(jCell)) {
                        throw new NotSupportedException("curved cells not supported!");
                    }

                    //qr.OutputQuadratureRuleAsVtpXML("NodesJ" + jCell + ".vtp");

                    var globTr = qr.CloneAs();
                    globTr.TransformLocal2Global(grd, jCell);
                    //globTr.OutputQuadratureRuleAsVtpXML("NodestransformedJ" + jCell + ".vtp");
                    

                    double metric_jCell = JacobiDet[jCell];
                    var WeightsGlobal_jCell = qr.Weights.CloneAs();
                    WeightsGlobal_jCell.Scale(metric_jCell);


                    ret = MultidimensionalArray.Create(globTr.NoOfNodes, globTr.SpatialDim + 1);


                    for (int n = 0; n < globTr.NoOfNodes; n++) {
                        for (int d=0; d < globTr.SpatialDim; d++) {
                            ret[n, d] = globTr.Nodes[n, d];
                        }
                        ret[n, globTr.SpatialDim] = WeightsGlobal_jCell[n];
                    }

                }
            }

            //if (ret == null)
            //    throw new ArgumentOutOfRangeException($"jCell{cellNo} could not be found");

            return ret;

        }

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
                    var globTr = qr.CloneAs();
                    globTr.TransformLocal2Global(grd, jCell);
                    globTr.OutputQuadratureRuleAsVtpXML("NodestransformedJ" + jCell + ".vtp");


                    //var NodesGlobal_jCell = NodesGlobal.ExtractSubArrayShallow(j, -1, -1);
                    //double metric_jCell = JacobiDet[jCell];
                    //var WeightsGlobal_jCell = qr.Weights.CloneAs();
                    //WeightsGlobal_jCell.Scale(metric_jCell);
                }
            }
            Console.WriteLine("Calculated the volume quadrature rule");
        }


        public void WriteSrfQuadRules(int deg) {

            var spcA = lsTrk.GetSpeciesId("A");

            int iLevSet = 0;

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spcA }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetLevelSetquadScheme(iLevSet, spcA, lsTrk.Regions.GetCutCellMask4LevSet(iLevSet));
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

