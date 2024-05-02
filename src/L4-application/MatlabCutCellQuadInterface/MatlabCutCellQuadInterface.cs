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


        /// <summary>
        /// the main purpose of this property is to guarantee that the assembly for all types are linked.
        /// </summary>
        static public Type[] ExplicitHooks {
            get {
                return new[] {
                    typeof(Foundation.Grid.Classic.GridData),
                    typeof(MatlabCutCellQuadInterface)
                };
            }

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
                Console.WriteLine("Test");
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

        public void SetLevelSets(int degree, _2D dummy) {

            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");
            levSet0.ProjectField(dummy);

            lsTrk = new LevelSetTracker(grd.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);
            PlotCurrentState(0);

        }

        public void PlotCurrentState(int superSampling=0) {
            Tecplot tecplot = new Tecplot(grd.GridData, (uint)superSampling);
            string path = Path.Combine(Path.GetFullPath("."), "plot_LS");
            double t = 0.0;
            DGField[] LevelSets = lsTrk.LevelSets.Select(s => (DGField)s).ToArray();

            tecplot.PlotFields(path, t, LevelSets);
        }


        /// <summary>
        /// Compile quadrature rules for the given degree and spieces
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="spec">Integer value for the phase:
        /// If level set < 0 then -1 and if level set > 0 then 1 </param>
        public void CompileQuadRules(int deg, int SpeciesId = -1) {
            var spcA = SpeciesId == -1 ? lsTrk.GetSpeciesId("A") : lsTrk.GetSpeciesId("B");

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spcA }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetVolumeQuadScheme(spcA);
            Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rules = scheme.Compile(grd.GridData, deg);

            if (SpeciesId == -1)
                rulesA = rules;
            else
                rulesB = rules;
        }

        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesA;
        Foundation.Quadrature.ICompositeQuadRule<Foundation.Quadrature.QuadRule> rulesB;

        public MultidimensionalArray GetQuadRules(int cellNo, int spec = -1) {
            ICompositeQuadRule<QuadRule> rules = spec == -1 ? rulesA : rulesB;

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

            if (ret == null)
                throw new ArgumentOutOfRangeException($"jCell{cellNo} could not be found");

            Console.WriteLine("Calculated the volume quadrature rule");
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


                    var NodesGlobal_jCell = NodesGlobal.ExtractSubArrayShallow(j, -1, -1);

                    double metric_jCell = JacobiDet[jCell];

                    var WeightsGlobal_jCell = qr.Weights.CloneAs();
                    WeightsGlobal_jCell.Scale(metric_jCell);
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

