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
        
        public void SetDomain(int Dim, double[] xNodes, double[] yNodes, double[] zNodes) {

            switch(Dim) {
                case 2: grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                    break; 
                case 3: grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes); 
                    break;
                default: 
                    throw new NotImplementedException();
            }
        }

        public delegate double Matlab2DFunctionDelegate(double x, double y);

        public delegate double Matlab3DFunctionDelegate(double x, double y, double z);

        LevelSetTracker lsTrk;

        public void SetLevelSets(int degree, _2D dummy) {

            Basis b = new Basis(grd, degree);
            var levSet0 = new LevelSet(b, "LevelSetField0");
            levSet0.ProjectField(dummy);

            lsTrk = new LevelSetTracker(grd.GridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, levSet0);
            lsTrk.UpdateTracker(0.0);

        }

        public void PlotCurrentState(int superSampling=0) {
            Tecplot tecplot = new Tecplot(grd.GridData, (uint)superSampling);
            string path = Path.Combine(Path.GetFullPath("."), "plot_");
            double t = 0.0;
            DGField[] LevelSets = lsTrk.LevelSets.Select(s => (DGField)s).ToArray();

            tecplot.PlotFields(path, t, LevelSets);
        }

        public void WriteVolQuadRules(int deg) {

            var spcA = lsTrk.GetSpeciesId("A");

            var metrics = lsTrk.GetXDGSpaceMetrics(new SpeciesId[] { spcA }, deg);
            var scheme = metrics.XQuadSchemeHelper.GetVolumeQuadScheme(spcA);
            var rule = scheme.Compile(grd.GridData, deg);

            var JacobiDet = grd.GridData.iGeomCells.JacobiDet;

            foreach (var pair in rule) {
                var qr = pair.Rule;
                var NodesGlobal = grd.GridData.GlobalNodes.GetValue_Cell(qr.Nodes, pair.Chunk.i0, pair.Chunk.Len);



                for(int jCell = pair.Chunk.i0; jCell < pair.Chunk.JE; jCell++) {
                    int j = jCell - pair.Chunk.i0;
                    if(!grd.GridData.Cells.IsCellAffineLinear(jCell)) {
                        throw new NotSupportedException("curved cells not supported!");
                    }
                    
                    var NodesGlobal_jCell = NodesGlobal.ExtractSubArrayShallow(j, -1, -1);
                    NodesGlobal_jCell.SaveToTextFile("NodesJ" + jCell + ".txt");


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

