using XESF.Fluxes;
using ApplicationWithIDT;
using XESF;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;

namespace XESTSF {
    public class XESTSFControl : IDTControl {
        public XESTSFControl():base() {
            base.quadOrderFunc = (int[] A, int[] B, int[] C) =>  Math.Abs(2*A.Max()) + Math.Abs(C.Max()) + Math.Max(this.LevelSetDegree,this.LevelSetTwoDegree);
            base.PartiallyFixLevelSetForSpaceTime = true;
        }
        public static void UpdateControl(XESTSFControl oldcontrol, double timestep, DGField[] previous_u,int _tNref=1, int MaxIterations=20) {
            var c = oldcontrol;
            c.previous_u= previous_u;
            c.NoOfTimesteps = MaxIterations;
            var grddat = (GridData)previous_u[0].GridDat;
            double tMin = grddat.GlobalBoundingBox.Max[grddat.SpatialDimension - 1];

            //change grid
            {
                double xMin = grddat.GlobalBoundingBox.Min[0];
                double xMax = grddat.GlobalBoundingBox.Max[0];
                
                double tMax = tMin+ timestep;
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, grddat.Cells.Count/c.tNref + 1);
                    double[] tNodes = GenericBlas.Linspace(tMin, tMax, _tNref + 1);
                    GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, tNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "SpaceTimeBoundary");
                    grid.DefineEdgeTags(delegate (double[] X) {
                        if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                            return 2;
                        } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                            return 1;
                        } else if(Math.Abs(X[1] - tMax) < 1e-14) { // top boundary
                            return 2;
                        } else { //bottom boundary
                            return 3;
                        }
                    });
                    return grid;
                };
                c.tNref = _tNref;
                //c.ProjectName=c.ProjectName + "Slab_t" + tMin+ "_to_t" + (tMin+timestep);
            }

        }
        public bool hasDirichletBoundary = false;
        public Func<double[], double[]> DirichletBoundaryFunc { get; set; }

        public override Type GetSolverType() {
            return typeof(XESTSFMain);
        }

        public int tNref =1;
        public XESFMain old_main;
        internal DGField[] previous_u;

        public string PointPath { get; set; }
        
        public ConvectiveBulkFluxes ConvectiveBulkFlux { get; set; } = ConvectiveBulkFluxes.OptimizedHLLC;

        public FluxVersion FluxVersion { get; set; } = FluxVersion.Optimized;

        public ConvectiveInterfaceFluxes ConvectiveInterfaceFlux_LsOne { get; set; } = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var;

        public ConvectiveInterfaceFluxes ConvectiveInterfaceFlux_LsTwo { get; set; } = ConvectiveInterfaceFluxes.OptimizedHLLCInterface;
        public int IVTimestepNumber { get; set; } = 0;
        public int StartDegree { get; set; } = 0;
        public double ExactEnthalpy { get; internal set; }
        public double NormfacAcoustic { get; internal set; } = 1000;
        public bool DoSpaceTimeSlabs { get; internal set; } = false;

        //public SensorTypes SensorType { get; internal set; };
        //public string SensorVariable { get; internal set; } = null;
        //public double SensorLimit { get; internal set; } = double.MinValue;
        //public ArtificialViscosityLawTypes ArtificialViscosityLawType { get; internal set; }
        //public DiffusiveBulkFluxes DiffusiveBulkFlux { get; internal set; }
        //public DiffusiveInterfaceFluxes DiffusiveInterfaceFlux { get; internal set; }

    }

}