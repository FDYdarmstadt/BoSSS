using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSERO_Solver {

    public class XNSERO_Control : XNSE_Control {

        public XNSERO_Control() {

        }

        /// <summary>
        /// ctor
        /// </summary>
        public XNSERO_Control(int degree, string projectName, string projectDescription = null, List<string> tags = null) : base() {
            ProjectName = projectName;
            SessionName = projectName;
            ProjectDescription = projectDescription;
            if(tags != null) {
                for(int i = 0; i < tags.Count(); i++) {
                    Tags.Add(tags[i]);
                }
            }
            SetDGdegree(degree);
            // Set default values to LevelSet (one could still overwrite those)
            InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), X => -1);
            Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.Prescribed;
            AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            UseImmersedBoundary = true;
            TimeSteppingScheme = TimeSteppingScheme.BDF2;
            NonlinearCouplingSolidFluid = true;
        }

        public bool ContainsSecondFluidSpecies = false;

        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            SetFieldOptions(p, Math.Max(12, p));
        }

        /// <summary>
        /// Set true during restart.
        /// </summary>
        public bool IsRestart = false;

        /// <summary>
        /// List of the boundary values at the domain boundary.
        /// </summary>
        readonly List<string> m_BoundaryValues = new List<string>();

        /// <summary>
        /// The position of all boundaries, independent of boundary condition.
        /// </summary>
        public double[][] BoundaryPositionPerDimension;

        /// <summary>
        /// The position of all boundaries, only for walls.
        /// </summary>
        public double[][] WallPositionPerDimension;

        /// <summary>
        /// True for periodic walls
        /// </summary>
        public bool[] BoundaryIsPeriodic;

        /// <summary>
        /// Max grid length
        /// </summary>
        public double MaxGridLength;

        /// <summary>
        /// Min grid length
        /// </summary>
        public double MinGridLength;

        public void SetSaveOptions(string dataBasePath = null, int savePeriod = 1) {
            if(dataBasePath != null) {
                savetodb = true;
                DbPath = dataBasePath;
                saveperiod = savePeriod;
            } else
                savetodb = false;
        }

        public void SetAddaptiveMeshRefinement(int amrLevel) {
            if(amrLevel == 0)
                return;
            AdaptiveMeshRefinement = true;
            RefinementLevel = amrLevel;
            AMR_startUpSweeps = amrLevel;
        }

        public void SetBoundaries(List<string> boundaryValues) {
            if(boundaryValues.Count() > 4)
                throw new NotImplementedException("max 4 boundary values");
            for(int i = 0; i < boundaryValues.Count(); i++) {
                AddBoundaryValue(boundaryValues[i]);
                m_BoundaryValues.Add(boundaryValues[i]);
            }
        }

        public void SetGrid(double lengthX, double lengthY, double cellsPerUnitLength, bool periodicX = false, bool periodicY = false) {
            MaxGridLength = 1 / cellsPerUnitLength;
            BoundaryPositionPerDimension = new double[2][];
            WallPositionPerDimension = new double[2][];
            WallPositionPerDimension[0] = new double[2];
            WallPositionPerDimension[1] = new double[2];
            BoundaryIsPeriodic = new bool[2];
            BoundaryPositionPerDimension[0] = new double[] { -lengthX / 2, lengthX / 2 };
            BoundaryPositionPerDimension[1] = new double[] { -lengthY / 2, lengthY / 2 };
            BoundaryIsPeriodic[0] = periodicX;
            BoundaryIsPeriodic[1] = periodicY;
            if(IsRestart)
                return;
            if(m_BoundaryValues.IsNullOrEmpty() && !BoundaryIsPeriodic[0] && !BoundaryIsPeriodic[1])
                SetBoundaries(new List<string> { "Wall" });
            GridFunc = delegate {
                int q = new int(); // #Cells in x-dircetion + 1
                int r = new int(); // #Cells in y-dircetion + 1

                q = (int)(cellsPerUnitLength * lengthX);
                r = (int)(cellsPerUnitLength * lengthY);

                double[] Xnodes = GenericBlas.Linspace(-lengthX / 2, lengthX / 2, q + 1);
                double[] Ynodes = GenericBlas.Linspace(-lengthY / 2, lengthY / 2, r + 1);

                Grid2D grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: periodicX, periodicY: periodicY);

                for(int i = 0; i < m_BoundaryValues.Count(); i++) {
                    byte iB = (byte)(i + 1);
                    grd.EdgeTagNames.Add(iB, m_BoundaryValues[i]);
                }

                if(m_BoundaryValues.Count() == 0 && periodicX == false && periodicY == false)
                    throw new Exception("Please specify boundaries before creating the grid");

                if(m_BoundaryValues.Count() == 1) {
                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 1;
                        if(m_BoundaryValues[0].Contains("wall") || m_BoundaryValues[0].Contains("Wall")) {
                            WallPositionPerDimension[0][0] = -lengthX / 2;
                            WallPositionPerDimension[0][1] = lengthX / 2;
                            WallPositionPerDimension[1][0] = -lengthY / 2;
                            WallPositionPerDimension[1][1] = lengthY / 2;
                        }
                        return et;
                    });
                } else {
                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[0] - (-lengthX / 2)) <= 1.0e-8) {
                            for(int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if(m_BoundaryValues[i].Contains("left") || m_BoundaryValues[i].Contains("Left")) {
                                    et = (byte)(i + 1);
                                    if(m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[0][0] = -lengthX / 2;
                                    }
                                }
                            }
                        }
                        if(Math.Abs(X[0] + (-lengthX / 2)) <= 1.0e-8) {
                            for(int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if(m_BoundaryValues[i].Contains("right") || m_BoundaryValues[i].Contains("Right")) {
                                    et = (byte)(i + 1);
                                    if(m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[0][1] = lengthX / 2;
                                    }
                                }
                            }
                        }

                        if(Math.Abs(X[1] - (-lengthY / 2)) <= 1.0e-8) {
                            for(int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if(m_BoundaryValues[i].Contains("lower") || m_BoundaryValues[i].Contains("Lower")) {
                                    et = (byte)(i + 1);
                                    if(m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[1][0] = -lengthY / 2;
                                    }
                                }
                            }
                        }
                        if(Math.Abs(X[1] + (-lengthY / 2)) <= 1.0e-8) {
                            for(int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if(m_BoundaryValues[i].Contains("upper") || m_BoundaryValues[i].Contains("Upper")) {
                                    et = (byte)(i + 1);
                                    if(m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[1][1] = lengthY / 2;
                                    }
                                }
                            }
                        }
                        Debug.Assert(et != 0);
                        return et;
                    });
                }
                Console.WriteLine("Cells:" + grd.NumberOfCells);
                return grd;
            };
        }


        public void SetTimesteps(double dt, int noOfTimesteps) {
            dtMax = dt;
            dtMin = dt;
            dtFixed = dt;
            Endtime = noOfTimesteps * dt;
            NoOfTimesteps = noOfTimesteps;
        }

        /// <summary>
        /// coefficient of restitution
        /// </summary>
        [DataMember]
        public double CoefficientOfRestitution = 1.0;

        /// <summary>
        /// Gravity acting on the particles, zero by default.
        /// </summary>
        [DataMember]
        private Vector Gravity = new Vector(0, 0);

        public void SetGravity(Vector Gravity) {
            this.Gravity = new Vector(Gravity);
            InitialValues_Evaluators.Add("GravityX#A", X => Gravity[0]);
            InitialValues_Evaluators.Add("GravityX#B", X => Gravity[0]);
            InitialValues_Evaluators.Add("GravityY#A", X => Gravity[1]);
            InitialValues_Evaluators.Add("GravityY#B", X => Gravity[1]);
        }

        public Vector GetGravity() => Gravity;

        /// <summary>
        /// See <see cref="LevelSetHandling"/>, Lie-Splitting with iterative coupling by default.
        /// </summary>
        [DataMember]
        public new LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

        public void SetParticles(List<Particle> ParticleList) {
            Particles = ParticleList.ToArray();
            // Initialize particle level-set
            double levelSet(double[] X) {
                double levelSetFunction = int.MinValue;
                for(int p = 0; p < Particles.Count(); p++) {
                    Particle currentParticle = Particles[p];
                    if(levelSetFunction < currentParticle.LevelSetFunction(X, 0))
                        levelSetFunction = currentParticle.LevelSetFunction(X, 0);
                }
                return levelSetFunction;
            }
            InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(1), levelSet);
            Option_LevelSetEvolution2 = Solution.LevelSetTools.LevelSetEvolution.RigidObject;
        }

        /// <summary>
        /// All particles in the FSI
        /// </summary>
        [DataMember]
        public Particle[] Particles { get; private set; }

        /// <summary>
        /// Fix all particles
        /// </summary>
        [DataMember]
        public bool fixPosition = false;

        /// <summary>
        /// For Collisions
        /// </summary>
        [DataMember]
        public double minDistanceThreshold = 0;
    }
}
