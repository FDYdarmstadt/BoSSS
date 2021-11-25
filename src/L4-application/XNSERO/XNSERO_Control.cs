using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
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

    [Serializable]
    public class XNSERO_Control : XNSE_Control {

        public XNSERO_Control() {
            //Debugger.Launch();
            base.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Set default values to LevelSet (one could still overwrite those)
            AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));
            Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.Prescribed;
            AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            TimeSteppingScheme = TimeSteppingScheme.BDF2;
            NonlinearCouplingSolidFluid = true;
            UseImmersedBoundary = true;
            CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public XNSERO_Control(int degree, string projectName, string projectDescription = null, List<string> tags = null) : this() {
            ProjectName = projectName;
            SessionName = projectName;
            ProjectDescription = projectDescription;
            if(tags != null) {
                for(int i = 0; i < tags.Count(); i++) {
                    Tags.Add(tags[i]);
                }
            }
            SetDGdegree(degree);
        }

        /// <summary>
        /// always true
        /// </summary>
        public override bool UseImmersedBoundary {
            get {
                return true;
            }
            set {
                if(value == false) {
                    Console.Error.WriteLine($"Immersed Boundary cannot be turned off for {typeof(XNSERO).Name}.");    
                }
            }
        }

        [DataMember]
        public bool ContainsSecondFluidSpecies = false;

        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            SetFieldOptions(p, Math.Max(4, 2 * p));

            FieldOptions.Add(VariableNames.Phoretic, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }

        /// <summary>
        /// Set true during restart.
        /// </summary>
        [DataMember]
        public bool IsRestart = false;

        /// <summary>
        /// List of the boundary values at the domain boundary.
        /// </summary>
        [DataMember]
        readonly List<string> m_BoundaryValues = new List<string>();

        /// <summary>
        /// The position of all boundaries, independent of boundary condition.
        /// </summary>
        [DataMember]
        public double[][] BoundaryPositionPerDimension;

        /// <summary>
        /// The position of all boundaries, only for walls.
        /// </summary>
        [DataMember]
        public double[][] WallPositionPerDimension;

        /// <summary>
        /// True for periodic walls
        /// </summary>
        [DataMember]
        public bool[] BoundaryIsPeriodic;

        /// <summary>
        /// Max grid length
        /// </summary>
        [DataMember]
        public double MaxGridLength;

        /// <summary>
        /// Min grid length
        /// </summary>
        [DataMember]
        public double MinGridLength;

        /// <summary>
        /// Dimension of the domain.
        /// </summary>
        [DataMember]
        public int SpatialDimension;

        /// <summary>
        /// Add database and the frequency of saves.
        /// </summary>
        /// <param name="dataBasePath"></param>
        /// <param name="savePeriod"></param>
        public void SetSaveOptions(string dataBasePath = null, int savePeriod = 1) {
            if(dataBasePath != null) {
                savetodb = true;
                DbPath = dataBasePath;
                saveperiod = savePeriod;
            } else
                savetodb = false;
        }

        /// <summary>
        /// Setup AMR Level at the level set.
        /// </summary>
        /// <param name="amrLevel"></param>
        public void SetAddaptiveMeshRefinement(int amrLevel) {
            if(amrLevel == 0)
                return;
            AdaptiveMeshRefinement = true;
            RefinementLevel = amrLevel;
            AMR_startUpSweeps = amrLevel;
        }

        /// <summary>
        /// Setup boundary values, e.g. wall, pressure_dirichlet...
        /// </summary>
        /// <param name="boundaryValues"></param>
        public void SetBoundaries(List<string> boundaryValues) {
            if(boundaryValues.Count() > 4)
                throw new NotImplementedException("max 4 boundary values");
            for(int i = 0; i < boundaryValues.Count(); i++) {
                AddBoundaryValue(boundaryValues[i]);
                m_BoundaryValues.Add(boundaryValues[i]);
            }
        }

        /// <summary>
        /// Set the grid in a rectangular domain. Currently only 2D!
        /// </summary>
        /// <param name="lengthX"></param>
        /// <param name="lengthY"></param>
        /// <param name="cellsPerUnitLength"></param>
        /// <param name="periodicX"></param>
        /// <param name="periodicY"></param>
        public void SetGrid(double lengthX, double lengthY, double cellsPerUnitLength, bool periodicX = false, bool periodicY = false, bool periodicParticles = false) {
            SpatialDimension = 2;
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
            if (periodicParticles) {
                BoundaryIsPeriodic[0] = periodicX;
                BoundaryIsPeriodic[1] = periodicY;
            }
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

        /// <summary>
        /// Set time-steps and length of the simulation.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="noOfTimesteps"></param>
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

        /// <summary>
        /// Set gravity for particles and fluid species.
        /// </summary>
        /// <param name="Gravity"></param>
        public void SetGravity(Vector Gravity) {
            this.Gravity = new Vector(Gravity);
            InitialValues_Evaluators.Add("GravityX#A", X => Gravity[0]);
            InitialValues_Evaluators.Add("GravityX#B", X => Gravity[0]);
            InitialValues_Evaluators.Add("GravityY#A", X => Gravity[1]);
            InitialValues_Evaluators.Add("GravityY#B", X => Gravity[1]);
        }

        /// <summary>
        /// Returns gravity vector.
        /// </summary>
        /// <returns></returns>
        public Vector GetGravity() => Gravity;

        /*
        /// <summary>
        /// See <see cref="LevelSetHandling"/>, Lie-Splitting with iterative coupling by default.
        /// </summary>
        [DataMember]
        public override LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
        */

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
        /// All particles
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

        /// <summary>
        /// type of appropriate solver
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XNSERO);
        }

        /// <summary>
        /// switch to turn the Phoretic Field on/off
        /// </summary>
        [DataMember]
        public bool UsePhoreticField = false;

        /// <summary>
        /// switch to turn the averaged equations on/off
        /// </summary>
        [DataMember]
        public bool UseAveragedEquations = false;
    }
}
