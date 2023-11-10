using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSERO_Solver {

    [Serializable]
    public class XNSERO_Control : XNSE_Control {

        public XNSERO_Control() {
            base.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            SetDefaultValues();

            void SetDefaultValues() {
                AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));
                Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.Prescribed;
                AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
                TimeSteppingScheme = TimeSteppingScheme.BDF2;
                NonlinearCouplingSolidFluid = true;
                UseImmersedBoundary = true;
                CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

                base.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
                base.NonLinearSolver.ConvergenceCriterion = 1.0e-8;

                Gravity = new(0, 0);
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        public XNSERO_Control(int degree, string projectName, string projectDescription = null, List<string> tags = null, bool IsRestart = false, string PathToOldSessionDir = "", int TimestepNoForRestart = 0) : this() {
            ProjectName = projectName;
            SessionName = projectName;
            ProjectDescription = projectDescription;
            if (tags != null) {
                for (int i = 0; i < tags.Count; i++) {
                    Tags.Add(tags[i]);
                }
            }
            SetDGdegree(degree);

            this.IsRestart = IsRestart;
            this.PathToOldSessionDir = PathToOldSessionDir;
            this.TimestepNoForRestart = TimestepNoForRestart;
        }

        /// <summary>
        /// type of appropriate solver
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XNSERO);
        }

        /// <summary>
        /// Add database and the frequency of saves.
        /// </summary>
        /// <param name="dataBasePath"></param>
        /// <param name="savePeriod"></param>
        public void SetSaveOptions(string dataBasePath = null, int savePeriod = 1) {
            if (dataBasePath != null) {
                savetodb = true;
                DbPath = dataBasePath;
                saveperiod = savePeriod;
            } else
                savetodb = false;
        }

        private readonly bool IsRestart;
        private readonly string PathToOldSessionDir;
        private readonly int TimestepNoForRestart;

        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            SetFieldOptions(p, Math.Max(4, 2 * p));
            //SetFieldOptions(p, p);
            FieldOptions.Add(VariableNames.Phoretic, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }

        /// <summary>
        /// always true
        /// </summary>
        public override bool UseImmersedBoundary {
            get {
                return true;
            }
            set {
                if (value == false) {
                    Console.Error.WriteLine($"Immersed Boundary cannot be turned off for {typeof(XNSERO).Name}.");
                }
            }
        }

        /// <summary>
        /// Flag if a second fluid species is used.
        /// </summary>
        [DataMember]
        public bool ContainsSecondFluidSpecies = false;

        /// <summary>
        /// List of the boundary values at the domain boundary.
        /// </summary>
        [DataMember]
        private readonly List<string> BoundaryValueList = new();

        /// <summary>
        /// The position of all boundaries, independent of boundary condition.
        /// </summary>
        [DataMember]
        public double[][] BoundaryPositionPerDimension { get; private set; }

        /// <summary>
        /// The position of all boundaries, only for walls.
        /// </summary>
        [DataMember]
        public double[][] WallPositionPerDimension { get; private set; }

        /// <summary>
        /// True for periodic walls
        /// </summary>
        [DataMember]
        public bool[] BoundaryIsPeriodic { get; private set; }

        /// <summary>
        /// Max grid length
        /// </summary>
        [DataMember]
        public double MaxGridLength { get; private set; }

        /// <summary>
        /// Setup AMR Level at the level set.
        /// </summary>
        /// <param name="MaxRefinementLevel"></param>
        public void SetAddaptiveMeshRefinement(int MaxRefinementLevel) {
            if (MaxRefinementLevel != 0) {
                AdaptiveMeshRefinement = true;
                RefinementLevel = MaxRefinementLevel;
                AMR_startUpSweeps = MaxRefinementLevel;
                activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = MaxRefinementLevel });
            }
        }

        /// <summary>
        /// All particles
        /// </summary>
        [DataMember]
        public Particle[] Particles { get; private set; }

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

        /// <summary>
        /// Coefficient of restitution for collision model.
        /// </summary>
        [DataMember]
        public double CoefficientOfRestitution = 1.0;

        /// <summary>
        /// Gravity acting on the particles, zero by default.
        /// </summary>
        [DataMember]
        public Vector Gravity { get; private set; }

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
        /// Setup boundary values, e.g. wall, pressure_dirichlet...
        /// </summary>
        /// <param name="boundaryValues"></param>
        public void SetBoundaries(List<string> boundaryValues) {
            if (boundaryValues.Count > 4)
                throw new NotImplementedException("max 4 boundary values");
            for (int i = 0; i < boundaryValues.Count; i++) {
                AddBoundaryValue(boundaryValues[i]);
                BoundaryValueList.Add(boundaryValues[i]);
            }
            //BoundaryValueList.Add("Pressure_Dirichlet_1");
            //BoundaryValueList.Add("Pressure_Dirichlet_2");
        }

        /// <summary>
        /// Set the grid in a rectangular domain. Currently only 2D!
        /// </summary>
        /// <param name="lengthX"></param>
        /// <param name="lengthY"></param>
        /// <param name="cellsPerUnitLength"></param>
        /// <param name="periodicX"></param>
        /// <param name="periodicY"></param>
        public void SetGrid2D(double lengthX, double lengthY, double cellsPerUnitLength, bool periodicX = false, bool periodicY = false) {
            MaxGridLength = 1 / (cellsPerUnitLength);
            Console.WriteLine("max " + 1 / cellsPerUnitLength + " min " + MaxGridLength);

            BoundaryPositionPerDimension = new double[2][];
            BoundaryPositionPerDimension[0] = new double[] { -lengthX / 2, lengthX / 2 };
            BoundaryPositionPerDimension[1] = new double[] { -lengthY / 2, lengthY / 2 };

            WallPositionPerDimension = new double[2][];
            WallPositionPerDimension[0] = new double[2];
            WallPositionPerDimension[1] = new double[2];

            BoundaryIsPeriodic = new bool[2];
            BoundaryIsPeriodic[0] = periodicX;
            BoundaryIsPeriodic[1] = periodicY;

            if (IsRestart)
                return;

            GridFunc = delegate {
                int q = (int)(cellsPerUnitLength * lengthX);
                int r = (int)(cellsPerUnitLength * lengthY);

                double[] Xnodes = GenericBlas.Linspace(-lengthX / 2, lengthX / 2, q + 1);
                double[] Ynodes = GenericBlas.Linspace(-lengthY / 2, lengthY / 2, r + 1);

                Grid2D grid = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: periodicX, periodicY: periodicY);

                for (int i = 0; i < BoundaryValueList.Count; i++) {
                    byte iB = (byte)(i + 1);
                    grid.EdgeTagNames.Add(iB, BoundaryValueList[i]);
                }

                if (BoundaryValueList.Count == 0 && periodicX == false && periodicY == false)
                    throw new Exception("Please specify boundaries before creating the grid");


                if (BoundaryValueList.Count == 1) {
                    grid.DefineEdgeTags(delegate (double[] X) {
                        byte et = 1;
                        if (BoundaryValueList[0].Contains("wall") || BoundaryValueList[0].Contains("Wall")) {
                            WallPositionPerDimension[0][0] = -lengthX / 2;
                            WallPositionPerDimension[0][1] = lengthX / 2;
                            WallPositionPerDimension[1][0] = -lengthY / 2;
                            WallPositionPerDimension[1][1] = lengthY / 2;
                        }
                        return et;
                    });
                } else {
                    grid.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[0] - (-lengthX / 2)) <= 1.0e-8) {
                            for (int i = 0; i < BoundaryValueList.Count; i++) {
                                if (BoundaryValueList[i].Contains("left") || BoundaryValueList[i].Contains("Left")) {
                                    if (BoundaryValueList[i].Contains("wall") || BoundaryValueList[i].Contains("Wall")) {
                                        WallPositionPerDimension[0][0] = -lengthX / 2;
                                    }
                                    return (byte)(i + 1);
                                }
                            }
                        }
                        if (Math.Abs(X[0] + (-lengthX / 2)) <= 1.0e-8) {
                            for (int i = 0; i < BoundaryValueList.Count; i++) {
                                if (BoundaryValueList[i].Contains("right") || BoundaryValueList[i].Contains("Right")) {
                                    if (BoundaryValueList[i].Contains("wall") || BoundaryValueList[i].Contains("Wall")) {
                                        WallPositionPerDimension[0][1] = lengthX / 2;
                                    }
                                    return (byte)(i + 1);
                                }
                            }
                        }
                        if (Math.Abs(X[1] - (-lengthY / 2)) <= 1.0e-8 && Math.Abs(X[0]) <= 2 / cellsPerUnitLength) {
                            for (int i = 0; i < BoundaryValueList.Count; i++) {
                                if (BoundaryValueList[i].Contains("_1")) {
                                    return (byte)(i + 1);
                                }
                            }
                        }
                        if (Math.Abs(X[1] + (-lengthY / 2)) <= 1.0e-8 && Math.Abs(X[0]) <= 2 / cellsPerUnitLength) {
                            for (int i = 0; i < BoundaryValueList.Count; i++) {
                                if (BoundaryValueList[i].Contains("_2")) {
                                    return (byte)(i + 1);
                                }
                            }
                        }
                        if (Math.Abs(X[1] - (-lengthY / 2)) <= 1.0e-8) {
                            for (int i = 0; i < BoundaryValueList.Count; i++) {
                                if (BoundaryValueList[i].Contains("lower") || BoundaryValueList[i].Contains("Lower")) {
                                    if (BoundaryValueList[i].Contains("wall") || BoundaryValueList[i].Contains("Wall")) {
                                        WallPositionPerDimension[1][0] = -lengthY / 2;
                                    }
                                    return (byte)(i + 1);
                                }
                            }
                        }
                        if (Math.Abs(X[1] + (-lengthY / 2)) <= 1.0e-8) {
                            for (int i = 0; i < BoundaryValueList.Count; i++) {
                                if (BoundaryValueList[i].Contains("upper") || BoundaryValueList[i].Contains("Upper")) {
                                    if (BoundaryValueList[i].Contains("wall") || BoundaryValueList[i].Contains("Wall")) {
                                        WallPositionPerDimension[1][1] = lengthY / 2;
                                    }
                                    return (byte)(i + 1);
                                }
                            }
                        }
                        Debug.Assert(et != 0);
                        return et;
                    });
                }
                Console.WriteLine("Cells:" + grid.NumberOfCells);
                if (BoundaryIsPeriodic[0]) {
                    WallPositionPerDimension[0][0] = 0;
                    WallPositionPerDimension[0][1] = 0;
                }
                if (BoundaryIsPeriodic[1]) {

                    WallPositionPerDimension[1][0] = 0;
                    WallPositionPerDimension[1][1] = 0;
                }
                Console.WriteLine("wall  positions" + new Vector(WallPositionPerDimension[0][0], WallPositionPerDimension[0][1]) + "   " + new Vector(WallPositionPerDimension[1][0], WallPositionPerDimension[1][1]));
                return grid;
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
        /// Initialize particle list
        /// </summary>
        /// <param name="ParticleList"></param>
        public void InitialiseParticles(List<Particle> ParticleList) {

            Particles = ParticleList.ToArray();
            if (IsRestart) {
                ParticleList = LoadParticlesOnRestart(PathToOldSessionDir, ParticleList, TimestepNoForRestart);
            } else {
                Particles = ParticleList.ToArray();
                //Console.WriteLine("Init");
                //Console.WriteLine(Particles[0].Motion.GetPosition(0));
                //Console.WriteLine(Particles[1].Motion.GetPosition(0));
                // Initialize particle level-set
                double levelSet(double[] X) {
                    double levelSetFunction = int.MinValue;
                    for (int p = 0; p < Particles.Length; p++) {
                        Particle currentParticle = Particles[p];
                        if (levelSetFunction < currentParticle.LevelSetFunction(X, 0))
                            levelSetFunction = currentParticle.LevelSetFunction(X, 0);
                    }
                    return levelSetFunction;
                }

                InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(1), levelSet);
                Option_LevelSetEvolution2 = Solution.LevelSetTools.LevelSetEvolution.RigidObject;

                Console.WriteLine("Simulation with " + Particles.Length + " particles");
            }
        }

        private static List<Particle> LoadParticlesOnRestart(string pathToOldSessionDir, List<Particle> ParticleList, int timestep = 0) {
            string pathToPhysicalData = Path.Combine(pathToOldSessionDir, "PhysicalData.txt");

            int historyLength = 3;
            string[] records = File.ReadAllLines(pathToPhysicalData);
            int timestepIndexOffset = 0;
            string lastLine = records[records.Length - 1];
            string[] lastLineFields = lastLine.Split(',');
            int lastTimestep = Convert.ToInt32(lastLineFields[0]);
            if (timestep != 0)
                lastTimestep = timestep;
            if (lastTimestep < historyLength + 1)
                throw new Exception("At least " + historyLength + " time-steps necessary for particle restart!");
            for (int r = 1; r < records.Length; r++) {// 0th line does not contain data
                string currentLine = records[r];
                string[] currentLineFields = currentLine.Split(',');
                if (lastTimestep == Convert.ToInt32(currentLineFields[0])) {
                    timestepIndexOffset = r;
                    break;
                }
            }

            for (int t = 0; t < historyLength; t++) {
                for (int p = 0; p < ParticleList.Count; p++) {
                    Particle currentParticle = ParticleList[p];
                    int index = timestepIndexOffset - ParticleList.Count * t + p;
                    string currentLine = records[index];
                    string[] currentLineFields = currentLine.Split(',');
                    double[] position = new double[2];
                    double[] translationalVelocity = new double[2];
                    double[] force = new double[2];
                    double[] duplicateDistance = new double[2];
                    double[] physicalData = currentLineFields.Select(eachElement => Convert.ToDouble(eachElement)).ToArray();
                    position[0] = Convert.ToDouble(currentLineFields[3]);
                    position[1] = Convert.ToDouble(currentLineFields[4]);
                    force[0] = Convert.ToDouble(currentLineFields[9]);
                    force[1] = Convert.ToDouble(currentLineFields[10]);
                    double angle = Convert.ToDouble(currentLineFields[5]) * 360 / (2 * Math.PI);
                    translationalVelocity[0] = Convert.ToDouble(currentLineFields[6]);
                    translationalVelocity[1] = Convert.ToDouble(currentLineFields[7]);
                    double angularVelocity = Convert.ToDouble(currentLineFields[8]);
                    double torque = Convert.ToDouble(currentLineFields[11]);
                    currentParticle.Motion.InitializeParticlePositionAndAngle(new double[] { physicalData[3], physicalData[4] }, physicalData[5] * 360 / (2 * Math.PI), historyLength, t);
                    currentParticle.Motion.InitializeParticleVelocity(new double[] { physicalData[6], physicalData[7] }, physicalData[8], historyLength, t);
                    double currentParticleMass = currentParticle.Motion.Density * currentParticle.Volume;
                    double[] transAcc = new double[] { physicalData[9] / currentParticleMass, physicalData[10] / currentParticleMass };
                    double rotAcc = physicalData[11] / currentParticle.MomentOfInertia;
                    currentParticle.Motion.InitializeParticleAcceleration(transAcc, rotAcc, historyLength, t);
                }
            }
            return ParticleList;
        }
    }
}
