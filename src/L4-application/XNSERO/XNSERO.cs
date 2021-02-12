using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;


namespace BoSSS.Application.XNSERO_Solver {

    /// <summary>
    /// eXtended Navier Stokes Equation plus Rigid Object (XNSERO) Solver
    /// Fluid-Rigid Body solver, based on XDG.
    /// - incompressible flows.
    /// - solid immersed boundaries.
    /// - passive rigid bodies, e.g. sand
    /// - static rigid bodies
    /// - active rigid bodies, e.g. bacteria
    /// </summary>
    /// <remarks>
    /// Development history:
    /// - Current (jan2021) Maintainers: Deußen, Kummer
    /// - successor of the old IBM+FSI solver <see cref="IBM_SolverMain"/> and <see cref="FSI_SolverMain"/>, 
    ///   which were mainly used for PhD thesis of D. Krause and B. Deußen and TRR146.
    /// - see also: Extended discontinuous Galerkin methods for two-phase flows: the spatial discretization, F. Kummer, IJNME 109 (2), 2017. 
    /// </remarks>
    public class XNSERO : XNSE<XNSERO_Control> {

        /// <summary>
        /// Main 
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            //InitMPI();
            //DeleteOldPlotFiles();
            //TestProgram.TestParticleInShearFlow();
            //throw new Exception("remove me");


            
            void KatastrophenPlot(DGField[] dGFields) {
                Tecplot.PlotFields(dGFields, "AgglomerationKatastrophe", 0.0, 3);
            }
            MultiphaseCellAgglomerator.Katastrophenplot = KatastrophenPlot;
            _Main(args, false, delegate () {
                var p = new XNSERO();
                return p;
            });//*/
        }

        /// <summary>
        /// An array of all particles (rigid objects). Particles are only added at the initialization of the simulation. 
        /// </summary>
        [DataMember]
        public Particle[] Particles { get; private set; }

        /// <summary>
        /// Spatial dimension
        /// </summary>
        [DataMember]
        private int SpatialDimension;

        /// <summary>
        /// FluidViscosity Phase A
        /// </summary>
        [DataMember]
        private double[] FluidViscosity => new double[] { Control.PhysicalParameters.mu_A, Control.PhysicalParameters.mu_B };

        /// <summary>
        /// The position of the (horizontal and vertical) boundaries.
        /// </summary>
        /// <remarks>
        /// First entry: vertical [0] or horizontal [1]
        /// Second entry: left/lower wall [0] or right/upper wall [1]
        /// </remarks>
        [DataMember]
        private double[][] BoundaryCoordinates => Control.BoundaryPositionPerDimension;

        /// <summary>
        /// Array with two entries (2D). [0] true: x-Periodic, [1] true: y-Periodic
        /// </summary>
        [DataMember]
        private bool[] IsPeriodic => Control.BoundaryIsPeriodic;

        /// <summary>
        /// Grid length parameter used as tolerance measurement for particles.
        /// </summary>
        [DataMember]
        private double MaxGridLength => Control.MaxGridLength;

        [DataMember]
        private bool ContainsSecondFluidSpecies => Control.ContainsSecondFluidSpecies;

        [DataMember]
        private Vector Gravity => Control.GetGravity();

        [DataMember]
        private bool AllParticlesFixed => Control.fixPosition;

        /// <summary>
        /// Provides information about the particle (rigid object) level set function to the level-set-updater.
        /// </summary>
        /// <param name="Basis">
        /// A basis for the level-set.
        /// </param>
        /// <param name="Name">
        /// The name of the level set.
        /// </param>
        /// <remarks>
        /// Tested by <see cref="TestProgram.TestRigidLevelSetProjection"/>
        /// </remarks>
        protected override RigidObjectLevelSet SetRigidLevelSet(Basis Basis, string Name) {
            if(Particles.IsNullOrEmpty())
                Particles = Control.Particles.ToArray();
            CreatePhysicalDataLogger();
            Func<double[], double, double>[] ParticleLevelSet = new Func<double[], double, double>[Particles.Length];
            for(int i = 0; i < ParticleLevelSet.Length; i++) {
                ParticleLevelSet[i] = Particles[i].LevelSetFunction;
            }
            return new RigidObjectLevelSet(ParticleLevelSet, MaxGridLength, null, Basis, Name);
        }

        /// <summary>
        /// Provides information about the evolution of the particle (rigid object) level set function to the level-set-updater.
        /// </summary>
        protected override RigidObjectLevelSetEvolver EvolveRigidLevelSet() {

            Func<double[], double, double>[] ParticleLevelSet = new Func<double[], double, double>[Particles.Length];
            for(int i = 0; i < ParticleLevelSet.Length; i++) {
                ParticleLevelSet[i] = Particles[i].LevelSetFunction;
            }
            return new RigidObjectLevelSetEvolver(ParticleLevelSet, MaxGridLength);
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystem(D, opFactory, lsUpdater);

            if(Control.UsePhoreticField) {
                opFactory.AddEquation(new Equations.PhoreticFieldBulk());
            }
        }


        /// <summary>
        /// Definition of the boundary condition on the immersed boundary, i.e. at the surface of the particles, <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// Override to customize.
        /// </summary>
        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            using(new FuncTrace()) {
                XNSE_OperatorConfiguration config = new XNSE_OperatorConfiguration(this.Control);
                for(int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new Equations.NSEROimmersedBoundary("A", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                    opFactory.AddEquation(new Equations.NSEROimmersedBoundary("B", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                }

                opFactory.AddEquation(new ImmersedBoundaryContinuity("A", "C", 1, config, D, LsTrk));
                opFactory.AddEquation(new ImmersedBoundaryContinuity("B", "C", 1, config, D, LsTrk));

                opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));
                opFactory.AddParameter((ParameterS)GetLevelSetActiveStress(1));

                if(Control.UsePhoreticField) {
                    opFactory.AddEquation(new Equations.ImmersedBoundaryPhoreticField(LsTrk));
                }
            }
        }

        /// <summary>
        /// Provides the parameter field levelSetVelocity, depending whether <paramref name="iLevSet"/> is the fluid-fluid interface or fluid-particle interface.
        /// </summary>
        /// <param name="iLevSet">
        /// The level-set index.
        /// </param>
        /// <returns></returns>
        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            using(new FuncTrace()) {
                SpatialDimension = GridData.SpatialDimension;
                if(IsFluidInterface(iLevSet)) {
                    ILevelSetParameter levelSetVelocity = new LevelSetVelocity(VariableNames.LevelSetCG, SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                    return levelSetVelocity;

                } else if(IsParticleInterface(iLevSet)) {
                    SetPeriodicityToParticles();
                    string[] fluidSpecies;
                    if(ContainsSecondFluidSpecies)
                        fluidSpecies = new string[] { "A", "B" };
                    else
                        fluidSpecies = new string[] { "A" };
                    ILevelSetParameter levelSetVelocity = new RigidObjectLevelSetVelocity(VariableNames.LevelSetCGidx(1), Particles, FluidViscosity, fluidSpecies, Gravity, Control.dtFixed, MaxGridLength);
                    return levelSetVelocity;

                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }
        }

        /// <summary>
        /// Fluid-Fluid
        /// </summary>
        /// <param name="iLevSet"></param>
        /// <returns></returns>
        private static bool IsFluidInterface(int iLevSet) {
            return iLevSet == 0;
        }

        /// <summary>
        /// Fluid-Particle
        /// </summary>
        /// <param name="iLevSet"></param>
        /// <returns></returns>
        private static bool IsParticleInterface(int iLevSet) {
            return iLevSet == 1;
        }

        /// <summary>
        /// Provide information about periodic boundaries to the particles. Does nothing if no periodic boundaries are defined.
        /// </summary>
        private void SetPeriodicityToParticles() {
            for(int d = 0; d < 2; d++) {
                if(IsPeriodic[d]) {
                    for(int p = 0; p < Particles.Length; p++) {
                        Particles[p].Motion.SetPeriodicBoundary(BoundaryCoordinates[d], d);
                    }
                }
            }
        }

        /// <summary>
        /// Provides the active stress at the surface of active particles as a parameter field. 
        /// Active stress is used to define a boundary condition for the velocity gradient.
        /// </summary>
        /// <param name="iLevSet"></param>
        /// <returns></returns>
        protected virtual ILevelSetParameter GetLevelSetActiveStress(int iLevSet) {
            ILevelSetParameter levelSetVelocity = new ActiveStress(VariableNames.LevelSetCGidx(iLevSet), Particles, MaxGridLength);
            return levelSetVelocity;
        }

        /// <summary>
        /// Update fluid variable fields and particle position and orientation angle.
        /// </summary>
        /// <param name="TimestepNo"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            dt = GetFixedTimestep();
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt}");
            InitializeParticlesNewTimestep(dt);
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

            if(!AllParticlesFixed) {
                CalculateCollision(Particles, dt);
                CalculateParticlePositionAndAngle(Particles, dt);
            }

            LogPhysicalData(phystime, TimestepNo);
            Console.WriteLine($"done with time step {TimestepNo}");
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            Console.WriteLine("RunTime per time-step" + elapsedTime);
            return dt;
        }


        bool initAddedDamping = true;

        /// <summary>
        /// Safes old values for the velocity of the particles and updates added damping tensors (if used).
        /// </summary>
        private void InitializeParticlesNewTimestep(double dt) {
            foreach(Particle p in Particles) {
                p.Motion.SaveVelocityOfPreviousTimestep();
                if(p.Motion.UseAddedDamping) {
                    if(initAddedDamping) {
                        double fluidViscosity = (FluidViscosity[0] + FluidViscosity[1]) / 2;
                        p.Motion.CalculateDampingTensor(p, LsTrk, fluidViscosity, 1, dt);
                        p.Motion.ExchangeAddedDampingTensors();
                        initAddedDamping = false;
                    }
                    p.Motion.UpdateDampingTensors();
                }
            }
        }

        /// <summary>
        /// Triggers the collision detection, which triggers the calculation of the collisions. 
        /// Note on parallelization: All particle operations are carried out on all processes, hence no communication is necessary.
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="dt">
        /// Time-step
        /// </param>
        /// <param name="DetermineOnlyOverlap">
        /// Set true if you are only interested in overlapping particles and not the actual distance between different particles, e.g. as check for the initialization of static particles. 
        /// </param>
        private void CalculateCollision(Particle[] Particles, double dt, bool DetermineOnlyOverlap = false) {
            using(new FuncTrace()) {
                foreach(Particle p in Particles) {
                    p.IsCollided = false;
                }
                ParticleCollision Collision = new ParticleCollision(MaxGridLength, 1, dt, ((XNSERO_Control)Control).WallPositionPerDimension, ((XNSERO_Control)Control).BoundaryIsPeriodic, 0, DetermineOnlyOverlap);
                Collision.Calculate(Particles);
            }
        }

        /// <summary>
        /// Calls the calculation of the particle position and orientation angle.
        /// </summary>
        /// <param name="Particles"></param>
        /// <param name="dt"></param>
        private void CalculateParticlePositionAndAngle(Particle[] Particles, double dt) {
            foreach(Particle p in Particles) {
                p.Motion.UpdateParticlePositionAndAngle(dt);
            }
        }

        /// <summary>
        /// Saves the physical data of all particles
        /// </summary>
        private TextWriter logPhysicalDataParticles;

        /// <summary>
        /// Creates a log file for the physical data of the particles. Only active if a database is specified.
        /// </summary>
        private void CreatePhysicalDataLogger() {
            if((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                logPhysicalDataParticles = DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", "time-step", "particle", "time", "posX", "posY", "angle", "velX", "velY", "rot", "fX", "fY", "T"));
            }
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Writes the physical data of the particles to a log file.
        /// </summary>
        /// <param name = phystime>
        /// </param>
        private void LogPhysicalData(double phystime, int timestepNo) {
            using(new FuncTrace()) {
                if((MPIRank == 0) && (logPhysicalDataParticles != null)) {
                    for(int p = 0; p < Particles.Length; p++) {
                        logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", timestepNo, p, phystime, Particles[p].Motion.GetPosition(0)[0], Particles[p].Motion.GetPosition(0)[1], Particles[p].Motion.GetAngle(0), Particles[p].Motion.GetTranslationalVelocity(0)[0], Particles[p].Motion.GetTranslationalVelocity(0)[1], Particles[p].Motion.GetRotationalVelocity(0), Particles[p].Motion.GetHydrodynamicForces(0)[0], Particles[p].Motion.GetHydrodynamicForces(0)[1], Particles[p].Motion.GetHydrodynamicTorque(0)));
                        logPhysicalDataParticles.Flush();
                    }
                }
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            base.AddMultigridConfigLevel(configsLevel);

            if(Control.UsePhoreticField) {
                int pVel = VelocityDegree();

                var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pVel },
                    //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                    mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Phoretic) }
                };
                configsLevel.Add(configPres);
            }
        }
    }
}
