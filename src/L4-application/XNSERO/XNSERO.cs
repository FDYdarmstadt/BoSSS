using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
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
    /// - Phoretic field is used by B. Liebchen (TRR146)
    /// </remarks>
    public class XNSERO : XNSE<XNSERO_Control> {
        /// <summary>
        /// Main 
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            //InitMPI();
            //TestProgram.TestRigidLevelSetProjection();
            //TestProgram.TestParticleParameter();
            //BoSSS.Application.XNSERO_Solver.TestProgram.TestRigidLevelSetProjection();
            //TestProgram.TestParticleInShearFlow_Phoretic();
            //Assert.IsFalse(true, "remove me");
            
            void KatastrophenPlot(DGField[] dGFields) {
                Tecplot.PlotFields(dGFields, "AgglomerationKatastrophe", 0.0, 3);
            }
            AgglomerationAlgorithm.Katastrophenplot = KatastrophenPlot;
            _Main(args, false, delegate () {
                var p = new XNSERO();
                return p;
            });
        }

        /// <summary>
        /// An array of all particles (rigid objects). Particles are only added at the initialization of the simulation. 
        /// </summary>
        public Particle[] Particles => Control.Particles;

        /// <summary>
        /// Spatial dimension
        /// </summary>
        private int SpatialDimension => GridData.SpatialDimension;

        /// <summary>
        /// Fluid viscosity
        /// </summary>
        private double[] FluidViscosity => new double[] { Control.PhysicalParameters.mu_A, Control.PhysicalParameters.mu_B };

        /// <summary>
        /// The position of the (horizontal and vertical) boundaries.
        /// </summary>
        /// <remarks>
        /// First entry: vertical [0] or horizontal [1]
        /// Second entry: left/lower wall [0] or right/upper wall [1]
        /// </remarks>
        private double[][] BoundaryCoordinates => Control.BoundaryPositionPerDimension;

        /// <summary>
        /// Array with two entries (2D). [0] true: x-Periodic, [1] true: y-Periodic
        /// </summary>
        private bool[] IsPeriodic => Control.BoundaryIsPeriodic;

        /// <summary>
        /// Grid length parameter used as tolerance measurement for particles.
        /// </summary>
        private double MaxGridLength => Control.MaxGridLength;

        /// <summary>
        /// Set this true in the control file to perform a simulation with two fluid species and the solid phase.
        /// </summary>
        private bool ContainsSecondFluidSpecies => Control.ContainsSecondFluidSpecies;

        /// <summary>
        /// Gravity vector, from control file.
        /// </summary>
        private Vector Gravity => Control.GetGravity();

        /// <summary>
        /// No position update, from control file.
        /// </summary>
        private bool AllParticlesFixed => Control.fixPosition;

        /// <summary>
        /// Coefficient of restitution for collisions.
        /// </summary>
        private double CoefficientOfRestitution => Control.CoefficientOfRestitution;
        
        /// <summary>
        /// Checks whether added damping tensors have been created. Only used if added damping is activated.
        /// </summary>
        private bool initAddedDamping = true;

        /// <summary>
        /// Additional helper methods. Checking the validity of calculated results and MPI consistency.
        /// </summary>
        private readonly Auxillary Auxillary = new();

        /// <summary>
        /// Saves the physical data of all particles
        /// </summary>
        private TextWriter LogParticleData;

        /// <summary>
        /// Prevents the creation of multiple log files
        /// </summary>
        private bool CreatedLogger = false;

        /// <summary>
        /// Particles are sorted into an R-tree to improve calculation of collisions
        /// </summary>
        private readonly RTree CollisionTree = new(2, 0.1);

        /// <summary>
        /// Prevents multiple initializations of <see cref="CollisionTree"/>
        /// </summary>
        private bool treeInitialized = false;

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
            var ParticleLevelSet = new Func<double[], double, double>[Particles.Length];
            for (int i = 0; i < ParticleLevelSet.Length; i++) {
                ParticleLevelSet[i] = Particles[i].LevelSetFunction;
            }
            return new RigidObjectLevelSet(ParticleLevelSet, MaxGridLength, null, Basis, Name);
        }

        /// <summary>
        /// Provides information about the evolution of the particle (rigid object) level set function to the level-set-updater.
        /// </summary>
        protected override RigidObjectLevelSetEvolver EvolveRigidLevelSet() {
            var ParticleLevelSet = new Func<double[], double, double>[Particles.Length];
            for (int i = 0; i < ParticleLevelSet.Length; i++) {
                ParticleLevelSet[i] = Particles[i].LevelSetFunction;
            }
            return new RigidObjectLevelSetEvolver(ParticleLevelSet, MaxGridLength);
        }

        /// <summary>
        /// Setup of the incompressible two-phase Navier-Stokes equation. If necessary adds the phoretic equations.
        /// </summary>
        /// <remarks>base: Navier Stokes, if(...): phoretic field</remarks>
        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            if (Control.UseAveragedEquations) {
                XNSE_OperatorConfiguration config = new(Control);

                // === mass equations === //
                DefineContinuityEquation(opFactory, config, D);

                // === momentum equations === //
                for (int d = 0; d < D; ++d) {
                    DefineMomentumEquation(opFactory, config, d, D);
                    // Add additional volume forces
                    if (config.isVolForce) {
                        var VolForceA = VolumeForce.CreateFrom("A", d, D, Control, Control.GetVolumeForce("A", d));
                        opFactory.AddParameter(VolForceA);
                        var VolForceB = VolumeForce.CreateFrom("B", d, D, Control, Control.GetVolumeForce("B", d));
                        opFactory.AddParameter(VolForceB);
                    }
                }
            } else {
                // NS 
                base.DefineSystem(D, opFactory, lsUpdater);
                // Phoretic field
                if (Control.UsePhoreticField) {
                    opFactory.AddEquation(new Equations.PhoreticFieldBulk());
                }
            }
        }

        /// <summary>
        /// Override this method to customize the assembly of the continuity equation
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="D">Spatial dimension (2 or 3)</param>
        protected override void DefineContinuityEquation(OperatorFactory opFactory, XNSE_OperatorConfiguration config, int D) {
            opFactory.AddEquation(new Continuity("A", config, D, boundaryMap));
            opFactory.AddEquation(new Continuity("B", config, D, boundaryMap));
            opFactory.AddEquation(new InterfaceContinuity("A", "B", config, D, config.isMatInt));
        }


        /// <summary>
        /// Definition of the boundary condition on the immersed boundary, i.e. at the surface of the particles, 
        /// <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// </summary>
        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            using (new FuncTrace()) {
                XNSE_OperatorConfiguration config = new(Control);

                for (int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new Equations.NSEROimmersedBoundary("A", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh, Control.UsePhoreticField, this.Particles));
                    opFactory.AddEquation(new Equations.NSEROimmersedBoundary("B", "C", 1, d, D, boundaryMap, LsTrk, config, config.isMovingMesh, Control.UsePhoreticField, this.Particles));
                }

                opFactory.AddEquation(new ImmersedBoundaryContinuity("A", "C", 1, config, D));
                opFactory.AddEquation(new ImmersedBoundaryContinuity("B", "C", 1, config, D));

                string[] fluidSpecies = CreateSpeciesArray(ContainsSecondFluidSpecies);

                RigidObjectLevelSetVelocity levelSetVelocity = new(VariableNames.LevelSetCGidx(1), Particles, FluidViscosity, fluidSpecies, Gravity, Control.dtFixed, MaxGridLength);
                opFactory.AddParameter(levelSetVelocity);

                Orientation OrientationVector = new(VariableNames.LevelSetCGidx(1), Particles, MaxGridLength);
                opFactory.AddParameter(OrientationVector);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCGidx(1), OrientationVector);

                if (Control.UsePhoreticField) {
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
            using (new FuncTrace()) {
                if (IsFluidInterface(iLevSet)) {
                    ILevelSetParameter levelSetVelocity = new LevelSetVelocity(VariableNames.LevelSetCG, SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                    return levelSetVelocity;

                } else if (IsParticleInterface(iLevSet)) {
                    SetPeriodicityToParticles();
                    string[] fluidSpecies = CreateSpeciesArray(ContainsSecondFluidSpecies);
                    ILevelSetParameter levelSetVelocity = new RigidObjectLevelSetVelocity(VariableNames.LevelSetCGidx(1), Particles, FluidViscosity, fluidSpecies, Gravity, Control.dtFixed, MaxGridLength);
                    return levelSetVelocity;

                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }
        }

        /// <summary>
        /// Creates an string array with either one species A or two species A and B.
        /// </summary>
        /// <returns></returns>
        private static string[] CreateSpeciesArray(bool ContainsSecondFluidSpecies) {
            string[] fluidSpecies;
            if (ContainsSecondFluidSpecies)
                fluidSpecies = new string[] { "A", "B" };
            else
                fluidSpecies = new string[] { "A" };
            return fluidSpecies;
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
            for (int d = 0; d < 2; d++) {
                if (IsPeriodic[d]) {
                    for (int p = 0; p < Particles.Length; p++) {
                        Particles[p].Motion.SetPeriodicBoundary(BoundaryCoordinates[d], d);
                    }
                }
            }
        }


        /// <summary>
        /// Update fluid variable fields and particle position and orientation angle.
        /// </summary>
        /// <param name="TimestepNo"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            dt = GetTimestep();
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt}");

            if (!CreatedLogger) {
                CreatedLogger = true;
                CreatePhysicalDataLogger();
            }

            InitializeParticlesNewTimestep(dt);
            if (TimestepNo - 1 % 10 == 0 || !treeInitialized) {
                CollisionTree.InitializeTree(Particles, dt);
                treeInitialized = true;
            } else {
                CollisionTree.UpdateTree(Particles, dt);
            }
            Auxillary.ParticleStateMPICheck(Particles, GridData, MPISize, TimestepNo);

            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

            CalculateCollision(Particles, dt);

            if (!AllParticlesFixed)
                CalculateParticlePositionAndAngle(Particles, dt);

            LogPhysicalData(phystime, TimestepNo);
            Console.WriteLine($"done with time step {TimestepNo}");
            return dt;
        }


        /// <summary>
        /// Safes old values for the velocity of the particles and updates added damping tensors (if used).
        /// </summary>
        private void InitializeParticlesNewTimestep(double dt) {
            foreach (Particle p in Particles) {
                p.Motion.SaveVelocityOfPreviousTimestep();
                p.LsTrk = LsTrk;
                if (p.Motion.UseAddedDamping) {
                    if (initAddedDamping) {
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
            using (new FuncTrace()) {
                int[][] potentialCollisionPartners = new int[Particles.Length][];
                for (int p = 0; p < Particles.Length; p++) {
                    Particles[p].IsCollided = false;
                }
                ParticleCollision Collision = new(MaxGridLength, CoefficientOfRestitution, dt, Control.WallPositionPerDimension, Control.BoundaryIsPeriodic, 0, DetermineOnlyOverlap);
                Collision.Calculate(Particles, CollisionTree);
            }
        }

        /// <summary>
        /// Calls the calculation of the particle position and orientation angle.
        /// </summary>
        /// <param name="Particles"></param>
        /// <param name="dt"></param>
        private static void CalculateParticlePositionAndAngle(Particle[] Particles, double dt) {
            foreach (Particle p in Particles) {
                p.Motion.UpdateParticlePositionAndAngle(dt);
            }
        }


        /// <summary>
        /// Creates a log file for the physical data of the particles. Only active if a database is specified.
        /// </summary>
        private void CreatePhysicalDataLogger() {
            if ((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                LogParticleData = DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                LogParticleData.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", "time-step", "particle", "time", "posX", "posY", "angle", "velX", "velY", "rot", "fX", "fY", "T"));
            }
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Writes the physical data of the particles to a log file.
        /// </summary>
        /// <param name = phystime>
        /// </param>
        private void LogPhysicalData(double phystime, int timestepNo) {
            using (new FuncTrace()) {
                if ((MPIRank == 0) && (LogParticleData != null)) {
                    for (int p = 0; p < Particles.Length; p++) {
                        LogParticleData.WriteLine($"{timestepNo},{p},{phystime},{Particles[p].Motion.GetPosition(0).x},{Particles[p].Motion.GetPosition(0).y},{Particles[p].Motion.GetAngle(0)},{Particles[p].Motion.GetTranslationalVelocity(0).x},{Particles[p].Motion.GetTranslationalVelocity(0).y},{Particles[p].Motion.GetRotationalVelocity(0)},{Particles[p].Motion.GetHydrodynamicForces(0).x},{Particles[p].Motion.GetHydrodynamicForces(0).y},{Particles[p].Motion.GetHydrodynamicTorque(0)}");
                        LogParticleData.Flush();
                    }
                }
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            base.AddMultigridConfigLevel(configsLevel, iLevel);

            if(Control.UsePhoreticField) {
                int pVel = VelocityDegree();

                var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pVel },
                    mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Phoretic) }
                };
                configsLevel.Add(configPres);
            }
        }
    }
}
