/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using FSI_Solver;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Application.FSI_Solver {
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain {

        /// <summary>
        /// Set the inital state of the simulation.
        /// </summary>
        protected override void SetInitial() {
            // Setup particles
            m_Particles = ((FSI_Control)this.Control).Particles;

            double phystimeInit = 0.0;
            UpdateLevelSetParticles(phystimeInit);

            // call base implementation
            base.SetInitial();
            int pCounter = 1;
            foreach (Particle p in m_Particles) {
                Console.WriteLine("===============================================================");
                Console.WriteLine("Particle properties, particle " + pCounter);
                Console.WriteLine("Type: " + p);
                Console.WriteLine("Density " + p.particleDensity);
                Console.WriteLine("Lengthscales #1 " + p.GetLengthScales()[0] + " #2 " + p.GetLengthScales()[1]);
            }
            if (!((FSI_Control)this.Control).pureDryCollisions) {
                Console.WriteLine("===============================================================");
                Console.WriteLine("Fluid properties: ");
                Console.WriteLine("Density: " + ((FSI_Control)this.Control).PhysicalParameters.rho_A);
                Console.WriteLine("Viscosity: " + ((FSI_Control)this.Control).PhysicalParameters.mu_A);
            }
            Console.WriteLine("===============================================================");
        }

        /// <summary>
        /// A list for all particles
        /// </summary>
        public IList<Particle> Particles {
            get {
                return m_Particles;
            }
        }
        List<Particle> m_Particles;

        /// <summary>
        /// The collision model
        /// </summary>
        private FSI_Solver.FSI_Control.CollisionModel CollisionModel {
            get {
                return ((FSI_Control)Control).collisionModel;
            }
        }

        /// <summary>
        /// Saves the physical data of all particles
        /// </summary>
        TextWriter logPhysicalDataParticles;

        bool CalculatedDampingTensors = false;
        readonly private FSI_Auxillary Auxillary = new FSI_Auxillary();

        //public static void MegaArschKakke2(DGField[] f) {
        //    int rank;
        //    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
        //    Tecplot.PlotFields(f, "MegaArschKakke-" + counter, 0.0, 0);
        //    counter++;
        //}

        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {
            //MultiphaseCellAgglomerator.Katastrophenplot = MegaArschKakke2;
            _Main(args, false, delegate () {
                var p = new FSI_SolverMain();
                return p;
            });
        }

        SinglePhaseField ParticleColor;
        SinglePhaseField LevelSetDistance;

        /// <summary>
        /// Create the colour and level set distance field. Necessary for the definition of the particle level set.
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();

            ParticleColor = new SinglePhaseField(new Basis(this.GridData, 0), "ParticleColor");
            m_RegisteredFields.Add(ParticleColor);
            m_IOFields.Add(ParticleColor);

            LevelSetDistance = new SinglePhaseField(new Basis(this.GridData, 0), "LevelSetDistance");
            m_RegisteredFields.Add(LevelSetDistance);
            m_IOFields.Add(LevelSetDistance);
        }

        /// <summary>
        /// Differentiate between splitting and coupled ansatz
        /// </summary>
        bool UseMovingMesh {
            get {
                switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                    case LevelSetHandling.Coupled_Once:
                    case LevelSetHandling.Coupled_Iterative:
                        return true;

                    case LevelSetHandling.LieSplitting:
                    case LevelSetHandling.StrangSplitting:
                    case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    case LevelSetHandling.None:
                        return false;

                    default:
                        throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)this.Control).Timestepper_LevelSetHandling);
                }
            }
        }

        /// <summary>
        /// Fully coupled LieSplitting?
        /// </summary>
        bool IsFullyCoupled {
            get {
                return ((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled;
            }
        }

        private double FluidViscosity => ((FSI_Control)Control).pureDryCollisions ? 0 : ((FSI_Control)Control).PhysicalParameters.mu_A;

        private double FluidDensity => ((FSI_Control)Control).pureDryCollisions ? 0 : ((FSI_Control)Control).PhysicalParameters.rho_A;

        private double HydrodynConvergenceCriterion => ((FSI_Control)Control).forceAndTorqueConvergenceCriterion;

        /// <summary>
        /// Creates Navier-Stokes and continuity eqution
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            if (IBM_Op != null)
                return;
            string[] CodNameSelected = new string[0];
            string[] DomNameSelected = new string[0];
            int spatialDim = GridData.SpatialDimension;

            // boundary conditions
            boundaryCondMap = new IncompressibleBoundaryCondMap(GridData, Control.BoundaryValues, PhysicsMode.Incompressible);

            // choose the operators
            NSEOperatorConfiguration IBM_Op_config = new NSEOperatorConfiguration {
                convection = Control.PhysicalParameters.IncludeConvection,
                continuity = true,
                Viscous = true,
                PressureGradient = true,
                Transport = true,
                CodBlocks = new bool[] { true, true },
                DomBlocks = new bool[] { true, true },
            };

            string[] CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, spatialDim)).Cat("div");
            string[] Params = ArrayTools.Cat(
                 VariableNames.Velocity0Vector(spatialDim),
                 VariableNames.Velocity0MeanVector(spatialDim));
            string[] DomName = ArrayTools.Cat(VariableNames.VelocityVector(spatialDim), VariableNames.Pressure);

            // selected part:
            if (IBM_Op_config.CodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, spatialDim));
            if (IBM_Op_config.CodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(spatialDim, 1));

            if (IBM_Op_config.DomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, spatialDim));
            if (IBM_Op_config.DomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(spatialDim, 1));

            IBM_Op = new XSpatialOperatorMk2(DomNameSelected, Params, CodNameSelected, (A, B, C) => HMForder, null);

            // Momentum equation
            // =============================
            // Convective part
            // =============================
            {
                if (IBM_Op_config.convection) {
                    for (int d = 0; d < spatialDim; d++) {

                        // The bulk
                        // -----------------------------
                        var comps = IBM_Op.EquationComponents[CodName[d]];
                        var convectionBulk = new Solution.NSECommon.LinearizedConvection(spatialDim, boundaryCondMap, d);
                        comps.Add(convectionBulk);

                        // Immersed boundary
                        // -----------------------------
                        if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                            var convectionAtIB = new BoSSS.Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, spatialDim, LsTrk,
                                this.Control.AdvancedDiscretizationOptions.LFFA, boundaryCondMap,
                                delegate (double[] X, double time) {
                                    throw new NotImplementedException("Currently not implemented for fixed motion");
                                },
                                this.Control.PhysicalParameters.rho_A,
                                UseMovingMesh);
                            comps.Add(convectionAtIB);
                        }
                        else {
                            var convectionAtIB = new BoSSS.Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, spatialDim, LsTrk,
                                Control.AdvancedDiscretizationOptions.LFFA, boundaryCondMap,
                                    delegate (double[] X, double time) {

                                        double[] result = new double[X.Length + 5];

                                        // Separating different boundary 
                                        // regions (for active particles)
                                        // -----------------------------
                                        foreach (Particle p in m_Particles) {
                                            // choose the correct particle at the current position
                                            bool containsParticle = m_Particles.Count == 1 ? true : p.Contains(X, GridData.iGeomCells.h_min.Min());
                                            p.CalculateRadialNormalVector(X, out double[] radialNormalVector);

                                            // active particles
                                            if (containsParticle && p.activeStress != 0) {
                                                double seperateBoundaryRegions = p.SeperateBoundaryRegions(X);
                                                result[0] = p.Motion.translationalVelocity[0][0];
                                                result[1] = p.Motion.translationalVelocity[0][1];
                                                result[2] = p.Motion.rotationalVelocity[0];
                                                result[3] = radialNormalVector[0];
                                                result[4] = radialNormalVector[1];
                                                result[5] = p.Motion.position[0].L2Distance(X);
                                                result[6] = -seperateBoundaryRegions;
                                            }
                                            // passive particles
                                            else if (containsParticle && p.activeStress == 0) {
                                                result[0] = p.Motion.translationalVelocity[0][0];
                                                result[1] = p.Motion.translationalVelocity[0][1];
                                                result[2] = p.Motion.rotationalVelocity[0];
                                                result[3] = radialNormalVector[0];
                                                result[4] = radialNormalVector[1];
                                                result[5] = p.Motion.position[0].L2Distance(X);
                                                result[6] = 0;
                                            }
                                        }
                                        return result;
                                    },
                                Control.PhysicalParameters.rho_A,
                                UseMovingMesh);
                            comps.Add(convectionAtIB);
                        }
                    }
                    U0MeanRequired = true;
                }
            }

            // Pressure part
            // =============================
            {
                if (IBM_Op_config.PressureGradient) {
                    for (int d = 0; d < spatialDim; d++) {
                        var comps = IBM_Op.EquationComponents[CodName[d]];

                        // The bulk
                        // -----------------------------
                        var pressureBulk = new PressureGradientLin_d(d, boundaryCondMap);
                        comps.Add(pressureBulk); // bulk component

                        // Immersed boundary
                        // -----------------------------
                        var pressureAtIB = new BoSSS.Solution.NSECommon.Operator.Pressure.ActivePressureAtIB(d, spatialDim, LsTrk);
                        comps.Add(pressureAtIB); // immersed boundary component

                        // if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                        if (this.Control.FixedStreamwisePeriodicBC) {
                            var presSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                            comps.Add(presSource);
                        }
                    }
                }
            }

            // Viscous part
            // =============================
            {
                if (IBM_Op_config.Viscous) {
                    for (int d = 0; d < spatialDim; d++) {
                        var comps = IBM_Op.EquationComponents[CodName[d]];
                        double penalty = this.Control.AdvancedDiscretizationOptions.PenaltySafety;

                        // The bulk
                        // -----------------------------
                        swipViscosity_Term1 viscousBulk = new swipViscosity_Term1(penalty, d, spatialDim, boundaryCondMap, ViscosityOption.ConstantViscosity, Control.PhysicalParameters.mu_A,
                            double.NaN, null);
                        comps.Add(viscousBulk);

                        // Immersed boundary
                        // -----------------------------
                        if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                            var viscousAtIB = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, spatialDim, LsTrk,
                                penalty, this.ComputePenaltyIB,
                                this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                                delegate (double[] X, double time) {
                                    throw new NotImplementedException("Currently not implemented for fixed motion");
                                });
                            comps.Add(viscousAtIB); // immersed boundary component
                        }
                        else {
                            Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB viscousAtIB = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, spatialDim, LsTrk,
                                penalty, this.ComputePenaltyIB,
                                this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                                delegate (double[] X, double time) {

                                    double[] result = new double[X.Length + 7];

                                    foreach (Particle p in m_Particles) {
                                        // which particle?
                                        bool containsParticle = m_Particles.Count == 1 ? true : p.Contains(X, GridData.iGeomCells.h_min.Min());
                                        p.CalculateRadialNormalVector(X, out double[] RadialNormalVector);

                                        // active particles
                                        if (containsParticle && p.activeStress != 0) {
                                            double seperateBoundaryRegions = p.SeperateBoundaryRegions(X);
                                            result[0] = p.Motion.translationalVelocity[0][0];
                                            result[1] = p.Motion.translationalVelocity[0][1];
                                            result[2] = p.Motion.rotationalVelocity[0];
                                            result[3] = RadialNormalVector[0];
                                            result[4] = RadialNormalVector[1];
                                            result[5] = p.Motion.position[0].L2Distance(X);
                                            result[6] = p.activeStress;
                                            result[7] = -seperateBoundaryRegions;
                                            result[8] = p.Motion.angle[0];
                                        }
                                        // passive particles
                                        else if (containsParticle && p.activeStress == 0) {
                                            result[0] = p.Motion.translationalVelocity[0][0];
                                            result[1] = p.Motion.translationalVelocity[0][1];
                                            result[2] = p.Motion.rotationalVelocity[0];
                                            result[3] = RadialNormalVector[0];
                                            result[4] = RadialNormalVector[1];
                                            result[5] = p.Motion.position[0].L2Distance(X);
                                            result[6] = 0;
                                            result[7] = 0;
                                            result[8] = p.Motion.angle[0];
                                        }
                                    }
                                    return result;
                                }
                             );
                            comps.Add(viscousAtIB); // immersed boundary component
                        }
                    }
                }
            }

            // Continuum equation
            // =============================
            {
                if (IBM_Op_config.continuity) {
                    for (int d = 0; d < spatialDim; d++) {
                        var src = new Divergence_DerivativeSource(d, spatialDim);
                        var flx = new Divergence_DerivativeSource_Flux(d, boundaryCondMap);
                        IBM_Op.EquationComponents["div"].Add(src);
                        IBM_Op.EquationComponents["div"].Add(flx);
                    }

                    if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                        var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(spatialDim, LsTrk, 1, delegate (double[] X, double time) {
                            throw new NotImplementedException("Currently not implemented for fixed motion");
                        });
                        IBM_Op.EquationComponents["div"].Add(divPen);  // immersed boundary component
                    }
                    else {
                        var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.ActiveDivergenceAtIB(spatialDim, LsTrk, 1,
                           delegate (double[] X, double time) {

                               double[] result = new double[X.Length + 4];

                               foreach (Particle p in m_Particles) {
                                   bool containsParticle = m_Particles.Count == 1 ? true : p.Contains(X, GridData.iGeomCells.h_min.Min());
                                   p.CalculateRadialNormalVector(X, out double[] RadialNormalVector);

                                   if (containsParticle) {
                                       result[0] = p.Motion.translationalVelocity[0][0];
                                       result[1] = p.Motion.translationalVelocity[0][1];
                                       result[2] = p.Motion.rotationalVelocity[0];
                                       result[3] = RadialNormalVector[0];
                                       result[4] = RadialNormalVector[1];
                                       result[5] = p.Motion.position[0].L2Distance(X);
                                       return result;
                                   }
                               }
                               return result;
                           });
                        IBM_Op.EquationComponents["div"].Add(divPen); // immersed boundary component 
                    }
                }
            }
            IBM_Op.Commit();

            // NSE or pure Stokes
            // =============================
            SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;
            if (this.Control.PhysicalParameters.IncludeConvection) {
                SpatialOp = SpatialOperatorType.Nonlinear;
            }

            // create timestepper, update level-set
            // =============================
            int bdfOrder;
            if (this.Control.Timestepper_Scheme == FSI_Control.TimesteppingScheme.CrankNicolson)
                bdfOrder = -1;
            else if (this.Control.Timestepper_Scheme == FSI_Control.TimesteppingScheme.ImplicitEuler)
                bdfOrder = 1;
            else if (this.Control.Timestepper_Scheme.ToString().StartsWith("BDF"))
                bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
            else
                throw new NotImplementedException("todo");

            MassMatrixShapeandDependence MassMatrixShape;
            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.Coupled_Once:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeDependent;
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    break;

                case LevelSetHandling.LieSplitting:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    if (!CalculatedDampingTensors) {
                        foreach (Particle p in m_Particles) {
                            if (p.Motion.m_AddedDampingCoefficient != -1) {
                                p.Motion.CalculateDampingTensor(p, LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                                Auxillary.ExchangeDampingTensors(m_Particles);
                            }
                        }
                    }
                    CalculatedDampingTensors = true;
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeDependent;
                    if (!CalculatedDampingTensors) {
                        foreach (Particle p in m_Particles) {
                            if (p.Motion.m_AddedDampingCoefficient != -1) {
                                p.Motion.CalculateDampingTensor(p, LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                                Auxillary.ExchangeDampingTensors(m_Particles);
                            }
                        }
                    }
                    CalculatedDampingTensors = true;
                    break;

                case LevelSetHandling.StrangSplitting:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeDependent;
                    break;

                case LevelSetHandling.None:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeDependent;
                    break;

                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)this.Control).Timestepper_LevelSetHandling);
            }

            m_BDF_Timestepper = new XdgBDFTimestepping(
                ArrayTools.Cat(this.Velocity, this.Pressure),
                ArrayTools.Cat(this.ResidualMomentum, this.ResidualContinuity),
                LsTrk,
                true,
                DelComputeOperatorMatrix, null, DelUpdateLevelset,
                bdfOrder,
                ((FSI_Control)this.Control).Timestepper_LevelSetHandling,
                MassMatrixShape,
                SpatialOp,
                MassScale,
                this.MultigridOperatorConfig, base.MultigridSequence,
                this.FluidSpecies, base.HMForder,
                this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold, true,
                this.Control.NonLinearSolver, this.Control.LinearSolver
                ) {
                m_ResLogger = base.ResLogger,
                m_ResidualNames = ArrayTools.Cat(this.ResidualMomentum.Select(f => f.Identification), this.ResidualContinuity.Identification),
                IterUnderrelax = ((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? ((FSI_Control)this.Control).LSunderrelax : 1.0,
                Config_LevelSetConvergenceCriterion = ((FSI_Control)this.Control).forceAndTorqueConvergenceCriterion,
                SessionPath = SessionPath,
                Timestepper_Init = Solution.Timestepping.TimeStepperInit.SingleInit
            };
        }

        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {

            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.None:
                    ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                    LevSet.ProjectField(Posfunction);
                    LsTrk.UpdateTracker();
                    break;

                case LevelSetHandling.Coupled_Once:
                    UpdateLevelSetParticles(phystime);
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    Console.WriteLine("WARNING: Coupled iterative solver is not tested!");
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                    CalculateHydrodynamicForces(m_Particles, dt);
                    UpdateLevelSetParticles(phystime);
                    foreach (Particle p in m_Particles) {
                        p.iteration_counter_P += 1;
                        p.forceAndTorque_convergence = ((FSI_Control)this.Control).forceAndTorqueConvergenceCriterion;
                    }
                    break;

                case LevelSetHandling.LieSplitting:
                    UpdateLevelSetParticles(phystime);
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    UpdateLevelSetParticles(phystime);
                    break;

                case LevelSetHandling.StrangSplitting:
                    UpdateLevelSetParticles(phystime);
                    break;

                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
            }

            // Forces and Torque residual
            /// <summary>
            /// Computes the Residual of the forces and torque acting from to fluid to the particle.
            /// </summary>
            double forces_PResidual;
            if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                double acc_force_P_x = 0;
                double acc_force_P_y = 0;
                double acc_torque_P = 0;
                double acc_force_P_x_old = 0;
                double acc_force_P_y_old = 0;
                double acc_torque_P_old = 0;
                double iterationCounter = -1;
                foreach (Particle p in m_Particles) {
                    // forces and torque of the previous iteration
                    acc_force_P_x_old += p.Motion.hydrodynamicForces[1][0];
                    acc_force_P_y_old += p.Motion.hydrodynamicForces[1][1];
                    acc_torque_P_old += p.Motion.hydrodynamicTorque[1];

                    // forces and torque of the current iteration
                    acc_force_P_x += p.Motion.hydrodynamicForces[0][0];
                    acc_force_P_y += p.Motion.hydrodynamicForces[0][1];
                    acc_torque_P += p.Motion.hydrodynamicTorque[0];
                    iterationCounter = p.iteration_counter_P;
                }
                // first iteration, to ensure at least two iterations per timestep
                if (iterationCounter == 0) {
                    forces_PResidual = 1;
                }
                // compute residual
                else {
                    forces_PResidual = Math.Sqrt((acc_force_P_x_old - acc_force_P_x).Pow2() + (acc_force_P_y_old - acc_force_P_y).Pow2() + (acc_torque_P_old - acc_torque_P).Pow2());
                }
                Console.WriteLine("Current forces_PResidual:   " + forces_PResidual);
            }
            // no iterative solver, no residual
            else {
                forces_PResidual = 0;
            }
            return forces_PResidual;
        }

        /// <summary>
        /// Array of all local cells with their specific color.
        /// </summary>
        private int[] cellColor = null;
        int[] globalParticleColor = null;
        /// <summary>
        /// Particle to Level-Set-Field 
        /// </summary>
        private void UpdateLevelSetParticles(double phystime) {
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate(LsTrk);
            // Step 1
            // Define an array with the respective cell colors
            // =======================================================
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            cellColor = cellColor == null ? InitializeColoring() : UpdateColoring();

            // Step 2
            // Delete the old level set
            // =======================================================
            DGLevSet.Current.Clear();

            // Step 3
            // Define level set per color
            // =======================================================
            CellMask allParticleMask = null;
            CellMask coloredCellMask = null;

            globalParticleColor = levelSetUpdate.DetermineGlobalParticleColor(GridData, cellColor, m_Particles);
            int[] _globalParticleColor = globalParticleColor.CloneAs();
            for (int p = 0; p < _globalParticleColor.Length; p++) {
                // Search for current colour on current process
                int currentColor = _globalParticleColor[p];
                bool processContainsCurrentColor = false;
                BitArray coloredCells = new BitArray(noOfLocalCells);
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (cellColor[j] == currentColor && currentColor != 0) {
                        processContainsCurrentColor = true;
                        coloredCells[j] = true;
                    }
                }

                if (processContainsCurrentColor) {
                    int[] particlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(_globalParticleColor, currentColor);
                    coloredCellMask = new CellMask(GridData, coloredCells);

                    // Save all colored cells of 
                    // any color in one cellmask
                    // -----------------------------
                    allParticleMask = allParticleMask == null ? coloredCellMask : allParticleMask.Union(coloredCellMask);

                    // Get particle level set
                    // -----------------------------
                    double phiComplete(double[] X, double t) {
                        // Generating the correct sign
                        double phi = Math.Pow(-1, particlesOfCurrentColor.Length - 1);
                        // Multiplication over all particle-level-sets within the current color
                        for (int pC = 0; pC < particlesOfCurrentColor.Length; pC++) {
                            Particle currentParticle = m_Particles[particlesOfCurrentColor[pC]];
                            phi *= currentParticle.Phi_P(X);
                            // Delete the particle within the current color from the particle color array
                            _globalParticleColor[particlesOfCurrentColor[pC]] = 0;
                        }
                        return phi;
                    }
                    // Set particle level set
                    // -----------------------------
                    SetLevelSet(phiComplete, coloredCellMask, phystime);
                }
            }

            // Step 4
            // Define level set of the remaining cells ("Fluid-Cells")
            // =======================================================
            double phiFluid(double[] X, double t) {
                return -1;
            }
            CellMask fluidCells = allParticleMask != null ? allParticleMask.Complement() : CellMask.GetFullMask(GridData);
            SetLevelSet(phiFluid, fluidCells, phystime);

            // Step 5
            // Smoothing
            // =======================================================
            PerformLevelSetSmoothing(allParticleMask, fluidCells, true);

            // Step 6
            // Update level set tracker
            // =======================================================
            LsTrk.UpdateTracker(__NearRegionWith: 2);
        }

        /// <summary>
        /// Set level set based on the function phi and the current cells
        /// </summary>
        private void SetLevelSet(Func<double[], double, double> phi, CellMask currentCells, double phystime) {
            ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(phi, phystime);
            DGLevSet.Current.Clear(currentCells);
            DGLevSet.Current.ProjectField(1.0, Function, new CellQuadratureScheme(UseDefaultFactories: true, domain: currentCells));
        }

        /// <summary>
        /// Update of <see cref="ParticleColor"/> and <see cref="LevelSetDistance"/>
        /// </summary>
        private int[] UpdateColoring() {
            // Step 1
            // Color all cells directlly related to the level set
            // ======================================================
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            ushort[] regionsCode = LsTrk.Regions.RegionsCode;
            int[] particleColor = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int[] particleColorExchange = particleColor.CloneAs();

            // No particles on current proc
            // -----------------------------
            if (m_Particles.Count == 0)
                return particleColor;

            particleColorExchange.MPIExchange(GridData);

            // Step 2
            // Color neighbour cells
            // =======================================================
            int neighbourSearchDepth = 2;
            for (int k = 0; k < neighbourSearchDepth; k++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (particleColorExchange[CellNeighbors[i]] != 0 && particleColorExchange[j] == 0) {
                            particleColor[j] = particleColorExchange[CellNeighbors[i]];
                        }
                    }
                }
                particleColorExchange = particleColor.CloneAs();
            }

            // Step 3
            // Communicate
            // =======================================================
            particleColorExchange.MPIExchange(GridData);

            // Step 4
            // Find neighbouring colours and recolour one of them
            // =======================================================
            // Find neighbouring cells with
            // different colours
            // -----------------------------
            int maxColor = particleColor.Max().MPIMax();
            int[,] colorToRecolorWith = new int[maxColor + 1, 2];
            for (int j = 0; j < noOfLocalCells; j++) {
                GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                for (int i = 0; i < CellNeighbors.Length; i++) {
                    if (particleColorExchange[CellNeighbors[i]] != particleColor[j] && particleColor[j] != 0 && particleColorExchange[CellNeighbors[i]] > 0) {
                        if (particleColorExchange[CellNeighbors[i]] < particleColor[j] || colorToRecolorWith[particleColor[j], 1] > particleColorExchange[CellNeighbors[i]]) {
                            colorToRecolorWith[particleColor[j], 0] = particleColor[j];
                            colorToRecolorWith[particleColor[j], 1] = particleColorExchange[CellNeighbors[i]];
                        }
                        if (particleColorExchange[CellNeighbors[i]] > particleColor[j]) {
                            if (colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 0] == 0 || colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 1] > particleColor[j]) {
                                colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 0] = particleColorExchange[CellNeighbors[i]];
                                colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 1] = particleColor[j];
                            }
                        }
                    }
                }
            }

            // Communicate
            // -----------------------------
            int[][,] GlobalColorToRecolorWith = colorToRecolorWith.MPIGatherO(0);
            GlobalColorToRecolorWith = GlobalColorToRecolorWith.MPIBroadcast(0);
            for (int m = 0; m < MPISize; m++) {
                for (int i = 0; i < maxColor + 1; i++) {
                    if (GlobalColorToRecolorWith[0][i, 1] == 0 || GlobalColorToRecolorWith[0][i, 1] > GlobalColorToRecolorWith[m][i, 1] && GlobalColorToRecolorWith[m][i, 1] != 0) {
                        GlobalColorToRecolorWith[0][i, 0] = GlobalColorToRecolorWith[m][i, 0];
                        GlobalColorToRecolorWith[0][i, 1] = GlobalColorToRecolorWith[m][i, 1];
                    }
                }
            }
            colorToRecolorWith = GlobalColorToRecolorWith[0];

            // Recolour
            // -----------------------------
            for (int i = maxColor; i > 0; i--) {
                if (colorToRecolorWith[i, 0] != 0) {
                    for (int j = 0; j < noOfLocalCells; j++) {
                        if (particleColor[j] == colorToRecolorWith[i, 0]) {
                            particleColor[j] = colorToRecolorWith[i, 1];
                        }
                    }
                }
            }

            for (int j = 0; j < noOfLocalCells; j++) {
                ParticleColor.SetMeanValue(j, particleColor[j]);
                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(regionsCode[j], 0));
            }
            return particleColor;
        }

        /// <summary>
        /// Initialization of <see cref="ParticleColor"/>  based on particle geometry
        /// </summary>
        private int[] InitializeColoring() {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cells = new int[J];
            List<int> coloredCells = new List<int>();
            for (int p = 0; p < m_Particles.Count; p++) {
                Particle currentParticle = m_Particles[p];
                double h_min = GridData.iGeomCells.h_min.Min() / 2;
                double h_max = GridData.iGeomCells.h_max.Max();
                for (int j = 0; j < J; j++) {
                    double[] center = GridData.iLogicalCells.GetCenter(j);
                    // Check for every cell whether their center is part of a particle or not (with tolerance sqrt(h_max^2+h_min^2))
                    if (currentParticle.Contains(center, h_min, h_min)) {
                        ParticleColor.SetMeanValue(j, p + 1);
                        coloredCells.Add(j);
                        cells[j] = p + 1;
                    }
                }
            }
            FixNeighbourColoring(cells);
            return cells;
        }

        /// <summary>
        /// Checks whether there are two different colours neighbouring each other. 
        /// </summary>
        /// <param name="coloredCells">
        /// All cells with their colour, uncoloured cells are set to zero.
        /// </param>
        private void FixNeighbourColoring(int[] coloredCells) {
            for (int i = 0; i < coloredCells.Length; i++) {
                if (coloredCells[i] != 0) {
                    GridData.GetCellNeighbours(i, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int j = 0; j < CellNeighbors.Length; j++) {
                        if (CellNeighbors[j] < coloredCells.Max() && coloredCells[i] != coloredCells[j] && coloredCells[CellNeighbors[j]] != 0) {
                            RecolorCells(coloredCells, coloredCells[i], coloredCells[CellNeighbors[j]]);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Recolours all cells with a specific colour.
        /// </summary>
        /// <param name="coloredCells">
        /// All cells with their colour, uncoloured cells are set to zero.
        /// </param>
        /// /// <param name="newColor">
        /// </param>
        /// /// <param name="oldColor">
        /// </param>
        private void RecolorCells(int[] coloredCells, int newColor, int oldColor) {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int i = 0; i < J; i++) {
                if (coloredCells[i] == oldColor) {
                    coloredCells[i] = newColor;
                    ParticleColor.SetMeanValue(i, newColor);
                }
            }
        }

        /// <summary>
        /// Calls the calculation of the hydrodyn. forces and torque in the particle.cs
        /// </summary>
        /// <param name="particles">
        /// </param>
        /// /// <param name="dt">
        /// </param>
        /// /// <param name="firstIteration">
        /// </param>
        private void CalculateHydrodynamicForces(List<Particle> particles, double dt, bool firstIteration = true) {
            // Note on MPI parallelization of particle solver:
            // ===============================================
            // - hydrodynamic forces are computed for each domain and added together;
            //   e.g. in the case of particles at MPI-boundaries each processor computes his part of the integral
            // - collisions are detected globally / collision forces are thus only computed on MPI rank 0
            // - finally, the forces on all particles are summed over all MPI processors and the result is stored on all processors (MPI_Allreduce)
            // - particle motion is computed for all particles simultaneously on all processors; every processor knows every particle
            // - since the particle solver is much cheaper than the flow solver, this "not-really parallel" approach may work up to a few hundreds of particles
            // ===============================================
            // Update forces
            // =============
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            for (int p = 0; p < particles.Count(); p++) {
                Particle currentParticle = particles[p];
                if (!((FSI_Control)Control).pureDryCollisions) { // && !currentParticle.skipForceIntegration) {
                    currentParticle.Motion.UpdateForcesAndTorque(Velocity, Pressure, LsTrk, currentParticle.CutCells_P(LsTrk), Control.PhysicalParameters.mu_A, Control.PhysicalParameters.rho_A, firstIteration, dt);

                }
                else {
                    currentParticle.Motion.hydrodynamicForces[0][0] = 0;
                    currentParticle.Motion.hydrodynamicForces[0][1] = 0;
                    currentParticle.Motion.hydrodynamicTorque[0] = 0;
                }
            }
            // MPISum over Forces moved to Particle.cs 
        }

        /// <summary>
        /// Calls the calculation of the Acceleration and the velocity
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="dt">
        /// The time step
        /// </param>
        /// <param name="FullyCoupled">
        /// Do you use FSI_Fully_Coupled?
        /// </param>
        /// <param name="IterationCounter">
        /// No of iterations
        /// </param>
        /// <param name="IncludeHydrodynamics"></param>
        internal void CalculateParticleVelocity(List<Particle> Particles, double dt, bool FullyCoupled, int IterationCounter, int TimestepInt, bool IncludeHydrodynamics = true) {
            foreach (Particle p in Particles) {
                p.iteration_counter_P = IterationCounter;
                if (IterationCounter == 0 && FullyCoupled) {
                    Console.WriteLine("Predicting forces for the next timestep...");
                    if (p.Motion.useAddedDamping) {
                        p.Motion.UpdateDampingTensors();
                    }
                    p.Motion.PredictForceAndTorque(TimestepInt);
                }
                //p.CalculateAcceleration(dt, FullyCoupled, IncludeHydrodynamics);
                //p.UpdateParticleVelocity(dt);
                p.Motion.UpdateParticleVelocity(dt);
            }
        }

        ///// <summary>
        ///// Calls the calculation of position
        ///// </summary>
        ///// <param name="Particles">
        ///// A list of all particles
        ///// </param>
        ///// <param name="dt">
        ///// The time step
        ///// </param>
        //internal void CalculateParticlePosition(List<Particle> Particles, double dt) {
        //    for (int p = 0; p < Particles.Count(); p++) {
        //        Particle currentParticle = Particles[p];
        //        currentParticle.CalculateParticlePosition(dt);
        //        currentParticle.CalculateParticleAngle(dt);
        //        currentParticle.collisionTimestep = 0;
        //    }
        //}

        /// <summary>
        /// runs solver one step?!
        /// </summary>
        /// <param name="TimestepInt">
        /// Timestep number
        /// </param>
        /// <param name="phystime">
        /// Physical time
        /// </param>
        /// <param name="dt">
        /// Timestep size
        /// </param>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {
                // init
                ResLogger.TimeStep = TimestepInt;
                dt = GetFixedTimestep();
                Console.WriteLine("Starting time-step " + TimestepInt + "...");
                // used later to check if there is exactly one push per timestep
                int OldPushCount = LsTrk.PushCount;

                // only particle motion & collisions, no flow solver
                // =================================================
                if (((FSI_Control)Control).pureDryCollisions) {
                    // in other branches, called by the BDF timestepper
                    // -------------------------------------------------
                    LsTrk.PushStacks();
                    DGLevSet.Push();

                    // physics
                    // -------------------------------------------------
                    foreach (Particle p in m_Particles) {
                        p.Motion.GetParticleDensity(p.particleDensity);
                        p.Motion.SaveHydrodynamicsOfPreviousTimestep();
                        p.Motion.SaveVelocityOfPreviousTimestep();
                        p.Motion.SavePositionAndAngleOfPreviousTimestep();
                    }
                    CalculateHydrodynamicForces(m_Particles, dt);
                    CalculateParticleVelocity(m_Particles, dt, IsFullyCoupled, 0, TimestepInt, false);
                    ResetCollisionState(m_Particles);
                    CalculateCollision(m_Particles, cellColor, dt);
                    foreach (Particle p in m_Particles) {
                        p.Motion.UpdateParticlePositionAndAngle(dt);
                        p.Motion.collisionTimestep = 0;
                    }
                    //CalculateParticlePosition(m_Particles, dt);
                    UpdateLevelSetParticles(phystime);
                    Auxillary.PrintResultToConsole(m_Particles, 0, 0, phystime, TimestepInt, out double MPIangularVelocity, out Test_Force);
                    // print and mpi check
                    // -------------------------------------------------
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                    //Auxillary.PrintResultToConsole(m_Particles, 0, 0, phystime, TimestepInt, out double MPIangularVelocity, out Test_Force);

                    // Save for NUnit Test
                    // -------------------------------------------------
                    base.QueryHandler.ValueQuery("C_Drag", 2 * Test_Force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                    base.QueryHandler.ValueQuery("C_Lift", 2 * Test_Force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
                    base.QueryHandler.ValueQuery("Angular_Velocity", MPIangularVelocity, true); // (TestCase FlowRotationalCoupling)

                }
                // particle motion & collisions plus flow solver
                // =================================================
                else {
                    if (((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.Coupled_Iterative) {
                        if (phystime == 0) { CreatePhysicalDataLogger(); }
                        int iterationCounter = 0;
                        double hydroDynForceTorqueResidual = 1e12;

                        //foreach (Particle p in m_Particles) {
                        //    p.Motion.SaveDataOfPreviousTimestep();
                        //}
                        foreach (Particle p in m_Particles) {
                            p.Motion.GetParticleDensity(p.particleDensity);
                            p.Motion.SaveHydrodynamicsOfPreviousTimestep();
                            p.Motion.activeStress = p.activeStress;
                        }
                        while (hydroDynForceTorqueResidual > HydrodynConvergenceCriterion) {
                            Auxillary.CheckForMaxIterations(iterationCounter, ((FSI_Control)Control).max_iterations_fully_coupled);
                            Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                            Auxillary.SaveOldParticleState(m_Particles, iterationCounter, ((FSI_Control)Control).forceAndTorqueConvergenceCriterion, IsFullyCoupled);
                            // actual physics
                            // -------------------------------------------------
                            if (iterationCounter != 0 || !IsFullyCoupled) { // in the first iteration of the fully coupled simulation the hydrodyn. forces are predicted by the particle.motion.cs
                                m_BDF_Timestepper.Solve(phystime, dt, false);
                                CalculateHydrodynamicForces(m_Particles, dt, false);
                            }
                            if (iterationCounter == 0) {
                                foreach (Particle p in m_Particles) {
                                    p.Motion.SaveVelocityOfPreviousTimestep();
                                }
                            }
                            CalculateParticleVelocity(m_Particles, dt, IsFullyCoupled, iterationCounter, TimestepInt);

                            // print
                            // -------------------------------------------------
                            if (IsFullyCoupled)
                                Auxillary.PrintResultToConsole(m_Particles, phystime, iterationCounter, out Test_Force);
                            else // not a fully coupled system? -> no iteration
                                break;

                            //residual
                            // -------------------------------------------------
                            hydroDynForceTorqueResidual = Auxillary.CalculateParticleResidual(m_Particles, ref iterationCounter);

                        }
                        foreach (Particle p in m_Particles) {
                            p.Motion.SavePositionAndAngleOfPreviousTimestep();
                        }
                        // collision
                        // -------------------------------------------------
                        ResetCollisionState(m_Particles);
                        CalculateCollision(m_Particles, cellColor, dt);

                        // particle position
                        // -------------------------------------------------
                        //CalculateParticlePosition(m_Particles, dt);
                        foreach (Particle p in m_Particles) {
                            p.Motion.UpdateParticlePositionAndAngle(dt);
                            p.Motion.collisionTimestep = 0;
                        }

                        // print
                        // -------------------------------------------------
                        Auxillary.PrintResultToConsole(m_Particles, FluidViscosity, FluidDensity, phystime, TimestepInt, out double Test_RotationalVelocity, out Test_Force);

                        // Save for NUnit Test
                        // -------------------------------------------------
                        base.QueryHandler.ValueQuery("C_Drag", 2 * Test_Force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                        base.QueryHandler.ValueQuery("C_Lift", 2 * Test_Force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
                        base.QueryHandler.ValueQuery("Angular_Velocity", Test_RotationalVelocity, true); // (TestCase FlowRotationalCoupling)

                        // level set tracker 
                        // -------------------------------------------------
                        if (IsFullyCoupled) {// in other branches, called by the BDF timestepper
                            LsTrk.IncreaseHistoryLength(1);
                            LsTrk.PushStacks();
                        }
                        LogPhysicalData(phystime);
                    }
                    else {// LevelSetHandling.Coupled_Iterative
                        foreach (Particle p in m_Particles) {
                            p.iteration_counter_P = -1;
                            p.forceAndTorque_convergence = ((FSI_Control)this.Control).forceAndTorqueConvergenceCriterion;
                        }
                        m_BDF_Timestepper.Solve(phystime, dt, false);
                    }
                }

                // finalize
                // ========
                if (LsTrk.PushCount - OldPushCount != 1) {
                    throw new ApplicationException("Illegal number of level-set push actions in time-step." + (LsTrk.PushCount - OldPushCount) + " It is important that LevelSetTracker.PushStacks() is called *exactly once per time-step*, at the beginning.");
                }
                ResLogger.NextTimestep(false);
                Console.WriteLine("Done with time-step " + TimestepInt + ".");
                return dt;
            }
        }

        /// <summary>
        /// Creates a log file for the physical data of the particles. Only active if a database is specified.
        /// </summary>
        private void CreatePhysicalDataLogger() {
            if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                logPhysicalDataParticles = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                string firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", "#Particle", "#Time", "Position X", "Position Y", "Angle", "Transl. Velocity X", "Transl. Velocity Y", "Rot. Velocity", "Force X", "Force Y", "Angular Momentum");
                logPhysicalDataParticles.WriteLine(firstline);
            }
        }

        /// <summary>
        /// Writes the physical data of the particles to a log file.
        /// </summary>
        /// <param name = phystime>
        /// </param>
        private void LogPhysicalData(double phystime) {
            if ((base.MPIRank == 0) && (logPhysicalDataParticles != null)) {
                for (int p = 0; p < m_Particles.Count(); p++) {
                    string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", p, phystime, m_Particles[p].Motion.position[0][0], m_Particles[p].Motion.position[0][1], m_Particles[p].Motion.angle[0], m_Particles[p].Motion.translationalVelocity[0][0], m_Particles[p].Motion.translationalVelocity[0][1], m_Particles[p].Motion.rotationalVelocity[0], m_Particles[p].Motion.hydrodynamicForces[0][0], m_Particles[p].Motion.hydrodynamicForces[0][1], m_Particles[p].Motion.hydrodynamicTorque[0]);
                    logPhysicalDataParticles.WriteLine(line);
                    logPhysicalDataParticles.Flush();
                }
            }
        }

        /// <summary>
        /// </summary>
        /// <param name = Particles>
        /// </param>
        internal void ResetCollisionState(List<Particle> Particles) {
            foreach (Particle p in Particles) {
                p.isCollided = false;
            }
        }

        // restart
        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {
            var tsi = new FSI_TimestepInfo(t, this.CurrentSessionInfo, timestepno, base.IOFields, m_Particles);
            SerialzeTester(tsi);
            return tsi;
        }

        /// <summary>
        /// Test the serialization of <see cref="FSI_TimestepInfo.Particles"/>
        /// </summary>
        private static void SerialzeTester(FSI_TimestepInfo b) {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Serialize
            };

            bool DebugSerialization = false;

            JsonReader GetJsonReader(Stream s) {
                if (DebugSerialization) {
                    return new JsonTextReader(new StreamReader(s));
                }
                else {
                    return new BsonReader(s);
                }
            }

            JsonWriter GetJsonWriter(Stream s) {
                if (DebugSerialization) {
                    return new JsonTextWriter(new StreamWriter(s));
                }
                else {
                    return new BsonWriter(s);
                }
            }


            byte[] buffer = null;
            using (var ms1 = new MemoryStream()) {
                using (var writer = GetJsonWriter(ms1)) {
                    formatter.Serialize(writer, b);
                    writer.Flush();
                    buffer = ms1.GetBuffer();
                    //writer.Close();
                }
            }

            FSI_TimestepInfo o;
            using (var ms2 = new MemoryStream(buffer)) {
                using (var reader = GetJsonReader(ms2)) {
                    o = formatter.Deserialize<FSI_TimestepInfo>(reader);
                    reader.Close();
                }
            }

            Debug.Assert(b.Particles.Length == o.Particles.Length);
            int L = b.Particles.Length;
            for (int l = 0; l < L; l++) { // loop over particles
                Debug.Assert(GenericBlas.L2Dist(b.Particles[l].Motion.position[0], o.Particles[l].Motion.position[0]) < 1e-13);
            }

        }

        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepNumber RestartFromDatabase(out double time) {

            // this sux, because the database API is totally fucked up
            var db = GetDatabase();
            Guid Rst_Tsid = base.GetRestartTimestepID();
            Guid Rst_SessionId = Control.RestartInfo.Item1;
            ISessionInfo session = db.Controller.GetSessionInfo(Rst_SessionId);

            var ArschInfo = ((DatabaseDriver)(base.DatabaseDriver)).LoadTimestepInfo<FSI_TimestepInfo>(Rst_Tsid, session, db);

            // init particles
            m_Particles = ArschInfo.Particles.ToList();
            UpdateLevelSetParticles(ArschInfo.PhysicalTime);

            // call base shit
            var R = base.RestartFromDatabase(out time);

            // return
            return R;
        }


        /// <summary>
        /// For restarting calculations, its important to reload old solutions if one uses a higher order method in time
        /// </summary>
        /// <param name="time"></param>
        /// <param name="timestep"></param>
        public override void PostRestart(double time, TimestepNumber timestep) {
            //var fsDriver = this.DatabaseDriver.FsDriver;
            //string pathToOldSessionDir = System.IO.Path.Combine(
            //    fsDriver.BasePath, "sessions", this.CurrentSessionInfo.RestartedFrom.ToString());
            //string pathToPhysicalData = System.IO.Path.Combine(pathToOldSessionDir,"PhysicalData.txt");
            //string[] records = File.ReadAllLines(pathToPhysicalData); 

            //string line1 = File.ReadLines(pathToPhysicalData).Skip(1).Take(1).First();
            //string line2 = File.ReadLines(pathToPhysicalData).Skip(2).Take(1).First();
            //string[] fields_line1 = line1.Split('\t');
            //string[] fields_line2 = line2.Split('\t');

            //Console.WriteLine("Line 1 " + fields_line1);

            //double dt = Convert.ToDouble(fields_line2[1]) - Convert.ToDouble(fields_line1[1]);

            //int idx_restartLine = Convert.ToInt32(time/dt + 1.0);
            //string restartLine = File.ReadLines(pathToPhysicalData).Skip(idx_restartLine-1).Take(1).First();
            //double[] values = Array.ConvertAll<string, double>(restartLine.Split('\t'), double.Parse);

            //if (time == values[1]+dt)
            //{
            //    Console.WriteLine("Restarting from time " + values[1]);
            //}

            //oldPosition[0] = values[7];
            //oldPosition[1] = values[8];
            //newTransVelocity[0] = values[4];
            //newTransVelocity[1] = values[5];
            //oldTransVelocity[0] = 0;
            //oldTransVelocity[1] = 0;
            //TransVelocityN2[0] = 0;
            //TransVelocityN2[1] = 0;
            //TransVelocityN3[0] = 0;
            //TransVelocityN3[1] = 0;
            //TransVelocityN4[0] = 0;
            //TransVelocityN4[1] = 0;
            //force[0] = values[2];
            //force[1] = values[3];

            //if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty))
            //{
            //    Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
            //    string firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#Timestep", "#Time", "DragForce", "LiftForce", "VelocityX", "VelocityY", "AngularVelocity", "xPosition", "yPosition", "ParticleRe");
            //    Log_DragAndLift.WriteLine(firstline);
            //    Log_DragAndLift.WriteLine(restartLine);
            //}
        }

        /// <summary>
        /// Triggers the collision detection, which triggers the calculation of the collisions
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="CellColor">
        /// All cells on the current process with their specific colour
        /// </param>
        /// <param name="dt">
        /// Timestep
        /// </param>
        private void CalculateCollision(List<Particle> Particles, int[] CellColor, double dt) {
            if (CollisionModel == FSI_Control.CollisionModel.NoCollisionModel)
                return;
            if (CollisionModel == FSI_Control.CollisionModel.RepulsiveForce)
                throw new NotImplementedException("Repulsive force model is currently unsupported, please use the momentum conservation model.");

            // Only particles with the same colour a close to each other, thus, we only test for collisions within those particles.
            // Determine colour.
            // =================================================
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate(LsTrk);
            //int[] globalParticleColor = levelSetUpdate.DetermineGlobalParticleColor(GridData, CellColor, Particles);
            int[] _GlobalParticleColor = globalParticleColor.CloneAs();
            for (int i = 0; i < _GlobalParticleColor.Length; i++) {
                int CurrentColor = _GlobalParticleColor[i];
                if (CurrentColor == 0)
                    continue;
                int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(_GlobalParticleColor, CurrentColor);

                // Multiple particles with the same colour, trigger collision detection
                // =================================================
                if (ParticlesOfCurrentColor.Length >= 1 && CurrentColor != 0) {
                    List<Particle> currentParticles = levelSetUpdate.GetParticleListOneColor(m_Particles, _GlobalParticleColor, CurrentColor);
                    FSI_Collision _Collision = new FSI_Collision(LsTrk, CurrentColor, FluidViscosity, FluidDensity, ((FSI_Control)Control).CoefficientOfRestitution, dt, LsTrk.GridDat.Cells.h_minGlobal);
                    _Collision.CalculateCollision(currentParticles, GridData, CellColor);
                }

                // Remove already investigated particles/colours from array
                // =================================================
                for (int j = 0; j < _GlobalParticleColor.Length; j++) {
                    if (_GlobalParticleColor[j] == CurrentColor)
                        _GlobalParticleColor[j] = 0;
                }
            }

            // Communicate
            // =================================================
            foreach (Particle p in m_Particles) {
                Collision_MPICommunication(p, MPISize);
            }
        }

        /// <summary>
        /// Ensures the communication between the processes after a collision. As collisions are triggered based on the (local) colouring of the cells only the owning process knows about them.
        /// Thus, it is necessary to inform the other processes about the collisions.
        /// </summary>
        /// <param name="currentParticle">
        /// The current particle.
        /// </param>
        /// <param name="MPISize">
        /// Number of mpi processes
        /// </param>
        private void Collision_MPICommunication(Particle currentParticle, int MPISize) {
            int noOfVars = 13;
            double[] isCollidedSend = new double[1];
            double[] sendSkipForceIntegration = new double[1];
            bool noCurrentCollision = true;

            // Did a collision take place on one of the processes?
            // ===================================================
            isCollidedSend[0] = currentParticle.isCollided ? 1 : 0;
            double[] isCollidedReceive = new double[MPISize];
            MPISendAndReceive(isCollidedSend, ref isCollidedReceive);

            for (int i = 0; i < isCollidedReceive.Length; i++) {
                // The particle is collided, thus, copy the data from
                // the owning process.
                // ===================================================
                if (isCollidedReceive[i] != 0) {
                    double[] dataSend = new double[noOfVars];
                    dataSend[0] = currentParticle.Motion.rotationalVelocity[0];
                    dataSend[1] = currentParticle.Motion.translationalVelocity[0][0];
                    dataSend[2] = currentParticle.Motion.translationalVelocity[0][1];
                    dataSend[3] = currentParticle.Motion.angle[0];
                    dataSend[4] = currentParticle.Motion.position[0][0];
                    dataSend[5] = currentParticle.Motion.position[0][1];
                    dataSend[6] = currentParticle.Motion.collisionTimestep;
                    dataSend[7] = currentParticle.Motion.rotationalVelocity[1];
                    dataSend[8] = currentParticle.Motion.translationalVelocity[1][0];
                    dataSend[9] = currentParticle.Motion.translationalVelocity[1][1];
                    dataSend[10] = currentParticle.Motion.angle[1];
                    dataSend[11] = currentParticle.Motion.position[1][0];
                    dataSend[12] = currentParticle.Motion.position[1][1];

                    double[] dataReceive = new double[noOfVars * MPISize];
                    MPISendAndReceive(dataSend, ref dataReceive);

                    currentParticle.Motion.rotationalVelocity[0] = dataReceive[0 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[0][0] = dataReceive[1 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[0][1] = dataReceive[2 + i * noOfVars];
                    currentParticle.Motion.angle[0] = dataReceive[3 + i * noOfVars];
                    currentParticle.Motion.position[0][0] = dataReceive[4 + i * noOfVars];
                    currentParticle.Motion.position[0][1] = dataReceive[5 + i * noOfVars];
                    currentParticle.Motion.collisionTimestep = dataReceive[6 + i * noOfVars];
                    currentParticle.Motion.rotationalVelocity[1] = dataReceive[7 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[1][0] = dataReceive[8 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[1][1] = dataReceive[9 + i * noOfVars];
                    currentParticle.Motion.angle[1] = dataReceive[10 + i * noOfVars];
                    currentParticle.Motion.position[1][0] = dataReceive[11 + i * noOfVars];
                    currentParticle.Motion.position[1][1] = dataReceive[12 + i * noOfVars];

                    currentParticle.isCollided = true;
                    noCurrentCollision = false;
                }
            }
            // nothing happend
            // ===================================================
            if (noCurrentCollision) {
                currentParticle.isCollided = false;
            }
        }

        /// <summary>
        /// MPI communication of an array.
        /// </summary>
        /// <param name="variableSend">
        /// The array to send.
        /// </param>
        /// <param name="variableReceive">
        /// An array of the data of all processes.
        /// </param>
        private void MPISendAndReceive(double[] variableSend, ref double[] variableReceive) {
            unsafe {
                fixed (double* pVariableSend = variableSend, pVariableReceive = variableReceive) {
                    csMPI.Raw.Allgather((IntPtr)pVariableSend, variableSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pVariableReceive, variableSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                }
            }
        }

        /// <summary>
        /// Adaptive mesh refinement
        /// </summary>
        /// <param name="TimestepNo">
        /// Currently unused.
        /// </param>
        /// <param name="newGrid">
        /// The adapted grid.
        /// </param>
        /// <param name="old2NewGrid">
        /// The correlation between old and new grid.
        /// </param>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            if (((FSI_Control)Control).AdaptiveMeshRefinement) {

                // Get the cells to refine and to coarse
                // ===================================================
                List<Tuple<int, CellMask>> AllCellsWithRefinementLevel = GetCellMaskWithRefinementLevels();
                bool AnyChange = GridRefinementController.ComputeGridChange((GridData)GridData, AllCellsWithRefinementLevel, out List<int> CellsToRefineList, out List<int[]> Coarsening);

                if (AnyChange) {
                    // Write stuff to console
                    // ===================================================
                    int[] consoleRefineCoarse = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                    int oldJ = this.GridData.CellPartitioning.TotalLength;
                    Console.WriteLine("       Refining " + consoleRefineCoarse[0] + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + consoleRefineCoarse[1] + " of " + oldJ + " cells");

                    // Adapt grid.
                    // ===================================================
                    newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);
                }
                else {
                    newGrid = null;
                    old2NewGrid = null;
                }
            }
            else {
                newGrid = null;
                old2NewGrid = null;
            }
        }

        private List<Tuple<int, CellMask>> GetCellMaskWithRefinementLevels() {
            int refinementLevel = ((FSI_Control)this.Control).RefinementLevel;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            double h_min = LsTrk.GridDat.Cells.h_minGlobal;
            double h_max = LsTrk.GridDat.Cells.h_maxGlobal;
            BitArray coarse = new BitArray(noOfLocalCells);
            BitArray fine = new BitArray(noOfLocalCells);
            foreach (Particle p in m_Particles) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (!fine[j])
                        fine[j] = p.Contains(new double[] { CellCenters[j, 0], CellCenters[j, 1] }, 4 * h_min);
                    if (LsTrk.Regions.IsSpeciesPresentInCell(LsTrk.GetSpeciesId("A"), j) && !coarse[j])
                        coarse[j] = p.Contains(new double[] { CellCenters[j, 0], CellCenters[j, 1] }, h_max / 2);
                }
            }

            CellMask coarseMask = new CellMask(GridData, coarse);
            CellMask fineMask = new CellMask(GridData, fine);
            int coarseRefinementLevel = refinementLevel > 2 ? refinementLevel / 2 : 1;
            if (refinementLevel - coarseRefinementLevel > coarseRefinementLevel)
                coarseRefinementLevel += 1;
            CellMask CutCells = LsTrk.Regions.GetCutCellMask();
            List<Tuple<int, CellMask>> AllCellsWithMaxRefineLevel = new List<Tuple<int, CellMask>>();
            AllCellsWithMaxRefineLevel.Add(new Tuple<int, CellMask>(refinementLevel, CutCells));
            AllCellsWithMaxRefineLevel.Add(new Tuple<int, CellMask>(refinementLevel, fineMask));
            AllCellsWithMaxRefineLevel.Add(new Tuple<int, CellMask>(coarseRefinementLevel, coarseMask));

            return AllCellsWithMaxRefineLevel;
        }
    }
}