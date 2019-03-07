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

using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LevelSetTools.Smoothing;
using BoSSS.Foundation.SpecFEM;
using MPI.Wrappers;
using BoSSS.Foundation.IO;
using System.Diagnostics;
using System.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Grid.Classic;
using static BoSSS.Application.FSI_Solver.FSI_Control;

namespace BoSSS.Application.FSI_Solver {
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain {
        double calculatedDampingTensors;
        #region start
        // =============================
        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {

            //            System.Threading.Thread.Sleep(5000);
            //            BoSSS.Application.FSI_Solver.TestProgram.Init();
            //            BoSSS.Application.FSI_Solver.TestProgram.TestFlowRotationalCoupling();
            //Debug.Assert(false);

            _Main(args, false, delegate () {
                var p = new FSI_SolverMain();
                return p;
            });
        }
        #endregion
        #region field instantiation
        // =============================
        /// <summary>
        /// Curvature; DG-polynomial degree should be 2 times the polynomial degree of <see cref="LevSet"/>.
        /// </summary>
        [InstantiateFromControlFile("Curvature", "Curvature", IOListOption.ControlFileDetermined)]
        public SinglePhaseField Curvature;
        #endregion

        #region Create equations and solvers
        // =============================

        bool UseMovingMesh {
            get {
                switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                    case LevelSetHandling.Coupled_Once:
                    case LevelSetHandling.Coupled_Iterative:
                        return true;

                    case LevelSetHandling.LieSplitting:
                    case LevelSetHandling.StrangSplitting:
                    case LevelSetHandling.None:
                        return false;

                    default:
                        throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)this.Control).Timestepper_LevelSetHandling);
                }
            }
        }


        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            #region Misc
            if (IBM_Op != null)
                return;


            string[] CodNameSelected = new string[0];
            string[] DomNameSelected = new string[0];

            int D = this.GridData.SpatialDimension;


            BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);
            #endregion

            #region operator
            #region Config
            int degU = this.Velocity[0].Basis.Degree;
            var IBM_Op_config = new NSEOperatorConfiguration {
                convection = this.Control.PhysicalParameters.IncludeConvection,
                continuity = true,
                Viscous = true,
                PressureGradient = true,
                Transport = true,
                CodBlocks = new bool[] { true, true },
                DomBlocks = new bool[] { true, true },
            };

            var CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D)).Cat("div");
            var Params = ArrayTools.Cat(
                 VariableNames.Velocity0Vector(D),
                 VariableNames.Velocity0MeanVector(D));
            var DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

            // selected part:
            if (IBM_Op_config.CodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, D));
            if (IBM_Op_config.CodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D, 1));

            if (IBM_Op_config.DomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, D));
            if (IBM_Op_config.DomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D, 1));

            IBM_Op = new XSpatialOperator(DomNameSelected, Params, CodNameSelected,
                (A, B, C) => this.HMForder
                );
            #endregion

            #region Momentum equation
            // =============================

            #region Convective part
            // =============================
            if (IBM_Op_config.convection) {
                for (int d = 0; d < D; d++) {

                    var comps = IBM_Op.EquationComponents[CodName[d]];

                    //var ConvBulk = new Solution.XNSECommon.Operator.Convection.ConvectionInBulk_LLF(D, BcMap, d, this.Control.PhysicalParameters.rho_A, 0, this.Control.AdvancedDiscretizationOptions.LFFA, this.Control.AdvancedDiscretizationOptions.LFFB, LsTrk);
                    var ConvBulk = new Solution.NSECommon.LinearizedConvection(D, BcMap, d);
                    //IBM_Op.OnIntegratingBulk += ConvBulk.SetParameter;
                    comps.Add(ConvBulk); // bulk component

                    //var ConvIB = new BoSSS.Solution.XNSECommon.Operator.Convection.ConvectionAtIB(d, D, LsTrk, IBM_Op_config.dntParams.LFFA, BcMap, uIBM, wIBM);

                    if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                        var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, D, LsTrk,
                            this.Control.AdvancedDiscretizationOptions.LFFA, BcMap,
                            delegate (double[] X, double time) {
                                throw new NotImplementedException("Currently not implemented for fixed motion");
                                //return new double[] { 0.0, 0.0 };
                            },
                            this.Control.PhysicalParameters.rho_A,
                            UseMovingMesh);
                        comps.Add(ConvIB); // immersed boundary component
                    } else {
                        var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, D, LsTrk,
                            this.Control.AdvancedDiscretizationOptions.LFFA, BcMap,
                                delegate (double[] X, double time) {

                                    double[] result = new double[X.Length + 3];

                                    foreach (Particle p in m_Particles) {
                                        // Separating different boundary regions (for active particles)
                                        double cos_theta;
                                        // The posterior side of the particle (Neumann boundary)
                                        if (Math.Cos(p.angleAtIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.angleAtIteration[0]) * (X[1] - p.positionAtIteration[0][1]) < 1e-8)// && Math.Cos(p.particleAnglePerIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.particleAnglePerIteration[0]) * (X[1] - p.positionAtIteration[0][1]) > -0.25)
                                        {
                                            cos_theta = (Math.Cos(p.angleAtIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.angleAtIteration[0]) * (X[1] - p.positionAtIteration[0][1])) / (Math.Sqrt((X[0] - p.positionAtIteration[0][0]).Pow2() + (X[1] - p.positionAtIteration[0][1]).Pow2()));
                                        }
                                        // The anterior side of the particle (Dirichlet boundary)
                                        else {
                                            cos_theta = 0;
                                        }

                                        // which particle?
                                        bool containsParticle;
                                        if (m_Particles.Count == 1) {
                                            containsParticle = true;
                                        } else { containsParticle = p.Contains(X, LsTrk); }

                                        // active particles
                                        if (containsParticle && p.activeParticle == true) {
                                            result[0] = p.transVelocityAtIteration[0][0];
                                            result[1] = p.transVelocityAtIteration[0][1];
                                            result[2] = p.rotationalVelocityAtIteration[0];
                                            result[3] = p.positionAtIteration[0].L2Distance(X);
                                            result[4] = -cos_theta;
                                            return result;
                                        }

                                        // passive particles
                                        else if (containsParticle && p.activeParticle == false) {
                                            result[0] = p.transVelocityAtIteration[0][0];
                                            result[1] = p.transVelocityAtIteration[0][1];
                                            result[2] = p.rotationalVelocityAtIteration[0];
                                            result[3] = p.positionAtIteration[0].L2Distance(X);
                                            result[4] = 0;
                                            return result;
                                        }
                                    }
                                    return result;
                                },
                            this.Control.PhysicalParameters.rho_A,
                            UseMovingMesh);
                        comps.Add(ConvIB); // immersed boundary component
                    }
                }
                this.U0MeanRequired = true;
            }
            #endregion

            #region Pressure part
            // =============================
            if (IBM_Op_config.PressureGradient) {
                for (int d = 0; d < D; d++) {
                    var comps = IBM_Op.EquationComponents[CodName[d]];
                    //var pres = new Solution.XNSECommon.Operator.Pressure.PressureInBulk(d, BcMap, 1, 0);
                    var pres = new PressureGradientLin_d(d, BcMap);
                    //IBM_Op.OnIntegratingBulk += pres.SetParameter;
                    comps.Add(pres); // bulk component

                    var presLs = new BoSSS.Solution.NSECommon.Operator.Pressure.ActivePressureAtIB(d, D, LsTrk);//no changes necessary for impl. of active particles
                    comps.Add(presLs); // immersed boundary component

                    // if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                    if (this.Control.FixedStreamwisePeriodicBC) {
                        var presSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                        comps.Add(presSource);
                    }
                }
            }
            #endregion

            #region Viscous part
            // =============================
            if (IBM_Op_config.Viscous) {
                for (int d = 0; d < D; d++) {
                    var comps = IBM_Op.EquationComponents[CodName[d]];
                    double penalty = this.Control.AdvancedDiscretizationOptions.PenaltySafety;


                    var Visc = new swipViscosity_Term1(penalty, d, D, BcMap, ViscosityOption.ConstantViscosity, this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A, double.NaN, null);

                    comps.Add(Visc);


                    //var Visc = new Solution.XNSECommon.Operator.Viscosity.ViscosityInBulk_GradUTerm(penalty, 1.0, BcMap, d, D, this.Control.PhysicalParameters.mu_A, 0, ViscosityImplementation.H);
                    //IBM_Op.OnIntegratingBulk += Visc.SetParameter;
                    //comps.Add(Visc); // bulk component GradUTerm

                    //delegate (double p, int i, int j, double[] cell) { return ComputePenalty(p, i, j, cell); });
                    //delegate (double p, int i, int j, double[] cell) { return ComputePenalty(p, i, j, cell); });
                    //FSI_Op.OnIntegratingBulk += Visc.SetParameter;                

                    if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                        var ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, D, LsTrk,
                            penalty, this.ComputePenaltyIB,
                            this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                            delegate (double[] X, double time) {
                                throw new NotImplementedException("Currently not implemented for fixed motion");
                                //return new double[] { 0.0, 0.0 };
                            });
                        comps.Add(ViscLs); // immersed boundary component

                    } else {
                        var ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, D, LsTrk,
                            penalty, this.ComputePenaltyIB,
                            this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                            delegate (double[] X, double time) {

                                double[] result = new double[X.Length + 5];

                                foreach (Particle p in m_Particles) {
                                    // which particle?
                                    bool containsParticle;
                                    if (m_Particles.Count == 1) {
                                        containsParticle = true;
                                    } else { containsParticle = p.Contains(X, LsTrk); }

                                    // active particles
                                    if (containsParticle && p.activeParticle == true) {
                                        // Separating different boundary regions (for active particles)
                                        double cos_theta;
                                        // The posterior side of the particle (Neumann boundary)
                                        if (Math.Cos(p.angleAtIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.angleAtIteration[0]) * (X[1] - p.positionAtIteration[0][1]) < 1e-8)// && Math.Cos(p.particleAnglePerIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.particleAnglePerIteration[0]) * (X[1] - p.positionAtIteration[0][1]) > -0.25)
                                        {
                                            cos_theta = (Math.Cos(p.angleAtIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.angleAtIteration[0]) * (X[1] - p.positionAtIteration[0][1])) / (Math.Sqrt((X[0] - p.positionAtIteration[0][0]).Pow2() + (X[1] - p.positionAtIteration[0][1]).Pow2()));
                                        }
                                        // The anterior side of the particle (Dirichlet boundary)
                                        else {
                                            cos_theta = 0;
                                        }
                                        result[0] = p.transVelocityAtIteration[0][0];
                                        result[1] = p.transVelocityAtIteration[0][1];
                                        result[2] = p.rotationalVelocityAtIteration[0];
                                        if (p is Particle_Sphere) {
                                            result[3] = ((Particle_Sphere)p).radius_P;
                                        } else {
                                            result[3] = p.positionAtIteration[0].L2Distance(X);
                                        }
                                        result[4] = p.active_stress_P;
                                        result[5] = -cos_theta;
                                        result[6] = p.angleAtIteration[0];
                                    }

                                    // passive particles
                                    else if (containsParticle && p.activeParticle == false) {
                                        result[0] = p.transVelocityAtIteration[0][0];
                                        result[1] = p.transVelocityAtIteration[0][1];
                                        result[2] = p.rotationalVelocityAtIteration[0];
                                        if (p is Particle_Sphere) {
                                            result[3] = ((Particle_Sphere)p).radius_P;
                                        } else {
                                            result[3] = p.positionAtIteration[0].L2Distance(X);
                                        }
                                        result[4] = 0;
                                        result[5] = 0;
                                        result[6] = p.angleAtIteration[0];
                                    }
                                }
                                return result;
                            }
                         );
                        comps.Add(ViscLs); // immersed boundary component
                    }
                }
            }
            #endregion
            #endregion

            #region Continuum equation
            // ==================
            if (IBM_Op_config.continuity) {
                for (int d = 0; d < D; d++) {
                    //var src = new Solution.XNSECommon.Operator.Continuity.DivergenceInBulk_Volume(d, D, 1, 0, 1, false);
                    //IBM_Op.OnIntegratingBulk += src.SetParameter;
                    //var flx = new Solution.XNSECommon.Operator.Continuity.DivergenceInBulk_Edge(d, BcMap, 1, 0, 1, false);
                    //IBM_Op.OnIntegratingBulk += flx.SetParameter;
                    var src = new Divergence_DerivativeSource(d, D);
                    //IBM_Op.OnIntegratingBulk += src.SetParameter;
                    var flx = new Divergence_DerivativeSource_Flux(d, BcMap);
                    IBM_Op.EquationComponents["div"].Add(src);
                    IBM_Op.EquationComponents["div"].Add(flx);

                }

                if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                    var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, LsTrk, 1, delegate (double[] X, double time) {
                        throw new NotImplementedException("Currently not implemented for fixed motion");
                        //return new double[] { 0.0, 0.0 };
                    });
                    IBM_Op.EquationComponents["div"].Add(divPen);  // immersed boundary component
                } else {
                    var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, LsTrk, 1,
                       delegate (double[] X, double time) {

                           double[] result = new double[X.Length + 2];

                           foreach (Particle p in m_Particles) {
                               bool containsParticle;
                               if (m_Particles.Count == 1) {
                                   containsParticle = true;
                               } else { containsParticle = p.Contains(X, LsTrk); }
                               if (containsParticle) {
                                   result[0] = p.transVelocityAtIteration[0][0];
                                   result[1] = p.transVelocityAtIteration[0][1];
                                   result[2] = p.rotationalVelocityAtIteration[0];
                                   if (p is Particle_Sphere) {
                                       result[3] = ((Particle_Sphere)p).radius_P;
                                   } else {
                                       result[3] = p.positionAtIteration[0].L2Distance(X);
                                   }
                                   return result;
                               }
                           }
                           return result;
                       });
                    IBM_Op.EquationComponents["div"].Add(divPen); // immersed boundary component 
                }
            }
            #endregion
            IBM_Op.Commit();
            #endregion

            #region NSE or pure Stokes
            // ------------------
            SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;
            if (this.Control.PhysicalParameters.IncludeConvection) {
                SpatialOp = SpatialOperatorType.Nonlinear;
            }
            #endregion

            #region create timestepper, update level-set
            // ------------------
            int bdfOrder;
            if (this.Control.Timestepper_Scheme == FSI_Control.TimesteppingScheme.CrankNicolson)
                bdfOrder = -1;
            //else if (this.Control.Timestepper_Scheme == IBM_Control.TimesteppingScheme.ExplicitEuler)
            //    bdfOrder = 0;
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
                    if (calculatedDampingTensors == 0)
                    {
                        foreach (Particle p in m_Particles)
                        {
                            p.CalculateDampingTensors(LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                        }
                        calculatedDampingTensors = 1;
                    }
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
                DelComputeOperatorMatrix, DelUpdateLevelset,
                bdfOrder,
                ((FSI_Control)this.Control).Timestepper_LevelSetHandling,
                MassMatrixShape,
                SpatialOp,
                MassScale,
                this.MultigridOperatorConfig, base.MultigridSequence,
                this.FluidSpecies, base.HMForder,
                this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold, true,
                this.Control.NonLinearSolver, this.Control.LinearSolver
                );
            m_BDF_Timestepper.m_ResLogger = base.ResLogger;
            m_BDF_Timestepper.m_ResidualNames = ArrayTools.Cat(this.ResidualMomentum.Select(f => f.Identification), this.ResidualContinuity.Identification);
            m_BDF_Timestepper.IterUnderrelax = ((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? ((FSI_Control)this.Control).LSunderrelax : 1.0;
            m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
            m_BDF_Timestepper.SessionPath = SessionPath;
            m_BDF_Timestepper.Timestepper_Init = Solution.Timestepping.TimeStepperInit.SingleInit;

        }

        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {
            #region Level-set handling
            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.None:
                    ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                    newTransVelocity[0] = (((FSI_Control)this.Control).transVelocityFunc[0])(phystime);
                    newTransVelocity[1] = (((FSI_Control)this.Control).transVelocityFunc[1])(phystime);
                    LevSet.ProjectField(Posfunction);
                    LsTrk.UpdateTracker();
                    break;

                case LevelSetHandling.Coupled_Once:
                    UpdateLevelSetParticles(dt);
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    UpdateForcesAndTorque(dt, phystime);
                    UpdateLevelSetParticles(dt);
                    foreach (Particle p in m_Particles) {
                        p.iteration_counter_P += 1;
                        p.forceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                    }
                    break;

                case LevelSetHandling.LieSplitting:
                    UpdateLevelSetParticles(dt);
                    break;

                case LevelSetHandling.StrangSplitting:
                    UpdateLevelSetParticles(dt);
                    break;

                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
            }
            #endregion
            #region Forces and Torque residual
            /// <summary>
            /// Computes the Residual of the forces and torque acting from to fluid to the particle.
            /// </summary>
            double forces_PResidual = 0;
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
                    acc_force_P_x_old += p.hydrodynForcesAtIteration[1][0];
                    acc_force_P_y_old += p.hydrodynForcesAtIteration[1][1];
                    acc_torque_P_old += p.hydrodynTorqueAtIteration[1];

                    // forces and torque of the current iteration
                    acc_force_P_x += p.hydrodynForcesAtIteration[0][0];
                    acc_force_P_y += p.hydrodynForcesAtIteration[0][1];
                    acc_torque_P += p.hydrodynTorqueAtIteration[0];
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
            #endregion
        }



        void UpdateLevelSetParticles(double dt) {
            // Call update methods
            foreach (Particle p in m_Particles) {
                p.ResetParticlePosition();
                p.UpdateDampingTensors();
                if (p.iteration_counter_P == 0 && ((FSI_Control)this.Control).splitting_fully_coupled == true)
                {
                    p.PredictTranslationalAccelaration();
                    p.PredictAngularAccelaration();
                    p.PredictTranslationalVelocity();
                    p.PredictAngularVelocity();
                }
                else
                {
                    p.UpdateAngularVelocity(dt, ((FSI_Control)this.Control).includeRotation);
                    if (((FSI_Control)this.Control).includeTranslation == true)
                    {
                        p.CalculateTranslationalAccelaration(dt, this.Control.PhysicalParameters.rho_A);
                        p.CalculateTranslationalVelocity(dt, this.Control.PhysicalParameters.rho_A);
                    }
                    p.ComputeParticleRe(this.Control.PhysicalParameters.mu_A);
                    p.UpdateParticlePosition(dt, this.Control.PhysicalParameters.rho_A);
                }
                
            }

            // Update phi complete
            Func<double[], double, double> phiComplete = delegate (double[] X, double t) {
                int exp = m_Particles.Count - 1;
                double ret = Math.Pow(-1, exp);
                for (int i = 0; i < m_Particles.Count; i++) {
                    ret *= m_Particles[i].phi_P(X, t);
                }
                return ret;
            };

            // Vectorize
            ScalarFunction function = NonVectorizedScalarFunction.Vectorize(phiComplete, hack_phystime);
            LevSet.ProjectField(function);
            DGLevSet.Current.ProjectField(function);
            LsTrk.UpdateTracker(__NearRegionWith: 2);
        }
        #endregion

        void UpdateForcesAndTorque(double dt, double phystime) {
            foreach (Particle p in m_Particles) {
                if (!((FSI_Control)this.Control).pureDryCollisions) {
                    p.UpdateForcesAndTorque(Velocity, Pressure, LsTrk, this.Control.PhysicalParameters.mu_A);
                }
                WallCollisionForces(p, LsTrk.GridDat.Cells.h_minGlobal);
            }

            double[] totalMomentum = new double[2] { 0, 0 };
            double[] totalKE = new double[3] { 0, 0, 0 };
            double xPos;
            double yPos;
            double ang;

            foreach (Particle p in m_Particles) {
                totalMomentum[0] += p.Mass_P * p.transVelocityAtIteration[0][0];
                totalMomentum[1] += p.Mass_P * p.transVelocityAtIteration[0][1];
                totalKE[0] += 0.5 * p.Mass_P * p.transVelocityAtIteration[0][0].Pow2();
                totalKE[1] += 0.5 * p.Mass_P * p.transVelocityAtIteration[0][1].Pow2();
                totalKE[2] += 0.5 * p.MomentOfInertia_P * p.rotationalVelocityAtIteration[0].Pow2();
            }

            Console.WriteLine("Total-Momentum in System:  " + Math.Sqrt(totalMomentum[0].Pow2() + totalMomentum[1].Pow2()));
            Console.WriteLine("Total-KineticEnergy in System:  " + (totalKE[0] + totalKE[1] + totalKE[2]));

            totalMomentumOld = Math.Sqrt(totalMomentum[0].Pow2() + totalMomentum[1].Pow2());

            if (m_Particles.Count > 1)
                UpdateCollisionForces(m_Particles, LsTrk.GridDat.Cells.h_minGlobal);

            force = m_Particles[0].hydrodynForcesAtIteration[0];
            torque = m_Particles[0].hydrodynTorqueAtIteration[0];

            xPos = m_Particles[0].positionAtIteration[0][0];
            yPos = m_Particles[0].positionAtIteration[0][1];
            ang = m_Particles[0].angleAtIteration[0];


            MPItransVelocity = m_Particles[0].transVelocityAtIteration[0];
            MPIangularVelocity = m_Particles[0].rotationalVelocityAtIteration[0];


            Console.WriteLine(newPosition[1].MPIMax());

            if ((base.MPIRank == 0) && (Log_DragAndLift != null)) {
                double drag = force[0];
                double lift = force[1];
                //string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", TimestepNo, phystime, m_Particles[0].positionAtIteration[0][0], m_Particles[0].positionAtIteration[0][1], m_Particles[0].particleAnglePerIteration[0], m_Particles[0].transVelocityAtIteration[0][0], m_Particles[0].transVelocityAtIteration[0][1], 0.0, (totalKE[0] + totalKE[1] + totalKE[2]), Math.Sqrt(totalMomentum[0].Pow2() + totalMomentum[1].Pow2()));
                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", phystime, m_Particles[0].positionAtIteration[0][0], m_Particles[0].positionAtIteration[0][1], m_Particles[0].angleAtIteration[0], m_Particles[0].transVelocityAtIteration[0][0], m_Particles[0].transVelocityAtIteration[0][1], 0.0, (totalKE[0] + totalKE[1] + totalKE[2]), Math.Sqrt(totalMomentum[0].Pow2() + totalMomentum[1].Pow2()));
                Log_DragAndLift.WriteLine(line);
                Log_DragAndLift.Flush();
            }

            oldAngularVelocity = newAngularVelocity;

            // Save for NUnit Test
            base.QueryHandler.ValueQuery("C_Drag", 2 * force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
            base.QueryHandler.ValueQuery("C_Lift", 2 * force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
            base.QueryHandler.ValueQuery("Angular_Velocity", MPIangularVelocity, true); // (TestCase FlowRotationalCoupling)


            Console.WriteLine("Drag Force:   {0}", force[0]);
            Console.WriteLine("Lift Force:   {0}", force[1]);
            Console.WriteLine("Torqe:   {0}", torque);
            Console.WriteLine("Transl VelocityX:   {0}", MPItransVelocity[0]);
            Console.WriteLine("Transl VelocityY:   {0}", MPItransVelocity[1]);
            Console.WriteLine("Angular Velocity:   {0}", MPIangularVelocity);
            Console.WriteLine("X-position:   {0}", xPos);
            Console.WriteLine("Y-position:   {0}", yPos);
            Console.WriteLine("Angle:   {0}", ang);
            Console.WriteLine();
            Console.WriteLine("=======================================================");
            Console.WriteLine();
        }
        #endregion

        #region Run solver one step
        /// <summary>
        /// Variables for FSI coupling
        /// </summary>
        double oldAngularVelocity,
            newAngularVelocity = 0.0, MPIangularVelocity;
        double[] TransVelocityN4 = new double[2];
        double[] TransVelocityN3 = new double[2];
        double[] TransVelocityN2 = new double[2];
        double[] oldTransVelocity = new double[2];
        double[] newTransVelocity = new double[2];
        double[] oldPosition = new double[2];
        double[] newPosition = new double[2];
        double[] oldforce = new double[2];
        double[] MPItransVelocity = new double[2];
        double[] MPIpos = new double[2];
        double totalMomentumOld = 0;

        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {
                
                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;

                base.ResLogger.TimeStep = TimestepInt;

                hack_phystime = phystime;
                dt = base.GetFixedTimestep();

                if (((FSI_Control)this.Control).pureDryCollisions) {
                    UpdateLevelSetParticles(dt);
                } else {
                    if (triggerOnlyCollisionProcedure) {
                        UpdateLevelSetParticles(dt);
                        triggerOnlyCollisionProcedure = false;
                        if (phystime == 0) {
                            if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                                Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                                string firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#Timestep", "#Time", "P0_PosX", "P0_PosY", "P0_angle", "P0_VelX", "P0_VelY", "xPosition", "TotalKineticEnergy", "TotalMomentum");
                                Log_DragAndLift.WriteLine(firstline);
                                if (m_Particles.Count > 1) {
                                    Log_DragAndLift_P1 = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData_P1", CurrentSessionInfo.ID);
                                    string firstline_P1 = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#Timestep", "#Time", "P1_PosX", "P1_PosY", "P1_angle", "P1_VelX", "P1_VelY", "xPosition", "TotalKineticEnergy", "TotalMomentum");
                                    Log_DragAndLift_P1.WriteLine(firstline_P1);
                                }
                            }


                        }

                        return dt;
                    } else if (((FSI_Control)this.Control).Timestepper_LevelSetHandling != LevelSetHandling.Coupled_Iterative) {
                        int iteration_counter = 0;
                        for (double posResidual_splitting = 1e12; posResidual_splitting > ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;)// && iteration_counter <= (this.Control).max_iterations_fully_coupled;)
                        {
                            foreach (Particle p in m_Particles) {
                                p.iteration_counter_P = iteration_counter;
                                p.forceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                            }
                            m_BDF_Timestepper.Solve(phystime, dt, false);
                            #region Get Drag and Lift Coefficiant
                            UpdateForcesAndTorque(dt, phystime);
                            double acc = 0;
                            foreach (Particle p in m_Particles) {
                                acc += (p.hydrodynForcesAtIteration[0][0] - p.hydrodynForcesAtIteration[1][0]).Pow2() + (p.hydrodynForcesAtIteration[0][1] - p.hydrodynForcesAtIteration[1][1]).Pow2() + (p.hydrodynTorqueAtIteration[0] - p.hydrodynTorqueAtIteration[1]).Pow2();
                            }
                            posResidual_splitting = Math.Sqrt(acc);
                            Console.WriteLine("Fully coupled system, number of iterations:  " + iteration_counter);
                            Console.WriteLine("Forces and torque residual: " + posResidual_splitting);
                            Console.WriteLine();
                            iteration_counter += 1;
                            if (((FSI_Control)this.Control).splitting_fully_coupled == false) {
                                break;
                            }
                            if (iteration_counter > ((FSI_Control)this.Control).max_iterations_fully_coupled) {
                                break;// throw new ApplicationException("no convergence in coupled iterative solver, number of iterations: " + iteration_counter);
                            }
                        }
                        if (phystime == 0) {
                            if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty) && iteration_counter == 0) {
                                Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                                string firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#Timestep", "#Time", "P0_PosX", "P0_PosY", "P0_angle", "P0_VelX", "P0_VelY", "xPosition", "TotalKineticEnergy", "TotalMomentum");
                                Log_DragAndLift.WriteLine(firstline);

                                if (m_Particles.Count > 1) {
                                    Log_DragAndLift_P1 = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData_P1", CurrentSessionInfo.ID);
                                    string firstline_P1 = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#Timestep", "#Time", "P1_PosX", "P1_PosY", "P1_angle", "P1_VelX", "P1_VelY", "xPosition", "TotalKineticEnergy", "TotalMomentum");
                                    Log_DragAndLift_P1.WriteLine(firstline_P1);
                                }
                            }
                        }
                    } else {
                        foreach (Particle p in m_Particles) {
                            p.iteration_counter_P = -1;
                            p.forceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                        }
                        m_BDF_Timestepper.Solve(phystime, dt, false);
                    }
                }


                this.ResLogger.NextTimestep(false);


                return dt;
            }
        }
        #endregion

        #region restart
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
        #endregion

        #region Initialize particles
        protected override void SetInitial() {
            // Setup particles
            m_Particles = ((FSI_Control)this.Control).Particles;
            hack_phystime = 0.0;
            UpdateLevelSetParticles(0.0);

            // call base implementation
            base.SetInitial();

            // Setup Collision Model
            m_collisionModel = ((FSI_Control)this.Control).collisionModel;


            foreach (Particle p in m_Particles) {
                p.m_collidedWithParticle = new bool[m_Particles.Count];
                p.m_collidedWithWall = new bool[4];
                p.m_closeInterfacePointTo = new double[m_Particles.Count][];

            }
        }

        List<Particle> m_Particles;
        bool collision = false;

        bool triggerOnlyCollisionProcedure = false;
        #endregion

        #region Collision models
        /// <summary>
        /// Update collisionforces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="particle0"></param>
        /// <param name="particle1"></param>
        public void UpdateCollisionForces(List<Particle> particles, double hmin) {


            if (m_collisionModel == FSI_Control.CollisionModel.NoCollisionModel) {
                return;
            }


            if (particles.Count < 2)
                return;

            // Most of the code resulted from old one, should be simplified soon
            for (int i = 0; i < particles.Count; i++) {
                for (int j = i + 1; j < particles.Count; j++) {
                    var particle0 = particles[i];
                    var particle1 = particles[j];

                    var particle0CutCells = particle0.cutCells_P(LsTrk);
                    var particle1CutCells = particle1.cutCells_P(LsTrk);

                    var particleCutCellArray_P0 = particle0CutCells.ItemEnum.ToArray();
                    var neighborCellsArray_P0 = particle0CutCells.AllNeighbourCells().ItemEnum.ToArray();
                    var allCellsArray_P0 = particleCutCellArray_P0.Concat(neighborCellsArray_P0).ToArray();
                    var allCells_P0 = new CellMask(GridData, neighborCellsArray_P0);

                    var particleCutCellArray_P1 = particle1CutCells.ItemEnum.ToArray();
                    var neighborCellsArray_P1 = particle1CutCells.AllNeighbourCells().ItemEnum.ToArray();
                    var allCellsArray_P1 = particleCutCellArray_P1.Concat(neighborCellsArray_P1).ToArray();
                    var allCells_P1 = new CellMask(GridData, neighborCellsArray_P1);

                    double distance = 1E20;
                    double[] distanceVec = new double[Grid.SpatialDimension];

                    var interSecMask = allCells_P0.Intersect(allCells_P1);

                    var p0intersect = interSecMask.AllNeighbourCells().Intersect(particle0CutCells);
                    var p1intersect = interSecMask.AllNeighbourCells().Intersect(particle1CutCells);

                    // If there is no element neighbour of both particle cut cells return
                    if (!interSecMask.IsEmpty) {

                        // All interface points at a specific subgrid containing all cut cells of one particle
                        var interfacePoints_P0 = BoSSS.Solution.XNSECommon.XNSEUtils.GetInterfacePoints(LsTrk, LevSet, new SubGrid(particle0CutCells));
                        var interfacePoints_P1 = BoSSS.Solution.XNSECommon.XNSEUtils.GetInterfacePoints(LsTrk, LevSet, new SubGrid(particle1CutCells));

                        var tempDistance = 0.0;
                        double[] tempPoint_P0 = new double[2] { 0.0, 0.0 };
                        double[] tempPoint_P1 = new double[2] { 0.0, 0.0 };

                        if (interfacePoints_P0 != null && interfacePoints_P1 != null) {

                            for (int f = 0; f < interfacePoints_P0.NoOfRows; f++) {
                                for (int g = 0; g < interfacePoints_P1.NoOfRows; g++) {
                                    tempDistance = Math.Sqrt((interfacePoints_P0.GetRow(f)[0] - interfacePoints_P1.GetRow(g)[0]).Pow2() + (interfacePoints_P0.GetRow(f)[1] - interfacePoints_P1.GetRow(g)[1]).Pow2());
                                    if (tempDistance < distance) {

                                        distanceVec = interfacePoints_P0.GetRow(f).CloneAs();
                                        distanceVec.AccV(-1, interfacePoints_P1.GetRow(g));
                                        tempPoint_P0 = interfacePoints_P0.GetRow(f);
                                        tempPoint_P1 = interfacePoints_P1.GetRow(g);
                                        distance = tempDistance;
                                    }
                                }
                            }
                        }

                        double realDistance = distance;
                        bool ForceCollision = false;

                        // Important to get normal vector if distance is overlapping in the next timestep
                        if (realDistance <= 0.0 && particle0.m_closeInterfacePointTo[m_Particles.IndexOf(particle1)] != null && particle1.m_closeInterfacePointTo[m_Particles.IndexOf(particle0)] != null) {
                            tempPoint_P0 = particle0.m_closeInterfacePointTo[m_Particles.IndexOf(particle1)];
                            tempPoint_P1 = particle1.m_closeInterfacePointTo[m_Particles.IndexOf(particle0)];
                            if (tempPoint_P0 == null || tempPoint_P1 == null)
                                Console.WriteLine("Overlap of particles which were not close in the previous timestep due to one timestep - this is just an output without assertion");
                            //    throw new ApplicationException("Overlap of particles which were not close in the previous timestep due to in one timestep");
                            distanceVec = tempPoint_P0.CloneAs();
                            distanceVec.AccV(-1, tempPoint_P1);
                            //ForceCollision = true;
                        }


                        particle0.m_closeInterfacePointTo[m_Particles.IndexOf(particle1)] = tempPoint_P0;
                        particle1.m_closeInterfacePointTo[m_Particles.IndexOf(particle0)] = tempPoint_P1;


                        double threshold = 2.5 * hmin;

                        double eps = threshold.Pow2() / 2; // Turek paper
                        double epsPrime = threshold / 2; // Turek paper

                        double[] collisionForce;

                        var massDifference = Math.Abs(this.Control.PhysicalParameters.rho_A - particle0.particleDensity);

                        Console.WriteLine("realDistance: " + realDistance);
                        Console.WriteLine("Threshold: " + threshold);
                        Console.WriteLine("hmin: " + hmin);


                        // test of Modell 2

                        switch (m_collisionModel) {

                            case (FSI_Solver.FSI_Control.CollisionModel.RepulsiveForce):
                                if ((realDistance <= threshold)) {
                                    distanceVec.ScaleV((threshold - realDistance).Pow2());
                                    distanceVec.ScaleV(1 / eps);

                                    Console.WriteLine("Strongly recommended to use conservation of momentum collision model. This one is highly experimental!!!!");

                                    collisionForce = distanceVec;
                                    var collisionForceP1 = collisionForce.CloneAs();
                                    collisionForce.ScaleV(-100.0);
                                    collisionForceP1.ScaleV(-100.0);
                                    particle0.hydrodynForcesAtIteration[0].AccV(-1, collisionForce);
                                    //particle0.hydrodynTorqueAtIteration[0] += 100 * (collisionForce[0] * (tempPoint_P0[0] - particle0.positionAtIteration[0][0]) + collisionForce[1] * (tempPoint_P0[1] - particle0.positionAtIteration[0][1]));
                                    particle1.hydrodynForcesAtIteration[0].AccV(1, collisionForceP1);
                                    //particle1.hydrodynTorqueAtIteration[0] += -100 * (collisionForceP1[0] * (tempPoint_P1[0] - particle1.positionAtIteration[0][0]) + collisionForceP1[1] * (tempPoint_P1[1] - particle1.positionAtIteration[0][1]));
                                    Console.WriteLine("Collision information: Particles coming close, force " + collisionForce.L2Norm());
                                    Console.WriteLine("Collision information: Particles coming close, torque " + particle1.hydrodynTorqueAtIteration[0]);

                                    if (realDistance <= 1.5 * hmin) {
                                        Console.WriteLine("Entering overlapping loop....");
                                        triggerOnlyCollisionProcedure = true;
                                    }

                                }
                                break;


                            case (FSI_Solver.FSI_Control.CollisionModel.MomentumConservation):

                                if (((realDistance <= threshold) || ForceCollision) && !particle0.m_collidedWithParticle[m_Particles.IndexOf(particle1)] && !particle1.m_collidedWithParticle[m_Particles.IndexOf(particle0)]) {



                                    // Bool if collided
                                    particle0.m_collidedWithParticle[m_Particles.IndexOf(particle1)] = true;
                                    particle1.m_collidedWithParticle[m_Particles.IndexOf(particle0)] = true;

                                    // Bool if force integration should be skipped
                                    particle0.skipForceIntegration = true;
                                    particle1.skipForceIntegration = true;

                                    //coefficient of restitution (e=0 pastic; e=1 elastic)
                                    double e = 1.0;

                                    //collision Nomal
                                    var normal = distanceVec.CloneAs();
                                    normal.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));

                                    double[] tangential = new double[] { -normal[1], normal[0] };


                                    //general definitions of normal and tangential components
                                    double collisionVn_P0 = particle0.transVelocityAtIteration[0][0] * normal[0] + particle0.transVelocityAtIteration[0][1] * normal[1];
                                    double collisionVt_P0 = particle0.transVelocityAtIteration[0][0] * tangential[0] + particle0.transVelocityAtIteration[0][1] * tangential[1];
                                    double collisionVn_P1 = particle1.transVelocityAtIteration[0][0] * normal[0] + particle1.transVelocityAtIteration[0][1] * normal[1];
                                    double collisionVt_P1 = particle1.transVelocityAtIteration[0][0] * tangential[0] + particle1.transVelocityAtIteration[0][1] * tangential[1];

                                    // exzentric collision
                                    // ----------------------------------------                                                                  
                                    tempPoint_P0.AccV(-1, particle0.positionAtIteration[0]);
                                    double a0 = (tempPoint_P0[0] * tangential[0] + tempPoint_P0[1] * tangential[1]);
                                    tempPoint_P1.AccV(-1, particle1.positionAtIteration[0]);
                                    double a1 = (tempPoint_P1[0] * tangential[0] + tempPoint_P1[1] * tangential[1]);

                                    // Fix for Sphere
                                    // ----------------------------------------  
                                    if (particle0 is Particle_Sphere)
                                        a0 = 0.0;
                                    if (particle1 is Particle_Sphere)
                                        a1 = 0.0;


                                    double Fx = (1 + e) * ((collisionVn_P0 - collisionVn_P1) / (1 / particle0.Mass_P + 1 / particle1.Mass_P + a0.Pow2() / particle0.MomentOfInertia_P + a1.Pow2() / particle1.MomentOfInertia_P));
                                    double Fxrot = (1 + e) * ((-a0 * particle0.rotationalVelocityAtIteration[0] + a1 * particle1.rotationalVelocityAtIteration[0]) / (1 / particle0.Mass_P + 1 / particle1.Mass_P + a0.Pow2() / particle0.MomentOfInertia_P + a1.Pow2() / particle1.MomentOfInertia_P));

                                    double tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / particle0.Mass_P;
                                    double tempCollisionVn_P1 = collisionVn_P1 + (Fx + Fxrot) / particle1.Mass_P;
                                    double tempCollisionVt_P0 = collisionVt_P0;
                                    double tempCollisionVt_P1 = collisionVt_P1;
                                    Console.WriteLine("a0:    " + a0 + "   Fx:    " + (-Fx) + "      Fxrot:    " + (-Fxrot));
                                    Console.WriteLine("a1:    " + a1 + "   Fx:    " + Fx + "      Fxrot:    " + Fxrot);
                                    particle0.rotationalVelocityAtIteration[0] = particle0.rotationalVelocityAtIteration[0] + a0 * (Fx + Fxrot) / particle0.MomentOfInertia_P;
                                    particle1.rotationalVelocityAtIteration[0] = particle1.rotationalVelocityAtIteration[0] - a1 * (Fx + Fxrot) / particle1.MomentOfInertia_P;


                                    // zentric collision
                                    // ----------------------------------------
                                    //double tempCollisionVn_P0 = (particle0.mass_P * collisionVn_P0 + particle1.mass_P * collisionVn_P1 + e * particle1.mass_P * (collisionVn_P1 - collisionVn_P0)) / (particle0.mass_P + particle1.mass_P);
                                    //double tempCollisionVt_P0 = collisionVt_P0;
                                    //double tempCollisionVn_P1 = (particle0.mass_P * collisionVn_P0 + particle1.mass_P * collisionVn_P1 + e * particle0.mass_P * (collisionVn_P0 - collisionVn_P1)) / (particle0.mass_P + particle1.mass_P);
                                    //double tempCollisionVt_P1 = collisionVt_P1;
                                    // ----------------------------------------


                                    particle0.transVelocityAtIteration[0] = new double[] { normal[0] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[0], normal[1] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[1] };
                                    particle1.transVelocityAtIteration[0] = new double[] { normal[0] * tempCollisionVn_P1 + tempCollisionVt_P1 * tangential[0], normal[1] * tempCollisionVn_P1 + tempCollisionVt_P1 * tangential[1] };
                                    //collided = true;

                                    //double contactForce = (1 + e)*(particle0.transVelocityAtIteration[0][0] - particle0.radius_P * particle0.rotationalVelocityAtIteration[0] - (particle1.transVelocityAtIteration[0][0] - particle1.radius_P * particle1.rotationalVelocityAtIteration[0])) / (1/particle0.mass_P+1/particle1.mass_P+particle0.radius_P.Pow2()/particle0.MomentOfInertia_P+particle1.radius_P.Pow2()/particle1.MomentOfInertia_P);
                                    //particle0.transVelocityAtIteration[0][0] -= contactForce / particle0.mass_P;
                                    //particle0.rotationalVelocityAtIteration[0] = particle0.rotationalVelocityAtIteration[0];
                                    //particle0.rotationalVelocityAtIteration[0] += particle0.radius_P * contactForce / particle0.MomentOfInertia_P;
                                    //particle1.transVelocityAtIteration[0][0] -= contactForce / particle1.mass_P;
                                    //particle1.rotationalVelocityAtIteration[0] = particle1.rotationalVelocityAtIteration[0];
                                    //particle1.rotationalVelocityAtIteration[0] += particle1.radius_P * contactForce / particle1.MomentOfInertia_P;

                                    if (realDistance <= 1.5 * hmin) {
                                        Console.WriteLine("Entering overlapping loop....");
                                        triggerOnlyCollisionProcedure = true;
                                    }
                                }

                                if (realDistance > 1.5 * hmin && particle0.m_collidedWithParticle[m_Particles.IndexOf(particle1)] && particle1.m_collidedWithParticle[m_Particles.IndexOf(particle0)]) {
                                    particle0.m_collidedWithParticle[m_Particles.IndexOf(particle1)] = false;
                                    particle1.m_collidedWithParticle[m_Particles.IndexOf(particle0)] = false;
                                    particle0.m_closeInterfacePointTo[m_Particles.IndexOf(particle1)] = null;
                                    particle1.m_closeInterfacePointTo[m_Particles.IndexOf(particle0)] = null;
                                }

                                ForceCollision = false;
                                break;


                            default:
                                throw new NotImplementedException("Collision model not available");
                        }
                    }
                }
            }
        }


        private FSI_Solver.FSI_Control.CollisionModel m_collisionModel;

        /// <summary>
        /// Calculation of collision forces between particle and wall
        /// </summary>
        /// <param name="particles"></param>
        /// <param name="hmin"></param>
        public void WallCollisionForces(Particle particle, double hmin) {

            if (m_collisionModel == FSI_Control.CollisionModel.NoCollisionModel) {
                return;
            }

            var particleCutCells = particle.cutCells_P(LsTrk);

            var particleCutCellArray = particleCutCells.ItemEnum.ToArray();
            var neighborCellsArray = particleCutCells.AllNeighbourCells().ItemEnum.ToArray();
            var allCellsArray = particleCutCellArray.Concat(neighborCellsArray).ToArray();
            var allCells = new CellMask(GridData, neighborCellsArray);

            collision = false;

            double distance = double.MaxValue;
            double[] distanceVec = new double[Grid.SpatialDimension];

            // All interface points at a specific subgrid containing all cut cells of one particle
            MultidimensionalArray interfacePoints = null;

            Console.WriteLine("ParticleCutCellCount:   " + particleCutCells.Count());

            var trafo = GridData.iGeomEdges.Edge2CellTrafos;

            SubGrid allCellsGrid = new SubGrid(allCells);

            double[] tempPoint = new double[2] { 0.0, 0.0 };

            foreach (int iEdge in allCellsGrid.BoundaryEdgesMask.ItemEnum) {

                // Collision forces have to act
                if (GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge)) {

                    if (interfacePoints == null)
                        interfacePoints = BoSSS.Solution.XNSECommon.XNSEUtils.GetInterfacePoints(LsTrk, LevSet, new SubGrid(particleCutCells));

                    collision = true;
                    var jCell = GridData.iGeomEdges.CellIndices[iEdge, 0];
                    int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

                    NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
                    NodeSet Nodes = refNodes.ElementAt(iKref);

                    var trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[iEdge, 0];
                    var transFormed = trafo[trafoIdx].Transform(Nodes);
                    var newVertices = transFormed.CloneAs();
                    GridData.TransformLocal2Global(transFormed, newVertices, jCell);
                    var tempDistance = 0.0;

                    for (int i = 0; i < interfacePoints.NoOfRows; i++) {
                        for (int j = 0; j < newVertices.NoOfRows; j++) {
                            tempDistance = Math.Sqrt((interfacePoints.GetRow(i)[0] - newVertices.GetRow(j)[0]).Pow2() + (interfacePoints.GetRow(i)[1] - newVertices.GetRow(j)[1]).Pow2());
                            if (tempDistance < distance) {
                                tempPoint = interfacePoints.GetRow(i);
                                distanceVec = interfacePoints.GetRow(i).CloneAs();
                                distanceVec.AccV(-1, newVertices.GetRow(j));
                                distance = tempDistance;
                            }

                        }


                    }
                }
            }

            double realDistance = distance;

            if (collision == false)
                return;


            Console.WriteLine("Closes Distance to wall is: " + distance);

            double threshold = 1.5 * hmin;

            double eps = threshold.Pow2() / 2; // Turek paper
            double epsPrime = threshold / 2; // Turek paper

            double[] collisionForce;


            switch (m_collisionModel) {

                case (FSI_Solver.FSI_Control.CollisionModel.RepulsiveForce):
                    if ((realDistance <= threshold)) {
                        Console.WriteLine("Strongly recommended to use conservation of momentum collision model. This one is highly experimental!!!!");
                        // Modell 1
                        distanceVec.ScaleV(1 / eps);
                        distanceVec.ScaleV(((threshold - realDistance).Abs()));

                        collisionForce = distanceVec;
                        collisionForce.ScaleV(100.0);

                        particle.hydrodynForcesAtIteration[0] = collisionForce;


                        return;
                    }


                    if (realDistance <= (1.5 * hmin)) {

                        distanceVec.ScaleV((threshold - realDistance).Abs());
                        distanceVec.ScaleV(1 / epsPrime);
                        collisionForce = distanceVec;

                        collisionForce.ScaleV(100.0);
                        particle.hydrodynForcesAtIteration[0].AccV(1, collisionForce);
                        particle.hydrodynTorqueAtIteration[0] -= (collisionForce[0] * (tempPoint[0] - particle.positionAtIteration[0][0]) + collisionForce[1] * (tempPoint[1] - particle.positionAtIteration[0][1]));
                        Console.WriteLine("Collision information: Wall overlapping, force X " + collisionForce[0]);
                        Console.WriteLine("Collision information: Wall overlapping, force Y " + collisionForce[1]);

                        if (realDistance <= 1.5 * hmin) {
                            Console.WriteLine("Entering wall overlapping loop....");
                            triggerOnlyCollisionProcedure = true;

                        }
                        return;
                    }
                    break;

                case (FSI_Solver.FSI_Control.CollisionModel.MomentumConservation):


                    if (realDistance <= (threshold) && !particle.m_collidedWithWall[0]) {

                        //coefficient of restitution (e=0 pastic; e=1 elastic)


                        double e = 1.0;

                        // Fully plastic for bottom wall
                        if (particle.positionAtIteration[0][1] < 0.5)
                            e = 0.0;


                        // if particle already collided with wall
                        particle.m_collidedWithWall[0] = true;

                        // Skip force integration for next timestep
                        particle.skipForceIntegration = false;

                        //collision Nomal
                        var normal = distanceVec.CloneAs();
                        normal.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));
                        double[] tangential = new double[] { -normal[1], normal[0] };


                        double collisionVn_P0 = particle.transVelocityAtIteration[0][0] * normal[0] + particle.transVelocityAtIteration[0][1] * normal[1];
                        double collisionVt_P0 = particle.transVelocityAtIteration[0][0] * tangential[0] + particle.transVelocityAtIteration[0][1] * tangential[1];


                        // exzentric collision
                        // ----------------------------------------
                        tempPoint.AccV(-1, particle.positionAtIteration[0]);
                        double a0 = (tempPoint[0] * tangential[0] + tempPoint[1] * tangential[1]);

                        if (particle is Particle_Sphere)
                            a0 = 0.0;


                        double Fx = (1 + e) * (collisionVn_P0) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);
                        double Fxrot = (1 + e) * (-a0 * particle.rotationalVelocityAtIteration[0]) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);

                        double tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / particle.Mass_P;
                        double tempCollisionVt_P0 = collisionVt_P0;

                        particle.rotationalVelocityAtIteration[0] = particle.rotationalVelocityAtIteration[0] + a0 * (Fx + Fxrot) / particle.MomentOfInertia_P;


                        particle.transVelocityAtIteration[0] = new double[] { normal[0] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[0], normal[1] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[1] };

                    }
                    if (realDistance > threshold && particle.m_collidedWithWall[0]) {
                        Console.WriteLine("Reset Wall");
                        particle.m_collidedWithWall[0] = false;
                    }
                    break;

                default:
                    throw new NotImplementedException("Collision model not available");
            }
        }
        #endregion

        #region Mesh refinement
        /// <summary>
        /// Very primitive refinement indicator, works on a LevelSet criterion.
        /// </summary>
        /// 
        int LevelIndicator(int j, int CurrentLevel) {
            var LevSetCells = LsTrk.Regions.GetCutCellMask();
            var LevSetNeighbours = LsTrk.Regions.GetNearFieldMask(1);
            //    var LevSetNeighboursNeighbours = LevSetNeighbours.AllNeighbourCells();
            ////    var LevSetNeighboursNeighboursNeighbours = LevSetNeighbours.AllNeighbourCells();

            int DesiredLevel_j = 0;
            //    if (LevSetCells.Contains(j))
            //    {
            //        if (CurrentLevel < ((FSI_Control)Control).RefinementLevel)
            //        {
            //            DesiredLevel_j = ((FSI_Control)Control).RefinementLevel;
            //        }

            //        //else if (((FSI_Control)this.Control).Timestepper_Mode != FSI_Control.TimesteppingMode.MovingMesh)
            //        //{
            //        //    double curv_max = 1.0 / (this.Control.maxCurvature * ((GridData)this.GridData).Cells.h_min[j]);
            //        //    double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
            //        //    double curv_thrshld = mean_curv;

            //        //    if (mean_curv > curv_max)
            //        //    {
            //        //        DesiredLevel_j = CurrentLevel + 1;
            //        //    }
            //        //    else if (mean_curv < (curv_max / 5))
            //        //    {

            //        //        DesiredLevel_j = CurrentLevel - 1;
            //        //    }
            //        //}
            //    }
            //else if (LevSetNeighbours.Contains(j) && ((FSI_Control)this.Control).Timestepper_Mode != FSI_Control.TimesteppingMode.MovingMesh)
            //{
            //    if (CurrentLevel < ((FSI_Control)Control).RefinementLevel)
            //    {
            //        DesiredLevel_j = CurrentLevel + 1;
            //    }
            //    else if (CurrentLevel > ((FSI_Control)Control).RefinementLevel)
            //    {
            //        DesiredLevel_j = CurrentLevel - 1;
            //    }
            //    else
            //    {
            //        DesiredLevel_j = CurrentLevel;
            //    }
            //}
            //else if (LevSetNeighboursNeighbours.Contains(j) && ((FSI_Control)this.Control).Timestepper_Mode != FSI_Control.TimesteppingMode.MovingMesh)
            //{
            //    if (CurrentLevel < ((FSI_Control)Control).RefinementLevel)
            //    {
            //        DesiredLevel_j = CurrentLevel + 1;
            //    }
            //    else if (CurrentLevel > ((FSI_Control)Control).RefinementLevel)
            //    {
            //        DesiredLevel_j = CurrentLevel - 1;
            //    }
            //    else
            //    {
            //        DesiredLevel_j = CurrentLevel;
            //    }
            //}
            //else if (DesiredLevel_j > 0)
            //{
            //    DesiredLevel_j = CurrentLevel - 1;
            //}

            //return DesiredLevel_j
            if (LevSetCells.Contains(j)) {
                DesiredLevel_j = ((FSI_Control)Control).RefinementLevel;
            } else if (LevSetNeighbours.Contains(j)) {
                DesiredLevel_j = ((FSI_Control)Control).RefinementLevel;
            }

            return DesiredLevel_j;
        }


        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {

            if (((FSI_Control)Control).AdaptiveMeshRefinement) {

                if (TimestepNo > 3 && TimestepNo % 3 != 0) {
                    newGrid = null;
                    old2NewGrid = null;
                    return;
                }

                // Check grid changes
                // ==================

                // compute curvature for levelindicator 
                CurvatureAlgorithms.CurvatureDriver(
                SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                CurvatureAlgorithms.FilterConfiguration.Default,
                this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                this.HMForder, this.DGLevSet.Current);

                CellMask CutCells = LsTrk.Regions.GetCutCellMask();
                CellMask CutCellNeighbors = LsTrk.Regions.GetNearFieldMask(1);
                var CutCellArray = CutCells.ItemEnum.ToArray();
                var CutCellNeighborsArray = CutCellNeighbors.ItemEnum.ToArray();
                var AllCells = CutCellArray.Concat(CutCellNeighborsArray).ToArray();
                var NoCoarseningcells = new CellMask(this.GridData, AllCells);

                // Only CutCells are NoCoarseningCells 
                bool AnyChange = GridRefinementController.ComputeGridChange((GridData)(this.GridData), CutCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if (AnyChange) {
                    int[] glb = (new int[] {
                    CellsToRefineList.Count,
                    Coarsening.Sum(L => L.Length),
                }).MPISum();
                    NoOfCellsToRefine = glb[0];
                    NoOfCellsToCoarsen = glb[1];
                }
                int oldJ = this.GridData.CellPartitioning.TotalLength;

                // Update Grid
                // ===========

                if (AnyChange) {

                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                    newGrid = ((GridData)(this.GridData)).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                    if (this.Control.savetodb == true) {
                        //Console.WriteLine("Save adaptive Mesh...");
                        //Console.WriteLine("GridGUID:   " + newGrid.GridGuid);
                        //DatabaseDriver.SaveGrid(newGrid, base.GetDatabase());
                        //Console.WriteLine("...done");
                    }
                } else {

                    Console.WriteLine("No changes in Grid");
                    newGrid = null;
                    old2NewGrid = null;
                }

                //debug = false;

            } else {

                newGrid = null;
                old2NewGrid = null;
            }
        }
        #endregion
    }
}




#endregion