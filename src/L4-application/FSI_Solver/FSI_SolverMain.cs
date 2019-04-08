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
using BoSSS.Solution.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.IO;
using System.Diagnostics;
using System.IO;
using ilPSP;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Grid.Classic;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using BoSSS.Foundation.Grid.RefElements;
using FSI_Solver;
using System.Collections;

namespace BoSSS.Application.FSI_Solver {
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain {

        static int counter = 0;

        public static void MegaArschKakke2(DGField[] f) {
            int rank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
            

            Tecplot.PlotFields(f, "MegaArschKakke-" + counter, 0.0, 2);
            counter++;
        }

        bool InitializedColor = false;

        double calculatedDampingTensors;
        // =============================
        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {
            MultiphaseCellAgglomerator.Katastrophenplot = MegaArschKakke2;

            _Main(args, false, delegate () {
                var p = new FSI_SolverMain();
                return p;
            });
        }

        /// <summary>
        /// Curvature; DG-polynomial degree should be 2 times the polynomial degree of <see cref="LevSet"/>.
        /// </summary>
        [InstantiateFromControlFile("Curvature", "Curvature", IOListOption.ControlFileDetermined)]
        public SinglePhaseField Curvature;

        SinglePhaseField ParticleColor;

        SinglePhaseField LevelSetDistance;

        protected override void CreateFields() {
            base.CreateFields();

            ParticleColor = new SinglePhaseField(new Basis(this.GridData, 0), "ParticleColor");
            m_RegisteredFields.Add(ParticleColor);
            m_IOFields.Add(ParticleColor);

            LevelSetDistance = new SinglePhaseField(new Basis(this.GridData, 0), "LevelSetDistance");
            m_RegisteredFields.Add(LevelSetDistance);
            m_IOFields.Add(LevelSetDistance);
        }

        // Create equations and solvers
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

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            if (IBM_Op != null)
                return;
            string[] CodNameSelected = new string[0];
            string[] DomNameSelected = new string[0];

            int D = this.GridData.SpatialDimension;

            BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

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

            // Momentum equation
            // =============================
            // Convective part
            // =============================
            {
                if (IBM_Op_config.convection)
                {
                    for (int d = 0; d < D; d++)
                    {

                        var comps = IBM_Op.EquationComponents[CodName[d]];

                        //var ConvBulk = new Solution.XNSECommon.Operator.Convection.ConvectionInBulk_LLF(D, BcMap, d, this.Control.PhysicalParameters.rho_A, 0, this.Control.AdvancedDiscretizationOptions.LFFA, this.Control.AdvancedDiscretizationOptions.LFFB, LsTrk);
                        var ConvBulk = new Solution.NSECommon.LinearizedConvection(D, BcMap, d);
                        //IBM_Op.OnIntegratingBulk += ConvBulk.SetParameter;
                        comps.Add(ConvBulk); // bulk component

                        //var ConvIB = new BoSSS.Solution.XNSECommon.Operator.Convection.ConvectionAtIB(d, D, LsTrk, IBM_Op_config.dntParams.LFFA, BcMap, uIBM, wIBM);

                        if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None)
                        {

                            var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, D, LsTrk,
                                this.Control.AdvancedDiscretizationOptions.LFFA, BcMap,
                                delegate (double[] X, double time)
                                {
                                    throw new NotImplementedException("Currently not implemented for fixed motion");
                                //return new double[] { 0.0, 0.0 };
                            },
                                this.Control.PhysicalParameters.rho_A,
                                UseMovingMesh);
                            comps.Add(ConvIB); // immersed boundary component
                        }
                        else
                        {
                            var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, D, LsTrk,
                                this.Control.AdvancedDiscretizationOptions.LFFA, BcMap,
                                    delegate (double[] X, double time)
                                    {

                                        double[] result = new double[X.Length + 3];

                                        foreach (Particle p in m_Particles)
                                        {
                                        // Separating different boundary regions (for active particles)
                                        double cos_theta;
                                        // The posterior side of the particle (Neumann boundary)
                                        if (Math.Cos(p.Angle[0]) * (X[0] - p.Position[0][0]) + Math.Sin(p.Angle[0]) * (X[1] - p.Position[0][1]) < 1e-8)// && Math.Cos(p.particleAnglePerIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.particleAnglePerIteration[0]) * (X[1] - p.positionAtIteration[0][1]) > -0.25)
                                        {
                                                cos_theta = (Math.Cos(p.Angle[0]) * (X[0] - p.Position[0][0]) + Math.Sin(p.Angle[0]) * (X[1] - p.Position[0][1])) / (Math.Sqrt((X[0] - p.Position[0][0]).Pow2() + (X[1] - p.Position[0][1]).Pow2()));
                                            }
                                        // The anterior side of the particle (Dirichlet boundary)
                                        else
                                            {
                                                cos_theta = 0;
                                            }

                                        // which particle?
                                        bool containsParticle;
                                            if (m_Particles.Count == 1)
                                            {
                                                containsParticle = true;
                                            }
                                            else { containsParticle = p.Contains(X, LsTrk); }

                                        // active particles
                                        if (containsParticle && p.ActiveParticle == true)
                                            {
                                                result[0] = p.TranslationalVelocity[0][0];
                                                result[1] = p.TranslationalVelocity[0][1];
                                                result[2] = p.RotationalVelocity[0];
                                                result[3] = p.Position[0].L2Distance(X);
                                                result[4] = -cos_theta;
                                                return result;
                                            }

                                        // passive particles
                                        else if (containsParticle && p.ActiveParticle == false)
                                            {
                                                result[0] = p.TranslationalVelocity[0][0];
                                                result[1] = p.TranslationalVelocity[0][1];
                                                result[2] = p.RotationalVelocity[0];
                                                result[3] = p.Position[0].L2Distance(X);
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
            }

            // Pressure part
            // =============================
            {
                if (IBM_Op_config.PressureGradient)
                {
                    for (int d = 0; d < D; d++)
                    {
                        var comps = IBM_Op.EquationComponents[CodName[d]];
                        //var pres = new Solution.XNSECommon.Operator.Pressure.PressureInBulk(d, BcMap, 1, 0);
                        var pres = new PressureGradientLin_d(d, BcMap);
                        //IBM_Op.OnIntegratingBulk += pres.SetParameter;
                        comps.Add(pres); // bulk component

                        var presLs = new BoSSS.Solution.NSECommon.Operator.Pressure.ActivePressureAtIB(d, D, LsTrk);//no changes necessary for impl. of active particles
                        comps.Add(presLs); // immersed boundary component

                        // if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                        if (this.Control.FixedStreamwisePeriodicBC)
                        {
                            var presSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                            comps.Add(presSource);
                        }
                    }
                }
            }

            // Viscous part
            // =============================
            {
                if (IBM_Op_config.Viscous)
                {
                    for (int d = 0; d < D; d++)
                    {
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

                        if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None)
                        {

                            var ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, D, LsTrk,
                                penalty, this.ComputePenaltyIB,
                                this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                                delegate (double[] X, double time)
                                {
                                    throw new NotImplementedException("Currently not implemented for fixed motion");
                                //return new double[] { 0.0, 0.0 };
                            });
                            comps.Add(ViscLs); // immersed boundary component

                        }
                        else
                        {
                            var ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, D, LsTrk,
                                penalty, this.ComputePenaltyIB,
                                this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                                delegate (double[] X, double time)
                                {

                                    double[] result = new double[X.Length + 5];

                                    foreach (Particle p in m_Particles)
                                    {
                                    // which particle?
                                    bool containsParticle;
                                        if (m_Particles.Count == 1)
                                        {
                                            containsParticle = true;
                                        }
                                        else { containsParticle = p.Contains(X, LsTrk); }

                                    // active particles
                                    if (containsParticle && p.ActiveParticle == true)
                                        {
                                        // Separating different boundary regions (for active particles)
                                        double cos_theta;
                                        // The posterior side of the particle (Neumann boundary)
                                        if (Math.Cos(p.Angle[0]) * (X[0] - p.Position[0][0]) + Math.Sin(p.Angle[0]) * (X[1] - p.Position[0][1]) < 1e-8)// && Math.Cos(p.particleAnglePerIteration[0]) * (X[0] - p.positionAtIteration[0][0]) + Math.Sin(p.particleAnglePerIteration[0]) * (X[1] - p.positionAtIteration[0][1]) > -0.25)
                                        {
                                                cos_theta = (Math.Cos(p.Angle[0]) * (X[0] - p.Position[0][0]) + Math.Sin(p.Angle[0]) * (X[1] - p.Position[0][1])) / (Math.Sqrt((X[0] - p.Position[0][0]).Pow2() + (X[1] - p.Position[0][1]).Pow2()));
                                            }
                                        // The anterior side of the particle (Dirichlet boundary)
                                        else
                                            {
                                                cos_theta = 0;
                                            }
                                            result[0] = p.TranslationalVelocity[0][0];
                                            result[1] = p.TranslationalVelocity[0][1];
                                            result[2] = p.RotationalVelocity[0];
                                            if (p is Particle_Sphere)
                                            {
                                                result[3] = ((Particle_Sphere)p).radius_P;
                                            }
                                            else
                                            {
                                                result[3] = p.Position[0].L2Distance(X);
                                            }
                                            result[4] = p.ActiveStress;
                                            result[5] = -cos_theta;
                                            result[6] = p.Angle[0];
                                        }

                                    // passive particles
                                    else if (containsParticle && p.ActiveParticle == false)
                                        {
                                            result[0] = p.TranslationalVelocity[0][0];
                                            result[1] = p.TranslationalVelocity[0][1];
                                            result[2] = p.RotationalVelocity[0];
                                            if (p is Particle_Sphere)
                                            {
                                                result[3] = ((Particle_Sphere)p).radius_P;
                                            }
                                            else
                                            {
                                                result[3] = p.Position[0].L2Distance(X);
                                            }
                                            result[4] = 0;
                                            result[5] = 0;
                                            result[6] = p.Angle[0];
                                        }
                                    }
                                    return result;
                                }
                             );
                            comps.Add(ViscLs); // immersed boundary component
                        }
                    }
                }
            }

            // Continuum equation
            // ==================
            {
                if (IBM_Op_config.continuity)
                {
                    for (int d = 0; d < D; d++)
                    {
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

                    if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None)
                    {

                        var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, LsTrk, 1, delegate (double[] X, double time)
                        {
                            throw new NotImplementedException("Currently not implemented for fixed motion");
                            //return new double[] { 0.0, 0.0 };
                        });
                        IBM_Op.EquationComponents["div"].Add(divPen);  // immersed boundary component
                    }
                    else
                    {
                        var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, LsTrk, 1,
                           delegate (double[] X, double time)
                           {

                               double[] result = new double[X.Length + 2];

                               foreach (Particle p in m_Particles)
                               {
                                   bool containsParticle;
                                   if (m_Particles.Count == 1)
                                   {
                                       containsParticle = true;
                                   }
                                   else { containsParticle = p.Contains(X, LsTrk); }
                                   if (containsParticle)
                                   {
                                       result[0] = p.TranslationalVelocity[0][0];
                                       result[1] = p.TranslationalVelocity[0][1];
                                       result[2] = p.RotationalVelocity[0];
                                       if (p is Particle_Sphere)
                                       {
                                           result[3] = ((Particle_Sphere)p).radius_P;
                                       }
                                       else
                                       {
                                           result[3] = p.Position[0].L2Distance(X);
                                       }
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
            // ------------------
            SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;
            if (this.Control.PhysicalParameters.IncludeConvection) {
                SpatialOp = SpatialOperatorType.Nonlinear;
            }

            // create timestepper, update level-set
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
                    if (calculatedDampingTensors == 0) {
                        foreach (Particle p in m_Particles) {
                            if (p.neglectAddedDamping == false) {
                                p.CalculateDampingTensor(LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                                ExchangeDampingTensors();
                            }
                        }
                    }
                    calculatedDampingTensors = 1;
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    if (calculatedDampingTensors == 0)
                    {
                        foreach (Particle p in m_Particles)
                        {
                            if (p.neglectAddedDamping == false)
                            {
                                p.CalculateDampingTensor(LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                            }
                        }
                    }
                    calculatedDampingTensors = 1;
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
                )
            {
                m_ResLogger = base.ResLogger,
                m_ResidualNames = ArrayTools.Cat(this.ResidualMomentum.Select(f => f.Identification), this.ResidualContinuity.Identification),
                IterUnderrelax = ((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? ((FSI_Control)this.Control).LSunderrelax : 1.0,
                Config_LevelSetConvergenceCriterion = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion,
                SessionPath = SessionPath,
                Timestepper_Init = Solution.Timestepping.TimeStepperInit.SingleInit
            };
        }

        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {
            
            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.None:
                    ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                    //newTransVelocity[0] = (((FSI_Control)this.Control).transVelocityFunc[0])(phystime);
                    //newTransVelocity[1] = (((FSI_Control)this.Control).transVelocityFunc[1])(phystime);
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
                        p.ForceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                    }
                    break;

                case LevelSetHandling.LieSplitting:
                    UpdateLevelSetParticles(dt);
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    UpdateLevelSetParticles(dt);
                    break;

                case LevelSetHandling.StrangSplitting:
                    UpdateLevelSetParticles(dt);
                    break;

                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
            }
        
            // Forces and Torque residual
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
                    acc_force_P_x_old += p.HydrodynamicForces[1][0];
                    acc_force_P_y_old += p.HydrodynamicForces[1][1];
                    acc_torque_P_old += p.HydrodynamicTorque[1];

                    // forces and torque of the current iteration
                    acc_force_P_x += p.HydrodynamicForces[0][0];
                    acc_force_P_y += p.HydrodynamicForces[0][1];
                    acc_torque_P += p.HydrodynamicTorque[0];
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
        /// Particle to Level-Set-Field 
        /// </summary>
        private void UpdateLevelSetParticles(double dt)
        {
            // Step 1
            // Define an array with the respective cell colors
            // ===============================================
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] CellColor = InitializedColor ? LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")] : InitializeColoring(J);

            // Step 2
            // Delete the old level set
            // ========================
            LevSet.Clear();

            // Step 3
            // Define level set per color
            // ==========================
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate();
            CellMask AgglParticleMask = null;
            List<int[]> ColoredCellsSorted = levelSetUpdate.ColoredCellsFindAndSort(CellColor);
            int[] ParticleColor = levelSetUpdate.FindParticleColor(GridData, m_Particles, ColoredCellsSorted);
            for (int p = 0; p < m_Particles.Count(); p++)
            {
                if (ParticleColor[p] != 0)
                {
                    int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(ParticleColor, ParticleColor[p]);
                    CellMask ColoredCellMask = levelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, ParticleColor[p], J);
                    AgglParticleMask = AgglParticleMask == null ? ColoredCellMask : AgglParticleMask.Union(ColoredCellMask);

                    double phiComplete(double[] X, double t)
                    {
                        // Generating the correct sign
                        // ===========================
                        double phi = Math.Pow(-1, ParticlesOfCurrentColor.Length - 1);

                        // Multiplication over all particle-level-sets within the current color
                        // ====================================================================
                        for (int pc = 0; pc < ParticlesOfCurrentColor.Length; pc++)
                        {
                            phi *= m_Particles[ParticlesOfCurrentColor[pc]].Phi_P(X);

                            // Delete all particles within the current color from the particle color array
                            // ===========================================================================
                            ParticleColor[ParticlesOfCurrentColor[pc]] = 0;
                        }
                        return phi;
                    }
                    SetLevelSet(phiComplete, ColoredCellMask, hack_phystime);
                }
            }

            // =======================================================
            // Step 4
            // Define level set of the remaining cells ("Fluid-Cells")
            // =======================================================
            double phiFluid(double[] X, double t)
            {
                return -1;
            }
            CellMask FluidCells = AgglParticleMask != null ? AgglParticleMask.Complement() : CellMask.GetFullMask(GridData);
            SetLevelSet(phiFluid, FluidCells, hack_phystime);

            // =======================================================
            // Step 5
            // Update level set tracker and coloring
            // =======================================================
            LsTrk.UpdateTracker(__NearRegionWith: 2);
            UpdateColoring();
        }

        /// <summary>
        /// Set level set based on the function phi and the current cells
        /// </summary>
        private void SetLevelSet(Func<double[], double, double> phi, CellMask CurrentCells, double phystime)
        {
            ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(phi, phystime);
            DGLevSet.Current.ProjectField(Function);
            LevSet.AccLaidBack(1.0, DGLevSet.Current, CurrentCells);
        }

        /// <summary>
        /// Update of <see cref="ParticleColor"/> and <see cref="LevelSetDistance"/>
        /// </summary>
        private void UpdateColoring()
        {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] PartCol = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            var rCode = LsTrk.Regions.RegionsCode;
            for (int j = 0; j < J; j++)
            {
                ParticleColor.SetMeanValue(j, PartCol[j]);

                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(rCode[j], 0));
            }
        }

        /// <summary>
        /// Initialization of <see cref="ParticleColor"/> 
        /// </summary>
        private int[] InitializeColoring(int J)
        {
            int[] Cells = new int[J];
            for (int p = 0; p < m_Particles.Count; p++)
            {
                double Hmin = Math.Sqrt(GridData.iGeomCells.GetCellVolume(0));
                double[] ParticlePos = m_Particles[p].Position[0];
                double ParticleAngle = m_Particles[p].Angle[0];
                double[] ParticleScales = m_Particles[p].GetLengthScales();
                double Upperedge = ParticlePos[1] + ParticleScales[1] * Math.Cos(ParticleAngle) + ParticleScales[0] * Math.Sin(ParticleAngle) + Hmin / 4;
                double Loweredge = ParticlePos[1] - ParticleScales[1] * Math.Cos(ParticleAngle) - ParticleScales[0] * Math.Sin(ParticleAngle) - Hmin / 4;
                double Leftedge = ParticlePos[0] - ParticleScales[0] * Math.Cos(ParticleAngle) - ParticleScales[1] * Math.Sin(ParticleAngle) - Hmin / 4;
                double Rightedge = ParticlePos[0] + ParticleScales[0] * Math.Cos(ParticleAngle) + ParticleScales[1] * Math.Sin(ParticleAngle) + Hmin / 4;
                for (int j = 0; j < J; j++)
                {
                    double[] center = GridData.iLogicalCells.GetCenter(j);
                    if (center[0] > Leftedge && center[0] < Rightedge && center[1] > Loweredge && center [1] < Upperedge)
                    {
                        ParticleColor.SetMeanValue(j, p + 1);
                        m_Particles[p].ParticleColoredCells.Add(new int[2] { j, p + 1 });
                        Cells[j] = p + 1;
                    } 
                }
            }
            CheckForNeighborColorsInit(Cells);
            InitializedColor = true;
            return Cells;
        }

        private void CheckForNeighborColorsInit(int[] ColoredCells)
        {
            for (int i = 0; i < ColoredCells.Length; i++)
            {
                if (ColoredCells[i] != 0)
                {
                    GridData.GetCellNeighbours(i, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out int[] ConnectingEntities);
                    for (int j = 0; j < CellNeighbors.Length; j++)
                    {
                        if (CellNeighbors[j] < ColoredCells.Max() && ColoredCells[i] != ColoredCells[j] && ColoredCells[CellNeighbors[j]] != 0)
                        {
                            RecolorCellsInit(ColoredCells, ColoredCells[i], ColoredCells[CellNeighbors[j]]);
                        }
                    }
                }
            }
        }

        private void RecolorCellsInit(int[] ColoredCells, int NewColor, int OldColor)
        {
            for (int i = 0; i < ColoredCells.Length; i++)
            {
                if(ColoredCells[i] == OldColor)
                {
                    ColoredCells[i] = NewColor;
                    ParticleColor.SetMeanValue(i, NewColor);
                }
            }
        }

        private void UpdateForcesAndTorque(double dt, double phystime) {
            //
            // Note on MPI parallelization of particle solver:
            // ===============================================
            //
            // - hydrodynamic forces are computed for each domain and added together;
            //   e.g. in the case of particles at MPI-boundaries each processor computes his part of the integral
            // - collisions are detected globally / collision forces are thus only computed on MPI rank 0
            // - finally, the forces on all particles are summed over all MPI processors and the result is stored on all processors (MPI_Allreduce)
            // - particle motion is computed for all particles simultaneously on all processors; every processor knows every particle
            // - since the particle solver is much cheaper than the flow solver, this "not-really parallel" approach may work up to a few hundreds of particles
            //

            // Initial check: is the motion state of the particles equal on all MPI processors?
            // ================================================================================
            int D = GridData.SpatialDimension;
            int NoOfParticles = m_Particles.Count;

            {
                // verify that we have the same number of particles on each processor
                int NoOfParticles_min = NoOfParticles.MPIMin();
                int NoOfParticles_max = NoOfParticles.MPIMax();
                if (NoOfParticles_max != NoOfParticles || NoOfParticles_max != NoOfParticles)
                    throw new ApplicationException("mismatch in number of MPI particles");

                // nor, compare those particles:
                int NoOfVars = (10 + D * 10); // variables per particle; size can be increased if more values should be compared
                double[] CheckSend = new double[NoOfParticles * NoOfVars]; 

                for (int iP = 0; iP < NoOfParticles; iP++) {
                    var P = m_Particles[iP];

                    // scalar values
                    CheckSend[iP * NoOfVars + 0] = P.Angle[0];
                    CheckSend[iP * NoOfVars + 1] = P.Angle[1];
                    //CheckSend[iP*NoOfVars + 2] = P.Area_P;
                    CheckSend[iP * NoOfVars + 3] = P.ClearSmallValues ? 1.0 : 0.0;
                    CheckSend[iP * NoOfVars + 4] = P.ForceAndTorque_convergence;
                    CheckSend[iP * NoOfVars + 5] = P.Mass_P;
                    CheckSend[iP * NoOfVars + 6] = P.particleDensity;
                    // todo: add more values here that might be relevant for the particle state;

                    // vector values
                    for (int d = 0; d < D; d++) {
                        int Offset = 10;
                        CheckSend[iP * NoOfVars + Offset + 0 * D + d] = P.Position[0][d];
                        CheckSend[iP * NoOfVars + Offset + 1 * D + d] = P.Position[1][d];
                        CheckSend[iP * NoOfVars + Offset + 2 * D + d] = P.TranslationalAcceleration[0][d];
                        CheckSend[iP * NoOfVars + Offset + 3 * D + d] = P.TranslationalAcceleration[1][d];
                        CheckSend[iP * NoOfVars + Offset + 4 * D + d] = P.TranslationalVelocity[0][d];
                        CheckSend[iP * NoOfVars + Offset + 5 * D + d] = P.TranslationalVelocity[1][d];
                        CheckSend[iP * NoOfVars + Offset + 6 * D + d] = P.HydrodynamicForces[0][d];
                        CheckSend[iP * NoOfVars + Offset + 7 * D + d] = P.HydrodynamicForces[1][d];
                        // todo: add more vector values here that might be relevant for the particle state;
                    }
                }

                int MPIsz = MPISize;
                double[] CheckReceive = new double[NoOfParticles * NoOfVars * MPIsz];
                unsafe {
                    fixed(double* pCheckSend = CheckSend, pCheckReceive = CheckReceive) {
                        csMPI.Raw.Allgather((IntPtr)pCheckSend, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                    }
                }

                for (int iP = 0; iP < NoOfParticles; iP++) {
                    for (int iVar = 0; iVar < NoOfVars; iVar++) {
                        // determine a tolerance...
                        int idx_l = 
                            iP * NoOfVars // particle index offset
                           + iVar; // variable index offset
                        double VarTol = Math.Abs(CheckSend[idx_l]) * 1.0e-10;

                        // compare
                        for (int r = 0; r < MPIsz; r++) {

                            int idx_g = CheckSend.Length * r // MPI index offset
                                + idx_l;

                            if (Math.Abs(CheckReceive[idx_g] - CheckSend[idx_l]) > VarTol)
                                throw new ApplicationException("Mismatch in particle state among MPI ranks. Index:  " + idx_l);
                        }
                        VarTol *= 1.0e-10;
                    }
                }
            }

            // Update forces
            // =============
            foreach (Particle p in m_Particles) {

                if (!((FSI_Control)Control).pureDryCollisions) {
                    p.UpdateForcesAndTorque(Velocity, Pressure, LsTrk, Control.PhysicalParameters.mu_A, dt, Control.PhysicalParameters.rho_A);
                }

                // wall collisions are computed on each processor
                WallCollisionForces(p, LsTrk.GridDat.Cells.h_minGlobal);
            }
            if (MPIRank == 0 && m_Particles.Count > 1) {
                // inter-particle collisions are computed only on rank 0
                UpdateCollisionForces(m_Particles, LsTrk.GridDat.Cells.h_minGlobal);
            }

            // Sum forces and moments over all MPI processors
            // ==============================================
            {
                // step 1: collect all variables that we need to sum up
                int NoOfVars = 1 + D * 1;
                double[] StateBuffer = new double[NoOfParticles * NoOfVars];

                for (int iP = 0; iP < NoOfParticles; iP++) {
                    var P = m_Particles[iP];
                    StateBuffer[NoOfVars * iP + 0] = 0;
                    StateBuffer[NoOfVars * iP + 0] = P.HydrodynamicTorque[0];
                    for(int d = 0; d < D; d++) {
                        int Offset = 1;
                        StateBuffer[NoOfVars * iP + Offset + 0 * D + d] = 0;
                        StateBuffer[NoOfVars * iP + Offset + 0*D + d] = P.HydrodynamicForces[0][d];
                    }
                }

                // step 2: sum over MPI processors
                // note: we want to sum all variables by a single MPI call, which is way more efficient
                // B. Deußen: a single call of MPISum() would only consider the first entry of StateBuffer, thus I implemented the loop over all entries
                double[] GlobalStateBuffer = new double[StateBuffer.Length];
                for (int i = 0; i < StateBuffer.Length; i++)
                {
                    GlobalStateBuffer[i] = 0;
                    GlobalStateBuffer[i] = StateBuffer[i].MPISum();
                }

                // step 3: write sum variables back 
                for (int iP = 0; iP < NoOfParticles; iP++) {
                    var P = m_Particles[iP];
                    P.HydrodynamicTorque[0] = 0;
                    P.HydrodynamicTorque[0] = GlobalStateBuffer[NoOfVars * iP + 0];
                    for(int d = 0; d < D; d++) {
                        int Offset = 1;
                        P.HydrodynamicForces[0][d] = 0;
                        P.HydrodynamicForces[0][d] = GlobalStateBuffer[NoOfVars * iP + Offset + 0 * D + d];
                    }
                }
            }
        }

        void ExchangeDampingTensors()
        {
            // Sum forces and moments over all MPI processors
            // ==============================================
            {
                // step 1: collect all variables that we need to sum up
                int NoOfParticles = m_Particles.Count;
                int NoOfVars = 3; //only for 2D at the moment
                double[,] StateBuffer = new double[NoOfParticles * NoOfVars, NoOfParticles * NoOfVars];
                for (int iP = 0; iP < NoOfParticles; iP++)
                {
                    var P = m_Particles[iP];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            StateBuffer[NoOfVars * iP + i, NoOfVars * iP + j] = P.AddedDampingTensor[i, j];
                        }
                    }

                }
                // step 2: sum over MPI processors
                // note: we want to sum all variables by a single MPI call, which is way more efficient
                // B. Deußen: a single call of MPISum() would only consider the first entry of StateBuffer, thus I implemented the loop over all entries
                double[,] GlobalStateBuffer = new double[NoOfParticles * NoOfVars, NoOfParticles * NoOfVars];
                for (int i = 0; i < NoOfParticles * NoOfVars; i++)
                {
                    for (int j = 0; j < NoOfParticles * NoOfVars; j++)
                    {
                        GlobalStateBuffer[i, j] = StateBuffer[i, j].MPISum();
                    }

                }
                // step 3: write sum variables back 
                for (int iP = 0; iP < NoOfParticles; iP++)
                {
                    var P = m_Particles[iP];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            StateBuffer[NoOfVars * iP + i, NoOfVars * iP + j] = P.AddedDampingTensor[i, j];
                            P.AddedDampingTensor[i, j] = GlobalStateBuffer[NoOfVars * iP + i, NoOfVars * iP + j];
                        }
                    }
                }
            }
        }

        void PrintResultToConsole(double phystime, double dt)
        {
            double[] TranslationalMomentum = new double[2] { 0, 0 };
            double RotationalMomentum = 0;
            double[] totalKE = new double[3] { 0, 0, 0 };
            double xPos;
            double yPos;
            double ang;

            foreach (Particle p in m_Particles)
            {
                double[] SingleParticleMomentum = p.CalculateParticleMomentum(dt);
                double[] SingleParticleKineticEnergy = p.CalculateParticleKineticEnergy(dt);
                TranslationalMomentum[0] += SingleParticleMomentum[0];
                TranslationalMomentum[1] += SingleParticleMomentum[1];
                RotationalMomentum += SingleParticleMomentum[SingleParticleMomentum.Length - 1];
                totalKE[0] += SingleParticleKineticEnergy[0];
                totalKE[1] += SingleParticleKineticEnergy[1];
                totalKE[2] += SingleParticleKineticEnergy[SingleParticleMomentum.Length - 1];
            }
            Console.WriteLine("Total momentum in system:  " + Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
            Console.WriteLine("Total kinetic energy in system:  " + (totalKE[0] + totalKE[1] + totalKE[2]));

            force = m_Particles[0].HydrodynamicForces[0];
            torque = m_Particles[0].HydrodynamicTorque[0];
            xPos = m_Particles[0].Position[0][0];
            yPos = m_Particles[0].Position[0][1];
            ang = m_Particles[0].Angle[0];
            var MPItransVelocity = m_Particles[0].TranslationalVelocity[0];
            var MPIangularVelocity = m_Particles[0].RotationalVelocity[0];
            

            //Console.WriteLine(newPosition[1].MPIMax());

            /*
            if ((base.MPIRank == 0) && (Log_DragAndLift != null)) {
                double drag = force[0];
                double lift = force[1];
                //string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", TimestepNo, phystime, m_Particles[0].positionAtIteration[0][0], m_Particles[0].positionAtIteration[0][1], m_Particles[0].particleAnglePerIteration[0], m_Particles[0].transVelocityAtIteration[0][0], m_Particles[0].transVelocityAtIteration[0][1], 0.0, (totalKE[0] + totalKE[1] + totalKE[2]), Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", phystime, m_Particles[0].Position[0][0], m_Particles[0].Position[0][1], m_Particles[0].Angle[0], m_Particles[0].TranslationalVelocity[0][0], m_Particles[0].TranslationalVelocity[0][1], 0.0, (totalKE[0] + totalKE[1] + totalKE[2]), Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
                Log_DragAndLift.WriteLine(line);
                Log_DragAndLift.Flush();
            }
            */

            // Save for NUnit Test
            base.QueryHandler.ValueQuery("C_Drag", 2 * force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
            base.QueryHandler.ValueQuery("C_Lift", 2 * force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
            base.QueryHandler.ValueQuery("Angular_Velocity", MPIangularVelocity, true); // (TestCase FlowRotationalCoupling)


            Console.WriteLine("Drag Force:   {0}", m_Particles[0].HydrodynamicForces[0][0]);
            Console.WriteLine("Lift Force:   {0}", m_Particles[0].HydrodynamicForces[0][1]);
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

        // <summary>
        // Variables for FSI coupling
        // </summary>
        //double MPIangularVelocity;
        //readonly double[] newTransVelocity = new double[2];
        //readonly double[] oldPosition = new double[2];
        //readonly double[] newPosition = new double[2];
        //readonly double[] oldforce = new double[2];
        //double[] MPItransVelocity = new double[2];

        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {
                
                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;

                base.ResLogger.TimeStep = TimestepInt;

                hack_phystime = phystime;
                dt = base.GetFixedTimestep();

                Console.WriteLine("Starting time-step " + TimestepInt + "...");

                int OldPushCount = LsTrk.PushCount; // used later to check if there is exactly one push per timestep
               

                if (((FSI_Control)this.Control).pureDryCollisions) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // only particle motion & collisions, no flow solver
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    // in other branches, called by the BDF timestepper
                    LsTrk.PushStacks();
                    DGLevSet.Push();
                    
                    UpdateForcesAndTorque(dt, phystime);
                    foreach (var p in m_Particles) {

                        p.CalculateAcceleration(dt, Control.PhysicalParameters.rho_A);

                        p.CalculateTranslationalVelocity(dt, Control.PhysicalParameters.rho_A, ((FSI_Control)this.Control).includeTranslation);
                        p.CalculateAngularVelocity(dt, ((FSI_Control)this.Control).includeRotation);

                        p.CalculateParticlePosition(dt, Control.PhysicalParameters.rho_A);
                        p.CalculateParticleAngle(dt);
                    }
                    UpdateLevelSetParticles(dt);
                    PrintResultToConsole(phystime, dt);

                } else {
                    if (triggerOnlyCollisionProcedure) {
                        UpdateLevelSetParticles(dt);
                        triggerOnlyCollisionProcedure = false;
                        /*
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
                        */

                        return dt;
                    } else if (((FSI_Control)this.Control).Timestepper_LevelSetHandling != LevelSetHandling.Coupled_Iterative) {
                        int iteration_counter = 0;
                        double posResidual_splitting = 1e12;
                        while (posResidual_splitting > ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion)
                        {
                            double[] ForcesOldSquared = new double[2];
                            double TorqueOldSquared = new double();
                            ForcesOldSquared[0] = 0;
                            ForcesOldSquared[1] = 0;
                            TorqueOldSquared = 0;

                            if (iteration_counter == 0 && ((FSI_Control)this.Control).splitting_fully_coupled == true)
                            {
                                foreach (Particle p in m_Particles)
                                {
                                    p.iteration_counter_P = iteration_counter;
                                    if (p.neglectAddedDamping == false && p.iteration_counter_P == 0)
                                    {
                                        p.UpdateDampingTensors();
                                        ExchangeDampingTensors();
                                    }
                                    p.PredictAcceleration();
                                    p.CalculateAngularVelocity(dt, ((FSI_Control)this.Control).includeRotation);
                                    p.CalculateTranslationalVelocity(dt, this.Control.PhysicalParameters.rho_A, ((FSI_Control)this.Control).includeTranslation);
                                    p.CalculateParticlePosition(dt, this.Control.PhysicalParameters.rho_A);
                                    p.CalculateParticleAngle(dt);
                                    p.ComputeParticleRe(this.Control.PhysicalParameters.mu_A);
                                }
                                posResidual_splitting = 1e12;
                            }
                            else
                            {
                                foreach (Particle p in m_Particles)
                                {
                                    p.iteration_counter_P = iteration_counter;
                                    p.ForceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                                    ForcesOldSquared[0] += p.HydrodynamicForces[0][0].Pow2();
                                    ForcesOldSquared[1] += p.HydrodynamicForces[0][1].Pow2();
                                    TorqueOldSquared += p.HydrodynamicTorque[0].Pow2();
                                    p.HydrodynamicForces[0][0] = 0;
                                    p.HydrodynamicForces[0][1] = 0;
                                    p.HydrodynamicTorque[0] = 0;
                                }

                                m_BDF_Timestepper.Solve(phystime, dt, false);

                                UpdateForcesAndTorque(dt, phystime);

                                foreach (Particle p in m_Particles)
                                {
                                    if (iteration_counter == 100)
                                    {
                                        p.PredictAccelerationWithinIteration();
                                    }
                                    else
                                    {
                                        p.CalculateAcceleration(dt, Control.PhysicalParameters.rho_A);
                                    }
                                    p.CalculateAngularVelocity(dt, ((FSI_Control)this.Control).includeRotation);
                                    p.CalculateTranslationalVelocity(dt, this.Control.PhysicalParameters.rho_A, ((FSI_Control)this.Control).includeTranslation);
                                    p.CalculateParticlePosition(dt, this.Control.PhysicalParameters.rho_A);
                                    p.CalculateParticleAngle(dt);
                                    p.ComputeParticleRe(this.Control.PhysicalParameters.mu_A);
                                }
                                double[] ForcesNewSquared = new double[2];
                                double TorqueNewSquared = new double();
                                foreach (Particle p in m_Particles)
                                {
                                    ForcesNewSquared[0] += p.HydrodynamicForces[0][0].Pow2();
                                    ForcesNewSquared[1] += p.HydrodynamicForces[0][1].Pow2();
                                    TorqueNewSquared += p.HydrodynamicTorque[0].Pow2();
                                }
                                posResidual_splitting = Math.Sqrt((Math.Sqrt(ForcesNewSquared[0]) - Math.Sqrt(ForcesOldSquared[0])).Pow2() + (Math.Sqrt(ForcesNewSquared[1]) - Math.Sqrt(ForcesOldSquared[1])).Pow2() + (Math.Sqrt(TorqueNewSquared) - Math.Sqrt(TorqueOldSquared)).Pow2());
                            }
                            PrintResultToConsole(phystime, dt);
                            //#region Get Drag and Lift Coefficiant
                            int PrintIteration = iteration_counter + 1;
                            Console.WriteLine("Fully coupled system, number of iterations:  " + PrintIteration);
                            Console.WriteLine("Forces and torque residual: " + posResidual_splitting);
                            Console.WriteLine();
                            iteration_counter += 1;
                            if (((FSI_Control)this.Control).splitting_fully_coupled == false) {
                                break;
                            }
                            if (iteration_counter > ((FSI_Control)this.Control).max_iterations_fully_coupled) {
                                throw new ApplicationException("no convergence in coupled iterative solver, number of iterations: " + iteration_counter);
                            }
                        }
                        if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled)
                        {
                            LsTrk.IncreaseHistoryLength(1);
                            LsTrk.PushStacks();
                        }
                        /*
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
                        */
                    }
                    else {
                        foreach (Particle p in m_Particles) {
                            p.iteration_counter_P = -1;
                            p.ForceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                        }
                        m_BDF_Timestepper.Solve(phystime, dt, false);
                    }
                }

                // finalize
                // ========
                if (LsTrk.PushCount - OldPushCount != 1) {
                    // To whom it concerns / who stumbles across this exception:
                    // It is important that LevelSetTracker.PushStacks() is called *exactly once per time-step*, at the beginning.
                    // Do not remove this check! Instead, remove any calls to 'PushStacks()' in subroutines of this method.
                    // Fk.
                    throw new ApplicationException("Illegal number of level-set push actions in time-step." + (LsTrk.PushCount - OldPushCount));
                }

                this.ResLogger.NextTimestep(false);

                Console.WriteLine("done with time-step.");
                return dt;
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
                } else {
                    return new BsonReader(s);
                }
            }

            JsonWriter GetJsonWriter(Stream s) {
                if (DebugSerialization) {
                    return new JsonTextWriter(new StreamWriter(s));
                } else {
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

            //Console.WriteLine(o.ToString());

            Debug.Assert(b.Particles.Length == o.Particles.Length);
            int L = b.Particles.Length;
            for (int l = 0; l < L; l++) { // loop over particles
                Debug.Assert(GenericBlas.L2Dist(b.Particles[l].Position[0], o.Particles[l].Position[0]) < 1e-13);
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
            hack_phystime = ArschInfo.PhysicalTime;
            UpdateLevelSetParticles(0.0);

            // call base shit
            var R = base.RestartFromDatabase(out time);



            foreach (Particle p in m_Particles) {
                p.m_collidedWithParticle = new bool[m_Particles.Count];
                p.m_collidedWithWall = new bool[4];
                p.m_closeInterfacePointTo = new double[m_Particles.Count][];

            }

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

        // Initialize particles
        protected override void SetInitial() {
            // Setup particles
            m_Particles = ((FSI_Control)this.Control).Particles;
            hack_phystime = 0.0;
            UpdateLevelSetParticles(0.0);

            // call base implementation
            base.SetInitial();

            foreach (Particle p in m_Particles) {
                p.m_collidedWithParticle = new bool[m_Particles.Count];
                p.m_collidedWithWall = new bool[4];
                p.m_closeInterfacePointTo = new double[m_Particles.Count][];

            }
        }

        public IList<Particle> Particles {
            get {
                return m_Particles;
            }
        }

        List<Particle> m_Particles;
        bool collision = false;

        bool triggerOnlyCollisionProcedure = false;

        // Collision models
        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
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
                    
                    var particle0CutCells = particle0.CutCells_P(LsTrk);
                    var particle1CutCells = particle1.CutCells_P(LsTrk);

                    //var neighborCellsArray_P0 = particle0CutCells.AllNeighbourCells().ItemEnum.ToArray();
                    //var allCells_P0 = new CellMask(GridData, neighborCellsArray_P0);
                    var allCells_P0 = particle0CutCells.AllNeighbourCells();

                    //var neighborCellsArray_P1 = particle1CutCells.AllNeighbourCells().ItemEnum.ToArray();
                    //var allCells_P1 = new CellMask(GridData, neighborCellsArray_P1);
                    var allCells_P1 = particle1CutCells.AllNeighbourCells();

                    double distance = 1E20;
                    double[] distanceVec = new double[Grid.SpatialDimension];

                    var interSecMask = allCells_P0.Intersect(allCells_P1);

                    var p0intersect = interSecMask.AllNeighbourCells().Intersect(particle0CutCells);
                    var p1intersect = interSecMask.AllNeighbourCells().Intersect(particle1CutCells);

                    // If there is no element neighbor of both particle cut cells return
                    if (!interSecMask.IsEmpty) {
                        ComputeCollissionModel(hmin, particle0, particle1, ref distance, ref distanceVec);
                    }
                }
            }
        }

        /// <summary>
        /// Update of particle state (velocity, force, etc.) for two particles where a collision is detected
        /// </summary>
        private void ComputeCollissionModel(double hmin, Particle particle0, Particle particle1, ref double distance, ref double[] distanceVec) {
            // All interface points at a specific subgrid containing all cut cells of one particle
            //var interfacePoints_P0 = XNSEUtils.GetInterfacePoints(LsTrk, LevSet, new SubGrid(particle0CutCells));
            //var interfacePoints_P1 = XNSEUtils.GetInterfacePoints(LsTrk, LevSet, new SubGrid(particle1CutCells));
            MultidimensionalArray interfacePoints_P0 = particle0.GetSurfacePoints(LsTrk, LevSet);
            MultidimensionalArray interfacePoints_P1 = particle1.GetSurfacePoints(LsTrk, LevSet);

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


                        collisionForce = distanceVec;
                        var collisionForceP1 = collisionForce.CloneAs();
                        collisionForce.ScaleV(-100.0);
                        collisionForceP1.ScaleV(-100.0);
                        particle0.HydrodynamicForces[0].AccV(-1, collisionForce);
                        //particle0.hydrodynTorqueAtIteration[0] += 100 * (collisionForce[0] * (tempPoint_P0[0] - particle0.positionAtIteration[0][0]) + collisionForce[1] * (tempPoint_P0[1] - particle0.positionAtIteration[0][1]));
                        particle1.HydrodynamicForces[0].AccV(1, collisionForceP1);
                        //particle1.hydrodynTorqueAtIteration[0] += -100 * (collisionForceP1[0] * (tempPoint_P1[0] - particle1.positionAtIteration[0][0]) + collisionForceP1[1] * (tempPoint_P1[1] - particle1.positionAtIteration[0][1]));
                        Console.WriteLine("Collision information: Particles coming close, force " + collisionForce.L2Norm());
                        Console.WriteLine("Collision information: Particles coming close, torque " + particle1.HydrodynamicTorque[0]);

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
                        double collisionVn_P0 = particle0.TranslationalVelocity[0][0] * normal[0] + particle0.TranslationalVelocity[0][1] * normal[1];
                        double collisionVt_P0 = particle0.TranslationalVelocity[0][0] * tangential[0] + particle0.TranslationalVelocity[0][1] * tangential[1];
                        double collisionVn_P1 = particle1.TranslationalVelocity[0][0] * normal[0] + particle1.TranslationalVelocity[0][1] * normal[1];
                        double collisionVt_P1 = particle1.TranslationalVelocity[0][0] * tangential[0] + particle1.TranslationalVelocity[0][1] * tangential[1];

                        // exzentric collision
                        // ----------------------------------------                                                                  
                        tempPoint_P0.AccV(-1, particle0.Position[0]);
                        double a0 = (tempPoint_P0[0] * tangential[0] + tempPoint_P0[1] * tangential[1]);
                        tempPoint_P1.AccV(-1, particle1.Position[0]);
                        double a1 = (tempPoint_P1[0] * tangential[0] + tempPoint_P1[1] * tangential[1]);




                        // Fix for Sphere
                        // ----------------------------------------  
                        if (particle0 is Particle_Sphere)
                            a0 = 0.0;
                        if (particle1 is Particle_Sphere)
                            a1 = 0.0;



                        double Fx = (1 + e) * ((collisionVn_P0 - collisionVn_P1) / (1 / particle0.Mass_P + 1 / particle1.Mass_P + a0.Pow2() / particle0.MomentOfInertia_P + a1.Pow2() / particle1.MomentOfInertia_P));
                        double Fxrot = (1 + e) * ((-a0 * particle0.RotationalVelocity[0] + a1 * particle1.RotationalVelocity[0]) / (1 / particle0.Mass_P + 1 / particle1.Mass_P + a0.Pow2() / particle0.MomentOfInertia_P + a1.Pow2() / particle1.MomentOfInertia_P));

                        double tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / particle0.Mass_P;
                        double tempCollisionVn_P1 = collisionVn_P1 + (Fx + Fxrot) / particle1.Mass_P;
                        double tempCollisionVt_P0 = collisionVt_P0;
                        double tempCollisionVt_P1 = collisionVt_P1;
                        Console.WriteLine("a0:    " + a0 + "   Fx:    " + (-Fx) + "      Fxrot:    " + (-Fxrot));
                        Console.WriteLine("a1:    " + a1 + "   Fx:    " + Fx + "      Fxrot:    " + Fxrot);
                        particle0.RotationalVelocity[0] = particle0.RotationalVelocity[0] + a0 * (Fx + Fxrot) / particle0.MomentOfInertia_P;
                        particle1.RotationalVelocity[0] = particle1.RotationalVelocity[0] - a1 * (Fx + Fxrot) / particle1.MomentOfInertia_P;



                        // zentric collision
                        // ----------------------------------------
                        //double tempCollisionVn_P0 = (particle0.mass_P * collisionVn_P0 + particle1.mass_P * collisionVn_P1 + e * particle1.mass_P * (collisionVn_P1 - collisionVn_P0)) / (particle0.mass_P + particle1.mass_P);
                        //double tempCollisionVt_P0 = collisionVt_P0;
                        //double tempCollisionVn_P1 = (particle0.mass_P * collisionVn_P0 + particle1.mass_P * collisionVn_P1 + e * particle0.mass_P * (collisionVn_P0 - collisionVn_P1)) / (particle0.mass_P + particle1.mass_P);
                        //double tempCollisionVt_P1 = collisionVt_P1;
                        // ----------------------------------------



                        particle0.TranslationalVelocity[0] = new double[] { normal[0] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[0], normal[1] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[1] };
                        particle1.TranslationalVelocity[0] = new double[] { normal[0] * tempCollisionVn_P1 + tempCollisionVt_P1 * tangential[0], normal[1] * tempCollisionVn_P1 + tempCollisionVt_P1 * tangential[1] };
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

        private FSI_Solver.FSI_Control.CollisionModel m_collisionModel {
            get {
                return ((FSI_Control)Control).collisionModel;
            }
        }

        /// <summary>
        /// Calculation of collision forces between particle and wall
        /// </summary>
        public void WallCollisionForces(Particle particle, double hmin) {

            if (m_collisionModel == FSI_Control.CollisionModel.NoCollisionModel) {
                return;
            }

            var particleCutCells = particle.CutCells_P(LsTrk);

            //var particleCutCellArray = particleCutCells.ItemEnum.ToArray();
            //var neighborCellsArray = particleCutCells.AllNeighbourCells().ItemEnum.ToArray();
            //var allCellsArray = particleCutCellArray.Concat(neighborCellsArray).ToArray();
            //var allCells = new CellMask(GridData, neighborCellsArray);
            var allCells = particleCutCells;

            collision = false;

            double distance = double.MaxValue;
            double[] distanceVec = new double[Grid.SpatialDimension];

            // All interface points at a specific subgrid containing all cut cells of one particle
            MultidimensionalArray interfacePoints = null;

            //Console.WriteLine("ParticleCutCellCount:   " + particleCutCells.Count());

            var trafo = GridData.iGeomEdges.Edge2CellTrafos;

            SubGrid allCellsGrid = new SubGrid(allCells);

            double[] tempPoint = new double[2] { 0.0, 0.0 };

            foreach (int iEdge in allCellsGrid.BoundaryEdgesMask.ItemEnum) {

                // Collision forces have to act
                if (GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge)) {

                    if (interfacePoints == null)
                        interfacePoints = particle.GetSurfacePoints(LsTrk, LevSet);

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


            Console.WriteLine("Closest Distance to wall is: " + distance);

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

                        particle.HydrodynamicForces[0] = collisionForce;


                        return;
                    }


                    if (realDistance <= (1.5 * hmin)) {

                        distanceVec.ScaleV((threshold - realDistance).Abs());
                        distanceVec.ScaleV(1 / epsPrime);
                        collisionForce = distanceVec;

                        collisionForce.ScaleV(100.0);
                        particle.HydrodynamicForces[0].AccV(1, collisionForce);
                        particle.HydrodynamicTorque[0] -= (collisionForce[0] * (tempPoint[0] - particle.Position[0][0]) + collisionForce[1] * (tempPoint[1] - particle.Position[0][1]));
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
                        if (particle.Position[0][1] < 0.5)
                            e = 0.0;


                        // if particle already collided with wall
                        particle.m_collidedWithWall[0] = true;

                        // Skip force integration for next timestep
                        particle.skipForceIntegration = false;

                        //collision Nomal
                        var normal = distanceVec.CloneAs();
                        normal.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));
                        double[] tangential = new double[] { -normal[1], normal[0] };


                        double collisionVn_P0 = particle.TranslationalVelocity[0][0] * normal[0] + particle.TranslationalVelocity[0][1] * normal[1];
                        double collisionVt_P0 = particle.TranslationalVelocity[0][0] * tangential[0] + particle.TranslationalVelocity[0][1] * tangential[1];


                        // exzentric collision
                        // ----------------------------------------
                        tempPoint.AccV(-1, particle.Position[0]);
                        double a0 = (tempPoint[0] * tangential[0] + tempPoint[1] * tangential[1]);

                        if (particle is Particle_Sphere)
                            a0 = 0.0;


                        double Fx = (1 + e) * (collisionVn_P0) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);
                        double Fxrot = (1 + e) * (-a0 * particle.RotationalVelocity[0]) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);

                        double tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / particle.Mass_P;
                        double tempCollisionVt_P0 = collisionVt_P0;

                        particle.RotationalVelocity[0] = particle.RotationalVelocity[0] + a0 * (Fx + Fxrot) / particle.MomentOfInertia_P;


                        particle.TranslationalVelocity[0] = new double[] { normal[0] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[0], normal[1] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[1] };

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


        /// <summary>
        /// Mesh refinement
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
                //var CutCellArray = CutCells.ItemEnum.ToArray();
                //var CutCellNeighborsArray = CutCellNeighbors.ItemEnum.ToArray();
                //var AllCells = CutCellArray.Concat(CutCellNeighborsArray).ToArray();
                //var NoCoarseningcells = new CellMask(this.GridData, AllCells);
                //var AllCells = CutCells.Union(CutCellNeighbors).ItemEnum.ToArray();

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
    }
}