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

namespace BoSSS.Application.FSI_Solver
{
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain
    {
        /// <summary>
        /// Set the inital state of the simulation.
        /// </summary>
        protected override void SetInitial()
        {
            // Setup particles
            m_Particles = ((FSI_Control)this.Control).Particles;
            hack_phystime = 0.0;
            UpdateLevelSetParticles();

            // call base implementation
            base.SetInitial();
        }

        public IList<Particle> Particles
        {
            get
            {
                return m_Particles;
            }
        }

        List<Particle> m_Particles;

        private FSI_Solver.FSI_Control.CollisionModel CollisionModel
        {
            get
            {
                return ((FSI_Control)Control).collisionModel;
            }
        }

        static int counter = 0;

        public static void AgglomerationFailDebugPlot(DGField[] f)
        {
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);


            Tecplot.PlotFields(f, "AgglomerationFailDebugPlot-" + counter, 0.0, 2);
            counter++;
        }

        bool CalculatedDampingTensors = false;
        readonly private FSI_Auxillary Auxillary = new FSI_Auxillary();
        readonly private FSI_LevelSetUpdate LevelSetUpdate = new FSI_LevelSetUpdate();

        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args)
        {
            //MultiphaseCellAgglomerator.Katastrophenplot = AgglomerationFailDebugPlot;
            //TestProgram.Init();
            //BoSSS.Application.FSI_Solver.TestProgram.SingleDryParticleAgainstWall(true);

            //Assert.IsTrue(false, "Remember to remove testcode!");

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

        private int IterationCounter = 0;

        protected override void CreateFields()
        {
            base.CreateFields();

            ParticleColor = new SinglePhaseField(new Basis(this.GridData, 0), "ParticleColor");
            m_RegisteredFields.Add(ParticleColor);
            m_IOFields.Add(ParticleColor);

            LevelSetDistance = new SinglePhaseField(new Basis(this.GridData, 0), "LevelSetDistance");
            m_RegisteredFields.Add(LevelSetDistance);
            m_IOFields.Add(LevelSetDistance);
        }

        // Create equations and solvers
        bool UseMovingMesh
        {
            get
            {
                switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling)
                {
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

        bool FullyCoupled
        {
            get
            {
                if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled)
                    return true;
                else
                    return false;
            }
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        {

            if (IBM_Op != null)
                return;
            string[] CodNameSelected = new string[0];
            string[] DomNameSelected = new string[0];

            int D = this.GridData.SpatialDimension;

            BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

            int degU = this.Velocity[0].Basis.Degree;
            NSEOperatorConfiguration IBM_Op_config = new NSEOperatorConfiguration
            {
                convection = this.Control.PhysicalParameters.IncludeConvection,
                continuity = true,
                Viscous = true,
                PressureGradient = true,
                Transport = true,
                CodBlocks = new bool[] { true, true },
                DomBlocks = new bool[] { true, true },
            };

            string[] CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D)).Cat("div");
            string[] Params = ArrayTools.Cat(
                 VariableNames.Velocity0Vector(D),
                 VariableNames.Velocity0MeanVector(D));
            string[] DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

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

                                        double[] result = new double[X.Length + 5];

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
                                            else { containsParticle = p.Contains(X, GridData.iGeomCells.h_min.Min()); }//, LsTrk

                                            FSI_Collision _FSI_Collision = new FSI_Collision();
                                            _FSI_Collision.CalculateRadialVector(p.Position[0], X, out _, out double RadialLength, out double[] RadialNormalVector);
                                            // active particles
                                            if (containsParticle && p.ActiveParticle == true)
                                            {
                                                result[0] = p.TranslationalVelocity[0][0];
                                                result[1] = p.TranslationalVelocity[0][1];
                                                result[2] = p.RotationalVelocity[0];
                                                result[3] = RadialNormalVector[0];
                                                result[4] = RadialNormalVector[1];
                                                result[5] = p.Position[0].L2Distance(X); //RadialLength;
                                                result[6] = -cos_theta;
                                                return result;
                                            }

                                            // passive particles
                                            else if (containsParticle && p.ActiveParticle == false)
                                            {
                                                result[0] = p.TranslationalVelocity[0][0];
                                                result[1] = p.TranslationalVelocity[0][1];
                                                result[2] = p.RotationalVelocity[0];
                                                result[3] = RadialNormalVector[0];
                                                result[4] = RadialNormalVector[1];
                                                result[5] = p.Position[0].L2Distance(X); //RadialLength;
                                                result[6] = 0;
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


                        swipViscosity_Term1 Visc = new swipViscosity_Term1(penalty, d, D, BcMap, ViscosityOption.ConstantViscosity, this.Control.PhysicalParameters.mu_A,// / this.Control.PhysicalParameters.rho_A,
                            double.NaN, null);

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
                            Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, D, LsTrk,
                                penalty, this.ComputePenaltyIB,
                                this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                                delegate (double[] X, double time)
                                {

                                    double[] result = new double[X.Length + 7];

                                    foreach (Particle p in m_Particles)
                                    {
                                        // which particle?
                                        bool containsParticle;
                                        if (m_Particles.Count == 1)
                                        {
                                            containsParticle = true;
                                        }
                                        else { containsParticle = p.Contains(X, GridData.iGeomCells.h_min.Min()); }//, LsTrk

                                        FSI_Collision _FSI_Collision = new FSI_Collision();
                                        _FSI_Collision.CalculateRadialVector(p.Position[0], X, out _, out double RadialLength, out double[] RadialNormalVector);
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
                                            result[3] = RadialNormalVector[0];
                                            result[4] = RadialNormalVector[1];
                                            result[5] = p.Position[0].L2Distance(X);
                                            result[6] = p.ActiveStress;
                                            result[7] = -cos_theta;
                                            result[8] = p.Angle[0];
                                        }

                                        // passive particles
                                        else if (containsParticle && p.ActiveParticle == false)
                                        {
                                            result[0] = p.TranslationalVelocity[0][0];
                                            result[1] = p.TranslationalVelocity[0][1];
                                            result[2] = p.RotationalVelocity[0];
                                            result[3] = RadialNormalVector[0];
                                            result[4] = RadialNormalVector[1];
                                            result[5] = p.Position[0].L2Distance(X);
                                            result[6] = 0;
                                            result[7] = 0;
                                            result[8] = p.Angle[0];
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
                        var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.ActiveDivergenceAtIB(D, LsTrk, 1,
                           delegate (double[] X, double time)
                           {

                               double[] result = new double[X.Length + 4];

                               foreach (Particle p in m_Particles)
                               {
                                   bool containsParticle;
                                   if (m_Particles.Count == 1)
                                   {
                                       containsParticle = true;
                                   }
                                   else { containsParticle = p.Contains(X, GridData.iGeomCells.h_min.Min()); }
                                   FSI_Collision _FSI_Collision = new FSI_Collision();
                                   _FSI_Collision.CalculateRadialVector(p.Position[0], X, out _, out double RadialLength, out double[] RadialNormalVector);
                                   if (containsParticle)
                                   {
                                       result[0] = p.TranslationalVelocity[0][0];
                                       result[1] = p.TranslationalVelocity[0][1];
                                       result[2] = p.RotationalVelocity[0];
                                       result[3] = RadialNormalVector[0];
                                       result[4] = RadialNormalVector[1];
                                       result[5] = p.Position[0].L2Distance(X);
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
            if (this.Control.PhysicalParameters.IncludeConvection)
            {
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
            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling)
            {
                case LevelSetHandling.Coupled_Once:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeDependent;
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    break;

                case LevelSetHandling.LieSplitting:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    if (!CalculatedDampingTensors)
                    {
                        foreach (Particle p in m_Particles)
                        {
                            if (p.UseAddedDamping)
                            {
                                p.CalculateDampingTensor(LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                                Auxillary.ExchangeDampingTensors(m_Particles);
                            }
                        }
                    }
                    CalculatedDampingTensors = true;
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    if (!CalculatedDampingTensors)
                    {
                        foreach (Particle p in m_Particles)
                        {
                            if (p.UseAddedDamping)
                            {
                                p.CalculateDampingTensor(LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
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

        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental)
        {

            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling)
            {
                case LevelSetHandling.None:
                    ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                    //newTransVelocity[0] = (((FSI_Control)this.Control).transVelocityFunc[0])(phystime);
                    //newTransVelocity[1] = (((FSI_Control)this.Control).transVelocityFunc[1])(phystime);
                    LevSet.ProjectField(Posfunction);
                    LsTrk.UpdateTracker();
                    break;

                case LevelSetHandling.Coupled_Once:
                    UpdateLevelSetParticles();
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    Console.WriteLine("WARNING: Coupled iterative solver is not tested!");
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                    CalculateHydrodynamicForces(m_Particles, dt);
                    UpdateLevelSetParticles();
                    foreach (Particle p in m_Particles)
                    {
                        p.iteration_counter_P += 1;
                        p.forceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                    }
                    break;

                case LevelSetHandling.LieSplitting:
                    UpdateLevelSetParticles();
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    UpdateLevelSetParticles();
                    break;

                case LevelSetHandling.StrangSplitting:
                    UpdateLevelSetParticles();
                    break;

                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
            }

            // Forces and Torque residual
            /// <summary>
            /// Computes the Residual of the forces and torque acting from to fluid to the particle.
            /// </summary>
            double forces_PResidual;
            if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative)
            {
                double acc_force_P_x = 0;
                double acc_force_P_y = 0;
                double acc_torque_P = 0;
                double acc_force_P_x_old = 0;
                double acc_force_P_y_old = 0;
                double acc_torque_P_old = 0;
                double iterationCounter = -1;
                foreach (Particle p in m_Particles)
                {
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
                if (iterationCounter == 0)
                {
                    forces_PResidual = 1;
                }
                // compute residual
                else
                {
                    forces_PResidual = Math.Sqrt((acc_force_P_x_old - acc_force_P_x).Pow2() + (acc_force_P_y_old - acc_force_P_y).Pow2() + (acc_torque_P_old - acc_torque_P).Pow2());
                }
                Console.WriteLine("Current forces_PResidual:   " + forces_PResidual);
            }
            // no iterative solver, no residual
            else
            {
                forces_PResidual = 0;
            }
            return forces_PResidual;
        }

        /// <summary>
        /// Array of all local cells with their specific color.
        /// </summary>
        private int[] CellColor = null;
        private CellMask ColoredCellMask = null;


        /// <summary>
        /// Particle to Level-Set-Field 
        /// </summary>
        private void UpdateLevelSetParticles()
        {
            // =======================================================
            // Step 1
            // Define an array with the respective cell colors
            // =======================================================
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            //CellColor = ((FSI_Control)Control).AdaptiveMeshRefinement ? InitializeColoring(J, GridData, ((FSI_Control)Control).AdaptiveMeshRefinement) : CellColor ?? InitializeColoring(J, GridData, ((FSI_Control)Control).AdaptiveMeshRefinement);
            CellColor = CellColor == null ? InitializeColoring(J, GridData, ((FSI_Control)Control).AdaptiveMeshRefinement) : LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")]; 

            // =======================================================
            // Step 2
            // Delete the old level set
            // =======================================================
            DGLevSet.Current.Clear();

            // =======================================================
            // Step 3
            // Define level set per color
            // =======================================================
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate();
            CellMask AgglParticleMask = null;
            List<int[]> ColoredCellsSorted = levelSetUpdate.ColoredCellsFindAndSort(CellColor);
            int[] ParticleColorArray = levelSetUpdate.FindParticleColor(GridData, m_Particles, ColoredCellsSorted);
            int NoOfParticles = ParticleColorArray.Length;
            int[] GlobalParticleColor = new int[NoOfParticles];
            double[] StateBuffer = new double[NoOfParticles];
            for (int i = 0; i < NoOfParticles; i++)
            {
                StateBuffer[i] = Convert.ToDouble(ParticleColorArray[i]);
            }
            double[] GlobalStateBuffer = StateBuffer.MPIMax();
            for (int i = 0; i < NoOfParticles; i++)
            {
                GlobalParticleColor[i] = Convert.ToInt32(GlobalStateBuffer[i]);
            }
            for (int p = 0; p < GlobalParticleColor.Length; p++)
            {
                int CurrentColor = GlobalParticleColor[p];
                bool ContainsCurrentColor = false;
                BitArray ColoredCells = new BitArray(J);
                List<int> ColoredCell_P = new List<int>();
                for (int j = 0; j < J; j++)
                {
                    if (CellColor[j] == CurrentColor && CurrentColor != 0)
                    {
                        ContainsCurrentColor = true;
                        ColoredCells[j] = true;
                        ColoredCell_P.Add(j);
                    }
                }
                if (ContainsCurrentColor)
                {
                    int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(GlobalParticleColor, CurrentColor);
                    //CellMask ColoredCellMask = levelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, CurrentColor, J, false);
                    ColoredCellMask = new CellMask(GridData, ColoredCells);
                    ColoredCellMask.Union(ColoredCellMask.AllNeighbourCells());

                    // Save all colored cells (of any color) in one mask
                    // =================================================
                    AgglParticleMask = AgglParticleMask == null ? ColoredCellMask : AgglParticleMask.Union(ColoredCellMask);

                    // Get particle level set
                    // ======================
                    double phiComplete(double[] X, double t)
                    {
                        // Generating the correct sign
                        // ===========================
                        double phi = Math.Pow(-1, ParticlesOfCurrentColor.Length - 1);
                        // Multiplication over all particle-level-sets within the current color
                        // ====================================================================
                        for (int pc = 0; pc < ParticlesOfCurrentColor.Length; pc++)
                        {
                            Particle Particle0 = m_Particles[ParticlesOfCurrentColor[pc]];
                            phi *= Particle0.Phi_P(X);
                            //phi = Math.Max(Particle0.Phi_P(X), phi);
                            Particle0.ParticleColor = CurrentColor;
                            Particle0.ParticleColoredCells = ColoredCell_P.ToArray();
                            // Delete all particles within the current color from the particle color array
                            // ===========================================================================
                            GlobalParticleColor[ParticlesOfCurrentColor[pc]] = 0;
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
            // Smoothing
            // =======================================================
            PerformLevelSetSmoothing(AgglParticleMask, FluidCells, true);

            // =======================================================
            // Step 6
            // Update level set tracker and coloring
            // =======================================================
            LsTrk.UpdateTracker(__NearRegionWith: 2);
            CellColor = UpdateColoring();
        }

        /// <summary>
        /// Set level set based on the function phi and the current cells
        /// </summary>
        private void SetLevelSet(Func<double[], double, double> phi, CellMask CurrentCells, double phystime)
        {
            ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(phi, phystime);
            DGLevSet.Current.Clear(CurrentCells);
            DGLevSet.Current.ProjectField(1.0, Function, new CellQuadratureScheme(UseDefaultFactories: true, domain: CurrentCells));
        }

        /// <summary>
        /// Update of <see cref="ParticleColor"/> and <see cref="LevelSetDistance"/>
        /// </summary>
        private int[] UpdateColoring()
        {
           //Debugger.Launch();
            // =======================================================
            // Step 1
            // Color all cells directlly related to the level set
            // =======================================================
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int Je = J + GridData.iLogicalCells.NoOfExternalCells;
            ushort[] rCode = LsTrk.Regions.RegionsCode;
            int[] PartCol = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int[] PartColEx = new int[Je];

            for (int j = 0; j < J; j++)
            {
                ParticleColor.SetMeanValue(j, PartCol[j]);
                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(rCode[j], 0));
                PartColEx[j] = PartCol[j];
            }
            PartColEx.MPIExchange(GridData);

            // =======================================================
            // Step 2
            // Color neighbour cells
            // =======================================================
            for (int j = 0; j < J; j++)
            {
                GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                for (int i = 0; i < CellNeighbors.Length; i++)
                {
                    if (PartColEx[CellNeighbors[i]] != 0 && PartColEx[j] == 0)
                    {
                        ParticleColor.SetMeanValue(j, PartColEx[CellNeighbors[i]]);
                        LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(rCode[j], 0));
                        PartCol[j] = PartColEx[CellNeighbors[i]];
                    }
                }
            }
            // rewrite to communication array
            for (int j = 0; j < J; j++)
            {
                PartColEx[j] = PartCol[j];
            }
            PartColEx.MPIExchange(GridData);

            List<int[]> ColoredCellsSorted2 = LevelSetUpdate.ColoredCellsFindAndSort(PartColEx);
            int[] ParticleColorArray = LevelSetUpdate.FindParticleColor(GridData, m_Particles, ColoredCellsSorted2);
            int NoOfParticles = ParticleColorArray.Length;
            int[] GlobalParticleColor = new int[NoOfParticles];
            double[] StateBuffer = new double[NoOfParticles];
            for (int i = 0; i < NoOfParticles; i++)
            {
                StateBuffer[i] = Convert.ToDouble(ParticleColorArray[i]);
            }
            double[] GlobalStateBuffer = StateBuffer.MPIMax();
            for (int i = 0; i < NoOfParticles; i++)
            {
                GlobalParticleColor[i] = Convert.ToInt32(GlobalStateBuffer[i]);
            }
            int[,] ColorToRecolorWith = new int[GlobalParticleColor.Max() + 1, 2];
            for (int j = 0; j < J; j++)
            {
                GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                for (int i = 0; i < CellNeighbors.Length; i++)
                {
                    if (PartColEx[CellNeighbors[i]] != PartCol[j] && PartCol[j] != 0 && PartColEx[CellNeighbors[i]] > 0)
                    {
                        if (PartColEx[CellNeighbors[i]] < PartCol[j] || ColorToRecolorWith[PartCol[j], 1] > PartColEx[CellNeighbors[i]])
                        {
                            ColorToRecolorWith[PartCol[j], 0] = PartCol[j];
                            ColorToRecolorWith[PartCol[j], 1] = PartColEx[CellNeighbors[i]];
                        }
                        if (PartColEx[CellNeighbors[i]] > PartCol[j])
                        {
                            if (ColorToRecolorWith[PartColEx[CellNeighbors[i]], 0] == 0 || ColorToRecolorWith[PartColEx[CellNeighbors[i]], 1] > PartCol[j])
                            {
                                ColorToRecolorWith[PartColEx[CellNeighbors[i]], 0] = PartColEx[CellNeighbors[i]];
                                ColorToRecolorWith[PartColEx[CellNeighbors[i]], 1] = PartCol[j];
                            }
                        }
                    }
                }
            }
            int[][,] GlobalColorToRecolorWith = ColorToRecolorWith.MPIGatherO(0);
            GlobalColorToRecolorWith = GlobalColorToRecolorWith.MPIBroadcast(0);
            for (int m = 0; m < MPISize; m++)
            {
                for (int i = 0; i <= GlobalParticleColor.Max(); i++)
                {
                    if (GlobalColorToRecolorWith[0][i, 1] == 0 || GlobalColorToRecolorWith[0][i, 1] > GlobalColorToRecolorWith[m][i, 1] && GlobalColorToRecolorWith[m][i, 1] != 0)
                    {
                        GlobalColorToRecolorWith[0][i, 0] = GlobalColorToRecolorWith[m][i, 0];
                        GlobalColorToRecolorWith[0][i, 1] = GlobalColorToRecolorWith[m][i, 1];
                    }
                }
            }
            ColorToRecolorWith = GlobalColorToRecolorWith[0];
            for (int i = GlobalParticleColor.Max(); i >= 0; i--)
            {
                if (ColorToRecolorWith[i, 0] != 0)
                {
                    for (int j = 0; j < J; j++)
                    {
                        if (PartCol[j] == ColorToRecolorWith[i, 0])
                        {
                            ParticleColor.SetMeanValue(j, ColorToRecolorWith[i, 1]);
                            LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(rCode[j], 0));
                            PartCol[j] = ColorToRecolorWith[i, 1];
                        }
                    }
                }
            }
            return PartCol;
        }
        
        /// <summary>
        /// Initialization of <see cref="ParticleColor"/> 
        /// </summary>
        private int[] InitializeColoring(int J, IGridData GridData, bool AdaptiveMeshRefinement)
        {
            int[] Cells = new int[J];
            List<int> ColoredCells = new List<int>();
            for (int p = 0; p < m_Particles.Count; p++)
            {
                Particle Particle = m_Particles[p];
                double h_min = GridData.iGeomCells.h_min.Min();
                double h_max = GridData.iGeomCells.h_max.Max();
                for (int j = 0; j < J; j++)
                {
                    double[] center = GridData.iLogicalCells.GetCenter(j);
                    if (Particle.Contains(center, h_max, h_max))
                    {
                        ParticleColor.SetMeanValue(j, p + 1);
                        ColoredCells.Add(j);
                        Cells[j] = p + 1;
                    }
                }
            }
            CheckForNeighborColorsInit(Cells, GridData);
            return Cells;
        }

        private void CheckForNeighborColorsInit(int[] ColoredCells, IGridData GridData)
        {
            for (int i = 0; i < ColoredCells.Length; i++)
            {
                if (ColoredCells[i] != 0)
                {
                    GridData.GetCellNeighbours(i, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int j = 0; j < CellNeighbors.Length; j++)
                    {
                        if (CellNeighbors[j] < ColoredCells.Max() && ColoredCells[i] != ColoredCells[j] && ColoredCells[CellNeighbors[j]] != 0)
                        {
                            RecolorCellsInit(ColoredCells, ColoredCells[i], ColoredCells[CellNeighbors[j]], GridData);
                        }
                    }
                }
            }
        }

        private void RecolorCellsInit(int[] ColoredCells, int NewColor, int OldColor, IGridData GridData)
        {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int i = 0; i < J; i++)
            {
                if (ColoredCells[i] == OldColor)
                {
                    ColoredCells[i] = NewColor;
                    ParticleColor.SetMeanValue(i, NewColor);
                }
            }
        }

        private void CalculateHydrodynamicForces(List<Particle> Particles, double dt, bool firstIteration = true)
        {
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
            // ===============================================
            // Update forces
            // =============
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            for (int p = 0; p < Particles.Count(); p++)
            {
                Particle CurrentParticle = Particles[p];
                if (!((FSI_Control)Control).pureDryCollisions && !CurrentParticle.skipForceIntegration)
                {
                    CurrentParticle.UpdateForcesAndTorque(Velocity, Pressure, LsTrk, Control.PhysicalParameters.mu_A, dt, Control.PhysicalParameters.rho_A, firstIteration);
                    
                }
                else
                {
                    CurrentParticle.HydrodynamicForces[0][0] = 0;
                    CurrentParticle.HydrodynamicForces[0][1] = 0;
                    CurrentParticle.HydrodynamicTorque[0] = 0;
                }
            }
            // MPISum over Forces moved to Particle.cs 
        }

        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt)
        {
            using (new FuncTrace())
            {

                ResLogger.TimeStep = TimestepInt;

                hack_phystime = phystime;
                dt = GetFixedTimestep();

                Console.WriteLine("Starting time-step " + TimestepInt + "...");

                // used later to check if there is exactly one push per timestep
                // =============================================================
                int OldPushCount = LsTrk.PushCount;

                // =================================================
                // only particle motion & collisions, no flow solver
                // =================================================
                if (((FSI_Control)Control).pureDryCollisions)
                {
                    // in other branches, called by the BDF timestepper
                    LsTrk.PushStacks();
                    DGLevSet.Push();
                    ResetCollisionState(m_Particles);
                    CalculateHydrodynamicForces(m_Particles, dt);
                    Auxillary.CalculateParticleVelocity(m_Particles, dt, ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled, 0, TimestepInt, false);
                    CalculateCollision(m_Particles, GridData, LsTrk, CellColor, dt);
                    Auxillary.CalculateParticlePosition(m_Particles, dt);
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                    UpdateLevelSetParticles();
                    Auxillary.PrintResultToConsole(m_Particles, 0, 0, phystime, TimestepInt, 0, true, out double MPIangularVelocity, out Test_Force);
                    // Save for NUnit Test
                    base.QueryHandler.ValueQuery("C_Drag", 2 * Test_Force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                    base.QueryHandler.ValueQuery("C_Lift", 2 * Test_Force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
                    base.QueryHandler.ValueQuery("Angular_Velocity", MPIangularVelocity, true); // (TestCase FlowRotationalCoupling)
                    
                }
                // =============================================
                // particle motion & collisions plus flow solver
                // =============================================
                else
                {
                    // Collision triggered, no calculation of hydrodynamics
                    // ====================================================
                    if (triggerOnlyCollisionProcedure)
                    {
                        UpdateLevelSetParticles();
                        triggerOnlyCollisionProcedure = false;
                        return dt;
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
                    }
                    else if (((FSI_Control)this.Control).Timestepper_LevelSetHandling != LevelSetHandling.Coupled_Iterative)
                    {
                        IterationCounter = 0;
                        double posResidual_splitting = 1e12;
                        while (posResidual_splitting > ((FSI_Control)Control).ForceAndTorque_ConvergenceCriterion)
                        {
                            Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                            Auxillary.SaveOldParticleState(m_Particles, IterationCounter, ((FSI_Control)Control).ForceAndTorque_ConvergenceCriterion, FullyCoupled);
                            if (IterationCounter != 0 || ((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.FSI_LieSplittingFullyCoupled)
                            {
                                m_BDF_Timestepper.Solve(phystime, dt, false);
                                CalculateHydrodynamicForces(m_Particles, dt, false);
                            }
                            Auxillary.CalculateParticleVelocity(m_Particles, dt, FullyCoupled, IterationCounter, TimestepInt);
                            if (IterationCounter != 100000 || ((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.FSI_LieSplittingFullyCoupled)
                                Auxillary.PrintResultToConsole(m_Particles, ((FSI_Control)Control).PhysicalParameters.mu_A, ((FSI_Control)Control).PhysicalParameters.rho_A, phystime, TimestepInt, IterationCounter, false, out double _, out Test_Force);
                            if (((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.FSI_LieSplittingFullyCoupled)
                                break;
                            Auxillary.CalculateParticleResidual(m_Particles, IterationCounter, ((FSI_Control)Control).max_iterations_fully_coupled, out posResidual_splitting, out IterationCounter);
                        }
                        ResetCollisionState(m_Particles);
                        CalculateCollision(m_Particles, GridData, LsTrk, CellColor, dt);
                        Auxillary.CalculateParticlePosition(m_Particles, dt);
                        Auxillary.PrintResultToConsole(m_Particles, ((FSI_Control)Control).PhysicalParameters.mu_A, ((FSI_Control)Control).PhysicalParameters.rho_A, phystime, TimestepInt, IterationCounter, true, out double Test_RotationalVelocity, out Test_Force);
                        // Save for NUnit Test
                        base.QueryHandler.ValueQuery("C_Drag", 2 * Test_Force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                        base.QueryHandler.ValueQuery("C_Lift", 2 * Test_Force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
                        base.QueryHandler.ValueQuery("Angular_Velocity", Test_RotationalVelocity, true); // (TestCase FlowRotationalCoupling)
                        if (((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled)
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
                    else
                    {
                        foreach (Particle p in m_Particles)
                        {
                            p.iteration_counter_P = -1;
                            p.forceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
                        }
                        m_BDF_Timestepper.Solve(phystime, dt, false);
                    }
                }

                // finalize
                // ========
                if (LsTrk.PushCount - OldPushCount != 1)
                {
                    // To whom it concerns / who stumbles across this exception:
                    // It is important that LevelSetTracker.PushStacks() is called *exactly once per time-step*, at the beginning.
                    // Do not remove this check! Instead, remove any calls to 'PushStacks()' in subroutines of this method.
                    // Fk.
                    throw new ApplicationException("Illegal number of level-set push actions in time-step." + (LsTrk.PushCount - OldPushCount));
                }

                ResLogger.NextTimestep(false);

                Console.WriteLine("done with time-step.");
                return dt;
            }
        }

        internal void ResetCollisionState(List<Particle> Particles)
        {
            foreach (Particle p in Particles)
            {
                p.skipForceIntegration = false;
                p.Collided = false;
            }
        }

        // restart
        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t)
        {
            var tsi = new FSI_TimestepInfo(t, this.CurrentSessionInfo, timestepno, base.IOFields, m_Particles);

            SerialzeTester(tsi);

            return tsi;
        }

        /// <summary>
        /// Test the serialization of <see cref="FSI_TimestepInfo.Particles"/>
        /// </summary>
        private static void SerialzeTester(FSI_TimestepInfo b)
        {
            JsonSerializer formatter = new JsonSerializer()
            {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Serialize
            };

            bool DebugSerialization = false;

            JsonReader GetJsonReader(Stream s)
            {
                if (DebugSerialization)
                {
                    return new JsonTextReader(new StreamReader(s));
                }
                else
                {
                    return new BsonReader(s);
                }
            }

            JsonWriter GetJsonWriter(Stream s)
            {
                if (DebugSerialization)
                {
                    return new JsonTextWriter(new StreamWriter(s));
                }
                else
                {
                    return new BsonWriter(s);
                }
            }


            byte[] buffer = null;
            using (var ms1 = new MemoryStream())
            {
                using (var writer = GetJsonWriter(ms1))
                {
                    formatter.Serialize(writer, b);
                    writer.Flush();
                    buffer = ms1.GetBuffer();
                    //writer.Close();
                }
            }

            FSI_TimestepInfo o;
            using (var ms2 = new MemoryStream(buffer))
            {
                using (var reader = GetJsonReader(ms2))
                {
                    o = formatter.Deserialize<FSI_TimestepInfo>(reader);
                    reader.Close();
                }
            }

            //Console.WriteLine(o.ToString());

            Debug.Assert(b.Particles.Length == o.Particles.Length);
            int L = b.Particles.Length;
            for (int l = 0; l < L; l++)
            { // loop over particles
                Debug.Assert(GenericBlas.L2Dist(b.Particles[l].Position[0], o.Particles[l].Position[0]) < 1e-13);
            }

        }

        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepNumber RestartFromDatabase(out double time)
        {

            // this sux, because the database API is totally fucked up
            var db = GetDatabase();
            Guid Rst_Tsid = base.GetRestartTimestepID();
            Guid Rst_SessionId = Control.RestartInfo.Item1;
            ISessionInfo session = db.Controller.GetSessionInfo(Rst_SessionId);

            var ArschInfo = ((DatabaseDriver)(base.DatabaseDriver)).LoadTimestepInfo<FSI_TimestepInfo>(Rst_Tsid, session, db);

            // init particles
            m_Particles = ArschInfo.Particles.ToList();
            hack_phystime = ArschInfo.PhysicalTime;
            UpdateLevelSetParticles();

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
        public override void PostRestart(double time, TimestepNumber timestep)
        {
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

        
            //PlotCurrentState(0.2, new TimestepNumber(2), 3);



        //bool collision = false;
            //PlotCurrentState(0.4, new TimestepNumber(3), 3);


            //PlotCurrentState(0.4, new TimestepNumber(4), 3);


        bool triggerOnlyCollisionProcedure = false;

        private void CalculateCollision(List<Particle> Particles, IGridData GridData, LevelSetTracker LsTrk, int[] CellColor, double dt)
        {
            if (CollisionModel == FSI_Control.CollisionModel.NoCollisionModel)
                return;

            if (CollisionModel == FSI_Control.CollisionModel.RepulsiveForce)
                throw new NotImplementedException("Repulsive force model is currently unsupported, please use the momentum conservation model.");
            FSI_Collision _Collision = new FSI_Collision(
                ((FSI_Control)Control).PhysicalParameters.mu_A,
                ((FSI_Control)Control).PhysicalParameters.rho_A, 
                ((FSI_Control)Control).CoefficientOfRestitution, 
                dt,
                ((GridData)GridData).Cells.h_minGlobal
                );
            _Collision.CalculateCollision(Particles, GridData, LsTrk, CellColor);
            foreach (Particle p in m_Particles)
            {
                _Collision.Collision_MPICommunication(m_Particles, p, MPISize);
            }
        }

        /// <summary>
        /// Mesh refinement
        /// Very primitive refinement indicator, works on a LevelSet criterion.
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel)
        {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate();
            CellMask LevSetCells = LsTrk.Regions.GetCutCellMask();
            int DesiredLevel_j = 0;
            if (LevSetCells.Contains(j))//ColoredCellMask != null &&
            {
                DesiredLevel_j = ((FSI_Control)this.Control).RefinementLevel;
            }
            else if (LevSetCells.AllNeighbourCells().Contains(j))
            {
                DesiredLevel_j = ((FSI_Control)this.Control).RefinementLevel - 1;
            }
            else
                DesiredLevel_j = 0;

            return DesiredLevel_j;
        }


        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid)
        {

            if (((FSI_Control)Control).AdaptiveMeshRefinement && IterationCounter == 0)
            {
                if (TimestepNo > 3 && TimestepNo % 3 != 0)
                {
                    newGrid = null;
                    old2NewGrid = null;
                    return;
                }

                // Check grid changes
                // ==================

                // compute curvature for levelindicator 
                //CurvatureAlgorithms.CurvatureDriver(
                //SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //CurvatureAlgorithms.FilterConfiguration.Default,
                //this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                //this.HMForder, this.DGLevSet.Current);

                CellMask CutCells = LsTrk.Regions.GetCutCellMask();
                CellMask CutCellNeighbors = LsTrk.Regions.GetNearFieldMask(1);
                CutCells = CutCells.Union(CutCellNeighbors);

                bool AnyChange = GridRefinementController.ComputeGridChange((GridData)GridData, CutCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if (AnyChange)
                {
                    int[] glb = (new int[] {

                    CellsToRefineList.Count,
                    Coarsening.Sum(L => L.Length),
                }).MPISum();
                    NoOfCellsToRefine = glb[0];
                    NoOfCellsToCoarsen = glb[1];
                }
                int oldJ = this.GridData.CellPartitioning.TotalLength;

                if (AnyChange)
                {
                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");
                    newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);
                }
                else
                {
                    newGrid = null;
                    old2NewGrid = null;
                }
            }
            else
            {
                newGrid = null;
                old2NewGrid = null;
            }
        }
    }
}