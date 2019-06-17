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
using BoSSS.Solution.LevelSetTools;
using NUnit.Framework;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Quadrature;

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

            foreach (Particle p in m_Particles)
            {
                p.m_collidedWithParticle = new bool[m_Particles.Count];
                p.m_collidedWithWall = new bool[4];
                p.m_closeInterfacePointTo = new double[m_Particles.Count][];
                p.ClosestPointToParticle = new double[m_Particles.Count, 2];
            }
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

        double calculatedDampingTensors;
        readonly private FSI_Auxillary Auxillary = new FSI_Auxillary();
        readonly private FSI_Collision Collision = new FSI_Collision();
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

        private int iteration_counter = 0;

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

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        {

            if (IBM_Op != null)
                return;
            string[] CodNameSelected = new string[0];
            string[] DomNameSelected = new string[0];

            int D = this.GridData.SpatialDimension;

            BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

            int degU = this.Velocity[0].Basis.Degree;
            var IBM_Op_config = new NSEOperatorConfiguration
            {
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
                                            else { containsParticle = p.Contains(X, LsTrk); }

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

                                    double[] result = new double[X.Length + 7];

                                    foreach (Particle p in m_Particles)
                                    {
                                        // which particle?
                                        bool containsParticle;
                                        if (m_Particles.Count == 1)
                                        {
                                            containsParticle = true;
                                        }
                                        else { containsParticle = p.Contains(X, LsTrk); }

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
                                   else { containsParticle = p.Contains(X, LsTrk); }
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
                    if (calculatedDampingTensors == 0)
                    {
                        foreach (Particle p in m_Particles)
                        {
                            if (p.neglectAddedDamping == false)
                            {
                                p.CalculateDampingTensor(LsTrk, ((FSI_Control)this.Control).PhysicalParameters.mu_A, ((FSI_Control)this.Control).PhysicalParameters.rho_A, ((FSI_Control)this.Control).dtMax);
                                Auxillary.ExchangeDampingTensors(m_Particles);
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
                                Auxillary.ExchangeDampingTensors(m_Particles);
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
                    UpdateForcesAndTorque(m_Particles, dt, 0);
                    UpdateLevelSetParticles();
                    foreach (Particle p in m_Particles)
                    {
                        p.iteration_counter_P += 1;
                        p.ForceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
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
            CellColor = ((FSI_Control)Control).AdaptiveMeshRefinement ? InitializeColoring(J, ((FSI_Control)Control).AdaptiveMeshRefinement) : CellColor ?? InitializeColoring(J, ((FSI_Control)Control).AdaptiveMeshRefinement);

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
                int CurrentParticle = p;
                bool ContainsCurrentColor = false;
                BitArray ColoredCells = new BitArray(J);
                for (int j = 0; j < J; j++)
                {
                    if (CellColor[j] == CurrentColor && CurrentColor != 0)
                    {
                        ContainsCurrentColor = true;
                        ColoredCells[j] = true;
                    }
                }
                if (ContainsCurrentColor)
                {
                    int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(GlobalParticleColor, CurrentColor);
                    //CellMask ColoredCellMask = levelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, CurrentColor, J, false);
                    CellMask ColoredCellMask = new CellMask(GridData, ColoredCells);
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
                            phi *= m_Particles[ParticlesOfCurrentColor[pc]].Phi_P(X);

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
            PerformLevelSetSmoothing(AgglParticleMask);
                        
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
            // =======================================================
            // Step 1
            // Color all cells directlly related to the level set
            // =======================================================
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int Je = J + GridData.iLogicalCells.NoOfExternalCells;
            int[] PartCol = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int ColorOffset = MPIRank * m_Particles.Count();
            List<int> ColoredCellsSorted = new List<int>();
            int ListIndex = 0;
            for (int CellID = 0; CellID < PartCol.Length; CellID++)
            {
                if (PartCol[CellID] != 0)
                {
                    ListIndex = 0;
                    if (ColoredCellsSorted.Count != 0)
                    {
                        while (ListIndex < ColoredCellsSorted.Count && PartCol[CellID] >= ColoredCellsSorted[ListIndex])
                        {
                            ListIndex += 1;
                        }
                    }
                    int temp = PartCol[CellID];
                    if (ListIndex == 0 || ColoredCellsSorted[ListIndex - 1] != temp)
                        ColoredCellsSorted.Insert(ListIndex, temp);
                }
            }
            int[] LocalColors = ColoredCellsSorted.ToArray();
            int ColorMax = PartCol.Max();
            List<int[]> LocalColorChunk = new List<int[]>();
            bool[] CheckArray = new bool[PartCol.Length];
            int index = 1;
            for (int j = 0; j < J; j++)
            {
                if (PartCol[j] != 0 && !CheckArray[j])
                {
                    int CurrentColor = PartCol[j];
                    CheckArray[j] = true;
                    LocalColorChunk.Add(new int[] { index, j });
                    for (int i = 0; i < LocalColorChunk.Count(); i++)
                    {
                        GridData.GetCellNeighbours(LocalColorChunk[i][1], GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                        for (int k = 0; k < CellNeighbors.Length; k++)
                        {
                            int NeighbourColor = PartCol[CellNeighbors[k]];
                            if (NeighbourColor == CurrentColor && !CheckArray[CellNeighbors[k]] && CellNeighbors[k] < J)
                            {
                                LocalColorChunk.Add(new int[] { index, CellNeighbors[k] });
                                CheckArray[CellNeighbors[k]] = true;
                            }
                        }
                    }
                    index += 1;
                }
            }
            for (int i = 0; i < LocalColorChunk.Count(); i++)
            {
                PartCol[LocalColorChunk[i][1]] = LocalColorChunk[i][0] + ColorOffset;
            }

            int[] PartColEx = new int[Je];
            var rCode = LsTrk.Regions.RegionsCode;
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
        private int[] InitializeColoring(int J, bool AdaptiveMeshRefinement)
        {
            int[] Cells = new int[J];
            for (int p = 0; p < m_Particles.Count; p++)
            {
                double Hmin = AdaptiveMeshRefinement ? Math.Sqrt(GridData.iGeomCells.GetCellVolume(0)) * 2 : Math.Sqrt(GridData.iGeomCells.GetCellVolume(0));
                double[] ParticlePos = m_Particles[p].Position[0];
                double ParticleAngle = m_Particles[p].Angle[0];
                double[] ParticleScales = m_Particles[p].GetLengthScales();
                double Upperedge = ParticlePos[1] + ParticleScales[1] * Math.Abs(Math.Cos(ParticleAngle)) + ParticleScales[0] * Math.Abs(Math.Sin(ParticleAngle)) + 1 * Hmin;
                double Loweredge = ParticlePos[1] - ParticleScales[1] * Math.Abs(Math.Cos(ParticleAngle)) - ParticleScales[0] * Math.Abs(Math.Sin(ParticleAngle)) - 1 * Hmin;
                double Leftedge = ParticlePos[0] - ParticleScales[0] * Math.Abs(Math.Cos(ParticleAngle)) - ParticleScales[1] * Math.Abs(Math.Sin(ParticleAngle)) - 1 * Hmin;
                double Rightedge = ParticlePos[0] + ParticleScales[0] * Math.Abs(Math.Cos(ParticleAngle)) + ParticleScales[1] * Math.Abs(Math.Sin(ParticleAngle)) + 1 * Hmin;
                for (int j = 0; j < J; j++)
                {
                    double[] center = GridData.iLogicalCells.GetCenter(j);
                    if (center[0] > Leftedge && center[0] < Rightedge && center[1] > Loweredge && center[1] < Upperedge)
                    {
                        ParticleColor.SetMeanValue(j, p + 1);
                        m_Particles[p].ParticleColoredCells.Add(new int[2] { j, p + 1 });
                        Cells[j] = p + 1;
                    }
                }
            }
            CheckForNeighborColorsInit(Cells);
            return Cells;
        }

        private void CheckForNeighborColorsInit(int[] ColoredCells)
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
                            RecolorCellsInit(ColoredCells, ColoredCells[i], ColoredCells[CellNeighbors[j]]);
                        }
                    }
                }
            }
        }

        private void RecolorCellsInit(int[] ColoredCells, int NewColor, int OldColor)
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

        private void UpdateForcesAndTorque(List<Particle> Particles, double dt, int iteration_counter)
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
            for (int p = 0; p < m_Particles.Count(); p++)
            {
                Particle CurrentParticle = m_Particles[p];
                if (!((FSI_Control)Control).pureDryCollisions && !CurrentParticle.skipForceIntegration)
                    CurrentParticle.UpdateForcesAndTorque(Velocity, Pressure, LsTrk, Control.PhysicalParameters.mu_A, dt, Control.PhysicalParameters.rho_A, ((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.FSI_LieSplittingFullyCoupled);
                else
                {
                    CurrentParticle.HydrodynamicForces[0][0] = 0;
                    CurrentParticle.HydrodynamicForces[0][1] = 0;
                    CurrentParticle.HydrodynamicTorque[0] = 0;
                }
                // wall collisions are computed on each processor
                //WallCollisionForces(CurrentParticle, p, LsTrk.GridDat.Cells.h_minGlobal);
                //Auxillary.Collision_MPICommunication(m_Particles, CurrentParticle, MPISize, true);
            }
            //if (m_Particles.Count > 1)
            //{
            //    UpdateCollisionForces(Particles, LsTrk.GridDat.Cells.h_minGlobal, dt, iteration_counter);
            //    foreach(Particle p in m_Particles)
            //    {
            //        Auxillary.Collision_MPICommunication(m_Particles, p, MPISize);
            //    }
                
            //}
            // MPISum over Forces moved to Particle.cs 
        }

        private void CalculateCollision(List<Particle> Particles, double Hmin, double dt, int iteration_counter)
        {
            UpdateCollisionForces(Particles, LsTrk.GridDat.Cells.h_minGlobal, dt, iteration_counter);
            foreach (Particle p in m_Particles)
            {
                Collision.Collision_MPICommunication(m_Particles, p, MPISize);
            }
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
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);

                    UpdateForcesAndTorque(m_Particles, dt, 0);
                    
                    foreach (Particle p in m_Particles)
                    {
                        p.CalculateAcceleration(dt, ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled, false);
                        p.UpdateParticleVelocity(dt);
                    }
                    foreach (Particle p in m_Particles)
                    {
                        if (p.skipForceIntegration)
                            p.skipForceIntegration = false;
                    }
                    if (m_Particles.Count() > 1)
                        CalculateCollision(m_Particles, LsTrk.GridDat.Cells.h_minGlobal, dt, iteration_counter);
                    for (int p = 0; p < m_Particles.Count(); p++)
                    {
                        Particle CurrentParticle = m_Particles[p];
                        WallCollisionForcesNew(CurrentParticle, dt);
                        Collision.Collision_MPICommunication(m_Particles, CurrentParticle, MPISize, true);
                    }
                    foreach (Particle p in m_Particles)
                    {
                        p.UpdateParticlePositionAndAngle(dt);
                    }
                    UpdateLevelSetParticles();
                    Auxillary.PrintResultToConsole(m_Particles, phystime, dt, 0, true, out double MPIangularVelocity, out force);
                    // Save for NUnit Test
                    base.QueryHandler.ValueQuery("C_Drag", 2 * force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                    base.QueryHandler.ValueQuery("C_Lift", 2 * force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
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
                        iteration_counter = 0;
                        double posResidual_splitting = 1e12;
                        while (posResidual_splitting > ((FSI_Control)Control).ForceAndTorque_ConvergenceCriterion)
                        {
                            double[] ForcesOldSquared = new double[2];
                            double TorqueOldSquared = new double();

                            if (iteration_counter != 0 || ((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.FSI_LieSplittingFullyCoupled)
                            {
                                Auxillary.SaveOldParticleState(m_Particles, iteration_counter, 2, ((FSI_Control)Control).ForceAndTorque_ConvergenceCriterion, ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled, out ForcesOldSquared, out TorqueOldSquared);
                                m_BDF_Timestepper.Solve(phystime, dt, false);
                                Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                                UpdateForcesAndTorque(m_Particles, dt, iteration_counter);
                            }

                            foreach (Particle p in m_Particles)
                            {
                                p.iteration_counter_P = iteration_counter;
                                Auxillary.UpdateParticleAccelerationAndDamping(p, iteration_counter, dt, ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled);
                                p.UpdateParticleVelocity(dt);
                            }

                            Auxillary.PrintResultToConsole(m_Particles, phystime, dt, iteration_counter, false, out double _, out force);

                            if (((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.FSI_LieSplittingFullyCoupled)
                                break;
                            Auxillary.CalculateParticleResidual(m_Particles, ForcesOldSquared, TorqueOldSquared, iteration_counter, ((FSI_Control)Control).max_iterations_fully_coupled, out posResidual_splitting, out iteration_counter);
                        }
                        foreach (Particle p in m_Particles)
                        {
                            if (p.skipForceIntegration)
                                p.skipForceIntegration = false;
                        }
                        if (m_Particles.Count() > 1)
                            CalculateCollision(m_Particles, LsTrk.GridDat.Cells.h_minGlobal, dt, iteration_counter);
                        for (int p = 0; p < m_Particles.Count(); p++)
                        {
                            Particle CurrentParticle = m_Particles[p];
                            WallCollisionForcesNew(CurrentParticle, p);
                            Collision.Collision_MPICommunication(m_Particles, CurrentParticle, MPISize, true);
                        }
                        foreach (Particle p in m_Particles)
                        {
                            p.UpdateParticlePositionAndAngle(dt);
                        }
                        Auxillary.PrintResultToConsole(m_Particles, phystime, dt, iteration_counter, true, out double MPIangularVelocity, out force);
                        // Save for NUnit Test
                        base.QueryHandler.ValueQuery("C_Drag", 2 * force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                        base.QueryHandler.ValueQuery("C_Lift", 2 * force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
                        base.QueryHandler.ValueQuery("Angular_Velocity", MPIangularVelocity, true); // (TestCase FlowRotationalCoupling)
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
                            p.ForceAndTorque_convergence = ((FSI_Control)this.Control).ForceAndTorque_ConvergenceCriterion;
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



            foreach (Particle p in m_Particles)
            {
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

        // Collision models
        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="particle0"></param>
        /// <param name="particle1"></param>
        public void UpdateCollisionForces(List<Particle> Particles, double hmin, double dt, int iteration_counter)
        {
            if (CollisionModel == FSI_Control.CollisionModel.NoCollisionModel)
                 return;

            if (Particles.Count < 2)
                return;

            LevelSetUpdate.DetermineGlobalParticleColor(GridData, CellColor, Particles, out int[] GlobalParticleColor);

            for (int i = 0; i < GlobalParticleColor.Length; i++)
            {
                int CurrentColor = GlobalParticleColor[i];
                int[] ParticlesOfCurrentColor = LevelSetUpdate.FindParticlesOneColor(GlobalParticleColor, CurrentColor);

                if (ParticlesOfCurrentColor.Length > 1 && CurrentColor != 0)
                {
                    int SpatialDim = 2;
                    MultidimensionalArray SaveTimeStepArray = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length);
                    MultidimensionalArray Distance = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length);
                    MultidimensionalArray DistanceVector = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length, SpatialDim);
                    MultidimensionalArray PositionDistanceVector = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length, SpatialDim);
                    MultidimensionalArray ClosestPoint_P0 = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length, SpatialDim);
                    MultidimensionalArray ClosestPoint_P1 = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length, SpatialDim);
                    double MaxDistance = 1e-4;
                    double AccDynamicTimestep = 0;
                    bool Overlapping = false;

                    while (AccDynamicTimestep < dt)
                    {
                        double MinDistance = double.MaxValue;
                        double SaveTimeStep = 0;
                        while (MinDistance > MaxDistance)
                        {
                            for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                            {
                                Collision.UpdateParticleState(Particles[ParticlesOfCurrentColor[p]], SaveTimeStep, 2);
                            }
                            Overlapping = false;
                            SaveTimeStep = double.MaxValue;
                            for (int p0 = 0; p0 < ParticlesOfCurrentColor.Length; p0++)
                            {
                                for (int p1 = p0 + 1; p1 < ParticlesOfCurrentColor.Length; p1++)
                                {
                                    Particle Particle0 = Particles[ParticlesOfCurrentColor[p0]];
                                    Particle Particle1 = Particles[ParticlesOfCurrentColor[p1]];
                                    Collision.ComputeMinimalDistance(Particle0, Particle1, LsTrk, out double temp_Distance, out MultidimensionalArray temp_DistanceVector, out MultidimensionalArray temp_ClosestPoint_p0, out MultidimensionalArray temp_ClosestPoint_p1, out bool temp_Overlapping);
                                    Collision.CalculateNormalAndTangentialVector(temp_DistanceVector.To1DArray(), out double[] NormalVector, out double[] TangentialVector);
                                    double temp_SaveTimeStep = Collision.DynamicTimestep(Particle0, Particle1, temp_ClosestPoint_p0.To1DArray(), temp_ClosestPoint_p1.To1DArray(), NormalVector, temp_Distance);
                                    Distance[p0, p1] = temp_Distance;
                                    SaveTimeStepArray[p0, p1] = temp_SaveTimeStep;
                                    DistanceVector.SetSubArray(temp_DistanceVector, new int[] { p0, p1, -1 });
                                    ClosestPoint_P0.SetSubArray(temp_ClosestPoint_p0, new int[] { p0, p1, -1 });
                                    ClosestPoint_P1.SetSubArray(temp_ClosestPoint_p1, new int[] { p0, p1, -1 });
                                    if (temp_SaveTimeStep < SaveTimeStep && temp_SaveTimeStep > 0)
                                    {
                                        SaveTimeStep = temp_SaveTimeStep;
                                        MinDistance = temp_Distance;
                                    } 
                                    if (temp_Overlapping)
                                    {
                                        SaveTimeStep = -dt;
                                        MinDistance = double.MaxValue;
                                        Overlapping = true;
                                        Console.WriteLine("Computing collision, Particle0Pos: " + Particle0.Position[0][0] + ", " + Particle0.Position[0][1] + "Particle0Pos1: " + Particle0.Position[1][0] + ", " + Particle0.Position[1][1]  + "Overlapping ? " + Overlapping);
                                    }
                                }
                            }
                            if (SaveTimeStep >= 0)
                                AccDynamicTimestep += SaveTimeStep;
                            if (AccDynamicTimestep > dt)
                                break;
                            Console.WriteLine("Computing collision, acc dynamic dt: " + AccDynamicTimestep + "Overlapping? " + Overlapping);
                        }

                        if (AccDynamicTimestep > dt)
                            break;
                        if (MinDistance > MaxDistance && !Overlapping)
                            break;
                        for (int p0 = 0; p0 < ParticlesOfCurrentColor.Length; p0++)
                        {
                            for (int p1 = p0 + 1; p1 < ParticlesOfCurrentColor.Length; p1++)
                            {
                                Particle Particle0 = Particles[ParticlesOfCurrentColor[p0]];
                                Particle Particle1 = Particles[ParticlesOfCurrentColor[p1]];
                                if (Distance[p0, p1] < MaxDistance && SaveTimeStepArray[p0, p1] > 0)
                                {
                                    double[] CurrentDistanceVector = DistanceVector.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                                    double[] CurrentClosestPoint_P0 = ClosestPoint_P0.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                                    double[] CurrentClosestPoint_P1 = ClosestPoint_P1.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                                    Console.WriteLine("Computing collision,MinDistance: " + MinDistance + "ParticlesOfCurrentColor.Length " + ParticlesOfCurrentColor.Length);
                                    Console.WriteLine("Computing collision, acc dynamic dt: " + AccDynamicTimestep + ", Distance: " + Distance[p0, p1] + ", SaveTimeStep " + SaveTimeStepArray[p0,p1]);
                                    Console.WriteLine("Particle0 Velocity: " + Particle0.TranslationalVelocity[0][0] + ", Particle0.TranslationalVelocity[0][1]: " + Particle0.TranslationalVelocity[0][1] + ", Particle0.TranslationalVelocity[1][0]: " + Particle0.TranslationalVelocity[1][0] + ", Particle0.TranslationalVelocity[1][1]: " + Particle0.TranslationalVelocity[1][1]);
                                    Collision.ComputeMomentumBalanceCollision(m_Particles, Particle0, Particle1, CurrentDistanceVector, CurrentClosestPoint_P0, CurrentClosestPoint_P1, ((FSI_Control)Control).CoefficientOfRestitution);
                                    Particle0.CollisionTimestep = AccDynamicTimestep;
                                    Particle1.CollisionTimestep = AccDynamicTimestep;
                                }
                                else
                                {
                                    Particle0.m_collidedWithParticle[m_Particles.IndexOf(Particle1)] = false;
                                    Particle1.m_collidedWithParticle[m_Particles.IndexOf(Particle0)] = false;
                                    Particle0.m_closeInterfacePointTo[m_Particles.IndexOf(Particle1)] = null;
                                    Particle1.m_closeInterfacePointTo[m_Particles.IndexOf(Particle0)] = null;
                                    triggerOnlyCollisionProcedure = false;
                                }
                            }
                        }
                        for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                        {
                            Collision.SumOverCollisionVelocities(Particles[ParticlesOfCurrentColor[p]]);
                        }
                    }
                }
                for (int j = 0; j < GlobalParticleColor.Length; j++)
                {
                    if (GlobalParticleColor[j] == CurrentColor)
                        GlobalParticleColor[j] = 0;
                }
            }
        }

        /// <summary>
        /// Update of particle state (velocity, force, etc.) for two particles where a collision is detected
        /// </summary>
        private void ComputeCollisionModel(double h_min, Particle Particle0, Particle Particle1, ref double Distance, ref double[] DistanceVector, double dt, int iteration_counter, out bool Collided)
        {
            // =======================================================
            // Step 0
            // Some Instantiations.
            // =======================================================
            FSI_Collision _FSI_Collision = new FSI_Collision();
            FSI_Auxillary _FSI_Auxillary = new FSI_Auxillary();
            int SpatialDim = DistanceVector.Length;
            bool ForceCollision = false;
            double[] ClosestPoint_P0 = new double[SpatialDim];
            double[] ClosestPoint_P1 = new double[SpatialDim];
            bool Overlapping = false;
            double collisionVn_P0 = 0;
            double collisionVt_P0 = 0;
            double collisionVn_P1 = 0;
            double collisionVt_P1 = 0;
            double[] NormalVector = new double[2];
            double[] TangentialVector = new double[2];
            // =======================================================
            // Step 1
            // Calculate the minimum distance between two particles.
            // =======================================================
            for (int i = 0; i < Particle0.NoOfSubParticles(); i++)
            {
                for (int j = 0; j < Particle1.NoOfSubParticles(); j++)
                {
                    Collision.GJK_DistanceAlgorithm(Particle0, i, Particle1, j, LsTrk, Particle0.Position[0], Particle1.Position[0], Particle0.Angle[0], Particle1.Angle[0], out double temp_Distance, out double[] temp_DistanceVector, out double[] temp_ClosestPoint_P0, out double[] temp_ClosestPoint_P1, out Overlapping);
                    if (Overlapping)
                        break;
                    if (temp_Distance < Distance)
                    {
                        Distance = temp_Distance;
                        DistanceVector = temp_DistanceVector;
                        ClosestPoint_P0 = temp_ClosestPoint_P0;
                        ClosestPoint_P1 = temp_ClosestPoint_P1;
                    }
                }
            }

            // Save closest points to particle.cs
            for (int d = 0; d < 2; d++)
            {
                Particle0.ClosestPointToParticle[m_Particles.IndexOf(Particle1), d] = ClosestPoint_P0[d];
                Particle1.ClosestPointToParticle[m_Particles.IndexOf(Particle0), d] = ClosestPoint_P1[d];
            }
            // =======================================================
            // Step 2
            // Project velocity on normal/tangential vector.
            // =======================================================
            _FSI_Collision.CalculateNormalAndTangentialVector(DistanceVector, out NormalVector, out TangentialVector);
            _FSI_Collision.ProjectVelocity(NormalVector, TangentialVector, Particle0.TranslationalVelocity[0], out collisionVn_P0, out collisionVt_P0);
            _FSI_Collision.ProjectVelocity(NormalVector, TangentialVector, Particle1.TranslationalVelocity[0], out collisionVn_P1, out collisionVt_P1);

            // =======================================================
            // Step 3
            // Calculate dynamic threshold.
            // =======================================================
            _FSI_Collision.CalculateDynamicCollisionThreshold(Particle0, Particle1, ClosestPoint_P0, ClosestPoint_P1, NormalVector, Distance, dt, out double Threshold);
            //Dynamic_dt = _FSI_Collision.DynamicTimestep(Particle0, Particle1, ClosestPoint_P0, ClosestPoint_P1, NormalVector, Distance, dt);

            if (double.IsNaN(collisionVn_P0) || double.IsInfinity(collisionVn_P0))
                throw new ArithmeticException("Error trying to update particle position. Value:  " + collisionVn_P0);
            // =======================================================
            // Step 4
            // Check whether the particles would collide 
            // with the velocities of the next timestep.
            // =======================================================
            //if (!Overlapping && Threshold == 0)
            //{
            //    double[] PointVelocity0 = new double[2];
            //    double[] PointVelocity1 = new double[2];
            //    double DetectCollisionVn_P0;
            //    double DetectCollisionVn_P1;
            //    bool Overlapping_NextTimestep = false;

            //    // Predict particle state of next timestep.
            //    _FSI_Collision.PredictParticleNextTimestep(Particle0, SpatialDim, dt, out double[] VirtualPosition0, out double[] VirtualVelocity0, out double VirtualAngle0, out double VirtualRotationalVelocity0);
            //    _FSI_Collision.PredictParticleNextTimestep(Particle1, SpatialDim, dt, out double[] VirtualPosition1, out double[] VirtualVelocity1, out double VirtualAngle1, out double VirtualRotationalVelocity1);

            //    // Calculate the minimum distance between the two particles.
            //    for (int i = 0; i < Particle0.NoOfSubParticles(); i++)
            //    {
            //        for (int j = 0; j < Particle1.NoOfSubParticles(); j++)
            //        {
            //            Collision.GJK_DistanceAlgorithm(Particle0, i, Particle1, j, LsTrk, Particle0.Position[0], Particle1.Position[0], Particle0.Angle[0], Particle1.Angle[0], out double temp_Distance, out double[] temp_DistanceVector, out double[] temp_ClosestPoint_P0, out double[] temp_ClosestPoint_P1, out Overlapping_NextTimestep);
            //            if (Overlapping_NextTimestep)
            //                break;
            //            if (temp_Distance < Distance)
            //            {
            //                Distance = temp_Distance;
            //                DistanceVector = temp_DistanceVector;
            //                ClosestPoint_P0 = temp_ClosestPoint_P0;
            //                ClosestPoint_P1 = temp_ClosestPoint_P1;
            //            }
            //        }
            //    }

            //    // Calculate dynamic threshold.
            //    _FSI_Collision.CalculateNormalAndTangentialVector(DistanceVector, out NormalVector, out TangentialVector);
            //    _FSI_Collision.CalculateRadialVector(VirtualPosition0, ClosestPoint_P0, out _, out double RadialLength0, out double[] RadialNormalVector0);
            //    _FSI_Collision.CalculateRadialVector(VirtualPosition1, ClosestPoint_P1, out _, out double RadialLength1, out double[] RadialNormalVector1);
            //    _FSI_Collision.TransformRotationalVelocity(VirtualRotationalVelocity0, RadialLength0, RadialNormalVector0, out double[] PointVelocityDueToRotation0);
            //    _FSI_Collision.TransformRotationalVelocity(VirtualRotationalVelocity1, RadialLength1, RadialNormalVector1, out double[] PointVelocityDueToRotation1);
            //    for (int d = 0; d < 2; d++)
            //    {
            //        PointVelocity0[d] = VirtualVelocity0[d] + PointVelocityDueToRotation0[d];
            //        PointVelocity1[d] = VirtualVelocity1[d] + PointVelocityDueToRotation1[d];
            //    }
            //    _FSI_Collision.ProjectVelocityOnVector(NormalVector, VirtualVelocity0, out DetectCollisionVn_P0);
            //    _FSI_Collision.ProjectVelocityOnVector(NormalVector, VirtualVelocity1, out DetectCollisionVn_P1);
            //    _FSI_Collision.ProjectVelocity(NormalVector, TangentialVector, VirtualVelocity0, out collisionVn_P0, out collisionVt_P0);
            //    _FSI_Collision.ProjectVelocity(NormalVector, TangentialVector, VirtualVelocity1, out collisionVn_P1, out collisionVt_P1);
            //    if (Distance <= Math.Abs((-DetectCollisionVn_P0 + DetectCollisionVn_P1) * dt))
            //        Threshold = Math.Abs((-DetectCollisionVn_P0 + DetectCollisionVn_P1) * dt);
            //    if (Overlapping_NextTimestep)
            //        Threshold = double.MaxValue;
            //}
            //Particle0.m_closeInterfacePointTo[m_Particles.IndexOf(Particle1)] = ClosestPoint_P0;
            //Particle1.m_closeInterfacePointTo[m_Particles.IndexOf(Particle0)] = ClosestPoint_P1;

            // =======================================================
            // Step 5
            // Emergency procedure if particles are overlapping.
            // =======================================================
            if (Overlapping)
            {
                // Set particle position to the last position
                Particle0.Angle[0] = Particle0.Angle[1];
                Particle1.Angle[0] = Particle1.Angle[1];
                Particle0.Position[0] = Particle0.Position[1].CloneAs();
                Particle1.Position[0] = Particle1.Position[1].CloneAs();

                // Calculate the minimum distance between the two particles.
                for (int i = 0; i < Particle0.NoOfSubParticles(); i++)
                {
                    for (int j = 0; j < Particle1.NoOfSubParticles(); j++)
                    {
                        Collision.GJK_DistanceAlgorithm(Particle0, i, Particle1, j, LsTrk, Particle0.Position[0], Particle1.Position[0], Particle0.Angle[0], Particle1.Angle[0], out double temp_Distance, out double[] temp_DistanceVector, out double[] temp_ClosestPoint_P0, out double[] temp_ClosestPoint_P1, out bool Overlapping_AfterReset);
                        if (Overlapping)
                            break;
                        if (temp_Distance < Distance)
                        {
                            Distance = temp_Distance;
                            DistanceVector = temp_DistanceVector;
                            ClosestPoint_P0 = temp_ClosestPoint_P0;
                            ClosestPoint_P1 = temp_ClosestPoint_P1;
                        }
                    }
                }

                // Ensure that the threshold is large enough
                Threshold = double.MaxValue;
                _FSI_Collision.CalculateNormalAndTangentialVector(DistanceVector, out NormalVector, out TangentialVector);
                _FSI_Collision.ProjectVelocity(NormalVector, TangentialVector, Particle0.TranslationalVelocity[0], out collisionVn_P0, out collisionVt_P0);
                _FSI_Collision.ProjectVelocity(NormalVector, TangentialVector, Particle1.TranslationalVelocity[0], out collisionVn_P1, out collisionVt_P1);
                ForceCollision = true;
            }
            if (double.IsNaN(collisionVn_P0) || double.IsInfinity(collisionVn_P0))
                throw new ArithmeticException("Error trying to update particle position. Value:  " + collisionVn_P0);
            double eps = Threshold.Pow2() / 2; // Turek paper
            double epsPrime = Threshold / 2; // Turek paper

            double[] collisionForce;

            Console.WriteLine("Distance: " + Distance);
            Console.WriteLine("Threshold: " + Threshold);
            Collided = false;

            // =======================================================
            // Step 6
            // Skip integration of hydrodynamic forces (maybe no 
            // longer necesarry, testing needed)
            // =======================================================
            if (Distance < Threshold)
            {
                Particle0.skipForceIntegration = true;
                Particle1.skipForceIntegration = true;
            }

            if (Distance > Threshold)
            {
                Particle0.m_collidedWithParticle[m_Particles.IndexOf(Particle1)] = false;
                Particle1.m_collidedWithParticle[m_Particles.IndexOf(Particle0)] = false;
                Particle0.m_closeInterfacePointTo[m_Particles.IndexOf(Particle1)] = null;
                Particle1.m_closeInterfacePointTo[m_Particles.IndexOf(Particle0)] = null;
                triggerOnlyCollisionProcedure = false;
                return;
            }

            // =======================================================
            // Step 7
            // Main collision procedure, only the momentum 
            // conservation model is fully supported.
            // =======================================================
            switch (CollisionModel)
            {
                case (FSI_Solver.FSI_Control.CollisionModel.RepulsiveForce):
                    if ((Distance <= Threshold))
                    {
                        DistanceVector.ScaleV((Threshold - Distance).Pow2());
                        DistanceVector.ScaleV(1 / eps);


                        collisionForce = DistanceVector;
                        var collisionForceP1 = collisionForce.CloneAs();
                        collisionForce.ScaleV(-100.0);
                        collisionForceP1.ScaleV(-100.0);
                        Particle0.HydrodynamicForces[0].AccV(-1, collisionForce);
                        //particle0.hydrodynTorqueAtIteration[0] += 100 * (collisionForce[0] * (tempPoint_P0[0] - particle0.positionAtIteration[0][0]) + collisionForce[1] * (tempPoint_P0[1] - particle0.positionAtIteration[0][1]));
                        Particle1.HydrodynamicForces[0].AccV(1, collisionForceP1);
                        //particle1.hydrodynTorqueAtIteration[0] += -100 * (collisionForceP1[0] * (tempPoint_P1[0] - particle1.positionAtIteration[0][0]) + collisionForceP1[1] * (tempPoint_P1[1] - particle1.positionAtIteration[0][1]));
                        Console.WriteLine("Collision information: Particles coming close, force " + collisionForce.L2Norm());
                        Console.WriteLine("Collision information: Particles coming close, torque " + Particle1.HydrodynamicTorque[0]);

                        if (Distance <= 1.5 * h_min)
                        {
                            Console.WriteLine("Entering overlapping loop....");
                            triggerOnlyCollisionProcedure = true;
                        }

                    }
                    throw new NotImplementedException("Please use the MomentumConservation model");
                    //break;


                case FSI_Control.CollisionModel.MomentumConservation:

                    if ((Distance <= Threshold || ForceCollision) && iteration_counter == 0)// && (!Particle0.m_collidedWithParticle[m_Particles.IndexOf(Particle1)] && !Particle1.m_collidedWithParticle[m_Particles.IndexOf(Particle0)] || iteration_counter != 0)))
                    {
                        // Bool if collided
                        Particle0.m_collidedWithParticle[m_Particles.IndexOf(Particle1)] = true;
                        Particle1.m_collidedWithParticle[m_Particles.IndexOf(Particle0)] = true;
                        Collided = true;

                        // Bool if force integration should be skipped
                        Particle0.skipForceIntegration = true;
                        Particle1.skipForceIntegration = true;

                        // coefficient of restitution (e=0 pastic; e=1 elastic)
                        double e = ((FSI_Control)Control).CoefficientOfRestitution;

                        // Calculate the position correction
                        double[] CollisionPositionCorrection0 = new double[SpatialDim];
                        double[] CollisionPositionCorrection1 = new double[SpatialDim];
                        for (int d=0; d< SpatialDim; d++)
                        {
                            double DistanceFraction0 = Math.Abs(collisionVn_P0) / (Math.Abs(collisionVn_P0) + Math.Abs(collisionVn_P1));
                            double DistanceFraction1 = Math.Abs(collisionVn_P1) / (Math.Abs(collisionVn_P0) + Math.Abs(collisionVn_P1));
                            CollisionPositionCorrection0[d] = -DistanceFraction0 * DistanceVector[d];
                            CollisionPositionCorrection1[d] = DistanceFraction1 * DistanceVector[d];
                        }
                        Particle0.CollisionPositionCorrection.Add(CollisionPositionCorrection0);
                        Particle1.CollisionPositionCorrection.Add(CollisionPositionCorrection1);
                        if (double.IsNaN(CollisionPositionCorrection0[0]) || double.IsInfinity(CollisionPositionCorrection0[0]))
                            throw new ArithmeticException("Error trying to update particle position. Value:  " + CollisionPositionCorrection0[0]);

                        // Calculate excentric parameter
                        double[] RadialDistance0 = Auxillary.VectorDiff(ClosestPoint_P0, Particle0.Position[0]);
                        double[] RadialDistance1 = Auxillary.VectorDiff(ClosestPoint_P1, Particle1.Position[0]);
                        double a0 = Particle0 is Particle_Sphere ? 0.0 : RadialDistance0[0] * TangentialVector[0] + RadialDistance0[1] * TangentialVector[1];
                        double a1 = Particle1 is Particle_Sphere ? 0.0 : RadialDistance1[0] * TangentialVector[0] + RadialDistance1[1] * TangentialVector[1];

                        // Calculate post collision velocities.
                        double Fx;
                        double Fxrot;
                        double tempCollisionVn_P0;
                        double tempCollisionVn_P1;
                        double tempCollisionRot_P0 = 0;
                        double tempCollisionRot_P1 = 0;
                        if (!Particle0.IncludeTranslation && !Particle0.IncludeRotation)
                        {
                            Fx = (1 + e) * ((collisionVn_P1) / (1 / Particle1.Mass_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                            Fxrot = (1 + e) * ((a1 * Particle1.RotationalVelocity[0]) / (1 / Particle1.Mass_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                            tempCollisionVn_P0 = collisionVn_P0;
                            tempCollisionVn_P1 = collisionVn_P1 + (Fx + Fxrot) / Particle1.Mass_P;
                            tempCollisionRot_P1 = Particle1.RotationalVelocity[0] - a1 * (Fx + Fxrot) / Particle1.MomentOfInertia_P;
                        }
                        else if (!Particle1.IncludeTranslation && !Particle1.IncludeRotation)
                        {
                            Fx = (1 + e) * ((collisionVn_P0) / (1 / Particle0.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P));
                            Fxrot = (1 + e) * ((-a0 * Particle0.RotationalVelocity[0]) / (1 / Particle0.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P));
                            tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / Particle0.Mass_P;
                            tempCollisionRot_P0 = Particle0.RotationalVelocity[0] + a0 * (Fx + Fxrot) / Particle0.MomentOfInertia_P;
                            tempCollisionVn_P1 = collisionVn_P1;
                        }
                        else
                        {
                            Fx = (1 + e) * ((collisionVn_P0 - collisionVn_P1) / (1 / Particle0.Mass_P + 1 / Particle1.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                            Fxrot = (1 + e) * ((-a0 * Particle0.RotationalVelocity[0] + a1 * Particle1.RotationalVelocity[0]) / (1 / Particle0.Mass_P + 1 / Particle1.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                            tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / Particle0.Mass_P;
                            tempCollisionRot_P0 = Particle0.RotationalVelocity[0] + a0 * (Fx + Fxrot) / Particle0.MomentOfInertia_P;
                            tempCollisionVn_P1 = collisionVn_P1 + (Fx + Fxrot) / Particle1.Mass_P;
                            tempCollisionRot_P1 = Particle1.RotationalVelocity[0] - a1 * (Fx + Fxrot) / Particle1.MomentOfInertia_P;
                        }

                        double tempCollisionVt_P0 = collisionVt_P0 * e;
                        double tempCollisionVt_P1 = collisionVt_P1 * e;
                        Console.WriteLine("a0:    " + a0 + "   Fx:    " + (-Fx) + "      Fxrot:    " + (-Fxrot));
                        Console.WriteLine("a1:    " + a1 + "   Fx:    " + Fx + "      Fxrot:    " + Fxrot);
                        
                        Particle0.CollisionNormal.Add(NormalVector);
                        Particle1.CollisionNormal.Add(NormalVector);
                        Particle0.CollisionTangential.Add(TangentialVector);
                        Particle1.CollisionTangential.Add(TangentialVector);
                        Particle0.CollisionRotationalVelocity.Add(tempCollisionRot_P0);
                        Particle1.CollisionRotationalVelocity.Add(tempCollisionRot_P1);
                        Particle0.CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn_P0, tempCollisionVt_P0 });
                        Particle1.CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn_P1, tempCollisionVt_P1 });
                        

                        for (int d = 0; d < 2; d++)
                        {
                            Particle0.TranslationalAcceleration[1][d] = 0;
                            Particle0.RotationalAcceleration[1] = 0;
                            Particle0.RotationalVelocity[1] = 0;
                            Particle1.TranslationalAcceleration[1][d] = 0;
                            Particle1.RotationalAcceleration[1] = 0;
                            Particle1.RotationalVelocity[1] = 0;
                        }
                    }
                    else
                    {
                        Particle0.m_collidedWithParticle[m_Particles.IndexOf(Particle1)] = false;
                        Particle1.m_collidedWithParticle[m_Particles.IndexOf(Particle0)] = false;
                        Particle0.m_closeInterfacePointTo[m_Particles.IndexOf(Particle1)] = null;
                        Particle1.m_closeInterfacePointTo[m_Particles.IndexOf(Particle0)] = null;
                    }
                    ForceCollision = false;
                    break;

                default:
                    throw new NotImplementedException("Collision model not available");
            }
        }

        public void GetWall(CellMask ParticleBoundaryCells, out double[,] WallPoints)
        {
            int SpatialDim = ParticleBoundaryCells.GridData.SpatialDimension;
            int NoOfMaxWallEdges = 4;
            WallPoints = new double[NoOfMaxWallEdges, SpatialDim];
            int[][] Cells2Edges = ((GridData)this.GridData).Cells.Cells2Edges;
            IList<Platform.LinAlg.AffineTrafo> trafo = GridData.iGeomEdges.Edge2CellTrafos;
            foreach (Chunk cnk in ParticleBoundaryCells)
            {
                for (int i = cnk.i0; i < cnk.JE; i++)
                {
                    foreach (int e in Cells2Edges[i])
                    {
                        int eId = (e < 0) ? -e - 1 : e - 1;
                        byte et = ((GridData)this.GridData).Edges.EdgeTags[eId];
                        if (GridData.EdgeTagNames[et].Contains("wall") || GridData.EdgeTagNames[et].Contains("Wall"))
                        {
                            int jCell = GridData.iGeomEdges.CellIndices[eId, 0];
                            int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

                            NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
                            NodeSet Nodes = refNodes.ElementAt(iKref);

                            int trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[eId, 0];
                            MultidimensionalArray transFormed = trafo[trafoIdx].Transform(Nodes);
                            MultidimensionalArray WallVerticies = transFormed.CloneAs();
                            GridData.TransformLocal2Global(transFormed, WallVerticies, jCell);
                            double[] WallPoint1 = WallVerticies.GetRow(0);
                            double[] WallPoint2 = WallVerticies.GetRow(1);
                            if (Math.Abs(WallPoint1[0] - WallPoint2[0]) < 1e-12)
                            {
                                if (WallPoints[0, 0] == 0 || Math.Abs(WallPoint1[0] - WallPoints[0, 0]) < 1e-12)
                                    WallPoints[0, 0] = WallPoint1[0];
                                else if (WallPoints[1, 0] == 0 || Math.Abs(WallPoint1[0] - WallPoints[1, 0]) < 1e-12)
                                    WallPoints[1, 0] = WallPoint1[0];
                                else
                                    throw new ArithmeticException("Error trying to get wall position. Please use horizontal/vertical boudaries");
                            }
                            if (Math.Abs(WallPoint1[1] - WallPoint2[1]) < 1e-12)
                            {
                                if (WallPoints[2, 1] == 0 || Math.Abs(WallPoint1[1] - WallPoints[2, 1]) < 1e-12)
                                    WallPoints[2, 1] = WallPoint1[1];
                                else if (WallPoints[3, 1] == 0 || Math.Abs(WallPoint1[1] - WallPoints[3, 1]) < 1e-12)
                                    WallPoints[3, 1] = WallPoint1[1];
                                else
                                    throw new ArithmeticException("Error trying to get wall position. Please use horizontal/vertical boudaries");
                            }
                        }
                    }
                }
            }
        }

        public void WallCollisionForcesNew(Particle particle, double dt)
        {
            if (CollisionModel == FSI_Control.CollisionModel.NoCollisionModel)
                return;

            int SpatialDim = GridData.SpatialDimension;
            // search for colored cells
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate();
            FSI_Collision _FSI_Collision = new FSI_Collision();
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            List<int[]> ColoredCellsSorted = levelSetUpdate.ColoredCellsFindAndSort(CellColor);
            int[] ParticleColorArray = levelSetUpdate.FindParticleColor(GridData, new List<Particle> { particle }, ColoredCellsSorted);
            CellMask ParticleCutCells = levelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, ParticleColorArray[0], J, false);
            // consider only colored boundary cells
            CellMask ParticleBoundaryCells = GridData.GetBoundaryCells().Intersect(ParticleCutCells);
            double Distance = double.MaxValue;
            double[] DistanceVec = new double[Grid.SpatialDimension];
            double[] ClosestPointParticle = new double[2];
            double[] ClosestPointWall = new double[2];
            GetWall(ParticleBoundaryCells, out double[,] WallPoints);
            double[] point0 = particle.Position[0].CloneAs();
            double[] point1 = point0.CloneAs();
            double[,] TempTranslationalVelocity = new double[WallPoints.GetLength(0), 2];
            double[] TempRotationalVelocity = new double[WallPoints.GetLength(0)];
            double[,] NormalVectors = new double[WallPoints.GetLength(0), 2];
            double[,] TangentialVectors = new double[WallPoints.GetLength(0), 2];
            int denominator = 0;
            for (int i = 0; i < WallPoints.GetLength(0); i++)
            {
                if (WallPoints[i, 0] != 0)
                    point1[0] = WallPoints[i, 0];
                else if (WallPoints[i, 1] != 0)
                    point1[1] = WallPoints[i, 1];
                else
                    continue;

                bool Overlapping = false;
                int test = particle.NoOfSubParticles();
                for (int j = 0; j < particle.NoOfSubParticles(); j++)
                {
                    Collision.GJK_DistanceAlgorithm(particle, j, null, 0, LsTrk, point0, point1, particle.Angle[0], 0, out double temp_Distance, out double[] temp_DistanceVector, out double[]  temp_ClosestPoint_P0, out double[]  temp_ClosestPoint_P1, out Overlapping); ;
                    if (Overlapping)
                        break;
                    if (temp_Distance < Distance)
                    {
                        Distance = temp_Distance;
                        DistanceVec = temp_DistanceVector;
                        ClosestPointParticle = temp_ClosestPoint_P0;
                        ClosestPointWall = temp_ClosestPoint_P1;
                    }
                }
                _FSI_Collision.CalculateNormalAndTangentialVector(DistanceVec, out double[] normal, out double[] tangential);
                _FSI_Collision.CalculateDynamicCollisionThreshold(particle, null, ClosestPointParticle, ClosestPointWall, normal, Distance, dt, out double threshold);
                _FSI_Collision.ProjectVelocity(normal, tangential, particle.TranslationalVelocity[0], out double collisionVn_P0, out double collisionVt_P0);

                if (Overlapping)// && particle0.m_closeInterfacePointTo[m_Particles.IndexOf(particle1)] != null && particle1.m_closeInterfacePointTo[m_Particles.IndexOf(particle0)] != null)
                {
                    particle.Position[0] = particle.Position[1].CloneAs();
                    particle.Angle[0] = particle.Angle[1];
                    point0 = particle.Position[0].CloneAs();
                    point1 = point0.CloneAs();
                    if (WallPoints[i, 0] != 0)
                        point1[0] = WallPoints[i, 0];
                    else if (WallPoints[i, 1] != 0)
                        point1[1] = WallPoints[i, 1];
                    for (int j = 0; j < particle.NoOfSubParticles(); i++)
                    {
                        Collision.GJK_DistanceAlgorithm(particle, j, null, 0, LsTrk, point0, point1, particle.Angle[0], 0, out Distance, out DistanceVec, out ClosestPointParticle, out ClosestPointWall, out Overlapping); ;
                    }
                    threshold = 1e20;
                    _FSI_Collision.CalculateNormalAndTangentialVector(DistanceVec, out normal, out tangential);
                    _FSI_Collision.ProjectVelocity(normal, tangential, particle.TranslationalVelocity[0], out collisionVn_P0, out collisionVt_P0);
                    threshold = 1e20;
                    _FSI_Collision.CalculateNormalAndTangentialVector(DistanceVec, out normal, out tangential);
                    _FSI_Collision.ProjectVelocity(normal, tangential, particle.TranslationalVelocity[0], out collisionVn_P0, out collisionVt_P0);
                }
                #region oldstuff
                //foreach (int iEdge in ColoredCellsGrid.BoundaryEdgesMask.ItemEnum)
                //{
                //    // Collision forces have to act
                //    if (GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge))
                //    {
                //        collision = true;
                //        int jCell = GridData.iGeomEdges.CellIndices[iEdge, 0];
                //        int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

                //        NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
                //        NodeSet Nodes = refNodes.ElementAt(iKref);

                //        int trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[iEdge, 0];
                //        MultidimensionalArray transFormed = trafo[trafoIdx].Transform(Nodes);
                //        MultidimensionalArray WallVerticies = transFormed.CloneAs();
                //        double[,] fsdf = WallVerticies.To2DArray();
                //        GridData.TransformLocal2Global(transFormed, WallVerticies, jCell);
                //        //double[] point0 = particle.Position[0].CloneAs();
                //        double[] point1 = WallVerticies.GetRow(0);
                //        double[] WallPoint1 = WallVerticies.GetRow(0);
                //        double[] WallPoint2 = WallVerticies.GetRow(1);
                //        _FSI_Auxillary.Wall_GJK_DistanceAlgorithm(particle, WallVerticies, LsTrk, point0, point1, 2, out double tempDistance, out double[] temp_distanceVec, out double[] temp_ClosestPointParticle, out double[] temp_ClosestPointWall, out Overlapping);
                //        if (Math.Abs(WallPoint1[0] - WallPoint2[0]) < 1e-12)
                //        {
                //            temp_ClosestPointWall[0] = WallPoint1[0];
                //            temp_ClosestPointWall[1] = temp_ClosestPointParticle[1];
                //            tempDistance = Math.Abs(temp_ClosestPointParticle[0] - temp_ClosestPointWall[0]);
                //            for (int d = 0; d < 2; d++)
                //            {
                //                temp_distanceVec[d] = temp_ClosestPointWall[d] - temp_ClosestPointParticle[d];
                //            }
                //        }
                //        if(Math.Abs(WallPoint1[1] - WallPoint2[1]) < 1e-12)
                //        {
                //            temp_ClosestPointWall[1] = WallPoint1[1];
                //            temp_ClosestPointWall[0] = temp_ClosestPointParticle[0];
                //            tempDistance = Math.Abs(temp_ClosestPointParticle[1] - temp_ClosestPointWall[1]);
                //            for (int d = 0; d < 2; d++)
                //            {
                //                temp_distanceVec[d] = temp_ClosestPointWall[d] - temp_ClosestPointParticle[d];
                //            }
                //        }
                //        if (tempDistance < distance)
                //        {
                //            distance = tempDistance;
                //            ClosestPointParticle = temp_ClosestPointParticle.CloneAs();
                //            ClosestPointWall = temp_ClosestPointWall.CloneAs();
                //            distanceVec = temp_distanceVec.CloneAs();
                //        }
                //    }
                //}
                //_FSI_Collision.FindNormalAndTangentialVector(distanceVec, out double[] normal, out double[] tangential);
                //_FSI_Collision.CalculateDynamicCollisionThreshold(particle, null, ClosestPointParticle, ClosestPointWall, normal, distance, dt, out double threshold);
                //_FSI_Collision.ProjectVelocity(normal, tangential, particle.TranslationalVelocity[0], out double collisionVn_P0, out double collisionVt_P0);
                //if (!Overlapping && threshold == 0)
                //{
                //    bool Overlapping_NextTimestep = false;
                //    _FSI_Collision.PredictParticleNextTimestep(particle, 2, dt, out double[] VirtualPosition0, out double[] VirtualVelocity0, out double VirtualAngle0, out double VirtualRotationalVelocity0);
                //    foreach (int iEdge in ColoredCellsGrid.BoundaryEdgesMask.ItemEnum)
                //    {
                //        // Collision forces have to act
                //        if (GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge))
                //        {
                //            collision = true;
                //            int jCell = GridData.iGeomEdges.CellIndices[iEdge, 0];
                //            int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

                //            NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
                //            NodeSet Nodes = refNodes.ElementAt(iKref);

                //            int trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[iEdge, 0];
                //            MultidimensionalArray transFormed = trafo[trafoIdx].Transform(Nodes);
                //            MultidimensionalArray WallVerticies = transFormed.CloneAs();
                //            GridData.TransformLocal2Global(transFormed, WallVerticies, jCell);
                //            double[] point0 = particle.Position[0].CloneAs();
                //            double[] point1 = WallVerticies.GetRow(0);
                //            _FSI_Auxillary.Wall_GJK_DistanceAlgorithm(particle, WallVerticies, LsTrk, point0, point1, 2, out double tempDistance, out distanceVec, out ClosestPointParticle, out ClosestPointWall, out Overlapping_NextTimestep);
                //            if (tempDistance < distance)
                //            {
                //                distance = tempDistance;
                //            }
                //        }
                //    }
                //    _FSI_Collision.FindNormalAndTangentialVector(distanceVec, out normal, out tangential);
                //    _FSI_Collision.FindRadialVector(VirtualPosition0, ClosestPointParticle, out _, out double RadialLength0, out double[] RadialNormalVector0);
                //    _FSI_Collision.TransformRotationalVelocity(VirtualRotationalVelocity0, RadialLength0, RadialNormalVector0, out double[] PointVelocityDueToRotation0);
                //    double[] PointVelocity0 = new double[2];
                //    for (int d = 0; d < 2; d++)
                //    {
                //        PointVelocity0[d] = VirtualVelocity0[d] + PointVelocityDueToRotation0[d];
                //    }
                //    _FSI_Collision.ProjectVelocityOnVector(normal, VirtualVelocity0, out double DetectCollisionVn_P0);
                //    _FSI_Collision.ProjectVelocity(normal, tangential, VirtualVelocity0, out collisionVn_P0, out collisionVt_P0);
                //    if (distance <= Math.Abs((-DetectCollisionVn_P0 + 0) * dt))
                //    {
                //        threshold = Math.Abs((-DetectCollisionVn_P0 + 0) * dt);
                //    }
                //    if (Overlapping_NextTimestep)
                //    {
                //        threshold = 1e20;
                //    }
                //}
                //if (Overlapping)// && particle0.m_closeInterfacePointTo[m_Particles.IndexOf(particle1)] != null && particle1.m_closeInterfacePointTo[m_Particles.IndexOf(particle0)] != null)
                //{
                //    particle.Position[0] = particle.Position[1].CloneAs();
                //    foreach (int iEdge in ColoredCellsGrid.BoundaryEdgesMask.ItemEnum)
                //    {
                //        // Collision forces have to act
                //        if (GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge))
                //        {
                //            collision = true;
                //            int jCell = GridData.iGeomEdges.CellIndices[iEdge, 0];
                //            int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

                //            NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
                //            NodeSet Nodes = refNodes.ElementAt(iKref);

                //            int trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[iEdge, 0];
                //            MultidimensionalArray transFormed = trafo[trafoIdx].Transform(Nodes);
                //            MultidimensionalArray WallVerticies = transFormed.CloneAs();
                //            GridData.TransformLocal2Global(transFormed, WallVerticies, jCell);
                //            double[] point0 = particle.Position[0].CloneAs();
                //            double[] point1 = WallVerticies.GetRow(0);
                //            _FSI_Auxillary.Wall_GJK_DistanceAlgorithm(particle, WallVerticies, LsTrk, point0, point1, 2, out double tempDistance, out distanceVec, out ClosestPointParticle, out ClosestPointWall, out bool _);
                //            if (tempDistance < distance)
                //            {
                //                distance = tempDistance;
                //            }
                //        }
                //    }
                //    threshold = 1e20;
                //    _FSI_Collision.FindNormalAndTangentialVector(distanceVec, out normal, out tangential);
                //    _FSI_Collision.ProjectVelocity(normal, tangential, particle.TranslationalVelocity[0], out collisionVn_P0, out collisionVt_P0);
                //}
                #endregion

                //if (collision == false)
                //{
                //    Console.WriteLine("Reset Wall");
                //    particle.m_collidedWithWall[0] = false;
                //    return;
                //}

                Console.WriteLine("Closest Distance to wall is: " + Distance);
                Console.WriteLine("Wall threshold: " + threshold);
                double eps = threshold.Pow2() / 2; // Turek paper
                double epsPrime = threshold / 2; // Turek paper

                //double[] collisionForce;
                switch (CollisionModel)
                {
                    case (FSI_Solver.FSI_Control.CollisionModel.RepulsiveForce):
                        //if ((Distance <= threshold))
                        //{
                        //    Console.WriteLine("Strongly recommended to use conservation of momentum collision model. This one is highly experimental!!!!");

                        //    // Modell 1
                        //    distanceVec.ScaleV(1 / eps);
                        //    distanceVec.ScaleV(((threshold - realDistance).Abs()));

                        //    collisionForce = distanceVec;
                        //    collisionForce.ScaleV(100.0);

                        //    particle.HydrodynamicForces[0] = collisionForce;
                        //    throw new NotImplementedException("The repulsive force model is not parallelized, please use the momentum conservation model.");
                        //    //return;
                        //}


                        //if (Distance <= (1.5 * hmin))
                        //{

                        //    distanceVec.ScaleV((threshold - realDistance).Abs());
                        //    distanceVec.ScaleV(1 / epsPrime);
                        //    collisionForce = distanceVec;

                        //    collisionForce.ScaleV(100.0);
                        //    particle.HydrodynamicForces[0].AccV(1, collisionForce);
                        //    particle.HydrodynamicTorque[0] -= (collisionForce[0] * (ClosestPointParticle[0] - particle.Position[0][0]) + collisionForce[1] * (ClosestPointParticle[1] - particle.Position[0][1]));
                        //    Console.WriteLine("Collision information: Wall overlapping, force X " + collisionForce[0]);
                        //    Console.WriteLine("Collision information: Wall overlapping, force Y " + collisionForce[1]);

                        //    if (realDistance <= 1.5 * hmin)
                        //    {
                        //        Console.WriteLine("Entering wall overlapping loop....");
                        //        triggerOnlyCollisionProcedure = true;

                        //    }
                        //    return;
                        //}
                        break;

                    case (FSI_Solver.FSI_Control.CollisionModel.MomentumConservation):

                        if (Distance <= (threshold) && !particle.m_collidedWithWall[0])
                        {
                            Console.WriteLine("I'm trying to calculate the wand collision.");
                            //coefficient of restitution (e=0 pastic; e=1 elastic)

                            double e = 1.0;

                            // Fully plastic for bottom wall
                            if (particle.Position[0][1] < 0.5 && ((FSI_Control)Control).LowerWallFullyPlastic)
                                e = 0.0;

                            // if particle already collided with wall
                            particle.m_collidedWithWall[0] = true;

                            // Skip force integration for next timestep
                            particle.skipForceIntegration = true;


                            collisionVn_P0 = particle.TranslationalVelocity[0][0] * normal[0] + particle.TranslationalVelocity[0][1] * normal[1];
                            collisionVt_P0 = particle.TranslationalVelocity[0][0] * tangential[0] + particle.TranslationalVelocity[0][1] * tangential[1];


                            // exzentric collision
                            // ----------------------------------------
                            ClosestPointParticle.AccV(-1, particle.Position[0]);
                            double a0 = 0;
                            if (particle is Particle_Sphere)
                                a0 = 0.0;
                            else
                                a0 = (ClosestPointParticle[0] * tangential[0] + ClosestPointParticle[1] * tangential[1]);
                            Console.WriteLine("a0: " + a0);



                            double Fx = (1 + e) * (collisionVn_P0) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);
                            Console.WriteLine("Fx: " + Fx);
                            double Fxrot = (1 + e) * (-a0 * particle.RotationalVelocity[0]) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);
                            Console.WriteLine("Fxrot: " + Fxrot);

                            double tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / particle.Mass_P;
                            Console.WriteLine("tempCollisionVn_P0: " + tempCollisionVn_P0);
                            double tempCollisionVt_P0 = collisionVt_P0;
                            Console.WriteLine("tempCollisionVt_P0: " + tempCollisionVt_P0);

                            //particle.RotationalVelocity[0] = particle.RotationalVelocity[0] + a0 * (Fx + Fxrot) / particle.MomentOfInertia_P;
                            TempRotationalVelocity[i] = particle.RotationalVelocity[0] + a0 * (Fx + Fxrot) / particle.MomentOfInertia_P;
                            for (int d = 0; d < SpatialDim; d++)
                            {
                                TempTranslationalVelocity[i, d] = normal[d] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[d];
                                NormalVectors[i, d] = normal[d];
                                TangentialVectors[i, d] = tangential[d];
                            }
                            denominator += 1;
                            //particle.TranslationalVelocity[0] = new double[] { normal[0] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[0], normal[1] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[1] };

                        }

                        if (Distance > threshold && particle.m_collidedWithWall[0])
                        {
                            Console.WriteLine("Reset Wall");
                            particle.m_collidedWithWall[0] = false;
                        }
                        break;

                    default:
                        throw new NotImplementedException("Collision model not available");
                }
            }
            if (denominator != 0)
            {
                particle.RotationalVelocity[0] = 0;
                for (int r = 0; r < TempRotationalVelocity.Length; r++)
                {
                    particle.RotationalVelocity[0] += TempRotationalVelocity[r];
                }

                double[] Normal = new double[SpatialDim];
                double[] Tangential = new double[SpatialDim];
                for (int t = 0; t < TempTranslationalVelocity.GetLength(0); t++)
                {
                    if (TempTranslationalVelocity[t, 0] == 0 && TempTranslationalVelocity[t, 1] == 0)
                        continue;
                    for (int d = 0; d < SpatialDim; d++)
                    {
                        Normal[d] += NormalVectors[t, d];
                        Tangential[d] += TangentialVectors[t, d];
                    }
                }
                Normal.ScaleV(1 / Math.Sqrt(Normal[0].Pow2() + Normal[1].Pow2()));
                Tangential.ScaleV(1 / Math.Sqrt(Tangential[0].Pow2() + Tangential[1].Pow2()));
                double temp_NormalVel = 0;
                double temp_TangentialVel = 0;
                for (int t = 0; t < TempTranslationalVelocity.GetLength(0); t++)
                {
                    temp_NormalVel += TempTranslationalVelocity[t, 0] * Normal[0] + TempTranslationalVelocity[t, 1] * Normal[1];
                    temp_TangentialVel += TempTranslationalVelocity[t, 0] * Tangential[0] + TempTranslationalVelocity[t, 1] * Tangential[1];
                }
                temp_NormalVel /= denominator;
                temp_TangentialVel /= denominator;
                for (int d = 0; d < SpatialDim; d++)
                {
                    particle.TranslationalVelocity[0][d] = Normal[d] * temp_NormalVel + Tangential[d] * temp_TangentialVel;
                }
            }
        }

        /// <summary>
        /// Calculation of collision forces between particle and wall
        /// </summary>
        //public void WallCollisionForces(Particle particle, int ParticleID, double hmin)
        //{
        //    if (CollisionModel == FSI_Control.CollisionModel.NoCollisionModel)
        //        return;

        //    int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
        //    FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate();
        //    List<int[]> ColoredCellsSorted = levelSetUpdate.ColoredCellsFindAndSort(CellColor);
        //    List<Particle> temp = new List<Particle> { particle };
        //    int[] ParticleColorArray = levelSetUpdate.FindParticleColor(GridData, temp, ColoredCellsSorted);
        //    CellMask particleCutCells = levelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, ParticleColorArray[0], J, false);

        //    //var particleCutCellArray = particleCutCells.ItemEnum.ToArray();
        //    //var neighborCellsArray = particleCutCells.AllNeighbourCells().ItemEnum.ToArray();
        //    //var allCellsArray = particleCutCellArray.Concat(neighborCellsArray).ToArray();
        //    //var allCells = new CellMask(GridData, neighborCellsArray);
        //    CellMask allCells = particleCutCells;

        //    collision = false;

        //    double distance = double.MaxValue;
        //    double[] distanceVec = new double[Grid.SpatialDimension];

        //    // All interface points at a specific subgrid containing all cut cells of one particle
        //    MultidimensionalArray interfacePoints = null;

        //    //Console.WriteLine("ParticleCutCellCount:   " + particleCutCells.Count());

        //    IList<Platform.LinAlg.AffineTrafo> trafo = GridData.iGeomEdges.Edge2CellTrafos;

        //    SubGrid allCellsGrid = new SubGrid(allCells);

        //    double[] tempPoint = new double[2] { 0.0, 0.0 };

        //    foreach (int iEdge in allCellsGrid.BoundaryEdgesMask.ItemEnum)
        //    {

        //        // Collision forces have to act
        //        if (GridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge))
        //        {

        //            if (interfacePoints == null)
        //                interfacePoints = particle.GetSurfacePoints(LsTrk, particle.Position[0], particle.Angle[0]);

        //            collision = true;
        //            int jCell = GridData.iGeomEdges.CellIndices[iEdge, 0];
        //            int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

        //            NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
        //            NodeSet Nodes = refNodes.ElementAt(iKref);

        //            int trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[iEdge, 0];
        //            MultidimensionalArray transFormed = trafo[trafoIdx].Transform(Nodes);
        //            MultidimensionalArray newVertices = transFormed.CloneAs();
        //            GridData.TransformLocal2Global(transFormed, newVertices, jCell);
        //            var tempDistance = 0.0;
                    
        //            for (int i = 0; i < interfacePoints.NoOfRows; i++)
        //            {
        //                for (int j = 0; j < newVertices.NoOfRows; j++)
        //                {
        //                    tempDistance = Math.Sqrt((interfacePoints.GetRow(i)[0] - newVertices.GetRow(j)[0]).Pow2() + (interfacePoints.GetRow(i)[1] - newVertices.GetRow(j)[1]).Pow2());
        //                    if (tempDistance < distance)
        //                    {
        //                        tempPoint = interfacePoints.GetRow(i);
        //                        distanceVec = interfacePoints.GetRow(i).CloneAs();
        //                        distanceVec.AccV(-1, newVertices.GetRow(j));
        //                        distance = tempDistance;
        //                    }

        //                }


        //            }
        //        }
        //    }

        //    double realDistance = distance;

        //    if (collision == false)
        //    {
        //        Console.WriteLine("Reset Wall");
        //        particle.m_collidedWithWall[0] = false;
        //        return;
        //    }


        //    Console.WriteLine("Closest Distance to wall is: " + distance);

        //    double threshold = 1.5 * hmin; // was 1.5 * hmin
        //    Console.WriteLine("threshold: " + threshold);
        //    double eps = threshold.Pow2() / 2; // Turek paper
        //    double epsPrime = threshold / 2; // Turek paper

        //    double[] collisionForce;


        //    switch (CollisionModel)
        //    {
        //        case (FSI_Solver.FSI_Control.CollisionModel.RepulsiveForce):
        //            if ((realDistance <= threshold))
        //            {
        //                Console.WriteLine("Strongly recommended to use conservation of momentum collision model. This one is highly experimental!!!!");
                        
        //                // Modell 1
        //                distanceVec.ScaleV(1 / eps);
        //                distanceVec.ScaleV(((threshold - realDistance).Abs()));

        //                collisionForce = distanceVec;
        //                collisionForce.ScaleV(100.0);

        //                particle.HydrodynamicForces[0] = collisionForce;
        //                throw new NotImplementedException("The repulsive force model is not parallelized, please use the momentum conservation model.");
        //                //return;
        //            }


        //            if (realDistance <= (1.5 * hmin))
        //            {

        //                distanceVec.ScaleV((threshold - realDistance).Abs());
        //                distanceVec.ScaleV(1 / epsPrime);
        //                collisionForce = distanceVec;

        //                collisionForce.ScaleV(100.0);
        //                particle.HydrodynamicForces[0].AccV(1, collisionForce);
        //                particle.HydrodynamicTorque[0] -= (collisionForce[0] * (tempPoint[0] - particle.Position[0][0]) + collisionForce[1] * (tempPoint[1] - particle.Position[0][1]));
        //                Console.WriteLine("Collision information: Wall overlapping, force X " + collisionForce[0]);
        //                Console.WriteLine("Collision information: Wall overlapping, force Y " + collisionForce[1]);

        //                if (realDistance <= 1.5 * hmin)
        //                {
        //                    Console.WriteLine("Entering wall overlapping loop....");
        //                    triggerOnlyCollisionProcedure = true;

        //                }
        //                return;
        //            }
        //            break;

        //        case (FSI_Solver.FSI_Control.CollisionModel.MomentumConservation):

        //            if (realDistance <= (threshold) && !particle.m_collidedWithWall[0])
        //            {
        //                Console.WriteLine("I'm trying to calculate the wand collision.");
        //                //coefficient of restitution (e=0 pastic; e=1 elastic)

        //                double e = 1.0;

        //                // Fully plastic for bottom wall
        //                 if (particle.Position[0][1] < 0.5 && ((FSI_Control)Control).LowerWallFullyPlastic)
        //                    e = 0.0;

        //                // if particle already collided with wall
        //                particle.m_collidedWithWall[0] = true;

        //                // Skip force integration for next timestep
        //                particle.skipForceIntegration = true;

        //                //collision Nomal
        //                var normal = distanceVec.CloneAs();
        //                normal.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));
        //                double[] tangential = new double[] { -normal[1], normal[0] };


        //                double collisionVn_P0 = particle.TranslationalVelocity[0][0] * normal[0] + particle.TranslationalVelocity[0][1] * normal[1];
        //                Console.WriteLine("collisionVn_P0: " + collisionVn_P0);
        //                double collisionVt_P0 = particle.TranslationalVelocity[0][0] * tangential[0] + particle.TranslationalVelocity[0][1] * tangential[1];
        //                Console.WriteLine("collisionVt_P0: " + collisionVt_P0);


        //                // exzentric collision
        //                // ----------------------------------------
        //                tempPoint.AccV(-1, particle.Position[0]);
        //                Console.WriteLine("tempPoint: " + tempPoint[0]);
        //                Console.WriteLine("tempPoint: " + tempPoint[1]);
        //                double a0 = (tempPoint[0] * tangential[0] + tempPoint[1] * tangential[1]);
        //                Console.WriteLine("a0: " + a0);

        //                if (particle is Particle_Sphere)
        //                    a0 = 0.0;

        //                double Fx = (1 + e) * (collisionVn_P0) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);
        //                Console.WriteLine("Fx: " + Fx);
        //                double Fxrot = (1 + e) * (-a0 * particle.RotationalVelocity[0]) / (1 / particle.Mass_P + a0.Pow2() / particle.MomentOfInertia_P);
        //                Console.WriteLine("Fxrot: " + Fxrot);

        //                double tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / particle.Mass_P;
        //                Console.WriteLine("tempCollisionVn_P0: " + tempCollisionVn_P0);
        //                double tempCollisionVt_P0 = collisionVt_P0;
        //                Console.WriteLine("tempCollisionVt_P0: " + tempCollisionVt_P0);

        //                particle.RotationalVelocity[0] = particle.RotationalVelocity[0] + a0 * (Fx + Fxrot) / particle.MomentOfInertia_P;

        //                particle.TranslationalVelocity[0] = new double[] { normal[0] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[0], normal[1] * tempCollisionVn_P0 + tempCollisionVt_P0 * tangential[1] };
                        
        //            }

        //            if (realDistance > threshold && particle.m_collidedWithWall[0])
        //            {
        //                Console.WriteLine("Reset Wall");
        //                particle.m_collidedWithWall[0] = false;
        //            }
        //            break;

        //        default:
        //            throw new NotImplementedException("Collision model not available");
        //    }

        //}


        /// <summary>
        /// Mesh refinement
        /// Very primitive refinement indicator, works on a LevelSet criterion.
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel)
        {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate();
            CellMask ColoredCellMask = null;
            List<int[]> ColoredCellsSorted = levelSetUpdate.ColoredCellsFindAndSort(CellColor);
            int[] ParticleColorArray = levelSetUpdate.FindParticleColor(GridData, m_Particles, ColoredCellsSorted);
            for (int p = 0; p < ParticleColorArray.Length; p++)
            {
                if (ParticleColorArray[p] != 0)
                {
                    ColoredCellMask = levelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, ParticleColorArray[p], J, false);
                }
            }
            CellMask LevSetCells = LsTrk.Regions.GetCutCellMask();
            //CellMask LevSetNeighbours = LsTrk.Regions.GetNearFieldMask(1);
            int DesiredLevel_j = 0;
            if (ColoredCellMask != null && LevSetCells.Contains(j))
            {
                DesiredLevel_j = ((FSI_Control)this.Control).RefinementLevel;
            }
            //else if (LevSetNeighbours.Contains(j))
            //{
            //    DesiredLevel_j = 1;
            //}

            return DesiredLevel_j;
        }


        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid)
        {

            if (((FSI_Control)Control).AdaptiveMeshRefinement && iteration_counter == 0)
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