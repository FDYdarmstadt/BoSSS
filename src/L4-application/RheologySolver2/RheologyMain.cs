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
using System.Diagnostics;
using System.Linq;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;

using BoSSS.Platform;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid.Aggregation;

using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.Gnuplot;

using MPI.Wrappers;
using NUnit.Framework;
using BoSSS.Foundation.SpecFEM;
using System.IO;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Solver for calculation with viscoelastic extra stresses using the Oldroyd B model or the upper convected Maxwell model (UCM).
    /// </summary>
    public class Rheology : BoSSS.Solution.Application<RheologyControl> {
        static void Main(string[] args) {
            //RheologyTestProgram.Init();
            //DeleteOldPlotFiles();
            //RheologyTestProgram.ScalingChannelTestStokesCondition(2, 1.0);
            //RheologyTestProgram.Cleanup();
            //Assert.IsTrue(false, "Testcode left in main routine.");

            Rheology._Main(args, false, () => new Rheology());
        }

        #region instantiation

        // Attributes for fields (Names), initialization of DG fields
        //==============================================================

        /// <summary>
        /// Velocity domain
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            null,
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorFieldHistory<SinglePhaseField> Velocity;

        /// <summary>
        /// Velocities codomain: Residuum in momentum equation
        /// </summary>
        [InstantiateFromControlFile(new string[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> ResidualMomentum;

        /// <summary>
        /// Pressure domain
        /// </summary>
        [InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField Pressure;

        /// <summary>
        /// Pressure codomain: Residuum in continuity equation
        /// </summary>
        [InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualConti;

        ///// <summary>
        ///// Extra stress domain (2D): StressXX
        ///// </summary>
        //[InstantiateFromControlFile("StressXX", null, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField StressXX;

        /// <summary>
        /// Extra stress domain (2D): StressXX
        /// </summary>
        [InstantiateFromControlFile(VariableNames.StressXX, null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXX;

        /// <summary>
        /// Extra stress domain (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile(VariableNames.StressXY, null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXY;

        /// <summary>
        /// Extra stress domain (2D): StressYY
        /// </summary>
        [InstantiateFromControlFile(VariableNames.StressYY, null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressYY;

        /// <summary>
        /// Extra stress codomain (2D): StressXX
        /// </summary>
        [InstantiateFromControlFile("ResidualStressXX", VariableNames.StressXX, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressXX;

        /// <summary>
        /// Extra stresses codomain (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("ResidualStressXY", VariableNames.StressXY, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressXY;

        /// <summary>
        /// Extra stresses codomain (2D): StressYY
        /// </summary>
        [InstantiateFromControlFile("ResidualStressYY", VariableNames.StressYY, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressYY;

        /// <summary>
        /// Extra stresses parameter (2D): StressXX
        /// </summary>
        [InstantiateFromControlFile("StressXXP", VariableNames.StressXX, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXXP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("StressXYP", VariableNames.StressXY, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXYP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("StressYYP", VariableNames.StressYY, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressYYP;

        // Parameters: Velocity Gradient
        VectorField<SinglePhaseField> VelocityXGradient;
        VectorField<SinglePhaseField> VelocityYGradient;


        //Parameters: external analytical velocity
        SinglePhaseField U;
        SinglePhaseField V;

        

        // Some initialisation of variables
        //============================================
        IncompressibleBoundaryCondMap BcMap;

        /// <summary>
        /// Current Weissenberg number,
        /// if the solver is used in Weissenberg-increment mode
        /// (<see cref="RheologyControl.RaiseWeissenberg"/>, <see cref="RheologyControl.WeissenbergIncrement"/>)
        /// </summary>
        public double currentWeissenberg {
            get {
                return ((double) XOP.UserDefinedValues["Weissenbergnumber"]);
            }
            set {
                double odVal;
                if(XOP.UserDefinedValues.ContainsKey("Weissenbergnumber"))
                    odVal = currentWeissenberg;
                else
                    odVal = double.NegativeInfinity;
                                    
                if(odVal != value)
                    Console.WriteLine("setting Weissenberg Number to " + value);
                XOP.UserDefinedValues["Weissenbergnumber"] = value;
            }
        }

        /// <summary>
        /// restart value for the Weissenberg-increment mode
        /// </summary>
        private double restartWeissenberg = 0.0;


        bool ChangeMesh = true;

        /// <summary>
        /// Spatial operator 
        /// </summary>
        DifferentialOperator XOP;

     
        /// <summary>
        /// Timestepping object
        /// </summary>
        protected XdgTimestepping m_Timestepper;

        // Persson sensor and artificial viscosity
        //=============================================
        /// <summary>
        /// Instance of Persson sensor
        /// </summary>
        protected PerssonSensor perssonsensor;

        /// <summary>
        /// Instance of artificial viscosity
        /// </summary>
        protected SinglePhaseField artificalViscosity;

        /// <summary>
        /// Instance of max value of artificial viscosity
        /// </summary>
        protected double artificialMaxViscosity;


        // Settings for calculation
        //===============================================
        /// <summary>
        /// Set true if Navier Stokes is solved, then the mean velocities as parameters for calculation of convective terms are needed
        /// </summary>
        protected bool U0MeanRequired {
            get {
                return (!this.Control.Stokes);
            }
        }

             CoordinateVector m_CurrentSolution = null;

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.Velocity.Current, this.Pressure, this.StressXX, this.StressXY, this.StressYY));
                }
                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(this.ResidualMomentum, this.ResidualConti, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY));
                }
                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// DG Field instantiation.
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();

            /*
            if (Control.UsePerssonSensor == true) {
                perssonsensor = new PerssonSensor(StressXX);
                this.IOFields.Add(perssonsensor.GetField());
            }
            if (Control.UseArtificialDiffusion == true) {
                artificalViscosity = new SinglePhaseField(new Basis(GridData, 1), "artificalViscosity");
                this.IOFields.Add(artificalViscosity);

            }
            */
        }


        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public override void DataBackupBeforeBalancing(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            m_Timestepper.DataBackupBeforeBalancing(L);
        }

        #endregion

        /// <summary>
        /// Initialize Calculation, Create Equations
        /// </summary>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            int D = this.GridData.SpatialDimension;

            if (XOP != null && L == null && Control.Weissenberg == 0.0)
                return;

            if (m_Timestepper != null) {
                if (L != null) {

                    m_Timestepper.DataRestoreAfterBalancing(L,
                        ArrayTools.Cat<DGField>(Velocity.Current, Pressure, StressXX, StressXY, StressYY),
                        ArrayTools.Cat<DGField>(ResidualMomentum, ResidualConti, ResidualStressXX, ResidualStressXY, ResidualStressYY),
                        null, this.MultigridSequence, XOP);
                    this.LsTrk = m_Timestepper.LsTrk;

                    m_CurrentSolution = null;
                    m_CurrentResidual = null;

                }
            } else {

                using (new FuncTrace()) {

                    D = this.GridData.SpatialDimension;
                    BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Viscoelastic);

                    string[] CodName = new string[] { "momX", "momY", "div", "constitutiveXX", "constitutiveXY", "constitutiveYY" };

                    string[] Params = new string[] { };

                    if (this.Control.useFDJacobianForOperatorMatrix == true) {
                        Params = ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector(), VariableNames.StressXXP, VariableNames.StressXYP, VariableNames.StressYYP, "artificialViscosity");
                    } else {
                        //Params = this.Control.UseArtificialDiffusion ? new string[] { "artificialViscosity" } : null;
                        Params = null;
                    }

                    string[] DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure, VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY);

                    XOP = new DifferentialOperator(DomName, Params, CodName, QuadOrderFunc.NonLinearWithoutParameters(2));

                    // Development switches to turn specific components on or off, 
                    // for the sake of iterative solver testing:
                    bool MomContinuitycoupling = true;
                    bool ConstitutiveEqs = true;
                                       

                    // Temporal operator
                    // ================================================================================

                    var tempOp = new ConstantTemporalOperator(XOP, 0.0);
                    double rho = 1; // this.Control.PhysicalParameters.rho_A;
                    tempOp.SetDiagonal(0, rho);
                    tempOp.SetDiagonal(1, rho);
                    tempOp.SetDiagonal(D + 1, 1);
                    tempOp.SetDiagonal(D + 2, 1);
                    tempOp.SetDiagonal(D + 3, 1);
                    XOP.TemporalOperator = tempOp;


                    // Momentum equation
                    // ================================================================================
                    for (int d = 0; d < D; d++) {
                        var comps = XOP.EquationComponents[CodName[d]];

                        // convective part:
                        if (this.Control.StokesConvection || this.Control.Stokes) {
                            Console.WriteLine("Using Stokes Equation - no convective term.");
                            
                        } else {
                            comps.Add(new LocalLaxFriedrichsConvection(D, BcMap, d, 1.0, null));
                        }


                        // pressure part:
                        var pres = new PressureGradientLin_d(d, BcMap);
                        comps.Add(pres);
                        //Console.WriteLine("!!!Warning!!!: no pressure grad.");

                        //if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                        if (this.Control.FixedStreamwisePeriodicBC) {
                            var pressSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                            comps.Add(pressSource);
                        }


                        // viscous part:
                        if (this.Control.beta < 0.0) {
                            throw new ArithmeticException("Illegal setting in control object: 'beta' is out of range, must be non-negative.");
                        }
                        if (this.Control.Reynolds <= 0.0) {
                            throw new ArithmeticException("Illegal setting in control object: 'Reynolds' is out of range, must be strictly positive.");
                        }
                        if (this.Control.beta > 0.0) {
                            var Visc = new SipViscosity_GradU(
                                this.Control.ViscousPenaltyScaling,
                                d,
                                D,
                                BcMap,
                                ViscosityOption.ConstantViscosityDimensionless,
                                reynolds: this.Control.Reynolds / this.Control.beta);
                            comps.Add(Visc);
                        }

                        if (ConstitutiveEqs) {
                            // extra stress divergence part:
                            comps.Add(new StressDivergence_Cockburn(d, BcMap, this.Control.Reynolds, this.Control.Penalty1, this.Control.Penalty2));
                        } else {
                            Console.WriteLine("!!!Warning!!!: stress divergence deactivated");
                        }
                    }


                    // Continuum equation
                    // ===============================================================================
                    if (MomContinuitycoupling) {
                        for (int d = 0; d < D; d++) {
                            XOP.EquationComponents["div"].Add(new Divergence_DerivativeSource(d, D));
                            XOP.EquationComponents["div"].Add(new Divergence_DerivativeSource_Flux(d, BcMap));

                            //Pressure stabilization for LDG
                            //var presStab = new PressureStabilization(this.Control.PresPenalty2, this.Control.Reynolds);
                            //Console.WriteLine("PresPenalty2 = " + this.Control.PresPenalty2);
                            //XOP.EquationComponents["div"].Add(presStab);
                        }
                    } else {
                        for (int d = 0; d < D; d++) {
                            Console.WriteLine("!!!Warning!!!: Continuity Equation deactivated.");
                            XOP.EquationComponents["div"].Add(new Idsource(VariableNames.Pressure));
                        }
                    }

                    // Constitutive equations
                    // ===============================================================================

                    // Identity part
                    XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Identity(0, this.Control.giesekusfactor, this.Control.Weissenberg, this.Control.beta));
                    XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Identity(1, this.Control.giesekusfactor, this.Control.Weissenberg, this.Control.beta));
                    XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Identity(2, this.Control.giesekusfactor, this.Control.Weissenberg, this.Control.beta));

                    if (ConstitutiveEqs) {
                        Console.WriteLine($"configuring Weissenberg number: {this.Control.Weissenberg:#.##e+00}");

                        //Convective part
                        XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Convective(0, BcMap, this.Control.Weissenberg, this.Control.useFDJacobianForOperatorMatrix, this.Control.alpha));
                        XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Convective(1, BcMap, this.Control.Weissenberg, this.Control.useFDJacobianForOperatorMatrix, this.Control.alpha));
                        XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Convective(2, BcMap, this.Control.Weissenberg, this.Control.useFDJacobianForOperatorMatrix, this.Control.alpha));

                        //Objective Part
                        XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Objective(0, BcMap, this.Control.Weissenberg, this.Control.StressPenalty, this.Control.useFDJacobianForOperatorMatrix));
                        XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Objective(1, BcMap, this.Control.Weissenberg, this.Control.StressPenalty, this.Control.useFDJacobianForOperatorMatrix));
                        XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Objective(2, BcMap, this.Control.Weissenberg, this.Control.StressPenalty, this.Control.useFDJacobianForOperatorMatrix));

                        // Viscous Part
                        XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Viscosity(0, BcMap, this.Control.beta, this.Control.Penalty1));
                        XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Viscosity(1, BcMap, this.Control.beta, this.Control.Penalty1));
                        XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Viscosity(2, BcMap, this.Control.beta, this.Control.Penalty1));

                        /*
                        // artificial diffusion part
                        if (this.Control.UseArtificialDiffusion == true) {
                            XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Diffusion(this.StressXX.Basis.Degree, Grid.SpatialDimension, ((GridData)GridData).Cells.cj, VariableNames.StressXX));
                            XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Diffusion(this.StressXY.Basis.Degree, Grid.SpatialDimension, ((GridData)GridData).Cells.cj, VariableNames.StressXY));
                            XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Diffusion(this.StressYY.Basis.Degree, Grid.SpatialDimension, ((GridData)GridData).Cells.cj, VariableNames.StressYY));
                        }
                        */
                    } else {
                        Console.WriteLine("!!!Warning!!!: Constitutive Equation deactivated.");
                    }

                    // source terms (gravity/manufactured solution)
                    // ============================================

                    // Gravity Source (default should be zero!)
                    if (Control.GravitySource == true) {
                        bool test = false;

                        if (this.Control.GravityX != null && this.Control.GravityY != null) {
                            XOP.EquationComponents["momX"].Add(new ConstantSource(this.Control.GravityX));
                            XOP.EquationComponents["momY"].Add(new ConstantSource(this.Control.GravityY));
                            test = true;
                        }

                        if (this.Control.GravityXX != null && this.Control.GravityXY != null && this.Control.GravityYY != null) {
                            XOP.EquationComponents["constitutiveXX"].Add(new ConstantSource(this.Control.GravityXX));
                            XOP.EquationComponents["constitutiveXY"].Add(new ConstantSource(this.Control.GravityXY));
                            XOP.EquationComponents["constitutiveYY"].Add(new ConstantSource(this.Control.GravityYY));
                            test = true;
                        }

                        if (this.Control.GravityDiv != null) {
                            XOP.EquationComponents["div"].Add(new ConstantSource(this.Control.GravityDiv));
                            test = true;
                        }

                        if (!test) {
                            throw new ApplicationException("Gravity is true, but no values set!");
                        }
                    }


                    // artificial viscosity
                    // ====================
                    /*
                    if(Control.UseArtificialDiffusion) {

                        XOP.ParameterUpdates.Add(ArtificialViscosityUpdate);
                        
                        XOP.ParameterFactories.Add(delegate (IReadOnlyDictionary<string, DGField> DomainVarFields) {
                            return new[] { ("artificialViscosity", this.artificalViscosity as DGField) };
                        });

                    }
                    */

                    // Solver-Controlled Homotopy
                    // ==========================

                    this.currentWeissenberg = Control.Weissenberg;
                    if(Control.RaiseWeissenberg) {
                        XOP.HomotopyUpdate.Add(delegate (double HomotopyScalar) {
                            if(HomotopyScalar < 0.0)
                                throw new ArgumentOutOfRangeException();
                            if(HomotopyScalar > 1.0)
                                throw new ArgumentOutOfRangeException();

                            this.currentWeissenberg = HomotopyScalar * Control.Weissenberg;
                        });
                    } else {
                        if(Control.Weissenberg > 0)
                            Console.WriteLine($"Warning: Homotopy is turned off -- for current Weissenberg number ({Control.Weissenberg}), the solver might not converge. (Set `RaiseWeissenberg` to ture))");
                    }


                    // Build spatial operator
                    // ======================

                    if (Control.useFDJacobianForOperatorMatrix)
                        XOP.LinearizationHint = LinearizationHint.FDJacobi;
                    else
                        XOP.LinearizationHint = LinearizationHint.GetJacobiOperator;

                    XOP.Commit();

                    //JacobiOp = XOP._GetJacobiOperator(2);

                    // create timestepper
                    //===============================================================

                    m_Timestepper = new XdgTimestepping(XOP, 
                        this.CurrentSolution.Fields, this.CurrentResidual.Fields, 
                        Control.TimeSteppingScheme, 
                        this.MultigridOperatorConfig,
                        Control.LinearSolver, Control.NonLinearSolver, null, this.QueryHandler);

                    m_Timestepper.RegisterResidualLogger(this.ResLogger);
                        


                }


            }
        }
                

        
        /// <summary>
        /// Depending on settings, computes either one timestep or a steady-state solution.
        /// </summary>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {
                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;

                //if (TimestepNo[0] > 1) {
                //    this.Control.RaiseWeissenberg = false;
                //}

                base.ResLogger.TimeStep = TimestepInt;

                dt = base.GetTimestep();


                int NoIncrementTimestep;

                Console.WriteLine("Instationary solve, timestep #{0}, dt = {1} ...", TimestepNo, dt);
                var overallstart = DateTime.Now;
                bool m_SkipSolveAndEvaluateResidual = this.Control.SkipSolveAndEvaluateResidual;

                /*
                bool SolverSuccess = true;

                if (Control.RaiseWeissenberg == true) {
                    


                    currentWeissenberg = restartWeissenberg;
                    restartWeissenberg = 0.0; // make sure the restart value is used only once
                    Console.WriteLine("current Weissenberg at " + currentWeissenberg);


                    if (Control.Weissenberg != 0.0) {

                        if (Control.WeissenbergIncrement != 0.0) {
                            NoIncrementTimestep = (int)(Control.Weissenberg / Control.WeissenbergIncrement);
                        } else {
                            throw new ArgumentException("Raise Weissenberg is turned on, but WeissenbergIncrement is zero!");
                        }

                    } else {
                        throw new ArgumentException("Raise Weissenberg is turned on, but aim Weissenberg is 0.0 (Newtonian)!");
                    }

                    for (int WeIncCounter = 0; WeIncCounter <= NoIncrementTimestep; WeIncCounter++) {

                        if (Control.UseArtificialDiffusion == true) {
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // USING Weissenberg increase, USING artificial viscosity
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            artificialMaxViscosity = 1.0;

                            for (int j = 0; j < 3; j++) {




                                if (Control.UsePerssonSensor == true) {
                                    perssonsensor.Update(StressXX);
                                } else {
                                    throw new ArgumentException("artificial viscosity is turned on, but Persson sensor is turned off!");
                                }

                                m_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                                //this.ResLogger.NextTimestep(false);

                                // this evaluation must later out of this loop. now here for comparing resluts with  
                                PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, WeIncCounter));
                                SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, WeIncCounter), phystime);

                                if (Control.Bodyforces == true) {
                                    if (Log != null) {
                                        WriteLogLine(TimestepNo.MajorNumber, phystime);
                                    } else {
                                        double[] force = IBMSolverUtils.GetForces_BoundaryFitted(Velocity.Current, StressXX, StressXY, StressYY, Pressure, m_Timestepper.LsTrk, 1 / Control.Reynolds, Control.beta, "Wall_cylinder");
                                        Console.WriteLine();
                                        Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                        Console.WriteLine();
                                    }
                                }

                                artificialMaxViscosity = artificialMaxViscosity - 0.5;
                            }
                        } else {
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // USING Weissenberg increase, NO artificial viscosity
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            m_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                            //this.ResLogger.NextTimestep(false);

                            // this evaluation must later out of this loop. now here for comparing results with  
                            //PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, WeIncCounter));
                            SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, WeIncCounter), phystime);

                            if (Control.Bodyforces == true) {
                                if (Log != null) {
                                    WriteLogLine(TimestepNo.MajorNumber, phystime);
                                } else {
                                    double[] force = IBMSolverUtils.GetForces_BoundaryFitted(Velocity.Current, StressXX, StressXY, StressYY, Pressure, m_Timestepper.LsTrk, 1 / Control.Reynolds, Control.beta, "Wall_cylinder");
                                    Console.WriteLine();
                                    Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                    Console.WriteLine();
                                }
                            }
                        }

                        ChangeMesh = Control.AdaptiveMeshRefinement;
                        while (ChangeMesh == true) {
                            this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                            perssonsensor.Update(StressXX);
                            PlotCurrentState(phystime, TimestepNo);
                            SaveToDatabase(TimestepNo, phystime);
                        }

                        if (currentWeissenberg < Control.Weissenberg) {
                            currentWeissenberg = currentWeissenberg + Control.WeissenbergIncrement;
                        } else {
                            WeIncCounter = int.MaxValue; // breaks the Weissenberg-increase loop.
                        }

                        currentWeissenberg = Math.Min(currentWeissenberg, Control.Weissenberg); // prevents any overshoot
                        Console.WriteLine("Raise Weissenberg number to " + currentWeissenberg);
                    }
                } else {

                    currentWeissenberg = Control.Weissenberg;

                    if (Control.UseArtificialDiffusion == true) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // USING artificial viscosity, but NO Weissenberg increase
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++


                        artificialMaxViscosity = 1.0;

                        for (int j = 0; j < 3; j++) {

                            if (Control.UsePerssonSensor == true) {
                                perssonsensor.Update(StressXX);
                            } else {
                                throw new ArgumentException("artificial viscosity is turned on, but Persson sensor is turned off!");
                            }

                            m_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                            // this evaluation must later out of this loop. now here for comparing resluts with  
                            //PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                            //SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                            if (Control.Bodyforces == true) {
                                if (Log != null) {
                                    WriteLogLine(TimestepNo.MajorNumber, phystime);
                                } else {
                                    double[] force = IBMSolverUtils.GetForces_BoundaryFitted(Velocity.Current, StressXX, StressXY, StressYY, Pressure, m_Timestepper.LsTrk, 1 / Control.Reynolds, Control.beta, "Wall_cylinder");
                                    Console.WriteLine();
                                    Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                    Console.WriteLine();
                                }
                            }

                            artificialMaxViscosity = artificialMaxViscosity - 0.5;
                        }

                        ChangeMesh = Control.AdaptiveMeshRefinement;
                        while (ChangeMesh == true) {
                            this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                        }
                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Using only timestepper, NO ADDITIONAL LOOP
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        SolverSuccess = m_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);
                                                

                        if (Control.Bodyforces == true) {
                            if (Log != null) {
                                WriteLogLine(TimestepNo.MajorNumber, phystime);
                            } else {
                                double[] force = IBMSolverUtils.GetForces_BoundaryFitted(Velocity.Current, StressXX, StressXY, StressYY, Pressure, m_Timestepper.LsTrk, 1 / Control.Reynolds, Control.beta, "Wall_cylinder");
                                Console.WriteLine();
                                Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                Console.WriteLine();
                            }
                        }
                    }
                }
                */

                bool SolverSuccess = m_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);


                if(Control.Bodyforces == true) {
                    if(Log != null) {
                        WriteLogLine(TimestepNo.MajorNumber, phystime);
                    } else {
                        double[] force = IBMSolverUtils.GetForces_BoundaryFitted(Velocity.Current, StressXX, StressXY, StressYY, Pressure, m_Timestepper.LsTrk, 1 / Control.Reynolds, Control.beta, "Wall_cylinder");
                        Console.WriteLine();
                        Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                        Console.WriteLine();
                    }
                }


                var overallstop = DateTime.Now;
                var overallduration = overallstop - overallstart;

                Console.WriteLine("Duration of this timestep: " + overallduration);

                if(!SolverSuccess)
                    base.CurrentSessionInfo.AddTag(SessionInfo.SOLVER_ERROR);


                this.ResLogger.NextTimestep(false);

                return dt;


            }
        }


        /// <summary>
        /// implements a <see cref="DelPartialParameterUpdate"/>
        /// </summary>
        /// <param name="DomainVarFields"></param>
        /// <param name="ParameterVarFields"></param>
        void ArtificialViscosityUpdate(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) { 
        //void ParameterUpdate(IEnumerable<DGField> CurrentState, IEnumerable<DGField> ParameterVar) {
            /*
            int D = this.GridData.SpatialDimension;

            var U0 = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());
            var Stress0 = new VectorField<SinglePhaseField>(CurrentState.Skip(D + 1).Take(3).Select(F => (SinglePhaseField)F).ToArray());

            if (this.U0MeanRequired) {

                SinglePhaseField[] __U0mean = ParameterVar.Skip(D).Take(D).Select(f => f as SinglePhaseField).ToArray();
                VectorField<SinglePhaseField> U0mean = new VectorField<SinglePhaseField>(__U0mean);

                U0mean.Clear();
                ComputeAverageU(U0, U0mean);

                SinglePhaseField[] __U0 = ParameterVar.Take(D).Select(f => f as SinglePhaseField).ToArray();
                Debug.Assert(ArrayTools.AreEqual(__U0, U0.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));
            } else {
                Debug.Assert(ParameterVar.Take(2 * D).Where(f => f != null).Count() == 0);
            }

            if (this.Control.SetParamsAnalyticalSol == false) {
                SinglePhaseField[] __VelocityXGradient = ParameterVar.Skip(2 * D).Take(D).Select(f => f as SinglePhaseField).ToArray();
                SinglePhaseField[] __VelocityYGradient = ParameterVar.Skip(3 * D).Take(D).Select(f => f as SinglePhaseField).ToArray();
                Debug.Assert(ArrayTools.AreEqual(__VelocityXGradient, VelocityXGradient.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));
                Debug.Assert(ArrayTools.AreEqual(__VelocityYGradient, VelocityYGradient.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));

                VelocityXGradient.Clear();
                VelocityXGradient.GradientByFlux(1.0, U0[0]);
                VelocityYGradient.Clear();
                VelocityYGradient.GradientByFlux(1.0, U0[1]);
            }
            */
            /*
            if (this.Control.UseArtificialDiffusion == true) {

                //SinglePhaseField __ArtificialViscosity = ParameterVar.Skip(5 * D + 1).Take(1).Select(f => f as SinglePhaseField).ToArray()[0];
                SinglePhaseField __ArtificialViscosity = ParameterVarFields["artificialViscosity"] as SinglePhaseField;
                if (!object.ReferenceEquals(this.artificalViscosity, __ArtificialViscosity))
                    throw new ApplicationException();

                ArtificialViscosity.ProjectArtificalViscosityToDGField(__ArtificialViscosity, perssonsensor, this.Control.SensorLimit, artificialMaxViscosity);
            }
            */
        }

        /// <summary>
        /// Only for testing / NUnit:
        /// checks whether the finite difference approximation of the Jacobian of <see cref="XOP"/>
        /// and the Jacobian operator (<see cref="IDifferentialOperator.GetJacobiOperator"/>)
        /// provide approximately the same matrix and affine vector.
        /// </summary>
        internal void CheckJacobian() {
            // Parameters
            DGField[] Params = null;
            /*
            if (this.Control.UseArtificialDiffusion) {
                Params = new[] { artificalViscosity };
            } else {
                Params = null;
            }
            */

            // initialize linearization point with random numbers
            var CurrentState = this.CurrentSolution.Fields.Select(f => f.CloneAs()).ToArray();
            var CurrentVector = new CoordinateVector(CurrentState);
            Random r = new Random(0); // seed of 0 guarantees the same random numbers on every run
            int L = CurrentVector.Count;
            for (int i = 0; i < L; i++) {
                CurrentVector[i] = r.NextDouble();
            }

            var domMap = CurrentVector.Mapping;
            var codMap = domMap;
            Assert.IsTrue(codMap.EqualsPartition(this.CurrentResidual.Mapping));

            // Finite Difference Linearization
            var FDbuilder = XOP.GetFDJacobianBuilder_(domMap, null, codMap, null);
            var JacobianFD = new BlockMsrMatrix(codMap, domMap);
            var AffineFD = new double[JacobianFD.NoOfRows];
            FDbuilder.ComputeMatrix(JacobianFD, AffineFD);

            // Jacobian Operator
            var JacobiOp = this.XOP._GetJacobiOperator(2);
            //var JacParams = JacobiOp.ParameterUpdate;
            var TmpParams = JacobiOp.InvokeParameterFactory(domMap.Fields);
            var map = new CoordinateMapping(CurrentState);
            var JacBuilder = JacobiOp.GetMatrixBuilder(map, TmpParams, map);
            //this.ParameterUpdate(CurrentState, TmpParams);
            //JacParams.PerformUpdate(CurrentState, TmpParams);
            JacobiOp.InvokeParameterUpdate(0.0, CurrentState, TmpParams);
            
            var JacobiDX = new BlockMsrMatrix(map);
            var AffineDX = new double[map.LocalLength];
            JacBuilder.ComputeMatrix(JacobiDX, AffineDX);

            // Comparison
            Console.WriteLine("Comparison of finite difference and direct Jacobian matrix:");
            var ErrMtx = JacobianFD.CloneAs();
            ErrMtx.Acc(-1.0, JacobiDX);
            double InfNorm_ErrMtx = ErrMtx.InfNorm();
            Console.WriteLine("  Jacobian Matrix Delta Norm: " + InfNorm_ErrMtx);

            var ErrAff = AffineFD.CloneAs();
            ErrAff.AccV(-1.0, AffineDX);
            double InfNorm_ErrAff = ErrAff.MPI_L2Norm();
            Console.WriteLine("  Affine Vector Delta Norm: " + InfNorm_ErrAff);

            // Error Threshold checks
            double DenomM = (JacobianFD.InfNorm(), JacobiDX.InfNorm()).Max();
            Assert.Less(InfNorm_ErrMtx / DenomM, 0.01, "Mismatch in between finite difference Jacobi matrix and direct Jacobi matrix");

            double DenomA = (CurrentVector.MPI_L2Norm(), AffineFD.MPI_L2Norm(), AffineDX.MPI_L2Norm()).Max();
            Assert.Less(InfNorm_ErrAff / DenomA, 0.01, "Mismatch in Affine Vector between finite difference Jacobi and direct Jacobi");
        }

        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        public MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int pVel = this.Velocity.Current[0].Basis.Degree;
                int pPrs = this.Pressure.Basis.Degree;
                int pStr = this.StressXX.Basis.Degree;
                int D = this.GridData.SpatialDimension;
                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration entry will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    int pVelLv = Math.Max(1, pVel - iLevel);
                    int pPreLv = Math.Max(0, pPrs - iLevel);
                    int pStrLv = Math.Max(1, pStr - iLevel);

                    /*
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[2];
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { pVelLv, pVelLv, pStrLv, pStrLv, pStrLv },
                            //mode = this.Control.VelocityBlockPrecondMode,
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            VarIndex = new int[] { 0, 1, 3, 4, 5 }
                        };
                    configs[iLevel][1] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                        mode = this.Control.PressureBlockPrecondMode,
                        VarIndex = new int[] { D }
                    };
                    //*/



                    /*
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[1];
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        mode = MultigridOperator.Mode.LeftInverse_DiagBlock,
                        VarIndex = new int[] { 0, 1, 2, 3, 4, 5 },
                        DegreeS = new int[] { pVel, pVel, pPrs, pStr, pStr, pStr }
                    };
                    //*/

                    
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 4];
                    
                    // configurations for velocity
                    for (int d = 0; d < D; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pVel) },
                            //DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                            //mode = this.Control.VelocityBlockPrecondMode,
                            //mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(0, pPrs) },
                        //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                        //mode = MultigridOperator.Mode.Eye,
                        mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                        VarIndex = new int[] { D }
                    };

                    
                    // configurations for stresses
                    for (int d = 3; d < 6; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pStr) },//DegreeS = new int[] { Math.Max(1, pStr - iLevel) },
                            //mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                            VarIndex = new int[] { d }
                        };
                    }
                    //*/
                }

                return configs;
            }
        }

        /// <summary>
        /// Plotting the current state
        /// </summary>
        protected override void PlotCurrentState(double physTime, Foundation.IO.TimestepNumber timestepNo, int superSampling = 0) {
            // Standard
            DGField[] myFields = ArrayTools.Cat<DGField>(Velocity.Current, ResidualMomentum, ResidualConti, Pressure, StressXX, StressXY, StressYY, ResidualStressXX, ResidualStressXY, ResidualStressYY); //, VelocityXGradient, VelocityYGradient, Gravity

            //Add sensor field only if Persson sensor exists
            if (perssonsensor != null) {
                myFields = ArrayTools.Cat<DGField>(myFields, perssonsensor.GetField());
            }

            //Add field only if artificial viscosity is turned on
            if (artificalViscosity != null) {
                myFields = ArrayTools.Cat<DGField>(myFields, artificalViscosity);
            }

            Tecplot.PlotFields(myFields, "Rheology-" + timestepNo.ToString(), physTime, superSampling);
        }

        /*
        /// <summary>
        /// Plotting the in interation callback
        /// </summary>
        protected void PlotOnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            // Standard
            DGField[] myFields = ArrayTools.Cat<DGField>(Velocity.Current, ResidualMomentum, ResidualConti, Pressure, StressXX, StressXY, StressYY, LevSet, ResidualStressXX, ResidualStressXY, ResidualStressYY); //, VelocityXGradient, VelocityYGradient, Gravity,

            //Add sensor field only if Persson sensor exists
            if (perssonsensor != null) {
                myFields = ArrayTools.Cat<DGField>(myFields, perssonsensor.GetField());
            }

            //Add field only if artificial viscosity is turned on
            if (artificalViscosity != null) {
                myFields = ArrayTools.Cat<DGField>(myFields, artificalViscosity);
            }

            Tecplot.PlotFields(myFields, "Rheology-" + iterIndex.ToString(), 0.0, 2);
        }
        */
        /// <summary>
        /// Initializing the DG fields
        /// </summary>
        protected override void SetInitial(double t) {
            int D = GridData.SpatialDimension;
            if (D != 2)
                throw new NotImplementedException("currently only support for 2 dimensions.");

            base.SetInitial(t);
            CreateEquationsAndSolvers(null);
            this.LsTrk = m_Timestepper.LsTrk;

            VelocityXGradient = new VectorField<SinglePhaseField>(D, Velocity.Current[0].Basis, "VelocityX_Gradient", SinglePhaseField.Factory);
            VelocityYGradient = new VectorField<SinglePhaseField>(D, Velocity.Current[1].Basis, "VelocityY_Gradient", SinglePhaseField.Factory);

            if (this.Control.SetParamsAnalyticalSol == true) {
                U = new SinglePhaseField(new Basis(this.GridData, Velocity.Current[0].Basis.Degree), "UAnalytical");
                V = new SinglePhaseField(new Basis(this.GridData, Velocity.Current[0].Basis.Degree), "VAnalytical");
                U.ProjectField(this.Control.VelFunctionU);
                V.ProjectField(this.Control.VelFunctionV);

                VelocityXGradient.Clear();
                VelocityXGradient.Gradient(1.0, U);
                VelocityYGradient.Clear();
                VelocityYGradient.Gradient(1.0, V);
            }

            Console.WriteLine("Total number of cells:    {0}", Grid.NumberOfCells);
            Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Mapping.TotalLength);

            if (this.CurrentSessionInfo.ID != Guid.Empty && base.MPIRank == 0) {
                InitLogFile(this.CurrentSessionInfo.ID);
            }

        }


        /// <summary>
        /// performs restart
        /// </summary>
        /// <param name="Time">
        /// on exit, the physical time associated with the field state
        /// </param>
        /// <param name="TimestepNo">
        /// on exit, the physical time associated with the field state
        /// </param>
        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);

            this.LsTrk.UpdateTracker(Time);
        }

        /// <summary>
        /// overriding the method to implement any user-specific tasks which
        /// should be carried out after a restart file has been loaded (e.g.,
        /// setting the correct time for a time-stepper)
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            VelocityXGradient = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, Velocity.Current[0].Basis, "VelocityX_Gradient", SinglePhaseField.Factory);
            VelocityYGradient = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, Velocity.Current[1].Basis, "VelocityY_Gradient", SinglePhaseField.Factory);
        }

       
        /// <summary>
        /// Computes the L2 Error of all Fields compared to exact solution specified in the control file
        /// </summary>
        protected void ComputeL2Error() {
            if (this.Control.ExSol_Velocity == null && this.Control.ExSol_Pressure == null && this.Control.ExSol_Stress == null) {
                // nothing to do
                return;
            }


            int D = this.GridData.SpatialDimension;

            int order = Velocity.Current[0].Basis.Degree * 2;

            // Velocity error
            // ===============================================
            if (this.Control.ExSol_Velocity != null) {
                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                for (int d = 0; d < D; d++) {
                    L2Error[d] = this.Velocity.Current[d].L2Error(this.Control.ExSol_Velocity[d].Vectorize(0.0), order);
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Console.WriteLine("L2err " + VariableNames.Velocity_d(d) + " is " + L2Error[d]);
                }
            }


            // pressure error
            // =============================================================
            if (this.Control.ExSol_Pressure != null) {

                double L2Error = 0;

                L2Error = this.Pressure.L2Error(this.Control.ExSol_Pressure.Vectorize(0.0), order - 1);
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
                Console.WriteLine("L2err " + VariableNames.Pressure + " is " + L2Error);
            }

            // Stress error
            // =============================================================
            if (this.Control.ExSol_Stress != null) {
                double[] L2Error = new double[3];

                L2Error[0] = this.StressXX.L2Error(this.Control.ExSol_Stress[0].Vectorize(0.0), order);
                L2Error[1] = this.StressXY.L2Error(this.Control.ExSol_Stress[1].Vectorize(0.0), order);
                L2Error[2] = this.StressYY.L2Error(this.Control.ExSol_Stress[2].Vectorize(0.0), order);

                base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressXX, L2Error[0], true);
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressXY, L2Error[1], true);
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressYY, L2Error[2], true);

                Console.WriteLine("L2err " + VariableNames.StressXX + " is " + L2Error[0]);
                Console.WriteLine("L2err " + VariableNames.StressXY + " is " + L2Error[1]);
                Console.WriteLine("L2err " + VariableNames.StressYY + " is " + L2Error[2]);
            }
        }



        /// <summary>
        /// Integration degree of HMF used throughout the application: this should ensure that
        /// only one HMF rule is created.
        /// </summary>
        public int HMForder {
            get {
                int VelDeg = this.Velocity.Current.Max(field => field.Basis.Degree);
                int Order = (VelDeg * (!this.Control.Stokes ? 3 : 2));
                Order += 2; // safety factor
                return Order;
            }
        }

        //ADAPTIVE MESH REFINEMENT
        //======================================================================

        /// <summary>
        /// refinement indicator
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel) {
            /*
            if (this.Control.UsePerssonSensor) {

                double maxVal = this.perssonsensor.GetValue(j);

                double[] coord = this.GridData.iGeomCells.GetCenter(j);

                //bound for perssonsensor should be around 1e-7 - 1e-8 that there is refinement behind the cylinder!
                double upperbound = this.Control.SensorLimit;
                double lowerbound = upperbound * 0.001;

                int DesiredLevel_j = CurrentLevel;

                if (maxVal != 0.0) {
                    if (maxVal > upperbound && DesiredLevel_j < this.Control.RefinementLevel) {

                        DesiredLevel_j = DesiredLevel_j + 1;

                    } else if (maxVal < lowerbound && DesiredLevel_j > 0) {
                        DesiredLevel_j = DesiredLevel_j - 1;
                    }
                } else {
                    if (Math.Abs(coord[0] - 10) < 2 && DesiredLevel_j < 2) // this.Control.RefinementLevel)
                        DesiredLevel_j = DesiredLevel_j + 1;
                }

                return DesiredLevel_j;

            } else {

                double celllength = ((GridData)GridData).Cells.cj[j];
                double maxVal = this.StressXX.GetMeanValue(j) / celllength;  // this.perssonsensor.GetValue(j);

                //bound for perssonsensor should be around 1e-7 - 1e-8 that there is refinement behind the cylinder!
                double upperbound = this.Control.SensorLimit / celllength;
                double lowerbound = -1 * this.Control.SensorLimit / celllength;

                int DesiredLevel_j = CurrentLevel;

                if (maxVal > upperbound && DesiredLevel_j < this.Control.RefinementLevel) {

                    DesiredLevel_j = DesiredLevel_j + 1;

                } else if (maxVal < lowerbound && DesiredLevel_j > 0) {
                    DesiredLevel_j = DesiredLevel_j - 1;
                }

                return DesiredLevel_j;
            }
            */
            return 0;
        }

        /// <summary>
        /// Adaptation of the current mesh.
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {

            if (this.Control.AdaptiveMeshRefinement) {
                GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData, null);
                bool AnyChange = gridRefinementController.ComputeGridChange(LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                ChangeMesh = AnyChange;
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
                long oldJ = this.GridData.CellPartitioning.TotalLength;

                // Update Grid
                // ===========

                if (AnyChange) {

                    //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 1 }), 2);
                    Console.WriteLine();
                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                    newGrid = ((GridData)(this.GridData)).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                    //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 2 }), 2);#

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }
            } else {

                newGrid = null;
                old2NewGrid = null;
            }
        }

        /// <summary>
        /// Appends the <see cref="currentWeissenberg"/> number to the timestep
        /// </summary>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {
            var Rtsi = new RheologyTimestepInfo(t, CurrentSessionInfo, timestepno, IOFields, currentWeissenberg);
            return Rtsi;
        }

        /// <summary>
        /// sets Weissenberg number from timestep-info 
        /// </summary>
        protected override void OnRestartTimestepInfo(TimestepInfo tsi) {
            if (this.Control.RaiseWeissenberg) {

                var Rtsi = tsi as RheologyTimestepInfo;
                if (Rtsi != null) {
                    Console.Write("Restoring Weissenberg number form database...  ");
                    Console.Write($" Weissenberg = {Rtsi.currentWeissenbergNumber}");
                    this.restartWeissenberg = Rtsi.currentWeissenbergNumber;
                    Console.WriteLine();
                } else {
                    Console.WriteLine($"No Weissenberg number contained in time-step; starting with pre-set.");
                }
            }
        }

        /// <summary>
        /// automatized analysis of condition number 
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis() {

            int[] varGroup_convDiff = new int[] { 0, 1 };
            int[] varGroup_Stokes = new int[] { 0, 1, 2 };
            int[] varGroup_Constitutive = new int[] { 3, 4, 5 };
            int[] varGroup_all = new int[] { 0, 1, 2, 3, 4, 5 };

            var res = m_Timestepper.TimesteppingBase.OperatorAnalysis(new[] {varGroup_convDiff, varGroup_Stokes, varGroup_Constitutive, varGroup_all });

            // filter only those results that we want;
            // this is a DG app, but it uses the LevelSetTracker; therefore, we want to filter analysis results for cut cells and only return uncut cells resutls
            var ret = new Dictionary<string, double>();
            foreach(var kv in res) {
                if(kv.Key.ToLowerInvariant().Contains("innercut") || kv.Key.ToLowerInvariant().Contains("bndycut")) {
                    // ignore
                } else {
                    ret.Add(kv.Key, kv.Value);
                }
            }

            return ret;
        }


        #region logging

        TextWriter Log;
        string header;

        public void InitLogFile(Guid sessionID) {

            if (this.Control.Bodyforces) {
                Log = base.DatabaseDriver.FsDriver.GetNewLog("BodyForces", sessionID);
                header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#timestep", "#time","Wi", "ForceX", "ForceY");
            }
        }

        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {

            double[] force = IBMSolverUtils.GetForces_BoundaryFitted(Velocity.Current, StressXX, StressXY, StressYY, Pressure, m_Timestepper.LsTrk, 1 / Control.Reynolds, Control.beta, "Wall_cylinder");
            string logline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, currentWeissenberg, force[0], force[1]);
            Log.WriteLine(logline);
            Log.Flush();
        }

        #endregion
    }  

    class Idsource : LinearSource {
        public Idsource(string _var) {
            m_var = _var;
        }

        string m_var;

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { m_var };
            }
        }

        protected override double Source(double[] x, double[] parameters, double[] U) {
            return U[0];
        }
    }

}


