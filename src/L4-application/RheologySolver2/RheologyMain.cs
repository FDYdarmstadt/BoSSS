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
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.Gnuplot;

using MPI.Wrappers;
using NUnit.Framework;
using BoSSS.Foundation.SpecFEM;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Solver for calculation with viscoelastic extra stresses using the Oldroyd B model or the upper convected Maxwell model (UCM).
    /// </summary>
<<<<<<< HEAD

=======
>>>>>>> root-public
    public class Rheology : BoSSS.Solution.Application<RheologyControl> {
        static void Main(string[] args) {
            
            Rheology._Main(args, false,
                delegate () {
                    var app = new Rheology();
                    return app;
                });
        }

<<<<<<< HEAD
        // Attributes for fields (Names), initialization of DG fields
        //==============================================================

        // Velocity
=======
        #region instantiation

        // Attributes for fields (Names), initialization of DG fields
        //==============================================================

        /// <summary>
        /// Velocity domain
        /// </summary>
>>>>>>> root-public
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            null,
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorFieldHistory<SinglePhaseField> Velocity;

<<<<<<< HEAD
        // Residuum in momentum equation
=======
        /// <summary>
        /// Velocities codomain: Residuum in momentum equation
        /// </summary>
>>>>>>> root-public
        [InstantiateFromControlFile(new string[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> ResidualMomentum;

<<<<<<< HEAD
        //Pressure
        [InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField Pressure;

        // Residuum in continuity equation
        [InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualConti;

        //Extra stresses domain(2D)
        [InstantiateFromControlFile("StressXX", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXX;

        [InstantiateFromControlFile("StressXY", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXY;

        [InstantiateFromControlFile("StressYY", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressYY;

        // Extra Stresses Codomain (2D)
        [InstantiateFromControlFile("ResidualStressXX", "StressXX", IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressXX;

        [InstantiateFromControlFile("ResidualStressXY", "StressXY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressXY;

        [InstantiateFromControlFile("ResidualStressYY", "StressYY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressYY;

        //Extra Stresses Parameter
        [InstantiateFromControlFile("StressXXP", "StressXX", IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXXP;

        [InstantiateFromControlFile("StressXYP", "StressXY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXYP;

        [InstantiateFromControlFile("StressYYP", "StressYY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressYYP;

        ////Velocity Gradient domain(2D)
        //[InstantiateFromControlFile("VelocityXGradientX", null, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField VelocityXGradientX;

        //[InstantiateFromControlFile("VelocityXGradientY", null, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField VelocityXGradientY;

        //[InstantiateFromControlFile("VelocityYGradientX", null, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField VelocityYGradientX;

        //[InstantiateFromControlFile("VelocityYGradientY", null, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField VelocityYGradientY;

        //// Velocity Gradient Codomain (2D)
        //[InstantiateFromControlFile("ResidualGradXX", "VelocityXGradientX", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ResidualGradXX;

        //[InstantiateFromControlFile("ResidualGradXY", "VelocityXGradientY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ResidualGradXY;

        //[InstantiateFromControlFile("ResidualGradYX", "VelocityYGradientX", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ResidualGradYX;

        //[InstantiateFromControlFile("ResidualGradYY", "VelocityYGradientY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ResidualGradYY;

        // Extra Source (e.g. Gravity)
=======
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

        /// <summary>
        /// Extra stress domain (2D): StressXX
        /// </summary>
        [InstantiateFromControlFile("StressXX", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXX;

        /// <summary>
        /// Extra stress domain (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("StressXY", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXY;

        /// <summary>
        /// Extra stress domain (2D): StressYY
        /// </summary>
        [InstantiateFromControlFile("StressYY", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressYY;

        /// <summary>
        /// Extra stress codomain (2D): StressXX
        /// </summary>
        [InstantiateFromControlFile("ResidualStressXX", "StressXX", IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressXX;

        /// <summary>
        /// Extra stresses codomain (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("ResidualStressXY", "StressXY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressXY;

        /// <summary>
        /// Extra stresses codomain (2D): StressYY
        /// </summary>
        [InstantiateFromControlFile("ResidualStressYY", "StressYY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualStressYY;

        /// <summary>
        /// Extra stresses parameter (2D): StressXX
        /// </summary>
        [InstantiateFromControlFile("StressXXP", "StressXX", IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXXP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("StressXYP", "StressXY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressXYP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        [InstantiateFromControlFile("StressYYP", "StressYY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField StressYYP;

        /// <summary>
        /// Extra source (e.g. gravity)
        /// </summary>
>>>>>>> root-public
        [InstantiateFromControlFile(new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
                    new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                    true, true,
                    IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Gravity;

<<<<<<< HEAD
        // Gravity source constitutive
        [InstantiateFromControlFile("GravityXX", "StressXX", IOListOption.ControlFileDetermined)]
        public SinglePhaseField GravityXX;

        [InstantiateFromControlFile("GravityXY", "StressXY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField GravityXY;

        [InstantiateFromControlFile("GravityYY", "StressYY", IOListOption.ControlFileDetermined)]
        public SinglePhaseField GravityYY;

        //Gravity source for divergence of u
        [InstantiateFromControlFile("GravityDiv", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        public SinglePhaseField GravityDiv;
=======
        //// Gravity source constitutive
        //[InstantiateFromControlFile("GravityXX", "StressXX", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityXX;

        //[InstantiateFromControlFile("GravityXY", "StressXY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityXY;

        //[InstantiateFromControlFile("GravityYY", "StressYY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityYY;

        ////Gravity source for divergence of u
        //[InstantiateFromControlFile("GravityDiv", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityDiv;
>>>>>>> root-public

        // Parameters: Velocity Gradient
        VectorField<SinglePhaseField> VelocityXGradient;
        VectorField<SinglePhaseField> VelocityYGradient;

        //Parameters: external analytical velocity
        SinglePhaseField U;
        SinglePhaseField V;

<<<<<<< HEAD
        // LEVEL-SET - Dummy, not needed for non-Level-set calculations
        //_______________________________________________________________________________________________

        // Level-Set tracker
=======
        // LEVEL-SET - not needed for non-Level-set calculations
        //_______________________________________________________________________________________________

        /// <summary>
        /// Levelset tracker
        /// </summary>
>>>>>>> root-public
        [LevelSetTracker("-:A +:B", 1)]
        public LevelSetTracker LevSetTrk;

        /// <summary>
        /// The  continuous level set field which defines the XDG space; 
<<<<<<< HEAD
        /// it is obtained from the projection of the discontinuous <see cref="DGLevSet"/> onto the 
=======
        /// it is obtained from the projection of the discontinuous DG Level set onto the 
>>>>>>> root-public
        /// continuous element space.
        /// </summary>
        [InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        public LevelSet LevSet;

<<<<<<< HEAD
=======
        /// <summary>
        /// Species which represents the flow domain.
        /// </summary>
>>>>>>> root-public
        protected SpeciesId[] FluidSpecies {
            get {
                return new SpeciesId[] { LsTrk.GetSpeciesId("A") }; // wir rechnen nur species A
            }
        }

<<<<<<< HEAD
=======
        /// <summary>
        /// Actual type of cut cell quadrature to use; If no XDG is used, resp. no cut cells are present,
        /// this setting has no effect.
        /// </summary>
>>>>>>> root-public
        protected XQuadFactoryHelper.MomentFittingVariants momentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.OneStepGauss;
        //_______________________________________________________________________________________________


<<<<<<< HEAD

        IncompressibleBoundaryCondMap BcMap;
        int D; // Spatial Dimension
        int PressureRefCellIndex; // Index of cell with reference Pressure
        double currentWeissenberg;
        bool ChangeMesh = true;


=======
        // Some initialisation of variables
        //============================================
        IncompressibleBoundaryCondMap BcMap;
        int D; // Spatial Dimension
        double currentWeissenberg;
        bool ChangeMesh = true;
        SpatialOperator XOP;
        CoordinateVector m_CurrentSolution = null;
        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// initialisation of BDF Timestepper
        /// </summary>
        protected XdgBDFTimestepping m_BDF_Timestepper;

        /// <summary>
        /// initialisation of spatial operator matrix analysis
        /// </summary>
        protected SpatialOperatorAnalysis SpatialOperatorAnalysis;


        // Persson sensor and artificial viscosity
        //=============================================
        /// <summary>
        /// initialisation of Persson sensor
        /// </summary>
        protected PerssonSensor perssonsensor;

        /// <summary>
        /// initialisation of artificial viscosity
        /// </summary>
        protected SinglePhaseField artificalViscosity;

        /// <summary>
        /// initialisation of max value of artificial viscosity
        /// </summary>
        protected double artificialMaxViscosity;


        // Settings for calculation
        //===============================================
        /// <summary>
        /// Set true if Navier Stokes is solved, then the mean velocities as parameters for calculation of convective terms are needed
        /// </summary>
>>>>>>> root-public
        protected bool U0MeanRequired {
            get {
                return (!this.Control.Stokes);
            }
        }

<<<<<<< HEAD
        SpatialOperator XOP;

        CoordinateVector m_CurrentSolution = null;

        public CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(
                                                             this.Velocity.Current,
                                                             this.Pressure,
                        this.StressXX, this.StressXY, this.StressYY)); //  this.VelocityXGradientX, this.VelocityXGradientY, this.VelocityYGradientX, this.VelocityYGradientY
=======
        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double rho = 1; // this.Control.PhysicalParameters.rho_A;

                int D = this.GridData.SpatialDimension;

                double[] _rho = new double[D + 4];
                _rho.SetAll(rho);
                //No MassMatrix for the pressure
                _rho[D] = 0;

                _rho[D + 1] = 1;
                _rho[D + 2] = 1;
                _rho[D + 3] = 1;
                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho);

                return R;
            }
        }

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.Velocity.Current, this.Pressure, this.StressXX, this.StressXY, this.StressYY));
>>>>>>> root-public
                }
                return m_CurrentSolution;
            }
        }

<<<<<<< HEAD
        CoordinateVector m_CurrentResidual = null;

        public CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(
                        this.ResidualMomentum,
                        this.ResidualConti,
                        this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY)); // this.ResidualGradXX, this. ResidualGradXY,this.ResidualGradYX, this.ResidualGradYY
=======
        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(this.ResidualMomentum, this.ResidualConti, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY));
>>>>>>> root-public
                }
                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// DG Field instantiation.
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();
            base.LsTrk = this.LevSetTrk;
            if(Control.CutCellQuadratureType != base.LsTrk.CutCellQuadratureType)
                throw new ApplicationException();
            if (Control.UsePerssonSensor == true) {
                perssonsensor = new PerssonSensor(StressXX);
                this.IOFields.Add(perssonsensor.GetField());
            }
            if (Control.UseArtificialDiffusion == true) {
                artificalViscosity = new SinglePhaseField(new Basis(GridData, 1), "artificalViscosity");
                this.IOFields.Add(artificalViscosity);

            }
        }

<<<<<<< HEAD
        protected XdgBDFTimestepping m_BDF_Timestepper;
        protected PerssonSensor perssonsensor;
        protected SinglePhaseField artificalViscosity;
        protected double artificialMaxViscosity;


=======
        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
>>>>>>> root-public
        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
        }

<<<<<<< HEAD
=======
        #endregion
>>>>>>> root-public

        /// <summary>
        /// Initialize Calculation, Create Equations
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            if (XOP != null && L == null && Control.Weissenberg == 0.0)
                return;

            if (m_BDF_Timestepper != null) {
                if (L != null) {

                    m_BDF_Timestepper.DataRestoreAfterBalancing(L,
                        ArrayTools.Cat<DGField>(Velocity.Current, Pressure, StressXX, StressXY, StressYY),
                        ArrayTools.Cat<DGField>(ResidualMomentum, ResidualConti, ResidualStressXX, ResidualStressXY, ResidualStressYY),
                        this.LsTrk, this.MultigridSequence);

                    m_CurrentSolution = null;
                    m_CurrentResidual = null;

                }
            } else {

                using (new FuncTrace()) {

                    D = this.GridData.SpatialDimension;
                    BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Viscoelastic);

<<<<<<< HEAD
                    string[] CodName = new string[] { "momX", "momY", "div", "constitutiveXX", "constitutiveXY", "constitutiveYY" }; //"velocitygradXX", "velocitygradXY", "velocitygradYX", "velocitygradYY"

                    string[] Params = ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D),
                        VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector(),
                        VariableNames.StressXXP, VariableNames.StressXYP, VariableNames.StressYYP, "artificialViscosity");

                    string[] DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure,
                        VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY); //VariableNames.VelocityXGradientX, VariableNames.VelocityXGradientY, VariableNames.VelocityYGradientX, VariableNames.VelocityYGradientY
=======
                    string[] CodName = new string[] { "momX", "momY", "div", "constitutiveXX", "constitutiveXY", "constitutiveYY" };

                    string[] Params = ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector(), VariableNames.StressXXP, VariableNames.StressXYP, VariableNames.StressYYP, "artificialViscosity");

                    string[] DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure, VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY);
>>>>>>> root-public

                    XOP = new SpatialOperator(DomName, Params, CodName, QuadOrderFunc.NonLinear(2));

                    // Momentum equation
                    //================================================================================
                    for (int d = 0; d < D; d++) {
                        var comps = XOP.EquationComponents[CodName[d]];

                        // convective part:
                        if (!this.Control.Stokes) {
                            comps.Add(new LinearizedConvection(D, BcMap, d));
                        }

                        // pressure part:
                        var pres = new PressureGradientLin_d(d, BcMap);
                        comps.Add(pres);

<<<<<<< HEAD

                        //if outflow boundary conditions are applied we need a pressure ref point
                        if (!this.BcMap.DirichletPressureBoundary) {
                            PressureRefCellIndex = SolverUtils.GetIndexOfPressureReferencePoint(new double[D], this.CurrentSolution.Mapping, D);
                        }

=======
>>>>>>> root-public
                        //if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                        if (this.Control.FixedStreamwisePeriodicBC) {
                            var pressSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                            comps.Add(pressSource);
                        }

                        // viscous part:
                        Type GridType = GridData.iGeomCells.RefElements[0].GetType();
                        double PenaltyBase;
                        int DegreeVelocity = this.Velocity.Current.Max(DGF => DGF.Basis.Degree);
                        if ((GridType == typeof(Triangle)) || (GridType == typeof(Tetra))) {
                            PenaltyBase = (DegreeVelocity + 1.0) * (DegreeVelocity + (double)D) / (double)D;
                        } else if ((GridType == typeof(Square)) || (GridType == typeof(Cube))) {
                            PenaltyBase = (DegreeVelocity + 1.0) * (DegreeVelocity + 1.0);
                        } else {
                            throw new NotImplementedException("Unknown RefElement");
                        }

                        if (this.Control.beta < 0.0) {
                            throw new ArithmeticException("Illegal setting in control object: 'beta' is out of range, must be non-negative.");
                        }
                        if (this.Control.Reynolds <= 0.0) {
                            throw new ArithmeticException("Illegal setting in control object: 'Reynolds' is out of range, must be strictly positive.");
                        }
                        if (this.Control.beta > 0.0) {
                            var Visc = new swipViscosity_Term1(
<<<<<<< HEAD
                                this.Control.ViscousPenaltyScaling * PenaltyBase, //GridData.Cells.cj, // / this.Control.beta
=======
                                this.Control.ViscousPenaltyScaling * PenaltyBase,
>>>>>>> root-public
                                d,
                                D,
                                BcMap,
                                ViscosityOption.ConstantViscosityDimensionless,
                                reynolds: this.Control.Reynolds / this.Control.beta);
                            comps.Add(Visc);
                        }

                        // extra stress divergence part:
<<<<<<< HEAD
                        //comps.Add(new StressDivergence_CentralDifference(d, BcMap, this.Control.Reynolds, this.Control.Penalty1, this.Control.Penalty2));
                        //comps.Add(new StressDivergence_Burman(d, BcMap, this.Control.Reynolds, this.Control.Penalty1, this.Control.Penalty2));
=======
>>>>>>> root-public
                        comps.Add(new StressDivergence_Cockburn(d, BcMap, this.Control.Reynolds, this.Control.Penalty1, this.Control.Penalty2));
                    }

                    // Continuum equation
                    // ===============================================================================
                    for (int d = 0; d < D; d++) {
                        XOP.EquationComponents["div"].Add(new Divergence_DerivativeSource(d, D));
                        XOP.EquationComponents["div"].Add(new Divergence_DerivativeSource_Flux(d, BcMap));

                        //Pressure stabilization for LDG
<<<<<<< HEAD
                        var presStab = new PressureStabilization(this.Control.PresPenalty2, this.Control.Reynolds);
                        XOP.EquationComponents["div"].Add(presStab);
=======
                        //var presStab = new PressureStabilization(this.Control.PresPenalty2, this.Control.Reynolds);
                        //XOP.EquationComponents["div"].Add(presStab);
>>>>>>> root-public
                    }

                    // Constitutive equations
                    // ===============================================================================

                    // Identity part
                    XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Identity(0));
                    XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Identity(1));
                    XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Identity(2));

                    //Convective part
                    XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_CellWiseForm(0, 0, BcMap, this.Control.Weissenberg, this.Control.alpha));
                    XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_CellWiseForm(0, 1, BcMap, this.Control.Weissenberg, this.Control.alpha));
                    XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_CellWiseForm(1, 1, BcMap, this.Control.Weissenberg, this.Control.alpha));

                    //Objective Part

                    // Objective PARAM ÜBERPRÜFEN!!!!!!
                    // GradU as params
                    XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Objective(0, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam, this.Control.StressPenalty));
                    XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Objective(1, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam, this.Control.StressPenalty));
                    XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Objective(2, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam, this.Control.StressPenalty));

                    //T as params
                    XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Objective_Tparam(0, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam));
                    XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Objective_Tparam(1, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam));
                    XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Objective_Tparam(2, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam));

<<<<<<< HEAD
                    //GradU and T as params
                    //XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Objective_allparam(0, BcMap, this.Control.Weissenberg));
                    //XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Objective_allparam(1, BcMap, this.Control.Weissenberg));
                    //XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Objective_allparam(2, BcMap, this.Control.Weissenberg));

                    //// GradU as new variable with L-GradU = 0
                    //XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Objective_withGrad(0, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam));
                    //XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Objective_withGrad(1, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam));
                    //XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Objective_withGrad(2, BcMap, this.Control.Weissenberg, this.Control.ObjectiveParam));

=======
>>>>>>> root-public
                    // Viscous Part
                    XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Viscosity(0, BcMap, this.Control.beta, this.Control.Penalty1));
                    XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Viscosity(1, BcMap, this.Control.beta, this.Control.Penalty1));
                    XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Viscosity(2, BcMap, this.Control.beta, this.Control.Penalty1));

                    // artificial diffusion part
                    if (this.Control.UseArtificialDiffusion == true) {
                        XOP.EquationComponents["constitutiveXX"].Add(new ConstitutiveEqns_Diffusion(this.StressXX.Basis.Degree, Grid.SpatialDimension, ((GridData)GridData).Cells.cj, VariableNames.StressXX));
                        XOP.EquationComponents["constitutiveXY"].Add(new ConstitutiveEqns_Diffusion(this.StressXY.Basis.Degree, Grid.SpatialDimension, ((GridData)GridData).Cells.cj, VariableNames.StressXY));
                        XOP.EquationComponents["constitutiveYY"].Add(new ConstitutiveEqns_Diffusion(this.StressYY.Basis.Degree, Grid.SpatialDimension, ((GridData)GridData).Cells.cj, VariableNames.StressYY));
                    }

<<<<<<< HEAD
                    // Velocity Gradient Equations
                    // ===============================================================================

                    //XOP.EquationComponents["velocitygradXX"].Add(new VelocityGradXX(0, BcMap));
                    //XOP.EquationComponents["velocitygradXY"].Add(new VelocityGradXX(1, BcMap));
                    //XOP.EquationComponents["velocitygradYX"].Add(new VelocityGradXX(2, BcMap));
                    //XOP.EquationComponents["velocitygradYY"].Add(new VelocityGradXX(3, BcMap));

                    //XOP.EquationComponents["velocitygradXX"].Add(new VelocityGrad_SU(0, BcMap));
                    //XOP.EquationComponents["velocitygradXY"].Add(new VelocityGrad_SU(1, BcMap));
                    //XOP.EquationComponents["velocitygradYX"].Add(new VelocityGrad_SU(2, BcMap));
                    //XOP.EquationComponents["velocitygradYY"].Add(new VelocityGrad_SU(3, BcMap));

                    // Build spatial operator
                    //================================================================================
                    XOP.Commit();
                    //XOP2.Commit();




                    // create timestepper
                    // ------------------

                    // LEVEL-SET - Dummy, not needed for non-Level-set calculations
                    //_______________________________________________________________________________________________
                    LevelSetHandling lsh = LevelSetHandling.None;
                    //_______________________________________________________________________________________________
=======
                    // Build spatial operator
                    XOP.Commit();


                    // create timestepper
                    //===============================================================

                    // level set - Not needed for non-Level-set calculations
                    LevelSetHandling lsh = LevelSetHandling.None;
>>>>>>> root-public


                    SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;

                    if (!this.Control.Stokes) {
                        SpatialOp = SpatialOperatorType.Nonlinear;
                    }

                    int bdfOrder;
                    if (this.Control.Timestepper_Scheme == RheologyControl.TimesteppingScheme.CrankNicolson)
                        bdfOrder = -1;
                    //else if (this.Control.Timestepper_Scheme == IBM_Control.TimesteppingScheme.ExplicitEuler)
                    //    bdfOrder = 0;
                    else if (this.Control.Timestepper_Scheme == RheologyControl.TimesteppingScheme.ImplicitEuler)
                        bdfOrder = 1;
                    else if (this.Control.Timestepper_Scheme.ToString().StartsWith("BDF"))
                        bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
                    else
<<<<<<< HEAD
                        throw new NotImplementedException("todo");

                    //PlotCurrentState(0, new TimestepNumber(new int[] { 0, 0 }), 2);

                    m_BDF_Timestepper = new XdgBDFTimestepping(ArrayTools.Cat(this.Velocity.Current, this.Pressure, this.StressXX, this.StressXY, this.StressYY), //this.VelocityXGradientX, this.VelocityXGradientY, this.VelocityYGradientX, this.VelocityYGradientY
                        ArrayTools.Cat(this.ResidualMomentum, this.ResidualConti, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY), //this.ResidualGradXX, this.ResidualGradXY, this.ResidualGradYX, this.ResidualGradYY
=======
                        throw new NotImplementedException("The chosen timestepper is not implemented!");


                    m_BDF_Timestepper = new XdgBDFTimestepping(ArrayTools.Cat(this.Velocity.Current, this.Pressure, this.StressXX, this.StressXY, this.StressYY),
                        ArrayTools.Cat(this.ResidualMomentum, this.ResidualConti, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY),
>>>>>>> root-public
                        LsTrk, false,
                        DelComputeOperatorMatrix, DelUpdateLevelset,
                        bdfOrder,
                        lsh,
<<<<<<< HEAD
                        MassMatrixShapeandDependence.IsIdentity,
=======
                        MassMatrixShapeandDependence.IsNonIdentity,
>>>>>>> root-public
                        SpatialOp,
                        MassScale,
                        this.MultigridOperatorConfig, base.MultigridSequence,
                        this.FluidSpecies, 1, // no hmf order required.
                        0, false); //HARDCODED AGGLOMERATION FACTOR -> NOT NEEDED FOR NON-LEVELSET
                    m_BDF_Timestepper.m_ResLogger = base.ResLogger;
                    m_BDF_Timestepper.m_ResidualNames = ArrayTools.Cat(this.ResidualMomentum.Select(f => f.Identification),
<<<<<<< HEAD
                        ResidualConti.Identification, ResidualStressXX.Identification, ResidualStressXY.Identification, ResidualStressYY.Identification); //ResidualGradXX.Identification, ResidualGradXY.Identification, ResidualGradYX.Identification, ResidualGradYY.Identification
=======
                        ResidualConti.Identification, ResidualStressXX.Identification, ResidualStressXY.Identification, ResidualStressYY.Identification);
>>>>>>> root-public
                }
                m_BDF_Timestepper.Config_linearSolver = this.Control.LinearSolver;
                m_BDF_Timestepper.Config_UnderRelax = this.Control.UnderRelax;
                m_BDF_Timestepper.Config_NonlinearSolver = this.Control.NonlinearMethod;
<<<<<<< HEAD
                m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;
=======
                //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;
>>>>>>> root-public
                //m_BDF_Timestepper.CustomIterationCallback += this.CoupledIterationCallback;

            }
        }

        bool solveVelocity = true;

        double VelocitySolver_ConvergenceCriterion = 1e-5;

        double StressSolver_ConvergenceCriterion = 1e-5;


        /// <summary>
        /// customizable callback routine for the handling of the coupled level-set iteration
        /// </summary>
        /// <param name="iterIndex"></param>
        /// <param name="currentSol"></param>
        /// <param name="currentRes"></param>
        /// <param name="Mgop"></param>
        protected void CoupledIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop)
        {

            var R = new CoordinateVector(this.CurrentSolution.Mapping.Fields.ToArray());
            Mgop.TransformRhsFrom(R, currentRes);
            int NF = R.Mapping.Fields.Count();

            double VelocityL2Res = 0.0;
            double StressL2Res = 0.0;

            for (int i = 0; i < NF; i++) {
                double L2Res = R.Mapping.Fields[i].L2Norm();
                if (i < 3) {
                    VelocityL2Res += L2Res;
                }
                else {
                    StressL2Res += L2Res;
                }
            }

<<<<<<< HEAD
            //double ResidualNorm = currentRes.L2NormPow2().MPISum().Sqrt();
            //if (this.Control.UsePerssonSensor == true && this.Control.UseArtificialDiffusion == true) {
            //    perssonsensor.Update(StressXX);
            //}

=======
>>>>>>> root-public
            if (solveVelocity && VelocityL2Res < this.VelocitySolver_ConvergenceCriterion) {
                this.solveVelocity = false;
            }
            else if (!solveVelocity && StressL2Res < this.StressSolver_ConvergenceCriterion) {
                this.solveVelocity = true;
            }

        }



        /// <summary>
        /// Computation of operator matrix used by the timestepper (<see cref="m_BDF_Timestepper"/>).
        /// </summary>
        protected virtual void DelComputeOperatorMatrix(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {
<<<<<<< HEAD
=======

>>>>>>> root-public
            BlockMsrMatrix OutputMatrix;
            double[] OutputAffine;

            // parameters...
            int D = this.LsTrk.GridDat.SpatialDimension;

<<<<<<< HEAD
=======

>>>>>>> root-public
            AssembleMatrix(out OutputMatrix, out OutputAffine, CurrentState, OpMatrix != null);
            if (OpMatrix != null) {
                OpMatrix.Clear();
                OpMatrix.Acc(1.0, OutputMatrix);
            }
<<<<<<< HEAD
            OpAffine.Clear();
            OpAffine.AccV(1.0, OutputAffine);

            if (Control.GravitySource == true)
            {
                // Gravity Source (default should be zero!)
=======

            OpAffine.Clear();
            OpAffine.AccV(1.0, OutputAffine);

            // Gravity Source (default should be zero!)
            if (Control.GravitySource == true)
            {
>>>>>>> root-public
                bool test = false;

                if (this.Control.GravityX != null && this.Control.GravityY != null)
                {
                    Gravity[0].ProjectField(this.Control.GravityX.Vectorize(0.0));
                    Gravity[1].ProjectField(this.Control.GravityY.Vectorize(0.0));
                    int[] MomEqIdx = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 0, 1);
                    OpAffine.AccV(-1.0, this.Gravity.CoordinateVector, MomEqIdx, default(int[]));
                    test = true;
                }

<<<<<<< HEAD
                if (this.Control.GravityXX != null && this.Control.GravityXY != null && this.Control.GravityYY != null)
                {
                    GravityXX.ProjectField(this.Control.GravityXX.Vectorize(0.0));
                    int[] ConstEqIdx1 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 3);
                    OpAffine.AccV(-1.0, this.GravityXX.CoordinateVector, ConstEqIdx1, default(int[]));

                    GravityXY.ProjectField(this.Control.GravityXY.Vectorize(0.0));
                    int[] ConstEqIdx2 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 4);
                    OpAffine.AccV(-1.0, this.GravityXY.CoordinateVector, ConstEqIdx2, default(int[]));

                    GravityYY.ProjectField(this.Control.GravityYY.Vectorize(0.0));
                    int[] ConstEqIdx3 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 5);
                    OpAffine.AccV(-1.0, this.GravityYY.CoordinateVector, ConstEqIdx3, default(int[]));
                    test = true;
                }

                if (this.Control.GravityDiv != null)
                {
                    GravityDiv.ProjectField(this.Control.GravityDiv.Vectorize(0.0));
                    int[] ContiEqIdx = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 2);
                    OpAffine.AccV(-1.0, this.GravityDiv.CoordinateVector, ContiEqIdx, default(int[]));
                    test = true;
                }
=======
                //if (this.Control.GravityXX != null && this.Control.GravityXY != null && this.Control.GravityYY != null)
                //{
                //    GravityXX.ProjectField(this.Control.GravityXX.Vectorize(0.0));
                //    int[] ConstEqIdx1 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 3);
                //    OpAffine.AccV(-1.0, this.GravityXX.CoordinateVector, ConstEqIdx1, default(int[]));

                //    GravityXY.ProjectField(this.Control.GravityXY.Vectorize(0.0));
                //    int[] ConstEqIdx2 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 4);
                //    OpAffine.AccV(-1.0, this.GravityXY.CoordinateVector, ConstEqIdx2, default(int[]));

                //    GravityYY.ProjectField(this.Control.GravityYY.Vectorize(0.0));
                //    int[] ConstEqIdx3 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 5);
                //    OpAffine.AccV(-1.0, this.GravityYY.CoordinateVector, ConstEqIdx3, default(int[]));
                //    test = true;
                //}

                //if (this.Control.GravityDiv != null)
                //{
                //    GravityDiv.ProjectField(this.Control.GravityDiv.Vectorize(0.0));
                //    int[] ContiEqIdx = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 2);
                //    OpAffine.AccV(-1.0, this.GravityDiv.CoordinateVector, ContiEqIdx, default(int[]));
                //    test = true;
                //}
>>>>>>> root-public

                if (!test)
                {
                    throw new ApplicationException("Gravity is true, but no values set!");
                }
            }
<<<<<<< HEAD

            // create mappings
            //var codMap = this.CurrentResidual.Mapping;
            //var domMap = this.CurrentSolution.Mapping;

            //double ResidualNorm = CurrentResidual.L2NormPow2().MPISum().Sqrt();
            //SinglePhaseField ResFieldvelX = (SinglePhaseField)codMap[0];

            //PlotCurrentState(phystime, 0, m_opt.SuperSampling);


=======
>>>>>>> root-public
        }

        /// <summary>
        /// Dummy function for level-set update, not used in this application, but required by the timestepper (<see cref="m_BDF_Timestepper"/>).
        /// </summary>
        public virtual double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {
            return 0.0;
        }


<<<<<<< HEAD
        int m_iterationCounter = 0;


        // Build and solve system
        //=================================================================
        /// <summary>
        /// Depending on settings <see cref="IBM_Control.Option_Timestepper"/>, computes either one timestep or a steady-state solution.
=======
        // Build and solve system
        //=================================================================
        /// <summary>
        /// Depending on settings, computes either one timestep or a steady-state solution.
>>>>>>> root-public
        /// </summary>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {

                if (this.Control.OperatorMatrixAnalysis == true) {
<<<<<<< HEAD
                    SpatialOperatorMatrixAnalysis(false, this.Control.AnalysisLevel);
=======
                    SpatialOperatorAnalysis.SpatialOperatorMatrixAnalysis(false, this.Control.AnalysisLevel);
>>>>>>> root-public
                }

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;

                base.ResLogger.TimeStep = TimestepInt;

                dt = base.GetFixedTimestep();


                int NoIncrementTimestep;

                Console.WriteLine("Instationary solve, timestep #{0}, dt = {1} ...", TimestepNo, dt);
                bool m_SkipSolveAndEvaluateResidual = this.Control.SkipSolveAndEvaluateResidual;
<<<<<<< HEAD
                //PlotCurrentState(phystime, TimestepNo);
=======
>>>>>>> root-public
                m_BDF_Timestepper.Config_SolverConvergenceCriterion = Control.ConvCrit;
                m_BDF_Timestepper.Config_MaxIterations = Control.MaxIter;
                m_BDF_Timestepper.Config_MinIterations = Control.MinIter;

                if (Control.RaiseWeissenberg == true) {

                    currentWeissenberg = 0.0;

                    if (Control.Weissenberg != 0.0) {

                        if (Control.WeissenbergIncrement != 0.0) {
                            NoIncrementTimestep = (int)(Control.Weissenberg / Control.WeissenbergIncrement);
                        } else {
                            throw new ArgumentException("Raise Weissenberg is turned on, but WeissenbergIncrement is zero!");
                        }

                    } else {
                        throw new ArgumentException("Raise Weissenberg is turned on, but aim Weissenberg is 0.0 (Newtonian)!");
                    }

                    for (int i = 0; i <= NoIncrementTimestep; i++) {

<<<<<<< HEAD
                        artificialMaxViscosity = 1.0;

                        for(int j = 0; j<3; j++) {

                            if (Control.UsePerssonSensor == true) {
                                perssonsensor.Update(StressXX);
                            }

=======
                        if(Control.UseArtificialDiffusion == true){

                        
                        artificialMaxViscosity = 1.0;

                            for (int j = 0; j < 3; j++) {

                                if (Control.UsePerssonSensor == true) {
                                    perssonsensor.Update(StressXX);
                                } else {
                                    throw new ArgumentException("artificial viscosity is turned on, but Persson sensor is turned off!");
                                }

                                m_BDF_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                                //this.ResLogger.NextTimestep(false);

                                // this evaluation must later out of this loop. now here for comparing resluts with  
                                PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                                SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                                if (Control.Bodyforces == true) {

                                    double[] force = IBMSolverUtils.GetForces_BoundaryFitted(VelocityXGradient, VelocityYGradient, StressXX, StressXY, StressYY, Pressure, LevSetTrk, 1 / Control.Reynolds, Control.beta);
                                    Console.WriteLine();
                                    Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                    Console.WriteLine();

                                }

                                artificialMaxViscosity = artificialMaxViscosity - 0.5;
                            }
                        } else {
>>>>>>> root-public
                            m_BDF_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                            //this.ResLogger.NextTimestep(false);

                            // this evaluation must later out of this loop. now here for comparing resluts with  
                            PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                            SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                            if (Control.Bodyforces == true) {

<<<<<<< HEAD
                                    double[] force = IBMSolverUtils.GetForces_BoundaryFitted(VelocityXGradient, VelocityYGradient, StressXX, StressXY, StressYY, Pressure, LevSetTrk, 1 / Control.Reynolds, Control.beta);
                                    Console.WriteLine();
                                    Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                    Console.WriteLine();

                            }

                            artificialMaxViscosity = artificialMaxViscosity - 0.5;
=======
                                double[] force = IBMSolverUtils.GetForces_BoundaryFitted(VelocityXGradient, VelocityYGradient, StressXX, StressXY, StressYY, Pressure, LevSetTrk, 1 / Control.Reynolds, Control.beta);
                                Console.WriteLine();
                                Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                Console.WriteLine();

                            }
>>>>>>> root-public
                        }

                        ChangeMesh = Control.AdaptiveMeshRefinement;
                        while (ChangeMesh == true) {
                            this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                            perssonsensor.Update(StressXX);
                            PlotCurrentState(phystime, TimestepNo);
                        }

                        if (currentWeissenberg < Control.Weissenberg) {
                            currentWeissenberg = currentWeissenberg + Control.WeissenbergIncrement;
                            Console.WriteLine();
                            Console.WriteLine("Raise Weissenberg number to " + currentWeissenberg);
                            Console.WriteLine();
                        }

                    }
                } else { 

                    currentWeissenberg = Control.Weissenberg;

<<<<<<< HEAD
                    artificialMaxViscosity = 1.0;

                    for (int j = 0; j < 3; j++) {

                        if (Control.UsePerssonSensor == true) {
                            perssonsensor.Update(StressXX);
                        }

                        m_BDF_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                        // this evaluation must later out of this loop. now here for comparing resluts with  
                        //PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                        //SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);
=======
                    if (Control.UseArtificialDiffusion == true) {
                        artificialMaxViscosity = 1.0;

                        for (int j = 0; j < 3; j++) {

                            if (Control.UsePerssonSensor == true) {
                                perssonsensor.Update(StressXX);
                            } else {
                                throw new ArgumentException("artificial viscosity is turned on, but Persson sensor is turned off!");
                            }

                            m_BDF_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                            // this evaluation must later out of this loop. now here for comparing resluts with  
                            //PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                            //SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                            if (Control.Bodyforces == true) {

                                double[] force = IBMSolverUtils.GetForces_BoundaryFitted(VelocityXGradient, VelocityYGradient, StressXX, StressXY, StressYY, Pressure, LevSetTrk, 1 / Control.Reynolds, Control.beta);
                                Console.WriteLine();
                                Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                                Console.WriteLine();

                            }

                            artificialMaxViscosity = artificialMaxViscosity - 0.5;
                        }

                        ChangeMesh = Control.AdaptiveMeshRefinement;
                        while (ChangeMesh == true) {
                            this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                        }
                    } else {

                        m_BDF_Timestepper.Solve(phystime, dt, m_SkipSolveAndEvaluateResidual);

                        // simple implicit Euler solve for debugging and excluding the bdf timestepper
                        //____________________________________________________________________________________
                        //Console.WriteLine("CAREFUL! Simple implicit Euler unsteady solve for Debugging!");
                        //var map = this.CurrentSolution.Mapping;
                        //var Mtx = new BlockMsrMatrix(map, map);
                        //double[] b = new double[map.LocalLength];
                        //this.DelComputeOperatorMatrix(Mtx, b, map, map.Fields.ToArray(), null, phystime + dt);

                        //double[] RHS = new double[map.LocalLength];
                        //RHS.AccV(-1, b);
                        //int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                        //int Np = this.Velocity.Current[0].Basis.Length;
                        //double oodt = 1.0 / dt;
                        //for (int j = 0; j < J; j++) { // loop over cells
                        //    for (int iVar = 0; iVar < 2; iVar++) { // loop over VelX, VelY
                        //        for (int n = 0; n < Np; n++) {
                        //            int iLoc = map.LocalUniqueCoordinateIndex(iVar, j, n);
                        //            RHS[iLoc] += oodt * this.CurrentSolution[iLoc];

                        //            int iGlob = map.GlobalUniqueCoordinateIndex(iVar, j, n);
                        //            Mtx[iGlob, iGlob] += oodt;
                        //        }
                        //    }
                        //}

                        //Mtx.Solve_Direct(this.CurrentSolution, RHS);
                        //____________________________________________________________________________________________
>>>>>>> root-public

                        if (Control.Bodyforces == true) {

                            double[] force = IBMSolverUtils.GetForces_BoundaryFitted(VelocityXGradient, VelocityYGradient, StressXX, StressXY, StressYY, Pressure, LevSetTrk, 1 / Control.Reynolds, Control.beta);
                            Console.WriteLine();
                            Console.WriteLine("Force in x:" + force[0] + ", force in y:" + force[1]);
                            Console.WriteLine();

                        }
<<<<<<< HEAD

                        artificialMaxViscosity = artificialMaxViscosity - 0.5;
                    }

                    ChangeMesh = Control.AdaptiveMeshRefinement;
                    while (ChangeMesh == true) {
                        this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                    }
                }



                //if (Control.RaiseWeissenberg == true)
                //{
                //    //int timeintervall = 20;
                //    //if (TimestepNo > timeintervall)
                //    //{
                //        if (Control.Weissenberg < 10)
                //        {
                //            var newWeissenberg = Control.Weissenberg + 0.1;
                //            Console.WriteLine();
                //            Console.WriteLine("Raise Weisenberg number from " + Control.Weissenberg + " to " + newWeissenberg);
                //            Console.WriteLine();
                //            Control.Weissenberg = newWeissenberg;
                //            m_BDF_Timestepper = null;
                //            CreateEquationsAndSolvers(null);
                            
                //       //     timeintervall += TimestepNo;
                //        }
                //   // }
                //}

=======
                    }
                }

>>>>>>> root-public
                if (Control.ComputeL2Error == true) {
                    this.ComputeL2Error();
                }

                this.ResLogger.NextTimestep(false);

                return dt;


            }
        }

<<<<<<< HEAD
        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double rho = 1; // this.Control.PhysicalParameters.rho_A;

                int D = this.GridData.SpatialDimension;

                double[] _rho = new double[D + 4];
                _rho.SetAll(rho);
                //No MassMatrix for the pressure
                _rho[D] = 0;

                _rho[D + 1] = 1;
                _rho[D + 2] = 1;
                _rho[D + 3] = 1;
                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho);

                return R;
            }
        }

        void ParameterUpdate(IEnumerable<DGField> CurrentState, IEnumerable<DGField> ParameterVar) {
            var U0 = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());
            var Stress0 = new VectorField<SinglePhaseField>(CurrentState.Skip(D+1).Take(3).Select(F => (SinglePhaseField)F).ToArray());

            //var Params = ArrayTools.Cat<DGField>(U0_U0mean, VelocityXGradient, VelocityYGradient, Stress0);
            //SinglePhaseField[] param_U0 = ParameterVar.Skip(D).Take(D).ToArray();

            //SinglePhaseField[] U0_U0mean;
=======

        void ParameterUpdate(IEnumerable<DGField> CurrentState, IEnumerable<DGField> ParameterVar) {

            var U0 = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());
            var Stress0 = new VectorField<SinglePhaseField>(CurrentState.Skip(D+1).Take(3).Select(F => (SinglePhaseField)F).ToArray());

>>>>>>> root-public
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
                SinglePhaseField[] __VelocityXGradient = ParameterVar.Skip(2*D).Take(D).Select(f => f as SinglePhaseField).ToArray();
                SinglePhaseField[] __VelocityYGradient = ParameterVar.Skip(3*D).Take(D).Select(f => f as SinglePhaseField).ToArray();
                Debug.Assert(ArrayTools.AreEqual(__VelocityXGradient, VelocityXGradient.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));
                Debug.Assert(ArrayTools.AreEqual(__VelocityYGradient, VelocityYGradient.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));

                VelocityXGradient.Clear();
                VelocityXGradient.GradientByFlux(1.0, U0[0]);
                VelocityYGradient.Clear();
                VelocityYGradient.GradientByFlux(1.0, U0[1]);
            }

            if (this.Control.UseArtificialDiffusion == true) {

                SinglePhaseField __ArtificialViscosity = ParameterVar.Skip(5 * D + 1).Take(1).Select(f => f as SinglePhaseField).ToArray()[0];
                if (!object.ReferenceEquals(this.artificalViscosity, __ArtificialViscosity))
                    throw new ApplicationException();

                ArtificialViscosity.ProjectArtificalViscosityToDGField(__ArtificialViscosity, perssonsensor, this.Control.SensorLimit, artificialMaxViscosity);
<<<<<<< HEAD

                //ContinuousDGField continuousField = new ContinuousDGField(artificalViscosity.Basis);
                //continuousField.ProjectDGField(1.0, __ArtificialViscosity);
                //__ArtificialViscosity.Clear();
                //continuousField.AccToDGField(1.0, __ArtificialViscosity);
            }
        }
        

        public void AssembleMatrix(out BlockMsrMatrix OpMatrix, out double[] OpAffine, DGField[] CurrentState, bool Linearization) {

            // check:
            D = this.GridData.SpatialDimension;

=======
            }
        }


        /// <summary>
        /// Computation of operator matrix to be used by DelComputeOperatorMatrix, the SpatialOperatorAnalysis and sone unit tests(<see cref="m_BDF_Timestepper"/>).
        /// </summary>
        public void AssembleMatrix(out BlockMsrMatrix OpMatrix, out double[] OpAffine, DGField[] CurrentState, bool Linearization) {

            D = this.GridData.SpatialDimension;
>>>>>>> root-public
            var U0 = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());
            var Stress0 = new VectorField<SinglePhaseField>(CurrentState.Skip(D+1).Take(3).Select(F => (SinglePhaseField)F).ToArray());

            if (U0.Count != D)
<<<<<<< HEAD
                throw new ArgumentException();

            if (Stress0.Count != D+1)
                throw new ArgumentException();

           

            // parameters...
=======
                throw new ArgumentException("Spatial dimesion and number of velocity parameter components does not match!");

            if (Stress0.Count != D+1)
                throw new ArgumentException("Spatial dimesion and number of stress parameter components does not match!");      


            // parameters
            //============================================================
>>>>>>> root-public
            SinglePhaseField[] U0_U0mean;
            if (this.U0MeanRequired) {
                Basis U0meanBasis = new Basis(GridData, 0);
                VectorField<SinglePhaseField> U0mean = new VectorField<SinglePhaseField>(D, U0meanBasis, "U0mean_", SinglePhaseField.Factory);
                U0mean.Clear();
                                
                U0_U0mean = ArrayTools.Cat<SinglePhaseField>(U0, U0mean);
            } else {
                U0_U0mean = new SinglePhaseField[2 * D];
            }

            var Params = ArrayTools.Cat<DGField>(U0_U0mean, VelocityXGradient, VelocityYGradient, Stress0, artificalViscosity);

<<<<<<< HEAD
            // create mappings
=======

            // create mappings
            //==========================================================
>>>>>>> root-public
            var codMap = this.CurrentResidual.Mapping;
            var domMap = this.CurrentSolution.Mapping;


<<<<<<< HEAD
            if (Linearization) {
                // +++++++++++++++++++++++++++++++++++++++
                // provide a linearization of the operator
                // +++++++++++++++++++++++++++++++++++++++

                bool useJacobianForOperatorMatrix = true;

                if (this.Control.NonlinearMethod == NonlinearSolverMethod.Picard)
                useJacobianForOperatorMatrix = false;
=======
            // provide a linearization of the operator
            //===========================================================
            if (Linearization) {

                bool useJacobianForOperatorMatrix = true;

                //if (this.Control.NonlinearMethod == NonlinearSolverMethod.Picard)
                    useJacobianForOperatorMatrix = false;
>>>>>>> root-public

                // create matrix and affine vector:
                OpMatrix = new BlockMsrMatrix(codMap, domMap);
                OpAffine = new double[codMap.LocalLength];

<<<<<<< HEAD
                // 'custom' Linearization 
                // ----------------------

=======

                // 'custom' Linearization 
>>>>>>> root-public
                if (!useJacobianForOperatorMatrix) {

                    var Mbuilder = XOP.GetMatrixBuilder(domMap, Params, codMap);
                    this.ParameterUpdate(domMap.Fields, Params);
                    Mbuilder.ComputeMatrix(OpMatrix, OpAffine);
                    Mbuilder.OperatorCoefficients.UserDefinedValues.Add("Weissenbergnumber", currentWeissenberg);

                } else {
                    // Finite Difference Linearization
<<<<<<< HEAD
                    // -------------------------------

                    //var CheckMatrix = new BlockMsrMatrix(codMap, domMap);
                    //var CheckAffine = new double[codMap.GlobalCount];



                    var FDbuilder = XOP.GetFDJacobianBuilder(domMap, Params, codMap, this.ParameterUpdate);
                    FDbuilder.ComputeMatrix(OpMatrix, OpAffine);
                    FDbuilder.OperatorCoefficients.UserDefinedValues.Add("Weissenbergnumber", currentWeissenberg);
                }

                //var ErrMatrix = OpMatrix.CloneAs();
                //var ErrAffine = OpAffine.CloneAs();
                //ErrMatrix.Acc(-1.0, CheckMatrix);
                //ErrAffine.AccV(-1.0, CheckAffine);
                //double LinfMtx = ErrMatrix.InfNorm();
                //double L2Aff = ErrAffine.L2NormPow2().MPISum().Sqrt();
                //Console.WriteLine("-----  Matrix/Affine delta norm {0} {1}", LinfMtx, L2Aff);

                //extract Matrix
                //Mbuilder.SaveToTextFileSparse("C:\\BoSSS-code\\internal\\src\\experimental\\L4-application\\RheologySolver2\\bin\\Release\\Mbuilder.txt");
                //FDbuilder.SaveToTextFileSparse("C:\\BoSSS-code\\internal\\src\\experimental\\L4-application\\RheologySolver2\\bin\\Release\\FDbuilder.txt");



                if (!this.BcMap.DirichletPressureBoundary)
                {
                    SolverUtils.SetRefPtPressure_Matrix(OpMatrix, this.PressureRefCellIndex);
                    SolverUtils.SetRefPtPressure_Rhs(OpAffine, PressureRefCellIndex, codMap.i0);
=======

                    var FDbuilder = XOP.GetFDJacobianBuilder(domMap, Params, codMap, this.ParameterUpdate);
                    FDbuilder.ComputeMatrix(OpMatrix, OpAffine);

                    // FDJacobian has (Mx +b) as RHS, for unsteady calc. we must subtract Mx for real affine Vector!
                    OpMatrix.SpMV(-1.0, new CoordinateVector(CurrentState), 1.0, OpAffine);

                    FDbuilder.OperatorCoefficients.UserDefinedValues.Add("Weissenbergnumber", currentWeissenberg);
                }

                // Set Pressure Reference Point
                //======================================================
                if (!this.BcMap.DirichletPressureBoundary) {
                    if (OpMatrix != null) {

                        IBMSolverUtils.SetPressureReferencePoint(
                            CurrentSolution.Mapping,
                            this.GridData.SpatialDimension,
                            this.LsTrk,
                            OpMatrix, OpAffine);
                        //OpMatrix.SaveToTextFileSparse("OpMatrix_3D");
                    } else {
                        IBMSolverUtils.SetPressureReferencePointResidual(
                            new CoordinateVector(CurrentState),
                            this.GridData.SpatialDimension,
                            this.LsTrk,
                            OpAffine);
                    }
>>>>>>> root-public
                }

                OpMatrix.CheckForNanOrInfM();
                OpAffine.CheckForNanOrInfV();
            } else {
<<<<<<< HEAD
                // +++++++++++++++++++++++++++++++++++++++
                // explicit evaluation of the operator
                // +++++++++++++++++++++++++++++++++++++++
=======

                // explicit evaluation of the operator
                //========================================================
>>>>>>> root-public

                OpMatrix = null;
                OpAffine = new double[codMap.LocalLength];
                var eval = XOP.GetEvaluatorEx(CurrentState, Params, codMap);
                this.ParameterUpdate(eval.DomainFields.Fields, Params);
                eval.OperatorCoefficients.UserDefinedValues.Add("Weissenbergnumber", currentWeissenberg);

                eval.Evaluate(1.0, 1.0, OpAffine);

            }

        }

        /// <summary>
<<<<<<< HEAD
        /// Computes condition number, etc. of the current system matrix.
        /// </summary>
        /// <param name="CheckAssertions"></param>
        /// <param name="AnalysisLevel">
        /// - equal 0: check that pressure gradient and velocity divergence are transpose
        /// - equal 1: in addition, positive definiteness test.
        /// - equal 2: in addition, check condition number and eigenvalues using MATLAB
        /// </param>
        public void SpatialOperatorMatrixAnalysis(bool CheckAssertions, int AnalysisLevel)
        {
            using (new FuncTrace())
            {
                int D = this.Grid.SpatialDimension;

                if (AnalysisLevel < 0 || AnalysisLevel > 2)
                    throw new ArgumentException();


                BlockMsrMatrix OpMatrix;
                double[] OpAffine;

                AssembleMatrix(out OpMatrix, out OpAffine, this.CurrentSolution.Mapping.ToArray(), true);


                // =============================
                // AnalysisLevel 0
                // =============================
                {
                    var OpMatrixT = OpMatrix.Transpose();

                    CoordinateVector TestVec = new CoordinateVector(this.CurrentSolution.Mapping.Fields.Select(f => f.CloneAs()).ToArray());

                    double testsumPos = 0.0;
                    double testsumNeg = 0.0;
                    for (int rnd_seed = 0; rnd_seed < 20; rnd_seed++)
                    {

                        // fill the pressure components of the test vector
                        TestVec.Clear();
                        Random rnd = new Random(rnd_seed);
                        DGField Pressack = TestVec.Mapping.Fields[D] as DGField;
                        int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++)
                        {
                            int N = Pressack.Basis.GetLength(j);

                            for (int n = 0; n < N; n++)
                                Pressack.Coordinates[j, n] = rnd.NextDouble();
                        }

                        // Gradient times P:
                        double[] R1 = new double[TestVec.Count];
                        OpMatrix.SpMV(1.0, TestVec, 0.0, R1);       // R1 = Grad * P
                        //Console.WriteLine("L2 of 'Grad * P': " + R1.L2Norm());

                        // transpose of Divergence times P: 
                        double[] R2 = new double[TestVec.Count];
                        OpMatrix.SpMV(1.0, TestVec, 0.0, R2);      // R2 = divT * P
                        //Console.WriteLine("L2 of 'divT * P': " + R2.L2Norm());

                        TestVec.Clear();
                        TestVec.Acc(1.0, R1);
                        TestVec.Acc(1.0, R2);


                        // analyze!
                        testsumNeg += GenericBlas.L2Dist(R1, R2);

                        R2.ScaleV(-1.0);
                        testsumPos += GenericBlas.L2Dist(R1, R2);

                    }

                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (+): " + testsumPos);
                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (-): " + testsumNeg);

                    if (CheckAssertions)
                        Assert.LessOrEqual(Math.Abs(testsumNeg), testsumPos * 1.0e-13);
                }


                // =============================
                // AnalysisLevel 1 and 2
                // =============================

                if (AnalysisLevel > 0)
                {
                    AggregationGridBasis[][] MgBasis = AggregationGridBasis.CreateSequence(this.MultigridSequence, this.CurrentSolution.Mapping.BasisS);
                    MultigridOperator mgOp = new MultigridOperator(MgBasis, this.CurrentSolution.Mapping,
                        OpMatrix, null, this.MultigridOperatorConfig);

                    // extract
                    ////////////

                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                    MsrMatrix DiffMatrix;
                    {
                        int[] VelVarIdx = D.ForLoop(d => d);

                        int[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int L = USubMatrixIdx_Row.Length;

                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(int[]), USubMatrixIdx_Col, default(int[]));

                        double DiffMatrix_sd = DiffMatrix.SymmetryDeviation();
                        Console.WriteLine("Diffusion assymetry:" + DiffMatrix_sd);
                    }

                    MsrMatrix SaddlePointMatrix;
                    {
                        int[] VelPVarIdx = new int[] { 0, 1, 2 };

                        int[] VelPSubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelPVarIdx);
                        int[] VelPSubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelPVarIdx);
                        int L = VelPSubMatrixIdx_Row.Length;

                        SaddlePointMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(SaddlePointMatrix, VelPSubMatrixIdx_Row, default(int[]), VelPSubMatrixIdx_Col, default(int[]));
                    }
                    //SaddlePointMatrix.SaveToTextFileSparse("C:\\Users\\kikker\\Documents\\MATLAB\\spm.txt");

                    MsrMatrix ConstitutiveMatrix;
                    {
                        int[] StressVarIdx = new int[] { 3, 4, 5 };

                        int[] StressSubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(StressVarIdx);
                        int[] StressSubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(StressVarIdx);
                        int L = StressSubMatrixIdx_Row.Length;

                        ConstitutiveMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(ConstitutiveMatrix, StressSubMatrixIdx_Row, default(int[]), StressSubMatrixIdx_Col, default(int[]));
                    }

                    // operator analysis
                    //////////////////////

                    bool posDef;
                    if (AnalysisLevel > 1)
                    {
                        // +++++++++++++++++++++++++++++++
                        // check condition number, etc
                        // +++++++++++++++++++++++++++++++

                        MultidimensionalArray ret = MultidimensionalArray.Create(1, 5);
                        Console.WriteLine("Calling MATLAB/Octave...");
                        using (BatchmodeConnector bmc = new BatchmodeConnector())
                        {
                            bmc.PutSparseMatrix(FullMatrix, "FullMatrix");
                            bmc.PutSparseMatrix(SaddlePointMatrix, "SaddlePointMatrix");
                            bmc.PutSparseMatrix(ConstitutiveMatrix, "ConstitutiveMatrix");
                            bmc.PutSparseMatrix(DiffMatrix, "DiffMatrix");

                            bmc.Cmd("DiffMatrix = 0.5*(DiffMatrix + DiffMatrix');");

                            bmc.Cmd("condNoFullMatrix = condest(FullMatrix);");
                            bmc.Cmd("condNoSaddlePointMatrix = condest(SaddlePointMatrix);");
                            bmc.Cmd("condNoConstitutiveMatrix = condest(ConstitutiveMatrix);");
                            bmc.Cmd("condNoDiffMatrix = condest(DiffMatrix);");

                            //bmc.Cmd("eigiMaxiSaddle = 1.0; % eigs(SaddlePointMatrix,1,'lm')");
                            //bmc.Cmd("eigiMiniSaddle = 1.0; % eigs(SaddlePointMatrix,1,'sm')");
                            //bmc.Cmd("eigiMaxiConst = 1.0; % eigs(ConstitutiveMatrix,1,'lm')");
                            //bmc.Cmd("eigiMiniConst = 1.0; % eigs(ConstitutiveMatrix,1,'sm')");
                            //bmc.Cmd("eigiMaxiDiff = 1.0; % eigs(DiffMatrix,1,'lm')");
                            //bmc.Cmd("eigiMiniDiff = 1.0; % eigs(DiffMatrix,1,'sm')");

                            bmc.Cmd("lasterr");
                            bmc.Cmd("[V,r]=chol(SaddlePointMatrix);");
                            bmc.Cmd("[V,r]=chol(ConstitutiveMatrix);");
                            bmc.Cmd("ret = [condNoFullMatrix, condNoSaddlePointMatrix, condNoConstitutiveMatrix, condNoDiffMatrix, r]"); //eigiMaxiSaddle, eigiMiniSaddle, eigiMaxiConst, eigiMiniConst, eigiMaxiDiff, eigiMiniDiff,
                            bmc.GetMatrix(ret, "ret");

                            bmc.Execute(false);
                        }

                        double condNoFullMatrix = ret[0, 0];
                        double condNoSaddlePMatrix = ret[0, 1];
                        double condNoConstitutiveMatrix = ret[0, 2];
                        double condNoDiffMatrix = ret[0, 3];
                        //double eigiMaxiSaddle = ret[0, 4];
                        //double eigiMiniSaddle = ret[0, 5];
                        //double eigiMaxiConst = ret[0, 6];
                        //double eigiMiniConst = ret[0, 7];
                        //double eigiMaxiDiff = ret[0, 8];
                        //double eigiMiniDiff = ret[0, 9];
                        posDef = ret[0, 4] == 0;

                        //Console.WriteLine("Eigenvalue range of saddle point matrix: {0} to {1}", eigiMiniSaddle, eigiMaxiSaddle);
                        //Console.WriteLine("Eigenvalue range of constitutive matrix: {0} to {1}", eigiMiniConst, eigiMaxiConst);
                        //Console.WriteLine("Eigenvalue range of diffusion matrix: {0} to {1}", eigiMiniDiff, eigiMaxiDiff);

                        Console.WriteLine("Condition number full operator: {0:0.####E-00}", condNoFullMatrix);
                        Console.WriteLine("Condition number saddle point operator: {0:0.####E-00}", condNoSaddlePMatrix);
                        Console.WriteLine("Condition number constitutive operator: {0:0.####E-00}", condNoConstitutiveMatrix);
                        Console.WriteLine("Condition number diffusion operator: {0:0.####E-00}", condNoDiffMatrix);

                        base.QueryHandler.ValueQuery("ConditionNumber", condNoFullMatrix);

                    }
                    else
                    {
                        // +++++++++++++++++++++++++++++++++++++++
                        // test only for positive definiteness
                        // +++++++++++++++++++++++++++++++++++++++

                        var SaddlePMatrixFull = SaddlePointMatrix.ToFullMatrixOnProc0();
                        var ConstMatrixFull = ConstitutiveMatrix.ToFullMatrixOnProc0();


                        posDef = true;
                        try
                        {
                            SaddlePMatrixFull.Cholesky();
                        }
                        catch (ArithmeticException)
                        {
                            posDef = false;
                        }

                        posDef = true;
                        try
                        {
                            ConstMatrixFull.Cholesky();
                        }
                        catch (ArithmeticException)
                        {
                            posDef = false;
                        }
                    }


                    double SaddlePSymm = SaddlePointMatrix.SymmetryDeviation();
                    Console.WriteLine("Symmetry deviation of saddle point matrix: " + SaddlePSymm);

                    if (posDef)
                        Console.WriteLine("Good news: Saddle point operator matrix seems to be positive definite.");
                    else
                        Console.WriteLine("WARNING: Saddle point operator matrix is not positive definite.");


                    double ConstSymm = ConstitutiveMatrix.SymmetryDeviation();
                    Console.WriteLine("Symmetry deviation of constitutive matrix: " + ConstSymm);

                    if (posDef)
                        Console.WriteLine("Good news: constitutive operator matrix seems to be positive definite.");
                    else
                        Console.WriteLine("WARNING: constitutive operator matrix is not positive definite.");

                    //if (CheckAssertions) {
                    //    if (Control.AdvancedDiscretizationOptions.ViscosityMode == ViscosityMode.FullySymmetric && Control.PhysicalParameters.IncludeConvection == false) {
                    //        Assert.IsTrue(posDef, "Positive definiteness test failed.");
                    //        double compVal = DiffMatrix.InfNorm() * 1e-13;
                    //        Assert.LessOrEqual(DiffSymm, compVal, "Diffusion matrix seems to be non-symmetric.");
                    //    }
                    //}
                }
            }
        }


        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
=======
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        public MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
>>>>>>> root-public
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
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 4];

                    // configurations for velocity
                    for (int d = 0; d < D; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            Degree = Math.Max(1, pVel - iLevel),
                            mode = this.Control.VelocityBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, pPrs - iLevel),
                        mode = this.Control.PressureBlockPrecondMode,
                        VarIndex = new int[] { D }
                    };

                    // configurations for stresses
                    for (int d = 3; d < 6; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            Degree = Math.Max(1, pStr - iLevel),
                            mode = this.Control.StressBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
<<<<<<< HEAD

                    //// configurations for velocitygradient
                    //for (int d = 6; d < 10; d++)
                    //{
                    //    configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig()
                    //    {
                    //        Degree = Math.Max(1, pStr - iLevel),
                    //        mode = this.Control.VelocityGradientBlockPrecondMode,
                    //        VarIndex = new int[] { d }
                    //    };
                    //}
=======
>>>>>>> root-public
                }

                return configs;
            }
        }

<<<<<<< HEAD
=======
        /// <summary>
        /// Plotting the current state
        /// </summary>
>>>>>>> root-public
        protected override void PlotCurrentState(double physTime, Foundation.IO.TimestepNumber timestepNo, int superSampling = 0) {
            // Standard
            DGField[] myFields = ArrayTools.Cat<DGField>(Velocity.Current, ResidualMomentum, ResidualConti, Pressure, StressXX, StressXY, StressYY, LevSet, ResidualStressXX, ResidualStressXY, ResidualStressYY); //, VelocityXGradient, VelocityYGradient, Gravity

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


<<<<<<< HEAD

=======
        /// <summary>
        /// Plotting the in interation callback
        /// </summary>
>>>>>>> root-public
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

            Tecplot.PlotFields( myFields, "Rheology-" + iterIndex.ToString(), 0.0, 2); 
        }

<<<<<<< HEAD
=======
        /// <summary>
        /// Initialising the DG fields
        /// </summary>
>>>>>>> root-public
        protected override void SetInitial() {
            base.SetInitial();
            this.LsTrk.UpdateTracker();
            CreateEquationsAndSolvers(null);
            m_BDF_Timestepper.SingleInit();
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

        }


<<<<<<< HEAD
=======
        /// <summary>
        /// performs restart
        /// </summary>
        /// <param name="Time">
        /// on exit, the physical time associated with the field state
        /// </param>
        /// <param name="TimestepNo">
        /// on exit, the physical time associated with the field state
        /// </param>
>>>>>>> root-public
        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);

            this.LsTrk.UpdateTracker();
        }

<<<<<<< HEAD

        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            //PlotCurrentState(0, new TimestepNumber(new int[] { 0, 10 }), 2);

            VelocityXGradient = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, Velocity.Current[0].Basis, "VelocityX_Gradient", SinglePhaseField.Factory);
            VelocityYGradient = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, Velocity.Current[1].Basis, "VelocityY_Gradient", SinglePhaseField.Factory);


            ////otherwise the DGfields are projected twice in refined cells!!!
            //foreach (var f in m_RegisteredFields) {
            //    f.Clear();
            //}


=======
        /// <summary>
        /// overriding the method to implement any user-specific tasks which
        /// should be carried out after a restart file has been loaded (e.g.,
        /// setting the correct time for a time-stepper)
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            VelocityXGradient = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, Velocity.Current[0].Basis, "VelocityX_Gradient", SinglePhaseField.Factory);
            VelocityYGradient = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, Velocity.Current[1].Basis, "VelocityY_Gradient", SinglePhaseField.Factory);
>>>>>>> root-public
        }


        /// <summary>
        /// Computes average velocity in case of Navier-Stokes Equations
        /// </summary>
        /// <param name="U0"></param>
        /// <param name="U0mean"></param>
        private void ComputeAverageU(VectorField<SinglePhaseField> U0, VectorField<SinglePhaseField> U0mean) {
            using(FuncTrace ft = new FuncTrace()) {
                var CC = this.LsTrk.Regions.GetCutCellMask();
                int D = this.LsTrk.GridDat.SpatialDimension;
                double minvol = Math.Pow(this.LsTrk.GridDat.Cells.h_minGlobal, D);

                int QuadDegree = this.HMForder;

                var qh = LsTrk.GetXDGSpaceMetrics(this.FluidSpecies, QuadDegree, 1).XQuadSchemeHelper;
                foreach(var Spc in this.FluidSpecies) { // loop over species...
                    //var Spc = this.LsTrk.GetSpeciesId("B"); {
                    // shadow fields
                    var U0_Spc = U0.ToArray();
                    var U0mean_Spc = U0mean.ToArray();


                    // normal cells:
                    for(int d = 0; d < D; d++) {
                        U0mean_Spc[d].AccLaidBack(1.0, U0_Spc[d], this.LsTrk.Regions.GetSpeciesMask(Spc));
                    }


                    // cut cells
                    var scheme = qh.GetVolumeQuadScheme(Spc, IntegrationDomain: this.LsTrk.Regions.GetCutCellMask());

                    var rule = scheme.Compile(this.LsTrk.GridDat, QuadDegree);
                    CellQuadrature.GetQuadrature(new int[] { D + 1 }, // vector components: ( avg_vel[0], ... , avg_vel[D-1], cell_volume )
                        this.LsTrk.GridDat,
                        rule,
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            EvalResult.Clear();
                            for(int d = 0; d < D; d++)
                                U0_Spc[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                            var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, D);
                            Vol.SetAll(1.0);
                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            for(int i = 0; i < Length; i++) {
                                int jCell = i + i0;

                                double Volume = ResultsOfIntegration[i, D];
                                if(Math.Abs(Volume) < minvol * 1.0e-12) {
                                    // keep current value
                                    // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                                } else {
                                    for(int d = 0; d < D; d++) {
                                        double IntVal = ResultsOfIntegration[i, d];
                                        U0mean_Spc[d].SetMeanValue(jCell, IntVal / Volume);
                                    }
                                }

                            }
<<<<<<< HEAD
                        }).Execute();

                }
                U0mean.ForEach(F => F.CheckForNanOrInf(true, true, true));

            }
        }




=======
                        }
                        ).Execute();
                }
                U0mean.ForEach(F => F.CheckForNanOrInf(true, true, true));
            }
        }

        /// <summary>
        /// Computes the L2 Error of all Fields compared to exact solution specified in the control file
        /// </summary>
>>>>>>> root-public
        protected void ComputeL2Error() {
            if (this.Control.ExSol_Velocity == null && this.Control.ExSol_Pressure == null && this.Control.ExSol_Stress == null) {
                // nothing to do
                return;
            }


            int D = this.GridData.SpatialDimension;

            int order = Velocity.Current[0].Basis.Degree * 2;
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(this.FluidSpecies, order).XQuadSchemeHelper;

            // Velocity error
<<<<<<< HEAD
            // ==============
=======
            // ===============================================
>>>>>>> root-public
            if (this.Control.ExSol_Velocity != null) {
                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                foreach (var spId in this.FluidSpecies) {
                    string spc = this.LsTrk.GetSpeciesName(spId);

                    L2Error_Species.Add(spc, new double[D]);

                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                    for (int d = 0; d < D; d++) {
                        L2Error_Species[spc][d] = this.Velocity.Current[d].L2Error(this.Control.ExSol_Velocity[d].Vectorize(0.0), order, scheme);
                        L2Error[d] += L2Error_Species[spc][d].Pow2();

                        base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);

                    }
                }
                L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

                for (int d = 0; d < D; d++) {
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Console.WriteLine("L2err " + VariableNames.Velocity_d(d) + " is " + L2Error[d]);                    
                }
            }


            // pressure error
<<<<<<< HEAD
            // ==============
            if (this.Control.ExSol_Pressure != null) {

                //// pass 1: mean value of pressure difference
                //double DiffInt = 0;
                //foreach (var spId in FluidSpecies) {

                //    string spc = this.LsTrk.GetSpeciesName(spId);
                //    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                //    var rule = scheme.Compile(this.GridData, order);

                //    DiffInt += this.Pressure.LxError(this.Control.ExSol_Pressure.Vectorize(0.0), (a, b) => (a - b), rule);
                //    //Volume +=  this.Pressure.GetSpeciesShadowField(spc).LxError(null, (a, b) => (1.0), rule);
                //}
                //double Volume2 = (new SubGrid(CellMask.GetFullMask(this.GridData))).Volume;
                //double PressureDiffMean = DiffInt / Volume2;


                double L2Error = 0;
                //Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                //foreach (var spId in this.FluidSpecies) {

                //    //SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                //    string spc = this.LsTrk.GetSpeciesName(spId);
                //    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                //    var rule = scheme.Compile(this.GridData, order);

                //    double IdV = this.Pressure.LxError(this.Control.ExSol_Pressure.Vectorize(0.0), (a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                //    L2Error += IdV;
                //    L2Error_Species.Add(spc, IdV.Sqrt());

                //    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
                //}


                //L2Error = L2Error.Sqrt();
=======
            // =============================================================
            if (this.Control.ExSol_Pressure != null) {

                double L2Error = 0;

>>>>>>> root-public
                L2Error = this.Pressure.L2Error(this.Control.ExSol_Pressure.Vectorize(0.0), order-1);
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
                Console.WriteLine("L2err " + VariableNames.Pressure + " is " + L2Error);
            }

            // Stress error
<<<<<<< HEAD
            // ==============
=======
            // =============================================================
>>>>>>> root-public
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
<<<<<<< HEAD

=======
>>>>>>> root-public
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
<<<<<<< HEAD
        //=============================================
=======
        //======================================================================
>>>>>>> root-public

        /// <summary>
        /// refinement indicator
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel){
<<<<<<< HEAD
            //if (j == 1) {
            //    return 1;

            //} else {
            //    return 0;
            //}
            if (this.Control.UsePerssonSensor) {

                double maxVal = this.perssonsensor.GetValue(j);
                //this.StressXX.GetExtremalValuesInCell(out double minVal, out double maxVal, j);
=======

            if (this.Control.UsePerssonSensor) {

                double maxVal = this.perssonsensor.GetValue(j);
>>>>>>> root-public

                //bound for perssonsensor should be around 1e-7 - 1e-8 that there is refinement behind the cylinder!
                double upperbound = this.Control.SensorLimit;
                double lowerbound = upperbound * 0.001;

                int DesiredLevel_j = CurrentLevel;

                if (maxVal > upperbound && DesiredLevel_j < this.Control.RefinementLevel) {

                    DesiredLevel_j = DesiredLevel_j + 1;

                } else if(maxVal < lowerbound && DesiredLevel_j > 0) {
                    DesiredLevel_j = DesiredLevel_j - 1;
                }

                return DesiredLevel_j;

            } else {
                throw new NotSupportedException("The Persson sensor is turned off. It is needed for adaptive mesh refinement!");
            }

        }

<<<<<<< HEAD
=======
        /// <summary>
        /// Adaptation of the current mesh.
        /// </summary>
>>>>>>> root-public
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {

            if (this.Control.AdaptiveMeshRefinement) {

                bool AnyChange = GridRefinementController.ComputeGridChange((GridData)(this.GridData), null, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
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
                int oldJ = this.GridData.CellPartitioning.TotalLength;

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
    }
}


