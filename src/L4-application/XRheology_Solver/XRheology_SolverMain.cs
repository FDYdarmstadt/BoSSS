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
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Utils;

using MPI.Wrappers;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.XRheology_Solver
{
    public class XRheology_SolverMain : BoSSS.Application.XNSE_Solver.XBase_Solver<XRheology_Control>
    {
        static void Main(string[] args) {
            XRheology_SolverMain._Main(args, false,
            delegate () {
                var app = new XRheology_SolverMain();
                return app;
            });
        }

        #region instantiation

        // Attributes for fields (Names), initialization of DG fields
        //==============================================================

        /// <summary>
        /// Velocity and related variables for the extended case, <see cref="XNSE_Control.UseXDG4Velocity"/> == false.
        /// </summary>
        VelocityRelatedVars<XDGField> XDGvelocity;

        /// <summary>
        /// Pressure domain
        /// </summary>
        XDGField Pressure;

        /// <summary>
        /// Pressure codomain: Residuum in continuity equation
        /// </summary>
        XDGField ResidualContinuity;

        /// <summary>
        /// Extra stress domain (2D): StressXX
        /// </summary>
        XDGField StressXX;

        /// <summary>
        /// Extra stress domain (2D): StressXY
        /// </summary>
        XDGField StressXY;

        /// <summary>
        /// Extra stress domain (2D): StressYY
        /// </summary>
        XDGField StressYY;

        /// <summary>
        /// Extra stress codomain (2D): StressXX
        /// </summary>
        XDGField ResidualStressXX;

        /// <summary>
        /// Extra stresses codomain (2D): StressXY
        /// </summary>
        XDGField ResidualStressXY;

        /// <summary>
        /// Extra stresses codomain (2D): StressYY
        /// </summary>
        XDGField ResidualStressYY;

        /// <summary>
        /// Extra stresses parameter (2D): StressXX
        /// </summary>
        XDGField StressXXP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        XDGField StressXYP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        XDGField StressYYP;

        /// <summary>
        /// Extra source (e.g. gravity)
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
                    new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                    true, true,
                    IOListOption.ControlFileDetermined)]
        public VectorField<XDGField> Gravity;

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

        // Parameters: Velocity Gradient
        VectorField<XDGField> VelocityXGradient;
        VectorField<XDGField> VelocityYGradient;


        //Parameters: external analytical velocity
        XDGField U;
        XDGField V;


        // Some initialisation of variables
        //============================================

        /// <summary>
        /// the spatial operator (momentum and continuity equation)
        /// </summary>
        XRheology_OperatorFactory XRheology_Operator;

        /// <summary>
        /// OperatorConfiguration for the <see cref="XRheology_Operator"/>
        /// </summary>
        XRheology_OperatorConfiguration XOpConfig;

        //not sure if needed....
        // ===================================================
        /// <summary>
        /// Current Velocity
        /// </summary>
        XDGField[] CurrentVel {
            get {
                return this.XDGvelocity.Velocity.ToArray();
            }
        }

        XDGField[] prevVel;

        /// <summary>
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        MassMatrixFactory MassFact;
        // ======================================================

        int D; // Spatial Dimension
        /// <summary>
        /// current Weissenberg number
        /// </summary>
        public double currentWeissenberg;
        bool ChangeMesh = true;

        /// <summary>
        /// initialisation of BDF Timestepper
        /// </summary>
        protected XdgBDFTimestepping m_BDF_Timestepper;


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
        protected bool U0MeanRequired {
            get {
                return (this.Control.PhysicalParameters.IncludeConvection);
            }
        }

        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double reynolds_A = this.Control.PhysicalParameters.reynolds_A,
                    reynolds_B = this.Control.PhysicalParameters.reynolds_B;

                int D = this.GridData.SpatialDimension;

                double[] _reynolds_A = new double[D + 1];
                _reynolds_A.SetAll(reynolds_A); // mass matrix in momentum equation
                _reynolds_A[D] = 0; // no  mass matrix for continuity equation
                double[] _reynolds_B = new double[D + 1];
                _reynolds_B.SetAll(reynolds_B); // mass matrix in momentum equation
                _reynolds_B[D] = 0; // no  mass matrix for continuity equation


                //double[] _rho = new double[D + 4];
                //_rho.SetAll(rho);
                ////No MassMatrix for the pressure
                //_rho[D] = 0;

                //_rho[D + 1] = 1;
                //_rho[D + 2] = 1;
                //_rho[D + 3] = 1;

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _reynolds_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), _reynolds_B);

                return R;
            }
        }

        CoordinateVector m_CurrentSolution = null;
        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(XDGvelocity.Velocity, this.Pressure, this.StressXX, this.StressXY, this.StressYY));
                }
                return m_CurrentSolution;
            }
        }

        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(XDGvelocity.ResidualMomentum, this.ResidualContinuity, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY));
                }
                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// DG Field instantiation.
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();
            int D = this.GridData.SpatialDimension;

            // LEVEL SET FIELDS
            DGLevSet = new ScalarFieldHistory<SinglePhaseField>(
                   new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions["Phi"].Degree), "PhiDG"));

            if (this.Control.FieldOptions["PhiDG"].Degree >= 0 && this.Control.FieldOptions["PhiDG"].Degree != this.DGLevSet.Current.Basis.Degree) {
                throw new ApplicationException("Specification of polynomial degree for 'PhiDG' is not supportet, since it is induced by polynomial degree of 'Phi'.");
            }

            this.LsTrk = new LevelSetTracker((GridData)this.GridData, base.Control.CutCellQuadratureType, base.Control.LS_TrackerWidth, new string[] { "A", "B" }, this.LevSet);
            base.RegisterField(this.LevSet);
            this.LevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(this.LevSet.Basis, "dPhi_dx[" + d + "]")));
            base.RegisterField(this.LevSetGradient);

            base.RegisterField(this.DGLevSet.Current);
            this.DGLevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(this.DGLevSet.Current.Basis, "dPhiDG_dx[" + d + "]")));
            base.RegisterField(this.DGLevSetGradient);


            //PRESSURE FIELD
            this.Pressure = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree), VariableNames.Pressure);
            base.RegisterField(this.Pressure);


            // CONTI FIELD
            this.ResidualContinuity = new XDGField(this.Pressure.Basis, "ResidualConti");
            base.RegisterField(this.ResidualContinuity);


            // ALL VELOCITY RELATED FIELDS
            this.XDGvelocity = new VelocityRelatedVars<XDGField>();
            InitFromAttributes.CreateFieldsAuto(this.XDGvelocity, this.GridData, base.Control.FieldOptions, base.Control.CutCellQuadratureType, base.IOFields, base.m_RegisteredFields);


            // ALL STRESS RELATED FIELDS
            this.StressXX = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), VariableNames.StressXX);
            base.RegisterField(this.StressXX);
            this.StressXY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), VariableNames.StressXY);
            base.RegisterField(this.StressXY);
            this.StressYY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), VariableNames.StressYY);
            base.RegisterField(this.StressYY);

            this.ResidualStressXX = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), "ResidualStressXX");
            base.RegisterField(this.ResidualStressXX);
            this.ResidualStressXY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), "ResidualStressXY");
            base.RegisterField(this.ResidualStressXY);
            this.ResidualStressYY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), "ResidualStressYY");
            base.RegisterField(this.ResidualStressYY);

            this.StressXXP = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), VariableNames.StressXXP);
            base.RegisterField(this.StressXXP);
            this.StressXYP = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), VariableNames.StressXYP);
            base.RegisterField(this.StressXYP);
            this.StressYYP = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), VariableNames.StressYYP);
            base.RegisterField(this.StressYYP);


            //PERSSON SENSOR FIELD
            if (Control.UsePerssonSensor == true) {
                perssonsensor = new PerssonSensor(StressXX);
                this.IOFields.Add(perssonsensor.GetField());
            }


            //ARTIFICIAL VISCOSITY FIELD
            if (Control.UseArtificialDiffusion == true) {
                artificalViscosity = new SinglePhaseField(new Basis(GridData, 1), "artificalViscosity");
                this.IOFields.Add(artificalViscosity);

            }
        }


        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
        }

        #endregion

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        {
            throw new NotImplementedException();
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0)
        {
            throw new NotImplementedException();
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt)
        {
            throw new NotImplementedException();
        }
    }
}
