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
using System.Runtime.Serialization;
using System.Linq;
using System.Text;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.NSECommon;

using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.Timestepping;
using static BoSSS.Application.XNSE_Solver.XNSE_Control;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// 
    /// </summary>
    [DataContract]
    [Serializable]
    public abstract class XBase_Control : AppControlSolver {

        /// <summary>
        /// An explicit expression of the Level-set over time.
        /// </summary>
        public Func<double[], double, double> Phi;

        /// <summary>
        /// Options for the initialization of the Fourier Level-set
        /// </summary>
        [DataMember]
        public FourierLevSetControl FourierLevSetControl;

        /// <summary>
        /// Enforce the level-set to be globally conervativ, by adding a constant to the level-set field
        /// </summary>
        public bool EnforceLevelSetConservation = false;

        /// <summary>
        /// Control Options for ExtVel
        /// </summary>
        public EllipticExtVelAlgoControl EllipticExtVelAlgoControl = new EllipticExtVelAlgoControl();

        /// <summary>
        /// See <see cref="ContinuityProjection"/>
        /// </summary>
        [DataMember]
        public ContinuityProjectionOption LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;

        /// <summary>
        /// Width of the narrow band.
        /// </summary>
        [DataMember]
        public int LS_TrackerWidth = 1;

        /// <summary>
        /// See <see cref="LevelSetEvolution"/>.
        /// </summary>
        [DataMember]
        public LevelSetEvolution Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

        /// <summary>
        /// if true, the jump condition for mass, momentum and energy will be checked
        /// </summary>
        [DataMember]
        public bool CheckJumpConditions = false;

        /// <summary>
        /// if true, the mass conservation and the surface changerate is checked
        /// </summary>
        [DataMember]
        public bool CheckInterfaceProps = false;

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public int ReInitPeriod = 0;

        /// <summary>
        /// Expert options regarding the spatial discretization.
        /// </summary>
        [DataMember]
        public DoNotTouchParameters AdvancedDiscretizationOptions = new DoNotTouchParameters() {
            PenaltySafety = 1.0,
            alpha = 1.0,
            ObjectiveParam = 1.0,
            //Penalty1 = { 0, 0 },
            Penalty2 = 1.0,
            //PresPenalty1 = {0,0},
            PresPenalty2 = 1.0,
            StressPenalty = 1.0,
            ViscosityMode = ViscosityMode.Viscoelastic

        };

        /// <summary>
        /// Viscosity, density and surface tension.
        /// </summary>
        [DataMember]
        public PhysicalParameters PhysicalParameters = new PhysicalParameters() {
            Material = true,
            IncludeConvection = false,
            reynolds_A = 1.0,
            reynolds_B = 1.0,
            mu_A = 1.0,
            mu_B = 1.0,
            rho_A = 1.0,
            rho_B = 1.0,
            Sigma = 0.0,
            Weissenberg_a = 0.0,
            Weissenberg_b = 0.0,
            beta_a = 0.0,
            beta_b = 0.0
        };

        /// <summary>
        /// See <see cref="InterfaceAveraging"/>
        /// </summary>
        public InterfaceAveraging InterAverage = InterfaceAveraging.density;

        /// <summary>
        /// average method for interface values
        /// </summary>
        public enum InterfaceAveraging {

            /// <summary>
            /// arithmetic mean
            /// </summary>
            mean,

            /// <summary>
            /// density weighted average
            /// </summary>
            density,

            /// <summary>
            /// viscosity weighted average
            /// </summary>
            viscosity

        }



        /// <summary>
        /// See <see cref="TimesteppingScheme"/>
        /// </summary>
        [DataMember]
        public TimesteppingScheme Timestepper_Scheme = TimesteppingScheme.ImplicitEuler;


        /*
        /// <summary>
        /// Timestepping schemes for the XdgTimestepper
        /// Either implicit timestepping using Backward-Differentiation-Formulas (BDF) formulas by <see cref="XdgBDFTimestepping"/> 
        /// or explicit/implicit using Runge-Kutta schemes <see cref="XdgRKTimestepping"/>
        /// </summary>
        public enum TimesteppingScheme {

            ImplicitEuler = 1,

            CrankNicolson = 2,

            BDF2 = 3,

            BDF3 = 4,

            BDF4 = 5,

            BDF5 = 6,

            BDF6 = 7,

            RK_ImplicitEuler = 201,

            RK_CrankNicolson = 202
        }
        */








        /// <summary>
        /// Data to be written in LogFile
        /// </summary>
        public enum LoggingValues {

            /// <summary>
            /// no data will be written
            /// </summary>
            None,

            /// <summary>
            /// for elemental test programm with line like interfaces
            /// </summary>
            LinelikeLS,

            /// <summary>
            /// for elemental test programm with circle like interfaces
            /// </summary>
            CirclelikeLS,

            /// <summary>
            /// for wavelike simulation as CapillaryWave, RT-Instability
            /// interface height (interface points)
            /// </summary>
            Wavelike,

            /// <summary>
            /// for the benchmark quantities of the Rising Bubble testcase
            /// </summary>
            RisingBubble,

            /// <summary>
            /// contact points and corresponding contact angle
            /// </summary>
            MovingContactLine,

            /// <summary>
            /// height of a rising capillary in a tube
            /// </summary>
            CapillaryHeight,

            /// <summary>
            /// Evaporative mass flux and speed of displacement 
            /// </summary>
            Evaporation,

            /// <summary>
            ///  half-radii and drop deformation for an elliptic shapes drop
            /// </summary>
            DropDeformation
        }

        /// <summary>
        /// See <see cref="LoggingValues"/>
        /// </summary>
        [DataMember]
        public LoggingValues LogValues = LoggingValues.None;

        [DataMember]
        public int LogPeriod = 1;

      


    }
}
