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
using System.Text;
using System.Threading.Tasks;

using ilPSP.Utils;

using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.EnergyCommon;


namespace BoSSS.Application.XNSE_Solver {

    public class XNSE_OperatorConfiguration : IEnergy_Configuration {

        public XNSE_OperatorConfiguration() {}

        public XNSE_OperatorConfiguration(XNSE_Control control) {


            Continuity = true;
            Viscous = true;
            PressureGradient = true;
            Transport = control.PhysicalParameters.IncludeConvection;
            CodBlocks = new bool[] { true, true };
            DomBlocks = new bool[] { true, true };
            dntParams = control.AdvancedDiscretizationOptions;
            physParams = control.PhysicalParameters;


            if(control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SemiImplicit)
                control.PhysicalParameters.mu_I = control.dtFixed * control.PhysicalParameters.Sigma;
            

            switch(control.Timestepper_LevelSetHandling) {
                case LevelSetHandling.Coupled_Once:
                    MovingMesh = true;
                    mmsd = MassMatrixShapeandDependence.IsTimeDependent;
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    MovingMesh = true;
                    mmsd = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    break;

                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.StrangSplitting:
                    MovingMesh = false;
                    mmsd = MassMatrixShapeandDependence.IsTimeDependent;
                    break;

                case LevelSetHandling.None:
                    MovingMesh = false;
                    mmsd = MassMatrixShapeandDependence.IsNonIdentity;
                    break;

                default:
                    throw new NotImplementedException();
            }


            kinEviscous = control.kinEViscousDiscretization;
            kinEpressure = control.kinEPressureDiscretization;
            withDissP = control.withDissipativePressure;

        }


        public PhysicalParameters physParams;

        //public ThermalParameters thermParams;

        /// <summary>
        /// advanced operator configuration
        /// </summary>
        public DoNotTouchParameters dntParams;

        /// <summary>
        /// Fucking Gravity
        /// </summary>
        public bool isGravity {
            get {
                return true;
            }
        }

        /// <summary>
        /// Controls the domain variables that the operator should contain. <br/>
        /// This controls only the formal operator shape, not the actual components.
        /// index: 0 -- Velocity components, 1 -- pressure component; 
        /// </summary>
        public bool[] DomBlocks = new bool[2];

        /// <summary>
        /// Controls the codomain variables that the operator should contain. <br/>
        ///This controls only the formal operator shape, not the actual components.
        /// index: 0 -- momentum equation, 1 -- continuity equation; 
        /// </summary>
        public bool[] CodBlocks = new bool[2];

        /// <summary>
        /// include transport operator
        /// </summary>
        public bool Transport;

        /// <summary>
        /// include viscous operator
        /// </summary>
        public bool Viscous;

        /// <summary>
        /// include pressure gradient
        /// </summary>
        public bool PressureGradient;

        /// <summary>
        /// Use the surface force term -- mostly only useful for manufactured solutions.
        /// </summary>
        public bool ArtificialsurfaceForce = false;

        /// <summary>
        /// include continuity equation
        /// </summary>
        public bool Continuity;

        // <summary>
        // include heat equation
        // </summary>
        //public bool Heat;

        // <summary>
        // include evaporation
        // </summary>
        //public bool Evaporation;

        /// <summary>
        /// true if the interface is a material interface
        /// </summary>
        public bool MatInt = true;

        /// <summary>
        /// Switch to turn velocity extension on/off.
        /// </summary>
        public bool UseXDG4Velocity = true;

        /// <summary>
        /// switch for moving mesh flux discretizations
        /// </summary>
        public bool MovingMesh;

        /// <summary>
        /// 
        /// </summary>
        public MassMatrixShapeandDependence mmsd;


        /// <summary>
        /// option for the discretization of the viscous kinetic energy source terms
        /// </summary>
        public KineticEnergyViscousSourceTerms kinEviscous;

        /// <summary>
        /// option for the discretization of the pressure kinetic energy source terms
        /// </summary>
        public KineticEnergyPressureSourceTerms kinEpressure;

        /// <summary>
        /// adds a locally discretized dissipation term regarding pressure
        /// </summary>
        public bool withDissP;


        // getter for interface
        // ====================

        public PhysicalParameters getPhysParams {
            get { return physParams; }
        }

        //public ThermalParameters getThermParams {
        //    get { return thermParams; }
        //}

        public DoNotTouchParameters getDntParams {
            get { return dntParams; }
        }


        public bool[] getDomBlocks {
            get { return DomBlocks; }
        }

        public bool[] getCodBlocks {
            get { return CodBlocks; }
        }

        public bool isTransport {
            get { return Transport; }
        }

        public bool isViscous {
            get { return Viscous; }
            set { Viscous = value; }
        }

        public bool isPressureGradient {
            get { return PressureGradient; }
        }

        public bool isContinuity {
            get { return Continuity; }
        }

        public bool isMovingMesh {
            get { return MovingMesh; }
        }

        public bool isMatInt {
            get { return MatInt; }
        }

        public virtual bool isPInterfaceSet {
            get { return false; }
        }


        public virtual KineticEnergyViscousSourceTerms getKinEviscousDiscretization {
            get { return kinEviscous; }
        }

        public virtual KineticEnergyPressureSourceTerms getKinEpressureDiscretization {
            get { return kinEpressure; }
        }

        public virtual bool withPressureDissipation {
            get { return withDissP; }
        }

    }


    public class XNSFE_OperatorConfiguration : XNSE_OperatorConfiguration, IXHeat_Configuration {


        public XNSFE_OperatorConfiguration(XNSE_Control control) 
            : base(control) {

            AgglomerationTreshold = control.AdvancedDiscretizationOptions.CellAgglomerationThreshold;

            thermParams = control.ThermalParameters;

            solveEnergy = control.solveKineticEnergyEquation;

            solveHeat = control.solveCoupledHeatEquation;
            Evaporation = (control.ThermalParameters.hVap > 0.0);
            if(control.prescribedMassflux_Evaluator != null)
                prescribedMassflux = control.prescribedMassflux_Evaluator;
            if (control.prescribedMassflux != null)
                prescribedMassflux = control.prescribedMassflux.Evaluate;
            MatInt = !Evaporation;

            int nBlocks = 2;
            //if (solveEnergy) {
            //    CodBlocks = new bool[3];
            //    DomBlocks = new bool[3];
            //    nBlocks = 3;
            //}

            if (solveHeat) {
                this.conductMode = control.conductMode;
                CodBlocks = (this.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) ? new bool[nBlocks + 1] : new bool[nBlocks + 2];
                DomBlocks = (this.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) ? new bool[nBlocks + 1] : new bool[nBlocks + 2];
            }

            CodBlocks.SetAll(true);
            DomBlocks.SetAll(true);
        }

        /// <summary>
        /// taken from <see cref="DoNotTouchParameters.CellAgglomerationThreshold"/>
        /// </summary>
        public double AgglomerationTreshold;

        /// <summary>
        /// 
        /// </summary>
        public ThermalParameters thermParams;

        /// <summary>
        /// include kinetic energy equation
        /// </summary>
        public bool solveEnergy;

        /// <summary>
        /// true if the heat equation is solved via the auxiliary heat flux formulation
        /// </summary>
        public ConductivityInSpeciesBulk.ConductivityMode conductMode;

        /// <summary>
        /// include heat equation
        /// </summary>
        public bool solveHeat;

        /// <summary>
        /// include transport operator
        /// </summary>
        public bool HeatTransport;

        /// <summary>
        /// use upwind discretization
        /// </summary>
        public bool HeatUpwinding = false;

        /// <summary>
        /// include evaporation
        /// </summary>
        public bool Evaporation;

        /// <summary>
        /// 
        /// </summary>
        public Func<double[], double, double> prescribedMassflux;

        public ThermalParameters getThermParams {
            get { return thermParams; }
        }

        public ConductivityInSpeciesBulk.ConductivityMode getConductMode {
            get { return conductMode; }
        }
        
        public bool isHeatTransport {
            get { return HeatTransport; }
        }

        public bool useUpwind {
            get { return HeatUpwinding;  }
        }

        public bool isEvaporation {
            get { return Evaporation; }
        }
    }
}
