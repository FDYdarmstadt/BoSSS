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

using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;


namespace BoSSS.Application.XNSE_Solver {

    public class XNSE_OperatorConfiguration : IXNSE_Configuration {

        public XNSE_OperatorConfiguration() {}

        public XNSE_OperatorConfiguration(XNSE_Control control) {

            if(control.FakePoisson) {
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("ACHTUNG: Fake-Poisson aktiviert!");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            }

            Continuity = true;
            Viscous = !control.FakePoisson;
            PressureGradient = true;
            Transport = !control.FakePoisson;
            CodBlocks = new bool[] { true, true };
            DomBlocks = new bool[] { true, true };
            dntParams = control.AdvancedDiscretizationOptions;
            physParams = control.PhysicalParameters;
            thermParams = control.ThermalParameters;
            UseXDG4Velocity = control.UseXDG4Velocity;

            if(control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SemiImplicit)
                control.PhysicalParameters.mu_I = control.dtFixed * control.PhysicalParameters.Sigma;

            Heat = control.solveCoupledHeatEquation;
            Evaporation = (control.ThermalParameters.hVap_A != 0.0 && control.ThermalParameters.hVap_B != 0.0);
            MatInt = !Evaporation;
            

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
        }


        public PhysicalParameters physParams;

        public ThermalParameters thermParams;

        /// <summary>
        /// advanced operator configuration
        /// </summary>
        public DoNotTouchParameters dntParams;

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

        /// <summary>
        /// include heat equation
        /// </summary>
        public bool Heat;

        /// <summary>
        /// include evaporation
        /// </summary>
        public bool Evaporation;

        /// <summary>
        /// true if the interface is a material interface
        /// </summary>
        public bool MatInt;

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


        // getter for interface
        // ====================

        public PhysicalParameters getPhysParams {
            get { return physParams; }
        }

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
    }

       
}
