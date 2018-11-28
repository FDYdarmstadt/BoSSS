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
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;



namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class FSI_Control : IBM_Solver.IBM_Control {

        /// <summary>
        /// Set true if translation of the particle should be induced by hydrodynamical forces.
        /// </summary>
        [DataMember]
        public bool includeTranslation = false;

        /// <summary>
        /// Set true if rotation of the particle should be indruced by hydrodynamical torque.
        /// </summary>
        [DataMember]
        public bool includeRotation = false;

        /// <summary>
        /// Function describing the fixed level-set movement
        /// </summary>
        [DataMember]
        public Func<double[], double, double> MovementFunc;

        /// <summary>
        /// Function describing the fixed level-set movement
        /// </summary>
        [DataMember]
        public Func<double, double>[] transVelocityFunc;

        /// <summary>
        /// Function describing the fixed level-set movement
        /// </summary>
        [DataMember]
        public Func<double, double>[] anglVelocityFunc;

        /// <summary>
        /// How should the level set be moved? Options: none, fixed, coupled
        /// </summary>
        [DataMember]
        public string LevelSetMovement = "none";

        /// <summary>
        /// 
        /// </summary>
        public enum TimesteppingMode {

            None = 0,

            /// <summary>
            /// 
            /// </summary>
            Splitting = 1,

            /// <summary>
            /// 
            /// </summary>
            MovingMesh = 2
        }
        [DataMember]
        public TimesteppingMode Timestepper_Mode = TimesteppingMode.Splitting;
        /*
        /// <summary>
        /// Function describing the boundary values at the level-set (VelocityX, VelocityY)
        /// </summary>
        public Func<double, double>[] BoundaryFunc;
<<<<<<< HEAD:public/src/L4-application/FSI_Solver/FSI_Control.cs
        */
=======

>>>>>>> 15057d440a0a91354d4b17eb02ab1e9b6945fd82:internal/src/experimental/L4-application/FSI_Solver/FSI_Control.cs
        [DataMember]
        public List<Particle> Particles;
       
        public enum CollisionModel {
            RepulsiveForce = 0,

            MomentumConservation = 1,

            NoCollisionModel = 2

        }
<<<<<<< HEAD:public/src/L4-application/FSI_Solver/FSI_Control.cs

        [DataMember]
        public CollisionModel collisionModel = CollisionModel.MomentumConservation;

        //public double particleMass;

        //public double particleRho;

=======
        [DataMember]
        public CollisionModel collisionModel = CollisionModel.MomentumConservation;

>>>>>>> 15057d440a0a91354d4b17eb02ab1e9b6945fd82:internal/src/experimental/L4-application/FSI_Solver/FSI_Control.cs
        [DataMember]
        public bool pureDryCollisions = false;

        /// <summary>
        /// Adds particle to particle list
        /// </summary>
        /// <param name="D"></param>
        /// <param name="HistoryLength"></param>
        /// <param name="start"></param>
        public void AddParticle(int D, int HistoryLength, double[] start) {
            this.Particles.Add(new Particle(D, HistoryLength, start));
        }


    }
}
