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

using BoSSS.Solution.XdgTimestepping;
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
        /// Set true if the coupling between fluid and particle should be calculated iterative, while using Lie-Splitting.
        /// </summary>
        [DataMember]
        public bool splitting_fully_coupled = false;

        /// <summary>
        /// Set true if the coupling between fluid and particle should be calculated iterative, while using Lie-Splitting.
        /// </summary>
        [DataMember]
        public int max_iterations_fully_coupled = 10000;

        /// <summary>
        /// Set true if translation of the particle should be induced by hydrodynamical forces.
        /// </summary>
        [DataMember]
        public bool instationarySolver = true;

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
        /// See <see cref="LevelSetHandling"/>
        /// </summary>
        [DataMember]
        public LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

        /// <summary>
        /// The termination criterion for fully coupled/implicit level-set evolution.
        /// </summary>
        [DataMember]
        public double ForceAndTorque_ConvergenceCriterion = 1.0e-6;

        /// <summary>
        /// underrelaxation of the level set movement in case of coupled iterative
        /// </summary>
        public double LSunderrelax = 1.0;

        /// <summary>
        /// desired minimum refinement level, 2 is minimum
        /// </summary>
        [DataMember]
        public int RefinementLevel = 2;


        /// <summary>
        /// reciprocal of the ratio between curvature and hmin
        /// </summary>
        [DataMember]
        public int maxCurvature = 2;

        public override void SetDGdegree(int k)
        {
            if (k < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");

            base.FieldOptions.Clear();
            this.AddFieldOption("Velocity*", k);
            this.AddFieldOption("Pressure", k - 1);
            this.AddFieldOption("PhiDG", k);
            this.AddFieldOption("Phi", k);
            this.AddFieldOption("Curvature", k);
        }

        ///// <summary>
        ///// How should the level set be moved? Options: none, fixed, coupled
        ///// </summary>
        //public string LevelSetMovement = "none";

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

        /// <summary>
        /// Function describing the boundary values at the level-set (VelocityX, VelocityY)
        /// </summary>
        public Func<double, double>[] BoundaryFunc;

        [DataMember]
        public List<Particle> Particles;
       
        public enum CollisionModel {
            RepulsiveForce = 0,

            MomentumConservation = 1,

            NoCollisionModel = 2

        }

        [DataMember]
        public CollisionModel collisionModel = CollisionModel.MomentumConservation;

        //public double particleMass;

        //public double particleRho;

        [DataMember]
        public bool pureDryCollisions = false;

        public override Type GetSolverType() {
            return typeof(FSI_SolverMain);
        }
    }
}
