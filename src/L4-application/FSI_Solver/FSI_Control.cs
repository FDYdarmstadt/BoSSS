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
        /// ctor
        /// </summary>
        public FSI_Control() : base() {
            this.Particles = new List<Particle>();
        }
        
        /// <summary>
        /// Set true if the coupling between fluid and particle should be calculated iterative, while using Lie-Splitting.
        /// </summary>
        public bool splitting_fully_coupled()
        {
            return Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled;
        }

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
        public bool LowerWallFullyPlastic = false;

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
        /// The termination criterion for fully coupled/implicit level-set evolution.
        /// </summary>
        [DataMember]
        public double ForceAndTorque_ConvergenceCriterion = 1.0e-6;

        /// <summary>
        /// under-relaxation of the level set movement in case of coupled iterative
        /// </summary>
        [DataMember]
        public double LSunderrelax = 1.0;

     

        /// <summary>
        /// Setting <see cref="Solution.Control.AppControl.FieldOptions"/>
        /// </summary>
        public override void SetDGdegree(int k) {
            if (k < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");
            int k_phiDG = Math.Max(2, k);
            int k_phi = Convert.ToInt32(Math.Pow(k_phiDG, 2));
            base.FieldOptions.Clear();
            this.AddFieldOption("Velocity*", k);
            this.AddFieldOption("Pressure", k - 1);
            this.AddFieldOption("PhiDG", k_phiDG);
            this.AddFieldOption("Phi", k_phi);
            this.AddFieldOption("Curvature", 2);
        }
        

        /// <summary>
        /// See <see cref="LevelSetHandling"/>
        /// </summary>
        [DataMember]
        public LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


        /// <summary>
        /// Function describing the boundary values at the level-set (VelocityX, VelocityY)
        /// </summary>
        public Func<double, double>[] BoundaryFunc;

        /// <summary>
        /// All particles in the FSI
        /// </summary>
        [DataMember]
        public List<Particle> Particles {
            get;
            set;
        }

        public enum CollisionModel {
            RepulsiveForce = 0,

            MomentumConservation = 1,

            NoCollisionModel = 2

        }

        [DataMember]
        public CollisionModel collisionModel = CollisionModel.MomentumConservation;

        /// <summary>
        /// if true the flow solver is turned off
        /// </summary>
        [DataMember]
        public bool pureDryCollisions = false;

        public override Type GetSolverType() {
            return typeof(FSI_SolverMain);
        }
    }
}
