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
using System.Runtime.Serialization;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using System.Collections;
using FSI_Solver;

namespace BoSSS.Application.FSI_Solver {

    /// <summary>
    /// Particle properties
    /// </summary>
    [DataContract]
    [Serializable]
    abstract public class Particle : ICloneable {

        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        protected Particle() {
            // noop
        }

        /// <summary>
        /// Constructor for an arbitrary particle to be implemented in the Particle_Shape classes.
        /// </summary>
        /// <param name="motionInit">
        /// Initializes the motion parameters of the particle (which model to use, whether it is a dry simulation etc.)
        /// </param>
        /// <param name="startPos">
        /// The initial position.
        /// </param>
        /// <param name="startAngl">
        /// The inital anlge.
        /// </param>
        /// <param name="activeStress">
        /// The active stress excerted on the fluid by the particle. Zero for passive particles.
        /// </param>
        /// <param name="startTransVelocity">
        /// The inital translational velocity.
        /// </param>
        /// <param name="startRotVelocity">
        /// The inital rotational velocity.
        /// </param>
        public Particle(ParticleMotionInit motionInit, double[] startPos, double startAngl = 0.0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) {
            SpatialDim = startPos.Length;
            ActiveStress = activeStress;
            Aux = new FSI_Auxillary();

            motionInit.CheckInput();
            Motion = motionInit.ParticleMotion;
            Motion.InitializeParticlePositionAndAngle(startPos, startAngl);
            Motion.InitializeParticleVelocity(startTransVelocity, startRotVelocity);
            particleDensity = Motion.Density;
        }

        [NonSerialized]
        protected readonly FSI_Auxillary Aux;
        [DataMember]
        private readonly double particleDensity;

        [DataMember]
        protected int SpatialDim { get; }

        /// <summary>
        /// Instantiate object for particle motion.
        /// </summary>
        [DataMember]
        public Motion Motion { get; private set; }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        [DataMember]
        protected double Mass_P => Area * particleDensity;

        /// <summary>
        /// Check whether any particles is collided with another particle
        /// </summary>
        [DataMember]
        public bool IsCollided { get; set; }

        /// <summary>
        /// No of convex (!) sub-particles (for GJK algorithm, distance calc)
        /// </summary>
        [DataMember]
        public virtual int NoOfSubParticles => 1;
        
        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double Eccentricity { get; private set; }

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public Vector ClosestPointToOtherObject = new Vector(2);

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public Vector ClosestPointOnOtherObjectToThis = new Vector(2);

        /// <summary>
        /// Active stress on the current particle.
        /// </summary>
        [DataMember]
        public double ActiveStress { get; private set; } = 0;

        /// <summary>
        /// Area of the current particle.
        /// </summary>
        [DataMember]
        public virtual double Area => throw new NotImplementedException();

        /// <summary>
        /// Level set function describing the particle.
        /// </summary>   
        public abstract double LevelSetFunction(double[] X);

        /// <summary>
        /// Circumference of the current particle.
        /// </summary>
        public abstract double Circumference { get; }

        /// <summary>
        /// Moment of inertia of the current particle.
        /// </summary>
        public abstract double MomentOfInertia { get; }

        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        public CellMask CutCells_P(LevelSetTracker LsTrk) {
            BitArray CellArray = new BitArray(LsTrk.GridDat.Cells.NoOfLocalUpdatedCells);
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            double h_min = LsTrk.GridDat.Cells.h_minGlobal;
            for (int i = 0; i < CellArray.Length; i++) {
                CellArray[i] = Contains(new Vector(CellCenters[i, 0], CellCenters[i, 1]), h_min);
            }
            CellMask CutCells = new CellMask(LsTrk.GridDat, CellArray, MaskType.Logical);
            return CutCells;
        }

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
        /// </summary>
        /// <param name="point"></param>
        public virtual bool Contains(Vector point, double tolerance = 0) => throw new NotImplementedException();

        /// <summary>
        /// Return the lengthscales of the particle (length and thickness)
        /// </summary>
        public virtual double[] GetLengthScales() => throw new NotImplementedException();

        /// <summary>
        /// Returns surface points (for distance calc)
        /// </summary>
        /// <param name="hMin"></param>
        public virtual MultidimensionalArray GetSurfacePoints(double hMin, double searchAngle, int subParticleID) => throw new NotImplementedException();

        /// <summary>
        /// Calculates the support point with an analytic formula (if applicable)
        /// </summary>
        /// <param name="supportVector"></param>
        public virtual Vector GetSupportPoint(Vector supportVector, int SubParticleID) => throw new NotImplementedException();

        /// <summary>
        /// Calculates the radial vector (SurfacePoint-ParticleReadOnlyPosition)
        /// </summary>
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialVector">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        internal void CalculateRadialVector(Vector SurfacePoint, out Vector RadialVector, out double RadialLength) {
            RadialVector = new Vector(SurfacePoint[0] - Motion.GetPosition(0)[0], SurfacePoint[1] - Motion.GetPosition(0)[1]);
            if (RadialVector.L2Norm() == 0)
                throw new ArithmeticException("The radial vector has no length");
            RadialLength = RadialVector.Abs();
            RadialVector.ScaleV(1 / RadialLength);
            Aux.TestArithmeticException(RadialVector, "particle radial vector");
            Aux.TestArithmeticException(RadialLength, "particle radial length");
        }

        /// <summary>
        /// Calculates the eccentricity of a collision
        /// </summary>
        internal void CalculateEccentricity() {
            CalculateRadialVector(ClosestPointToOtherObject, out Vector RadialVector, out _);
            Eccentricity = RadialVector * Motion.GetLastCollisionTangentialVector();
            Aux.TestArithmeticException(Eccentricity, "particle eccentricity");
        }

        /// <summary>
        /// clone, not implemented
        /// </summary>
        public virtual object Clone() => throw new NotImplementedException("Currently cloning of a particle is not available");
    }
}

