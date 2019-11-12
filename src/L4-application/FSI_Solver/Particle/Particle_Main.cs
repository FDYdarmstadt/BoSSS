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
using System.Linq;
using System.Runtime.Serialization;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using System.Collections;
using FSI_Solver;
using System.Collections.Generic;

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
            m_MotionInit = motionInit;
            Aux = new FSI_Auxillary();

            m_MotionInit.CheckInput();
            Motion = m_MotionInit.ParticleMotion;
            Motion.InitializeParticlePositionAndAngle(startPos, startAngl);
            Motion.InitializeParticleVelocity(startTransVelocity, startRotVelocity);
            particleDensity = Motion.Density;
        }

        [NonSerialized]
        protected readonly FSI_Auxillary Aux;
        [DataMember]
        private readonly ParticleMotionInit m_MotionInit;
        [DataMember]
        private readonly double particleDensity;

        [DataMember]
        protected int SpatialDim { get; }

        /// <summary>
        /// Instantiate object for particle motion.
        /// </summary>
        [DataMember]
        public Motion_Wet Motion { get; private set; } = new Motion_Wet(gravity: new double[] { 0, 9.81 }, density: 1);

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
        public double[] ClosestPointToOtherObject { get; set; }

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double[] ClosestPointOnOtherObjectToThis { get; set; }

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
        /// Necessary for active particles. Returns 0 for the non active boundary region and a number between 0 and 1 for the active region.
        /// </summary>
        public double SeperateBoundaryRegions(double[] X) => Math.Cos(Motion.GetAngle(0)) * (X[0] - Motion.GetPosition(0)[0]) + Math.Sin(Motion.GetAngle(0)) * (X[1] - Motion.GetPosition(0)[1]) < 1e-8
            ? (Math.Cos(Motion.GetAngle(0)) * (X[0] - Motion.GetPosition(0)[0]) + Math.Sin(Motion.GetAngle(0)) * (X[1] - Motion.GetPosition(0)[1])) / Math.Sqrt((X[0] - Motion.GetPosition(0)[0]).Pow2() + (X[1] - Motion.GetPosition(0)[1]).Pow2())
            : 0;

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
            double h_max = LsTrk.GridDat.Cells.h_maxGlobal;

            for (int i = 0; i < CellArray.Length; i++) {
                CellArray[i] = Contains(new double[] { CellCenters[i, 0], CellCenters[i, 1] }, h_min, h_max, false);
            }
            CellMask CutCells = new CellMask(LsTrk.GridDat, CellArray, MaskType.Logical);
            CutCells = CutCells.Intersect(LsTrk.Regions.GetCutCellMask());
            return CutCells;
        }

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
        /// </summary>
        /// <param name="point"></param>
        /// <param name="h_min"></param>
        /// <param name="h_max"></param>
        /// <param name="WithoutTolerance"></param>
        public virtual bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) => throw new NotImplementedException();

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
        /// <param name="Vector"></param>
        public virtual double[] GetSupportPoint(double[] Vector, int SubParticleID) => throw new NotImplementedException();

        /// <summary>
        /// Calculates the radial vector (SurfacePoint-ParticleReadOnlyPosition)
        /// </summary>
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialVector">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        internal void CalculateRadialVector(double[] SurfacePoint, out double[] RadialVector, out double RadialLength) {
            RadialVector = new double[] { SurfacePoint[0] - Motion.GetPosition(0)[0], SurfacePoint[1] - Motion.GetPosition(0)[1] };
            if (RadialVector.L2Norm() == 0)
                throw new ArithmeticException("The given vector has no length");
            RadialLength = RadialVector.L2Norm();
            RadialVector.ScaleV(1 / RadialLength);
            Aux.TestArithmeticException(RadialVector, "particle radial vector");
            Aux.TestArithmeticException(RadialLength, "particle radial length");
        }

        /// <summary>
        /// Calculates the vector normal to the radial vector.
        /// </summary>
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialNormalVector">
        /// </param>
        internal void CalculateRadialNormalVector(double[] SurfacePoint, out double[] RadialNormalVector) {
            RadialNormalVector = new double[] { SurfacePoint[1] - Motion.GetPosition(0)[1], -SurfacePoint[0] + Motion.GetPosition(0)[0] };
            RadialNormalVector.ScaleV(1 / RadialNormalVector.L2Norm());
            Aux.TestArithmeticException(RadialNormalVector, "particle vector normal to radial vector");
        }

        /// <summary>
        /// Calculates the eccentricity of a collision
        /// </summary>
        internal void CalculateEccentricity() {
            CalculateRadialVector(ClosestPointToOtherObject, out double[] RadialVector, out _);
            Eccentricity = RadialVector[0] * Motion.GetLastCollisionTangentialVector()[0] + RadialVector[1] * Motion.GetLastCollisionTangentialVector()[1];
            Aux.TestArithmeticException(Eccentricity, "particle eccentricity");
        }

        /// <summary>
        /// clone, not implemented
        /// </summary>
        public virtual object Clone() => throw new NotImplementedException("Currently cloning of a particle is not available");
    }
}

