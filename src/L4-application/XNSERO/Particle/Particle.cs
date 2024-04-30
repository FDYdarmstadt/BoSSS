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
using BoSSS.Foundation.Grid;
using System.Collections;
using System.Linq;
using ilPSP.Tracing;
using System.Text.Json.Serialization;

namespace BoSSS.Application.XNSERO_Solver {

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
        /// <param name="motion">
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
        public Particle(IMotion motion, double[] startPos, double startAngl = 0.0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) {
            SpatialDim = startPos.Length;
            ActiveStress = activeStress;
            Aux = new Auxillary();
            this.Motion = motion.CloneAs() ?? throw new ArgumentNullException("Missing definition of particle motion");
            this.Motion.InitializeParticlePositionAndAngle(startPos, startAngl);
            this.Motion.InitializeParticleVelocity(startTransVelocity, startRotVelocity);
            Density = this.Motion.Density;
        }

        [DataMember]
        public bool IncludeRotation = true;
        [DataMember]
        public bool IncludeTranslation = true;

        [DataMember]
        public double phoreticActivity = 0;

        /// <summary>
        /// Level set tracker of the solid level set.
        /// </summary>
        [JsonIgnore]
        public LevelSetTracker LsTrk { get; set; }
                
        /// <summary>
        /// Provides additional methods for testing and checking variables.
        /// </summary>
        [NonSerialized]
        protected Auxillary Aux;

        /// <summary>
        /// The density of the particle.
        /// </summary>
        [DataMember]
        private readonly double Density;

        /// <summary>
        /// The spatial dimension.
        /// </summary>
        [DataMember]
        protected int SpatialDim;

        /// <summary>
        /// Instantiate object for particle motion.
        /// </summary>
        [DataMember]
        public IMotion Motion { get; private set; }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        [DataMember]
        public double Mass => Volume * Density;

        /// <summary>
        /// Check whether any particles is collided with another particle
        /// </summary>
        [DataMember]
        public bool IsCollided { get; set; }

        /// <summary>
        /// No of convex (!) sub-particles (for GJK algorithm, distance calculation)
        /// </summary>
        [DataMember]
        public virtual int NoOfSubParticles => 1;
        
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
        /// Area of the current particle. Implementation within the specific particle shapes.
        /// </summary>
        [DataMember]
        public virtual double Volume => throw new NotImplementedException();
        
        /// <summary>
        /// Circumference of the current particle.
        /// </summary>
        public abstract double Circumference { get; }

        /// <summary>
        /// Moment of inertia of the current particle.
        /// </summary>
        public abstract double MomentOfInertia { get; }

        /// <summary>
        /// Level set function describing the particle. Adds also particle level set over periodic boundaries.
        /// </summary>   
        public double LevelSetFunction(double[] X, double GridLength) {
            double levelSet = ParticleLevelSetFunction(X, Motion.GetPosition());
            for (int i = 0; i < Motion.OriginInVirtualPeriodicDomain.Count; i++) {
                Vector virtualPosition = Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition();
                if (Motion.IsInsideOfPeriodicDomain(virtualPosition, (GridLength * 2 + GetLengthScales().Max())))
                    levelSet = Math.Max(levelSet, ParticleLevelSetFunction(X, Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition()));
            }
            return levelSet;
        }

        /// <summary>
        /// The level set function of the particle, depending on the shape of the particle.
        /// </summary>
        /// <param name="X">A point within the computational domain.</param>
        /// <param name="Postion">The position of the particle.</param>
        /// <returns></returns>
        protected virtual double ParticleLevelSetFunction(double[] X, Vector Postion) => throw new NotImplementedException();

        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        public CellMask ParticleCutCells(LevelSetTracker LsTrk, CellMask AllCutCells) {
            BitArray CellArray = AllCutCells.GetBitMask();
            BitArray ContainArray = new(CellArray.Length);
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            double h = 1;// LsTrk.GridDat.Cells.h_maxGlobal;
            for (int i = 0; i < CellArray.Length; i++) {
                if (CellArray[i]) {
                    ContainArray[i] = Contains(new Vector(CellCenters[i, 0], CellCenters[i, 1]),h);
                }
            }
            CellMask CutCells = new(LsTrk.GridDat, ContainArray, MaskType.Logical);
            return CutCells;
        }

        /// <summary>
        /// Checks whether a point is inside of a particle. Checks also for particles at periodic boundaries.
        /// </summary>
        /// <param name="Point">The point to be checked.</param>
        /// <param name="Tolerance">A tolerance parameter</param>
        /// <returns></returns>
        public bool Contains(Vector Point, double Tolerance = 0) {
            bool contains = ParticleContains(Point, Motion.GetPosition(), Tolerance);
            if (!contains) {
                for (int i = 0; i < Motion.OriginInVirtualPeriodicDomain.Count; i++) {
                    Vector virtualPosition = Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition();
                    if (Motion.IsInsideOfPeriodicDomain(virtualPosition, Tolerance + GetLengthScales().Max()))
                        contains = ParticleContains(Point, virtualPosition, Tolerance);
                }
            }
            return contains;
        }

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not. Implementation depends on particle shape.
        /// </summary>
        /// <param name="point"></param>
        protected virtual bool ParticleContains(Vector point, Vector Position, double tolerance = 0) => throw new NotImplementedException();

        /// <summary>
        /// Return the length-scales of the particle (length and thickness)
        /// </summary>
        public virtual double[] GetLengthScales() => throw new NotImplementedException();

        /// <summary>
        /// Returns surface points (for distance calc)
        /// </summary>
        /// <param name="hMin"></param>
        /// <param name="searchAngle"></param>
        /// <param name="subParticleID">between 0 and <see cref="NoOfSubParticles"/>, i guess</param>
        public virtual MultidimensionalArray GetSurfacePoints(double hMin, double searchAngle, int subParticleID) => throw new NotImplementedException();

        /// <summary>
        /// Calculates the support point with an analytic formula (if applicable), else it uses a binary search.
        /// </summary>
        /// <param name="supportVector"></param>
        /// <param name="SubParticleID">between 0 and <see cref="NoOfSubParticles"/>, i guess</param>
        /// <param name="Position"></param>
        public virtual Vector GetSupportPoint(Vector supportVector, Vector Position, Vector Angle, int SubParticleID, double tolerance = 0) {
            int spatialDim = Position.Dim;
            Vector currentSupportPoint = new(spatialDim);
            if (spatialDim != 2)
                throw new NotImplementedException("Calculation of support point only implemented in 2D");
            double angle = Angle[0]; // hardcoded 2D
            Vector particleDirection = new(Math.Cos(angle), Math.Sin(angle));
            double crossProductDirectionSupportVector = particleDirection[0] * supportVector[1] - particleDirection[1] * supportVector[0];
            double searchStartAngle = (1 - Math.Sign(crossProductDirectionSupportVector)) * Math.PI / 2 + Math.Acos((supportVector * particleDirection) / supportVector.L2Norm());
            double L = searchStartAngle - Math.PI;
            double R = searchStartAngle + Math.PI;
            while (L < R && Math.Abs(L - R) > 1e-15) {
                searchStartAngle = (L + R) / 2;
                double dAngle = 1e-8;
                MultidimensionalArray SurfacePoints = GetSurfacePoints(dAngle, searchStartAngle, SubParticleID);
                Vector RightNeighbour = new(spatialDim);
                Vector LeftNeighbour = new(spatialDim);
                for (int d = 0; d < spatialDim; d++) {
                    currentSupportPoint[d] = SurfacePoints[1, d];
                    LeftNeighbour[d] = SurfacePoints[0, d];
                    RightNeighbour[d] = SurfacePoints[2, d];
                }
                if ((currentSupportPoint * supportVector) > (RightNeighbour * supportVector) && (currentSupportPoint * supportVector) > (LeftNeighbour * supportVector))
                    break; // The current temp_supportPoint is the actual support point.
                else if ((RightNeighbour * supportVector) > (LeftNeighbour * supportVector))
                    L = searchStartAngle; // Search on the right side of the current point.
                else
                    R = searchStartAngle; // Search on the left side.
            }
            currentSupportPoint.Acc(Position);
            Aux.TestArithmeticException(currentSupportPoint, "supportPoint");
            return currentSupportPoint;
        }

        /// <summary>
        /// Calculates the radial vector (SurfacePoint-ParticleReadOnlyPosition)
        /// </summary>
        /// <param name="Point">
        /// </param>
        internal Vector CalculateRadialVector(Vector Point) {
            Aux = new Auxillary();
            Vector RadialVector = new(Point[0] - Motion.GetPosition(0)[0], Point[1] - Motion.GetPosition(0)[1]);
            if(RadialVector.L2Norm() > GetLengthScales().Max()) {//Point is in a different virtual domain (Periodic bndy only):
                for (int i = 0; i < Motion.OriginInVirtualPeriodicDomain.Count; i++) {
                    Vector virtualPosition = Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition();
                    Vector tempRadialVector = new(Point[0] - virtualPosition[0], Point[1] - virtualPosition[1]);
                    if (tempRadialVector.L2Norm() < RadialVector.L2Norm())
                        RadialVector = new Vector(tempRadialVector);
                }
            }
            if (RadialVector.L2Norm() == 0)
                throw new ArithmeticException("The radial vector has no length. Surface point: " + Point + " Position: " + Motion.GetPosition(0));
            Aux.TestArithmeticException(RadialVector, "particle radial vector");
            return RadialVector;
        }

        /// <summary>
        /// Calculates the eccentricity of a collision
        /// </summary>
        internal double CalculateEccentricity(Vector normalVector, Vector ClosestPoint) {
            Vector radialVector = CalculateRadialVector(ClosestPoint);
            Vector normalRadialVector = new(-radialVector[1], radialVector[0]);
            return normalRadialVector * normalVector;
        }

        internal double CalculateSecondOrderEccentricity(Vector NormalVector, Vector ClosestPoint) {
            Vector radialVector = CalculateRadialVector(ClosestPoint);
            double firstCrossProduct2D = radialVector[0] * NormalVector[1] - radialVector[1] * NormalVector[0];
            Vector secondCrossProduct2D = new(-firstCrossProduct2D * radialVector[1], firstCrossProduct2D * radialVector[0]);
            return secondCrossProduct2D * NormalVector;
        }

        /// <summary>
        /// clone, not implemented
        /// </summary>
        public virtual object Clone() => throw new NotImplementedException("Currently cloning of this type of particle is not available");
    }
}

