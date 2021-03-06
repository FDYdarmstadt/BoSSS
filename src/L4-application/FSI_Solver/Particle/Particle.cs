﻿/* =======================================================================
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
using System.Collections.Generic;
using System.Linq;

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
        public Particle(InitializeMotion motionInit, double[] startPos, double startAngl = 0.0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) {
            SpatialDim = startPos.Length;
            ActiveStress = activeStress;
            Aux = new FSIAuxillary();

            if (motionInit != null) {
                motionInit.CheckInput();
                Motion = motionInit.ParticleMotion;
                Motion.InitializeParticlePositionAndAngle(startPos, startAngl);
                Motion.InitializeParticleVelocity(startTransVelocity, startRotVelocity);
                particleDensity = Motion.Density;
                MotionInitializer = motionInit;
            }
        }

        

        [DataMember]
        public bool IsMaster = true;

        [DataMember]
        public int[] MasterDuplicateIDs = new int[4];

        public void SetMaster(Motion newMotionType) {
            newMotionType.TransferDataFromOtherParticle(Motion);
            Motion = newMotionType;
            IsMaster = true;
        }

        public void SetDuplicate() {
            IsMaster = false;
        }

        public void SetDuplicateHierachy(int[] hierachy) {
            MasterDuplicateIDs = hierachy.CloneAs();
        }
        
        [NonSerialized]
        protected FSIAuxillary Aux;
        [DataMember]
        private readonly double particleDensity;

        [DataMember]
        protected int SpatialDim;

        /// <summary>
        /// Instantiate object for particle motion.
        /// </summary>
        [DataMember]
        public Motion Motion { get; private set; }

        [DataMember]
        public InitializeMotion MotionInitializer;

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
        /// Circumference of the current particle.
        /// </summary>
        public abstract double Circumference { get; }

        /// <summary>
        /// Moment of inertia of the current particle.
        /// </summary>
        public abstract double MomentOfInertia { get; }

        /// <summary>
        /// Level set function describing the particle.
        /// </summary>   
        public double LevelSetFunction(double[] X, double GridLength) {
            double levelSet = ParticleLevelSetFunction(X, Motion.GetPosition());
            for (int i = 0; i < Motion.OriginInVirtualPeriodicDomain.Count(); i++) {
                Vector virtualPosition = Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition();
                if (Motion.IsInsideOfPeriodicDomain(virtualPosition, (3 * GridLength + GetLengthScales().Max())))
                    levelSet = Math.Max(levelSet, ParticleLevelSetFunction(X, Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition()));
            }
            return levelSet;
        }

        protected virtual double ParticleLevelSetFunction(double[] X, Vector Postion) => throw new NotImplementedException();

        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        public CellMask ParticleCutCells(LevelSetTracker LsTrk, CellMask AllCutCells) {
            BitArray CellArray = new BitArray(LsTrk.GridDat.Cells.NoOfLocalUpdatedCells);
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            var h_min = LsTrk.Regions.GetCutCellSubGrid().h_minSubGrd;
            for (int i = 0; i < CellArray.Length; i++) {
                CellArray[i] = Contains(new Vector(CellCenters[i, 0], CellCenters[i, 1]), 2 * h_min);
            }
            CellMask CutCells = new CellMask(LsTrk.GridDat, CellArray, MaskType.Logical);
            CutCells = CutCells.Intersect(AllCutCells);
            return CutCells;
        }

        public bool Contains(Vector Point, double Tolerance = 0) {
            bool contains = ParticleContains(Point, Motion.GetPosition(), Tolerance);
            if (!contains) {
                for (int i = 0; i < Motion.OriginInVirtualPeriodicDomain.Count(); i++) {
                    Vector virtualPosition = Motion.OriginInVirtualPeriodicDomain[i] + Motion.GetPosition();
                    if (Motion.IsInsideOfPeriodicDomain(virtualPosition, Tolerance + GetLengthScales().Max()))
                        contains = ParticleContains(Point, virtualPosition, Tolerance);
                }
            }
            return contains;
        }

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
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
        public virtual MultidimensionalArray GetSurfacePoints(double hMin, double searchAngle, int subParticleID) => throw new NotImplementedException();

        /// <summary>
        /// Calculates the support point with an analytic formula (if applicable)
        /// </summary>
        /// <param name="supportVector"></param>
        public virtual Vector GetSupportPoint(Vector supportVector, Vector Position, int SubParticleID) {
            int spatialDim = Position.Dim;
            Vector currentSupportPoint = new Vector(spatialDim);
            double angle = Motion.GetAngle(0);
            Vector particleDirection = new Vector(Math.Cos(angle), Math.Sin(angle));
            double crossProductDirectionSupportVector = particleDirection[0] * supportVector[1] - particleDirection[1] * supportVector[0];
            double searchStartAngle = (1 - Math.Sign(crossProductDirectionSupportVector)) * Math.PI / 2 + Math.Acos((supportVector * particleDirection) / supportVector.L2Norm());
            double L = searchStartAngle - Math.PI;
            double R = searchStartAngle + Math.PI;
            while (L < R && Math.Abs(L - R) > 1e-15) {
                searchStartAngle = (L + R) / 2;
                double dAngle = 1e-8;
                MultidimensionalArray SurfacePoints = GetSurfacePoints(dAngle, searchStartAngle, SubParticleID);
                Vector RightNeighbour = new Vector(spatialDim);
                Vector LeftNeighbour = new Vector(spatialDim);
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
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialVector">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        internal Vector CalculateRadialVector(Vector SurfacePoint) {
            Aux = new FSIAuxillary();
            Vector RadialVector = new Vector(SurfacePoint[0] - Motion.GetPosition(0)[0], SurfacePoint[1] - Motion.GetPosition(0)[1]);
            if (RadialVector.L2Norm() == 0)
                throw new ArithmeticException("The radial vector has no length. Surface point: " + SurfacePoint + " Position: " + Motion.GetPosition(0));
            Aux.TestArithmeticException(RadialVector, "particle radial vector");
            return RadialVector;
        }

        /// <summary>
        /// Calculates the eccentricity of a collision
        /// </summary>
        internal double CalculateEccentricity(Vector normalVector, Vector ClosestPoint) {
            Vector radialVector = CalculateRadialVector(ClosestPoint);
            Vector normalRadialVector = new Vector(-radialVector[1], radialVector[0]);
            return normalRadialVector * normalVector;
        }

        internal double CalculateSecondOrderEccentricity(Vector NormalVector, Vector ClosestPoint) {
            Vector radialVector = CalculateRadialVector(ClosestPoint);
            double firstCrossProduct2D = radialVector[0] * NormalVector[1] - radialVector[1] * NormalVector[0];
            Vector secondCrossProduct2D = new Vector(-firstCrossProduct2D * radialVector[1], firstCrossProduct2D * radialVector[0]);
            return secondCrossProduct2D * NormalVector;
        }

        /// <summary>
        /// clone, not implemented
        /// </summary>
        public virtual object Clone() => throw new NotImplementedException("Currently cloning of this type of particle is not available");
    }
}

