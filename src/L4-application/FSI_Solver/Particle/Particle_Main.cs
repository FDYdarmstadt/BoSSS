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
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using System.Collections;
using FSI_Solver;

namespace BoSSS.Application.FSI_Solver {

    /// <summary>
    /// Particle properties (for disk shape and spherical particles only).
    /// </summary>
    [DataContract]
    [Serializable]
    abstract public class Particle : ICloneable {

        /// <summary>
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        protected Particle() {
            // noop
        }

        protected static int spatialDim = 2;

        public Particle(ParticleMotionInit motionInit, double[] startPos = null, double startAngl = 0.0, double[] startTransVelocity = null, double startRotVelocity = 0) {
            spatialDim = startPos.Length;
            if (startPos == null) {
                startPos = new double[spatialDim];
            }

            m_MotionInit = motionInit;
            m_MotionInit.CheckInput();
            Motion = m_MotionInit.GetParticleMotion();
            Motion.InitializeParticlePositionAndAngle(startPos, startAngl);
            Motion.InitializeParticleVelocity(startTransVelocity, startRotVelocity);
        }

        /// <summary>
        /// Initialize particle motion.
        /// </summary>
        [DataMember]
        public ParticleMotionInit m_MotionInit;

        /// <summary>
        /// Instantiate object for particle motion.
        /// </summary>
        [DataMember]
        public Motion_Wet Motion = new Motion_Wet(new double[] { 0, 9.81 });

        /// <summary>
        /// Check whether any particles is collided with another particle
        /// </summary>
        public bool isCollided;

        /// <summary>
        /// Number of iterations
        /// </summary>
        [DataMember]
        public int iteration_counter_P = 0;

        /// <summary>
        /// Number of iterations
        /// </summary>
        [DataMember]
        public double ForceTorqueResidual;

        virtual internal int NoOfSubParticles() { return 1; }

        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double particleDensity = 1;

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double eccentricity;

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> CollisionTranslationalVelocity = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> collisionNormalVector = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> collisionTangentialVector = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double[] closestPointToOtherObject = new double[spatialDim];

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double[] closestPointOnOtherObjectToThis = new double[spatialDim];

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double[] TotalCollisionPositionCorrection = new double[spatialDim];

        /// <summary>
        /// The angular velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double> CollisionRotationalVelocity = new List<double>();

        /// <summary>
        /// Level set function describing the particle.
        /// </summary>       
        public abstract double Phi_P(double[] X);

        /// <summary>
        /// Convergence criterion for the calculation of the Forces and Torque
        /// </summary>
        [DataMember]
        public double forceAndTorque_convergence = 1e-8;

        /// <summary>
        /// Active stress on the current particle.
        /// </summary>
        public double activeStress = 0;

        /// <summary>
        /// Area of the current particle.
        /// </summary>
        virtual public double Area_P() {
            throw new NotImplementedException("");
        }

        /// <summary>
        /// Necessary for active particles. Returns 0 for the non active boundary region and a number between 0 and 1 for the active region.
        /// </summary>
        internal double SeperateBoundaryRegions(double[] X) {
            double seperateBoundaryRegions;
            // The posterior side of the particle 
            if (Math.Cos(Motion.angle[0]) * (X[0] - Motion.position[0][0]) + Math.Sin(Motion.angle[0]) * (X[1] - Motion.position[0][1]) < 1e-8) {
                seperateBoundaryRegions = (Math.Cos(Motion.angle[0]) * (X[0] - Motion.position[0][0]) + Math.Sin(Motion.angle[0]) * (X[1] - Motion.position[0][1]));
                seperateBoundaryRegions /= Math.Sqrt((X[0] - Motion.position[0][0]).Pow2() + (X[1] - Motion.position[0][1]).Pow2());
            }
            // The anterior side of the particle 
            else {
                seperateBoundaryRegions = 0;
            }
            return seperateBoundaryRegions;
        }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        public double Mass_P {
            get {
                Aux.TestArithmeticException(Area_P(), "particle area");
                Aux.TestArithmeticException(particleDensity, "particle density");
                return Area_P() * particleDensity;
            }
        }

        /// <summary>
        /// Circumference of the current particle.
        /// </summary>
        abstract protected double Circumference_P { get; }

        /// <summary>
        /// Moment of inertia of the current particle.
        /// </summary>
        abstract public double MomentOfInertia_P { get; }

        [NonSerialized]
        readonly internal FSI_Auxillary Aux = new FSI_Auxillary();

        /// <summary>
        /// clone, not implemented
        /// </summary>
        virtual public object Clone() {
            throw new NotImplementedException("Currently cloning of a particle is not available");
        }

        public double[] CalculateParticleMomentum() {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++) {
                temp[d] = Mass_P * Motion.translationalVelocity[0][d];
            }
            temp[spatialDim] = MomentOfInertia_P * Motion.rotationalVelocity[0];
            return temp;
        }

        public double[] CalculateParticleKineticEnergy() {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++) {
                temp[d] = 0.5 * Mass_P * Motion.translationalVelocity[0][d].Pow2();
            }
            temp[spatialDim] = 0.5 * MomentOfInertia_P * Motion.rotationalVelocity[0].Pow2();
            return temp;
        }

        /// <summary>
        /// Calculating the particle reynolds number
        /// </summary>
        public double ComputeParticleRe(double ViscosityFluid) {
            return Motion.translationalVelocity[0].L2Norm() * GetLengthScales().Max() / ViscosityFluid;
        }

        public double ComputeParticleSt(double ViscosityFluid, double DensityFluid) {
            return ComputeParticleRe(ViscosityFluid) * particleDensity / (9 * DensityFluid);
        }

        public double ComputeParticleRe(double ViscosityFluid, double[] relativeVelocity) {
            return relativeVelocity.L2Norm() * GetLengthScales().Max() / ViscosityFluid;
        }

        public double ComputeParticleSt(double ViscosityFluid, double DensityFluid, double[] relativeVelocity) {
            return ComputeParticleRe(ViscosityFluid, relativeVelocity) * particleDensity / (9 * DensityFluid);
        }

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
        /// <returns></returns>
        public abstract bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false);

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public abstract bool ParticleInternalCell(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false);

        virtual public double[] GetLengthScales() {
            throw new NotImplementedException();
        }

        virtual public MultidimensionalArray GetSurfacePoints(double hMin) {
            throw new NotImplementedException();
        }

        virtual public void GetSupportPoint(int SpatialDim, double[] Vector, out double[] SupportPoint) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates the radial vector (SurfacePoint-ParticlePosition)
        /// </summary>
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialVector">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        internal void CalculateRadialVector(double[] SurfacePoint, out double[] RadialVector, out double RadialLength) {
            RadialVector = new double[] { SurfacePoint[0] - Motion.position[0][0], SurfacePoint[1] - Motion.position[0][1] };
            RadialLength = RadialVector.L2Norm();
            RadialVector.ScaleV(1 / RadialLength);
            Aux.TestArithmeticException(RadialVector, "particle radial vector");
            Aux.TestArithmeticException(RadialLength, "particle radial length");
        }

        internal void CalculateRadialNormalVector(double[] SurfacePoint, out double[] RadialNormalVector) {
            RadialNormalVector = new double[] { SurfacePoint[1] - Motion.position[0][1], -SurfacePoint[0] + Motion.position[0][0] };
            RadialNormalVector.ScaleV(1 / RadialNormalVector.L2Norm());
            Aux.TestArithmeticException(RadialNormalVector, "particle vector normal to radial vector");
        }

        internal void CalculateEccentricity() {
            CalculateRadialVector(closestPointToOtherObject, out double[] RadialVector, out _);
            double[] tangentialVector = collisionTangentialVector.Last();
            eccentricity = RadialVector[0] * tangentialVector[0] + RadialVector[1] * tangentialVector[1];
            Aux.TestArithmeticException(eccentricity, "particle eccentricity");
        }
    }
}

