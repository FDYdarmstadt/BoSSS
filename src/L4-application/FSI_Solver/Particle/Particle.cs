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
using System.Runtime.Serialization;
using BoSSS.Solution.Control;
using System.Runtime.InteropServices;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using MPI.Wrappers;
using NUnit.Framework;
using FSI_Solver;
using System.Collections;

namespace BoSSS.Application.FSI_Solver
{

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
        protected Particle()
        {
            // noop
        }

        private const int historyLength = 4;
        protected static int spatialDim = 2;

        public Particle(int Dim, double[] startPos = null, double startAngl = 0.0) {
            
            spatialDim = Dim;

            // Particle history
            // =============================   
            for (int i = 0; i < historyLength; i++)
            {
                position.Add(new double[Dim]);
                angle.Add(new double());
                translationalVelocity.Add(new double[Dim]);
                translationalAcceleration.Add(new double[Dim]);
                rotationalVelocity.Add(new double());
                rotationalAcceleration.Add(new double());
                HydrodynamicForces.Add(new double[Dim]);
                HydrodynamicTorque.Add(new double());
            }

            // ============================= 
            if (startPos == null) {
                startPos = new double[Dim];
            }
            position[0] = startPos;
            position[1] = startPos;
            //From degree to radiant
            angle[0] = StartingAngle = startAngl * 2 * Math.PI / 360;
            angle[1] = startAngl * 2 * Math.PI / 360;
        }
        
        /// <summary>
        /// Set true if translation of the particle should be induced by hydrodynamical forces.
        /// </summary>
        [DataMember]
        public bool IncludeTranslation = true;

        /// <summary>
        /// Set true if rotation of the particle should be induced by hydrodynamical torque.
        /// </summary>
        [DataMember]
        public bool IncludeRotation = true;

        /// <summary>
        /// Check whether any particles is collided with another particle
        /// </summary>
        public bool Collided;
        
        /// <summary>
        /// Skip calculation of hydrodynamic force and Torque if particles are too close -----> to be tested whether it is still necessary
        /// </summary>
        [DataMember]
        public bool skipForceIntegration = false;

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

        /// <summary>
        /// Constant Forces and Torque underrelaxation?
        /// </summary>
        [DataMember]
        public bool useAddaptiveUnderrelaxation = false;

        /// <summary>
        /// Defines the order of the underrelaxation factor
        /// </summary>
        [DataMember]
        public int underrelaxationFT_exponent = 0;

        /// <summary>
        /// Underrelaxation factor
        /// </summary>
        [DataMember]
        public double underrelaxation_factor = 1;

        /// <summary>
        /// Set true if you want to delete all values of the Forces anf Torque smaller than convergenceCriterion*1e-2
        /// </summary>
        [DataMember]
        public bool clearSmallValues = false;

        /// <summary>
        /// The color of the particle.
        /// </summary>
        public int ParticleColor = new int();

        /// <summary>
        /// Colored cells of this particle. 0: CellID, 1: Color
        /// </summary>
        public int[] ParticleColoredCells;

        /// <summary>
        /// Set false if you want to include the effects of added damping
        /// </summary>
        [DataMember]
        public bool UseAddedDamping = false;

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        [DataMember]
        public double[,] addedDampingTensor = new double[6, 6];

        /// <summary>
        /// AddedDampingCoefficient
        /// </summary>
        [DataMember]
        public double addedDampingCoefficient = 1;
        
        virtual internal int NoOfSubParticles() { return 1; }

        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double particleDensity = 1;
        
        /// <summary>
        /// The position (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> position = new List<double[]>();
        
        /// <summary>
        /// The angle (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> angle = new List<double>();

        ///// <summary>
        /// The angle (center of mass) of the particle at the starting point.
        /// </summary>
        [DataMember]
        private readonly double StartingAngle = new double();
        
        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> translationalVelocity = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double[] PreCollisionVelocity;

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
        public List<double[]> CollisionNormalVector = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> CollisionTangentialVector = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double CollisionTimestep = new double();

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
        public double CollisionPreviousTimestep = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public double[] TotalCollisionPositionCorrection = new double[spatialDim];

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> rotationalVelocity = new List<double>();

        /// <summary>
        /// The angular velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double> CollisionRotationalVelocity = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> translationalAcceleration = new List<double[]>();
        
        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> rotationalAcceleration = new List<double>();
        
        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> HydrodynamicForces = new List<double[]>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public double[] forcesPrevIteration = new double[spatialDim];

        /// <summary>
        /// The Torque acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> HydrodynamicTorque = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public double torquePrevIteration = new double();


        /// <summary>
        /// Level set function describing the particle.
        /// </summary>       
        public abstract double Phi_P(double[] X);

        /// <summary>
        /// Sets the gravity in vertical direction, default is 0.0
        /// </summary>
        [DataMember]
        public double GravityVertical = 0.0;

        /// <summary>
        /// Set true if the particle should be an active particle, i.e. self driven
        /// </summary>
        [DataMember]
        public bool ActiveParticle = false;

        /// <summary>
        /// Convergence criterion for the calculation of the Forces and Torque
        /// </summary>
        [DataMember]
        public double forceAndTorque_convergence = 1e-8;

        /// <summary>
        /// Active stress on the current particle.
        /// </summary>
        public double ActiveStress = 0;

        /// <summary>
        /// Active velocity (alternative to active stress) on the current particle.
        /// </summary>
        [DataMember]
        public double ActiveVelocity;

        /// <summary>
        /// Area of the current particle.
        /// </summary>
        abstract public double Area_P { get; }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        public double Mass_P
        {
            get
            {
                if (Area_P <= 0.0 || double.IsNaN(Area_P) || double.IsInfinity(Area_P))
                    throw new ArithmeticException("Particle volume/area is " + Area_P);
                if (particleDensity <= 0.0 || double.IsNaN(particleDensity) || double.IsInfinity(particleDensity))
                    throw new ArithmeticException("Particle volume/area is " + particleDensity);

                return Area_P * particleDensity;
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
        readonly internal ParticleAuxillary Aux = new ParticleAuxillary();
        [NonSerialized]
        readonly private ParticleForceIntegration ForceIntegration = new ParticleForceIntegration();
        [NonSerialized]
        readonly private ParticleAddedDamping AddedDamping = new ParticleAddedDamping();
        [NonSerialized]
        readonly private ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
        [NonSerialized]
        readonly private ParticleAcceleration Acceleration = new ParticleAcceleration();

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public void CalculateParticlePosition(double dt)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(position);
            }

            if (spatialDim != 2 && spatialDim != 3)
                throw new NotSupportedException("Unknown particle dimension: SpatialDim = " + spatialDim);

            int ClearAcceleartion = CollisionTimestep != 0 ? 0 : 1;
            if (IncludeTranslation == true) {
                for (int d = 0; d < spatialDim; d++)
                {
                    position[0][d] = position[1][d] + ( translationalVelocity[0][d] + ClearAcceleartion * (4 * translationalVelocity[1][d] + translationalVelocity[2][d])) * (dt - CollisionTimestep) / 6;
                    if (double.IsNaN(position[0][d]) || double.IsInfinity(position[0][d]))
                        throw new ArithmeticException("Error trying to update particle position. Value:  " + position[0][d]);
                }
            }
            else
            {
                for (int d = 0; d < spatialDim; d++)
                {
                    position[0][d] = position[1][d];
                    translationalAcceleration[0][d] = 0;
                    translationalVelocity[0][d] = 0;
                }
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public void CalculateParticleAngle(double dt)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveValueOfLastTimestep(angle);
            }

            if (spatialDim != 2)
                throw new NotSupportedException("Unknown particle dimension: SpatialDim = " + spatialDim);

            int ClearAcceleartion = CollisionTimestep != 0 ? 0 : 1;
            if (IncludeRotation == true)
            {
                angle[0] = angle[1] + (rotationalVelocity[0] + ClearAcceleartion * (4 * rotationalVelocity[1] + rotationalVelocity[2])) * (dt - CollisionTimestep) / 6;
                if (double.IsNaN(angle[0]) || double.IsInfinity(angle[0]))
                    throw new ArithmeticException("Error trying to update particle angle. Value:  " + angle[0]);
            }
            else
            {
                angle[0] = angle[1];
                rotationalAcceleration[0] = 0;
                rotationalVelocity[0] = 0;
            }
        }

        /// <summary>
        /// Predict the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public void PredictAcceleration()
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(translationalAcceleration);
                Aux.SaveValueOfLastTimestep(rotationalAcceleration);
                //Aux.SaveMultidimValueOfLastTimestep(HydrodynamicForces);
                //Aux.SaveValueOfLastTimestep(HydrodynamicTorque);
            }

            for (int d = 0; d < spatialDim; d++)
            {
                translationalAcceleration[0][d] = (translationalAcceleration[1][d] + 4 * translationalAcceleration[2][d] + translationalAcceleration[3][d]) / 8;
                if (Math.Abs(translationalAcceleration[0][d]) < 1e-20)
                    translationalAcceleration[0][d] = 0;
            }

            rotationalAcceleration[0] = (rotationalAcceleration[1] + 4 * rotationalAcceleration[2] + rotationalAcceleration[3]) / 8;
            if (Math.Abs(rotationalAcceleration[0]) < 1e-20)
                rotationalAcceleration[0] = 0;
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public void CalculateAcceleration(double dt, bool FullyCoupled, bool IncludeHydrodynamics)
        {
            if (iteration_counter_P == 0 || FullyCoupled == false)
            {
                Aux.SaveMultidimValueOfLastTimestep(translationalAcceleration);
                Aux.SaveValueOfLastTimestep(rotationalAcceleration);
            }

            // Include gravitiy for dry simulations
            if (!Collided && !IncludeHydrodynamics)
            {
                HydrodynamicForces[0][1] += GravityVertical * Mass_P;
            }

            double[,] CoefficientMatrix = Acceleration.CalculateCoefficientMatrix(addedDampingTensor, Mass_P, MomentOfInertia_P, dt, addedDampingCoefficient);
            double Denominator = Acceleration.CalculateDenominator(CoefficientMatrix);

            if (IncludeTranslation) 
                translationalAcceleration[0] = Acceleration.Translational(CoefficientMatrix, Denominator, HydrodynamicForces[0], HydrodynamicTorque[0]);

            for (int d = 0; d < spatialDim; d++)
            {
                if (Math.Abs(translationalAcceleration[0][d]) < 1e-20 || IncludeTranslation == false)
                    translationalAcceleration[0][d] = 0;
            }

            if (IncludeRotation)
                rotationalAcceleration[0] = Acceleration.Rotational(CoefficientMatrix, Denominator, HydrodynamicForces[0], HydrodynamicTorque[0]);
            if (Math.Abs(rotationalAcceleration[0]) < 1e-20 || IncludeRotation == false)
                rotationalAcceleration[0] = 0;
        }

        internal void UpdateParticleVelocity(double dt)
        {
            CalculateTranslationalVelocity(dt);
            CalculateAngularVelocity(dt);
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public void CalculateTranslationalVelocity(double dt)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(translationalVelocity);
            }

            double[] tempActiveVelcotiy = new double[2];

            int ClearAcceleartion = CollisionTimestep != 0 ? 0 : 1;

            if (this.IncludeTranslation == false)
            {
                for (int d = 0; d < spatialDim; d++)
                {
                    translationalVelocity[0][d] = 0;
                }
            }
            else if (ActiveVelocity != 0)
            {
                tempActiveVelcotiy[0] = Math.Cos(angle[0]) * ActiveVelocity;
                tempActiveVelcotiy[1] = Math.Sin(angle[0]) * ActiveVelocity;
                for (int d = 0; d < spatialDim; d++)
                {
                    translationalVelocity[0][d] = tempActiveVelcotiy[d];
                    if (double.IsNaN(translationalVelocity[0][d]) || double.IsInfinity(translationalVelocity[0][d]))
                        throw new ArithmeticException("Error trying to calculate particle velocity Value:  " + translationalVelocity[0][d]);
                }
            }
            else
            {
                for (int d = 0; d < spatialDim; d++) {

                    //TranslationalVelocity[0][d] = TranslationalVelocity[1][d] + (TranslationalAcceleration[1][d] + TranslationalAcceleration[0][d]) * dt / 2;
                    translationalVelocity[0][d] = translationalVelocity[1][d] + (translationalAcceleration[0][d] + ClearAcceleartion * (4 * translationalAcceleration[1][d] + translationalAcceleration[2][d])) * dt / 6;
                    if (double.IsNaN(translationalVelocity[0][d]) || double.IsInfinity(translationalVelocity[0][d]))
                        throw new ArithmeticException("Error trying to calculate particle velocity Value:  " + translationalVelocity[0][d]);
                }
            }
        }
        
        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public void CalculateAngularVelocity(double dt)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveValueOfLastTimestep(rotationalVelocity);
            }

            int ClearAcceleartion = CollisionTimestep != 0 ? 0 : 1;
            if (this.IncludeRotation == false) {
                rotationalVelocity[0] = 0;
                return;
            }
            else
            {
                rotationalVelocity[0] = rotationalVelocity[1] + dt * (rotationalAcceleration[0] + ClearAcceleartion * (4 * rotationalAcceleration[1] + rotationalAcceleration[2])) / 6;
                if (double.IsNaN(rotationalVelocity[0]) || double.IsInfinity(rotationalVelocity[0]))
                    throw new ArithmeticException("Error trying to calculate particle angluar velocity. Value:  " + rotationalVelocity[0]);
            }
        }
        
        /// <summary>
        /// clone, not implemented
        /// </summary>
        virtual public object Clone() {
            throw new NotImplementedException("Currently cloning of a particle is not available");
        }

        /// <summary>
        /// Calculate tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public void CalculateDampingTensor(Particle particle, LevelSetTracker LsTrk, double muA, double rhoA, double dt)
        {
            addedDampingTensor = AddedDamping.IntegrationOverLevelSet(particle, LsTrk, muA, rhoA, dt, position[0]);
        }

        /// <summary>
        /// Update in every timestep tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public void UpdateDampingTensors()
        {
            addedDampingTensor = AddedDamping.RotateTensor(angle[0], StartingAngle, addedDampingTensor);
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public void PredictForceAndTorque(int TimestepInt)
        {
            if (TimestepInt == 1)
            {
                HydrodynamicForces[0][0] = 20 * Math.Cos(angle[0]) * ActiveStress * Circumference_P;
                HydrodynamicForces[0][1] = 20 * Math.Sin(angle[0]) * ActiveStress * Circumference_P + GravityVertical * Mass_P;
            }
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(translationalAcceleration);
                Aux.SaveValueOfLastTimestep(rotationalAcceleration);
                Aux.SaveMultidimValueOfLastTimestep(HydrodynamicForces);
                Aux.SaveValueOfLastTimestep(HydrodynamicTorque);
            }
            for (int d = 0; d < spatialDim; d++)
            {
                HydrodynamicForces[0][d] = (HydrodynamicForces[1][d] + 4 * HydrodynamicForces[2][d] + HydrodynamicForces[3][d]) / 6;
                if (Math.Abs(HydrodynamicForces[0][d]) < 1e-20)
                    HydrodynamicForces[0][d] = 0;
            }
            HydrodynamicTorque[0] = (HydrodynamicTorque[1] + 4 * HydrodynamicTorque[2] + HydrodynamicTorque[3]) / 6;
            if (Math.Abs(HydrodynamicTorque[0]) < 1e-20)
                HydrodynamicTorque[0] = 0;
        }

        public void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, double muA, double dt, double fluidDensity, bool firstIteration)
        {
            CalculateHydrodynamicForces(U, P, LsTrk, muA, dt, fluidDensity, out double[] tempForces);
            CalculateHydrodynamicTorque(U, P, LsTrk, muA, dt, out double tempTorque);
            HydrodynamicsPostprocessing(tempForces, tempTorque, firstIteration);
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        private void CalculateHydrodynamicForces(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, double muA, double dt, double fluidDensity, out double[] tempForces)
        { 
            HydrodynamicForces[0][0] = 0;
            HydrodynamicForces[0][1] = 0;
            tempForces = new double[spatialDim];

            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P; 

            if (IncludeTranslation)
                ForcesIntegration(UA, pA, LsTrk, RequiredOrder, muA, out tempForces);
            Force_MPISum(ref tempForces);
            CalculateGravitationalForce(ref tempForces, fluidDensity);
            if (UseAddedDamping)
                ForceAddedDamping(ref tempForces, dt);
        }

        private void ForcesIntegration(SinglePhaseField[] UA, ConventionalDGField pA, LevelSetTracker LsTrk, int RequiredOrder, double FluidViscosity, out double[] tempForces)
        {
            tempForces = new double[spatialDim];
            double[] IntegrationForces = tempForces.CloneAs();
            for (int d = 0; d < spatialDim; d++)
            {
                void ErrFunc(int CurrentCellID, int Length, NodeSet Ns, MultidimensionalArray result)
                {

                    int K = result.GetLength(1);
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Length, K, spatialDim, spatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Length, K);
                    MultidimensionalArray Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, CurrentCellID, Length);
                    for (int i = 0; i < spatialDim; i++)
                    {
                        UA[i].EvaluateGradient(CurrentCellID, Length, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    pA.Evaluate(CurrentCellID, Length, Ns, pARes);
                    for (int j = 0; j < Length; j++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            result[j, k] = ForceIntegration.CalculateStressTensor(Grad_UARes, pARes, Normals, FluidViscosity, k, j, spatialDim, d);
                        }
                    }
                }
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, CutCells_P(LsTrk));
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat, cqs.Compile(LsTrk.GridDat, RequiredOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult)
                    {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                    {
                        IntegrationForces[d] = ForceTorqueSummationWithNeumaierArray(IntegrationForces[d], ResultsOfIntegration, Length);
                    }
                ).Execute();
            }
            tempForces = IntegrationForces.CloneAs();
        }

        private void Force_MPISum(ref double[] forces)
        {
            double[] stateBuffer = forces.CloneAs();
            double[] globalStateBuffer = stateBuffer.MPISum();
            for (int d = 0; d < spatialDim; d++)
            {
                forces[d] = globalStateBuffer[d];
            }
        }

        private void CalculateGravitationalForce(ref double[] Forces, double fluidDensity)
        {
            Forces[1] += (particleDensity - fluidDensity) * Area_P * GravityVertical;
        }

        private void ForceAddedDamping(ref double[] forces, double dt)
        {
            for (int d = 0; d < spatialDim; d++)
            {
                forces[d] += addedDampingCoefficient * dt * (addedDampingTensor[0, d] * translationalAcceleration[0][0] + addedDampingTensor[1, d] * translationalAcceleration[0][1] + addedDampingTensor[d, 2] * rotationalAcceleration[0]);
            }
        }

        private void CalculateHydrodynamicTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, double muA, double dt, out double tempTorque)
        {
            HydrodynamicTorque[0] = 0;
            tempTorque = new double();
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;

            if (IncludeRotation)
                TorqueIntegration(UA, pA, LsTrk, RequiredOrder, muA, out tempTorque);

            Torque_MPISum(ref tempTorque);

            if (UseAddedDamping)
                TorqueAddedDamping(ref tempTorque, dt);
        }

        private void TorqueIntegration(SinglePhaseField[] UA, ConventionalDGField pA, LevelSetTracker LsTrk, int RequiredOrder, double FluidViscosity, out double tempTorque)
        {
            double IntegrationTorque = new double();
            void ErrFunc2(int j0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                int K = result.GetLength(1);
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, spatialDim, spatialDim); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                for (int i = 0; i < spatialDim; i++)
                {
                    UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }
                pA.Evaluate(j0, Len, Ns, pARes);
                for (int j = 0; j < Len; j++)
                {
                    MultidimensionalArray Ns_Global = Ns.CloneAs();
                    LsTrk.GridDat.TransformLocal2Global(Ns, Ns_Global, j0 + j);
                    for (int k = 0; k < K; k++)
                    {
                        result[j, k] = ForceIntegration.CalculateTorqueFromStressTensor2D(Grad_UARes, pARes, Normals, Ns_Global, FluidViscosity, k, j, position[0]);
                    }
                }
            }
            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, CutCells_P(LsTrk));
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat, cqs2.Compile(LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult)
                {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                {
                    IntegrationTorque = ForceTorqueSummationWithNeumaierArray(IntegrationTorque, ResultsOfIntegration, Length);
                }
            ).Execute();
            tempTorque = IntegrationTorque;
        }

        private void Torque_MPISum(ref double torque)
        {
            double stateBuffer = torque;
            double globalStateBuffer = stateBuffer.MPISum();
            torque = globalStateBuffer;
        }

        private void TorqueAddedDamping(ref double torque, double dt)
        {
            torque += addedDampingCoefficient * dt * (addedDampingTensor[2, 0] * translationalAcceleration[0][0] + addedDampingTensor[2, 1] * translationalAcceleration[0][1] + addedDampingTensor[2, 2] * rotationalAcceleration[0]);
        }

        private void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque, bool firstIteration)
        {
            if (!firstIteration)
            {
                Underrelaxation.CalculateAverageForces(tempForces, tempTorque, GetLengthScales().Max(), out double averagedForces);
                Underrelaxation.Forces(ref tempForces, forcesPrevIteration, forceAndTorque_convergence, underrelaxation_factor, useAddaptiveUnderrelaxation, averagedForces);
                Underrelaxation.Torque(ref tempTorque, torquePrevIteration, forceAndTorque_convergence, underrelaxation_factor, useAddaptiveUnderrelaxation, averagedForces);
            }
            ForceClearSmallValues(tempForces);
            TorqueClearSmallValues(tempTorque);
            ForcesCheckPlausibility();
            TorqueCheckPlausibility();
        }

        private void ForceClearSmallValues(double[] tempForces)
        {
            for (int d = 0; d < spatialDim; d++)
            {
                HydrodynamicForces[0][d] = 0;
                if (Math.Abs(tempForces[d]) > forceAndTorque_convergence * 1e-2 || !clearSmallValues)
                    HydrodynamicForces[0][d] = tempForces[d];
            }
        }

        private void TorqueClearSmallValues(double tempTorque)
        {
            HydrodynamicTorque[0] = 0;
            if (Math.Abs(tempTorque) > forceAndTorque_convergence * 1e-2 || !clearSmallValues)
                HydrodynamicTorque[0] = tempTorque;
        }

        private void ForcesCheckPlausibility()
        {
            if (double.IsNaN(HydrodynamicForces[0][0]) || double.IsInfinity(HydrodynamicForces[0][0]))
                throw new ArithmeticException("Error trying to calculate hydrodynamic forces (x). Value:  " + HydrodynamicForces[0][0]);
            if (double.IsNaN(HydrodynamicForces[0][1]) || double.IsInfinity(HydrodynamicForces[0][1]))
                throw new ArithmeticException("Error trying to calculate hydrodynamic forces (y). Value:  " + HydrodynamicForces[0][1]);
        }

        private void TorqueCheckPlausibility()
        {
            if (double.IsNaN(HydrodynamicTorque[0]) || double.IsInfinity(HydrodynamicTorque[0]))
                throw new ArithmeticException("Error trying to calculate hydrodynamic torque. Value:  " + HydrodynamicTorque[0]);
        }

        public double[] CalculateParticleMomentum()
        {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++)
            {
                temp[d] = Mass_P * translationalVelocity[0][d];
            }
            temp[spatialDim] = MomentOfInertia_P * rotationalVelocity[0];
            return temp;
        }

        public double[] CalculateParticleKineticEnergy()
        {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++)
            {
                temp[d] = 0.5 * Mass_P * translationalVelocity[0][d].Pow2();
            }
            temp[spatialDim] = 0.5 * MomentOfInertia_P * rotationalVelocity[0].Pow2();
            return temp;
        }
        
        /// <summary>
        /// Calculating the particle reynolds number
        /// </summary>
        public double ComputeParticleRe(double ViscosityFluid)
        {
            return translationalVelocity[0].L2Norm() * GetLengthScales().Max() / ViscosityFluid;
        }
        
        public double ComputeParticleSt(double ViscosityFluid, double DensityFluid)
        {
            return ComputeParticleRe(ViscosityFluid) * particleDensity / (9 * DensityFluid);
        }

        public double ComputeParticleRe(double ViscosityFluid, double[] relativeVelocity)
        {
            return relativeVelocity.L2Norm() * GetLengthScales().Max() / ViscosityFluid;
        }

        public double ComputeParticleSt(double ViscosityFluid, double DensityFluid, double[] relativeVelocity)
        {
            return ComputeParticleRe(ViscosityFluid, relativeVelocity) * particleDensity / (9 * DensityFluid);
        }

        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        public CellMask CutCells_P(LevelSetTracker LsTrk)
        {
            BitArray CellArray = new BitArray(LsTrk.GridDat.Cells.NoOfLocalUpdatedCells);
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            double h_min = LsTrk.GridDat.Cells.h_minGlobal;
            double h_max = LsTrk.GridDat.Cells.h_maxGlobal;

            for (int i = 0; i < CellArray.Length; i++)
            {
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

        virtual public double[] GetLengthScales()
        {
            throw new NotImplementedException();
        }

        virtual public MultidimensionalArray GetSurfacePoints(double hMin)
        {
            throw new NotImplementedException();
        }

        virtual public void GetSupportPoint(int SpatialDim, double[] Vector, out double[] SupportPoint)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates the radial vector (SurfacePoint-ParticlePosition)
        /// </summary>
        /// <param name="ParticlePosition">
        /// </param>
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialVector">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        /// <param name="RadialNormalVector">
        /// Vector normal to the radial vector.
        /// </param>
        internal void CalculateRadialVector(double[] SurfacePoint, out double[] RadialVector, out double RadialLength)
        {
            RadialVector = new double[] { SurfacePoint[0] - position[0][0], SurfacePoint[1] - position[0][1] };
            RadialLength = RadialVector.L2Norm();
            RadialVector.ScaleV(1 / RadialLength);
        }

        internal void CalculateRadialNormalVector(double[] SurfacePoint, out double[] RadialNormalVector)
        {
            RadialNormalVector = new double[] { SurfacePoint[1] - position[0][1], -SurfacePoint[0] + position[0][0] };
            RadialNormalVector.ScaleV(1 / RadialNormalVector.L2Norm());
        }

        internal void CalculateEccentricity()
        {
            CalculateRadialVector(closestPointToOtherObject, out double[] RadialVector, out _);
            double[] tangentialVector = CollisionTangentialVector.Last();
            eccentricity = RadialVector[0] * tangentialVector[0] + RadialVector[1] * tangentialVector[1];
        }

        internal void CalculateNormalAndTangentialVelocity()
        {
            double[] Velocity = translationalVelocity[0];
            double[] NormalVector = CollisionNormalVector.Last();
            double[] TangentialVector = CollisionTangentialVector.Last();
            PreCollisionVelocity = new double[] { Velocity[0] * NormalVector[0] + Velocity[1] * NormalVector[1], Velocity[0] * TangentialVector[0] + Velocity[1] * TangentialVector[1] };
        }

        /// <summary>
        /// This method performs the Neumaier algorithm form the sum of the entries of an array.
        /// </summary>
        /// <param name="ResultVariable">
        /// The variable where  the sum will be saved.
        /// </param>
        /// <param name="Summands">
        /// The array of summands
        /// </param>
        /// <param name="Length">
        /// The number of summands.
        /// </param>
        private double ForceTorqueSummationWithNeumaierArray(double ResultVariable, MultidimensionalArray Summands, double Length)
        {
            double sum = ResultVariable;
            double naiveSum;
            double c = 0.0;
            for (int i = 0; i < Length; i++)
            {
                naiveSum = sum + Summands[i, 0];
                if (Math.Abs(sum) >= Math.Abs(Summands[i, 0]))
                {
                    c += (sum - naiveSum) + Summands[i, 0];
                }
                else
                {
                    c += (Summands[i, 0] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            return sum + c;
        }
    }
}

