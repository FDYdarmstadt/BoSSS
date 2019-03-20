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

namespace BoSSS.Application.FSI_Solver
{

    /// <summary>
    /// Particle properties (for disk shape and spherical particles only).
    /// </summary>
    [DataContract]
    [Serializable]
    abstract public class Particle : ICloneable {

        #region particle init
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        protected Particle()
        {
            // noop

        }
        
        public Particle(int Dim, double[] startPos = null, double startAngl = 0.0) {
            
            m_HistoryLength = 4;
            m_Dim = Dim;

            #region Particle history
            // =============================   
            for (int i = 0; i < 4; i++) {
                positionAtTimestep.Add(new double[Dim]);
                angleAtTimestep.Add(new double());
                transVelocityAtTimestep.Add(new double[Dim]);
                transAccelerationAtTimestep.Add(new double[Dim]);
                rotationalVelocityAtTimestep.Add(new double());
                rotationalAccelarationAtTimestep.Add(new double());
                hydrodynForcesAtTimestep.Add(new double[Dim]);
                hydrodynTorqueAtTimestep.Add(new double());
            }
            #endregion

            #region Initial values
            // ============================= 
            if (startPos == null) {
                startPos = new double[Dim];
            }
            positionAtTimestep[0] = startPos;
            positionAtTimestep[1] = startPos;
            //From degree to radiant
            angleAtTimestep[0] = StartingAngle = startAngl * 2 * Math.PI / 360;
            angleAtTimestep[1] = startAngl * 2 * Math.PI / 360;

            //UpdateLevelSetFunction();
            #endregion
        }
        #endregion


        #region Collision parameters
        /// <summary>
        /// Check whether any particles is collided with another particle
        /// </summary>
        public bool[] m_collidedWithParticle;

        /// <summary>
        /// Check whether any particles is collided with the wall
        /// </summary>
        public bool[] m_collidedWithWall;

        public double[][] m_closeInterfacePointTo;

        /// <summary>
        /// Skip calculation of hydrodynamic force and Torque if particles are too close
        /// </summary>
        [DataMember]
        public bool skipForceIntegration = false;
        #endregion

        #region Iteration parameters
        /// <summary>
        /// Number of iterations
        /// </summary>
        [DataMember]
        public int iteration_counter_P = 0;

        /// <summary>
        /// Constant Forces and Torque underrelaxation?
        /// </summary>
        [DataMember]
        public bool AddaptiveUnderrelaxation = false;

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
        public bool ClearSmallValues = false;
        #endregion

        #region Misc parameters
        /// <summary>
        /// Length of history for time, velocity, position etc.
        /// </summary>
        readonly int m_HistoryLength;
        

        /// <summary>
        /// Length of history for time, velocity, position etc.
        /// </summary>
        [DataMember]
        public bool neglectAddedDamping = true;

        double[,] AddedDampingTensor = new double[6, 6];
        
        #endregion

        #region Geometric parameters
        /// <summary>
        /// Spatial Dimension of the particle 
        /// </summary>
        [DataMember]
        readonly int m_Dim;
        #endregion

        #region Virtual force model parameter
        ///// <summary>
        ///// needed for second velocity model
        ///// </summary>
        //public double C_v = 0.5;

        ///// <summary>
        ///// needed for second velocity model, obsolete?
        ///// </summary>
        //public double velResidual_ConvergenceCriterion = 1e-6;

        ///// <summary>
        ///// needed for second velocity model, obsolete?
        ///// </summary>
        //public double MaxParticleVelIterations = 10000;

        //private int vel_iteration_counter;
        #endregion

        #region Physical parameters
        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double particleDensity;

        double asdf = 0;

        ///// <summary>
        ///// The position (center of mass) of the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double[]> positionAtTimestep = new List<double[]>();

        /// <summary>
        /// The position (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> positionAtTimestep = new List<double[]>();

        ///// <summary>
        ///// The angle (center of mass) of the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double> angleAtTimestep = new List<double>();

        /// <summary>
        /// The angle (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> angleAtTimestep = new List<double>();

        ///// <summary>
        /// The angle (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public double StartingAngle = new double();
        
        ///// The translational velocity of the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double[]> transVelocityAtTimestep = new List<double[]>();

        /// <summary>

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> transVelocityAtTimestep = new List<double[]>();

        ///// <summary>
        ///// The angular velocity of the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double> rotationalVelocityAtTimestep = new List<double>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> rotationalVelocityAtTimestep = new List<double>();

        ///// <summary>
        ///// The translational velocity of the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double[]> transAccelerationAtTimestep = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> transAccelerationAtTimestep = new List<double[]>();

        ///// <summary>
        ///// The angular velocity of the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double> rotationalAccelarationAtTimestep = new List<double>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> rotationalAccelarationAtTimestep = new List<double>();

        ///// <summary>
        ///// The force acting on the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double[]> hydrodynForcesAtTimestep = new List<double[]>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> hydrodynForcesAtTimestep = new List<double[]>();
        
        ///// <summary>
        ///// The Torque acting on the particle in the current iteration.
        ///// </summary>
        //[DataMember]
        //public List<double> hydrodynTorqueAtTimestep = new List<double>();

        /// <summary>
        /// The Torque acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> hydrodynTorqueAtTimestep = new List<double>();

        /// <summary>
        /// Level set function describing the particle.
        /// </summary>       
        public abstract double phi_P(double[] X);

        /// <summary>
        /// Sets the gravity in vertical direction, default is 0.0
        /// </summary>
        [DataMember]
        public double gravityVertical = 0.0;

        /// <summary>
        /// Set true if the particle should be an active particle, i.e. self driven
        /// </summary>
        [DataMember]
        public bool activeParticle = false;

        /// <summary>
        /// Convergence criterion for the calculation of the Forces and Torque
        /// </summary>
        [DataMember]
        public double forceAndTorque_convergence = 1e-8;

        /// <summary>
        /// Active stress on the current particle.
        /// </summary>
        public double active_stress_P;
        
        /// <summary>
        /// Area of the current particle.
        /// </summary>
        abstract public double Area_P {
            get;
        }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        public double Mass_P {
            get {
                double a = Area_P;
                if (a <= 0.0 || double.IsNaN(a) || double.IsInfinity(a))
                    throw new ArithmeticException("Particle volume/area is " + a);
                return Area_P * particleDensity;
            }
        }

        /// <summary>
        /// Circumference of the current particle.
        /// </summary>
        abstract public double Circumference_P {
            get;
        }

        /// <summary>
        /// Moment of inertia of the current particle.
        /// </summary>
        abstract public double MomentOfInertia_P {
            get;
        }
        #endregion

        #region Administrative tasks

        [NonSerialized]
        ParticleAuxillary Aux = new ParticleAuxillary();
        [NonSerialized]
        ParticleForceIntegration ForceIntegration = new ParticleForceIntegration();
        [NonSerialized]
        ParticleAddedDamping AddedDamping = new ParticleAddedDamping();
        [NonSerialized]
        ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
        [NonSerialized]
        ParticleAcceleration Acceleration = new ParticleAcceleration();
        #endregion
        
        /// <summary>
        /// Move particle with current velocity
        /// </summary>
        /// <param name="dt"></param>
        public void CalculateParticlePosition(double dt, double rho_Fluid)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(positionAtTimestep);
            }

            if (m_Dim != 2 && m_Dim != 3)
                throw new NotSupportedException("Unknown particle dimension: m_Dim = " + m_Dim);

            for (int d = 0; d < m_Dim; d++)
            {
                positionAtTimestep[0][d] = positionAtTimestep[1][d] + transVelocityAtTimestep[1][d] * dt + (transAccelerationAtTimestep[1][d] + transAccelerationAtTimestep[0][d]) * dt.Pow2() / 4;
                if (double.IsNaN(positionAtTimestep[0][d]) || double.IsInfinity(positionAtTimestep[0][d]))
                    throw new ArithmeticException("Error trying to update particle position");
            }
        }

        public void CalculateParticleAngle(double dt)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveValueOfLastTimestep(angleAtTimestep);
            }

            double tempAngle = angleAtTimestep[1] + rotationalVelocityAtTimestep[1] * dt + dt.Pow2() * (rotationalAccelarationAtTimestep[1] + rotationalAccelarationAtTimestep[0]) / 4;
            if (double.IsNaN(tempAngle) || double.IsInfinity(tempAngle))
                throw new ArithmeticException("Error trying to update particle angle");
            angleAtTimestep[0] = tempAngle;
        }
        
        //abstract public void UpdateLevelSetFunction();

        //public void PredictTranslationalAccelaration()
        //{
        //    if (iteration_counter_P == 0)
        //    {
        //        Aux.SaveMultidimValueOfLastTimestep(transAccelerationAtTimestep);
        //    }
        //    double[] temp = new double[m_Dim];
        //    for (int d = 0; d < m_Dim; d++)
        //    {
        //        temp[d] = (0.8 * transAccelerationAtTimestep[1][d] + 0.2 * transAccelerationAtTimestep[2][d]);
        //        if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
        //            throw new ArithmeticException("Error trying to predict particle acceleration");
        //    }
        //    Aux.SaveMultidimValueToList(transAccelerationAtTimestep, temp);
        //    transAccelerationAtTimestep[0] = transAccelerationAtTimestep[0];
        //}

        public void CalculateAcceleration(double dt, double fluidDensity, double addedDampingCoeff = 1)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(transAccelerationAtTimestep);
                Aux.SaveValueOfLastTimestep(rotationalAccelarationAtTimestep);
            }

            double[] tempAccTrans = new double[2];
            double[,] CoefficientMatrix = Acceleration.CalculateCoefficients(AddedDampingTensor, Mass_P, MomentOfInertia_P, dt);
            double Denominator = Acceleration.CalculateDenominator(CoefficientMatrix);

            transAccelerationAtTimestep[0] = Acceleration.Translational(CoefficientMatrix, Denominator, hydrodynForcesAtTimestep[0], hydrodynTorqueAtTimestep[0]);
            for (int i = 0; i< m_Dim; i++)
            {
                if (double.IsNaN(transAccelerationAtTimestep[0][i]) || double.IsInfinity(transAccelerationAtTimestep[0][i]))
                    throw new ArithmeticException("Error trying to calculate particle acceleration");
            }

            rotationalAccelarationAtTimestep[0] = Acceleration.Rotational(CoefficientMatrix, Denominator, hydrodynForcesAtTimestep[0], hydrodynTorqueAtTimestep[0]);
            if (double.IsNaN(rotationalAccelarationAtTimestep[0]) || double.IsInfinity(rotationalAccelarationAtTimestep[0]))
                throw new ArithmeticException("Error trying to calculate particle rotational acceleration");
        }

        //public void PredictTranslationalVelocity()
        //{
        //    if (iteration_counter_P == 0)
        //    {
        //        Aux.SaveMultidimValueOfLastTimestep(transVelocityAtTimestep);
        //    }
        //    double[] temp = new double[m_Dim];
        //    for (int d = 0; d < m_Dim; d++)
        //    {
        //        temp[d] = (0.8 * transVelocityAtTimestep[1][d] + 0.2 * transVelocityAtTimestep[2][d]);
        //        if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
        //            throw new ArithmeticException("Error trying to predict particle velocity");
        //    }
        //    Aux.SaveMultidimValueToList(transVelocityAtTimestep, temp);
        //    transVelocityAtTimestep[0] = transVelocityAtTimestep[0];
        //}

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public void CalculateTranslationalVelocity(double dt, double fluidDensity)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(transVelocityAtTimestep);
            }

            double[] temp = new double[m_Dim];
            for (int d = 0; d < m_Dim; d++)
            {
                transVelocityAtTimestep[0][d] = transVelocityAtTimestep[1][d] + (transAccelerationAtTimestep[1][d] + transAccelerationAtTimestep[0][d]) * dt / 2;
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to calculate particle velocity");
            }
            return;
        }

        public void VirtualForceModel(double dt, double fluidDensity)
        {


            double[] temp = new double[2];
            double[] old_temp = new double[2];
            double[] tempForces = new double[2];
            double massDifference = (particleDensity - fluidDensity) * (Area_P);

            #region virtual force model
            // Virtual force model (Schwarz et al. - 2015 A temporal discretization scheme to compute the motion of light particles in viscous flows by an immersed boundary")
            // =============================
            double[] f_vTemp = new double[2];
            double[] f_vNew = new double[2];
            double[] f_vOld = new double[2];
            double[] k_1 = new double[2];
            double[] k_2 = new double[2];
            double[] k_3 = new double[2];
            double[] C_v_mod = new double[2];
            double[] tempForceNew = new double[2];
            double[] tempForceOld = new double[2];
            //C_v_mod[0] = C_v;
            //C_v_mod[1] = C_v;
            double[] c_a = new double[2];
            double[] c_u = new double[2];
            //double vel_iteration_counter = 0;
            double[] test = new double[2];
            // 2nd order Adam Bashford
            //for (int i = 0; i < 2; i++)
            //{
            //    dt = 1e-3;
            //    C_v_mod[i] = 0.1;
            //    c_a[i] = (C_v_mod[i] * fluidDensity) / (particleDensity + C_v_mod[i] * fluidDensity);
            //    c_u[i] = 1 / (Area_P * (particleDensity + C_v_mod[i] * particleDensity));
            //    f_vNew[i] = c_a[i] * (3 * transVelocityAtTimestep[0][i] - 4 * transVelocityAtTimestep[1][i] + transVelocityAtTimestep[2][i]) / (2 * dt);
            //    f_vOld[i] = c_a[i] * (3 * transVelocityAtTimestep[1][i] - 4 * transVelocityAtTimestep[2][i] + transVelocityAtTimestep[3][i]) / (2 * dt);
            //    tempForceNew[i] = (hydrodynForcesAtTimestep[0][i] + massDifference * gravity[i]) * (c_u[i]) + f_vNew[i];
            //    tempForceOld[i] = (hydrodynForcesAtTimestep[1][i] + massDifference * gravity[i]) * (c_u[i]) + f_vOld[i];
            //    temp[i] = transVelocityAtTimestep[0][i] + (3 * tempForceNew[i] - tempForceOld[i]) * dt / 2;
            //}

            // implicit Adams Moulton (modified)
            //for (double velResidual = 1; velResidual > velResidual_ConvergenceCriterion;)
            //{
            //    dt = 1e-3;
            //    for (int i = 0; i < 2; i++)
            //    {
            //        gravity[0] = 0;
            //        if (includeGravity == true)
            //        {
            //            gravity[1] = -9.81;
            //        }
            //        C_v_mod[i] = 300;// * Math.Abs(forces_P[0][i] / (forces_P[0][i] + forces_P[1][i] + 1e-30));
            //        c_a[i] = (C_v_mod[i] * fluidDensity) / (particleDensity + C_v_mod[i] * fluidDensity);
            //        c_u[i] = 1 / (Area_P * (particleDensity + C_v_mod[i] * particleDensity));
            //        f_vTemp[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * temp[i] - 18 * transVelocityAtTimestep[0][i] + 9 * transVelocityAtTimestep[1][i] - 2 * transVelocityAtTimestep[2][i]) / (8 * dt);
            //        f_vNew[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * transVelocityAtTimestep[0][i] - 18 * transVelocityAtTimestep[1][i] + 9 * transVelocityAtTimestep[2][i] - 2 * transVelocityAtTimestep[3][i]) / (6 * dt);
            //        f_vOld[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * transVelocityAtTimestep[1][i] - 18 * transVelocityAtTimestep[2][i] + 9 * transVelocityAtTimestep[3][i] - 2 * transVelocityAtTimestep[4][i]) / (6 * dt);
            //        tempForces[i] = (hydrodynForcesAtTimestep[0][i] + massDifference * gravity[i]) * (c_u[i]) + f_vTemp[i];
            //        tempForceNew[i] = (hydrodynForcesAtTimestep[1][i] + massDifference * gravity[i]) * (c_u[i]) + f_vNew[i];
            //        tempForceOld[i] = (hydrodynForcesAtTimestep[2][i] + massDifference * gravity[i]) * (c_u[i]) + f_vOld[i];
            //        old_temp[i] = temp[i];
            //        temp[i] = previous_vel[i] + (1 * tempForces[i] + 4 * tempForceNew[i] + 1 * tempForceOld[i]) * dt / 6;
            //    }
            //    vel_iteration_counter += 1;
            //    if (vel_iteration_counter == MaxParticleVelIterations)
            //    {
            //        throw new ApplicationException("no convergence in particle velocity calculation");
            //    }
            //    velResidual = Math.Sqrt((temp[0] - old_temp[0]).Pow2() + (temp[1] - old_temp[1]).Pow2());

            //    Console.WriteLine("Current velResidual:  " + velResidual);
            //}
            //Console.WriteLine("Number of Iterations for translational velocity calculation:  " + vel_iteration_counter);
            //Console.WriteLine("C_v_mod:  " + C_v_mod[0]);
            #endregion
        }//unused

        //public void PredictAngularAcceleration() {
        //    if (iteration_counter_P == 0) {
        //        Aux.SaveValueOfLastTimestep(rotationalAccelarationAtTimestep);
        //    }
        //    double temp = (0.8 * rotationalAccelarationAtTimestep[1] + 0.2 * rotationalAccelarationAtTimestep[2]);
        //    if (double.IsNaN(temp) || double.IsInfinity(temp))
        //        throw new ArithmeticException("Error trying to predict particle angluar acceleration");
        //    Aux.SaveValueToList(rotationalAccelarationAtTimestep, temp);
        //}

        //public void CalculateAngularAcceleration(double dt, double addedDampingCoeff = 1)
        //{
        //    double MomentofInertia_m = MomentOfInertia_P;// + addedDampingCoeff * dt * AddedDampingTensor[0, 0];
        //    Aux.SaveValueToList(rotationalAccelarationAtIteration, hydrodynTorqueAtIteration[0] / MomentofInertia_m);
        //}

        //public void PredictAngularVelocity() {
        //    if (iteration_counter_P == 0) {
        //        Aux.SaveValueOfLastTimestep(rotationalVelocityAtTimestep);
        //    }
        //    double temp = (0.8 * rotationalVelocityAtTimestep[1] + 0.2 * rotationalVelocityAtTimestep[2]);
        //    if (double.IsNaN(temp) || double.IsInfinity(temp))
        //        throw new ArithmeticException("Error trying to predict particle angluar velocity");
        //    Aux.SaveValueToList(rotationalVelocityAtTimestep, temp);
        //}

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public void CalculateAngularVelocity(double dt, int noOfSubtimesteps = 1)
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveValueOfLastTimestep(rotationalVelocityAtTimestep);
            }

            // no rotation
            // =============================
            //if (includeRotation == false) {
            //    rotationalVelocityAtTimestep.Insert(0, 0.0);
            //    rotationalVelocityAtTimestep.Remove(rotationalVelocityAtTimestep.Last());
            //    return;
            //}

            double newAngularVelocity = 0;
            //double oldAngularVelocity = new double();
            double subtimestep;
            noOfSubtimesteps = 1;
            subtimestep = dt / noOfSubtimesteps;
            for (int i = 1; i <= noOfSubtimesteps; i++) {
                newAngularVelocity = rotationalVelocityAtTimestep[1] + dt * (rotationalAccelarationAtTimestep[1] + rotationalAccelarationAtTimestep[0]) / 2;
            }
            if (double.IsNaN(newAngularVelocity) || double.IsInfinity(newAngularVelocity))
                throw new ArithmeticException("Error trying to calculate particle angluar velocity");
        }
        
        #region Cloning
        /// <summary>
        /// clone
        /// </summary>
        virtual public object Clone() {
            throw new NotImplementedException("Currently cloning of a particle is not available");
        }
        #endregion

        /// <summary>
        /// some length scale 
        /// </summary>
        abstract protected double AverageDistance { get; }

        public void CalculateDampingTensor(LevelSetTracker LsTrk, double muA, double rhoA, double dt)
        {
            AddedDampingTensor = AddedDamping.IntegrationOverLevelSet(LsTrk, muA, rhoA, dt, positionAtTimestep[0], CutCells_P(LsTrk));
        }

        public void UpdateDampingTensors()
        {
            AddedDampingTensor = AddedDamping.RotateTensor(angleAtTimestep[0], StartingAngle, AddedDampingTensor);
        }

        #region Update forces and torque
        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, double muA, double dt, double fluidDensity) {

            if (skipForceIntegration) {
                skipForceIntegration = false;
                return;
            }
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(hydrodynForcesAtTimestep);
                Aux.SaveValueOfLastTimestep(hydrodynTorqueAtTimestep);
            }

            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            int spatialDim = LsTrk.GridDat.SpatialDimension;
            #region Force
            double[] Forces = new double[spatialDim];
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = null;
            pA = P;
            for (int d = 0; d < spatialDim; d++)
            {
                void ErrFunc(int j0, int NumberOfCells, NodeSet Ns, MultidimensionalArray result)
                {
                    
                    int NumberOfNodes = result.GetLength(1);
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(NumberOfCells, NumberOfNodes, spatialDim, spatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(NumberOfCells, NumberOfNodes);
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, NumberOfCells);
                    for (int i = 0; i < spatialDim; i++) {
                        UA[i].EvaluateGradient(j0, NumberOfCells, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    pA.Evaluate(j0, NumberOfCells, Ns, pARes);
                    for (int j = 0; j < NumberOfCells; j++) {
                        for (int k = 0; k < NumberOfNodes; k++) {
                            result[j, k] = ForceIntegration.CalculateStressTensor(Grad_UARes, pARes, Normals, muA, k, j, m_Dim, d);
                        }
                    }
                }
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, this.CutCells_P(LsTrk));
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, RequiredOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        Forces[d] = Aux.ForceTorqueSummationWithNeumaierArray(Forces[d], ResultsOfIntegration, Length);
                    }
                ).Execute();
            }

            // add gravity
            {
                 Forces[1] =+ (particleDensity - fluidDensity) * Area_P * gravityVertical;
            }

            #endregion

            #region Torque
            double Torque = 0;
            void ErrFunc2(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, spatialDim, spatialDim); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                // Evaluate tangential velocity to level-set surface
                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                for (int i = 0; i < spatialDim; i++) {
                    UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }
                MultidimensionalArray tempArray = Ns.CloneAs();
                LsTrk.GridDat.TransformLocal2Global(Ns, tempArray, j0);
                pA.Evaluate(j0, Len, Ns, pARes);
                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {
                        result[j, k] = ForceIntegration.CalculateTorqueFromStressTensor2D(Grad_UARes, pARes, Normals, tempArray, muA, k, j, positionAtTimestep[0]);
                    }
                }
            }
            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, this.CutCells_P(LsTrk));
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs2.Compile(LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    Torque = Aux.ForceTorqueSummationWithNeumaierArray(Torque, ResultsOfIntegration, Length);
                }

            ).Execute();
            #endregion
            double beta = 1;
            if (neglectAddedDamping == false) {
                Forces[0] = Forces[0] + beta * dt * (AddedDampingTensor[0, 0] * transAccelerationAtTimestep[0][0] + AddedDampingTensor[1, 0] * transAccelerationAtTimestep[0][1] + AddedDampingTensor[0, 2] * rotationalAccelarationAtTimestep[0]);
                Forces[1] = Forces[1] + beta * dt * (AddedDampingTensor[0, 1] * transAccelerationAtTimestep[0][0] + AddedDampingTensor[1, 1] * transAccelerationAtTimestep[0][1] + AddedDampingTensor[1, 2] * rotationalAccelarationAtTimestep[0]);
                Torque = Torque + beta * dt * (AddedDampingTensor[2, 0] * transAccelerationAtTimestep[0][0] + AddedDampingTensor[2, 1] * transAccelerationAtTimestep[0][1] + AddedDampingTensor[2, 2] * rotationalAccelarationAtTimestep[0]);
            }
            if (iteration_counter_P == 0 && asdf == 0) {
                Console.WriteLine("First iteration of the current timestep, all relaxation factors are set to 1");
                for (int d = 0; d < spatialDim; d++) {
                    for (int t = 0; t < m_HistoryLength; t++) {
                        hydrodynForcesAtTimestep[t][d] = hydrodynForcesAtTimestep[1][d];
                        hydrodynTorqueAtTimestep[t] = hydrodynTorqueAtTimestep[1];
                    }
                    if (Math.Abs(Forces[d]) < forceAndTorque_convergence * 1e-2 && ClearSmallValues == true) {
                        Forces[d] = 0;
                    }
                }
                if (Math.Abs(Torque) < forceAndTorque_convergence * 1e-2 && ClearSmallValues == true) {
                    Torque = 0;
                }
                asdf = 100;
            }
            else if (iteration_counter_P == 100)
            {
                Console.WriteLine("No convergence after 100 iterations, I will try to restart");
                double ForceSummation = 0;
                for (int d = 0; d < spatialDim; d++)
                {
                    ForceSummation += hydrodynForcesAtTimestep[1][d]; 
                }
                hydrodynForcesAtTimestep[0][0] = ForceSummation * Math.Cos(angleAtTimestep[1]);
                hydrodynForcesAtTimestep[0][1] = ForceSummation * Math.Sin(angleAtTimestep[1]);
                hydrodynTorqueAtTimestep[0] = 0;
            }
            else
            {
                double[] RelaxatedForceAndTorque = Underrelaxation.RelaxatedForcesAndTorque(Forces, Torque, hydrodynForcesAtTimestep[0], hydrodynTorqueAtTimestep[0], forceAndTorque_convergence, underrelaxation_factor, ClearSmallValues, AddaptiveUnderrelaxation);
                for (int d = 0; d < m_Dim; d++) {
                    hydrodynForcesAtTimestep[0][d] = RelaxatedForceAndTorque[d];
                }
                hydrodynTorqueAtTimestep[0] = RelaxatedForceAndTorque[m_Dim];
            }
        }
        #endregion

        #region Particle reynolds number
        /// <summary>
        /// Calculating the particle reynolds number according to paper Turek and testcase ParticleUnderGravity
        /// </summary>
        abstract public double ComputeParticleRe(double ViscosityFluid);
        #endregion

        #region Cut cells
        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        abstract public CellMask CutCells_P(LevelSetTracker LsTrk);

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        abstract public bool Contains(double[] point, LevelSetTracker LsTrk); 
        #endregion
        
    }
}

