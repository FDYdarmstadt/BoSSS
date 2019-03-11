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

        }
        
        public Particle(int Dim, double[] startPos = null, double startAngl = 0.0) {
            
            m_HistoryLength = 4;
            m_Dim = Dim;

            #region Particle history
            // =============================   
            for (int i = 0; i < m_HistoryLength; i++) {
                positionAtIteration.Add(new double[Dim]);
                angleAtIteration.Add(new double());
                transVelocityAtIteration.Add(new double[Dim]);
                transAccelerationAtIteration.Add(new double[Dim]);
                rotationalVelocityAtIteration.Add(new double());
                rotationalAccelarationAtIteration.Add(new double());
                hydrodynForcesAtIteration.Add(new double[Dim]);
                hydrodynTorqueAtIteration.Add(new double());
            }
            for (int i = 0; i < 4; i++)
            {
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
            if (startPos == null)
            {
                if (Dim == 2)
                {
                    startPos = new double[] { 0.0, 0.0 };
                }
                else
                {
                    startPos = new double[] { 0.0, 0.0, 0.0 };
                }
            }
            positionAtTimestep[0] = startPos;
            positionAtTimestep[1] = startPos;
            //From degree to radiant
            angleAtTimestep[0] = startAngl * 2 * Math.PI / 360;
            angleAtTimestep[1] = startAngl * 2 * Math.PI / 360;

            UpdateLevelSetFunction();
            #endregion
        }
        #endregion

        #region Particle parameter
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
        public int underrelaxation_factor = 1;

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
        int m_HistoryLength;
        

        /// <summary>
        /// Length of history for time, velocity, position etc.
        /// </summary>
        [DataMember]
        public bool neglectAddedDamping = true;

        double[,] addedDampingTensorVV = new double[2, 2];
        double[,] addedDampingTensorVW = new double[2, 2];
        double[,] addedDampingTensorWV = new double[2, 2];
        double[,] addedDampingTensorWW = new double[2, 2];

        
        #endregion

        #region Geometric parameters
        /// <summary>
        /// Spatial Dimension of the particle 
        /// </summary>
        [DataMember]
        int m_Dim;
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
        
        /// <summary>
        /// The position (center of mass) of the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double[]> positionAtIteration = new List<double[]>();

        /// <summary>
        /// The position (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> positionAtTimestep = new List<double[]>();

        /// <summary>
        /// The angle (center of mass) of the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double> angleAtIteration = new List<double>();

        /// <summary>
        /// The angle (center of mass) of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> angleAtTimestep = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double[]> transVelocityAtIteration = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> transVelocityAtTimestep = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double> rotationalVelocityAtIteration = new List<double>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> rotationalVelocityAtTimestep = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double[]> transAccelerationAtIteration = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> transAccelerationAtTimestep = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double> rotationalAccelarationAtIteration = new List<double>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> rotationalAccelarationAtTimestep = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double[]> hydrodynForcesAtIteration = new List<double[]>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double[]> hydrodynForcesAtTimestep = new List<double[]>();
        
        /// <summary>
        /// The Torque acting on the particle in the current iteration.
        /// </summary>
        [DataMember]
        public List<double> hydrodynTorqueAtIteration = new List<double>();

        /// <summary>
        /// The Torque acting on the particle in the current time step.
        /// </summary>
        [DataMember]
        public List<double> hydrodynTorqueAtTimestep = new List<double>();

        /// <summary>
        /// Level set function describing the particle.
        /// </summary>       
        public Func<double[], double, double> phi_P;

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
        [DataMember]
        abstract public double Circumference_P {
            get;
        }

        /// <summary>
        /// Moment of inertia of the current particle.
        /// </summary>
        [DataMember]
        abstract public double MomentOfInertia_P {
            get;
        }
        #endregion
        #endregion

        #region Administrative tasks
        ParticleAuxillary Aux = new ParticleAuxillary();
        ParticlePhysics Physics = new ParticlePhysics();
        ParticleAddedDamping AddedDamping = new ParticleAddedDamping();
        ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
        #region obsolete
        ///// <summary>
        ///// Clean all Particle iteration histories until a certain length, obsolete?
        ///// </summary>
        ///// <param name="length"></param>
        //public void CleanHistoryIter() {
        //    if (positionAtIteration.Count > m_HistoryLength)
        //    {
        //        for (int j = positionAtIteration.Count; j > m_HistoryLength; j--)
        //        {
        //            int tempPos = j - 1;
        //            positionAtIteration.RemoveAt(tempPos);
        //            transVelocityAtIteration.RemoveAt(tempPos);
        //            hydrodynForcesAtIteration.RemoveAt(tempPos);
        //            rotationalVelocityAtIteration.RemoveAt(tempPos);
        //            hydrodynTorqueAtIteration.RemoveAt(tempPos);
        //        }
        //    }
        //}
        ///// <summary>
        ///// Clean all Particle histories until a certain length, obsolete?
        ///// </summary>
        ///// <param name="length"></param>
        //public void CleanHistory()
        //{
        //    if (positionAtTimestep.Count > 4)
        //    {
        //        for (int j = positionAtTimestep.Count; j > 4; j--)
        //        {
        //            int tempPos = j - 1;
        //            positionAtTimestep.RemoveAt(tempPos);
        //            angleAtTimestep.RemoveAt(tempPos);
        //            transVelocityAtTimestep.RemoveAt(tempPos);
        //            hydrodynForcesAtTimestep.RemoveAt(tempPos);
        //            rotationalVelocityAtTimestep.RemoveAt(tempPos);
        //            hydrodynTorqueAtTimestep.RemoveAt(tempPos);
        //        }
        //    }
        //}
        #endregion
        #endregion

        #region Move particle with current velocity
        /// <summary>
        /// Move particle with current velocity
        /// </summary>
        /// <param name="dt"></param>
        public void CalculateParticlePosition(double dt, double rho_Fluid) {

            if (m_Dim != 2 && m_Dim != 3)
                throw new NotSupportedException("Unknown particle dimension: m_Dim = " + m_Dim);
            for (int d = 0; d < m_Dim; d++)
            {
                //gravity[d] = 0;
                //gravity[1] = gravityVertical;
                //double massDifference = (particleDensity - fluidDensity) * (Area_P);
                //double tempForces = (hydrodynForcesAtIteration[0][d] + hydrodynForcesAtTimestep[1][d]) / 2;
                //tempPos[d] = particlePositionPerTimestep[1][d] + transVelocityAtTimestep[1][d] * dt + 0.5 * dt * dt * (tempForces + massDifference * gravity[d]) / Mass_P;
                positionAtIteration[0][d] = positionAtTimestep[1][d] + transVelocityAtTimestep[1][d] * dt + (transAccelerationAtTimestep[1][d] + transAccelerationAtIteration[0][d]) * dt.Pow2() / 4;
                if (double.IsNaN(positionAtIteration[0][d]) || double.IsInfinity(positionAtIteration[0][d]))
                    throw new ArithmeticException("Error trying to update particle position");
            }
            positionAtTimestep[0] = positionAtIteration[0];
        }

        public void CalculateParticleAngle(double dt)
        {
            angleAtIteration[0] = angleAtTimestep[1] + rotationalVelocityAtTimestep[1] + dt * (rotationalAccelarationAtTimestep[1] + rotationalAccelarationAtIteration[0]) / 2;
            if (double.IsNaN(angleAtIteration[0]) || double.IsInfinity(angleAtIteration[0]))
                throw new ArithmeticException("Error trying to update particle angle");
            angleAtTimestep[0] = angleAtIteration[0];
        }

        public void ResetParticlePosition()
        {
            if (iteration_counter_P == 0)
            {
                Aux.SaveMultidimValueOfLastTimestep(positionAtTimestep);
                Aux.SaveValueOfLastTimestep(angleAtTimestep);
            }

            Aux.SaveMultidimValueToList(positionAtIteration, positionAtIteration[0], 1);
            positionAtIteration[0] = positionAtTimestep[1];

            Aux.SaveValueToList(angleAtIteration, angleAtIteration[0], 1);
            angleAtIteration[0] = angleAtTimestep[1];

            UpdateLevelSetFunction();
        }
        #endregion
        
        abstract public void UpdateLevelSetFunction();

        public void PredictTranslationalAccelaration()
        {
            double[] temp = new double[m_Dim];
            for (int d = 0; d < m_Dim; d++)
            {
                temp[d] = (transAccelerationAtTimestep[0][d] + transAccelerationAtTimestep[1][d]) / 2;
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to predict particle acceleration");
            }
            Aux.SaveMultidimValueToList(transAccelerationAtIteration, temp);
            transAccelerationAtTimestep[0] = transAccelerationAtIteration[0];
        }

        public void CalculateTranslationalAcceleration(double dt, double fluidDensity, double addedDampingCoeff = 1)
        {
            double D1 = Mass_P + addedDampingCoeff * dt * addedDampingTensorVV[0, 0];
            double D2 = Mass_P + addedDampingCoeff * dt * addedDampingTensorVV[1, 1];
            double D3 = addedDampingCoeff * dt * addedDampingTensorVV[1, 0];
            double D4 = addedDampingCoeff * dt * addedDampingTensorVV[0, 1];
            double[] tempAcc = new double[2];

            tempAcc[0] = (hydrodynForcesAtIteration[0][0] - (D3 * D4 * hydrodynForcesAtIteration[0][0] - D1 * D2 * hydrodynForcesAtIteration[0][1]) / (D3 * D4 - D1 * D2)) / D1;
            if (double.IsNaN(tempAcc[0]) || double.IsInfinity(tempAcc[0]))
                throw new ArithmeticException("Error trying to calculate particle acceleration");

            tempAcc[1] = (D3 * hydrodynForcesAtIteration[0][0] - D1 * hydrodynForcesAtIteration[0][1]) / (D3 * D4 - D1 * D2);
            if (double.IsNaN(tempAcc[1]) || double.IsInfinity(tempAcc[1]))
                throw new ArithmeticException("Error trying to calculate particle acceleration");

            Aux.SaveMultidimValueToList(transAccelerationAtIteration, tempAcc);
            transAccelerationAtTimestep[0] = transAccelerationAtIteration[0];
        }

        public void PredictTranslationalVelocity()
        {
            double[] temp = new double[m_Dim];
            for (int d = 0; d < m_Dim; d++)
            {
                temp[d] = (transVelocityAtTimestep[0][d] + transVelocityAtTimestep[1][d]) / 2;
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to predict particle velocity");
            }
            Aux.SaveMultidimValueToList(transVelocityAtIteration, temp);
            transVelocityAtTimestep[0] = transVelocityAtIteration[0];
        }

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
                temp[d] = transVelocityAtTimestep[1][d] + (transAccelerationAtTimestep[1][d] + transAccelerationAtIteration[0][d]) * dt / 2;
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to calculate particle velocity");
            }
            Aux.SaveMultidimValueToList(transVelocityAtIteration, temp);
            transVelocityAtTimestep[0] = transVelocityAtIteration[0];
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
            //    f_vNew[i] = c_a[i] * (3 * transVelocityAtIteration[0][i] - 4 * transVelocityAtIteration[1][i] + transVelocityAtIteration[2][i]) / (2 * dt);
            //    f_vOld[i] = c_a[i] * (3 * transVelocityAtIteration[1][i] - 4 * transVelocityAtIteration[2][i] + transVelocityAtIteration[3][i]) / (2 * dt);
            //    tempForceNew[i] = (hydrodynForcesAtIteration[0][i] + massDifference * gravity[i]) * (c_u[i]) + f_vNew[i];
            //    tempForceOld[i] = (hydrodynForcesAtIteration[1][i] + massDifference * gravity[i]) * (c_u[i]) + f_vOld[i];
            //    temp[i] = transVelocityAtIteration[0][i] + (3 * tempForceNew[i] - tempForceOld[i]) * dt / 2;
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
            //        f_vTemp[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * temp[i] - 18 * transVelocityAtIteration[0][i] + 9 * transVelocityAtIteration[1][i] - 2 * transVelocityAtIteration[2][i]) / (8 * dt);
            //        f_vNew[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * transVelocityAtIteration[0][i] - 18 * transVelocityAtIteration[1][i] + 9 * transVelocityAtIteration[2][i] - 2 * transVelocityAtIteration[3][i]) / (6 * dt);
            //        f_vOld[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * transVelocityAtIteration[1][i] - 18 * transVelocityAtIteration[2][i] + 9 * transVelocityAtIteration[3][i] - 2 * transVelocityAtIteration[4][i]) / (6 * dt);
            //        tempForces[i] = (hydrodynForcesAtIteration[0][i] + massDifference * gravity[i]) * (c_u[i]) + f_vTemp[i];
            //        tempForceNew[i] = (hydrodynForcesAtIteration[1][i] + massDifference * gravity[i]) * (c_u[i]) + f_vNew[i];
            //        tempForceOld[i] = (hydrodynForcesAtIteration[2][i] + massDifference * gravity[i]) * (c_u[i]) + f_vOld[i];
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

        public void PredictAngularAcceleration()
        {
            double temp = (rotationalAccelarationAtTimestep[0] + rotationalAccelarationAtTimestep[1]) / 2;
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to predict particle angluar acceleration");
            Aux.SaveValueToList(rotationalAccelarationAtIteration, temp);
            rotationalAccelarationAtTimestep[0] = rotationalAccelarationAtIteration[0];
        }

        public void CalculateAngularAcceleration(double dt, double addedDampingCoeff = 1)
        {
            double MomentofInertia_m = MomentOfInertia_P + addedDampingCoeff * dt * addedDampingTensorVV[0, 0];
            Aux.SaveValueToList(rotationalAccelarationAtIteration, hydrodynTorqueAtIteration[0] / MomentofInertia_m);
        }

        public void PredictAngularVelocity()
        {
            double temp = (rotationalVelocityAtTimestep[0] + rotationalVelocityAtTimestep[1]) / 2;
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to predict particle angluar velocity");
            Aux.SaveValueToList(rotationalVelocityAtIteration, temp);
            rotationalVelocityAtTimestep[0] = rotationalVelocityAtIteration[0];
        }
        
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
            //    rotationalVelocityAtIteration.Insert(0, 0.0);
            //    rotationalVelocityAtIteration.Remove(rotationalVelocityAtIteration.Last());
            //    return;
            //}

            double newAngularVelocity = 0;
            //double oldAngularVelocity = new double();
            double subtimestep;
            noOfSubtimesteps = 1;
            subtimestep = dt / noOfSubtimesteps;
            for (int i = 1; i <= noOfSubtimesteps; i++)
            {
                newAngularVelocity = rotationalVelocityAtTimestep[1] + dt * (rotationalAccelarationAtTimestep[1] + rotationalAccelarationAtIteration[0]) / 2;
            }
            if (double.IsNaN(newAngularVelocity) || double.IsInfinity(newAngularVelocity))
                throw new ArithmeticException("Error trying to calculate particle angluar velocity");
            Aux.SaveValueToList(rotationalVelocityAtIteration, newAngularVelocity);
            rotationalVelocityAtTimestep[0] = rotationalVelocityAtIteration[0];
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
        abstract protected double averageDistance { get; }
        
        public void CalculateDampingTensors(LevelSetTracker LsTrk, double muA, double rhoA, double dt)
        {
            addedDampingTensorVV = AddedDamping.IntegrationOverLevelSet(0, LsTrk, muA, rhoA, dt, positionAtIteration[0], cutCells_P(LsTrk), neglectAddedDamping);
            addedDampingTensorVW = AddedDamping.IntegrationOverLevelSet(1, LsTrk, muA, rhoA, dt, positionAtIteration[0], cutCells_P(LsTrk), neglectAddedDamping);
            addedDampingTensorWV = AddedDamping.IntegrationOverLevelSet(2, LsTrk, muA, rhoA, dt, positionAtIteration[0], cutCells_P(LsTrk), neglectAddedDamping);
            addedDampingTensorWW = AddedDamping.IntegrationOverLevelSet(3, LsTrk, muA, rhoA, dt, positionAtIteration[0], cutCells_P(LsTrk), neglectAddedDamping);
        }

        public void UpdateDampingTensors()
        {
            addedDampingTensorVV = AddedDamping.RotateTensor(angleAtIteration[0], addedDampingTensorVV);
            addedDampingTensorVW = AddedDamping.RotateTensor(angleAtIteration[0], addedDampingTensorVW);
            addedDampingTensorWV = AddedDamping.RotateTensor(angleAtIteration[0], addedDampingTensorWV);
            addedDampingTensorWW = AddedDamping.RotateTensor(angleAtIteration[0], addedDampingTensorWW);
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

            int spatialDim = LsTrk.GridDat.SpatialDimension;

            var UA = U.ToArray();

            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            ConventionalDGField pA = null;
            pA = P;

            #region Force
            double[] Forces = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                ScalarFunctionEx ErrFunc = delegate (int j0, int NumberOfCells, NodeSet Ns, MultidimensionalArray result) {
                    int NumberOfNodes = result.GetLength(1); 
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(NumberOfCells, NumberOfNodes, spatialDim, spatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(NumberOfCells, NumberOfNodes);
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, NumberOfCells);
                    for (int i = 0; i < spatialDim; i++) {
                        UA[i].EvaluateGradient(j0, NumberOfCells, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    pA.Evaluate(j0, NumberOfCells, Ns, pARes);
                    for (int j = 0; j < NumberOfCells; j++)
                    {
                        for (int k = 0; k < NumberOfNodes; k++)
                        {
                            result[j, k] = Physics.CalculateStressTensor(Grad_UARes, pARes, Normals, muA, k, j, m_Dim, d);
                        }
                    }
                };
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, this.cutCells_P(LsTrk));
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
            #endregion

            #region Torque
            double Torque = 0;
            ScalarFunctionEx ErrFunc2 = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
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
                        result[j, k] = Physics.CalculateTorqueFromStressTensor2D(Grad_UARes, pARes, Normals, tempArray, muA, k, j, positionAtIteration[0]);
                    }
                }
            };
            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, this.cutCells_P(LsTrk));
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
            Forces[0] = Forces[0] + addedDampingTensorVV[0, 0] * beta * transAccelerationAtIteration[0][0] * dt + addedDampingTensorVV[1, 0] * beta * transAccelerationAtIteration[0][1] * dt;
            Forces[1] = Forces[1] + addedDampingTensorVV[0, 1] * beta * transAccelerationAtIteration[0][0] * dt + addedDampingTensorVV[1, 1] * beta * transAccelerationAtIteration[0][1] * dt + (particleDensity - fluidDensity) * Area_P * gravityVertical;
            Torque = Torque - beta * dt * addedDampingTensorWW[0, 0] * rotationalAccelarationAtIteration[0];
            if (iteration_counter_P == 0)
            {
                Console.WriteLine("First iteration of the current timestep, all relaxation factors are set to 1");
                for (int d = 0; d < spatialDim; d++)
                {
                    for (int t = 0; t < m_HistoryLength; t++)
                    {
                        hydrodynForcesAtIteration[t][d] = hydrodynForcesAtTimestep[1][d];
                        hydrodynTorqueAtIteration[t] = hydrodynTorqueAtTimestep[1];
                    }
                    if (Math.Abs(Forces[d]) < forceAndTorque_convergence * 1e-2 && ClearSmallValues == true)
                    {
                        Forces[d] = 0;
                    }
                }
                if (Math.Abs(Torque) < forceAndTorque_convergence * 1e-2 && ClearSmallValues == true)
                {
                    Torque = 0;
                }
            }
            else
            {
                double[] RelaxatedForceAndTorque = Underrelaxation.RelaxatedForcesAndTorque(Forces, Torque, hydrodynForcesAtIteration[0], hydrodynTorqueAtIteration[0], forceAndTorque_convergence, underrelaxation_factor, ClearSmallValues, AddaptiveUnderrelaxation);
                for (int d = 0; d < m_Dim; d++)
                {
                    Forces[d] = RelaxatedForceAndTorque[d];
                }
                Torque = RelaxatedForceAndTorque[m_Dim];
            }
            Aux.SaveMultidimValueToList(hydrodynForcesAtIteration, Forces);
            Aux.SaveValueToList(hydrodynTorqueAtIteration, Torque);
            hydrodynForcesAtTimestep[0] = hydrodynForcesAtIteration[0];
            hydrodynTorqueAtTimestep[0] = hydrodynTorqueAtIteration[0];
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
        abstract public CellMask cutCells_P(LevelSetTracker LsTrk); 

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        abstract public bool Contains(double[] point, LevelSetTracker LsTrk); 
        #endregion
        
    }
}

