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
        
        public Particle(int Dim, int HistoryLength, double[] startPos = null, double startAngl = 0.0) {
            
            m_HistoryLength = HistoryLength;
            m_Dim = Dim;

            #region Particle history
            // =============================   
            for (int i = 0; i < HistoryLength; i++) {
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
        public bool underrelaxationFT_constant = true;

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
        ParticleAuxillary aux = new ParticleAuxillary();
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
                aux.SaveMultidimValueOfLastTimestep(positionAtTimestep);
                aux.SaveValueOfLastTimestep(angleAtTimestep);
            }

            aux.SaveMultidimValueToList(positionAtIteration, positionAtIteration[0], 1);
            positionAtIteration[0] = positionAtTimestep[1];

            aux.SaveValueToList(angleAtIteration, angleAtIteration[0], 1);
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
                temp[d] = 2 * transAccelerationAtTimestep[0][d] - transAccelerationAtTimestep[1][d];
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to predict particle acceleration");
            }
            aux.SaveMultidimValueToList(transAccelerationAtIteration, temp);
            transAccelerationAtTimestep[0] = transAccelerationAtIteration[0];
        }

        public void CalculateTranslationalAcceleration(double dt, double fluidDensity, double addedDampingCoeff = 1)
        {
            double D1 = Mass_P + addedDampingCoeff * dt * addedDampingTensorVV[0, 0];
            double D2 = Mass_P + addedDampingCoeff * dt * addedDampingTensorVV[1, 1];
            double D3 = addedDampingCoeff * dt * addedDampingTensorVV[1, 0];
            double D4 = addedDampingCoeff * dt * addedDampingTensorVV[0, 1];
            double[] tempAcc = new double[2];

            tempAcc[0] = ((hydrodynForcesAtIteration[0][0] + hydrodynForcesAtTimestep[1][0]) / 2 - D4 * (D3 * hydrodynForcesAtIteration[0][0] - D1 * hydrodynForcesAtIteration[0][1]) / (D3 * D4 - D1 * D2)) / D1;
            if (double.IsNaN(tempAcc[0]) || double.IsInfinity(tempAcc[0]))
                throw new ArithmeticException("Error trying to calculate particle acceleration");

            tempAcc[1] = (D3 * hydrodynForcesAtIteration[0][0] - D1 * hydrodynForcesAtIteration[0][1]) / (D3 * D4 - D1 * D2);
            if (double.IsNaN(tempAcc[1]) || double.IsInfinity(tempAcc[1]))
                throw new ArithmeticException("Error trying to calculate particle acceleration");

            aux.SaveMultidimValueToList(transAccelerationAtIteration, tempAcc);
            transAccelerationAtTimestep[0] = transAccelerationAtIteration[0];
        }

        public void PredictTranslationalVelocity()
        {
            double[] temp = new double[m_Dim];
            for (int d = 0; d < m_Dim; d++)
            {
                temp[d] = 2 * transVelocityAtTimestep[0][d] - transVelocityAtTimestep[1][d];
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to predict particle velocity");
            }
            aux.SaveMultidimValueToList(transVelocityAtIteration, temp);
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
                aux.SaveMultidimValueOfLastTimestep(transVelocityAtTimestep);
            }

            double[] temp = new double[m_Dim];
            for (int d = 0; d < m_Dim; d++)
            {
                temp[d] = transVelocityAtTimestep[1][d] + (transAccelerationAtTimestep[1][d] + transAccelerationAtIteration[0][d]) * dt / 2;
                if (double.IsNaN(temp[d]) || double.IsInfinity(temp[d]))
                    throw new ArithmeticException("Error trying to calculate particle velocity");
            }
            aux.SaveMultidimValueToList(transVelocityAtIteration, temp);
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
            double temp = 2 * rotationalAccelarationAtTimestep[0] - rotationalAccelarationAtTimestep[1];
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to predict particle angluar acceleration");
            aux.SaveValueToList(rotationalAccelarationAtIteration, temp);
            rotationalAccelarationAtTimestep[0] = rotationalAccelarationAtIteration[0];
        }

        public void CalculateAngularAcceleration(double dt, double addedDampingCoeff = 1)
        {
            double MomentofInertia_m = MomentOfInertia_P + addedDampingCoeff * dt * addedDampingTensorVV[0, 0];
            aux.SaveValueToList(rotationalAccelarationAtIteration, hydrodynTorqueAtIteration[0] / MomentofInertia_m);
        }

        public void PredictAngularVelocity()
        {
            double temp = 2 * rotationalVelocityAtTimestep[0] - rotationalVelocityAtTimestep[1];
            if (double.IsNaN(temp) || double.IsInfinity(temp))
                throw new ArithmeticException("Error trying to predict particle angluar velocity");
            aux.SaveValueToList(rotationalVelocityAtIteration, temp);
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
                aux.SaveValueOfLastTimestep(rotationalVelocityAtTimestep);
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
            aux.SaveValueToList(rotationalVelocityAtIteration, newAngularVelocity);
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
            if (neglectAddedDamping == true)
            {
                return;
            }
            int D = LsTrk.GridDat.SpatialDimension;
            double alpha = 0.5;
            int RequiredOrder = 2;
            for (int i = 0; i < 4; i++)
            {
                for (int d1 = 0; d1 < D; d1++)
                {
                    for (int d2 = 0; d2 < D; d2++)
                    {
                        ScalarFunctionEx evalfD = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result)
                        {
                            int K = result.GetLength(1);
                            // Normal vector
                            var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);


                            if (LsTrk.GridDat.SpatialDimension == 2)
                            {
                                for (int j = 0; j < Len; j++)
                                {
                                    double dh = Math.Sqrt(LsTrk.GridDat.iGeomCells.GetCellVolume(j)); //just an approx, needs to be revisited
                                    double delta = dh * Math.Sqrt(rhoA) / (Math.Sqrt(alpha * muA * dt));
                                    double dn = dh / (1 - Math.Exp(-delta));

                                    for (int k = 0; k < K; k++)
                                    {
                                        double[] R = new double[D];
                                        R[0] = 1;// Ns[0] - particlePositionPerTimestep[0][0];
                                        R[1] = 1;//Ns[1] - particlePositionPerTimestep[0][1];
                                        switch (i)
                                        {
                                            case 0:
                                                result[j, k] = d1 == d2 ? (1 - Normals[j, k, d1] * Normals[j, k, d2]) * muA / dn : (0 - Normals[j, k, d1] * Normals[j, k, d2]) * muA / dn;
                                                break;
                                            case 1:
                                                result[j, k] = 0;
                                                break;
                                            case 2:
                                                result[j, k] = 0;
                                                break;
                                            case 3:
                                                result[j, k] = d1 == d2 ? -(R[1 - d1] * R[1 - d2]) * muA / dn : R[1 - d1] * R[1 - d2] * muA / dn;
                                                break;
                                        }
                                    }
                                }
                            }
                        };

                        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                        //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );

                        //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
                        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, this.cutCells_P(LsTrk));

                        CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                            cqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult)
                            {
                                evalfD(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                            },
                            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                            {
                                for (int l = 0; l < Length; l++)
                                {
                                    switch (i)
                                    {
                                        case 0:
                                            addedDampingTensorVV[d1, d2] += ResultsOfIntegration[l, 0];
                                            break;
                                        case 1:
                                            addedDampingTensorVW[d1, d2] += ResultsOfIntegration[l, 0];
                                            break;
                                        case 2:
                                            addedDampingTensorWV[d1, d2] += ResultsOfIntegration[l, 0];
                                            break;
                                        case 3:
                                            addedDampingTensorWW[d1, d2] += ResultsOfIntegration[l, 0];
                                            break;
                                    }
                                }
                            }
                        ).Execute();
                    }
                }
            }
        }

        public void UpdateDampingTensors()
        {
            // form rotation matrix R=EpEp^T, where Ep is the matrix of the principle axis of inertia
            // symmetry axis are always axis of inertia:
            double[,] Ep = new double[2, 2];
            Ep[0, 0] = Math.Cos(angleAtIteration[0]);
            Ep[1, 0] = Math.Sin(angleAtIteration[0]);
            Ep[0, 1] = -Math.Sin(angleAtIteration[0]);
            Ep[1, 1] = Math.Cos(angleAtIteration[0]);

            double[,] R = new double[2, 2];
            R[0, 0] = Ep[0, 0].Pow2() + Ep[0, 1] * Ep[1, 0];
            R[1, 0] = Ep[0, 0] * Ep[0, 1] + Ep[0, 1] * Ep[1, 1];
            R[0, 1] = Ep[1, 0] * Ep[0, 0] + Ep[1, 1] * Ep[1, 0];
            R[1, 1] = Ep[1, 0] * Ep[0, 1] + Ep[1, 1].Pow2();

            addedDampingTensorVV[0, 0] = R[0, 0].Pow2() * addedDampingTensorVV[0, 0] + R[1, 0] * R[0, 0] * addedDampingTensorVV[0, 1] + R[0, 0] * R[1, 0] * addedDampingTensorVV[1, 0] + R[1, 0].Pow2() * addedDampingTensorVV[1, 1];
            addedDampingTensorVV[1, 0] = R[0, 0] * R[0, 1] * addedDampingTensorVV[0, 0] + R[1, 0] * R[0, 1] * addedDampingTensorVV[0, 1] + R[0, 0] * R[1, 1] * addedDampingTensorVV[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorVV[1, 1];
            addedDampingTensorVV[0, 1] = R[0, 0] * R[0, 1] * addedDampingTensorVV[0, 0] + R[0, 0] * R[1, 1] * addedDampingTensorVV[0, 1] + R[0, 0] * R[0, 1] * addedDampingTensorVV[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorVV[1, 1];
            addedDampingTensorVV[1, 1] = R[0, 1].Pow2() * addedDampingTensorVV[0, 0] + R[1, 1] * R[0, 1] * addedDampingTensorVV[0, 1] + R[0, 1] * R[1, 1] * addedDampingTensorVV[1, 0] + R[1, 1].Pow2() * addedDampingTensorVV[1, 1];

            addedDampingTensorVW[0, 0] = R[0, 0].Pow2() * addedDampingTensorVW[0, 0] + R[1, 0] * R[0, 0] * addedDampingTensorVW[0, 1] + R[0, 0] * R[1, 0] * addedDampingTensorVW[1, 0] + R[1, 0].Pow2() * addedDampingTensorVW[1, 1];
            addedDampingTensorVW[1, 0] = R[0, 0] * R[0, 1] * addedDampingTensorVW[0, 0] + R[1, 0] * R[0, 1] * addedDampingTensorVW[0, 1] + R[0, 0] * R[1, 1] * addedDampingTensorVW[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorVW[1, 1];
            addedDampingTensorVW[0, 1] = R[0, 0] * R[0, 1] * addedDampingTensorVW[0, 0] + R[0, 0] * R[1, 1] * addedDampingTensorVW[0, 1] + R[0, 0] * R[0, 1] * addedDampingTensorVW[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorVW[1, 1];
            addedDampingTensorVW[1, 1] = R[0, 1].Pow2() * addedDampingTensorVW[0, 0] + R[1, 1] * R[0, 1] * addedDampingTensorVW[0, 1] + R[0, 1] * R[1, 1] * addedDampingTensorVW[1, 0] + R[1, 1].Pow2() * addedDampingTensorVW[1, 1];

            addedDampingTensorWV[0, 0] = R[0, 0].Pow2() * addedDampingTensorWV[0, 0] + R[1, 0] * R[0, 0] * addedDampingTensorWV[0, 1] + R[0, 0] * R[1, 0] * addedDampingTensorWV[1, 0] + R[1, 0].Pow2() * addedDampingTensorWV[1, 1];
            addedDampingTensorWV[1, 0] = R[0, 0] * R[0, 1] * addedDampingTensorWV[0, 0] + R[1, 0] * R[0, 1] * addedDampingTensorWV[0, 1] + R[0, 0] * R[1, 1] * addedDampingTensorWV[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorWV[1, 1];
            addedDampingTensorWV[0, 1] = R[0, 0] * R[0, 1] * addedDampingTensorWV[0, 0] + R[0, 0] * R[1, 1] * addedDampingTensorWV[0, 1] + R[0, 0] * R[0, 1] * addedDampingTensorWV[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorWV[1, 1];
            addedDampingTensorWV[1, 1] = R[0, 1].Pow2() * addedDampingTensorWV[0, 0] + R[1, 1] * R[0, 1] * addedDampingTensorWV[0, 1] + R[0, 1] * R[1, 1] * addedDampingTensorWV[1, 0] + R[1, 1].Pow2() * addedDampingTensorWV[1, 1];

            addedDampingTensorWW[0, 0] = R[0, 0].Pow2() * addedDampingTensorWW[0, 0] + R[1, 0] * R[0, 0] * addedDampingTensorWW[0, 1] + R[0, 0] * R[1, 0] * addedDampingTensorWW[1, 0] + R[1, 0].Pow2() * addedDampingTensorWW[1, 1];
            addedDampingTensorWW[1, 0] = R[0, 0] * R[0, 1] * addedDampingTensorWW[0, 0] + R[1, 0] * R[0, 1] * addedDampingTensorWW[0, 1] + R[0, 0] * R[1, 1] * addedDampingTensorWW[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorWW[1, 1];
            addedDampingTensorWW[0, 1] = R[0, 0] * R[0, 1] * addedDampingTensorWW[0, 0] + R[0, 0] * R[1, 1] * addedDampingTensorWW[0, 1] + R[0, 0] * R[0, 1] * addedDampingTensorWW[1, 0] + R[1, 0] * R[1, 1] * addedDampingTensorWW[1, 1];
            addedDampingTensorWW[1, 1] = R[0, 1].Pow2() * addedDampingTensorWW[0, 0] + R[1, 1] * R[0, 1] * addedDampingTensorWW[0, 1] + R[0, 1] * R[1, 1] * addedDampingTensorWW[1, 0] + R[1, 1].Pow2() * addedDampingTensorWW[1, 1];

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
                aux.SaveMultidimValueOfLastTimestep(hydrodynForcesAtTimestep);
                aux.SaveValueOfLastTimestep(hydrodynTorqueAtTimestep);
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
                ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of Nodes
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, spatialDim, spatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                    // Evaluate tangential velocity to level-set surface
                    // =============================
                    // Normal vector
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                    // Velocity
                    for (int i = 0; i < spatialDim; i++) {
                        UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    // Pressure
                    pA.Evaluate(j0, Len, Ns, pARes);

                    if (LsTrk.GridDat.SpatialDimension == 2) {
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                // Defining variables
                                double acc = 0.0;
                                double[] SummandsVelGradient = new double[3];
                                double SummandsPressure;
                                switch (d) {
                                    case 0:
                                        SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                        SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        SummandsPressure = pARes[j, k] * Normals[j, k, 0];
                                        acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                                        break;

                                    case 1:
                                        SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                        SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                        SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                        SummandsPressure = pARes[j, k] * Normals[j, k, 1];
                                        acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                                        break;
                                    default:
                                        throw new NotImplementedException();
                                }
                                result[j, k] = acc;
                            }
                        }
                    }
                    else {
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;
                                double[] SummandsVelGradient = new double[5];
                                double SummandsPressure;
                                switch (d) {
                                    case 0:
                                        SummandsPressure = pARes[j, k] * Normals[j, k, 0];
                                        SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                        SummandsVelGradient[1] = -Grad_UARes[j, k, 0, 2] * Normals[j, k, 2];
                                        SummandsVelGradient[2] = -Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        SummandsVelGradient[3] = -Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 0] * Normals[j, k, 2];
                                        acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                                        break;
                                    case 1:
                                        SummandsPressure = pARes[j, k] * Normals[j, k, 1];
                                        SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                        SummandsVelGradient[1] = -Grad_UARes[j, k, 1, 2] * Normals[j, k, 2];
                                        SummandsVelGradient[2] = -Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                        SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                        SummandsVelGradient[4] = -Grad_UARes[j, k, 2, 1] * Normals[j, k, 2];
                                        acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                                        break;
                                    case 2:
                                        SummandsPressure = pARes[j, k] * Normals[j, k, 2];
                                        SummandsVelGradient[0] = -2 * Grad_UARes[j, k, 2, 2] * Normals[j, k, 2];
                                        SummandsVelGradient[1] = -Grad_UARes[j, k, 2, 0] * Normals[j, k, 0];
                                        SummandsVelGradient[2] = -Grad_UARes[j, k, 2, 1] * Normals[j, k, 1];
                                        SummandsVelGradient[3] = -Grad_UARes[j, k, 0, 2] * Normals[j, k, 0];
                                        SummandsVelGradient[4] = -Grad_UARes[j, k, 1, 2] * Normals[j, k, 1];
                                        acc += aux.SummationWithNeumaier(SummandsVelGradient, SummandsPressure, muA);
                                        break;
                                    default:
                                        throw new NotImplementedException();
                                }

                                result[j, k] = acc;
                            }
                        }
                    }

                };

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );

                //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, this.cutCells_P(LsTrk));

                double forceNaiveSum = 0.0;
                double forceSum = 0.0;
                double forceC = 0.0;
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                        {
                            //Forces[d] += ResultsOfIntegration[i, 0];
                            forceNaiveSum = forceSum + ResultsOfIntegration[i, 0];
                            if (Math.Abs(forceSum) >= Math.Abs(ResultsOfIntegration[i, 0]))
                            {
                                forceC += (forceSum - forceNaiveSum) + ResultsOfIntegration[i, 0];
                            }
                            else
                            {
                                forceC += (ResultsOfIntegration[i, 0] - forceNaiveSum) + forceSum;
                            }
                            forceSum = forceNaiveSum;
                        }
                        Forces[d] = forceSum + forceC;
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

                //var trafo = LsTrk.GridDat.Edges.Edge2CellTrafos;
                //var trafoIdx = LsTrk.GridDat.TransformLocal2Global(Ns)
                //var transFormed = trafo[trafoIdx].Transform(Nodes);
                //var newVertices = transFormed.CloneAs();
                //GridData.TransformLocal2Global(transFormed, newVertices, jCell);


                MultidimensionalArray tempArray = Ns.CloneAs();

                LsTrk.GridDat.TransformLocal2Global(Ns, tempArray, j0);

                pA.Evaluate(j0, Len, Ns, pARes);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        double[] integrand = new double[4];
                        double[] integrand2 = new double[4];
                        double naiveSum = 0.0;
                        double c = 0.0;
                        double sum = 0.0;
                        double sum2 = 0.0;
                        double naiveSum2 = 0.0;
                        double c2 = 0.0;

                        // Calculate the Torque around a circular particle with a given radius (Paper Wan and Turek 2005)
                        integrand[0] = -2 * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                        integrand[1] = -Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                        integrand[2] = -Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                        integrand[3] = (pARes[j, k] * Normals[j, k, 0]);
                        sum = integrand[0];
                        for (int i = 1; i < integrand.Length - 1; i++)
                        {
                            naiveSum = sum + integrand[i];
                            if (Math.Abs(sum) >= integrand[i])
                            {
                                c += (sum - naiveSum) + integrand[i];
                            }
                            else
                            {
                                c += (integrand[i] - naiveSum) + sum;
                            }
                            sum = naiveSum;
                        }
                        sum *= muA;
                        c *= muA;
                        naiveSum = sum + integrand[3];
                        if (Math.Abs(sum) >= integrand[3])
                        {
                            c += (sum - naiveSum) + integrand[3];
                        }
                        else
                        {
                            c += (integrand[3] - naiveSum) + sum;
                        }
                        sum *= -Normals[j, k, 1] * (this.positionAtIteration[0][1] - tempArray[k, 1]).Abs();
                        c *= -Normals[j, k, 1] * (this.positionAtIteration[0][1] - tempArray[k, 1]).Abs();

                        integrand2[0] = -2 * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                        integrand2[1] = -Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                        integrand2[2] = -Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                        integrand2[3] = pARes[j, k] * Normals[j, k, 1];
                        sum2 = integrand2[0];
                        for (int i = 1; i < integrand2.Length - 1; i++)
                        {
                            naiveSum2 = sum2 + integrand2[i];
                            if (Math.Abs(sum2) >= integrand2[i])
                            {
                                c2 += (sum2 - naiveSum2) + integrand2[i];
                            }
                            else
                            {
                                c2 += (integrand2[i] - naiveSum2) + sum2;
                            }
                            sum2 = naiveSum2;
                        }
                        sum2 *= muA;
                        c2 *= muA;
                        naiveSum2 = sum2 + integrand2[3];
                        if (Math.Abs(sum2) >= integrand2[3])
                        {
                            c2 += (sum2 - naiveSum2) + integrand2[3];
                        }
                        else
                        {
                            c2 += (integrand2[3] - naiveSum2) + sum2;
                        }
                        sum2 *= Normals[j, k, 0] * (this.positionAtIteration[0][0] - tempArray[k, 0]).Abs();
                        c2 *= Normals[j, k, 0] * (this.positionAtIteration[0][0] - tempArray[k, 0]).Abs();
                        sum += sum2;
                        c += c2;

                        result[j, k] = sum + c;
                    }
                }
            };

            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );
            //CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, this.cutCells_P(LsTrk));
            double torqueNaiveSum = 0.0;
            double torqueSum = 0.0;
            double torqueC = 0.0;
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs2.Compile(LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    //for (int i = 0; i < Length; i++)
                    //{
                    //    Torque += ResultsOfIntegration[i, 0];
                    //}
                    for (int i = 0; i < Length; i++)
                    {
                        torqueNaiveSum = torqueSum + ResultsOfIntegration[i, 0];
                        if (Math.Abs(torqueSum) >= Math.Abs(ResultsOfIntegration[i, 0]))
                        {
                            torqueC += (torqueSum - torqueNaiveSum) + ResultsOfIntegration[i, 0];
                        }
                        else
                        {
                            torqueC += (ResultsOfIntegration[i, 0] - torqueNaiveSum) + torqueSum;
                        }
                        torqueSum = torqueNaiveSum;
                    }
                    Torque = torqueSum + torqueC;
                }

            ).Execute();

            // determine underrelaxation factor (URF)
            // =============================
            double[] ForcesUnderrelaxation = new double[spatialDim];
            double TorqueUnderrelaxation;

            double averageForce = aux.CalculateAverageForces(Forces, Torque, averageDistance);
            if (iteration_counter_P == 0)
            {
                for (int k = 0; k < spatialDim; k++)
                {
                    ForcesUnderrelaxation[k] = 1;
                    for (int t = 0; t < m_HistoryLength; t++)
                    {
                        hydrodynForcesAtIteration[t][k] = hydrodynForcesAtTimestep[1][k];
                        hydrodynTorqueAtIteration[t] = hydrodynTorqueAtTimestep[1];
                    }
                }
                TorqueUnderrelaxation = 1;
            }
            //// restart iteration
            //// =============================
            //else if ((iteration_counter_P - 1) / 100 % 2 == 1 && Math.Sqrt((Forces[0] - hydrodynForcesAtIteration[1][0]).Pow2()+ (Forces[1] - hydrodynForcesAtIteration[1][1]).Pow2()+ (Torque-hydrodynTorqueAtIteration[1]).Pow2()) > 100 * forceAndTorque_convergence)
            //{
            //    Forces[0] = 0.00125 * Circumference_P * active_stress_P.Pow2() * Math.Cos(angleAtIteration[0]) / (muA);
            //    Forces[1] = 0.00125 * Circumference_P * active_stress_P.Pow2() * Math.Sin(angleAtIteration[0]) / muA;
            //    Torque = 0;
            //}
            // constant predefined URF
            // =============================
            else if (underrelaxationFT_constant == true)
            {
                for (int k = 0; k < spatialDim; k++)
                {
                    ForcesUnderrelaxation[k] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                }
                TorqueUnderrelaxation = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
            }
            // calculation of URF for adaptive underrelaxation
            // =============================
            else 
            {
                ForcesUnderrelaxation = aux.CalculateAdaptiveForceUnderrelaxation(Forces, hydrodynForcesAtIteration[0], averageForce, forceAndTorque_convergence, underrelaxation_factor);
                TorqueUnderrelaxation = aux.CalculateAdaptiveTorqueUnderrelaxation(Torque, hydrodynTorqueAtIteration[0], averageForce, forceAndTorque_convergence, underrelaxation_factor);
            }
            Console.WriteLine("ForcesUnderrelaxation[0]  " + ForcesUnderrelaxation[0] + ", ForcesUnderrelaxation[1]: " + ForcesUnderrelaxation[1] + ", TorqueUnderrelaxation " + TorqueUnderrelaxation);
            Console.WriteLine("tempfForces[0]  " + Forces[0] + ", temp_Forces[1]: " + Forces[1] + ", tempTorque " + Torque);

            // calculation of Forces and Torque with underrelaxation
            // =============================
            int beta = 1;
            Forces[0] = Forces[0] - addedDampingTensorVV[0, 0] * beta * transAccelerationAtIteration[0][0] * dt - addedDampingTensorVV[1, 0] * beta * transAccelerationAtIteration[0][1] * dt;
            Forces[1] = Forces[1] - addedDampingTensorVV[0, 1] * beta * transAccelerationAtIteration[0][0] * dt - addedDampingTensorVV[1, 1] * beta * transAccelerationAtIteration[0][1] * dt + (particleDensity - fluidDensity) * Area_P * gravityVertical;
            Torque = Torque - beta * dt * addedDampingTensorWW[0, 0] * rotationalAccelarationAtIteration[0];
            Forces = aux.RelaxatedForce(ForcesUnderrelaxation, Forces, hydrodynForcesAtIteration[0], ClearSmallValues, forceAndTorque_convergence);
            Torque = aux.RelaxatedTorque(TorqueUnderrelaxation, Torque, hydrodynTorqueAtIteration[0], ClearSmallValues, forceAndTorque_convergence);
            aux.SaveMultidimValueToList(hydrodynForcesAtIteration, Forces);
            aux.SaveValueToList(hydrodynTorqueAtIteration, Torque);
            hydrodynForcesAtTimestep[0] = hydrodynForcesAtIteration[0];
            hydrodynTorqueAtTimestep[0] = hydrodynTorqueAtIteration[0];

            #endregion
        }
        #endregion

        #region Particle reynolds number
        /// <summary>
        /// Calculating the particle reynolds number according to paper Turek and testcase ParticleUnderGravity
        /// </summary>
        /// <param name="currentVelocity"></param>
        /// <param name="Radius"></param>
        /// <param name="particleDensity"></param>
        /// <param name="viscosity"></param>
        /// <returns></returns>
        abstract public double ComputeParticleRe(double mu_Fluid);
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

