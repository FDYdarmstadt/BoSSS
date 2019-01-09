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

namespace BoSSS.Application.FSI_Solver
{

    /// <summary>
    /// Particle properties (for disk shape and spherical particles only).
    /// </summary>
    [DataContract]
    [Serializable]
    public class Particle : ICloneable {

        #region particle
        public Particle(int Dim, int HistoryLength, double[] startPos = null, double startAngl = 0.0, ParticleShape shape = ParticleShape.spherical) {

            
            m_Dim = Dim;
            m_HistoryLength = HistoryLength;
            m_shape = shape;


            #region Particle history
            // =============================   
            for (int i = 0; i < HistoryLength; i++) {
                currentIterPos_P.Add(new double[Dim]);
                currentIterAng_P.Add(new double());
                currentIterVel_P.Add(new double[Dim]);
                currentIterRot_P.Add(new double());
                currentIterForces_P.Add(new double[Dim]);
                temporalForces_P.Add(new double[Dim]);
                currentIterTorque_P.Add(new double());
                temporalTorque_P.Add(new double());
            }
            for (int i = 0; i < 4; i++)
            {
                currentTimePos_P.Add(new double[Dim]);
                currentTimeAng_P.Add(new double());
                currentTimeVel_P.Add(new double[Dim]);
                currentTimeRot_P.Add(new double());
                currentTimeForces_P.Add(new double[Dim]);
                currentTimeTorque_P.Add(new double());
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
            currentTimePos_P[0] = startPos;
            currentTimePos_P[1] = startPos;
            //From degree to radiant
            currentTimeAng_P[0] = startAngl * 2 * Math.PI / 360;
            currentTimeAng_P[1] = startAngl * 2 * Math.PI / 360;
            //currentIterVel_P[0][0] = 2e-8;

            UpdateLevelSetFunction();
            #endregion
        }
        #endregion

        #region Particle parameter
        /// <summary>
        /// empty constructor important for serialization
        /// </summary>
        private Particle() {
        }

        /// <summary>
        /// Dimension of the particles (Disk or Sphere)
        /// </summary>
        int m_Dim;

        public bool[] m_collidedWithParticle;
        public bool[] m_collidedWithWall;
        public double[][] m_closeInterfacePointTo;


        public double[] tempPos_P = new double[2];
        public double tempAng_P;
        public int iteration_counter_P = 0;
        public bool underrelaxationFT_constant = true;
        public int underrelaxationFT_exponent = 0;
        public double underrelaxation_factor = 1;
        public int underrelaxationFT_exponent_min = 0;
        public int underrelaxationFT_exponent_max = 0;

        /// <summary>
        /// Shape of the particle
        /// </summary>
        public ParticleShape m_shape;

        /// <summary>
        /// Skip calculation of hydrodynamic force and torque if particles are too close
        /// </summary>
        public bool skipForceIntegration = false;

        /// <summary>
        /// Length of history for time, velocity, position etc.
        /// </summary>
        int m_HistoryLength;

        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double rho_P;

        /// <summary>
        /// Radius of the particle.
        /// </summary>
        [DataMember]
        public double radius_P;

        /// <summary>
        /// Lenght of an elliptic particle.
        /// </summary>
        [DataMember]
        public double length_P;

        /// <summary>
        /// Thickness of an elliptic particle.
        /// </summary>
        [DataMember]
        public double thickness_P;

        /// <summary>
        /// Current position of the particle (the center of mass)
        /// </summary>
        [DataMember]
        public List<double[]> currentIterPos_P = new List<double[]>();

        /// <summary>
        /// Current position of the particle (the center of mass)
        /// </summary>
        [DataMember]
        public List<double[]> currentTimePos_P = new List<double[]>();

        /// <summary>
        /// Current angle of the particle (the center of mass)
        /// </summary>
        [DataMember]
        public List<double> currentIterAng_P = new List<double>();

        /// <summary>
        /// Current angle of the particle (the center of mass)
        /// </summary>
        [DataMember]
        public List<double> currentTimeAng_P = new List<double>();

        /// <summary>
        /// Current translational Velocity of the particle
        /// </summary>
        [DataMember]
        public List<double[]> currentIterVel_P = new List<double[]>();

        /// <summary>
        /// Current translational Velocity of the particle
        /// </summary>
        [DataMember]
        public List<double[]> currentTimeVel_P = new List<double[]>();

        /// <summary>
        /// Current rotational Velocity of the particle
        /// </summary>
        [DataMember]
        public List<double> currentIterRot_P = new List<double>();

        /// <summary>
        /// Current rotational Velocity of the particle
        /// </summary>
        [DataMember]
        public List<double> currentTimeRot_P = new List<double>();

        /// <summary>
        /// Force acting on the particle
        /// </summary>
        [DataMember]
        public List<double[]> currentIterForces_P = new List<double[]>();

        /// <summary>
        /// Force acting on the particle
        /// </summary>
        [DataMember]
        public List<double[]> currentTimeForces_P = new List<double[]>();

        /// <summary>
        /// Force acting on the particle
        /// </summary>
        [DataMember]
        public List<double[]> temporalForces_P = new List<double[]>();

        /// <summary>
        /// Torque acting on the particle
        /// </summary>
        [DataMember]
        public List<double> currentIterTorque_P = new List<double>();

        /// <summary>
        /// Torque acting on the particle
        /// </summary>
        [DataMember]
        public List<double> currentTimeTorque_P = new List<double>();

        /// <summary>
        /// Torque acting on the particle
        /// </summary>
        [DataMember]
        public List<double> temporalTorque_P = new List<double>();

        /// <summary>
        /// Level set function describing the particle
        /// </summary>       
        public Func<double[], double, double> phi_P;

        /// <summary>
        /// Set true if the particle should be an active particle, i.e. self driven
        /// </summary>
        [DataMember]
        public bool active_P = false;

        /// <summary>
        /// Set true if the particle should be an active particle, i.e. self driven
        /// </summary>
        [DataMember]
        public bool includeGravity = true;

        /// <summary>
        /// Active stress on the current particle
        /// </summary>
        public double stress_magnitude_P;

        /// <summary>
        /// heaviside function depending on arclength 
        /// </summary>
        public int H_P;

        /// <summary>
        /// heaviside function depending on arclength 
        /// </summary>
        public double C_v = 0.5;

        /// <summary>
        /// heaviside function depending on arclength 
        /// </summary>
        public double velResidual_ConvergenceCriterion = 1e-6;

        /// <summary>
        /// heaviside function depending on arclength 
        /// </summary>
        public double forceAndTorque_convergence = 1e-6;

        /// <summary>
        /// heaviside function depending on arclength 
        /// </summary>
        public double MaxParticleVelIterations = 100;

        /// <summary>
        /// Active stress on the current particle
        /// </summary>
        virtual public double active_stress_P 
        {
            get
            {
                double stress;
                switch (m_shape)
                {
                    case ParticleShape.spherical:
                        stress = 2 * Math.PI * radius_P * stress_magnitude_P;
                        break;

                    case ParticleShape.elliptic:
                        //Approximation formula for circumference according to Ramanujan
                        double circumference;
                        circumference = Math.PI * ((length_P + thickness_P) + (3 * (length_P - thickness_P).Pow2()) / (10 * (length_P + thickness_P) + Math.Sqrt(length_P.Pow2() + 14 * length_P * thickness_P + thickness_P.Pow2()))); 
                        stress = 0.5 * circumference * stress_magnitude_P;
                        break;

                    case ParticleShape.squircle:
                        stress = 2 * Math.PI * 3.708 * stress_magnitude_P;
                        break;

                    default:

                        throw new NotImplementedException("");
                }

                return stress;
            }

        }
        
        /// <summary>
        /// Mass of the current particle
        /// </summary>
        [DataMember]
        public double mass_P {
            get
            {
                return area_P * rho_P;
            }

        }

        /// <summary>
        /// Area of the current particle
        /// </summary>
        [DataMember]
        virtual public double area_P {
            get {
                double area;
                switch (m_shape) {
                    case ParticleShape.spherical:
                        area = Math.PI * radius_P * radius_P;
                        break;

                    case ParticleShape.elliptic:
                        area = length_P * thickness_P * Math.PI * radius_P;
                        break;

                    case ParticleShape.hippopede:
                        // not correct mass
                        area = Math.PI * radius_P * radius_P;
                        break;

                    case ParticleShape.bean:
                        // not correct mass
                        area = Math.PI * radius_P * radius_P;
                        break;

                    case ParticleShape.squircle:
                        area = 3.708 * radius_P.Pow2();
                        break;

                    default:

                        throw new NotImplementedException("");
                }

                return area;
            }

        }

        /// <summary>
        /// Moment of inertia of the current particle
        /// </summary>
        [DataMember]
        virtual public double MomentOfInertia_P {
            get {
                double moment;
                switch (m_shape) {
                    case ParticleShape.spherical:
                        moment = (1 / 2.0) * (mass_P * radius_P * radius_P);
                        break;

                    case ParticleShape.elliptic:
                        moment = (1 / 4.0) * (mass_P * (length_P * length_P + thickness_P * thickness_P) * radius_P * radius_P);
                        break;

                    case ParticleShape.hippopede:
                        // not correct moment of inertia
                        moment = (1 / 2.0) * (mass_P * radius_P * radius_P);
                        break;

                    case ParticleShape.bean:
                        // not correct moment of inertia
                        moment = (1 / 2.0) * (mass_P * radius_P * radius_P);
                        break;

                    case ParticleShape.squircle:
                        // not correct moment of inertia
                        moment = (1 / 2.0) * (mass_P * radius_P * radius_P);
                        break;


                    default:

                        throw new NotImplementedException("");
                }

                return moment;
            }
        }
        #endregion

        #region Particle history
        /// <summary>
        /// Clean all Particle histories until a certain length
        /// </summary>
        /// <param name="length"></param>
        // iteration
        public void CleanHistoryIter() {
            if (currentIterPos_P.Count > m_HistoryLength)
            {
                for (int j = currentIterPos_P.Count; j > m_HistoryLength; j--)
                {
                    int tempPos = j - 1;
                    currentIterPos_P.RemoveAt(tempPos);
                    currentIterVel_P.RemoveAt(tempPos);
                    currentIterForces_P.RemoveAt(tempPos);
                    currentIterRot_P.RemoveAt(tempPos);
                    currentIterTorque_P.RemoveAt(tempPos);
                    temporalForces_P.RemoveAt(tempPos);
                    temporalTorque_P.RemoveAt(tempPos);
                }
            }
        }
        // time
        public void CleanHistory()
        {
            if (currentTimePos_P.Count > 4)
            {
                for (int j = currentTimePos_P.Count; j > 4; j--)
                {
                    int tempPos = j - 1;
                    currentTimePos_P.RemoveAt(tempPos);
                    currentTimeAng_P.RemoveAt(tempPos);
                    currentTimeVel_P.RemoveAt(tempPos);
                    currentTimeForces_P.RemoveAt(tempPos);
                    currentTimeRot_P.RemoveAt(tempPos);
                    currentTimeTorque_P.RemoveAt(tempPos);
                }
            }
        }
        #endregion

        #region Move particle with current velocity
        /// <summary>
        /// Move particle with current velocity
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticlePosition(double dt) {
            
            // Position
            var tempPos = currentIterPos_P[0].CloneAs();
            tempPos.AccV(dt, currentIterVel_P[0]);

            //currentIterPos_P.Insert(0, tempPos);
            //currentIterPos_P.Remove(currentIterPos_P.Last());
            currentIterPos_P[0] = tempPos;

            // Angle
            var tempAng = currentIterAng_P[0];
            tempAng += dt * currentIterRot_P[0];
            //currentIterAng_P.Insert(0, tempAng);
            //currentIterAng_P.Remove(currentIterAng_P.Last());
            currentIterAng_P[0] = tempAng;

            currentTimePos_P[0] = currentIterPos_P[0];
            currentTimeAng_P[0] = currentIterAng_P[0];

            //Console.WriteLine("Current angle speed is " + currentIterRot_P[0]*360/Math.PI);
            //Console.WriteLine("Current angle is " + currentIterAng_P[0] + " rad");

            UpdateLevelSetFunction();
        }
        public void ResetParticlePosition()
        {
            // save position of the last timestep
            if (iteration_counter_P == 1)
            {
                currentTimePos_P[1] = currentTimePos_P[0];
                currentTimeAng_P[1] = currentTimeAng_P[0];
            }

            // Position
            currentIterPos_P.Insert(1, currentIterPos_P[0]);
            currentIterPos_P.Remove(currentIterPos_P.Last());
            currentIterPos_P[0] = currentTimePos_P[1];

            // Angle
            currentIterAng_P.Insert(1, currentTimeAng_P[0]);
            currentIterAng_P.Remove(currentIterAng_P.Last());
            currentIterAng_P[0] = currentTimeAng_P[1];

            UpdateLevelSetFunction();
        }
        #endregion

        #region Update Level-set
        public void UpdateLevelSetFunction() {

            double a;
            double b;
            double alpha;

            switch (m_shape) {

                case ParticleShape.spherical:
                    phi_P = (X, t) => -(X[0] - currentIterPos_P[0][0]).Pow2() + -(X[1] - currentIterPos_P[0][1]).Pow2() + radius_P.Pow2();
                    break;

                case ParticleShape.elliptic:
                    alpha = -(currentIterAng_P[0]);
                    a = 3.0;//no longer necessary
                    b = 1.0;//no longer necessary
                    phi_P = (X, t) => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radius_P.Pow2();
                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radius_P.Pow2();
                    b = 1.0 * radius_P.Pow2();
                    alpha = -(currentIterAng_P[0]);
                    phi_P = (X, t) => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow2() - b * ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2());
                    break;


                case ParticleShape.bean:
                    a = 3.0 * radius_P.Pow2();
                    b = 1.0 * radius_P.Pow2();
                    alpha = -(currentIterAng_P[0]);
                    phi_P = (X, t) => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2());
                    break;

                case ParticleShape.squircle:
                    alpha = -(currentIterAng_P[0]);
                    phi_P = (X, t) => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(4) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(4)) - radius_P.Pow(4));
                    break;

                default:
                    throw new NotImplementedException("Shape is not implemented yet");
            }
        }
        #endregion

        #region Update translational velocity
        /// <summary>
        /// Calculate the new translational velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="D"></param>
        /// <param name="dt"></param>
        /// <param name="oldTransVelocity"></param>
        /// <param name="particleMass"></param>
        /// <param name="force"></param>
        /// <returns></returns>
        public void UpdateTransVelocity(double dt, double rho_Fluid = 1, bool includeTranslation = true, bool LiftAndDrag = true) {

            // save velocity of the last timestep
            if (iteration_counter_P == 1)
            {
                currentTimeVel_P[1] = currentIterVel_P[0];
            }

            // no translation
            // =============================
            if (includeTranslation == false)
            {
                currentIterVel_P[0][0] = 0.0;
                currentIterVel_P[0][1] = 0.0;
                return;
            }

            double[] temp = new double[2];
            double[] old_temp = new double[2];
            double[] tempForceTemp = new double[2];
            double[] tempForceNew = new double[2];
            double[] tempForceOld = new double[2];
            double[] gravity = new double[2];
            double massDifference = (rho_P - rho_Fluid) * (area_P);
            double[] previous_vel = new double[2];

            if (iteration_counter_P == 1)
            {
                previous_vel = currentIterVel_P[0];
            }
            
            // gravity
            // =============================
            gravity[0] = 0;
            if (includeGravity == true)
            {
                gravity[1] = -9.81e0;
            }
            else
            {
                gravity[1] = 0;
            }

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
            C_v_mod[0] = C_v;
            C_v_mod[1] = C_v;
            double[] c_a = new double[2];
            double[] c_u = new double[2];
            double vel_iteration_counter = 0;
            double[] test = new double[2];
            // 2nd order Adam Bashford
            //for (int i = 0; i < 2; i++)
            //{
            //    f_vNew[i] = c_a * (3 * currentIterVel_P[0][i] - 4 * currentIterVel_P[1][i] + currentIterVel_P[2][i]) / (2 * dt);
            //    f_vOld[i] = c_a * (3 * currentIterVel_P[1][i] - 4 * currentIterVel_P[2][i] + currentIterVel_P[3][i]) / (2 * dt);
            //    tempForceNew[i] = (forces_P[0][i] + massDifference * gravity[i]) * (c_u) + f_vNew[i];
            //    tempForceOld[i] = (forces_P[1][i] + massDifference * gravity[i]) * (c_u) + f_vOld[i];
            //    temp[i] = currentIterVel_P[0][i] + (3 * tempForceNew[i] - tempForceOld[i]) * dt / 2;
            //}

            // implicit Adams Moulton (modified)
            //for (double velResidual = 1; velResidual > velResidual_ConvergenceCriterion;)
            //{
            //    for (int i = 0; i < 2; i++)
            //    {
            //        C_v_mod[i] = 0.1;// * Math.Abs(forces_P[0][i] / (forces_P[0][i] + forces_P[1][i] + 1e-30));
            //        c_a[i] = (C_v_mod[i] * rho_Fluid) / (rho_P + C_v_mod[i] * rho_Fluid);
            //        c_u[i] = 1 / (area_P * (rho_P + C_v_mod[i] * rho_P));
            //        f_vTemp[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * temp[i] - 18 * currentIterVel_P[0][i] + 9 * currentIterVel_P[1][i] - 2 * currentIterVel_P[2][i]) / (8 * dt);
            //        f_vNew[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * currentIterVel_P[0][i] - 18 * currentIterVel_P[1][i] + 9 * currentIterVel_P[2][i] - 2 * currentIterVel_P[3][i]) / (6 * dt);
            //        f_vOld[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * currentIterVel_P[1][i] - 18 * currentIterVel_P[2][i] + 9 * currentIterVel_P[3][i] - 2 * currentIterVel_P[4][i]) / (6 * dt);
            //        tempForceTemp[i] = (forces_P[0][i] + massDifference * gravity[i]) * (c_u[i]) + f_vTemp[i];
            //        tempForceNew[i] = (forces_P[1][i] + massDifference * gravity[i]) * (c_u[i]) + f_vNew[i];
            //        tempForceOld[i] = (forces_P[2][i] + massDifference * gravity[i]) * (c_u[i]) + f_vOld[i];
            //        old_temp[i] = temp[i];
            //        temp[i] = previous_vel[i] + (1 * tempForceTemp[i] + 4 * tempForceNew[i] + 1 * tempForceOld[i]) * dt / 6;
            //    }
            //    vel_iteration_counter += 1;
            //    if (vel_iteration_counter == MaxParticleVelIterations)
            //    {
            //        throw new ApplicationException("no convergence in particle velocity calculation");
            //    }
            //    velResidual = Math.Sqrt((temp[0] - old_temp[0]).Pow2() + (temp[1] - old_temp[1]).Pow2());

            //    //Console.WriteLine("Current velResidual:  " + velResidual);
            //}
            //Console.WriteLine("Number of Iterations for translational velocity calculation:  " + vel_iteration_counter);
            //Console.WriteLine("C_v_mod:  " + C_v_mod[0]);
            #endregion

            //Crank Nicolson
            // =============================
            tempForceNew[0] = (currentIterForces_P[0][0] + currentTimeForces_P[1][0]) / 2 + massDifference * gravity[0];
            tempForceNew[1] = (currentIterForces_P[0][1] + currentTimeForces_P[1][1]) / 2 + massDifference * gravity[1];
            //temp.SetV(currentIterVel_P[0], 1);
            //temp.AccV(dt / mass_P, tempForceNew);
            temp[0] = currentTimeVel_P[1][0] + dt * tempForceNew[0] / mass_P;
            temp[1] = currentTimeVel_P[1][1] + dt * tempForceNew[1] / mass_P;

            // Save new velocity
            // =============================
            currentIterVel_P.Insert(0, temp);
            currentIterVel_P.Remove(currentIterVel_P.Last());
            currentTimeVel_P[0] = currentIterVel_P[0];

            return;
        }
        #endregion

        #region Update angular velocity
        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="D"></param>
        /// <param name="dt"></param>
        /// <param name="oldAngularVelocity"></param>
        /// <param name="Radius"></param>
        /// <param name="ParticleMass"></param>
        /// <param name="Torque"></param>
        /// <returns></returns>
        public void UpdateAngularVelocity(double dt, double rho_Fluid = 1, bool includeRotation = true, int noOfSubtimesteps = 1) {

            // save rotation of the last timestep
            if (iteration_counter_P == 1)
            {
                currentTimeRot_P[1] = currentTimeRot_P[0];
            }

            // no rotation
            // =============================
            if (includeRotation == false) {
                currentIterRot_P.Insert(0, 0.0);
                currentIterRot_P.Remove(currentIterRot_P.Last());
                return;
            }
            double temp = 0;
            double old_temp = 0;
            double tempMomentTemp = 0;
            double tempMomentNew = 0;
            double tempMomentOld = 0;
            double m_vTemp;
            double c_a = (C_v * rho_Fluid) / (rho_P + C_v * rho_Fluid);
            double c_u = 1 / (rho_P + C_v * rho_Fluid * MomentOfInertia_P);
            double newAngularVelocity = 0;
            double oldAngularVelocity = new double();
            double subtimestep;
            double rot_iteration_counter = 0;
            noOfSubtimesteps = 1;
            subtimestep = dt / noOfSubtimesteps;
            double previous_rot = new double();
            if (iteration_counter_P == 1)
            {
                previous_rot = currentIterRot_P[0];
            }

            // Benjamin Stuff
            //for (int i = 1; i <= noOfSubtimesteps; i++) {
            //    newAngularVelocity = c_a * (3 * currentIterRot_P[0] - 4 * currentIterRot_P[1] + currentIterRot_P[2]) / 2 + 0.5 * c_u * dt * (3 * currentIterTorque_P[0] - currentIterTorque_P[1]); // for 2D
            //    oldAngularVelocity = newAngularVelocity;
            //}

            for (int i = 1; i <= noOfSubtimesteps; i++) {
                //newAngularVelocity = currentIterRot_P[0] + (dt / MomentOfInertia_P) * (currentIterTorque_P[0] + currentIterTorque_P[1]) / 2; // for 2D
                newAngularVelocity = currentTimeRot_P[1] + (dt / MomentOfInertia_P) * (currentIterTorque_P[0] + currentTimeTorque_P[1]) / 2; // for 2D

                oldAngularVelocity = newAngularVelocity;

            }
            #region Müll
            //for (double rotResidual = 1; rotResidual > velResidual_ConvergenceCriterion;) {
            //    m_vTemp = c_a * (1 * temp + 3 * currentIterRot_P[0] + 3 * currentIterRot_P[1] + 1 * currentIterRot_P[2]) / (8 * dt);
            //    tempMomentTemp = (currentIterTorque_P[0]) * (c_u) + m_vTemp;
            //    tempMomentNew = (currentIterTorque_P[1]);
            //    tempMomentOld = (currentIterTorque_P[2]);
            //    old_temp = temp;
            //    temp = currentIterRot_P[0] + (1 * tempMomentTemp + 4 * tempMomentNew + 1 * tempMomentOld) * dt / 6;
            //    rot_iteration_counter += 1;
            //    if (rot_iteration_counter == MaxParticleVelIterations) {
            //        throw new ApplicationException("no convergence in particle velocity calculation");
            //    }
            //    rotResidual = Math.Sqrt((temp - old_temp).Pow2());
            //    //Console.WriteLine("Current velResidual:  " + velResidual);
            //}
            //Console.WriteLine("Number of Iterations for angular velocity calculation:  " + rot_iteration_counter);
            //if (Math.Abs(temp) < 1e-9)
            //{
            //    temp = 0;
            //}
            #endregion
            currentIterRot_P.Insert(0, newAngularVelocity);
            currentIterRot_P.Remove(currentIterRot_P.Last());
            currentTimeRot_P[0] = currentIterRot_P[0];
        }
        #endregion
        
        #region Cloning
        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            throw new NotImplementedException("Currently cloning of a particle is not available");
        }
        #endregion

        #region Update forces and torque
        /// <summary>
        /// Update forces and torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P,
            LevelSetTracker LsTrk,
            double muA) {

            if (skipForceIntegration) {
                skipForceIntegration = false;
                return;
            }
            // save forces and torque of the last timestep
            if (iteration_counter_P == 1)
            {
                //currentTimeForces_P[1] = currentTimeForces_P[0];
                //currentTimeTorque_P[1] = currentTimeTorque_P[0];
                currentTimeForces_P.Insert(1, currentTimeForces_P[0]);
                currentTimeForces_P.Remove(currentTimeForces_P.Last());
                currentTimeTorque_P.Insert(1, currentTimeTorque_P[0]);
                currentTimeTorque_P.Remove(currentTimeTorque_P.Last());
            }
            int D = LsTrk.GridDat.SpatialDimension;
            // var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UA = U.ToArray();

            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            //int RequiredOrder = LsTrk.GetXQuadFactoryHelper(momentFittingVariant).GetCachedSurfaceOrders(0).Max();
            //Console.WriteLine("Order reduction: {0} -> {1}", _RequiredOrder, RequiredOrder);

            //if (RequiredOrder > agg.HMForder)
            //    throw new ArgumentException();

            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);


            ConventionalDGField pA = null;

            //pA = P.GetSpeciesShadowField("A");
            pA = P;

            #region Force
            double[] forces = new double[D];
            for (int d = 0; d < D; d++) {
                ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D); ;
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                    // Evaluate tangential velocity to level-set surface
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);


                    for (int i = 0; i < D; i++) {
                        UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }

                    pA.Evaluate(j0, Len, Ns, pARes);

                    if (LsTrk.GridDat.SpatialDimension == 2) {

                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;
                                double scale = Normals[j, k, 0] * Math.Cos(currentIterAng_P[0]) + Normals[j, k, 1] * Math.Sin(currentIterAng_P[0]);
                                // pressure
                                switch (d) {
                                    case 0:
                                        acc += (pARes[j, k]) * Normals[j, k, 0];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0]; 
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        break;
                                    case 1:
                                        acc += (pARes[j, k]) * Normals[j, k, 1];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                        break;
                                    default:
                                        throw new NotImplementedException();
                                }

                                result[j, k] = acc;
                            }
                        }
                    } else {
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;

                                // pressure
                                switch (d) {
                                    case 0:
                                        acc += pARes[j, k] * Normals[j, k, 0];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 2] * Normals[j, k, 2];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 0] * Normals[j, k, 2];
                                        break;
                                    case 1:
                                        acc += pARes[j, k] * Normals[j, k, 1];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 2] * Normals[j, k, 2];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 1] * Normals[j, k, 2];
                                        break;
                                    case 2:
                                        acc += pARes[j, k] * Normals[j, k, 2];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 2, 2] * Normals[j, k, 2];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 2] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 2] * Normals[j, k, 1];
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


                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            forces[d] += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }
            #endregion

            #region Torque
            double torque = 0;
            ScalarFunctionEx ErrFunc2 = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                // Evaluate tangential velocity to level-set surface
                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);

                for (int i = 0; i < D; i++) {
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
                        
                        double acc = 0.0;
                        double acc2 = 0.0;
                        // Calculate the torque around a circular particle with a given radius (Paper Wan and Turek 2005)

                        acc += (pARes[j, k] * Normals[j, k, 0]);
                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                        //acc *= -Normals[j, k, 1] * this.radius_P;
                        acc *= -Normals[j, k, 1] * (this.currentIterPos_P[0][1] - tempArray[k, 1]).Abs();


                        acc2 += pARes[j, k] * Normals[j, k, 1];
                        acc2 -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                        acc2 -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                        acc2 -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                        //acc2 *= Normals[j, k, 0] * this.radius_P;
                        acc2 *= Normals[j, k, 0] * (this.currentIterPos_P[0][0] - tempArray[k, 0]).Abs();

                        result[j, k] = acc + acc2;

                    }


                }

            };

            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );
            //CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, this.cutCells_P(LsTrk));

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs2.Compile(LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        torque += ResultsOfIntegration[i, 0];
                    }
                }

            ).Execute();

            // determine underrelaxation factor (URF)
            // =============================
            temporalForces_P.Insert(0, forces);
            temporalForces_P.Remove(temporalForces_P.Last());
            temporalTorque_P.Insert(0, torque);
            temporalTorque_P.Remove(temporalTorque_P.Last());
            double underrelaxationFT = 1.0;
            int[] pre_under_R = new int[D];
            double[] temp_underR = new double[D + 1];
            double[] mu = new double[D + 1];
            double[] sigmaSquared = new double[D + 1];
            double[] muTemp = new double[D + 1];
            double[] sigmaSquaredTemp = new double[D + 1];
            for (int k = 0; k < D + 1; k++)
            {
                temp_underR[k] = underrelaxation_factor;
            }
            // first iteration, set URF to 1
            if (iteration_counter_P == 1)
            {
                for (int k = 0; k < D; k++)
                {
                    temp_underR[k] = 1;
                    for (int t = 0; t < m_HistoryLength; t++)
                    {
                        currentIterForces_P[t][k] = currentTimeForces_P[1][k];
                        currentIterTorque_P[t] = currentTimeTorque_P[1];
                    }
                }
                temp_underR[D] = 1;
            }
            // constant predefined URF
            else if (underrelaxationFT_constant == true)
            {
                for (int k = 0; k < D + 1; k++)
                {
                    temp_underR[k] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                }
            }
            // calculation of URF for adaptive underrelaxation
            else if (underrelaxationFT_constant == false)
            {
                // forces
                bool underrelaxation_ok = false;
                for (int j = 0; j < D; j++)
                {
                    underrelaxation_ok = false;
                    temp_underR[j] = underrelaxation_factor;
                    pre_under_R[j] = underrelaxationFT_exponent;
                    underrelaxationFT_exponent = 0;
                    mu[j] = 0;
                    muTemp[j] = 0;
                    for (int t = 0; t < m_HistoryLength; t++)
                    {
                        mu[j] += (currentIterForces_P[t][j]) / (m_HistoryLength);
                        muTemp[j] += (temporalForces_P[t][j]) / (m_HistoryLength);
                    }
                    sigmaSquared[j] = 0;
                    sigmaSquaredTemp[j] = 0;
                    for (int t = 0; t < m_HistoryLength; t++)
                    {
                        sigmaSquared[j] += Math.Pow((currentIterForces_P[t][j] - mu[j]) / (m_HistoryLength), 2);
                        sigmaSquaredTemp[j] += Math.Pow((temporalForces_P[t][j] - muTemp[j]) / (m_HistoryLength), 2);
                    }
                    for (int i = 0; underrelaxation_ok == false; i++)
                    {
                        if (Math.Abs(2 * temp_underR[j] * forces[j]) > Math.Abs(currentIterForces_P[0][j]) && temp_underR[j] > forceAndTorque_convergence * 100)
                        {
                            underrelaxationFT_exponent -= 1;
                            temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        }
                        //if (Math.Abs((1 - temp_underR[j]) * (forces[j] - currentIterForces_P[0][j])) < sigmaSquared[j] && temp_underR[j] > forceAndTorque_convergence * 1000)
                        //{
                        //    underrelaxationFT_exponent -= 1;
                        //    temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        //}
                        //else if (Math.Abs((1 - temp_underR[j]) * (forces[j] - currentIterForces_P[0][j])) < mu[j])
                        //{
                        //    underrelaxationFT_exponent -= 2;
                        //    temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        //}
                        //else if (Math.Abs((1 - temp_underR[j]) * (forces[j] - currentIterForces_P[0][j])) < 2 * mu[j])
                        //{
                        //    underrelaxationFT_exponent -= 2;
                        //    temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        //}
                        //if (Math.Abs(temp_underR[j] * (forces[j] - currentIterForces_P[0][j])) > Math.Abs(temp_underR[j] * (forces[j] - mu[j])) && temp_underR[j] > forceAndTorque_convergence * 1000)
                        //{
                        //    underrelaxationFT_exponent -= 1;
                        //    temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        //}
                        else
                        {
                            underrelaxation_ok = true;
                            if (Math.Abs(mu[j] - muTemp[j]) > Math.Abs(temp_underR[j] * forces[j]))
                            {
                                underrelaxationFT_exponent += 1;
                                temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                            }
                            if (underrelaxationFT_exponent >= underrelaxationFT_exponent_min && iteration_counter_P > 30)
                            {
                                underrelaxationFT_exponent = underrelaxationFT_exponent_min;
                                temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                            }
                            if (underrelaxationFT_exponent > 0)
                            {
                                underrelaxationFT_exponent = 0;
                                temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                            }
                        }
                    }
                }
                // torque
                underrelaxation_ok = false;
                temp_underR[D] = underrelaxation_factor;
                underrelaxationFT_exponent = 0;
                mu[D] = 0;
                muTemp[D] = 0;
                for (int t = 0; t < m_HistoryLength; t++)
                {
                    mu[D] += (currentIterTorque_P[t]) / (m_HistoryLength);
                    muTemp[D] += (temporalTorque_P[t]) / (m_HistoryLength);
                }
                sigmaSquared[D] = 0;
                sigmaSquaredTemp[D] = 0;
                for (int t = 0; t < m_HistoryLength; t++)
                {
                    sigmaSquared[D] += Math.Pow((currentIterTorque_P[t] - mu[D]) / (m_HistoryLength), 2);
                    sigmaSquaredTemp[D] += Math.Pow((temporalTorque_P[t] - muTemp[D]) / (m_HistoryLength), 2);
                }
                for (int i = 0; underrelaxation_ok == false; i++)
                {
                    if (Math.Abs(temp_underR[D] * torque) > Math.Abs(currentIterTorque_P[0]) && temp_underR[D] > forceAndTorque_convergence * 1000)
                    {
                        underrelaxationFT_exponent -= 1;
                        temp_underR[D] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                    }
                    else
                    {
                        underrelaxation_ok = true;
                        if (mu[D] - muTemp[D] > Math.Abs(temp_underR[D] * torque))
                        {
                            underrelaxationFT_exponent = underrelaxationFT_exponent_min;
                            temp_underR[D] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        }
                        else if (underrelaxationFT_exponent > 0)
                        {
                            underrelaxationFT_exponent = 0;
                            temp_underR[D] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        }
                    }
                }
            }
            Console.WriteLine("tempunderR[0]  " + temp_underR[0] + ", temp_underR[1]: " + temp_underR[1] + ", temp_underR[D] " + temp_underR[D]);
            Console.WriteLine("tempfForces[0]  " + forces[0] + ", temp_Forces[1]: " + forces[1] + ", tempTorque " + torque);

            // calculation of forces and torque with underrelaxation
            // =============================
            // forces
            double[] forces_underR = new double[D];
            //temp_underR[D] = 0;
            for (int i = 0; i < D; i++)
            {
                forces_underR[i] = temp_underR[i] * forces[i] + (1 - temp_underR[i]) * currentIterForces_P[0][i];
                //forces_underR[i] = (forces[i] + 257 * currentIterForces_P[0][i] + 27 * currentIterForces_P[1][i] + 272 * currentIterForces_P[2][i] + 27 * currentIterForces_P[3][i] + 216 * currentIterForces_P[4][i] + 41 * currentIterForces_P[5][i]) / 840;
            }
            // torque
            double torque_underR = temp_underR[D] * torque + (1 - temp_underR[D]) * currentIterTorque_P[0];
            // update forces and torque
            this.currentIterForces_P.Insert(0, forces_underR);
            currentIterForces_P.Remove(currentIterForces_P.Last());
            this.currentIterTorque_P.Remove(currentIterTorque_P.Last());
            this.currentIterTorque_P.Insert(0, torque_underR);
            currentTimeForces_P[0] = currentIterForces_P[0];
            currentTimeTorque_P[0] = currentIterTorque_P[0];

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
        public double ComputeParticleRe(double mu_Fluid) {
            double particleReynolds = 0;

            particleReynolds = Math.Sqrt(currentIterVel_P[0][0] * currentIterVel_P[0][0] + currentIterVel_P[0][1] * currentIterVel_P[0][1]) * 2 * radius_P * rho_P / mu_Fluid;

            return particleReynolds;
        }
        #endregion

        #region Cut cells
        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        virtual public CellMask cutCells_P(LevelSetTracker LsTrk) {

            // tolerance is very important
            var radiusTolerance = radius_P + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(currentIterAng_P[0]);

            switch (m_shape) {
                case ParticleShape.spherical:
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => (-(X[0] - currentIterPos_P[0][0]).Pow2() + -(X[1] - currentIterPos_P[0][1]).Pow2() + radiusTolerance.Pow2()) > 0);
                    break;

                case ParticleShape.elliptic:
                    double a = 3.0;
                    double b = 1.0;
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0);

                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow2() - b * ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0);
                    break;

                case ParticleShape.bean:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0);
                    break;

                case ParticleShape.squircle:
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(4) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(4)) - radiusTolerance.Pow(4)) > 0);
                    break;


                default:
                    throw new NotImplementedException("Shape is not implemented yet");

            }


            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }

        /// <summary>
        /// Gives a bool whether the particle contains a certain point or not
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        virtual public bool Contains(double[] point,LevelSetTracker LsTrk) {

            // only for squared cells
            double radiusTolerance = radius_P+2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());


            switch (m_shape) {
                case ParticleShape.spherical:
                    var distance = point.L2Distance(currentIterPos_P[0]);

                    if (distance < (radiusTolerance)) {
                        
                        return true;
                    }
                    break;

                //case ParticleShape.spherical:
                //    if (((point[0] - currentIterPos_P[0][0]).Pow2() + -(point[1] - currentIterPos_P[0][1]).Pow2() + radiusTolerance.Pow2()) > 0) {
                //        return true;
                //    }
                //    break;

                case ParticleShape.elliptic:
                    double a = 3.0;
                    double b = 1.0;
                    if (-((((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow2()) / length_P.Pow2()) + -(((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0) { 
                        return true;
                    }
                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    if (-((((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow(2) + ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow(2)).Pow2() - length_P * ((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow2() - thickness_P * ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow2()) > 0)
                        return true;
                    break;

                case ParticleShape.bean:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    if (-((((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow(2) + ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow(2)).Pow2() - a * ((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow(3) - b * ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow2()) > 0)
                        return true;
                    break;

                case ParticleShape.squircle:
                    if (-((((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow(4) + ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow(4)) - radiusTolerance.Pow(4)) > 0)
                        return true;
                    break;


                default:
                    throw new NotImplementedException("Shape is not implemented yet");
            }
            return false;
        }
        #endregion

        #region Particle shape
        public enum ParticleShape {
            spherical = 0,

            elliptic = 1,

            hippopede = 2,

            bean = 3,

            squircle = 4
        }
        #endregion
    }

}

