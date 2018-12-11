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

namespace BoSSS.Application.FSI_Solver {

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
                currentPos_P.Add(new double[Dim]);
                currentAng_P.Add(new double());
                vel_P.Add(new double[Dim]);
                rot_P.Add(new double());
                forces_P.Add(new double[Dim]);
                torque_P.Add(new double());
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
            currentPos_P[0] = startPos;
            //From degree to radiant
            currentAng_P[0] = startAngl * 2 * Math.PI / 360;
            //vel_P[0][0] = 2e-8;

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
        public int underrelaxationFT_exponent = 1;
        public double underrelaxation_factor = 1;

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
        public List<double[]> currentPos_P = new List<double[]>();

        /// <summary>
        /// Current angle of the particle (the center of mass)
        /// </summary>
        [DataMember]
        public List<double> currentAng_P = new List<double>();

        /// <summary>
        /// Current translational Velocity of the particle
        /// </summary>
        [DataMember]
        public List<double[]> vel_P = new List<double[]>();

        /// <summary>
        /// Current rotational Velocity of the particle
        /// </summary>
        [DataMember]
        public List<double> rot_P = new List<double>();

        /// <summary>
        /// Force acting on the particle
        /// </summary>
        [DataMember]
        public List<double[]> forces_P = new List<double[]>();

        /// <summary>
        /// Torque acting on the particle
        /// </summary>
        [DataMember]
        public List<double> torque_P = new List<double>();

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
        public double MaxParticleVelIterations = 100;

        /// <summary>
        /// Active stress on the current particle
        /// </summary>
        public double active_stress_P 
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
            get {
                double mass;
                switch (m_shape) {
                    case ParticleShape.spherical:
                        mass = area_P * rho_P;
                        break;

                    case ParticleShape.elliptic:
                        mass = area_P * rho_P;
                        break;

                    case ParticleShape.hippopede:
                        mass = area_P * rho_P;
                        break;

                    case ParticleShape.bean:
                        mass = area_P * rho_P;
                        break;

                    case ParticleShape.squircle:
                        mass = area_P * rho_P;
                        break;

                    default:

                        throw new NotImplementedException("");
                }

                return mass;
            }

        }

        /// <summary>
        /// Area of the current particle
        /// </summary>
        [DataMember]
        public double area_P {
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
        public double MomentOfInertia_P {
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
        public void CleanHistory() {
            if (currentPos_P.Count > m_HistoryLength) {

                for (int j = currentPos_P.Count; j > m_HistoryLength; j--) {

                    int tempPos = j - 1;

                    currentPos_P.RemoveAt(tempPos);
                    vel_P.RemoveAt(tempPos);
                    forces_P.RemoveAt(tempPos);
                    rot_P.RemoveAt(tempPos);
                    torque_P.RemoveAt(tempPos);

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
            var tempPos = currentPos_P[0].CloneAs();
            tempPos.AccV(dt, vel_P[0]);

            currentPos_P.Insert(0, tempPos);
            currentPos_P.Remove(currentPos_P.Last());

            // Angle
            var tempAng = currentAng_P[0];
            tempAng += dt * rot_P[0];
            currentAng_P.Insert(0, tempAng);
            currentAng_P.Remove(currentAng_P.Last());

            //Console.WriteLine("Current angle speed is " + rot_P[0]*360/Math.PI);
            //Console.WriteLine("Current angle is " + currentAng_P[0] + " rad");

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
                    phi_P = (X, t) => -(X[0] - currentPos_P[0][0]).Pow2() + -(X[1] - currentPos_P[0][1]).Pow2() + radius_P.Pow2();
                    break;

                case ParticleShape.elliptic:
                    alpha = -(currentAng_P[0]);
                    a = 3.0;//no longer necessary
                    b = 1.0;//no longer necessary
                    phi_P = (X, t) => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radius_P.Pow2();
                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radius_P.Pow2();
                    b = 1.0 * radius_P.Pow2();
                    alpha = -(currentAng_P[0]);
                    phi_P = (X, t) => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2() - b * ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2());
                    break;


                case ParticleShape.bean:
                    a = 3.0 * radius_P.Pow2();
                    b = 1.0 * radius_P.Pow2();
                    alpha = -(currentAng_P[0]);
                    phi_P = (X, t) => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2());
                    break;

                case ParticleShape.squircle:
                    alpha = -(currentAng_P[0]);
                    phi_P = (X, t) => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(4) + ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(4)) - radius_P.Pow(4));
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

            // no translation
            // =============================
            if (includeTranslation == false)
            {
                vel_P[0][0] = 0.0;
                vel_P[0][1] = 0.0;
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

            if (iteration_counter_P == 0)
            {
                previous_vel = vel_P[0];
            }
            Console.WriteLine("Previous Velocity:  " + previous_vel[0]);
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
            //    f_vNew[i] = c_a * (3 * vel_P[0][i] - 4 * vel_P[1][i] + vel_P[2][i]) / (2 * dt);
            //    f_vOld[i] = c_a * (3 * vel_P[1][i] - 4 * vel_P[2][i] + vel_P[3][i]) / (2 * dt);
            //    tempForceNew[i] = (forces_P[0][i] + massDifference * gravity[i]) * (c_u) + f_vNew[i];
            //    tempForceOld[i] = (forces_P[1][i] + massDifference * gravity[i]) * (c_u) + f_vOld[i];
            //    temp[i] = vel_P[0][i] + (3 * tempForceNew[i] - tempForceOld[i]) * dt / 2;
            //}

            // implicit Adams Moulton (modified)
            //for (double velResidual = 1; velResidual > velResidual_ConvergenceCriterion;)
            //{
            //    for (int i = 0; i < 2; i++)
            //    {
            //        C_v_mod[i] = 0.1;// * Math.Abs(forces_P[0][i] / (forces_P[0][i] + forces_P[1][i] + 1e-30));
            //        c_a[i] = (C_v_mod[i] * rho_Fluid) / (rho_P + C_v_mod[i] * rho_Fluid);
            //        c_u[i] = 1 / (area_P * (rho_P + C_v_mod[i] * rho_P));
            //        f_vTemp[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * temp[i] - 18 * vel_P[0][i] + 9 * vel_P[1][i] - 2 * vel_P[2][i]) / (8 * dt);
            //        f_vNew[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * vel_P[0][i] - 18 * vel_P[1][i] + 9 * vel_P[2][i] - 2 * vel_P[3][i]) / (6 * dt);
            //        f_vOld[i] = (C_v_mod[i]) / (1 + C_v_mod[i]) * (11 * vel_P[1][i] - 18 * vel_P[2][i] + 9 * vel_P[3][i] - 2 * vel_P[4][i]) / (6 * dt);
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

            //Crank Nicolson
            // =============================
            tempForceNew[0] = (forces_P[0][0] + forces_P[1][0]) / 2+ massDifference * gravity[0];
            tempForceNew[1] = (forces_P[1][1] + forces_P[0][1]) / 2 + massDifference * gravity[1];
            //temp.SetV(vel_P[0], 1);
            //temp.AccV(dt / mass_P, tempForceNew);
            temp[0] = previous_vel[0] + dt / mass_P * tempForceNew[0];
            temp[1] = previous_vel[1] + dt / mass_P * tempForceNew[1];
            
            // Save new velocity
            // =============================

            vel_P.Insert(0, temp);
            vel_P.Remove(vel_P.Last());

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

            // no rotation
            // =============================
            if (includeRotation == false) {
                rot_P.Insert(0, 0.0);
                rot_P.Remove(rot_P.Last());
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

            // Benjamin Stuff
            //for (int i = 1; i <= noOfSubtimesteps; i++) {
            //    newAngularVelocity = c_a * (3 * rot_P[0] - 4 * rot_P[1] + rot_P[2]) / 2 + 0.5 * c_u * dt * (3 * torque_P[0] - torque_P[1]); // for 2D
            //    oldAngularVelocity = newAngularVelocity;
            //}

            for (int i = 1; i <= noOfSubtimesteps; i++) {
                newAngularVelocity = rot_P[0] + (dt / MomentOfInertia_P) * (torque_P[0] + torque_P[1]); // for 2D

                oldAngularVelocity = newAngularVelocity;

            }

            //for (double rotResidual = 1; rotResidual > velResidual_ConvergenceCriterion;) {
            //    m_vTemp = c_a * (1 * temp + 3 * rot_P[0] + 3 * rot_P[1] + 1 * rot_P[2]) / (8 * dt);
            //    tempMomentTemp = (torque_P[0]) * (c_u) + m_vTemp;
            //    tempMomentNew = (torque_P[1]);
            //    tempMomentOld = (torque_P[2]);
            //    old_temp = temp;
            //    temp = rot_P[0] + (1 * tempMomentTemp + 4 * tempMomentNew + 1 * tempMomentOld) * dt / 6;
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
            rot_P.Insert(0, newAngularVelocity);
            rot_P.Remove(rot_P.Last());
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
                        acc *= -Normals[j, k, 1] * (this.currentPos_P[0][1] - tempArray[k, 1]).Abs();


                        acc2 += pARes[j, k] * Normals[j, k, 1];
                        acc2 -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                        acc2 -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                        acc2 -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                        //acc2 *= Normals[j, k, 0] * this.radius_P;
                        acc2 *= Normals[j, k, 0] * (this.currentPos_P[0][0] - tempArray[k, 0]).Abs();

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
            
            double underrelaxationFT = 1.0;
            double[] temp_underR = new double[D + 1];
            for (int k = 0; k < D+1; k++)
            {
                temp_underR[k] = underrelaxation_factor;
            }
            if (iteration_counter_P == 0)
            {
                underrelaxationFT = 1;
            }
            else if (underrelaxationFT_constant == true)
            {
                underrelaxationFT = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
            }
            else if (underrelaxationFT_constant == false)
            {
                //double[] temp_underR = new double[D + 1];
                bool underrelaxation_ok = false;
                underrelaxationFT_exponent = 1;
                for (int j = 0; j < D; j++)
                {
                    underrelaxation_ok = false;
                    temp_underR[j] = underrelaxation_factor;
                    for (int i = 0; underrelaxation_ok == false; i++)
                    {
                        if (Math.Abs(temp_underR[j] * forces[j]) > Math.Abs(forces_P[0][j]))
                        {
                            underrelaxationFT_exponent -= 1;
                            temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        }
                        else
                        {
                            underrelaxation_ok = true;
                            if (underrelaxationFT_exponent > -0)
                            {
                                underrelaxationFT_exponent = -0;
                                temp_underR[j] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                            }
                        }
                    }
                }
                underrelaxation_ok = false;
                temp_underR[D] = underrelaxation_factor;
                for (int i = 0; underrelaxation_ok == false; i++)
                {
                    if (Math.Abs(temp_underR[D] * torque) > Math.Abs(torque_P[0]))
                    {
                        underrelaxationFT_exponent -= 1;
                        temp_underR[D] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                    }
                    else
                    {
                        underrelaxation_ok = true;
                        if (underrelaxationFT_exponent > -0)
                        {
                            underrelaxationFT_exponent = -0;
                            temp_underR[D] = underrelaxation_factor * Math.Pow(10, underrelaxationFT_exponent);
                        }
                    }
                }
            }

            double[] forces_underR = new double[D];
            for (int i = 0; i < D; i++)
            {
                forces_underR[i] = temp_underR[i] * forces[i] + (1 - temp_underR[i]) * forces_P[0][i];
            }
            double torque_underR = temp_underR[D] * torque + (1 - temp_underR[D]) * torque_P[0];
            this.forces_P.Insert(0, forces_underR);
            forces_P.Remove(forces_P.Last());
            this.torque_P.Remove(torque_P.Last());
            this.torque_P.Insert(0, torque_underR);

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

            particleReynolds = Math.Sqrt(vel_P[0][0] * vel_P[0][0] + vel_P[0][1] * vel_P[0][1]) * 2 * radius_P * rho_P / mu_Fluid;

            return particleReynolds;
        }
        #endregion

        #region Cut cells
        /// <summary>
        /// get cut cells describing the boundary of this particle
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        public CellMask cutCells_P(LevelSetTracker LsTrk) {

            // tolerance is very important
            var radiusTolerance = radius_P + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(currentAng_P[0]);

            switch (m_shape) {
                case ParticleShape.spherical:
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => (-(X[0] - currentPos_P[0][0]).Pow2() + -(X[1] - currentPos_P[0][1]).Pow2() + radiusTolerance.Pow2()) > 0);
                    break;

                case ParticleShape.elliptic:
                    double a = 3.0;
                    double b = 1.0;
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0);

                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2() - b * ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0);
                    break;

                case ParticleShape.bean:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0);
                    break;

                case ParticleShape.squircle:
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(4) + ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(4)) - radiusTolerance.Pow(4)) > 0);
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
        public bool Contains(double[] point,LevelSetTracker LsTrk) {

            // only for squared cells
            double radiusTolerance = radius_P+2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());


            switch (m_shape) {
                case ParticleShape.spherical:
                    var distance = point.L2Distance(currentPos_P[0]);

                    if (distance < (radiusTolerance)) {
                        
                        return true;
                    }
                    break;

                //case ParticleShape.spherical:
                //    if (((point[0] - currentPos_P[0][0]).Pow2() + -(point[1] - currentPos_P[0][1]).Pow2() + radiusTolerance.Pow2()) > 0) {
                //        return true;
                //    }
                //    break;

                case ParticleShape.elliptic:
                    double a = 3.0;
                    double b = 1.0;
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(currentAng_P[0]) - (point[1] - currentPos_P[0][1]) * Math.Sin(currentAng_P[0])).Pow2()) / length_P.Pow2()) + -(((point[0] - currentPos_P[0][0]) * Math.Sin(currentAng_P[0]) + (point[1] - currentPos_P[0][1]) * Math.Cos(currentAng_P[0])).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0) { 
                        return true;
                    }
                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(currentAng_P[0]) - (point[1] - currentPos_P[0][1]) * Math.Sin(currentAng_P[0])).Pow(2) + ((point[0] - currentPos_P[0][0]) * Math.Sin(currentAng_P[0]) + (point[1] - currentPos_P[0][1]) * Math.Cos(currentAng_P[0])).Pow(2)).Pow2() - length_P * ((point[0] - currentPos_P[0][0]) * Math.Cos(currentAng_P[0]) - (point[1] - currentPos_P[0][1]) * Math.Sin(currentAng_P[0])).Pow2() - thickness_P * ((point[0] - currentPos_P[0][0]) * Math.Sin(currentAng_P[0]) + (point[1] - currentPos_P[0][1]) * Math.Cos(currentAng_P[0])).Pow2()) > 0)
                        return true;
                    break;

                case ParticleShape.bean:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(currentAng_P[0]) - (point[1] - currentPos_P[0][1]) * Math.Sin(currentAng_P[0])).Pow(2) + ((point[0] - currentPos_P[0][0]) * Math.Sin(currentAng_P[0]) + (point[1] - currentPos_P[0][1]) * Math.Cos(currentAng_P[0])).Pow(2)).Pow2() - a * ((point[0] - currentPos_P[0][0]) * Math.Cos(currentAng_P[0]) - (point[1] - currentPos_P[0][1]) * Math.Sin(currentAng_P[0])).Pow(3) - b * ((point[0] - currentPos_P[0][0]) * Math.Sin(currentAng_P[0]) + (point[1] - currentPos_P[0][1]) * Math.Cos(currentAng_P[0])).Pow2()) > 0)
                        return true;
                    break;

                case ParticleShape.squircle:
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(currentAng_P[0]) - (point[1] - currentPos_P[0][1]) * Math.Sin(currentAng_P[0])).Pow(4) + ((point[0] - currentPos_P[0][0]) * Math.Sin(currentAng_P[0]) + (point[1] - currentPos_P[0][1]) * Math.Cos(currentAng_P[0])).Pow(4)) - radiusTolerance.Pow(4)) > 0)
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

