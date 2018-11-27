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
using ilPSP.Utils;

namespace BoSSS.Application.FSI_Solver {


    public static class IBMMover
    {

        /// <summary>
        /// Calculate the new translational velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="D"></param>
        /// <param name="dt"></param>
        /// <param name="oldTransVelocity"></param>
        /// <param name="particleMass"></param>
        /// <param name="force"></param>
        /// <returns></returns>
        //public static double[] GetTransVelocity(double phystime, double dt, double[] oldTransVelocity, double radius, double particleRho, double[] force, double[] oldforce, bool includeTranslation)
        public static double[] GetTransVelocity(double phystime, double dt, double[] oldTransVelocity, double[] TransVelocityN2, double[] TransVelocityN3, double[] TransVelocityN4, double radius, double particleRho, double[] force, double[] oldforce, bool includeTranslation)
        {
            double[] temp = new double[2];
           
            if (includeTranslation == false)
            {

                temp[0] = 0;
                temp[1] = 0;

                return temp;
            }
            // Uncommented parts are for introducing virtual mass force in order to stabilize calculations for light spheres. Approach taken from
            // "Schwarz et al. - 2015 A temporal discretization scheme to compute the motion of light particles in viscous flows by an immersed boundary", paper is saved in 
            // Z:\Stange\HiWi\Arbeitspakete

            double[] tempForceNew = new double[2];
            double[] tempForceOld = new double[2];
            //double[] tempForce = new double[2];
            double[] gravity = new double[2];
            gravity[0] = 0;
            gravity[1] = -981;
            double C_v = 0.5;
            // In case we want to quickly adapt to all possible fluid densities
            double rho_A = 1;
            double[] f_vNew = new double[2];
            double[] f_vOld = new double[2];
            double c_a = (C_v * rho_A) / (particleRho + C_v * rho_A);
            
            double particleMass = particleRho * Math.PI * radius * radius;
            double massDifference = (particleRho - rho_A) * (Math.PI * radius * radius);

            for (int i = 0; i < 2; i++)
            {
                tempForceNew[i] = force[i] + massDifference * gravity[i];
                tempForceOld[i] = oldforce[i] + massDifference * gravity[i];
                f_vNew[i] = c_a * (3 * oldTransVelocity[i] - 4 * TransVelocityN2[i] + TransVelocityN3[i]) / dt;
                f_vOld[i] = c_a * (3 * TransVelocityN2[i] - 4 * TransVelocityN3[i] + TransVelocityN4[i]) / dt;
                tempForceNew[i] = tempForceNew[i] / ((Math.PI * radius * radius) * (particleRho + C_v * rho_A)) + f_vNew[i];
                tempForceOld[i] = tempForceOld[i] / ((Math.PI * radius * radius) * (particleRho + C_v * rho_A)) + f_vOld[i];
                temp[i] = oldTransVelocity[i] + (3 * tempForceNew[i] - tempForceOld[i]) * dt / 2;
            }
            //double massDifference = Math.Abs((Math.PI * radius * radius) - particleMass);


            //tempForce[0] = 0.5 * (oldforce[0] + force[0]);
            //tempForce[1] = 0.5 * (oldforce[1] + force[1]);
            //tempForce.SetV(force, 0.5);
            //tempForce.SetV(oldforce,0.5);
            //temp.SetV(oldTransVelocity, 1);

            //tempForce.AccV(massDifference, gravity);

            //temp.AccV(dt / particleMass, tempForce);

            return temp;
        }

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
        public static double GetAngularVelocity(double dt, double oldAngularVelocity, double Radius, double ParticleMass, double Torque, double oldTorque, bool includeRotation)
        {
            if (includeRotation == false)
                return 0;

            double newAngularVelocity = new double();
            int noOfSubtimesteps;
            double subtimestep;

            noOfSubtimesteps = 1;
            subtimestep = dt / noOfSubtimesteps;

            for (int i = 1; i <= noOfSubtimesteps; i++)
            {         
                 newAngularVelocity = oldAngularVelocity + (dt / (ParticleMass * Radius * Radius)) * (Torque+oldTorque); // for 2D

                oldAngularVelocity = newAngularVelocity;

            }

            return newAngularVelocity;
        }

        /// <summary>
        /// Get new postion of the particle with given Velocities.
        /// </summary>
        /// <param name="D"></param>
        /// <param name="dt"></param>
        /// <param name="oldVelocity"></param>
        /// <returns></returns>
        public static double[] MoveCircularParticle(double dt, double[] currentVelocity, double[] oldPosition)
        {

            oldPosition.AccV(dt, currentVelocity);

            return oldPosition;
        }


        /// <summary>
        /// Calculating the particle reynolds number according to paper Turek and testcase ParticleUnderGravity
        /// </summary>
        /// <param name="currentVelocity"></param>
        /// <param name="Radius"></param>
        /// <param name="particleDensity"></param>
        /// <param name="viscosity"></param>
        /// <returns></returns>
        public static double ComputeParticleRe(double[] currentVelocity, double Radius, double pDensity)
        {
            double particleReynolds = 0;

            particleReynolds = Math.Sqrt(currentVelocity[0]* currentVelocity[0]+ currentVelocity[1] * currentVelocity[1])*2*Radius*pDensity/0.1;

            return particleReynolds;
        }
    }
}
