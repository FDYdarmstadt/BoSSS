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

        public Particle(int Dim =2, int HistoryLength=4, double[] startPos = null, double startAngl = 0.0, ParticleShape shape = ParticleShape.spherical) {

            if (startPos == null) {
                if (Dim == 2) {
                    startPos = new double[] { 0.0, 0.0 };
                } else {
                    startPos = new double[] { 0.0, 0.0, 0.0 };
                }
            }

            m_Dim = Dim;
            m_HistoryLength = HistoryLength;
            m_shape = shape;

            for (int i = 0; i < HistoryLength; i++) {
                currentPos_P.Add(new double[Dim]);
                currentAng_P.Add(new double());
                vel_P.Add(new double[Dim]);
                rot_P.Add(new double());
                forces_P.Add(new double[Dim]);
                torque_P.Add(new double());
            }

            currentPos_P[0] = startPos;
            //From degree to radiant
            currentAng_P[0] = startAngl * 2 * Math.PI / 360;

            UpdateLevelSetFunction();

        }

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
        /// Initial position of the particle (the center of mass)
        /// </summary>
        [DataMember]
        public List<double[]> currentPos_P = new List<double[]>();

        /// <summary>
        /// Initial angle of the particle (the center of mass)
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
                        area = (3.0 * radius_P) * (radius_P* 1.0) * Math.PI;
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
                        moment = (1 / 4.0) * (mass_P * (3.0 * 3.0 + 1.0 * 1.0) * radius_P * radius_P);
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

        public void UpdateLevelSetFunction() {

            double a;
            double b;
            double alpha;


            // For angle and shift positions
            // for X[0] choose ((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha))
            // for X[1] choose ((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha))

            switch (m_shape) {

                case ParticleShape.spherical:
                    phi_P = (X, t) => -(X[0] - currentPos_P[0][0]).Pow2() + -(X[1] - currentPos_P[0][1]).Pow2() + radius_P.Pow2();
                    break;

                case ParticleShape.elliptic:
                    alpha = -(currentAng_P[0]);
                    a = 3.0;
                    b = 1.0;
                    phi_P = (X, t) => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / a.Pow2()) + -(((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2() / b.Pow2()) + radius_P.Pow2();
                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radius_P.Pow2();
                    b = 1.0 * radius_P.Pow2();
                    alpha = -(currentAng_P[0]);
                    // hippopede
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
        public void UpdateTransVelocity(double dt, double rho_Fluid = 1, bool includeTranslation = true, bool LiftAndDrag = true) {

            double[] temp = new double[2];

            if (includeTranslation == false) {

                vel_P[0][0] = 0.0;
                vel_P[0][1] = 0.0;

                return;
            }
            // Uncommented parts are for introducing virtual mass force in order to stabilize calculations for light spheres. Approach taken from
            // "Schwarz et al. - 2015 A temporal discretization scheme to compute the motion of light particles in viscous flows by an immersed boundary", paper is saved in 
            // Z:\Stange\HiWi\Arbeitspakete

            double[] tempForceNew = new double[2];
            double[] tempForceOld = new double[2];
            //double[] tempForce = new double[2];
            double[] gravity = new double[2];
            gravity[0] = 0;
            gravity[1] = -98.0;
            double C_v = 0.5;
            // In case we want to quickly adapt to all possible fluid densities
            double[] f_vNew = new double[2];
            double[] f_vOld = new double[2];
            double c_a = (C_v * rho_Fluid) / (rho_P + C_v * rho_Fluid);


            double massDifference = (rho_P - rho_Fluid) * (area_P);

            //Console.WriteLine("Mass difference:    " + massDifference);

            //Console.WriteLine("Mass difference times gravity:    " + (massDifference*gravity[1]));

            // VIRTUAL MASS MODEL
            //for (int i = 0; i < 2; i++)
            //{
            //    tempforcenew[i] = forces_p[0][i] + massdifference * gravity[i];
            //    tempforceold[i] = forces_p[1][i] + massdifference * gravity[i];
            //    f_vnew[i] = c_a * (3 * vel_p[0][i] - 4 * vel_p[1][i] + vel_p[2][i]) / dt;
            //    f_vold[i] = c_a * (3 * vel_p[1][i] - 4 * vel_p[2][i] + vel_p[3][i]) / dt;
            //    tempforcenew[i] = tempforcenew[i] / ((math.pi * radius_p * radius_p) * (rho_p + c_v * rho_fluid)) + f_vnew[i];
            //    tempforceold[i] = tempforceold[i] / ((math.pi * radius_p * radius_p) * (rho_p + c_v * rho_fluid)) + f_vold[i];
            //    temp[i] = vel_p[0][i] + (3 * tempforcenew[i] - tempforceold[i]) * dt / 2;
            //}

            

            // modell 1
            tempForceNew[0] = 0.5 * (forces_P[1][0] + forces_P[0][0]) + massDifference * gravity[0];
            tempForceNew[1] = 0.5 * (forces_P[1][1] + forces_P[0][1]) + massDifference * gravity[1];
            temp.SetV(vel_P[0], 1);
            temp.AccV(dt / mass_P, tempForceNew);


            vel_P.Insert(0, temp);

            vel_P.Remove(vel_P.Last());

            return;
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
        public void UpdateAngularVelocity(double dt, bool includeRotation = true, int noOfSubtimesteps = 1) {
            if (includeRotation == false) {
                rot_P.Insert(0, 0.0);
                rot_P.Remove(rot_P.Last());
                return;
            }

            double newAngularVelocity = new double();
            double oldAngularVelocity = new double();
            double subtimestep;

            noOfSubtimesteps = 1;
            subtimestep = dt / noOfSubtimesteps;

            for (int i = 1; i <= noOfSubtimesteps; i++) {
                newAngularVelocity = rot_P[0] + (dt / MomentOfInertia_P) * (torque_P[0] + torque_P[1]); // for 2D

                oldAngularVelocity = newAngularVelocity;

            }

            rot_P.Insert(0, newAngularVelocity);

            rot_P.Remove(rot_P.Last());
        }

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            throw new NotImplementedException("Currently cloning of a particle is not available");
        }

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
                                        acc += pARes[j, k] * Normals[j, k, 0];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        break;
                                    case 1:
                                        acc += pARes[j, k] * Normals[j, k, 1];
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


            this.forces_P.Insert(0, forces);
            forces_P.Remove(forces_P.Last());
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

                        acc += pARes[j, k] * Normals[j, k, 0];
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
            this.torque_P.Remove(torque_P.Last());
            this.torque_P.Insert(0, torque);

            #endregion

        }

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
                    cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / a.Pow2()) + -(((X[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2() / b.Pow2()) + radiusTolerance.Pow2() > 0);

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


            double alpha = -(currentAng_P[0]);

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
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (point[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / a.Pow2()) + -(((point[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (point[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2() / b.Pow2()) + radiusTolerance.Pow2() > 0) {
                        return true;
                    }
                    break;

                case ParticleShape.hippopede:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (point[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((point[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (point[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((point[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (point[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow2() - b * ((point[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (point[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0)
                        return true;
                    break;

                case ParticleShape.bean:
                    a = 4.0 * radiusTolerance.Pow2();
                    b = 1.0 * radiusTolerance.Pow2();
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (point[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((point[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (point[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((point[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (point[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((point[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (point[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0)
                        return true;
                    break;

                case ParticleShape.squircle:
                    if (-((((point[0] - currentPos_P[0][0]) * Math.Cos(alpha) - (point[1] - currentPos_P[0][1]) * Math.Sin(alpha)).Pow(4) + ((point[0] - currentPos_P[0][0]) * Math.Sin(alpha) + (point[1] - currentPos_P[0][1]) * Math.Cos(alpha)).Pow(4)) - radiusTolerance.Pow(4)) > 0)
                        return true;
                    break;


                default:
                    throw new NotImplementedException("Shape is not implemented yet");
            }
            return false;
        }

        public enum ParticleShape {
            spherical = 0,

            elliptic = 1,

            hippopede = 2,

            bean = 3,

            squircle = 4
        }

    }
}

