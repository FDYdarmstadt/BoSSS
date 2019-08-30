using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class ParticleMotion {

        private const int historyLength = 4;
        protected static int spatialDim = 2;
        [NonSerialized]
        readonly internal FSI_Auxillary Aux = new FSI_Auxillary();
        [NonSerialized]
        readonly private ParticleForceIntegration ForceIntegration = new ParticleForceIntegration();

        public ParticleMotion(double[] gravity) {
            m_Gravity = gravity;
            for (int i = 0; i < historyLength; i++) {
                position.Add(new double[spatialDim]);
                angle.Add(new double());
                translationalVelocity.Add(new double[spatialDim]);
                translationalAcceleration.Add(new double[spatialDim]);
                rotationalVelocity.Add(new double());
                rotationalAcceleration.Add(new double());
                hydrodynamicForces.Add(new double[spatialDim]);
                hydrodynamicTorque.Add(new double());
            }
        }

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public double[] m_Gravity;

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public double particleMass = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public double particleMomentOfInertia = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> position = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> angle = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> translationalVelocity = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> rotationalVelocity = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> translationalAcceleration = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> rotationalAcceleration = new List<double>();


        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public List<double[]> hydrodynamicForces = new List<double[]>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double[] forcesPrevIteration = new double[spatialDim];

        /// <summary>
        /// The Torque acting on the particle in the current time step.
        /// </summary>
        public List<double> hydrodynamicTorque = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double torquePrevIteration = new double();

        /// <summary>
        /// Saves position, angle, acceleration and velocity of the last timestep
        /// </summary>
        public void SaveDataOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(position);
            Aux.SaveValueOfLastTimestep(angle);
            Aux.SaveMultidimValueOfLastTimestep(translationalVelocity);
            Aux.SaveValueOfLastTimestep(rotationalVelocity);
            Aux.SaveMultidimValueOfLastTimestep(translationalAcceleration);
            Aux.SaveValueOfLastTimestep(rotationalAcceleration);
            Aux.SaveMultidimValueOfLastTimestep(hydrodynamicForces);
            Aux.SaveValueOfLastTimestep(hydrodynamicTorque);
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public virtual void CalculateParticlePosition(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + (translationalVelocity[0][d] + 4 * translationalVelocity[1][d] + translationalVelocity[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        public virtual void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + translationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public virtual void CalculateParticleAngle(double dt) {
            angle[0] = angle[1] + (rotationalVelocity[0] + 4 * rotationalVelocity[1] + rotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        public virtual void CalculateParticleAngle(double dt, double collisionTimestep) {
            angle[0] = angle[1] + rotationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public virtual void CalculateAcceleration(double dt, double[] force, double torque) {
            // Translation
            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = force[d] / particleMass;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");
            // Rotation
            rotationalAcceleration[0] = torque / particleMomentOfInertia;
            Aux.TestArithmeticException(rotationalAcceleration[0], "particle rotational acceleration");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        public virtual void CalculateTranslationalVelocity(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + (translationalAcceleration[0][d] + 4 * translationalAcceleration[1][d] + translationalAcceleration[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        public virtual void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + translationalAcceleration[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        public virtual void CalculateAngularVelocity(double dt) {
            rotationalVelocity[0] = rotationalVelocity[1] + (rotationalAcceleration[0] + 4 * rotationalAcceleration[1] + rotationalAcceleration[2]) * dt / 6;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        public virtual void CalculateAngularVelocity(double dt, double collisionTimestep) {
            rotationalVelocity[0] = rotationalVelocity[1] + rotationalAcceleration[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public virtual void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt = 0) {
            double[] tempForces = CalculateHydrodynamicForces(U, P, LsTrk, CutCells_P, muA, relativeParticleMass);
            double tempTorque = CalculateHydrodynamicTorque(U, P, LsTrk, CutCells_P, muA);
            HydrodynamicsPostprocessing(tempForces, tempTorque);
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        protected virtual double[] CalculateHydrodynamicForces(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt = 0) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double[] tempForces = ForcesIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Force_MPISum(ref tempForces);
            for (int d = 0; d < spatialDim; d++) {
                tempForces[d] += relativeParticleMass * m_Gravity[d];
            }
            return tempForces;
        }

        protected double[] ForcesIntegration(SinglePhaseField[] UA, ConventionalDGField pA, LevelSetTracker LsTrk, CellMask CutCells_P, int RequiredOrder, double FluidViscosity) {
            double[] tempForces = new double[spatialDim];
            double[] IntegrationForces = tempForces.CloneAs();
            for (int d = 0; d < spatialDim; d++) {
                void ErrFunc(int CurrentCellID, int Length, NodeSet Ns, MultidimensionalArray result) {

                    int K = result.GetLength(1);
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Length, K, spatialDim, spatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Length, K);
                    MultidimensionalArray Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, CurrentCellID, Length);
                    for (int i = 0; i < spatialDim; i++) {
                        UA[i].EvaluateGradient(CurrentCellID, Length, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    pA.Evaluate(CurrentCellID, Length, Ns, pARes);
                    for (int j = 0; j < Length; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = ForceIntegration.CalculateStressTensor(Grad_UARes, pARes, Normals, FluidViscosity, k, j, spatialDim, d);
                        }
                    }
                }
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, CutCells_P);
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat, cqs.Compile(LsTrk.GridDat, RequiredOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        IntegrationForces[d] = ForceTorqueSummationWithNeumaierArray(IntegrationForces[d], ResultsOfIntegration, Length);
                    }
                ).Execute();
            }
            return tempForces = IntegrationForces.CloneAs();
        }

        protected void Force_MPISum(ref double[] forces) {
            double[] stateBuffer = forces.CloneAs();
            double[] globalStateBuffer = stateBuffer.MPISum();
            for (int d = 0; d < spatialDim; d++) {
                forces[d] = globalStateBuffer[d];
            }
        }

        protected virtual double CalculateHydrodynamicTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double dt = 0) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double tempTorque = TorqueIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Torque_MPISum(ref tempTorque);
            return tempTorque;
        }

        protected double TorqueIntegration(SinglePhaseField[] UA, ConventionalDGField pA, LevelSetTracker LsTrk, CellMask CutCells_P, int RequiredOrder, double FluidViscosity) {
            double tempTorque = new double();
            void ErrFunc2(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1);
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, spatialDim, spatialDim); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                for (int i = 0; i < spatialDim; i++) {
                    UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }
                pA.Evaluate(j0, Len, Ns, pARes);
                for (int j = 0; j < Len; j++) {
                    MultidimensionalArray Ns_Global = Ns.CloneAs();
                    LsTrk.GridDat.TransformLocal2Global(Ns, Ns_Global, j0 + j);
                    for (int k = 0; k < K; k++) {
                        result[j, k] = ForceIntegration.CalculateTorqueFromStressTensor2D(Grad_UARes, pARes, Normals, Ns_Global, FluidViscosity, k, j, position[0]);
                    }
                }
            }
            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, CutCells_P);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat, cqs2.Compile(LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    tempTorque = ForceTorqueSummationWithNeumaierArray(tempTorque, ResultsOfIntegration, Length);
                }
            ).Execute();
            return tempTorque;
        }

        protected void Torque_MPISum(ref double torque) {
            double stateBuffer = torque;
            double globalStateBuffer = stateBuffer.MPISum();
            torque = globalStateBuffer;
        }

        protected virtual void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque) {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = tempForces[d];
            }
            hydrodynamicTorque[0] = tempTorque;
            Aux.TestArithmeticException(hydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(hydrodynamicTorque[0], "hydrodynamic torque");
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
        private double ForceTorqueSummationWithNeumaierArray(double ResultVariable, MultidimensionalArray Summands, double Length) {
            double sum = ResultVariable;
            double naiveSum;
            double c = 0.0;
            for (int i = 0; i < Length; i++) {
                naiveSum = sum + Summands[i, 0];
                if (Math.Abs(sum) >= Math.Abs(Summands[i, 0])) {
                    c += (sum - naiveSum) + Summands[i, 0];
                }
                else {
                    c += (Summands[i, 0] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            return sum + c;
        }
    }
}
