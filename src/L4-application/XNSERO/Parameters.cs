using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSERO_Solver {

    /// <summary>
    /// Level set velocity of a rigid object.
    /// </summary>
    /// <remarks>
    /// \boldsymbol{u}_{LevelSet}=\boldsymbol{u}_trans+\boldsymbol{\omega}\times\boldsymbol{r}
    /// </remarks>
    public class RigidObjectLevelSetVelocity : ParameterS, ILevelSetParameter {

        public RigidObjectLevelSetVelocity(string levelSetName, Particle[] Particles, double[] FluidViscosity, string[] FluidSpecies, Vector Gravity, double TimeStep, double GridToleranceParam) : base() {
            ParticleHydrodynamics = new ParticleHydrodynamics(2);
            ParticleHydrodynamics.SaveHydrodynamicOfPreviousTimestep(Particles);
            SpatialDimension = Particles[0].Motion.GetPosition().Dim;
            this.Particles = Particles;
            this.FluidSpecies = FluidSpecies;
            this.TimeStep = TimeStep;
            this.GridToleranceParam = GridToleranceParam;
            m_ParameterNames = VariableNames.AsLevelSetVariable(levelSetName, VariableNames.VelocityVector(SpatialDimension)).ToArray();
            this.FluidViscosity = FluidViscosity;
            this.Gravity = Gravity;
        }

        private readonly string[] FluidSpecies;
        private readonly int SpatialDimension;
        private readonly Particle[] Particles;
        private readonly ParticleHydrodynamics ParticleHydrodynamics;
        private readonly double TimeStep;
        private readonly double GridToleranceParam;
        private readonly string[] m_ParameterNames;
        private readonly double[] FluidViscosity;
        private readonly Vector Gravity;
        private int HistoryPosition = 0;

        public override IList<string> ParameterNames {
            get {
                return m_ParameterNames;
            }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override DelPartialParameterUpdate Update {
            get {
                return InternalParameterUpdate;
            }
        }

        void InternalParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                HistoryPosition = 0;
                XDGField Pressure = (XDGField)DomainVarFields[VariableNames.Pressure];
                XDGField[] Velocity = SpatialDimension.ForLoop(d => (XDGField)DomainVarFields[VariableNames.Velocity_d(d)]);
                ParticleHydrodynamics.SaveHydrodynamicOfPreviousIteration(Particles);
                ParticleHydrodynamicsIntegration hydrodynamicsIntegration = new ParticleHydrodynamicsIntegration(SpatialDimension, Velocity, Pressure, Particles[0].LsTrk, FluidViscosity);
                ParticleHydrodynamics.CalculateHydrodynamics(Particles, hydrodynamicsIntegration, FluidSpecies, Gravity, TimeStep);
                CalculateParticleVelocity(Particles, TimeStep);
                DGField[] levelSetVelocity = new ConventionalDGField[SpatialDimension];
                LevelSetTracker LsTrk = Particles[0].LsTrk;
                for (int d = 0; d < SpatialDimension; d++) {
                    levelSetVelocity[d] = ParameterVarFields[ParameterNames[d]];
                    levelSetVelocity[d].Clear();
                    levelSetVelocity[d].ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1); 
                        for (int j = 0; j < Len; j++) {
                            MultidimensionalArray globCoord = MultidimensionalArray.Create(K, SpatialDimension);
                            double[] cellCenter = LsTrk.GridDat.Cells.GetCenter(j + j0);
                            LsTrk.GridDat.TransformLocal2Global(NS, globCoord, j + j0);
                            for (int k = 0; k < K; k++) {
                                double velocityLocal = 0.0;
                                for (int p = 0; p < Particles.Length; p++) {
                                    Particle particle = Particles[p];
                                    if (particle.Contains(cellCenter, GridToleranceParam)) {
                                        double[] X = new double[] { globCoord[k, 0], globCoord[k, 1] };
                                        Vector radialVector = particle.CalculateRadialVector(X);
                                        switch (d) {
                                            case 0:
                                            velocityLocal = particle.Motion.GetTranslationalVelocity(HistoryPosition)[0] - particle.Motion.GetRotationalVelocity(HistoryPosition) * radialVector[1];
                                            break;
                                            case 1:
                                            velocityLocal = particle.Motion.GetTranslationalVelocity(HistoryPosition)[1] + particle.Motion.GetRotationalVelocity(HistoryPosition) * radialVector[0];
                                            break;
                                            default:
                                            throw new NotImplementedException("Rigid object solver only for 2D");
                                        }
                                    }
                                }
                                result[j, k] = velocityLocal;
                            }
                        }
                    });
                }
            }
        }

        private void CalculateParticleVelocity(Particle[] Particles, double dt) {
            using (new FuncTrace()) {
                foreach (Particle p in Particles) {
                    p.Motion.UpdateParticleVelocity(dt);
                }
            }
        }

        /// <summary>
        /// Velocity at rigid level set during level set parameter update
        /// </summary>
        /// <param name="levelSet"></param>
        /// <param name="time"></param>
        /// <param name="DomainVarFields"></param>
        /// <param name="ParameterVarFields"></param>
        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                foreach (Particle p in Particles) {
                    p.LsTrk = levelSet.Tracker;
                }
                ParticleHydrodynamics.SaveHydrodynamicOfPreviousTimestep(Particles);
                HistoryPosition = 1;
                InternalParameterUpdate(time, DomainVarFields, ParameterVarFields);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[SpatialDimension];
            var bv = DomainVarFields[VariableNames.VelocityX].Basis;
            var b = new Basis(bv.GridDat, Math.Max(6, bv.Degree));

            for (int d = 0; d < SpatialDimension; ++d) {
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(b, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }


    public class ActiveStress : ParameterS, ILevelSetParameter {

        public ActiveStress(string levelSetName, Particle[] Particles, double GridToleranceParam) : base() {
            D = Particles[0].Motion.GetPosition().Dim;
            this.Particles = Particles;
            this.GridToleranceParam = GridToleranceParam;
            m_ParameterNames = VariableNames.AsLevelSetVariable(levelSetName, VariableNames.SurfaceForceVector(D)).ToArray();
        }

        private readonly int D;
        private readonly Particle[] Particles;
        private readonly double GridToleranceParam;
        private readonly string[] m_ParameterNames;

        public override IList<string> ParameterNames {
            get {
                return m_ParameterNames;
            }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override DelPartialParameterUpdate Update {
            get {
                return InternalParameterUpdate;
            }
        }

        void InternalParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                DGField[] activeStress = new ConventionalDGField[D];
                LevelSetTracker LsTrk = Particles[0].LsTrk;
                for (int d = 0; d < D; d++) {
                    activeStress[d] = ParameterVarFields[ParameterNames[d]];
                    activeStress[d].Clear();
                    activeStress[d].ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1); // No nof Nodes
                        MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                        MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                        var Normals = LsTrk.DataHistories[1].Current.GetLevelSetNormals(NS, j0, Len);

                        for (int j = 0; j < Len; j++) {
                            MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                            double[] cellCenter = LsTrk.GridDat.Cells.GetCenter(j + j0);
                            LsTrk.GridDat.TransformLocal2Global(NS, globCoord, j+ j0);
                            for (int k = 0; k < K; k++) {
                                double activeStressLocal = 0.0;
                                for (int p = 0; p < Particles.Length; p++) {
                                    Particle particle = Particles[p];
                                    if (particle.Contains(cellCenter, GridToleranceParam)) {
                                        double activeStressMagnitude = particle.ActiveStress;
                                        double angle = particle.Motion.GetAngle(0);
                                        Vector orientation = new Vector(Math.Cos(angle), Math.Sin(angle));
                                        Vector orientationNormal = new Vector(-orientation[1], orientation[0]);
                                        Vector normalVector = new Vector(Normals[j, k, 0], Normals[j, k, 1]);
                                        if (orientation * normalVector <= 0)
                                            activeStressLocal = 0;
                                        else {
                                            //activeStressMagnitude *= Math.Pow(orientationNormal * normalVector, 2);
                                            //switch (d) {
                                            //    case 0:
                                            //    activeStressLocal = orientationNormal * normalVector > 0 ? -activeStressMagnitude * normalVector[1] : activeStressMagnitude * normalVector[1];
                                            //    break;
                                            //    case 1:
                                            //    activeStressLocal = orientationNormal * normalVector > 0 ? (activeStressMagnitude * normalVector[0]) : -activeStressMagnitude * normalVector[0];
                                            //    break;
                                            //    default:
                                            //    throw new NotImplementedException("Rigid object solver only for 2D");
                                            //}
                                            activeStressLocal = -activeStressMagnitude * orientation[d];
                                        }
                                    }
                                }
                                result[j, k] = activeStressLocal;
                            }
                        }
                    });
                }
            }
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            InternalParameterUpdate(time, DomainVarFields, ParameterVarFields);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var stresses = new (string, DGField)[D];
            var velocityBasis = DomainVarFields[VariableNames.VelocityX].Basis;
            var basis = new Basis(velocityBasis.GridDat, velocityBasis.Degree);

            for (int d = 0; d < D; ++d) {
                string paramName = ParameterNames[d];
                DGField stress = new SinglePhaseField(basis, paramName);
                stresses[d] = (paramName, stress);
            }
            return stresses;
        }
    }
}
