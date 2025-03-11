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
        private int CurrentDimension = 0;
        double OldTime = -1;

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

        private double VelocityFunction(double[] X, double t) {
            double VelocityFunction = 0;
            for (int p = 0; p < Particles.Length; p++) {
                Particle particle = Particles[p];
                if (particle.Contains(X, 2 * GridToleranceParam)) {
                    Vector radialVector = particle.CalculateRadialVector(X);
                    VelocityFunction = CurrentDimension switch {
                        0 => particle.Motion.GetTranslationalVelocity(0)[0] - particle.Motion.GetRotationalVelocity(0) * radialVector[1],
                        1 => particle.Motion.GetTranslationalVelocity(0)[1] + particle.Motion.GetRotationalVelocity(0) * radialVector[0],
                        _ => throw new NotImplementedException("Rigid object solver only for 2D"),
                    };
                }
            }
            return VelocityFunction;
        }

        private void InternalParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                XDGField Pressure = (XDGField)DomainVarFields[VariableNames.Pressure];
                XDGField[] Velocity = SpatialDimension.ForLoop(d => (XDGField)DomainVarFields[VariableNames.Velocity_d(d)]);
                LevelSetTracker LsTrk = Particles[0].LsTrk;//???? Besserer Weg?
                CellMask AllCutCells = Particles[0].LsTrk.Regions.GetCutCellMask();

                ParticleHydrodynamics.SaveHydrodynamicOfPreviousIteration(Particles);
                ParticleHydrodynamicsIntegration hydrodynamicsIntegration = new(SpatialDimension, Velocity, Pressure, LsTrk, FluidViscosity);
                ParticleHydrodynamics.CalculateHydrodynamics(Particles, AllCutCells, hydrodynamicsIntegration, FluidSpecies, Gravity, TimeStep);
                UpdateVelocity();
                ProjectFieldVel(t, ParameterVarFields, AllCutCells);
            }
        }

        private void ProjectFieldVel(double t, IReadOnlyDictionary<string, DGField> ParameterVarFields, CellMask AllCutCells) {
            DGField[] levelSetVelocity = new ConventionalDGField[SpatialDimension];
            for (int d = 0; d < SpatialDimension; d++) {
                CurrentDimension = d;
                ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(VelocityFunction, t);
                levelSetVelocity[d] = ParameterVarFields[ParameterNames[d]];
                levelSetVelocity[d].Clear();
                levelSetVelocity[d].ProjectField(1.0, Function, new CellQuadratureScheme(true, AllCutCells));
            }
        }

        private void UpdateVelocity() {
            foreach (Particle p in Particles) {
                p.Motion.UpdateParticleVelocity(TimeStep);
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
                if (time > OldTime) {
                    foreach (Particle p in Particles) {
                        p.LsTrk = levelSet.Tracker;
                    }
                    ParticleHydrodynamics.SaveHydrodynamicOfPreviousTimestep(Particles);
                    OldTime = time;
                }
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[SpatialDimension];
            var bv = DomainVarFields[VariableNames.VelocityX].Basis;
            var b = new Basis(bv.GridDat, bv.Degree);

            for (int d = 0; d < SpatialDimension; ++d) {
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(b, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }


    /// <summary>
    /// Sets orientation vector of particles as parameter DGField
    /// </summary>
    public class Orientation : ParameterS, ILevelSetParameter {

        private readonly IList<string> parameterNames;
        private readonly Particle[] Particles;
        private readonly double Tolerance;
        private int CurrentDimension = 0;
        private readonly int D;
        private double OldTime = -1;

        public Orientation(string levelSetName, Particle[] Particles, double GridToleranceParam) : base() {
            D = Particles[0].Motion.GetPosition(0).Count;
            parameterNames = VariableNames.AsLevelSetVariable(levelSetName, VariableNames.OrientationVector(D));
            this.Particles = Particles;
            this.Tolerance = GridToleranceParam;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => parameterNames;

        double OrientationFunction(double[] X, double t) {
            double OrientationFunction = 0;
            for (int p = 0; p < Particles.Length; p++) {
                if (Particles[p].Contains(X, 2 * Tolerance) && Particles[p].ActiveStress != 0) {
                    OrientationFunction = CurrentDimension switch {
                        0 => Math.Cos(Particles[p].Motion.GetAngle(1)),
                        1 => Math.Sin(Particles[p].Motion.GetAngle(1)),
                        _ => throw new NotImplementedException("Rigid object solver only for 2D"),
                    };
                }
            }
            return OrientationFunction;
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if (time > OldTime) {
                DGField[] Orientation = new SinglePhaseField[D];
                for (int d = 0; d < D; ++d) {
                    Orientation[d] = ParameterVarFields[parameterNames[d]];
                    CurrentDimension = d;
                    ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(OrientationFunction, time);
                    Orientation[d].Clear();
                    Orientation[d].ProjectField(1.0, Function, new CellQuadratureScheme(true, Particles[0].LsTrk.Regions.GetCutCellMask()));
                }
                OldTime = time;
            } 
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var orientation = new (string, DGField)[D];
            var velocityBasis = DomainVarFields[VariableNames.VelocityX].Basis;
            var basis = new Basis(velocityBasis.GridDat, velocityBasis.Degree);
            for (int d = 0; d < D; ++d) {
                string paramName = ParameterNames[d];
                DGField orientationElement = new SinglePhaseField(basis, paramName);
                orientation[d] = (paramName, orientationElement);
            }
            return orientation;
        }
    }

    public class PhoreticActivity : ParameterS, ILevelSetParameter {

        public PhoreticActivity(string levelSetName, Particle[] Particles, double GridToleranceParam) : base() {
            D = Particles[0].Motion.GetPosition().Dim;
            this.Particles = Particles;
            this.GridToleranceParam = GridToleranceParam;
            m_ParameterNames = new string[1];
            m_ParameterNames[0] = VariableNames.AsLevelSetVariable(levelSetName, VariableNames.Phoretic);
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
                DGField[] activeStress = new ConventionalDGField[1];
                LevelSetTracker LsTrk = Particles[0].LsTrk;
                activeStress[0] = ParameterVarFields[ParameterNames[0]];
                activeStress[0].Clear();
                activeStress[0].ProjectField(1.0,
                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes

                        var Normals = LsTrk.DataHistories[1].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {
                        MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                        double[] cellCenter = LsTrk.GridDat.Cells.GetCenter(j + j0);
                        LsTrk.GridDat.TransformLocal2Global(NS, globCoord, j + j0);
                        for (int k = 0; k < K; k++) {
                            double activeStressLocal = 0.0;
                            for (int p = 0; p < Particles.Length; p++) {
                                Particle particle = Particles[p];
                                if (particle.Contains(cellCenter, GridToleranceParam)) {
                                    double activeStressMagnitude = particle.phoreticActivity;
                                    double angle = particle.Motion.GetAngle(0);
                                    Vector orientation = new Vector(Math.Cos(angle), Math.Sin(angle));
                                    Vector orientationNormal = new Vector(-orientation[1], orientation[0]);
                                    Vector normalVector = new Vector(Normals[j, k, 0], Normals[j, k, 1]);
                                    if (orientation * normalVector <= 0)
                                        activeStressLocal = 0;
                                    else {
                                        activeStressLocal = activeStressMagnitude;
                                    }
                                }
                            }
                            result[j, k] = activeStressLocal;
                        }
                    }
                });
            }
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            InternalParameterUpdate(time, DomainVarFields, ParameterVarFields);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var phoreticBoundary = new (string, DGField)[1];
            var velocityBasis = DomainVarFields[VariableNames.VelocityX].Basis;
            var basis = new Basis(velocityBasis.GridDat, 12);
            string paramName = ParameterNames[0];
            DGField stress = new SinglePhaseField(basis, paramName);
            phoreticBoundary[0] = (paramName, stress);
            return phoreticBoundary;
        }
    }
}
