using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Linearization point for the Lax-Friedrichs flux
    /// </summary>
    /// <remarks>
    /// This is only for the ad-hoc linearization used by <see cref="XNSECommon.Operator.Convection.ConvectionInBulk_LLF"/>;
    /// the automatic differentiation of the Newton method will make this obsolete soon.
    /// </remarks>
    public class Velocity0 : ParameterS {
        int D;

        public Velocity0(int D) {
            this.D = D;
        }

        public override DelParameterFactory Factory => Velocity0Factory;

        public override IList<string> ParameterNames => BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D);

        /// <summary>
        /// Returns the <see cref="BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector"/>
        /// as a **shallow copy** of
        /// <see cref="BoSSS.Solution.NSECommon.VariableNames.VelocityVector"/> (i.e. the domain var);
        /// Thereby, the <see cref="Update"/> does not have to do anything.
        /// </summary>
        (string, DGField)[] Velocity0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocity0 = new (string, DGField)[D];
            for (int d = 0; d < D; ++d) {
                string velocityname = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d];
                DGField velocity = DomainVarFields[velocityname];
                string paramName = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d];
                velocity0[d] = (paramName, velocity);
            }
            return velocity0;
        }


        void InternalParameterUpdate(double t, IReadOnlyDictionary<string,DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            for (int d = 0; d < D; ++d) {
                string velocityname = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d];
                DGField velocity = DomainVarFields[velocityname];
                
                string paramName = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d];
                DGField velocity0 = ParameterVarFields[paramName];

                if(!object.ReferenceEquals(velocity, velocity0))
                    throw new ApplicationException("Messed up parameter system; Velocit0 is supposed to be a shallow copy of domain var Velocity"); 
            }
        }

        public override DelPartialParameterUpdate Update {
            get {
                return InternalParameterUpdate;
            }
        }
    }

    public class Velocity0Prescribed : ParameterS {
        int degree;

        IDictionary<string, Func<double[], double, double>> initial;

        string[] names;

        LevelSetTracker LsTrk;

        public override DelParameterFactory Factory => Velocity0PrescribedFactory;

        public override DelPartialParameterUpdate Update => ParameterUpdate;

        public Velocity0Prescribed(LevelSetTracker LsTrk, int d, int D, IDictionary<string, Func<double[], double, double>> initial, int degree) {
            this.degree = degree;
            this.LsTrk = LsTrk;


            names = new string[1];
            string velocity = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d];
            names[0] = velocity;
            this.initial = initial;
        }

        public static Velocity0Prescribed CreateFrom(LevelSetTracker LsTrk, int d, int D, AppControl control) {

            string velocity = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d];

            IDictionary<string, Func<double[], double, double>> initial = new Dictionary<string, Func<double[], double, double>>();
            foreach (string species in LsTrk.SpeciesNames) {
                string velocityOfSpecies = velocity + "#" + species;
                Func<double[], double, double> initialVelocity;
                if (control.InitialValues_Evaluators_TimeDep.TryGetValue(velocityOfSpecies, out Func<double[], double, double> initialValue_TimeDep)) {
                    initialVelocity = initialValue_TimeDep;
                } else if (control.InitialValues_Evaluators.TryGetValue(velocityOfSpecies, out Func<double[], double> initialValue)) {
                    initialVelocity = (X,t) => initialValue(X);
                } else {
                    initialVelocity = (X,t) => 0.0;
                }
                initial.Add(species, initialVelocity);
            }

            int velocityDegree;
            if (control.FieldOptions.TryGetValue(velocity, out FieldOpts opts)) {
                velocityDegree = opts.Degree;
            } else if (control.FieldOptions.TryGetValue("Velocity*", out FieldOpts velOpts)) {
                velocityDegree = velOpts.Degree;
            } else {
                velocityDegree = -1;
            }
            return new Velocity0Prescribed(LsTrk, d, D, initial, velocityDegree);
        }

        public override IList<string> ParameterNames => names;

        (string, DGField)[] Velocity0PrescribedFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            XDGBasis basis = new XDGBasis(LsTrk, degree != -1 ? degree : DomainVarFields.First().Value.Basis.Degree);
            XDGField velocity = new XDGField(basis, names[0]);

            foreach (var species in LsTrk.SpeciesNames)
                velocity.GetSpeciesShadowField(species).ProjectField(initial[species].Convert_Xt2X(0.0));

            return new (string, DGField)[] { (names[0], velocity) };
        }

        double timeOfLastUpdate = 0.0;
        public void ParameterUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if (timeOfLastUpdate != time) {
                foreach (var species in LsTrk.SpeciesNames) {
                    timeOfLastUpdate = time;
                    ((XDGField)ParameterVarFields[names[0]]).GetSpeciesShadowField(species).Clear();
                    ((XDGField)ParameterVarFields[names[0]]).GetSpeciesShadowField(species).ProjectField(X => initial[species](X, time));
                }
            }
        }


    }

    public class Velocity0MeanPrescribed : Velocity0Mean {
        public Velocity0MeanPrescribed(int D, LevelSetTracker LsTrk, int cutCellQuadOrder) :
            base(D, LsTrk, cutCellQuadOrder) { }

        protected override void Velocity0MeanUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            for (int d = 0; d < D; ++d) {
                foreach (string speciesName in SpeciesNames) {
                    XDGField paramMeanVelocity = (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]];
                    DGField speciesParam = paramMeanVelocity.GetSpeciesShadowField(speciesName);

                    XDGField velocity = (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d]];
                    DGField speciesVelocity = velocity.GetSpeciesShadowField(speciesName);

                    //Uncut
                    speciesParam.SetMeanValueTo(speciesVelocity);

                    //Cut
                    CellMask cutCells = regions.GetSpeciesMask(speciesName);
                    SpeciesId speciesId = speciesMap[speciesName];
                    CellQuadratureScheme scheme = schemeHelper.GetVolumeQuadScheme(speciesId, IntegrationDomain: cutCells);
                    SetMeanValueToMeanOf(speciesParam, speciesVelocity, minvol, cutCellQuadOrder, scheme);
                }
            }
        }

    }

    /// <summary>
    /// Cell-wise mean value, required for the for the localized Lax-Friedrichs flux <see cref="XNSECommon.Operator.Convection.ConvectionInBulk_LLF"/>,
    /// to have a constant Eigenvalue along an edge.
    /// </summary>
    public class Velocity0Mean : ParameterS, ILevelSetParameter {
        protected int D;

        protected int cutCellQuadOrder;

        protected LevelSetTracker LsTrk;

        public Velocity0Mean(int D, LevelSetTracker LsTrk, int cutCellQuadOrder) {
            this.D = D;
            this.cutCellQuadOrder = cutCellQuadOrder;
            this.LsTrk = LsTrk;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D);

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            this.LsTrk = ((XDGBasis)DomainVarFields.First().Value.Basis).Tracker;
            var velocity0Mean = new (string, DGField)[D];
            for (int d = 0; d < D; ++d) {
                XDGBasis U0meanBasis = new XDGBasis(LsTrk, 0);
                string paramName = BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d];
                XDGField U0mean = new XDGField(U0meanBasis, paramName);
                velocity0Mean[d] = (paramName, U0mean);
            }
            return velocity0Mean;
        }

        protected IList<string> SpeciesNames;

        protected LevelSetTracker.LevelSetRegions regions;

        protected IDictionary<string, SpeciesId> speciesMap;

        protected XQuadSchemeHelper schemeHelper;

        protected double minvol;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            LevelSetTracker lsTrkr = levelSet.Tracker;
            SpeciesNames = lsTrkr.SpeciesNames;
            regions = lsTrkr.Regions;
            IList<SpeciesId> speciesIds = lsTrkr.SpeciesIdS;
            schemeHelper = lsTrkr.GetXDGSpaceMetrics(speciesIds.ToArray(), cutCellQuadOrder).XQuadSchemeHelper;
            minvol = Math.Pow(lsTrkr.GridDat.Cells.h_minGlobal, D);

            speciesMap = new Dictionary<string, SpeciesId>(SpeciesNames.Count);
            foreach (string name in SpeciesNames) {
                speciesMap.Add(name, lsTrkr.GetSpeciesId(name));
            }
        }

        public override DelPartialParameterUpdate Update {
            get {
                return Velocity0MeanUpdate;
            }
        }

        protected virtual void Velocity0MeanUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(new FuncTrace()) {
                for(int d = 0; d < D; ++d) {
                    foreach(string speciesName in SpeciesNames) {
                        XDGField paramMeanVelocity = (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]];
                        DGField speciesParam = paramMeanVelocity.GetSpeciesShadowField(speciesName);

                        XDGField velocity = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d]];
                        DGField speciesVelocity = velocity.GetSpeciesShadowField(speciesName);

                        //Uncut
                        speciesParam.SetMeanValueTo(speciesVelocity);

                        //Cut
                        CellMask cutCells = regions.GetSpeciesMask(speciesName);
                        SpeciesId speciesId = speciesMap[speciesName];
                        CellQuadratureScheme scheme = schemeHelper.GetVolumeQuadScheme(speciesId, IntegrationDomain: cutCells);
                        SetMeanValueToMeanOf(speciesParam, speciesVelocity, minvol, cutCellQuadOrder, scheme);
                    }
                }
            }
        }

        protected static void SetMeanValueToMeanOf(DGField target, DGField source, double minvol, int order, CellQuadratureScheme scheme) {
            //Cut
            int D = source.GridDat.SpatialDimension;
            var rule = scheme.Compile(source.GridDat, order);
            CellQuadrature.GetQuadrature(new int[] { 2 },
                source.GridDat,
                rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.Clear();
                    source.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, 1);
                    Vol.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        int jCell = i + i0;

                        double Volume = ResultsOfIntegration[i, 1];
                        if (Math.Abs(Volume) < minvol * 1.0e-12) {
                            // keep current value
                            // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                        } else {
                            double IntVal = ResultsOfIntegration[i, 0];
                            target.SetMeanValue(jCell, IntVal / Volume);
                        }

                    }
                }).Execute();
        }
    }

    /// <summary>
    /// Gravity Parameter, note a few specials:
    /// 1.  When setting gravity in control file, set an acceleration. 
    ///     However in this class the acceleration is multiplied by the density of the respective phase, thus representing a gravitational volume force.
    /// 2.  The sign of this volume force is opposite to that of the acceleration. 
    ///     This is due to the implementation, where the volume force is brought to the left-hand-side of the Navier-Stokes-Equations <see cref="Solution.XNSECommon.Operator.MultiPhaseSource"/>
    /// </summary>
    public class Gravity : ParameterS {
        int degree;

        //Func<double[], double, double> initial;
        ScalarFunctionTimeDep initial;
        double rho;

        string[] names;

        public override DelParameterFactory Factory => ParameterFactory;

        public override DelPartialParameterUpdate Update => ParameterUpdate;

        public Gravity(string species, int d, int D, ScalarFunctionTimeDep initial, double rho, int degree) {
            this.degree = degree;

            names = new string[1];
            string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
            names[0] = gravity + "#" + species;
            this.initial = initial;
            this.rho = rho;
            //this.initial = (X, t) => -initial(X, t) * rho;
        }

        public static Gravity CreateFrom(string species, int d, int D, AppControl control, double rho, ScalarFunctionTimeDep gravityFunc) {
            string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
            string gravityOfSpecies = gravity + "#" + species;
            //Func<double[], double, double> initialGravity;
            //if (gravityFunc != null) {
            //    initialGravity = gravityFunc;
            //} else {
            //    initialGravity = (X, t) => 0.0;
            //}

            int gravityDegree;
            if (control.FieldOptions.TryGetValue(gravity, out FieldOpts opts)) {
                gravityDegree = Math.Max(0, opts.Degree);                
            } else if (control.FieldOptions.TryGetValue("Velocity*", out FieldOpts velOpts)) {
                gravityDegree = velOpts.Degree;
            } else {
                gravityDegree = 0;
            }

            return new Gravity(species, d, D, gravityFunc, rho, gravityDegree);
        }

        public override IList<string> ParameterNames => names;

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            DGField gravity = new SinglePhaseField(basis, names[0]);
            gravity.Clear();
            if(initial != null)
                gravity.ProjectField(-rho, initial.SetTime(0.0));
            //gravity.ProjectField(X => initial(X, 0.0));
            
            return new (string, DGField)[] { (names[0], gravity) };
        }

        public void ParameterUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if (initial != null) {
                DGField gravity = ParameterVarFields[names[0]];
                gravity.Clear();
                gravity.ProjectField(-rho, initial.SetTime(time));
            }
        }
    }

    public class VolumeForce : ParameterS {
        int degree;

        ScalarFunctionTimeDep initial;

        string[] names;

        public override DelParameterFactory Factory => ParameterFactory;

        public override DelPartialParameterUpdate Update => ParameterUpdate;

        public VolumeForce(string species, int d, int D, ScalarFunctionTimeDep initial, int degree) {
            this.degree = degree;

            names = new string[1];
            string volforce = BoSSS.Solution.NSECommon.VariableNames.VolumeForceVector(D)[d];
            names[0] = volforce + "#" + species;
            this.initial = initial;
        }

        public static VolumeForce CreateFrom(string species, int d, int D, AppControl control, ScalarFunctionTimeDep VolForceFunc) {
            string gravity = BoSSS.Solution.NSECommon.VariableNames.VolumeForce_d(d);
            string gravityOfSpecies = gravity + "#" + species;


            int volForceDegree;
            if (control.FieldOptions.TryGetValue(gravityOfSpecies, out FieldOpts opts)) {
                volForceDegree = Math.Max(0, opts.Degree);
            } else if (control.FieldOptions.TryGetValue("Velocity*", out FieldOpts velOpts)) {
                volForceDegree = velOpts.Degree;
            } else {
                volForceDegree = 0;
            }
            return new VolumeForce(species, d, D, VolForceFunc, volForceDegree);
        }

        public override IList<string> ParameterNames => names;

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            DGField volforce = new SinglePhaseField(basis, names[0]);
            if (initial != null)
                volforce.ProjectField(-1, initial.SetTime(0.0));
            return new (string, DGField)[] { (names[0], volforce) };
        }

        public void ParameterUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if (initial != null) {
                DGField volforce = ParameterVarFields[names[0]];
                volforce.Clear();
                volforce.ProjectField(-1, initial.SetTime(time));
            }
        }
    }


    /// <summary>
    /// Computation of normals from the level-set;
    /// - computed normals are typically **not of unit length**, i.e. the vectors must be normalized before use!
    /// - computed from broken derivatives, i.e. un-filtered
    /// </summary>
    public class Normals : ParameterS, ILevelSetParameter {

        int D;
        int degree;
        IList<string> parameterNames;

        public Normals(int D, int degree) {
            this.D = D;
            this.degree = degree;
            parameterNames = BoSSS.Solution.NSECommon.VariableNames.NormalVector(D);
        }

        public Normals(int D, int degree, string levelSetName) : this(D, degree){
            parameterNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable( 
                levelSetName, BoSSS.Solution.NSECommon.VariableNames.NormalVector(D));
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => parameterNames;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            LevelSet Phi = levelSet.CGLevelSet;
            DGField[] Normals = new SinglePhaseField[D];
            for (int i = 0; i < D; ++i) {
                Normals[i] = ParameterVarFields[parameterNames[i]];
            }
            VectorField<DGField> normalVector = new VectorField<DGField>(Normals);
            normalVector.Clear();
            normalVector.Gradient(1.0, Phi);
        }

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis basis = new Basis(gridData, degree);
            VectorField<SinglePhaseField> Normals = new VectorField<SinglePhaseField>(D, basis, parameterNames[0], SinglePhaseField.Factory);

            (string, DGField)[] normals = new (string, DGField)[D];
            for (int d = 0; d < D; ++d) {
                normals[d] = (parameterNames[d], Normals[d]);
            }
            return normals;
        }
    }

    public class BeltramiGradient : ILevelSetParameter {
        DoNotTouchParameters AdvancedDiscretizationOptions;

        string[] parameters;

        int degree;

        public BeltramiGradient(int D, DoNotTouchParameters AdvancedDiscretizationOptions, int degree) {
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D);
            this.AdvancedDiscretizationOptions = AdvancedDiscretizationOptions;
            this.degree = degree;
        }

        public IList<string> ParameterNames => parameters;

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            int paramCount = parameters.Length;
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[parameters.Length];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis basis = new Basis(gridData, degree);
            for (int i = 0; i < paramCount; ++i) {
                fields[i] = (parameters[i], new SinglePhaseField(basis, parameters[i]));
            }
            return fields;
        }

        public void LevelSetParameterUpdate(
           DualLevelSet phaseInterface,
           double time,
           IReadOnlyDictionary<string, DGField> DomainVarFields,
           IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            //dgLevSetGradient and update curvature
            VectorField<SinglePhaseField> filtLevSetGradient;
            CurvatureAlgorithms.LaplaceBeltramiDriver(
                AdvancedDiscretizationOptions.SST_isotropicMode,
                AdvancedDiscretizationOptions.FilterConfiguration,
                out filtLevSetGradient,
                phaseInterface.Tracker,
                phaseInterface.DGLevelSet);
            for (int i = 0; i < parameters.Length; ++i) {
                ParameterVarFields[parameters[i]].Clear();
                ParameterVarFields[parameters[i]].Acc(1.0, filtLevSetGradient[i]);
            }
        }
    }

    public class BeltramiGradientAndCurvature : ParameterS, ILevelSetParameter {
        DoNotTouchParameters AdvancedDiscretizationOptions;

        int m_HMForder;

        string[] lsParameters;

        int gradientDegree;

        int curvatureDegree;

        public BeltramiGradientAndCurvature(int curvatureDegree, int gradientDegree, int m_HMForder, DoNotTouchParameters AdvancedDiscretizationOptions, int D) {
            lsParameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            this.AdvancedDiscretizationOptions = AdvancedDiscretizationOptions;
            this.m_HMForder = m_HMForder;
            this.gradientDegree = gradientDegree;
            this.curvatureDegree = curvatureDegree;
        }

        IList<string> ILevelSetParameter.ParameterNames => lsParameters;

        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };

        public override DelParameterFactory Factory => CurvatureFactory;

        (string ParameterName, DGField ParamField)[] CurvatureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            //Curvature
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[1];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            string curvatureName = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            fields[0] = (curvatureName, new SinglePhaseField(curvatureBasis, curvatureName));
            return fields;
        }

        (string ParameterName, DGField ParamField)[] ILevelSetParameter.ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            int paramCount = lsParameters.Length;
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[lsParameters.Length];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            
            Basis basis = new Basis(gridData, gradientDegree);
            for (int i = 0; i < paramCount - 1; ++i) {
                fields[i] = (lsParameters[i], new SinglePhaseField(basis, lsParameters[i]));
            }
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            fields[2] = (lsParameters[2], new SinglePhaseField(curvatureBasis, lsParameters[2]));
            return fields;
        }

        public void LevelSetParameterUpdate(
           DualLevelSet phaseInterface,
           double time,
           IReadOnlyDictionary<string, DGField> DomainVarFields,
           IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];
            VectorField<SinglePhaseField> filtLevSetGradient;
            CurvatureAlgorithms.CurvatureDriver(
                AdvancedDiscretizationOptions.SST_isotropicMode,
                AdvancedDiscretizationOptions.FilterConfiguration,
                Curvature,
                out filtLevSetGradient,
                phaseInterface.Tracker,
                m_HMForder,
                phaseInterface.DGLevelSet);
            for (int i = 0; i < lsParameters.Length - 1; ++i) {
                ParameterVarFields[lsParameters[i]].Clear();
                if(filtLevSetGradient != null)
                    ParameterVarFields[lsParameters[i]].Acc(1.0, filtLevSetGradient[i]);
            }
        }
    }

    public class MaxSigma : ParameterS, ILevelSetParameter {
        PhysicalParameters physParams;

        DoNotTouchParameters dntParams;

        int cutCellQuadOrder;

        double dt;

        public MaxSigma(PhysicalParameters physParams, DoNotTouchParameters dntParams, int cutCellQuadOrder, double dt) {
            this.physParams = physParams;
            this.dntParams = dntParams;
            this.cutCellQuadOrder = cutCellQuadOrder;
            this.dt = dt;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        string[] name = new string[] { BoSSS.Solution.NSECommon.VariableNames.MaxSigma };

        public override IList<string> ParameterNames => name;

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            string name = BoSSS.Solution.NSECommon.VariableNames.MaxSigma;
            IGridData grid = DomainVarFields.FirstOrDefault().Value.GridDat;
            Basis constant = new Basis(grid, 0);
            SinglePhaseField sigmaField = new SinglePhaseField(constant, name);
            return new (string ParameterName, DGField ParamField)[] { (name, sigmaField) };
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            DGField sigmaMax = ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.MaxSigma];
            LevelSetTracker lsTrkr = levelSet.Tracker;

            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengths =
                lsTrkr.GetXDGSpaceMetrics(lsTrkr.SpeciesIdS.ToArray(), cutCellQuadOrder).CutCellMetrics.InterfaceArea;

            foreach (Chunk cnk in lsTrkr.Regions.GetCutCellMask()) {
                for (int i = cnk.i0; i < cnk.JE; i++) {
                    double ILen = InterfaceLengths.ElementAt(0).Value[i];
                    //ILen /= LevSet_Deg;
                    double sigmaILen_Max = (this.physParams.rho_A + this.physParams.rho_B)
                           * Math.Pow(ILen, 3) / (2 * Math.PI * dt.Pow2());

                    if (dntParams.SetSurfaceTensionMaxValue && (physParams.Sigma > sigmaILen_Max)) {
                        sigmaMax.SetMeanValue(i, sigmaILen_Max * 0.5);
                        //Console.WriteLine("set new sigma value: {0}; {1}", sigmaILen_Max, sigmaILen_Max/physParams.Sigma);
                    } else {
                        sigmaMax.SetMeanValue(i, this.physParams.Sigma * 0.5);
                    }
                }
            }
        }

    }
}
