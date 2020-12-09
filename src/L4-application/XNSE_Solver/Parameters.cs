using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.XSolver;

namespace BoSSS.Application.XNSE_Solver
{
    class Velocity0 : Parameter
    {
        int D;

        public Velocity0(int D)
        {
            this.D = D;
        }

        public override DelParameterFactory Factory => Velocity0Factory;

        public override IList<string> ParameterNames => BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D);

        (string, DGField)[] Velocity0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            var velocity0 = new (string, DGField)[D];
            for(int d = 0; d < D; ++d)
            {
                string velocityname = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d];
                DGField velocity = DomainVarFields[velocityname];
                string paramName = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d];
                velocity0[d] = (paramName, velocity);
            }
            return velocity0;
        }
    }

    class Velocity0Prescribed : Parameter {
        int degree;

        IDictionary<string, Func<double[], double>> initial;

        string[] names;

        LevelSetTracker LsTrk;

        public override DelParameterFactory Factory => Velocity0PrescribedFactory;

        public Velocity0Prescribed(LevelSetTracker LsTrk, int d,int D, IDictionary<string, Func<double[], double>> initial, int degree) {
            this.degree = degree;
            this.LsTrk = LsTrk;


            names = new string[1];
            string velocity = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d];
            names[0] = velocity;
            this.initial = initial;
        }

        public static Velocity0Prescribed CreateFrom(LevelSetTracker LsTrk, int d, int D, AppControl control) {
            
            string velocity = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d];

            IDictionary<string, Func<double[], double>> initial = new Dictionary<string, Func<double[], double>>();
            foreach (string species in LsTrk.SpeciesNames) {
                string velocityOfSpecies = velocity + "#" + species;
                Func<double[], double> initialVelocity;
                if (control.InitialValues_Evaluators.TryGetValue(velocityOfSpecies, out Func<double[], double> initialValue)) {
                    initialVelocity = initialValue;
                } else {
                    initialVelocity = X => 0.0;
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

            foreach(var species in LsTrk.SpeciesNames)
                velocity.GetSpeciesShadowField(species).ProjectField(initial[species]);

            return new (string, DGField)[] { (names[0], velocity) };
        } 


    }

    class Velocity0MeanPrescribed : Velocity0Mean {
        public Velocity0MeanPrescribed(int D, LevelSetTracker LsTrk, int cutCellQuadOrder) :
            base(D, LsTrk, cutCellQuadOrder) { }        

        protected override void Velocity0MeanUpdate(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
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

    class Velocity0Mean : Parameter, ILevelSetParameter
    {
        protected int D;

        protected int cutCellQuadOrder;

        protected LevelSetTracker LsTrk;

        public Velocity0Mean(int D, LevelSetTracker LsTrk, int cutCellQuadOrder)
        {
            this.D = D;
            Update = Velocity0MeanUpdate;
            this.cutCellQuadOrder = cutCellQuadOrder;
            this.LsTrk = LsTrk;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D);

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            var velocity0Mean = new (string, DGField)[D];
            for(int d = 0; d < D; ++d)
            {
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
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            LevelSetTracker lsTrkr = levelSet.Tracker;
            SpeciesNames = lsTrkr.SpeciesNames;
            regions = lsTrkr.Regions;
            IList<SpeciesId> speciesIds = lsTrkr.SpeciesIdS;
            schemeHelper = lsTrkr.GetXDGSpaceMetrics(speciesIds.ToArray(), cutCellQuadOrder).XQuadSchemeHelper;
            minvol = Math.Pow(lsTrkr.GridDat.Cells.h_minGlobal, D);

            speciesMap = new Dictionary<string, SpeciesId>(SpeciesNames.Count);
            foreach (string name in SpeciesNames)
            {
                speciesMap.Add(name, lsTrkr.GetSpeciesId(name));
            }
        }

        protected virtual void Velocity0MeanUpdate(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            for(int d = 0; d < D; ++d)
            {
                foreach(string speciesName in SpeciesNames)
                {
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

        protected static void SetMeanValueToMeanOf(DGField target, DGField source, double minvol, int order, CellQuadratureScheme scheme)
        {
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
                    for (int i = 0; i < Length; i++)
                    {
                        int jCell = i + i0;

                        double Volume = ResultsOfIntegration[i, 1];
                        if (Math.Abs(Volume) < minvol * 1.0e-12)
                        {
                            // keep current value
                            // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                        }
                        else
                        {
                            double IntVal = ResultsOfIntegration[i, 0];
                            target.SetMeanValue( jCell, IntVal / Volume);
                        }

                    }
                }).Execute();
        }
    }

    class Gravity : Parameter
    {
        int degree;

        Func<double[], double> initial;

        string[] names;

        public override DelParameterFactory Factory => GravityFactory;

        public Gravity(string species, int d, int D, Func<double[], double> initial, double rho, int degree)
        {
            this.degree = degree;

            names = new string[1];
            string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
            names[0] = gravity + "#" + species;
            this.initial = X => -initial(X) * rho;
        }

        public static Gravity CreateFrom(string species, int d, int D, AppControl control, double rho)
        {
            string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
            string gravityOfSpecies = gravity + "#" + species;
            Func<double[], double> initialGravity;
            if (control.InitialValues_Evaluators.TryGetValue(gravityOfSpecies, out Func<double[], double> initialValue))
            {
                initialGravity = initialValue;
            }
            else
            {
                initialGravity = X => 0.0;
            }

            int gravityDegree;
            if (control.FieldOptions.TryGetValue(gravityOfSpecies, out FieldOpts opts))
            {
                gravityDegree = opts.Degree;
            }
            else if(control.FieldOptions.TryGetValue("Velocity*", out FieldOpts velOpts))
            {
                gravityDegree = velOpts.Degree;
            }
            else
            {
                gravityDegree = 0;
            }
            return new Gravity(species, d, D, initialGravity, rho, gravityDegree);
        }

        public override IList<string> ParameterNames => names;

        (string, DGField)[] GravityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            DGField gravity = new SinglePhaseField(basis, names[0]);
            gravity.ProjectField(initial);
            return new (string, DGField)[] { (names[0], gravity) };
        }
    }

    class Normals : Parameter, ILevelSetParameter
    {
        int D;

        int degree;

        public Normals(int D, int degree)
        {
            this.D = D;
            this.degree = degree;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => BoSSS.Solution.NSECommon.VariableNames.NormalVector(D);

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            LevelSet Phi = levelSet.DGLevelSet;
            DGField[] Normals = new SinglePhaseField[D];
            for (int i = 0; i < D; ++i)
            {
                Normals[i] = ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[i]];
            }
            VectorField<DGField> normalVector = new VectorField<DGField>(Normals);
            Normals.Clear();
            normalVector.Gradient(1.0, Phi);
        }

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis basis = new Basis(gridData, degree);
            VectorField<SinglePhaseField> Normals = new VectorField<SinglePhaseField>(D, basis, SinglePhaseField.Factory);

            (string, DGField)[] normals = new (string, DGField)[D];
            for(int d = 0; d <D; ++d)
            {
                normals[d] = (BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[d], Normals[d] );
            }
            return normals;
        }
    }

    class BeltramiGradient : ILevelSetParameter
    {
        DoNotTouchParameters AdvancedDiscretizationOptions;

        string[] parameters;

        int degree;

        public BeltramiGradient(int D, DoNotTouchParameters AdvancedDiscretizationOptions, int degree)
        {
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D);
            this.AdvancedDiscretizationOptions = AdvancedDiscretizationOptions;
            this.degree = degree;
        }

        public IList<string> ParameterNames => parameters;

        public static BeltramiGradient CreateFrom(XNSE_Control control, string levelSetName, int D)
        {
            string levelSet = levelSetName;
            int levelSetDegree;
            if (control.FieldOptions.TryGetValue(levelSet, out FieldOpts lsOpts))
            {
                var levelSetSource = control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource;
                levelSetDegree = (levelSetSource == CurvatureAlgorithms.LevelSetSource.fromDG) ? lsOpts.Degree : lsOpts.Degree + 1;
            }
            else
            {
                levelSetDegree = 1;
            }

            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new BeltramiGradient(D, AdvancedDiscretizationOptions, levelSetDegree);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            int paramCount = parameters.Length;
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[parameters.Length];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis basis = new Basis(gridData, degree);
            for (int i = 0; i < paramCount; ++i)
            {
                fields[i] = (parameters[i], new SinglePhaseField(basis, parameters[i]));
            }
            return fields;
        }

        public void LevelSetParameterUpdate(
           DualLevelSet phaseInterface,
           double time,
           IReadOnlyDictionary<string, DGField> DomainVarFields,
           IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            //dgLevSetGradient and update curvature
            VectorField<SinglePhaseField> filtLevSetGradient;
            CurvatureAlgorithms.LaplaceBeltramiDriver(
                AdvancedDiscretizationOptions.SST_isotropicMode,
                AdvancedDiscretizationOptions.FilterConfiguration,
                out filtLevSetGradient,
                phaseInterface.Tracker,
                phaseInterface.DGLevelSet);
            for(int i = 0; i < parameters.Length; ++i)
            {
                ParameterVarFields[parameters[i]].Clear();
                ParameterVarFields[parameters[i]].Acc(1.0, filtLevSetGradient[i]);
            }
        }
    }

    class BeltramiGradientAndCurvature : Parameter, ILevelSetParameter
    {
        DoNotTouchParameters AdvancedDiscretizationOptions;

        int m_HMForder;

        string[] lsParameters;

        int gradientDegree;

        int curvatureDegree;

        public BeltramiGradientAndCurvature(int curvatureDegree, int gradientDegree, int m_HMForder, DoNotTouchParameters AdvancedDiscretizationOptions, int D)
        {
            lsParameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            this.AdvancedDiscretizationOptions = AdvancedDiscretizationOptions;
            this.m_HMForder = m_HMForder;
            this.gradientDegree = gradientDegree;
            this.curvatureDegree = curvatureDegree;
        }

        IList<string> ILevelSetParameter.ParameterNames => lsParameters;

        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };

        public override DelParameterFactory Factory => CurvatureFactory;

        (string ParameterName, DGField ParamField)[] CurvatureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            //Curvature
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[1];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            string curvatureName = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            fields[0] = (curvatureName, new SinglePhaseField(curvatureBasis, curvatureName));
            return fields;
        }

        public static BeltramiGradientAndCurvature CreateFrom(XNSE_Control control, string levelSetName, int m_HMForder, int D)
        {
            string curvature = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            int curvatureDegree;
            if (control.FieldOptions.TryGetValue(curvature, out FieldOpts opts))
            {
                curvatureDegree = opts.Degree;
            }
            else
            {
                curvatureDegree = 1;
            }
            string levelSet = levelSetName;
            int levelSetDegree;
            if (control.FieldOptions.TryGetValue(levelSet, out FieldOpts lsOpts))
            {
                var levelSetSource = control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource;
                levelSetDegree = (levelSetSource == CurvatureAlgorithms.LevelSetSource.fromDG) ? lsOpts.Degree : lsOpts.Degree + 1;
            }
            else
            {
                levelSetDegree = 1;
            }
            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new BeltramiGradientAndCurvature(curvatureDegree, levelSetDegree, m_HMForder, AdvancedDiscretizationOptions, D);
        }

        (string ParameterName, DGField ParamField)[] ILevelSetParameter.ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            int paramCount = lsParameters.Length;
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[lsParameters.Length];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            //Basis basis = new Basis(gridData, gradientDegree);
            for (int i = 0; i < paramCount; ++i)
            {
                if (i == 2) {
                    Basis basis = new Basis(gridData, curvatureDegree);
                    fields[i] = (lsParameters[i], new SinglePhaseField(basis, lsParameters[i]));
                } else {
                    Basis basis = new Basis(gridData, gradientDegree);
                    fields[i] = (lsParameters[i], new SinglePhaseField(basis, lsParameters[i]));
                }
            }
            return fields;
        }

        public void LevelSetParameterUpdate(
           DualLevelSet phaseInterface,
           double time,
           IReadOnlyDictionary<string, DGField> DomainVarFields,
           IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
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
            for (int i = 0; i < lsParameters.Length - 1; ++i)
            {
                ParameterVarFields[lsParameters[i]].Clear();
                ParameterVarFields[lsParameters[i]].Acc(1.0, filtLevSetGradient[i]);
            }
        }
    }

    class MaxSigma : Parameter, ILevelSetParameter
    {
        PhysicalParameters physParams;

        DoNotTouchParameters dntParams;

        int cutCellQuadOrder;

        double dt;

        public MaxSigma(PhysicalParameters physParams, DoNotTouchParameters dntParams, int cutCellQuadOrder, double dt)
        {
            this.physParams = physParams;
            this.dntParams = dntParams;
            this.cutCellQuadOrder = cutCellQuadOrder;
            this.dt = dt;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        string[] name = new string[] { BoSSS.Solution.NSECommon.VariableNames.MaxSigma };

        public override IList<string> ParameterNames => name;

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            string name = BoSSS.Solution.NSECommon.VariableNames.MaxSigma;
            IGridData grid = DomainVarFields.FirstOrDefault().Value.GridDat;
            Basis constant = new Basis(grid, 0);
            SinglePhaseField sigmaField = new SinglePhaseField(constant, name);
            return new (string ParameterName, DGField ParamField)[] { (name, sigmaField) };
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, 
            IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            DGField sigmaMax = ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.MaxSigma];
            LevelSetTracker lsTrkr = levelSet.Tracker;

            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengths =
                lsTrkr.GetXDGSpaceMetrics(lsTrkr.SpeciesIdS.ToArray(), cutCellQuadOrder).CutCellMetrics.InterfaceArea;

            foreach (Chunk cnk in lsTrkr.Regions.GetCutCellMask())
            {
                for (int i = cnk.i0; i < cnk.JE; i++)
                {
                    double ILen = InterfaceLengths.ElementAt(0).Value[i];
                    //ILen /= LevSet_Deg;
                    double sigmaILen_Max = (this.physParams.rho_A + this.physParams.rho_B)
                           * Math.Pow(ILen, 3) / (2 * Math.PI * dt.Pow2());

                    if (dntParams.SetSurfaceTensionMaxValue && (physParams.Sigma > sigmaILen_Max))
                    {
                        sigmaMax.SetMeanValue(i, sigmaILen_Max * 0.5);
                        //Console.WriteLine("set new sigma value: {0}; {1}", sigmaILen_Max, sigmaILen_Max/physParams.Sigma);
                    }
                    else
                    {
                        sigmaMax.SetMeanValue(i, this.physParams.Sigma * 0.5);
                    }
                }
            }
        }

    }

    class LevelSetVelocity : ILevelSetParameter
    {
        protected int D;

        protected IList<string> parameters;

        protected int degree;

        public IList<string> ParameterNames => parameters;

        protected XNSE_Control.InterfaceVelocityAveraging averagingMode;

        protected PhysicalParameters physicalParameters;

        public LevelSetVelocity(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters)
        {
            this.averagingMode = averagingMode;
            this.physicalParameters = physicalParameters;
            this.D = D;
            this.degree = degree;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
        }

        public virtual void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            //Mean Velocity
            XDGField[] EvoVelocity; // = new XDGField[]
            try {
                EvoVelocity = new XDGField[]
                {
                    (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX],
                    (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityY],
                };
            } catch {
                Console.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                EvoVelocity = new XDGField[]
                {
                    (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0X],
                    (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0Y],
                };
            }


            int D = EvoVelocity.Length;
            DGField[] meanVelocity;

            meanVelocity = new ConventionalDGField[D];

            double rho_A = physicalParameters.rho_A, rho_B = physicalParameters.rho_B;
            double mu_A = physicalParameters.mu_A, mu_B = physicalParameters.mu_B;
            LevelSetTracker lsTrkr = levelSet.Tracker;
            CellMask CC = lsTrkr.Regions.GetCutCellMask4LevSet(0);
            CellMask Neg = lsTrkr.Regions.GetLevelSetWing(0, -1).VolumeMask;
            CellMask Pos = lsTrkr.Regions.GetLevelSetWing(0, +1).VolumeMask;
            CellMask posNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
            CellMask negNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

            for (int d = 0; d < D; d++)
            {
                Basis b = EvoVelocity[d].Basis.NonX_Basis;
                meanVelocity[d] = ParameterVarFields[ParameterNames[d]];


                foreach (string spc in lsTrkr.SpeciesNames)
                {
                    double rhoSpc;
                    double muSpc;
                    switch (spc)
                    {
                        case "A": rhoSpc = rho_A; muSpc = mu_A; break;
                        case "B": rhoSpc = rho_B; muSpc = mu_B; break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }

                    double scale = 1.0;
                    switch (averagingMode)
                    {
                        case XNSE_Control.InterfaceVelocityAveraging.mean:
                            {
                                scale = 0.5;
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.density:
                            {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.viscosity:
                            {
                                scale = muSpc / (mu_A + mu_B);
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.phaseA:
                            {
                                scale = (spc == "A") ? 1.0 : 0.0;
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.phaseB:
                            {
                                scale = (spc == "B") ? 1.0 : 0.0;
                                break;
                            }
                    }

                    meanVelocity[d].Acc(scale, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), CC);
                    switch (spc)
                    {
                        //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                        case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), negNear); break;
                        case "B": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), posNear); break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }
                }
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            var velocties = new (string, DGField)[D];
            for (int d = 0; d < D; ++d)
            {
                Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(basis, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }

    class LevelSetVelocityEvaporative : LevelSetVelocity {

        ThermalParameters thermalParameters;
        XNSFE_OperatorConfiguration config;

        public LevelSetVelocityEvaporative(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters, XNSFE_OperatorConfiguration config) : base(levelSetName, D, degree, averagingMode, physicalParameters) {
            this.thermalParameters = config.getThermParams;
            this.config = config;
        }

        public override void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            // Evaporative part
            BitArray EvapMicroRegion = new BitArray(levelSet.Tracker.GridDat.Cells.Count);  //this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();

            double kA = 0, kB = 0, rhoA = 0, rhoB = 0;

            foreach(var spc in levelSet.Tracker.SpeciesNames) {
                switch (spc) {
                    case "A": { kA = thermalParameters.k_A; rhoA = thermalParameters.rho_A; break; }
                    case "B": { kB = thermalParameters.k_B; rhoB = thermalParameters.rho_B; break; }
                    default: { throw new ArgumentException("unknown species"); }
                }
            }

            XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];
            for (int d = 0; d < D; d++) {
                DGField evapVelocity = ParameterVarFields[ParameterNames[d]];

                int order = evapVelocity.Basis.Degree * evapVelocity.Basis.Degree + 2;
                evapVelocity.Clear();
                evapVelocity.ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1); // No nof Nodes

                       MultidimensionalArray VelA = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray VelB = MultidimensionalArray.Create(Len, K, D);

                       for (int dd = 0; dd < D; dd++) {
                           ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, VelA.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                           ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, VelB.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                       }

                       MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                       temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                       temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                       MultidimensionalArray HeatFluxA_Res = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray HeatFluxB_Res = MultidimensionalArray.Create(Len, K, D);

                       for (int dd = 0; dd < D; dd++) {
                           ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, HeatFluxA_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                           ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, HeatFluxB_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                       }


                       MultidimensionalArray TempA_Res = MultidimensionalArray.Create(Len, K);
                       MultidimensionalArray TempB_Res = MultidimensionalArray.Create(Len, K);
                       //MultidimensionalArray Curv_Res = MultidimensionalArray.Create(Len, K);
                       //MultidimensionalArray Pdisp_Res = MultidimensionalArray.Create(Len, K);

                       temperature.GetSpeciesShadowField("A").Evaluate(j0, Len, NS, TempA_Res);
                       temperature.GetSpeciesShadowField("B").Evaluate(j0, Len, NS, TempB_Res);
                       //ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature].Evaluate(j0, Len, NS, Curv_Res);
                       //this.DisjoiningPressure.Evaluate(j0, Len, NS, Pdisp_Res);

                       var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                       for (int j = 0; j < Len; j++) {

                           MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                           levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);

                           for (int k = 0; k < K; k++) {

                               double qEvap = 0.0;
                               if (EvapMicroRegion[j]) {
                                   throw new NotImplementedException("Check consistency for micro regions");
                                   // micro region
                                   //double Tsat = this.Control.ThermalParameters.T_sat;
                                   //double pc = this.Control.ThermalParameters.pc;
                                   //double pc0 = (pc < 0.0) ? this.Control.PhysicalParameters.Sigma * Curv_Res[j, k] + Pdisp_Res[j, k] : pc;
                                   //double f = this.Control.ThermalParameters.fc;
                                   //double R = this.Control.ThermalParameters.Rc;
                                   //if (this.Control.ThermalParameters.hVap_A > 0) {
                                   //    hVap = this.Control.ThermalParameters.hVap_A;
                                   //    rho_l = this.Control.PhysicalParameters.rho_A;
                                   //    rho_v = this.Control.PhysicalParameters.rho_B;
                                   //    double TintMin = Tsat * (1 + (pc0 / (hVap * rho_l)));
                                   //    double Rint = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rho_v * hVap.Pow2());
                                   //    if (TempA_Res[j, k] > TintMin)
                                   //        qEvap = -(TempA_Res[j, k] - TintMin) / Rint;
                                   //} else {
                                   //    hVap = -this.Control.ThermalParameters.hVap_A;
                                   //    rho_l = this.Control.PhysicalParameters.rho_B;
                                   //    rho_v = this.Control.PhysicalParameters.rho_A;
                                   //    double TintMin = Tsat * (1 + (pc0 / (hVap * rho_l)));
                                   //    double Rint = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rho_v * hVap.Pow2());
                                   //    if (TempB_Res[j, k] > TintMin)
                                   //        qEvap = (TempB_Res[j, k] - TintMin) / Rint;
                                   //}

                               } else {
                                   //macro region
                                   for (int dd = 0; dd < D; dd++) {
                                       qEvap += (HeatFluxB_Res[j, k, dd] - HeatFluxA_Res[j, k, dd]) * Normals[j, k, dd];
                                   }
                               }
                               //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
                               double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
                               double mEvap = (config.prescribedMassflux != null) ? config.prescribedMassflux(globX, time) : qEvap / thermalParameters.hVap; // mass flux
                                                                                                                                                             //double mEvap = qEvap / this.Control.ThermalParameters.hVap;
                                                                                                                                                             //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
                                                                                                                                                             //Console.WriteLine("prescribedMassFlux = {0}", this.XOpConfig.prescribedMassflux(globX, hack_Phystime));

                               double sNeg = VelA[j, k, d] + mEvap * (1 / rhoA) * Normals[j, k, d];
                               //Console.WriteLine("sNeg = {0}", sNeg);
                               double sPos = VelB[j, k, d] + mEvap * (1 / rhoB) * Normals[j, k, d];
                               //Console.WriteLine("sPos = {0}", sPos);

                               result[j, k] = (rhoA * sNeg + rhoB * sPos) / (rhoA + rhoB);   // density averaged evap velocity 
                           }
                       }
                   }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));

            } 

        }
        
    }

    class LevelSetVelocityGeneralNonMaterial : LevelSetVelocity {

        public LevelSetVelocityGeneralNonMaterial(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters) : base(levelSetName, D, degree, averagingMode, physicalParameters) {
        }

        public override void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            double rhoA = 0, rhoB = 0;

            foreach (var spc in levelSet.Tracker.SpeciesNames) {
                switch (spc) {
                    case "A": { rhoA = physicalParameters.rho_A; break; }
                    case "B": { rhoB = physicalParameters.rho_B; break; }
                    default: { throw new ArgumentException("unknown species"); }
                }
            }

            for (int d = 0; d < D; d++) {
                DGField interfaceVelocity = ParameterVarFields[ParameterNames[d]];

                int order = interfaceVelocity.Basis.Degree * interfaceVelocity.Basis.Degree + 2;
                interfaceVelocity.Clear();
                interfaceVelocity.ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1); // No nof Nodes

                       MultidimensionalArray VelA = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray VelB = MultidimensionalArray.Create(Len, K, D);

                       for (int dd = 0; dd < D; dd++) {
                           ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, VelA.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                           ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, VelB.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                       }

                       for (int j = 0; j < Len; j++) {

                           for (int k = 0; k < K; k++) {     
                               result[j, k] = (rhoA * VelA[j, k, d] - rhoB * VelB[j, k, d]) / (rhoA - rhoB);   // interface velocity for arbitrary mass flux
                           }
                       }
                   }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));

            }

        }

    }

    class MassFluxExtension_Evaporation : Parameter, ILevelSetParameter {

        XNSFE_OperatorConfiguration config;
        DualLevelSet levelSet;
        double time;
        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension };

        public override DelParameterFactory Factory => ParameterFactory;

        public MassFluxExtension_Evaporation(XNSFE_OperatorConfiguration config) {
            this.config = config;
            Update = MassFluxExtension_Evaporation_Update;
        }

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var massfluxext = new (string, DGField)[1];
            string paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;
            Basis b = new Basis(DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
            DGField MassFluxExtension = new SinglePhaseField(b, paramName);
            massfluxext[0] = (paramName, MassFluxExtension);
            return massfluxext;
        }

        public void MassFluxExtension_Evaporation_Update(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var thermalParams = config.getThermParams;
            double kA = 0, kB = 0;
            foreach (var spc in levelSet.Tracker.SpeciesNames) {
                switch (spc) {
                    case "A": { kA = thermalParams.k_A; break; }
                    case "B": { kB = thermalParams.k_B; break; }
                    default: { throw new ArgumentException("unknown species"); }
                }
            }
            var paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;

            SinglePhaseField MassFluxField = (SinglePhaseField)ParameterVarFields[paramName];
            int order = MassFluxField.Basis.Degree * MassFluxField.Basis.Degree + 2;

            XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];

            int D = levelSet.Tracker.GridDat.SpatialDimension;
            MassFluxField.Clear();
            MassFluxField.ProjectField(1.0,
                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes


                    MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                    temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                    temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                    //MultidimensionalArray HeatFluxA_Res = MultidimensionalArray.Create(Len, K, D);
                    //MultidimensionalArray HeatFluxB_Res = MultidimensionalArray.Create(Len, K, D);

                    //for (int dd = 0; dd < D; dd++) {
                    //    ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, HeatFluxA_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                    //    ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, HeatFluxB_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                    //}

                    var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {

                        MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                        levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);

                        for (int k = 0; k < K; k++) {

                            double qEvap = 0.0;
                            //macro region
                            for (int dd = 0; dd < D; dd++) {
                                qEvap += ((-kB) * GradTempB_Res[j, k, dd] - (-kA) * GradTempA_Res[j, k, dd]) * Normals[j, k, dd];
                                //qEvap += (HeatFluxB_Res[j, k, dd] - HeatFluxA_Res[j, k, dd]) * Normals[j, k, dd];
                            }

                            //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
                            double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
                            double mEvap = (config.prescribedMassflux != null) ? config.prescribedMassflux(globX, time) : qEvap / thermalParams.hVap; // mass flux

                            //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
                            result[j, k] = mEvap;
                        }
                    }
                }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            this.levelSet = levelSet;
            this.time = time;

            //var thermalParams = config.getThermParams;
            //double kA, kB;
            //foreach (var spc in levelSet.Tracker.SpeciesNames) {
            //    switch (spc) {
            //        case "A": { kA = thermalParams.k_A; break; }
            //        case "B": { kB = thermalParams.k_B; break; }
            //        default: { throw new ArgumentException("unknown species"); }
            //    }
            //}
            //var paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;
            //SinglePhaseField MassFluxExtensionField = (SinglePhaseField)ParameterVarFields[paramName];
            //SinglePhaseField MassFluxField = new SinglePhaseField(MassFluxExtensionField.Basis);
            //int order = MassFluxField.Basis.Degree * MassFluxField.Basis.Degree + 2;

            //XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];

            //int D = levelSet.Tracker.GridDat.SpatialDimension;
            //MassFluxField.Clear();
            //MassFluxField.ProjectField(1.0,
            //    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //        int K = result.GetLength(1); // No nof Nodes
                    

            //        MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
            //        MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

            //        temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
            //        temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

            //        MultidimensionalArray HeatFluxA_Res = MultidimensionalArray.Create(Len, K, D);
            //        MultidimensionalArray HeatFluxB_Res = MultidimensionalArray.Create(Len, K, D);
                    
            //        for (int dd = 0; dd < D; dd++) {
            //            ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, HeatFluxA_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
            //            ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, HeatFluxB_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
            //        }
                    

            //        var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

            //        for (int j = 0; j < Len; j++) {

            //            MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
            //            levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);

            //            for (int k = 0; k < K; k++) {

            //            double qEvap = 0.0;
            //            //macro region
            //            for (int dd = 0; dd < D; dd++) {                                
            //                    qEvap += (HeatFluxB_Res[j, k, dd] - HeatFluxA_Res[j, k, dd]) * Normals[j, k, dd];                                
            //            }

            //            //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
            //            double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
            //            double mEvap = (config.prescribedMassflux != null) ? config.prescribedMassflux(globX, time) : qEvap /thermalParams.hVap; // mass flux

            //            //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
            //            result[j, k] = mEvap;
            //            }
            //        }
            //    }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));     


            

            //SubGrid CCgrid = levelSet.Tracker.Regions.GetCutCellSubGrid();
            //CellMask CC = levelSet.Tracker.Regions.GetCutCellMask();
            //CellMask NEAR = levelSet.Tracker.Regions.GetNearFieldMask(1);
            //int J = levelSet.Tracker.GridDat.Cells.NoOfLocalUpdatedCells;
            //double[][] MassFluxMin = new double[1][];
            //double[][] MassFluxMax = new double[1][];
            //MassFluxMin[0] = new double[J];
            //MassFluxMax[0] = new double[J];

            //VectorField<SinglePhaseField> DGLevSetGradient =  new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(levelSet.DGLevelSet.Basis)));
            //DGLevSetGradient.Clear();
            //DGLevSetGradient.Gradient(1.0, levelSet.DGLevelSet);
            //NarrowMarchingBand.ConstructExtVel_PDE(levelSet.Tracker, CCgrid, new SinglePhaseField[] { MassFluxExtensionField }, new SinglePhaseField[] { MassFluxField },
            //    levelSet.DGLevelSet, DGLevSetGradient, MassFluxMin, MassFluxMax, order);

            //var marcher = new FastMarchReinit(levelSet.DGLevelSet.Basis);
            ////timestepnumber??
            //marcher.ConstructExtension(levelSet.DGLevelSet, NEAR.Except(CC), CC, new SinglePhaseField[] { MassFluxExtensionField },
            //    MassFluxMin, MassFluxMax, DGLevSetGradient, (int)time);

            //MassFluxExtensionField.CheckForNanOrInf(true, true, true);

            
        }
    }

    class Temperature0 : Parameter {
        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.Temperature0 };

        public override DelParameterFactory Factory => Temperature0Factory;

        public Temperature0() {
        }

        (string, DGField)[] Temperature0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var temperature0 = new (string, DGField)[1];
            
            string temperaturename = BoSSS.Solution.NSECommon.VariableNames.Temperature;
            DGField temperature = DomainVarFields[temperaturename];
            string paramName = BoSSS.Solution.NSECommon.VariableNames.Temperature0;
            temperature0[0] = (paramName, temperature);
            
            return temperature0;
        }
    }

    class HeatFlux0 : Parameter {

        int D;
        string[] parameters;
        LevelSetTracker lstrk;
        Dictionary<string, double> k_spc;
        public override IList<string> ParameterNames => parameters;

        public override DelParameterFactory Factory => HeatFlux0Factory;

        public HeatFlux0(int D, LevelSetTracker lstrk, XNSFE_OperatorConfiguration config) {
            this.D = D;
            this.lstrk = lstrk;
            parameters = BoSSS.Solution.NSECommon.VariableNames.HeatFlux0Vector(D);
            k_spc = new Dictionary<string, double>();
            var thermParams = config.getThermParams;
            foreach(var spc in lstrk.SpeciesNames) {
                switch (spc) {
                    case "A": { k_spc.Add(spc, thermParams.k_A); break; }
                    case "B": { k_spc.Add(spc, thermParams.k_B); break; }
                    default: { throw new ArgumentException("Unknown species.");}                    
                }
            }     
            if(config.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP)
                Update = HeatFlux0Update;
        }

        (string, DGField)[] HeatFlux0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var heatflux0 = new (string, DGField)[D];

            for (int d = 0; d < D; d++) {

                string heatfluxname = BoSSS.Solution.NSECommon.VariableNames.HeatFluxVectorComponent(d);
                string paramName = BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(d);

                XDGField heatflux;
                if (DomainVarFields.ContainsKey(heatfluxname)) {
                    // if the heatflux is a domainvariable, we can directly use him, no special update necessary
                    heatflux = (XDGField)DomainVarFields[heatfluxname];
                } else {
                    heatflux = new XDGField((XDGBasis)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature].Basis, paramName);
                }
                heatflux0[d] = (paramName, heatflux);
            }          

            return heatflux0;
        }
        private void HeatFlux0Update(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var varname = BoSSS.Solution.NSECommon.VariableNames.Temperature;
            DGField temperaure = DomainVarFields[varname];

            SubGrid Sgrd = lstrk.Regions.GetCutCellSubGrid();

            for (int d = 0; d < D; d++) {
                foreach(var spc in lstrk.SpeciesNames) {
                    
                    DGField temperaure_Spc = ((temperaure as XDGField).GetSpeciesShadowField(spc));

                    double aSpc = (k_spc != null && k_spc.ContainsKey(spc)) ? k_spc[spc] : 1.0;

                    var paramname = BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(d);
                    DGField heatflux = ParameterVarFields[paramname];
                    heatflux.Clear();
                    // q=-k*grad(T)
                    (heatflux as XDGField).GetSpeciesShadowField(spc).DerivativeByFlux(-aSpc, temperaure_Spc, d, Sgrd);
                }
            }            
        }

    }
}

