using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

    class Velocity0Mean : Parameter, ILevelSetParameter
    {
        int D;

        int cutCellQuadOrder;

        LevelSetTracker LsTrk;

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

        IList<string> SpeciesNames;

        LevelSetTracker.LevelSetRegions regions;

        IDictionary<string, SpeciesId> speciesMap;

        XQuadSchemeHelper schemeHelper;

        double minvol;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, LevelSetTracker lsTrkr, double time, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
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

        void Velocity0MeanUpdate(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
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

        static void SetMeanValueToMeanOf(DGField target, DGField source, double minvol, int order, CellQuadratureScheme scheme)
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

        public void LevelSetParameterUpdate(DualLevelSet levelSet, LevelSetTracker lsTrkr, double time, IReadOnlyDictionary<string, DGField> ParameterVarFields)
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
           LevelSetTracker lsTrkr,
           double time,
           IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            //dgLevSetGradient and update curvature
            VectorField<SinglePhaseField> filtLevSetGradient;
            CurvatureAlgorithms.LaplaceBeltramiDriver(
                AdvancedDiscretizationOptions.SST_isotropicMode,
                AdvancedDiscretizationOptions.FilterConfiguration,
                out filtLevSetGradient,
                lsTrkr,
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

        public BeltramiGradientAndCurvature(int curvatureDegree, int degree, int m_HMForder, DoNotTouchParameters AdvancedDiscretizationOptions, int D)
        {
            lsParameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            this.AdvancedDiscretizationOptions = AdvancedDiscretizationOptions;
            this.m_HMForder = m_HMForder;
            this.gradientDegree = degree;
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
                curvatureDegree = opts.Degree * 2;
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
           LevelSetTracker lsTrkr,
           double time,
           IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];
            VectorField<SinglePhaseField> filtLevSetGradient;
            CurvatureAlgorithms.CurvatureDriver(
                AdvancedDiscretizationOptions.SST_isotropicMode,
                AdvancedDiscretizationOptions.FilterConfiguration,
                Curvature,
                out filtLevSetGradient,
                lsTrkr,
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

        public void LevelSetParameterUpdate(DualLevelSet levelSet, LevelSetTracker lsTrkr, double time, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            DGField sigmaMax = ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.MaxSigma];

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
}

