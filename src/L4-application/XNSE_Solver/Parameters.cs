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
            Factory = Velocity0Factory;
        }

        public override IList<string> Names => BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D);

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

    class Velocity0Mean : Parameter
    {
        int D;

        LevelSetTracker LsTrk;

        int cutCellQuadOrder;

        public Velocity0Mean(int D, LevelSetTracker LsTrk, int cutCellQuadOrder)
        {
            this.D = D;
            this.LsTrk = LsTrk;
            Factory = Velocity0MeanFactory;
            Update = Velocity0MeanUpdate;
            this.cutCellQuadOrder = cutCellQuadOrder;
        }

        public override IList<string> Names => BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D);

        (string, DGField)[] Velocity0MeanFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
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

        void Velocity0MeanUpdate(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            for(int d = 0; d < D; ++d)
            {
                foreach(string speciesName in LsTrk.SpeciesNames)
                {
                    XDGField paramMeanVelocity = (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]];
                    DGField speciesParam = paramMeanVelocity.GetSpeciesShadowField(speciesName);

                    XDGField velocity = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d]];
                    DGField speciesVelocity = velocity.GetSpeciesShadowField(speciesName);

                    //Uncut
                    speciesParam.SetMeanValueTo(speciesVelocity);

                    //Cut
                    CellMask cutCells = LsTrk.Regions.GetSpeciesMask(speciesName);
                    SpeciesId[ ] speciesIds = LsTrk.SpeciesIdS.ToArray();
                    XQuadSchemeHelper qh = LsTrk.GetXDGSpaceMetrics(speciesIds, cutCellQuadOrder).XQuadSchemeHelper;

                    SpeciesId speciesId = LsTrk.GetSpeciesId(speciesName);
                    double minvol = Math.Pow(this.LsTrk.GridDat.Cells.h_minGlobal, D);
                    CellQuadratureScheme scheme = qh.GetVolumeQuadScheme(speciesId, IntegrationDomain: cutCells);

                    SetMeanValueToMeanOf(speciesParam, speciesVelocity, minvol, cutCellQuadOrder, scheme);
                }
            }
        }

        void SetMeanValueToMeanOf(DGField target, DGField source, double minvol, int order, CellQuadratureScheme scheme)
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

        public Gravity(string species, int d, int D, Func<double[], double> initial, double rho, int degree)
        {
            this.degree = degree;
            Factory = GravityFactory;

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

        public override IList<string> Names => names;

        (string, DGField)[] GravityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            DGField gravity = new SinglePhaseField(basis, names[0]);
            gravity.ProjectField(initial);
            return new (string, DGField)[] { (names[0], gravity) };
        }
    }

    class Normals : Parameter
    {
        LevelSetTracker LsTrk;

        int D;

        public Normals(int D, LevelSetTracker LsTrk)
        {
            this.LsTrk = LsTrk;
            this.D = D;
            Factory = NormalFactory;
            Update = NormalUpdate;
        }

        public override IList<string> Names => BoSSS.Solution.NSECommon.VariableNames.NormalVector(D);

        (string, DGField)[] NormalFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            LevelSet Phi = (LevelSet)(LsTrk.LevelSets[0]);
            VectorField<SinglePhaseField> Normals = new VectorField<SinglePhaseField>(D, Phi.Basis, SinglePhaseField.Factory);
            Normals.Gradient(1.0, Phi);

            (string, DGField)[] normals = new (string, DGField)[D];
            for(int d = 0; d <D; ++d)
            {
                normals[d] = (BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[d], Normals[d] );
            }
            return normals;
        }

        void NormalUpdate(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            LevelSet Phi = (LevelSet)(LsTrk.LevelSets[0]);
            DGField[] Normals = new SinglePhaseField[D];
            for (int i = 0; i < D; ++i)
            {
                Normals[i] = ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[i]];
            }
            VectorField<DGField> normalVector = new VectorField<DGField>(Normals);
            Normals.Clear();
            normalVector.Gradient(1.0, Phi);
        }
    }
    
    class Curvature : Parameter
    {
        int D;

        LevelSetTracker LsTrk;

        int degree;

        int m_HMForder;

        bool solveKineticEnergyEquation;

        bool isEvaporation;

        DoNotTouchParameters AdvancedDiscretizationOptions;

        public Curvature(LevelSetTracker LsTrkr, int degree, int m_HMForder, bool isEvaporation, bool solveKineticEnergyEquation, DoNotTouchParameters AdvancedDiscretizationOptions)
        {
            this.LsTrk = LsTrkr;
            this.degree = degree;
            this.m_HMForder = m_HMForder;
            this.isEvaporation = isEvaporation;
            this.solveKineticEnergyEquation = solveKineticEnergyEquation;
            this.AdvancedDiscretizationOptions = AdvancedDiscretizationOptions;

            Factory = CurvatureFactory;
        }

        public static Curvature CreateFrom(XNSE_Control control, XNSFE_OperatorConfiguration xopConfig, LevelSetTracker LsTrkr, int m_HMForder)
        {
            string curvature = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            int degree;
            if (control.FieldOptions.TryGetValue(curvature, out FieldOpts opts))
            {
                degree = opts.Degree;
            }
            else
            {
                degree = 0;
            }
            bool solveKineticEnergyEquation = control.solveKineticEnergyEquation;

            bool isEvaporation = xopConfig.isEvaporation;

            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new Curvature(LsTrkr, degree, m_HMForder, isEvaporation, solveKineticEnergyEquation, AdvancedDiscretizationOptions);
        }

        public override IList<string> Names => new string[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };

        (string, DGField)[] CurvatureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            string name = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            SinglePhaseField curvature = new SinglePhaseField(basis, name);

            return new (string, DGField)[]
            {
                (name, curvature)
            };
        }
    }
}
