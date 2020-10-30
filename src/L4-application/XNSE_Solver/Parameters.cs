using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    class Velocity0 : IParameter
    {
        int D; 

        public Velocity0(int D)
        {
            this.D = D;
        }

        public IList<string> Names => BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D);

        public DelParameterFactory Factory => Velocity0Factory;

        public DelPartialParameterUpdate Update => Velocity0Update;

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

        void Velocity0Update(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            Console.WriteLine("Update Velocity0");
        }
    }

    class Velocity0Mean : IParameter
    {
        int D;

        LevelSetTracker LsTrk;

        public Velocity0Mean(int D, LevelSetTracker LsTrk)
        {
            this.D = D;
            this.LsTrk = LsTrk;
        }

        public IList<string> Names => BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D);

        public DelParameterFactory Factory => Velocity0MeanFactory;

        public DelPartialParameterUpdate Update => Velocity0MeanUpdate;

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

                    speciesParam.SetMeanValue(speciesVelocity);
                }
            }
            Console.WriteLine("Update Velocity0Mean");
            //Tecplot.PlotFields(new DGField[] { paramMeanVelocity, velocity }, "params-", 0, 3);
        }
    }

    class Gravity : IParameter
    {
        int degree;

        Func<double[], double> initial;

        string[] names;

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
            else
            {
                gravityDegree = 0;
            }
            return new Gravity(species, d, D, initialGravity, rho, gravityDegree);
        }

        public IList<string> Names => names;

        public DelParameterFactory Factory => GravityFactory;

        public DelPartialParameterUpdate Update => null;

        (string, DGField)[] GravityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            DGField gravity = new SinglePhaseField(basis, names[0]);
            gravity.ProjectField(initial);
            return new (string, DGField)[] { (names[0], gravity) };
        }
    }

    class Normals : IParameter
    {
        LevelSetTracker LsTrk;

        int D;

        public Normals(int D, LevelSetTracker LsTrk)
        {
            this.LsTrk = LsTrk;
            this.D = D;
        }

        public IList<string> Names => BoSSS.Solution.NSECommon.VariableNames.NormalVector(D);

        public DelParameterFactory Factory => NormalFactory;

        public DelPartialParameterUpdate Update => NormalUpdate;

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
            normalVector.Gradient(1.0, Phi);
        }
    }
    
    /*
    class Curvature : IParameter
    {
        
        public Curvature()
        {
            CurvatureAlgorithms.CurvatureDriver(
                this.Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                this.Control.AdvancedDiscretizationOptions.FilterConfiguration,
                this.Curvature, out filtLevSetGradient, this.LsTrk,
                this.m_HMForder,
                this.DGLevSet.Current);
        }

        public static Curvature CreateFrom(AppControl control)
        {

        }

        public IList<string> Names => new string[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };

        public DelParameterFactory Factory => CurvatureFactory;

        public DelPartialParameterUpdate Update => throw new NotImplementedException();

        (string, DGField)[] CurvatureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            LevelSet Phi = (LevelSet)(LsTrk.LevelSets[0]);
            VectorField<SinglePhaseField> Normals = new VectorField<SinglePhaseField>(D, Phi.Basis, SinglePhaseField.Factory);
            Normals.Gradient(1.0, Phi);

            (string, DGField)[] normals = new (string, DGField)[D];
            for (int d = 0; d < D; ++d)
            {
                normals[d] = (BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)[d], Normals[d]);
            }
            return normals;
        }

    }
    */
}
