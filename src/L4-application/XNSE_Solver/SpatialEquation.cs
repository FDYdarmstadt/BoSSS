using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {

    struct Parameter
    {
        public string Name;

        public DelParameterFactory Factory;

        public DelPartialParameterUpdate Update;
    }

    abstract class SpatialEquation {

        public SpatialEquation() 
        {
            Components = new LinkedList<IEquationComponent>();
        }

        public LinkedList<IEquationComponent> Components { get; set; }

        public abstract string CodomainName { get; }

        public string[] VariableNames { get; private set; }

        public Parameter[] Parameters { get; private set; }

        public void AddVariableNames(params string[] names)
        {
            if(VariableNames == null)
            {
                VariableNames = names;
            }
            else
            {
                VariableNames = VariableNames.Cat(names);
            }
        }
        
        public void AddParameter(string name, DelParameterFactory factory = null, DelPartialParameterUpdate update = null)
        {
            Parameter parameter = new Parameter
            {
                Name = name,
                Factory = factory,
                Update = update
            };
            if (Parameters == null)
            {
                Parameters = new Parameter[] { parameter };
            }
            else
            {
                Parameters = Parameters.Cat(parameter);
            }
        }

        public void AddComponent(IEquationComponent component) 
        {
            Components.AddLast(component);
        }
    }

    abstract class SurfaceEquation : SpatialEquation 
    {
        public abstract string FirstSpeciesName { get; }

        public abstract string SecondSpeciesName { get; }

        public LinkedList<IEquationComponent> SurfaceComponents { get; private set; }

        public SurfaceEquation() {
            SurfaceComponents = new LinkedList<IEquationComponent>();
        }

        public void AddSurfaceComponent(IEquationComponent surfaceComponent) {
            SurfaceComponents.AddLast(surfaceComponent);
        }
    }

    abstract class BulkEquation : SpatialEquation 
    {
        public abstract string SpeciesName { get; }

        public abstract double MassScale { get; }

        public LinkedList<IEquationComponent> GhostComponents { get; private set; }

        public BulkEquation() {
            GhostComponents = new LinkedList<IEquationComponent>();
        }

        public void AddGhostComponent(IEquationComponent ghostComponent) {
            GhostComponents.AddLast(ghostComponent);
        }
    }

    class SystemOfEquations 
    {
        protected LinkedList<SurfaceEquation> interfaceEquations;

        protected LinkedList<BulkEquation> bulkEquations;

        public SystemOfEquations() {
            interfaceEquations = new LinkedList<SurfaceEquation>();
            bulkEquations = new LinkedList<BulkEquation>();
        }

        public void AddEquation(SurfaceEquation A) {
            interfaceEquations.AddLast(A);
        }

        public void AddEquation(BulkEquation A) {
            bulkEquations.AddLast(A);
        }

        public virtual XSpatialOperatorMk2 GetSpatialOperator(Func<int[], int[], int[], int> QuadOrderFunc) {
            var spatialOperator = CreateSpatialOperator(QuadOrderFunc);
            AddEquationComponents(spatialOperator);
            AddSurfaceEquationComponents(spatialOperator);
            AddGhostEquationComponents(spatialOperator);
            AddTemporalOperator(spatialOperator);
            AddParameterDelegates(spatialOperator);

            return spatialOperator;
        }

        XSpatialOperatorMk2 CreateSpatialOperator(Func<int[], int[], int[], int> QuadOrderFunc) {
            string[] domainVars = ExtractDomainVarsFromEquations();
            string[] codomainVars = ExtractCoDomainVarsFromEquations();
            string[] parameters = ExtractParametersFromEquations();
            string[] species = ExtractSpeciesFromEquations();

            var spatialOperator = new XSpatialOperatorMk2(
                domainVars,
                parameters,
                codomainVars,
                QuadOrderFunc,
                species);
            return spatialOperator;
        }

        void AddEquationComponents(XSpatialOperatorMk2 spatialOperator) {
            foreach(SpatialEquation equation in interfaceEquations) {
                foreach(IEquationComponent component in equation.Components) {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
            foreach(SpatialEquation equation in bulkEquations) {
                foreach(IEquationComponent component in equation.Components) {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
        }

        void AddSurfaceEquationComponents(XSpatialOperatorMk2 spatialOperator) {
            foreach(SurfaceEquation equation in interfaceEquations) {
                if(equation.SurfaceComponents != null) {
                    foreach(IEquationComponent component in equation.SurfaceComponents) {
                        spatialOperator.SurfaceElementOperator.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
            }
        }

        void AddGhostEquationComponents(XSpatialOperatorMk2 spatialOperator) {
            foreach(BulkEquation equation in bulkEquations) {
                if(equation.GhostComponents != null) {
                    foreach(IEquationComponent component in equation.GhostComponents) {
                        spatialOperator.GhostEdgesOperator.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
            }
        }

        string[] ExtractDomainVarsFromEquations() {
            LinkedList<string> domainVars = new LinkedList<string>();
            foreach(SpatialEquation equation in bulkEquations) {
                if(domainVars.Count == 0) 
                {
                    domainVars.AddRange(equation.VariableNames);
                } 
                else 
                {
                    // a real, sorted merging would be nicer, but it is really tricky to implement
                    foreach(string n in equation.VariableNames)
                    {
                        if (!domainVars.Contains(n))
                        {
                            domainVars.AddLast(n);
                        }
                    }
                }
            }
            return domainVars.ToArray();
        }

        string[] ExtractCoDomainVarsFromEquations() {
            LinkedList<string> coDomainVars = new LinkedList<string>();
            foreach(BulkEquation equation in bulkEquations) {
                string coDomainVar = equation.CodomainName;
                if(!coDomainVars.Contains(coDomainVar)) {
                    coDomainVars.AddLast(coDomainVar);
                }
            }
            return coDomainVars.ToArray();
        }

        string[] ExtractParametersFromEquations() {
            LinkedList<string> parameterNames = new LinkedList<string>();
            foreach(SpatialEquation equation in interfaceEquations) 
            {
                if(equation.Parameters != null)
                {
                    foreach (Parameter parameter in equation.Parameters)
                    {
                        if (!parameterNames.Contains(parameter.Name))
                        {
                            parameterNames.AddLast(parameter.Name);
                        }
                    }
                }
            }
            foreach(SpatialEquation equation in bulkEquations) 
            {
                if (equation.Parameters != null)
                {
                    foreach (Parameter parameter in equation.Parameters)
                    {
                        if (!parameterNames.Contains(parameter.Name))
                        {
                            parameterNames.AddLast(parameter.Name);
                        }
                    }
                }
            }
            return VariableNames.Velocity0Vector(2).Cat(VariableNames.Velocity0MeanVector(2));
        }

        string[] ExtractSpeciesFromEquations() {
            LinkedList<string> species = new LinkedList<string>();
            foreach(BulkEquation equation in bulkEquations) {
                string componentSpecies = equation.SpeciesName;
                if(componentSpecies != null && !species.Contains(componentSpecies)) {
                    species.AddLast(componentSpecies);
                }
            }
            return species.ToArray();
        }

        void AddTemporalOperator(XSpatialOperatorMk2 spatialOperator) 
        {
            (string, double[])[] diagonal = GetMassDiagonal(spatialOperator);
            spatialOperator.TemporalOperator = new ConstantXTemporalOperator(spatialOperator, diagonal);
        }

        (string, double[])[] GetMassDiagonal(XSpatialOperatorMk2 spatialOperator) 
        {
            StringArrayDictionary<StringArrayDictionary<double>> diag = new StringArrayDictionary<StringArrayDictionary<double>>(spatialOperator.Species);
            foreach(BulkEquation equation in bulkEquations) 
            {
                if(diag.TryGetValue(equation.SpeciesName, out StringArrayDictionary<double> scales)) 
                {
                    scales.Add(equation.CodomainName, equation.MassScale);
                } 
                else 
                {
                    scales = new StringArrayDictionary<double>(spatialOperator.CodomainVar);
                    scales.Add(equation.CodomainName, equation.MassScale);
                    diag.Add(equation.SpeciesName, scales);
                }
            }
            (string, double[])[] diagonal = ToTupleArray(diag);
            return diagonal;
        }

        static (string, double[])[] ToTupleArray(StringArrayDictionary<StringArrayDictionary<double>> d) {
            (string, double[])[] diagonal = new (string, double[])[d.Length];
            for(int i = 0; i < diagonal.Length; ++i) 
            {
                diagonal[i] = (d.Keys[i], d.Values[i].Values);
            }
            return diagonal;
        }

        class StringArrayDictionary<T> {
            bool[] isInitialized;

            public string[] Keys;

            public T[] Values;

            public int Length;

            public StringArrayDictionary(ICollection<string> keys) {
                Length = keys.Count;
                isInitialized = new bool[Length];
                Values = new T[Length];
                Keys = new string[Length];
                int i = 0;
                foreach(string stringKey in keys) {
                    isInitialized[i] = false;
                    Keys[i] = stringKey;
                    ++i;
                }
            }

            public void Add(string key, T value) {
                int pointer = ValuePointer(key);
                if(isInitialized[pointer]) {
                    throw new Exception("already initialized");
                }
                Values[pointer] = value;
                isInitialized[pointer] = true;
            }

            int ValuePointer(string key) {
                for(int i = 0; i < Length; ++i) {
                    if(Keys[i] == key) {
                        return i;
                    }
                }
                throw new Exception("key not found");
            }

            public bool TryGetValue(string key, out T value) {
                for(int i = 0; i < Length; ++i) {
                    if(Keys[i] == key) {
                        if(isInitialized[i]) {
                            value = Values[i];
                            return true;
                        } else {
                            break;
                        }
                    }
                }
                value = default(T);
                return false;
            }
        }

        void AddParameterDelegates(XSpatialOperatorMk2 spatialOperator)
        {
            foreach (SpatialEquation equation in interfaceEquations)
            {
                if (equation.Parameters != null)
                {
                    foreach (Parameter parameter in equation.Parameters)
                    {
                        if(parameter.Factory != null)
                        {
                            spatialOperator.ParameterFactories.Add(parameter.Factory);
                        }
                        if(parameter.Update != null)
                        {
                            spatialOperator.ParameterUpdates.Add(parameter.Update);
                        }
                    }
                }
            }
            foreach (SpatialEquation equation in bulkEquations)
            {
                if (equation.Parameters != null)
                {
                    foreach (Parameter parameter in equation.Parameters)
                    {
                        if (parameter.Factory != null)
                        {
                            spatialOperator.ParameterFactories.Add(parameter.Factory);
                        }
                        if (parameter.Update != null)
                        {
                            spatialOperator.ParameterUpdates.Add(parameter.Update);
                        }
                    }
                }
            }
            
        }
    }
}
