using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    abstract class SpatialEquation
    {
        public SpatialEquation()
        {
            Components = new LinkedList<IEquationComponent>();
        }

        public LinkedList<IEquationComponent> Components { get; set; }

        public string CodomainName { get; set; }

        public void AddComponent(IEquationComponent component)
        {
            Components.AddLast(component);
        }
    }

    abstract class SurfaceEquation : SpatialEquation
    {
        public string FirstSpeciesName { get; set; }
        
        public string SecondSpeciesName { get; set; }

        public LinkedList<IEquationComponent> SurfaceComponents { get; set; }

        public SurfaceEquation()
        {
            SurfaceComponents = new LinkedList<IEquationComponent>();
        }

        public void AddSurfaceComponent(IEquationComponent surfaceComponent)
        {
            SurfaceComponents.AddLast(surfaceComponent);
        }
    }

    abstract class BulkEquation : SpatialEquation
    {
        public string SpeciesName { get; set; }

        public double MassScale { get; set; }

        public LinkedList<IEquationComponent> GhostComponents { get; set; }

        public BulkEquation()
        {
            GhostComponents = new LinkedList<IEquationComponent>();
        }

        public void AddGhostComponent(IEquationComponent ghostComponent)
        {
            GhostComponents.AddLast(ghostComponent);
        }
    }

    class SystemOfEquations
    {
        protected LinkedList<SurfaceEquation> interfaceEquations;
        
        protected LinkedList<BulkEquation> bulkEquations;

        public SystemOfEquations()
        {
            interfaceEquations = new LinkedList<SurfaceEquation>();
            bulkEquations = new LinkedList<BulkEquation>();
        }

        public void AddEquation(SurfaceEquation A)
        {
            interfaceEquations.AddLast(A);
        }

        public void AddEquation(BulkEquation A)
        {
            bulkEquations.AddLast(A);
        }

        public virtual XSpatialOperatorMk2 GetSpatialOperator(int quadratureOrder)
        {
            var spatialOperator = CreateSpatialOperator(quadratureOrder);
            AddEquationComponents(spatialOperator);
            AddSurfaceEquationComponents(spatialOperator);
            AddGhostEquationComponents(spatialOperator);
            AddTemporalOperator(spatialOperator);
            
            return spatialOperator;
        }

        XSpatialOperatorMk2 CreateSpatialOperator(int quadratureOrder)
        {
            string[] domainVars = ExtractDomainVarsFromEquations();
            string[] codomainVars = ExtractCoDomainVarsFromEquations();
            string[] parameters = ExtractParametersFromEquations();
            string[] species = ExtractSpeciesFromEquations();

            var spatialOperator = new XSpatialOperatorMk2(
                domainVars,
                parameters,
                codomainVars,
                (A, B, C) => quadratureOrder,
                species);
            return spatialOperator;
        }

        void AddEquationComponents(XSpatialOperatorMk2 spatialOperator)
        {
            foreach(SpatialEquation equation in interfaceEquations)
            {
                foreach (IEquationComponent component in equation.Components)
                {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
            foreach (SpatialEquation equation in bulkEquations)
            {
                foreach (IEquationComponent component in equation.Components)
                {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
        }

        void AddSurfaceEquationComponents(XSpatialOperatorMk2 spatialOperator)
        {
            foreach (SurfaceEquation equation in interfaceEquations)
            {
                if(equation.SurfaceComponents != null)
                {
                    foreach (IEquationComponent component in equation.SurfaceComponents)
                    {
                        spatialOperator.SurfaceElementOperator.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
            }
        }

        void AddGhostEquationComponents(XSpatialOperatorMk2 spatialOperator)
        {
            foreach (BulkEquation equation in bulkEquations)
            {
                if (equation.GhostComponents != null)
                {
                    foreach (IEquationComponent component in equation.GhostComponents)
                    {
                        spatialOperator.GhostEdgesOperator.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
            }
        }

        string[] ExtractDomainVarsFromEquations()
        {
            LinkedList<string> domainVars = new LinkedList<string>();
            foreach (SpatialEquation equation in bulkEquations)
            {
                foreach (IEquationComponent component in equation.Components)
                {
                    IList<string> equationVariables = component.ArgumentOrdering;
                    foreach (string equationVariable in equationVariables)
                    {
                        if (!domainVars.Contains(equationVariable))
                        {
                            domainVars.AddLast(equationVariable);
                        }
                    }
                }
            }
            return domainVars.ToArray();
        }

        string[] ExtractCoDomainVarsFromEquations()
        {
            LinkedList<string> coDomainVars = new LinkedList<string>();
            foreach (BulkEquation equation in bulkEquations)
            {
                string coDomainVar = equation.CodomainName;
                if (!coDomainVars.Contains(coDomainVar))
                {
                    coDomainVars.AddLast(coDomainVar);
                }
            }
            return coDomainVars.ToArray();
        }

        string[] ExtractParametersFromEquations()
        {
            LinkedList<string> parameters = new LinkedList<string>();
            foreach (SpatialEquation equation in interfaceEquations)
            {
                foreach (IEquationComponent component in equation.Components)
                {
                    IList<string> equationParameters = component.ParameterOrdering;
                    if (equationParameters != null)
                    {
                        foreach (string equationparameter in equationParameters)
                        {
                            if (!parameters.Contains(equationparameter))
                            {
                                parameters.AddLast(equationparameter);
                            }
                        }
                    }
                }
            }
            foreach (SpatialEquation equation in bulkEquations)
            {
                foreach (IEquationComponent component in equation.Components)
                {
                    IList<string> equationParameters = component.ParameterOrdering;
                    if (equationParameters != null)
                    {
                        foreach (string equationparameter in equationParameters)
                        {
                            if (!parameters.Contains(equationparameter))
                            {
                                parameters.AddLast(equationparameter);
                            }
                        }
                    }
                }
            }
            return parameters.ToArray();
        }

        string[] ExtractSpeciesFromEquations()
        {
            LinkedList<string> species = new LinkedList<string>();
            foreach (BulkEquation equation in bulkEquations)
            {
                string componentSpecies = equation.SpeciesName;
                if (componentSpecies != null && !species.Contains(componentSpecies))
                {
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
            foreach (BulkEquation equation in bulkEquations)
            {
                if (diag.TryGetValue(equation.SpeciesName, out StringArrayDictionary<double> scales))
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

        static (string, double[])[] ToTupleArray(StringArrayDictionary<StringArrayDictionary<double>> d)
        {
            (string, double[])[] diagonal = new (string, double[])[d.Length];
            for(int i = 0; i < diagonal.Length;  ++i)
            {
                diagonal[i] = (d.Keys[i], d.Values[i].Values);
            }
            return diagonal;
        }

        class StringArrayDictionary<T>
        {
            bool[] isInitialized;

            public string[] Keys;

            public T[] Values;

            public int Length;

            public StringArrayDictionary(ICollection<string> keys)
            {
                Length = keys.Count;
                isInitialized = new bool[Length];
                Values = new T[Length];
                Keys = new string[Length];
                int i = 0;
                foreach (string stringKey in keys)
                {
                    isInitialized[i] = false;
                    Keys[i] = stringKey;
                    ++i;
                }
            }

            public void Add(string key, T value)
            {
                int pointer = ValuePointer(key);
                if (isInitialized[pointer])
                {
                    throw new Exception("already initialized");
                }
                Values[pointer] = value;
                isInitialized[pointer] = true;
            }

            int ValuePointer(string key)
            {
                for (int i = 0; i < Length; ++i)
                {
                    if (Keys[i] == key)
                    {
                        return i;
                    }
                }
                throw new Exception("key not found");
            }

            public bool TryGetValue(string key, out T value)
            {
                for (int i = 0; i < Length; ++i)
                {
                    if (Keys[i] == key)
                    {
                        if (isInitialized[i])
                        {
                            value = Values[i];
                            return true;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                value = default(T);
                return false;
            }
        }
    }
}
