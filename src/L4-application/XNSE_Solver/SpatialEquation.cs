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

namespace BoSSS.Application.XNSE_Solver 
{
    abstract class SpatialEquation {

        public SpatialEquation() 
        {
            Components = new LinkedList<IEquationComponent>();
        }

        public LinkedList<IEquationComponent> Components { get; set; }

        public abstract string CodomainName { get; }

        public string[] VariableNames { get; private set; }

        public string[] Parameters { get; private set; }

        public void AddVariableNames(params string[] names)
        {
            if(VariableNames == null)
            {
                VariableNames = names;
            }
            else
            {
                foreach(string name in names)
                {
                    if (!VariableNames.Contains(name))
                    {
                        VariableNames = VariableNames.Cat(name);
                    }
                }
            }
        }
        
        public void AddParameter(params string[] names)
        {
            if (Parameters == null)
            {
                Parameters = names;
            }
            else
            {
                foreach (string name in names)
                {
                    if (!Parameters.Contains(name))
                    {
                        Parameters = Parameters.Cat(name);
                    }
                }
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
        public LinkedList<SpatialEquation> SpatialEquations;

        public LinkedList<SurfaceEquation> InterfaceEquations;

        public LinkedList<BulkEquation> BulkEquations;

        public SystemOfEquations() {
            InterfaceEquations = new LinkedList<SurfaceEquation>();
            BulkEquations = new LinkedList<BulkEquation>();
            SpatialEquations = new LinkedList<SpatialEquation>();
        }

        public void AddEquation(SurfaceEquation A) {
            InterfaceEquations.AddLast(A);
            SpatialEquations.AddLast(A);
        }

        public void AddEquation(BulkEquation A) {
            BulkEquations.AddLast(A);
            SpatialEquations.AddLast(A);
        }

        public void AddEquation(SpatialEquation A)
        {
            SpatialEquations.AddLast(A);
        }

        public string[] DomainVars() {
            LinkedList<string> domainVars = new LinkedList<string>();
            foreach(SpatialEquation equation in BulkEquations) {
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

        public string[] CoDomainVars() {
            LinkedList<string> coDomainVars = new LinkedList<string>();
            foreach(BulkEquation equation in BulkEquations) {
                string coDomainVar = equation.CodomainName;
                if(!coDomainVars.Contains(coDomainVar)) {
                    coDomainVars.AddLast(coDomainVar);
                }
            }
            return coDomainVars.ToArray();
        }

        public string[] Parameters() {
            LinkedList<string> parameterNames = new LinkedList<string>();
            foreach(SpatialEquation equation in SpatialEquations) 
            {
                if(equation.Parameters != null)
                {
                    foreach (string parameter in equation.Parameters)
                    {
                        if (!parameterNames.Contains(parameter))
                        {
                            parameterNames.AddLast(parameter);
                        }
                    }
                }
            }
            return parameterNames.ToArray();
        }

        public string[] Species() {
            LinkedList<string> species = new LinkedList<string>();
            foreach(BulkEquation equation in BulkEquations) {
                string componentSpecies = equation.SpeciesName;
                if(componentSpecies != null && !species.Contains(componentSpecies)) {
                    species.AddLast(componentSpecies);
                }
            }
            return species.ToArray();
        }

        public (string, double[])[] MassDiagonal() 
        {
            StringArrayDictionary<StringArrayDictionary<double>> diag = new StringArrayDictionary<StringArrayDictionary<double>>(Species());
            foreach(BulkEquation equation in BulkEquations) 
            {
                if(diag.TryGetValue(equation.SpeciesName, out StringArrayDictionary<double> scales)) 
                {
                    scales.Add(equation.CodomainName, equation.MassScale);
                } 
                else 
                {
                    scales = new StringArrayDictionary<double>(CoDomainVars());
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
    }
}
