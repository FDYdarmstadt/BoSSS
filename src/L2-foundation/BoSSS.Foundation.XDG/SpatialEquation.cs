using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Foundation.XDG.OperatorFactory {
    /// <summary>
    /// Base class for all equations of the operator factory.
    /// A spatial equation is an equation for one codomain. 
    /// Implement this class, if you have a single phase equation. 
    /// </summary>
    public abstract class SpatialEquation : ICodomainEquation<SpatialEquation> {
        /// <summary>
        /// Empty spatial equation.
        /// </summary>
        public SpatialEquation() {
            Components = new LinkedList<IEquationComponent>();
        }

        /// <summary>
        /// Equation components of this equation. All components belong to a single codomain.
        /// </summary>
        public LinkedList<IEquationComponent> Components { get; set; }

        /// <summary>
        /// Single codomain name of this equation. A codomain a name for the row in a system of equations.
        /// </summary>
        public abstract string CodomainName { get; }

        /// <summary>
        /// All names of variables in equation. Must match variables of components.
        /// </summary>
        public string[] VariableNames { get; private set; }

        /// <summary>
        /// All names of parameters in equation. Must match parameters of components.
        /// </summary>
        public string[] Parameters { get; private set; }

        /// <summary>
        /// All names of coefficients. Must match coefficients of components.
        /// </summary>
        public string[] Coefficients { get; private set; }

        /// <summary>
        /// Add name of variable. Will only add, if name was not already added.
        /// </summary>
        /// <param name="names"></param>
        public void AddVariableNames(params string[] names) {
            if(VariableNames == null) {
                VariableNames = names;
            } else {
                foreach(string name in names) {
                    if(!VariableNames.Contains(name)) {
                        VariableNames = VariableNames.Cat(name);
                    }
                }
            }
        }

        /// <summary>
        /// Add name of parameter. Will only add, if name was not already added.
        /// </summary>
        /// <param name="names"></param>
        public void AddParameter(params string[] names) {
            if(Parameters == null) {
                Parameters = names;
            } else {
                foreach(string name in names) {
                    if(!Parameters.Contains(name)) {
                        Parameters = Parameters.Cat(name);
                    }
                }
            }
        }

        /// <summary>
        /// Add name of coefficient. Will only add, if name was not already added.
        /// </summary>
        /// <param name="names"></param>
        public void AddCoefficient(params string[] names) {
            if(Coefficients == null) {
                Coefficients = names;
            } else {
                foreach(string name in names) {
                    if(!Coefficients.Contains(name)) {
                        Coefficients = Coefficients.Cat(name);
                    }
                }
            }
        }

        /// <summary>
        /// Add component to this equation. The variable/parameter/coefficient names must be added separately.
        /// </summary>
        /// <param name="component">Equation component for this equation's codomain</param>
        public void AddComponent(IEquationComponent component) {
            Components.AddLast(component);
        }

        public void Combine(SpatialEquation other) {
            if(other.CodomainName != this.CodomainName) {
                throw new Exception("Can only add eqations of same codomain name");
            }
            AddVariableNames(other.VariableNames);
            if(other.Parameters != null) {
                AddParameter(other.Parameters);
            }
            if(other.Coefficients != null) {
                AddCoefficient(other.Coefficients);
            }
            Components.AddRange(other.Components);
        }

        public bool EqualCodomain(SpatialEquation other) {
            return CodomainName == other.CodomainName;
        }
    }

    /// <summary>
    /// XDG Equations on a surface. 
    /// Surface is between to species.
    /// </summary>
    public abstract class SurfaceEquation : SpatialEquation, ICodomainEquation<SurfaceEquation>{
        /// <summary>
        /// Empty surface equation.
        /// </summary>
        public SurfaceEquation() {
            SurfaceComponents = new LinkedList<IEquationComponent>();
            ContactLineComponents = new LinkedList<IEquationComponent>();
        }

        /// <summary>
        /// First/Negative species (with respect to level-set)
        /// Order not important.
        /// </summary>
        public abstract string FirstSpeciesName { get; }

        /// <summary>
        /// Second/Positive species (with respect to level-set)
        /// Order not important.
        /// </summary>
        public abstract string SecondSpeciesName { get; }

        /// <summary>
        /// Specialized surface components that are part of a special spatial operator in the general spatial operator.
        /// Generally, surface equations are also collected in <see cref="SpatialEquation.Components"/>.
        /// </summary>
        public LinkedList<IEquationComponent> SurfaceComponents { get; private set; }

        /// <summary>
        /// Specialized contact line components that are part of a special spatial operator in the general spatial operator.
        /// </summary>
        public LinkedList<IEquationComponent> ContactLineComponents { get; private set; }

        /// <summary>
        /// Add Specialized surface components that will be part of a special spatial operator.
        /// Generally, surface equations are also added via <see cref="SpatialEquation.AddComponent(IEquationComponent)"/>.
        /// </summary>
        /// <param name="surfaceComponent"> Specialized surface component </param>
        public void AddSurfaceComponent(IEquationComponent surfaceComponent) {
            SurfaceComponents.AddLast(surfaceComponent);
        }

        /// <summary>
        /// Add Specialized surface components that will be part of a special spatial operator.
        /// Generally, surface equations are also added via <see cref="SpatialEquation.AddComponent(IEquationComponent)"/>.
        /// </summary>
        /// <param name="surfaceComponent"> Specialized surface component </param>
        public void AddContactLineComponent(IEquationComponent surfaceComponent) {
            ContactLineComponents.AddLast(surfaceComponent);
        }

        public void Combine(SurfaceEquation other) {
            if(other.FirstSpeciesName != FirstSpeciesName || other.SecondSpeciesName != SecondSpeciesName) {
                throw new Exception("Species Names do not match");
            }
            base.Combine(other);
            SurfaceComponents.AddRange(other.SurfaceComponents);
            ContactLineComponents.AddRange(other.ContactLineComponents);
        }

        public bool EqualCodomain(SurfaceEquation other) {
            bool isEqual = CodomainName == other.CodomainName;
            isEqual &= (FirstSpeciesName == other.FirstSpeciesName);
            isEqual &= (SecondSpeciesName == other.SecondSpeciesName);
            return isEqual;
        }
    }

    /// <summary>
    /// XDG Equations for bulk phase.
    /// </summary>
    public abstract class BulkEquation : SpatialEquation, ICodomainEquation<BulkEquation> {
        
        /// <summary>
        /// Name of species for which equation is valid.
        /// </summary>
        public abstract string SpeciesName { get; }

        /// <summary>
        /// This is required, when the whole equation is scaled, e.g. by density. 
        /// In time stepping, the respective entries of the mass matrix will be scaled by this factor.
        /// </summary>
        public abstract double MassScale { get; }

        /// <summary>
        /// Special component of an extra ghost spatial operator.
        /// </summary>
        public LinkedList<IEquationComponent> GhostComponents { get; private set; }

        /// <summary>
        /// Empty bulk equation.
        /// </summary>
        public BulkEquation() {
            GhostComponents = new LinkedList<IEquationComponent>();
        }

        /// <summary>
        /// Add a specialized ghost component.
        /// </summary>
        /// <param name="ghostComponent"></param>
        public void AddGhostComponent(IEquationComponent ghostComponent) {
            GhostComponents.AddLast(ghostComponent);
        }

        public void Combine(BulkEquation other) {
            if (other.SpeciesName != SpeciesName || other.MassScale != MassScale) {
                throw new Exception("species name or mass scale do not match");
            };
            base.Combine(other);
            GhostComponents.AddRange(other.GhostComponents);
        }

        public bool EqualCodomain(BulkEquation other) {
            bool isEqual = CodomainName == other.CodomainName;
            isEqual &= (SpeciesName == other.SpeciesName);
            return isEqual;
        }
    }


    interface ICodomainEquation<T> {
        void Combine(T other);

        bool EqualCodomain(T other);
    }

    /// <summary>
    /// List of Equations. There is only one equation for each codomain species pair.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    class EquationList<T> : IEnumerable<T> where T : ICodomainEquation<T> {
        LinkedList<T> equations;

        public EquationList() {
            equations = new LinkedList<T>();
        }

        /// <summary>
        ///  
        /// </summary>
        /// <param name="equation">
        /// If the codomain of equation is already present, equation will be combined with the present entry.
        /// </param>
        public void Add(T equation) {
            LinkedListNode<T> node = equations.First;
            while(node != null) {
                T q = node.Value;
                if(q.EqualCodomain(equation)) {
                    T combo = Combine(equation, q);
                    node.Value = combo;
                    return;
                }
                node = node.Next;
            }
            equations.AddLast(equation);
        }

        public IEnumerator<T> GetEnumerator() {
            return equations.GetEnumerator();
        }

        T Combine(T a, T b) {
            a.Combine(b);
            return a;
        }

        IEnumerator IEnumerable.GetEnumerator() {
            return equations.GetEnumerator();
        }
    }


    /// <summary>
    /// Accumulation of equations
    /// Finds the names of all required parameters and coefficients
    /// </summary>
    class SystemOfEquations {
        public EquationList<SpatialEquation> SpatialEquations;

        public EquationList<SurfaceEquation> InterfaceEquations;

        public EquationList<BulkEquation> BulkEquations;

        /// <summary>
        /// Create empty System
        /// </summary>
        public SystemOfEquations() {
            InterfaceEquations = new EquationList<SurfaceEquation>();
            BulkEquations = new EquationList<BulkEquation>();
            SpatialEquations = new EquationList<SpatialEquation>();
        }

        IEnumerable<SpatialEquation> AllEquations() {
            foreach (SpatialEquation equation in SpatialEquations) {
                yield return equation;
            }
            foreach (SpatialEquation equation in BulkEquations) {
                yield return equation;
            }
            foreach (SpatialEquation equation in InterfaceEquations) {
                yield return equation;
            }
        }

        public void AddEquation(SurfaceEquation A) {
            InterfaceEquations.Add(A);
        }

        public void AddEquation(BulkEquation A) {
            BulkEquations.Add(A);
        }

        public void AddEquation(SpatialEquation A) {
            SpatialEquations.Add(A);
        }

        public string[] DomainVars() {
            LinkedList<string> domainVars = new LinkedList<string>();
            foreach(SpatialEquation equation in AllEquations()) {
                if(equation.VariableNames != null) {
                    if(domainVars.Count == 0) {
                        domainVars.AddRange(equation.VariableNames);
                    } else {
                        // a real, sorted merging would be nicer, but it is really tricky to implement
                        foreach(string n in equation.VariableNames) {
                            if(!domainVars.Contains(n)) {
                                domainVars.AddLast(n);
                            }
                        }
                    }
                }
            }
            return domainVars.ToArray();
        }

        public string[] CoDomainVars() {
            LinkedList<string> coDomainVars = new LinkedList<string>();
            foreach(SpatialEquation equation in AllEquations()) {
                string coDomainVar = equation.CodomainName;
                if(!coDomainVars.Contains(coDomainVar)) {
                    coDomainVars.AddLast(coDomainVar);
                }
            }
            return coDomainVars.ToArray();
        }

        public string[] Parameters() {
            LinkedList<string> parameterNames = new LinkedList<string>();
            foreach(SpatialEquation equation in AllEquations()) {
                if(equation.Parameters != null) {
                    foreach(string parameter in equation.Parameters) {
                        if(!parameterNames.Contains(parameter)) {
                            parameterNames.AddLast(parameter);
                        }
                    }
                }
            }
            return parameterNames.ToArray();
        }

        public string[] Coefficients() {
            LinkedList<string> coefficientsNames = new LinkedList<string>();
            foreach(SpatialEquation equation in AllEquations()) {
                if(equation.Coefficients != null) {
                    foreach(string coefficient in equation.Coefficients) {
                        if(!coefficientsNames.Contains(coefficient)) {
                            coefficientsNames.AddLast(coefficient);
                        }
                    }
                }
            }
            return coefficientsNames.ToArray();
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

        public (string, double[])[] MassDiagonal() {
            StringArrayDictionary<StringArrayDictionary<double>> diag = new StringArrayDictionary<StringArrayDictionary<double>>(Species());
            foreach(BulkEquation equation in BulkEquations) {
                if(diag.TryGetValue(equation.SpeciesName, out StringArrayDictionary<double> scales)) {
                    scales.Add(equation.CodomainName, equation.MassScale);
                } else {
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
            for(int i = 0; i < diagonal.Length; ++i) {
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
                    throw new Exception("value has already been initialized");
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
