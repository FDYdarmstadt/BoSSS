using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Foundation.XDG.OperatorFactory {
    /// <summary>
    /// Factory to create a spatial operator.
    /// Usage:
    /// 1. Create System by:
    ///    - Adding Equations: <see cref="AddEquation(BulkEquation)"/> (and variants)
    ///    - Adding Parameters (Parameters are fields, e.g. spatially dependent viscosity): <see cref="AddParameter(ParameterS)"/>
    ///    - Adding Coefficients (Coefficients are single numbers, e.g. the Reynolds number): <see cref="AddCoefficient(Coefficient)"/>
    /// 2. Create spatial operator by calling GetSpatialOperator 
    /// </summary>
    public class OperatorFactory {
        SystemOfEquations eqSystem;

        class ParameterList {

            public void VerifyList(IEnumerable<string> parameterNames) {
                foreach(var p in this.parameters) {
                    foreach(string OtherPname in p.ParameterNames) {
                        if(!parameterNames.Contains(OtherPname)) {
                            //throw new ArgumentException($"Smells like configuration error: an updater for parameter {OtherPname} was added, but none of the equations seem to need this parameter.");
                            //Console.WriteLine($"Warning: smells like configuration error: an updater for parameter {OtherPname} was added, but none of the equations seem to need this parameter.");
                        }
                    }
                }

                foreach(var OtherPname in parameterNames) {
                    bool bfound = false;
                    foreach(var p in this.parameters) {
                        if(p.ParameterNames.Contains(OtherPname))  {
                            bfound = true;
                            break;
                        }
                    }
                    if(!bfound)
                        //Console.Error.WriteLine($"Smells like configuration error: an parameter {OtherPname} is specified by the operator, but none of the updaters seems to care about this parameter.");
                        throw new ArgumentException($"Smells like configuration error: an parameter {OtherPname} is specified by the operator, but none of the updaters seems to care about this parameter.");
                }

            }

            List<ParameterS> parameters;

            public ParameterList(int capacity = 10) {
                parameters = new List<ParameterS>(capacity);
            }

            public void AddParameter(ParameterS parameter) {
                parameters.Add(parameter);
            }

            public ICollection<DelParameterFactory> Factories(IList<string> names) {
                LinkedList<string> nameList = new LinkedList<string>(names);
                LinkedList<DelParameterFactory> parameterFactories = new LinkedList<DelParameterFactory>();
                //Find parameters and remove all found parameters from list;

                while(nameList.Count > 0) {
                    string name = nameList.First.Value;
                    nameList.RemoveFirst();
                    //Find currentName
                    for(int i = 0; i < parameters.Count; ++i) {
                        ParameterS parameter = parameters[i];
                        if(parameter.ParameterNames.Contains(name)) {
                            if(parameter.Factory != null) {
                                parameterFactories.AddLast(parameter.Factory);
                            }
                            foreach(string otherParamName in parameter.ParameterNames) {
                                nameList.Remove(otherParamName);
                            }
                            break;
                        }
                    }
                }
                return parameterFactories;
            }

            public ICollection<DelPartialParameterUpdate> ParameterUpdates(IList<string> names) {
                LinkedList<string> nameList = new LinkedList<string>(names);
                LinkedList<DelPartialParameterUpdate> parameterUpdates = new LinkedList<DelPartialParameterUpdate>();

                //Find parameters and remove all found parameters from list;
                while(nameList.Count > 0) {
                    string name = nameList.First.Value;
                    nameList.RemoveFirst();
                    //Find currentName
                    for(int i = 0; i < parameters.Count; ++i) {
                        ParameterS parameter = parameters[i];
                        if(parameter.ParameterNames.Contains(name)) {
                            if(parameter.Update != null) {
                                parameterUpdates.AddLast(parameter.Update);
                            }
                            foreach(string otherParamName in parameter.ParameterNames) {
                                nameList.Remove(otherParamName);
                            }
                            break;
                        }
                    }
                }
                return parameterUpdates;
            }
        }


        ParameterList parameters;

        CoefficientsList coefficients;

        /// <summary>
        /// Default constructor 
        /// </summary>
        public OperatorFactory() {
            eqSystem = new SystemOfEquations();
            parameters = new ParameterList();
            coefficients = new CoefficientsList();
        }

        /// <summary>
        /// Add an spatial equation to this factory. 
        /// </summary>
        /// <param name="equation"></param>
        public void AddEquation(SpatialEquation equation) {
            eqSystem.AddEquation(equation);
        }

        /// <summary>
        /// Adds a surface equation to the factory.
        /// </summary>
        /// <param name="equation"></param>
        public void AddEquation(SurfaceEquation equation) {
            eqSystem.AddEquation(equation);
        }

        /// <summary>
        /// Adds a bulk equation to the factory.
        /// </summary>
        /// <param name="equation"></param>
        public void AddEquation(BulkEquation equation) {
            eqSystem.AddEquation(equation);
        }

        /// <summary>
        /// Add parameter to factory. Only parameters, that are present in at least one equation will be used.
        /// </summary>
        /// <param name="parameter"></param>
        public void AddParameter(ParameterS parameter) {
            parameters.AddParameter(parameter);
        }

        /// <summary>
        /// Add coefficient to factory. Only coefficients, that are present in at least one equation will be used.
        /// </summary>
        /// <param name="coefficient"></param>
        public void AddCoefficient(Coefficient coefficient) {
            coefficients.AddCoefficient(coefficient);
        }

        /// <summary>
        /// Creates Spatial operator
        /// </summary>
        /// <param name="quadOrder">Quadrature Order of regular cells</param>
        /// <returns>Configured spatial operator. Not committed.</returns>
        public XDifferentialOperatorMk2 GetSpatialOperator(int quadOrder) {
            int QuadOrderFunc(int[] DomvarDegs, int[] ParamDegs, int[] CodvarDegs) {
                return quadOrder;
            }
            XDifferentialOperatorMk2 XOP = GetSpatialOperator(QuadOrderFunc);
            return XOP;
        }

        /// <summary>
        /// Creates Spatial operator
        /// </summary>
        /// <param name="QuadOrderFunc">Function Mapping from Domain Variable Degrees, 
        /// Parameter Degrees and CoDomain Variable Degrees to the Quadrature Order </param>
        /// <returns>Configured spatial operator. Not committed.</returns>
        public XDifferentialOperatorMk2 GetSpatialOperator(Func<int[], int[], int[], int> QuadOrderFunc) {
            var XOP = CreateSpatialOperator(QuadOrderFunc);
            AddEquationComponents(XOP);
            AddSurfaceEquationComponents(XOP);
            AddGhostEquationComponents(XOP);
            AddTemporalOperator(XOP);
            AddParameterDelegates(XOP);
            AddCoefficients(XOP);

            this.parameters.VerifyList(XOP.ParameterVar);

            return XOP;
        }

        XDifferentialOperatorMk2 CreateSpatialOperator(Func<int[], int[], int[], int> QuadOrderFunc) {
            string[] domainVars = eqSystem.DomainVars();
            string[] codomainVars = eqSystem.CoDomainVars();
            string[] parameters = eqSystem.Parameters();
            string[] species = eqSystem.Species();

            var spatialOperator = new XDifferentialOperatorMk2(
                domainVars,
                parameters,
                codomainVars,
                QuadOrderFunc,
                species);
            return spatialOperator;
        }

        void AddTemporalOperator(XDifferentialOperatorMk2 spatialOperator) {
            (string, double[])[] diagonal = eqSystem.MassDiagonal();
            spatialOperator.TemporalOperator = new ConstantXTemporalOperator(spatialOperator, diagonal);
        }

        void AddParameterDelegates(XDifferentialOperatorMk2 spatialOperator) {
            ICollection<DelParameterFactory> factories = parameters.Factories(spatialOperator.ParameterVar);
            foreach(DelParameterFactory factory in factories) {
                spatialOperator.ParameterFactories.Add(factory);
            }
            ICollection<DelPartialParameterUpdate> updates = parameters.ParameterUpdates(spatialOperator.ParameterVar);
            foreach(DelPartialParameterUpdate update in updates) {
                spatialOperator.ParameterUpdates.Add(update);
            }
        }

        void AddEquationComponents(XDifferentialOperatorMk2 spatialOperator) {
            foreach(SpatialEquation equation in eqSystem.SpatialEquations) {
                foreach(IEquationComponent component in equation.Components) {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
            foreach (SpatialEquation equation in eqSystem.BulkEquations) {
                foreach (IEquationComponent component in equation.Components) {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
            foreach (SpatialEquation equation in eqSystem.InterfaceEquations) {
                foreach (IEquationComponent component in equation.Components) {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
        }

        void AddSurfaceEquationComponents(XDifferentialOperatorMk2 spatialOperator) {
            foreach(SurfaceEquation equation in eqSystem.InterfaceEquations) {
                if(equation.SurfaceComponents != null) {
                    foreach(IEquationComponent component in equation.SurfaceComponents) {
                        spatialOperator.SurfaceElementOperator_Ls0.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
                if(equation.ContactLineComponents != null) {
                    foreach (IEquationComponent component in equation.ContactLineComponents) {
                        spatialOperator.ContactLineOperator_Ls0.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
            }
        }

        void AddGhostEquationComponents(XDifferentialOperatorMk2 spatialOperator) {
            foreach(BulkEquation equation in eqSystem.BulkEquations) {
                if(equation.GhostComponents != null) {
                    foreach(IEquationComponent component in equation.GhostComponents) {
                        spatialOperator.GhostEdgesOperator.EquationComponents[equation.CodomainName].Add(component);
                    }
                }
            }
        }

        void AddCoefficients(XDifferentialOperatorMk2 spatialOperator) {
            spatialOperator.OperatorCoefficientsProvider = delegate (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
                var r = Coefficients(lstrk, spc, quadOrder, TrackerHistoryIdx, time);
                r.HomotopyValue = spatialOperator.CurrentHomotopyValue;
                return r;
            };
        }

        CoefficientSet Coefficients(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            var r = new CoefficientSet() {
                GrdDat = lstrk.GridDat
            };
            var g = lstrk.GridDat;
            if(g is Foundation.Grid.Classic.GridData cgdat) {
                r.CellLengthScales = cgdat.Cells.CellLengthScale;
                r.EdgeLengthScales = cgdat.Edges.h_min_Edge;

            } else {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + g.GetType().FullName);
            }

            r.SpeciesSubGrdMask = lstrk.Regions.GetSpeciesSubGrid(spc).VolumeMask.GetBitMaskWithExternal();

            string[] coeffs = eqSystem.Coefficients();

            ICollection<DelCoefficientFactory> factories = coefficients.Factories(coeffs);
            foreach(DelCoefficientFactory factory in factories) {
                var ttt = factory(lstrk, spc, quadOrder, TrackerHistoryIdx, time);
                foreach(var tt in ttt)
                    r.UserDefinedValues[tt.CoefficientName] = tt.CoefficientValue;
            }

            string[] actualcoeffs = r.UserDefinedValues.Keys.ToArray();
            foreach(string coeff in coeffs)
                if(Array.IndexOf(actualcoeffs, coeff) < 0) throw new ApplicationException("configuration error in spatial differential operator; some equation component depends on coefficient  \""
                                    + coeff
                                    + "\", but this name is not a member of the UserDefinedValues list.");            
            return r;
        }
    }
}
