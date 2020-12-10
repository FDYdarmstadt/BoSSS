using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Foundation.XDG.OperatorFactory {

    public class OperatorFactory
    {
        SystemOfEquations eqSystem;

        ParameterList parameters;

        CoefficientsList coefficients;

        public OperatorFactory()
        {
            eqSystem = new SystemOfEquations();
            parameters = new ParameterList();
            coefficients = new CoefficientsList();
        }

        public void AddEquation(SpatialEquation equation)
        {
            eqSystem.AddEquation(equation);
        }

        public void AddEquation(SurfaceEquation equation)
        {
            eqSystem.AddEquation(equation);
        }

        public void AddEquation(BulkEquation equation)
        {
            eqSystem.AddEquation(equation);
        }

        public void AddParameter(Parameter parameter)
        {
            parameters.AddParameter(parameter);
        }        

        public void AddCoefficient(Coefficient coefficient)
        {
            coefficients.AddCoefficient(coefficient);
        }

        public XSpatialOperatorMk2 GetSpatialOperator(int quadOrder)
        {
            int QuadOrderFunc(int[] DomvarDegs, int[] ParamDegs, int[] CodvarDegs)
            {
                return quadOrder;
            }
            XSpatialOperatorMk2 XOP =  GetSpatialOperator(QuadOrderFunc);
            return XOP;
        }

        public XSpatialOperatorMk2 GetSpatialOperator(Func<int[], int[], int[], int> QuadOrderFunc)
        {
            var XOP = CreateSpatialOperator(QuadOrderFunc);
            AddEquationComponents(XOP);
            AddSurfaceEquationComponents(XOP);
            AddGhostEquationComponents(XOP);
            AddTemporalOperator(XOP);
            AddParameterDelegates(XOP);
            AddCoefficients(XOP);
            return XOP;
        }

        XSpatialOperatorMk2 CreateSpatialOperator(Func<int[], int[], int[], int> QuadOrderFunc)
        {
            string[] domainVars = eqSystem.DomainVars();
            string[] codomainVars = eqSystem.CoDomainVars();
            string[] parameters = eqSystem.Parameters();
            string[] species = eqSystem.Species();

            var spatialOperator = new XSpatialOperatorMk2(
                domainVars,
                parameters,
                codomainVars,
                QuadOrderFunc,
                species);
            return spatialOperator;
        }

        void AddTemporalOperator(XSpatialOperatorMk2 spatialOperator)
        {
            (string, double[])[] diagonal = eqSystem.MassDiagonal();
            spatialOperator.TemporalOperator = new ConstantXTemporalOperator(spatialOperator, diagonal);
        }

        void AddParameterDelegates(XSpatialOperatorMk2 spatialOperator)
        {
            ICollection<DelParameterFactory> factories = parameters.Factories(spatialOperator.ParameterVar);
            foreach(DelParameterFactory factory in factories)
            {
                spatialOperator.ParameterFactories.Add(factory);
            }
            ICollection<DelPartialParameterUpdate> updates = parameters.ParameterUpdates(spatialOperator.ParameterVar);
            foreach (DelPartialParameterUpdate update in updates)
            {
                spatialOperator.ParameterUpdates.Add(update);
            }
        }

        void AddEquationComponents(XSpatialOperatorMk2 spatialOperator)
        {
            foreach (SpatialEquation equation in eqSystem.SpatialEquations)
            {
                foreach (IEquationComponent component in equation.Components)
                {
                    spatialOperator.EquationComponents[equation.CodomainName].Add(component);
                }
            }
        }

        void AddSurfaceEquationComponents(XSpatialOperatorMk2 spatialOperator)
        {
            foreach (SurfaceEquation equation in eqSystem.InterfaceEquations)
            {
                if (equation.SurfaceComponents != null)
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
            foreach (BulkEquation equation in eqSystem.BulkEquations)
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

        void AddCoefficients(XSpatialOperatorMk2 spatialOperator)
        {
            spatialOperator.OperatorCoefficientsProvider = Coefficients;
        }

        CoefficientSet Coefficients(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time)
        {
            var r = new CoefficientSet()
            {
                GrdDat = lstrk.GridDat
            };
            var g = lstrk.GridDat;
            if (g is Foundation.Grid.Classic.GridData cgdat)
            {
                r.CellLengthScales = cgdat.Cells.CellLengthScale;
                r.EdgeLengthScales = cgdat.Edges.h_min_Edge;

            }
            else
            {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + g.GetType().FullName);
            }

            string[] coeffs = eqSystem.Coefficients();

            ICollection<DelCoefficientFactory> factories = coefficients.Factories(coeffs);            
            foreach (DelCoefficientFactory factory in factories) {
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
