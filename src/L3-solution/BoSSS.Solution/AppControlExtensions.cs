using BoSSS.Foundation;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;


namespace BoSSS.Solution.Control {
    
    
    /// <summary>
    /// 
    /// </summary>
    public static class AppControlExtensions {

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_Evaluators"/> or <see cref="AppControl.ExactSolutions_Evaluators"/>
        /// </summary>
        static public Func<double[], double> GetValueOrDefault(this IDictionary<string, Func<double[], double>> collection, string FieldName) {
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return (_ => 0.0);
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_Evaluators"/> or <see cref="AppControl.ExactSolutions_Evaluators"/>, Multiphase version
        /// </summary>
        static public Func<double[], double> GetValueOrDefault(this IDictionary<string, Func<double[], double>> collection, string FieldName, string SpeciesName) {
            if ( collection.TryGetValue(FieldName + "#" + SpeciesName, out var result) )
                return result;
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return (_ => 0.0);
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_Evaluators_TimeDep"/> or <see cref="AppControl.ExactSolutions_Evaluators_TimeDep"/> or <see cref="AppControl.BoundaryValueCollection.Evaluators"/>
        /// </summary>
        static public Func<double[], double, double> GetValueOrDefault(this IDictionary<string, Func<double[], double, double>> collection, string FieldName) {
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return ((X, t) => 0.0);
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_Evaluators_TimeDep"/> or <see cref="AppControl.ExactSolutions_Evaluators_TimeDep"/> or <see cref="AppControl.BoundaryValueCollection.Evaluators"/>, Multiphase version
        /// </summary>
        static public Func<double[], double, double> GetValueOrDefault(this IDictionary<string, Func<double[], double, double>> collection, string FieldName, string SpeciesName) {
            if ( collection.TryGetValue(FieldName + "#" + SpeciesName, out var result) )
                return result;
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return ((X, t) => 0.0);
        }

        static void DefaultScalarFunctionTimeDep(MultidimensionalArray input, double time, MultidimensionalArray output) {
            output.Clear();
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_EvaluatorsVec"/> or <see cref="AppControl.InitialValues_EvaluatorsVec"/>
        /// </summary>
        static public ScalarFunctionTimeDep GetValueOrDefault(this IDictionary<string, ScalarFunctionTimeDep> collection, string FieldName) {
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return DefaultScalarFunctionTimeDep;
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_EvaluatorsVec"/> or <see cref="AppControl.InitialValues_EvaluatorsVec"/>, Multiphase version
        /// </summary>
        static public ScalarFunctionTimeDep GetValueOrDefault(this IDictionary<string, ScalarFunctionTimeDep> collection, string FieldName, string SpeciesName) {
            if ( collection.TryGetValue(FieldName + "#" + SpeciesName, out var result) )
                return result;
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return DefaultScalarFunctionTimeDep;
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues"/> or <see cref="AppControl.ExactSolutions"/> or <see cref="AppControl.BoundaryValueCollection.Values"/>
        /// </summary>
        static public IBoundaryAndInitialData GetValueOrDefault(this IDictionary<string, IBoundaryAndInitialData> collection, string FieldName) {
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return new ConstantValue(0.0);
        }

        /// <summary>
        /// Tries to obtain an value from <see cref="AppControl.InitialValues_Evaluators"/> or <see cref="AppControl.ExactSolutions_Evaluators"/>, Multiphase version
        /// </summary>
        static public IBoundaryAndInitialData GetValueOrDefault(this IDictionary<string, IBoundaryAndInitialData> collection, string FieldName, string SpeciesName) {
            if ( collection.TryGetValue(FieldName + "#" + SpeciesName, out var result) )
                return result;
            if ( collection.TryGetValue(FieldName, out var result2) )
                return result2;
            return new ConstantValue(0.0);
        }


    }
}
