using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using Microsoft.SqlServer.Server;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    /// <summary>
    /// Factory for the allocation of storage for storing the coefficients for this component
    /// </summary>
    /// <returns>
    /// a pair, containing:
    /// - coefficient name
    /// - an object to store the respective coefficient value(s)
    /// </returns>
    public delegate (string CoefficientName, object CoefficientValue)[] DelCoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time);

    /// <summary>
    /// For now <see cref="DelCoefficientFactory"/> is called repeatedly,
    /// maybe this will be changed in the future to align coefficient handling to Parameterhandling
    /// </summary>
    //public delegate void DelCoefficientUpdate(IReadOnlyDictionary<string, object> Coefficients, double time);

    abstract class Coefficient
    {
        public abstract IList<string> CoefficientsNames { get;}

        public abstract DelCoefficientFactory Factory { get; }

        //public DelCoefficientUpdate Update;
    }

    class CoefficientsList {

        List<Coefficient> coefficients;

        public CoefficientsList(int capacity = 10)
        {
            coefficients = new List<Coefficient>(capacity);
        }

        public void AddCoefficient(Coefficient coefficient)
        {
            coefficients.Add(coefficient);
        }

        public ICollection<DelCoefficientFactory> Factories(IList<string> names) {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelCoefficientFactory> coefficientFactories = new LinkedList<DelCoefficientFactory>();
            //Find parameters and remove all found parameters from list;

            while (nameList.Count > 0) {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for (int i = 0; i < coefficients.Count; ++i) {
                    Coefficient coefficient = coefficients[i];
                    if (coefficient.CoefficientsNames.Contains(name)) {
                        if (coefficient.Factory != null) {
                            coefficientFactories.AddLast(coefficient.Factory);
                        }
                        foreach (string otherCoefficientName in coefficient.CoefficientsNames) {
                            nameList.Remove(otherCoefficientName);
                        }
                        break;
                    }
                }
            }
            return coefficientFactories;
        }

        //public ICollection<DelCoefficientUpdate> CoefficientUpdates(IList<string> names) {
        //    LinkedList<string> nameList = new LinkedList<string>(names);
        //    LinkedList<DelCoefficientUpdate> coefficientUpdates = new LinkedList<DelCoefficientUpdate>();

        //    //Find parameters and remove all found parameters from list;
        //    while (nameList.Count > 0) {
        //        string name = nameList.First.Value;
        //        nameList.RemoveFirst();
        //        //Find currentName
        //        for (int i = 0; i < coefficients.Count; ++i) {
        //            Coefficient coefficient = coefficients[i];
        //            if (coefficient.CoefficientsNames.Contains(name)) {
        //                if (coefficient.Update != null) {
        //                    coefficientUpdates.AddLast(coefficient.Update);
        //                }
        //                foreach (string otherCoefficientName in coefficient.CoefficientsNames) {
        //                    nameList.Remove(otherCoefficientName);
        //                }
        //                break;
        //            }
        //        }
        //    }
        //    return coefficientUpdates;
        //}

    }
}
