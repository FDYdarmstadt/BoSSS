/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using ilPSP;
using Mono.CSharp;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Runtime.Serialization;
using System.Linq;

namespace BoSSS.Solution.Control {

    /// <summary>
    /// This class encapsulates the representation of a mathematical formula (dependent on space and time)
    /// \f[
    ///   (\vec{x},t) \mapsto f(\vec{x},t)
    ///   \ \text{ resp. } \
    ///   (\vec{x}) \mapsto f(\vec{x})
    /// \f]
    /// which is used to provide boundary or initial values.
    /// The mathematical expression is compiled from C#-code on the fly.
    /// In contrast to delegates, this is representation of mathematical formulas is serializeable.
    /// </summary>
    [Serializable]
    [DataContract]
    public class Formula : IBoundaryAndInitialData {

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="code">
        /// C#-code which represents the mathematical formula,
        ///  - if <paramref name="TimeDep"/> is true, a `Func<double[], double, double>`, representing an expression \f$ (\vec{x},t) \mapsto f(\vec{x},t) \f$
        ///  - if <paramref name="TimeDep"/> is false, a `Func<double[], double>`, representing an expression \f$ (\vec{x}) \mapsto f(\vec{x}) \f$
        /// </param>
        /// <param name="AdditionalPrefixCode">
        /// Optional, additional C#-statements, e.g. auxiliary definitions, which is entered before <paramref name="code"/>.
        /// </param>
        /// <param name="TimeDep">
        /// Whether the function is time dependent or no, see <paramref name="TimeDep"/>.
        /// </param>
        public Formula(string code, bool TimeDep, string AdditionalPrefixCode = "") {
            m_Code = code;
            m_TimeDep = TimeDep;
            m_AdditionalPrefixCode = AdditionalPrefixCode;
            Compile();
        }

        [DataMember]
        bool m_TimeDep;

        [DataMember]
        string m_Code;

        [DataMember]
        string m_AdditionalPrefixCode;

        [NonSerialized]
        Func<double[], double, double> m_Xt_Del;

        [NonSerialized]
        Func<double[], double> m_X__Del;



        void Compile() {
            if (m_Xt_Del == null && m_X__Del == null) {
                using (var err = new StringWriter()) {
                    var Settings = new CompilerSettings();
#if DEBUG
                    Settings.Optimize = false;
#else
                    Settings.Optimize = false;
#endif
                    var Printer = new StreamReportPrinter(err);


                    CompilerContext cmpCont = new CompilerContext(Settings, Printer);


                    Evaluator eval = new Evaluator(cmpCont);
                    eval.InteractiveBaseClass = typeof(Object);

                    Assembly[] allAssis = BoSSS.Solution.Application.GetAllAssemblies().ToArray();

                    foreach (var assi in allAssis) {
                        eval.ReferenceAssembly(assi);
                    }

                    eval.Compile(@"
                                using System;
                                using System.Collections.Generic;
                                using System.Linq;
                                using ilPSP;
                                using ilPSP.Utils;
                                using BoSSS.Platform;
                                using BoSSS.Platform.Utils;
                                using BoSSS.Foundation;
                                using BoSSS.Foundation.Grid;
                                using BoSSS.Foundation.IO;
                                using BoSSS.Solution;
                                using BoSSS.Solution.Utils;
                      ");

                    if (!m_AdditionalPrefixCode.IsEmptyOrWhite()) {
                        try {
                            object figdi;
                            bool dummy;
                            eval.Evaluate(m_AdditionalPrefixCode, out figdi, out dummy);
                        } catch (Exception e) {
                            throw new AggregateException(e.GetType().Name + " during the interpretation of code snippet '"
                                + m_AdditionalPrefixCode + "'" + err.NewLine + "Error(s): " + err.NewLine + err.ToString(),
                                e);
                        }
                    }

                    object formula;
                    try {
                        string Prefix = m_TimeDep ? "Func<double[], double, double>" : "Func<double[], double>";
                        Prefix = Prefix + " myfunc = ";
                        object result;
                        bool result_set;
                        string ans = eval.Evaluate(Prefix + m_Code + ";", out result, out result_set);
                        formula = eval.Evaluate("myfunc;");
                    } catch (Exception e) {
                        throw new AggregateException(e.GetType().Name + " during the interpretation of code snippet '"
                            + m_Code + "'" + err.NewLine + "Error(s): " + err.NewLine + err.ToString(),
                            e);
                    }

                    if (formula != null && cmpCont.Report.Errors == 0) {
                        if (formula is Func<double[], double, double>) {
                            m_Xt_Del = (Func<double[], double, double>)formula;
                            return;
                        }

                        if (formula is Func<double[], double>) {
                            m_X__Del = (Func<double[], double>)formula;
                            return;
                        }
                    }

                    throw new ArgumentException("Unable to cast result of code snippet '" + m_Code + " to a valid expression (Func<double[],double,double> or Func<double[],double>)." + err.NewLine + "Error(s): " + err.NewLine + err.ToString());
                }
            }
        }

        /// <summary>
        /// Evaluates the function.
        /// </summary>
        /// <param name="X">global/physical coordinates.</param>
        /// <param name="t">Physical time.</param>
        /// <returns>Function value.</returns>
        public double Evaluate(double[] X, double t) {
            Compile();
            Debug.Assert((m_Xt_Del == null) != (m_X__Del == null));

            if (m_Xt_Del != null) {
                return m_Xt_Del(X, t);
            } else if (m_X__Del != null) {
                return m_X__Del(X);
            } else {
                throw new ApplicationException();
            }
        }

    }
}
