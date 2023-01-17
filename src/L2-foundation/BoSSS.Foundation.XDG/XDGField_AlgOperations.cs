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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;

namespace BoSSS.Foundation.XDG {

    partial class XDGField {

        virtual public void ProjectFunction(double alpha, Func<Vector, double[], int, double> f, params XDGField[] U) {

            IList<string> species = U[0].Basis.Tracker.SpeciesNames;
            string[] Dom = new string[U.Length];
            for (int i = 0; i < Dom.Length; i++)
                Dom[i] = "_" + i;

            string[] Cod = new string[] { "res" };

            XSpatialOperatorMk2 src = new XSpatialOperatorMk2(0.1, species.ToArray<string>());
            src.EquationComponents[Cod[0]].Add(new ProjectFunctionSource("A", f, Dom));
            src.EquationComponents[Cod[0]].Add(new ProjectFunctionSource("B", f, Dom));
            src.Commit();

            var ev = src.GetEvaluatorEx(
                new CoordinateMapping(U), null, this.Mapping);

            ev.Evaluate(alpha, 1.0, this.CoordinateVector);
        }

        class ProjectFunctionSource : ISpeciesFilter, IVolumeForm {
            string species;
            Func<Vector, double[], int, double> f;
            string[] arguments;


            public TermActivationFlags VolTerms => TermActivationFlags.V;

            public IList<string> ArgumentOrdering => arguments;

            public IList<string> ParameterOrdering => null;

            public string ValidSpecies => species;

            public ProjectFunctionSource(string species, Func<Vector, double[], int, double> f,string[] Dom) {
                this.species = species;
                this.f = f;
                arguments = Dom;
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                return f(cpv.Xglobal, U, cpv.jCell) * V;
            }
        }
    }
}