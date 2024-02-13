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
using ilPSP.Utils;

namespace BoSSS.Foundation.XDG {

    partial class XDGField {
        /// <summary>
        /// Projects a function onto an XDG Field, using all species
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <param name="U"></param>
        virtual public void ProjectFunctionXDG(double alpha, Func<Vector, double[], int, double> f, params XDGField[] U)
        {
            ProjectFunctionXDG(alpha, f, this.Basis.Tracker.SpeciesIdS, U);
        }
        /// <summary>
        /// Projects a function onto an XDG Field, only on some prescribed species
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <param name="speciesToEvaluateIDs"></param>
        /// <param name="U"></param>
        virtual public void ProjectFunctionXDG(double alpha, Func<Vector, double[], int, double> f,  IList<SpeciesId> speciesToEvaluateIDs, params XDGField[] U) {
            var LsTrk = this.Basis.Tracker;
            IList<string> speciesToEvaluateNames = new List<string>(); 
            foreach (SpeciesId spcID in speciesToEvaluateIDs)
            {
                speciesToEvaluateNames.Add(LsTrk.GetSpeciesName(spcID));
            }
            string[] Dom = new string[U.Length];
            for (int i = 0; i < Dom.Length; i++)
                Dom[i] = "_" + i;

            string[] Cod = new string[] { "res" };

            XDifferentialOperatorMk2 src = new XDifferentialOperatorMk2(0.1, speciesToEvaluateNames.ToArray<string>());

            foreach (string spc in speciesToEvaluateNames)
            {
                src.EquationComponents[Cod[0]].Add(new ProjectFunctionSource(spc, f, Dom));
            }
            src.FluxesAreNOTMultithreadSafe = false;
            src.Commit();

            var ev = src.GetEvaluatorEx(
                new CoordinateMapping(U), null, this.Mapping);

            ev.Evaluate(alpha, 1.0, this.CoordinateVector);

            //Next, in the cut cells, we need to multiply the CoordinateVector with the Inverse Mass Matrix
            
            int deg = this.Basis.Degree;
            int N = this.CoordinateVector.Length;
            double[] tmp = new double[N];

            //chose a quadOrder
            var AvailOrders = LsTrk.GetCachedOrders().Where(order => order >= 2 * deg);
            int order2Pick = AvailOrders.Any() ? AvailOrders.Min() : 2 * deg;

            var MMF = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, order2Pick).MassMatrixFactory;
            var MM = MMF.GetMassMatrix(this.Mapping, inverse: true);
            
            MM.SpMV(1.0, this.CoordinateVector, 0.0, tmp);

            this.CoordinateVector.SetV(tmp);
            ////helper Function
            //double[] GetCoords(int j, SpeciesId spc, XDGField field)
            //{
            //    int iSpc = LsTrk.Regions.GetSpeciesIndex(spc, j);
            //    if (iSpc < 0)
            //        return null;
            //    else
            //        return field.Coordinates.GetRowPart(j, iSpc * N, N);
            //}
            //// loop over every species
            //foreach (SpeciesId speciesId in LsTrk.SpeciesIdS)
            //{
            //    var speciesMask = LsTrk.Regions.GetSpeciesMask(speciesId);
            //    var cutCellMask = LsTrk.Regions.GetCutCellMask().Intersect(speciesMask);

            //    var MMF = LsTrk.GetXDGSpaceMetrics(speciesId, order2Pick).MassMatrixFactory;
            //    var MMblox = MMF.GetMassMatrixBlocks(this.Basis.NonX_Basis, speciesId);

            //    double[] tmp = new double[N];
            //    for (int iSub = 0; iSub < MMblox.jSub2jCell.Length; iSub++) // loop over every cut cell
            //    {
            //        int jCell = MMblox.jSub2jCell[iSub];
            //        double[] CoordsFTT = GetCoords(jCell, speciesId, this);
            //        if (CoordsFTT == null)
            //            continue; // species not present in cell; no contribution.

            //        var MM_j = MMblox.MassMatrixBlocks.ExtractSubArrayShallow(new int[] { iSub, 0, 0 }, new int[] { iSub - 1, N - 1, N - 1 });

            //        MM_j.GEMV(1.0, CoordsFTT, 0.0, tmp);
            //        double denominator = CoordsFTT.InnerProd(tmp);


            //        //if(LsTrk.GetSpeciesName(speciesId) == "A") {
            //        //    PlotCurrentState(CurrentStepNo + "0000" + 2*jCell);
            //        //} else {
            //        //    PlotCurrentState(CurrentStepNo + "0000" + 2 * jCell+1);
            //        //}
            //    }
            //}
        }

        class ProjectFunctionSource : ISpeciesFilter, IVolumeForm {
            string species;
            Func<Vector, double[], int, double> f;
            string[] arguments;


            public TermActivationFlags VolTerms => TermActivationFlags.UxV;

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