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
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Connectors.Matlab;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.Grid;

namespace BoSSS.Application.SipPoisson {

    /// <summary>
    /// predefined control-objects
    /// </summary>
    static public class SipHardcodedControl
    {
        /// <summary>
        /// Poisson Equation on a (-1,1)x(-1,1), Dirichlet everywhere
        /// </summary>
        public static SipControl RegularSquare(int xRes = 5, int yRes = 5, int deg = 2) {

            Func<double[], double> exSol = 
                X => -Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5);
            Func<double[], double> exRhs = 
                X => (Math.PI * Math.PI * 0.5 * Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5)); // == - /\ exSol

            var R = new SipControl();
            R.ProjectName = "ipPoison/square";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 4 });
            R.InitialValues_Evaluators.Add("RHS", exRhs);
            R.InitialValues_Evaluators.Add("Tex", exSol);
            R.ExactSolution_provided = true;
            R.SuppressExceptionPrompt = true;
            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T", exSol);
            R.NoOfSolverRuns = 1;

            R.AdaptiveMeshRefinement = true;
            R.NoOfTimesteps = 1;

            R.ImmediatePlotPeriod = 1;
            R.SuperSampling = 2;

            R.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-1, 1, xRes);
                double[] yNodes = GenericBlas.Linspace(-1, 1, yRes);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret = 1;
                    return ret;
                });
                return grd;
            };
            return R;
        }

        public static SipControl VoronoiSquare(int res = 50, int deg = 2)
        {
            var R = new SipControl
            {
                ProjectName = "SipPoisson-Voronoi",
                SessionName = "testrun",
                ImmediatePlotPeriod = 1,
                SuperSampling = 2,
                savetodb = false,
                ExactSolution_provided = false,

            };
            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg * 2 });
            R.InitialValues_Evaluators.Add("RHS", X => X[0] * X[0]);
            R.InitialValues_Evaluators.Add("Tex", X => X[0]);

            Func<double[], double> dirichletBoundary =
                X => 0.0;
            //  X => Math.Pow(X[0], 2) + Math.Pow(X[1], 2);
            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T", dirichletBoundary);

            R.GridFunc = delegate (){
                Vector[] DomainBndyPolygon = new[] {
                    new Vector(-1,1),
                    new Vector(1,1),
                    new Vector(1,-1),
                    new Vector(-1,-1)
                };
                AggregationGrid grid;
                grid = BoSSS.Foundation.Grid.Voronoi.VoronoiGrid2D.Polygonal(DomainBndyPolygon, 10, res);
                grid.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grid.DefineEdgeTags( (double[] X) => (byte)1);
                return grid;
            };
            return R;
        }

        public static SipControl JumpingSquare(int xRes = 5, int yRes = 5, int deg = 3)
        {
            var R = new SipControl();
            R.ProjectName = "ipPoison/square";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 4 });
            R.InitialValues_Evaluators.Add("RHS", X => 0.0);
            R.InitialValues_Evaluators.Add("Tex", X => 0.1);
            R.ExactSolution_provided = true;
            R.SuppressExceptionPrompt = true;
            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T", 
                X => { 
                    if (X[0] > 0.9999) 
                    {
                        return 0.1;
                    } 
                    else 
                    { 
                        return 0.1; 
                    }
                });
            R.NoOfSolverRuns = 1;

            R.AdaptiveMeshRefinement = true;
            R.NoOfTimesteps = 1;

            R.ImmediatePlotPeriod = 1;
            R.SuperSampling = 2;

            R.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-1, 1, xRes);
                double[] yNodes = GenericBlas.Linspace(-1, 1, yRes);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret = 1;
                    return ret;
                });
                return grd;
            };
            return R;
        }
    }
}
