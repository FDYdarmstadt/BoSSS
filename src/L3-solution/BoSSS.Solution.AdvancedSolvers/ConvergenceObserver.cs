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
using ilPSP.LinSolvers;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Platform;
using MPI.Wrappers;
using BoSSS.Solution.Gnuplot;
using System.IO;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// A helper class to analyze and visualize the convergence of an iterative solver algorithm.
    /// - writes tecplot files containing the residuals, error (if the solution is known), <see cref="TecplotOut"/>.
    /// - records convergence trends for different grid and \f$ p \f$-levels.
    /// </summary>
    public class ConvergenceObserver {

        public static double last_resReduction;
        public static double last_errReduction;

        /// <summary>
        /// Performs a mode decay analysis (<see cref="Waterfall(bool, bool, int)"/>) on this solver.
        /// </summary>
        public static Dictionary<string,Plot2Ddata> WaterfallAnalysis(ISolverWithCallback linearSolver, MultigridOperator mgOperator, BlockMsrMatrix MassMatrix) {
            int L = mgOperator.BaseGridProblemMapping.LocalLength;

            Console.WriteLine($"Waterfall analysis for {mgOperator.GridData.CellPartitioning.TotalLength} cells, {L} DOFs...");

            var RHS = new double[L];
            var exSol = new double[L];

            ConvergenceObserver co = new ConvergenceObserver(mgOperator, MassMatrix, exSol);
            var bkup = linearSolver.IterationCallback;
            linearSolver.IterationCallback = co.IterationCallback;

            var pc = linearSolver as IProgrammableTermination;
            var termBkup = pc?.TerminationCriterion;
            //if (pc != null) {
            //    pc.TerminationCriterion = (iter, R0_l2, R_l2) => {
            //        return (R_l2 / R0_l2 > 1e-8) && (iter < 1500);
            //    };
            //}

            // use a random init for intial guess; if RHS = 0, then the exact solution is 0, so the error for approximate solution x_i in iteration i is (x_i - 0).
            double[] x0 = GenericBlas.RandomVec(L, 0);

            // execute solver 
            linearSolver.Init(mgOperator);
            mgOperator.UseSolver(linearSolver, x0, RHS);

            // reset and return
            linearSolver.IterationCallback = bkup;
            if(termBkup != null) {
                pc.TerminationCriterion = termBkup;
            }
            //var p = co.PlotIterationTrend(true, false, true, true);

            double resReduction = Math.Exp(Math.Log(co.Iter0ResidualNorm / co.LastResidualNorm) / co.NumberOfIterations);
            double errReduction = Math.Exp(Math.Log(co.Iter0SolNorm / co.LastSolNorm) / co.NumberOfIterations);
            last_resReduction = resReduction;
            last_errReduction = errReduction;

            Console.WriteLine("Residual reduction:   " + Math.Pow(co.LastResidualNorm / co.Iter0ResidualNorm, -1));
            Console.WriteLine("Error    reduction:   " + Math.Pow(co.LastSolNorm / co.Iter0SolNorm, -1));
            Console.WriteLine("Residual reduction per iteration:   " + resReduction);
            Console.WriteLine("Error    reduction per iteration:   " + errReduction);
            
            
            var p = co.Waterfall(true, true, 100);
            return p;
        }



        /// <summary>
        /// ctor
        /// </summary>
        public ConvergenceObserver(MultigridOperator muop, BlockMsrMatrix MassMatrix, double[] __ExactSolution) {
            Setup(muop, MassMatrix, __ExactSolution);
        }

        /// <summary>
        /// another constructor
        /// </summary>
        public ConvergenceObserver(MultigridOperator muop, BlockMsrMatrix MassMatrix, double[] __ExactSolution, ISolverFactory SF) {
            m_SF = SF;
            Setup(muop, MassMatrix, __ExactSolution);
        }

        private void Setup(MultigridOperator muop, BlockMsrMatrix MassMatrix, double[] __ExactSolution) {
            if (__ExactSolution != null) {
                if (__ExactSolution.Length != muop.BaseGridProblemMapping.LocalLength)
                    throw new ArgumentException();
            }
            this.SolverOperator = muop;

            List<AggregationGridBasis[]> aggBasisSeq = new List<AggregationGridBasis[]>();
            for (var mo = muop; mo != null; mo = mo.CoarserLevel) {
                aggBasisSeq.Add(mo.Mapping.AggBasis);
            }

            this.ExactSolution = __ExactSolution;
            int[] Degrees = muop.BaseGridProblemMapping.BasisS.Select(b => b.Degree).ToArray();

            BlockMsrMatrix DummyOpMatrix = new BlockMsrMatrix(muop.BaseGridProblemMapping, muop.BaseGridProblemMapping);
            DummyOpMatrix.AccEyeSp(123);


            MultigridOperator.ChangeOfBasisConfig[][] config = new MultigridOperator.ChangeOfBasisConfig[1][];
            config[0] = new MultigridOperator.ChangeOfBasisConfig[muop.BaseGridProblemMapping.BasisS.Count];
            for (int iVar = 0; iVar < config[0].Length; iVar++) {
                config[0][iVar] = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { Degrees[iVar] },
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { iVar }
                };
            }

            //this.DecompositionOperator = muop; this.DecompositionOperator_IsOrthonormal = false;
            this.DecompositionOperator = new MultigridOperator(aggBasisSeq, muop.BaseGridProblemMapping, DummyOpMatrix, MassMatrix, config, muop.AbstractOperator);
            this.DecompositionOperator_IsOrthonormal = true;

            ResNormTrend = new Dictionary<(int, int, int), List<double>>();
            ErrNormTrend = new Dictionary<(int, int, int), List<double>>();
            for (var mgop = this.DecompositionOperator; mgop != null; mgop = mgop.CoarserLevel) {
                int[] _Degrees = mgop.Mapping.DgDegree;
                for (int iVar = 0; iVar < _Degrees.Length; iVar++) {
                    for (int p = 0; p <= _Degrees[iVar]; p++) {
                        ResNormTrend.Add((mgop.LevelIndex, iVar, p), new List<double>());
                        ErrNormTrend.Add((mgop.LevelIndex, iVar, p), new List<double>());
                    }
                }
            }
        }

        bool DecompositionOperator_IsOrthonormal {
            get;
            set;
        }

        private ISolverFactory m_SF;

        /// <summary>
        /// Used to compute an orthonormal decomposition 
        /// </summary>
        public MultigridOperator DecompositionOperator {
            get;
            private set;
        }
                
        
        MultigridOperator SolverOperator;

        double[] ExactSolution;

        /// <summary>
        /// L2-norm of the residual per variable, multigrid level and DG polynomial degree;
        /// - 1st index: multigrid level index
        /// - 2nd index: variable index
        /// - 3rd index: polynomial degree
        /// </summary>
        Dictionary<(int MGlevel, int iVar, int deg), List<double>> ResNormTrend;

        /// <summary>
        /// L2-norm of the error per variable, multigrid level and DG polynomial degree;
        /// - 1st index: multigrid level index
        /// - 2nd index: variable index
        /// - 3rd index: polynomial degree 
        /// </summary>
        Dictionary<(int MGlevel, int iVar, int deg), List<double>> ErrNormTrend;


        /// <summary>
        /// Visualization of data from <see cref="WriteTrendToTable"/> in gnuplot
        /// </summary>
        public Plot2Ddata PlotIterationTrend(bool ErrorOrResidual, bool SepVars, bool SepPoly, bool SepLev) {
            var Ret = new Plot2Ddata();

            Ret.Title = ErrorOrResidual ? "\"Error trend\"" : "\"Residual trend\"";
            Ret.LogY = true;

            string[] Titels;
            MultidimensionalArray ConvTrendData;
            WriteTrendToTable(ErrorOrResidual, SepVars, SepPoly, SepLev, out Titels, out ConvTrendData);
            double[] IterNo = ConvTrendData.GetLength(0).ForLoop(i => ((double)i));

            for(int iCol = 0; iCol < Titels.Length; iCol++) {
                var g = new Plot2Ddata.XYvalues(Titels[iCol], IterNo, ConvTrendData.GetColumn(iCol));
                g.Format = new PlotFormat(lineColor: ((LineColors)(iCol + 1)), Style: Styles.Lines);

            }

            return Ret;
        }

        /// <summary>
        /// Visualization of error data from <see cref="WriteTrendToTable"/> in gnuplot
        /// </summary>
        public Dictionary<string, Plot2Ddata> Waterfall(bool IncludeError, bool IncludeResidual, int NoOfIter = int.MaxValue) {


            var RetDict = new Dictionary<string, Plot2Ddata>();

            //string[] Titels;
            //MultidimensionalArray ConvTrendData;
            //WriteTrendToTable(ErrorOrResidual, SepVars, SepPoly, SepLev, out var Titels, out ConvTrendData);
            //var AllData = ErrorOrResidual ? this.ErrNormTrend : this.ResNormTrend;

            int i = 0;
            foreach (var AllData in new[] { this.ErrNormTrend, this.ResNormTrend }) {
                i++;

                int MglMax = AllData.Keys.Max(tt => tt.MGlevel);
                int VarMax = AllData.Keys.Max(tt => tt.iVar);
                int MaxIter = AllData.First().Value.Count - 1;

                for (int iVar = 0; iVar < VarMax + 1; iVar++) {
                    if (!RetDict.TryGetValue("Waterfall_of_Var" + iVar, out Plot2Ddata Ret)) {
                        Ret = new Plot2Ddata();
                        RetDict.Add("Waterfall_of_Var" + iVar, Ret);
                    }
                    Ret.Title = "Waterfall of Variable " + iVar;

                    int DegMax = AllData.Keys.Where(v => v.iVar == iVar).Max(tt => tt.deg); // pressure may have DG-1
                    var WaterfallData = new List<double[]>();
                    var Row = new List<double>();
                    var xCoords = new List<double>();

                    for (int iIter = 0; iIter <= Math.Min(NoOfIter, MaxIter); iIter++) {
                        Row.Clear();
                        for (int iLv = MglMax; iLv >= 0; iLv--) {
                            for (int p = 0; p <= DegMax; p++) {
                                double Acc = 0;

                                foreach (var kv in AllData) {
                                    if (kv.Key.deg == p && kv.Key.MGlevel == iLv && kv.Key.iVar == iVar)
                                        Acc += kv.Value[iIter].Pow2();
                                }

                                Acc = Acc.Sqrt();
                                Row.Add(Acc);

                                if (iIter == 0) {
                                    double xCoord = -iLv + MglMax + 1.0 + (p) * 0.1;
                                    xCoords.Add(xCoord);
                                }
                            }
                        }

                        WaterfallData.Add(Row.ToArray());
                    }

                    Ret.LogY = true;
                    Ret.ShowLegend = false;

                    for (int iIter = 1; iIter <= Math.Min(NoOfIter, MaxIter); iIter++) {
                        var PlotRow = WaterfallData[iIter];

                        var g = new Plot2Ddata.XYvalues("iter" + iIter, xCoords, PlotRow);
                        if (i == 1) {
                            // error
                            g.Format.PointType = PointTypes.Circle;
                            g.Format.LineColor = LineColors.Black;
                        } else {
                            // residual
                            g.Format.PointType = PointTypes.OpenBox;
                            g.Format.LineColor = LineColors.Red;
                        }
                        Ret.AddDataGroup(g);

                    }
                    //RetDict.Add("Waterfall_of_Var" + iVar, Ret);
                }
            }
            return RetDict;
        }


        /// <summary>
        /// Writes the table obtained through <see cref="WriteTrendToTable(bool, bool, bool, bool, out string[], out MultidimensionalArray)"/> into a CSV file.
        /// </summary>
        public void WriteTrendToCSV(bool ErrorOrResidual, bool SepVars, bool SepPoly, bool SepLev, string name) {
            string[] Titels;
            MultidimensionalArray ConvTrendData;
            WriteTrendToTable(ErrorOrResidual, SepVars, SepPoly, SepLev, out Titels, out ConvTrendData);


            using (StreamWriter stw = new StreamWriter(name)) {
                stw.WriteLine(Titels.CatStrings(" "));
                stw.Flush();
                ConvTrendData.SaveToStream(stw);
                stw.Flush();
            }
        }
        
        
        /// <summary>
        /// provides data collected during a solver run in tabular form
        /// </summary>
        /// <param name="ErrorOrResidual">
        /// - true: output is error against exact solution (if provided) during each iteration
        /// - false: output is residual against exact solution (if provided) during each iteration
        /// </param>
        /// <param name="SepVars">
        /// - true: separate column for each variable
        /// - false:  l2-norm over all variables
        /// </param>
        /// <param name="SepPoly">
        /// - true:   separate column for each polynomial degree
        /// - false:  l2-norm over all polynomial degrees
        /// </param>
        /// <param name="SepLev">
        /// - true:   separate column for each multi-grid level
        /// - false:  l2-norm over all multi-grid levels
        /// </param>
        /// <param name="Titels">
        /// Column names/titles
        /// </param>
        /// <param name="ConvTrendData">
        /// data table:
        /// - columns: correspond to <paramref name="Titels"/>
        /// - rows: solver iterations
        /// </param>
        public void WriteTrendToTable(bool ErrorOrResidual, bool SepVars, bool SepPoly, bool SepLev, out string[] Titels, out MultidimensionalArray ConvTrendData) {

            var data = ErrorOrResidual ? ErrNormTrend : ResNormTrend;

            List<string> titleS = new List<string>();
            int iCol = 0;
            Dictionary<(int MGlevel, int iVar, int deg), int> ColumnIndex = new Dictionary<(int, int, int), int>();
            foreach (var kv in data) {
                int iLevel = kv.Key.Item1;
                int pDG = kv.Key.Item3;
                int iVar = kv.Key.Item2;

                string title;
                {
                    if(SepVars)
                        title = $"Var{iVar}";
                    else
                        title = "";

                    if(SepLev) {
                        if(title.Length > 0)
                            title = title + ",";
                        title = title + $"Mglv{iLevel}";
                    }

                    if(SepPoly) {
                        if(title.Length > 0)
                            title = title + ",";
                        title = title + $"p={pDG}";
                    }

                    if(title.Length > 0)
                        title = "_" + title;

                    if(ErrorOrResidual)
                        title = "Err" + title;
                    else
                        title = "Res" + title;
                }
                /*
                if (SepPoly == false && SepLev == false) {
                    title = string.Format("(var#{0})", iVar);
                } else if (SepPoly == false && SepLev == true) {
                    title = string.Format("(var#{0},mg.lev.{1})", iVar, iLevel);
                } else if (SepPoly == true && SepLev == false) {
                    title = string.Format("(var#{0},p={1})", iVar, pDG);
                } else if (SepPoly == true && SepLev == true) {
                    title = string.Format("(var#{0},mg.lev.{1},p={2})", iVar, iLevel, pDG);
                } else {
                    throw new ApplicationException();
                }
                */

                int iColKv = titleS.IndexOf(title);
                if (iColKv < 0) {
                    titleS.Add(title);
                    iColKv = iCol;
                    iCol++;
                }
                ColumnIndex.Add(kv.Key, iColKv);
            }
            Titels = titleS.ToArray();
            ConvTrendData = MultidimensionalArray.Create(data.First().Value.Count, titleS.Count);

            foreach (var kv in data) { // over all data columns
                
                double[] ColumnSquared = kv.Value.Select(val => val.Pow2()).ToArray();

                if (ConvTrendData.GetLength(0) != ColumnSquared.Length)
                    throw new ApplicationException();

                int iColkv = ColumnIndex[kv.Key]; // output column index in which we sum up the column
                ConvTrendData.ExtractSubArrayShallow(-1, iColkv).AccVector(1.0, ColumnSquared);

            }
            ConvTrendData.ApplyAll(x => Math.Sqrt(x));
        }

        /// <summary>
        /// Decomposition of some solution vector <paramref name="vec"/> into the different multigrid levels.
        /// </summary>
        /// <param name="vec">
        /// approximate solution
        /// </param>
        /// <param name="plotName">
        /// name of the tecplot file
        /// </param>
        public void PlotDecomposition<V>(V vec, string plotName)
            where V : IList<double> //
        {
            int L0 = this.DecompositionOperator.Mapping.LocalLength;
            double[] vec0 = new double[L0];
            this.DecompositionOperator.TransformSolInto(vec, vec0);
            var Decomp = this.OrthonormalMultigridDecomposition(vec0, false);

            List<DGField> DecompFields = new List<DGField>();

            MultigridOperator op = this.DecompositionOperator;
            for (int iLevel = 0; iLevel < Decomp.Count; iLevel++) {
                double[] vec_i = Decomp[iLevel];

                var DecompVec = this.InitProblemDGFields("Level" + iLevel);

                MultigridOperator opi = op;
                for (int k = iLevel; k > 0; k--) {
                    int Lk1 = Decomp[k - 1].Length;
                    double[] vec_i1 = new double[Lk1];
                    opi.Prolongate(1.0, vec_i1, 0.0, vec_i);
                    vec_i = vec_i1;
                    opi = opi.FinerLevel;
                }

                this.DecompositionOperator.TransformSolFrom(DecompVec, vec_i);

                DecompFields.AddRange(DecompVec.Mapping.Fields);

                op = op.CoarserLevel;
            }

            Tecplot.Tecplot.PlotFields(DecompFields, plotName, 0.0, 3);
        }


        /// <summary>
        /// Decomposition of a certain vector into its frequencies on all mesh levels.
        /// </summary>
        /// <param name="Vec"></param>
        /// <param name="decompose">
        /// - true: orthogonality (of the respective DG/XDG) representation across all mesh levels:
        ///         the high level modes contain only contributions which cannot be represented on the lower levels.
        /// - false: the higher level meshes include also lower frequencies;
        /// </param>
        /// <returns>
        /// - one vector per multigrid level
        /// </returns>
        public IList<double[]> OrthonormalMultigridDecomposition(double[] Vec, bool decompose = true) {
            // vector length on level 0
            int L0 = DecompositionOperator.Mapping.LocalLength;
            if (Vec.Length != L0)
                throw new ArgumentException("Mismatch in vector length.", "Vec");


            List<double[]> OrthoVecs = new List<double[]>();
            OrthoVecs.Add(Vec.CloneAs());

            
            for (var mgop = this.DecompositionOperator.CoarserLevel; mgop != null; mgop = mgop.CoarserLevel) {
                int L = mgop.Mapping.LocalLength;
                int iLevel = mgop.LevelIndex;
                OrthoVecs.Add(new double[L]);

                mgop.Restrict(OrthoVecs[iLevel - 1], OrthoVecs[iLevel]);


                if(decompose) {
                    double L2higher = OrthoVecs[iLevel - 1].MPI_L2NormPow2();
                    mgop.Prolongate(-1.0, OrthoVecs[iLevel - 1], 1.0, OrthoVecs[iLevel]);

                    double L2higher2 = OrthoVecs[iLevel - 1].MPI_L2NormPow2();
                    double L2lower = OrthoVecs[iLevel].MPI_L2NormPow2();
                    double L2comb = L2higher2 + L2lower;

                }
            }


            return OrthoVecs;
        }

        /// <summary>
        /// Basis filename for the Tecplot output.
        /// </summary>
        public string TecplotOut = null;

        /// <summary>
        /// L2 norm of residual in last/most recent iteration
        /// </summary>
        public double LastResidualNorm;

        /// <summary>
        /// L2 norm of solution (= error against 0 RHS) in last/most recent iteration
        /// </summary>
        public double LastSolNorm;


        /// <summary>
        /// L2 norm of residual in first iteration
        /// </summary>
        public double Iter0ResidualNorm;

        /// <summary>
        /// L2 norm of solution (= error against 0 RHS) in first recent iteration
        /// </summary>
        public double Iter0SolNorm;


        private bool firstLog = false;

        /// <summary>
        /// Highest iteration number
        /// </summary>
        public int NumberOfIterations = 0;

        /// <summary>
        /// Callback routine, see <see cref="ISolverWithCallback.IterationCallback"/> or <see cref="NonlinearSolver.IterationCallback"/>.
        /// </summary>
        public void IterationCallback(int iter, double[] xI, double[] rI, MultigridOperator mgOp) {
            if (xI.Length != SolverOperator.Mapping.LocalLength)
                throw new ArgumentException();
            if (rI.Length != SolverOperator.Mapping.LocalLength)
                throw new ArgumentException();

            NumberOfIterations = Math.Max(NumberOfIterations, iter);

            int Lorg = SolverOperator.BaseGridProblemMapping.LocalLength;

            // transform residual and solution back onto the original grid
            // ==========================================================

            double[] Res_Org = new double[Lorg];
            double[] Sol_Org = new double[Lorg];

            SolverOperator.TransformRhsFrom(Res_Org, rI);
            SolverOperator.TransformSolFrom(Sol_Org, xI);

            double[] Err_Org = Sol_Org.CloneAs();
            Err_Org.AccV(-1.0, this.ExactSolution);

            if (TecplotOut != null) {
                var ErrVec = InitProblemDGFields("Err");
                var ResVec = InitProblemDGFields("Res");
                var SolVec = InitProblemDGFields("Sol");

                ErrVec.SetV(Err_Org);
                ResVec.SetV(Res_Org);
                SolVec.SetV(Sol_Org);
                List<DGField> ErrResSol = new List<DGField>();
                ErrResSol.AddRange(ErrVec.Mapping.Fields);
                ErrResSol.AddRange(ResVec.Mapping.Fields);
                ErrResSol.AddRange(SolVec.Mapping.Fields);

                Tecplot.Tecplot.PlotFields(ErrResSol, TecplotOut + "." + iter, iter, 4);

                PlotDecomposition(xI, TecplotOut + "-sol-decomp." + iter);
                PlotDecomposition(rI, TecplotOut + "-res-decomp." + iter);
            }

            // Console out
            // ===========
            double l2_RES = rI.L2NormPow2().MPISum().Sqrt();
            double l2_ERR = Err_Org.L2NormPow2().MPISum().Sqrt();
            Console.WriteLine("Iter: {0}\tRes: {1:0.##E-00}\tErr: {2:0.##E-00}", iter, l2_RES, l2_ERR);

            LastResidualNorm = l2_RES;
            LastSolNorm = l2_ERR;
            if(firstLog == false) {
                Iter0ResidualNorm = l2_RES;
                Iter0SolNorm = l2_ERR;
                firstLog = true;
            }


            // decompose error and residual into orthonormal vectors
            // =====================================================


            int L0 = DecompositionOperator.Mapping.LocalLength;
            double[] Err_0 = new double[L0], Res_0 = new double[L0];
            DecompositionOperator.TransformSolInto(Err_Org, Err_0);
            DecompositionOperator.TransformRhsInto(Res_Org, Res_0, false);

            IList<double[]> Err_OrthoLevels = OrthonormalMultigridDecomposition(Err_0, decompose:true);
            IList<double[]> Res_OrthoLevels = OrthonormalMultigridDecomposition(Res_0, decompose:true);


            // compute L2 norms on each level
            // ==============================
            for (var mgop = this.DecompositionOperator; mgop != null; mgop = mgop.CoarserLevel) {
                int[] _Degrees = mgop.Mapping.DgDegree;

                double[] Resi = Res_OrthoLevels[mgop.LevelIndex];
                double[] Errr = Err_OrthoLevels[mgop.LevelIndex];
                int JAGG = mgop.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;


                for (int iVar = 0; iVar < _Degrees.Length; iVar++) { // loop over variables...
                    for (int p = 0; p <= _Degrees[iVar]; p++) { // loop over degrees...
                        List<double> ResNorm = this.ResNormTrend[(mgop.LevelIndex, iVar, p)];
                        List<double> ErrNorm = this.ErrNormTrend[(mgop.LevelIndex, iVar, p)];

                        double ResNormAcc = 0.0;
                        double ErrNormAcc = 0.0;

                        for (int jagg = 0; jagg < JAGG; jagg++) { // sum of var 'iVar' of all modes with degree 'p' over all cells
                            int[] NN = mgop.Mapping.AggBasis[iVar].ModeIndexForDegree(jagg, p, _Degrees[iVar]);

                            foreach (int n in NN) {
                                int idx = mgop.Mapping.LocalUniqueIndex(iVar, jagg, n);

                                ResNormAcc += Resi[idx].Pow2();
                                ErrNormAcc += Errr[idx].Pow2();
                            }
                        }

                        ResNorm.Add(ResNormAcc.Sqrt());
                        ErrNorm.Add(ErrNormAcc.Sqrt());
                    }
                }
            }
        }

        

        /// <summary>
        /// 
        /// </summary>
        /// <param name="FsDriver"></param>
        /// <param name="SI"></param>
        public void WriteTrendToSession(IFileSystemDriver FsDriver, SessionInfo SI) {
            this.WriteTrendToTable(false, true, true, true, out string[] columns, out MultidimensionalArray table);

            int MPIrank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);

            if ((MPIrank == 0) && (SI.ID != Guid.Empty)) {
                var LogRes = FsDriver.GetNewLog("ResTrend", SI.ID);
                foreach (var col in columns) LogRes.Write(col + "\t");
                int nocol = columns.Length;
                int norow = table.GetLength(0);
                Debug.Assert(nocol == table.GetLength(1));
                LogRes.WriteLine();
                for (int iRow = 0; iRow < norow; iRow++) {
                    for (int iCol = 0; iCol < nocol; iCol++) {
                        LogRes.Write(table[iRow, iCol] + "\t");
                    }
                    LogRes.WriteLine();
                }
                LogRes.Flush();
            }
        }

        private int[] Iterationcounter {
            get {
                /*
                if (m_SF.GetIterationcounter == null)
                    throw new ArgumentNullException("switch verbose mode on for the solver you like to plot! Iterationcounter is null!");
                Debug.Assert(m_SF.GetIterationcounter.Length == 6);
                return m_SF.GetIterationcounter;
                */
                throw new NotImplementedException();
            }
        }

        /*
        public void ResItCallbackAtAll(int iter, double[] xI, double[] rI, MultigridOperator mgOp) {

            var Ptr_mgOp = mgOp;
            int iLevel = mgOp.LevelIndex;
            double[] vec_i = rI;

            CoordinateVector DecompVec = this.InitProblemDGFields("Res");

            for (int k = iLevel; k > 0; k--) {
                double[] vec_i1 = new double[Ptr_mgOp.FinerLevel.Mapping.LocalLength];
                Ptr_mgOp.Prolongate(1.0, vec_i1, 0.0, vec_i);
                vec_i = vec_i1;
                Ptr_mgOp = Ptr_mgOp.FinerLevel;
            }

            Ptr_mgOp.TransformSolFrom(DecompVec, vec_i);
            //DecompVec.AccV(1, vec_i);

            string plotName = TecplotOut + "Res-decomp" + "."+Iterationcounter[3] + "."+ ItWithinMGCycle();


            Tecplot.Tecplot.PlotFields(DecompVec.Mapping.Fields, plotName, 0.0, 3);
            //DecomposedDGFields.AddRange(DecompVec.Mapping.Fields);
        }
        */

        private int CurrentMLevel_down=0;
        private int CurrentMLevel_up = 0;

        private bool IsDownstep(int currentIt)
        {
            bool stepdown = false;
            if (CurrentMLevel_down - currentIt == -1)
            {
                stepdown = true;
            } 
            if (ItWithinMGCycle() == 1)
            {
                stepdown = true;
            }

            CurrentMLevel_down = currentIt;
            return stepdown;
        }

        private bool IsUpstep (int currentIt)
        {
            bool stepup = false;
            if (CurrentMLevel_up - currentIt == +1)
            {
                stepup = true;
            }

            CurrentMLevel_up = currentIt;
            return stepup;
        }

        private int CurrentMCycle = 0;
        private int MG_internal_counter = 0;

        private int ItWithinMGCycle() {
                if(Iterationcounter[3]> CurrentMCycle)
                {
                    CurrentMCycle = Iterationcounter[3];
                    MG_internal_counter = Iterationcounter[5]-1;
                }
            return Iterationcounter[5] - MG_internal_counter;
        }

        /*
        public void ResItCallbackAtDownstep(int iter, double[] xI, double[] rI, MultigridOperator mgOp)
        {

            if (IsDownstep(mgOp.LevelIndex)||IsUpstep(mgOp.LevelIndex))
            {
                var Ptr_mgOp = mgOp;
                int iLevel = mgOp.LevelIndex;
                double[] vec_i = rI;

                CoordinateVector DecompVec = this.InitProblemDGFields("Res");

                for (int k = iLevel; k > 0; k--)
                {
                    double[] vec_i1 = new double[Ptr_mgOp.FinerLevel.Mapping.LocalLength];
                    Ptr_mgOp.Prolongate(1.0, vec_i1, 0.0, vec_i);
                    vec_i = vec_i1;
                    Ptr_mgOp = Ptr_mgOp.FinerLevel;
                }

                Ptr_mgOp.TransformSolFrom(DecompVec, vec_i);
                //DecompVec.AccV(1, vec_i);

                string plotName = TecplotOut + "Res-decomp" + "." + Iterationcounter[3] + "." + ItWithinMGCycle();
                Tecplot.Tecplot.PlotFields(DecompVec.Mapping.Fields, plotName, 0.0, 3);
            }
        }
        */

        private CoordinateVector InitProblemDGFields(string NamePrefix) {
            Basis[] BS = this.SolverOperator.BaseGridProblemMapping.BasisS.ToArray();
            DGField[] Fields = new DGField[BS.Length];

            for(int iFld = 0; iFld < BS.Length; iFld++) {
                var name = string.Format("{0}_var_{1}", NamePrefix, iFld);
                if(BS[iFld] is BoSSS.Foundation.XDG.XDGBasis) {
                    Fields[iFld] = new BoSSS.Foundation.XDG.XDGField((BoSSS.Foundation.XDG.XDGBasis)BS[iFld], name);
                } else {
                    Fields[iFld] = new SinglePhaseField(BS[iFld], name);
                }
            }

            return new CoordinateVector(Fields);
        }

        // extract the Fields from the solution, Resample them equally spaced and ready to use in an fft
        public void Resample(int iterIndex, double[] currentSol, string component) {
            var Mgop = this.SolverOperator;

            if (Mgop.GridData.SpatialDimension == 2 && Mgop.LevelIndex == 0) {
                MultidimensionalArray SamplePoints;

                GridData GD = (GridData)Mgop.Mapping.AggGrid.AncestorGrid;

                BoundingBox BB = GD.GlobalBoundingBox;

                double xDist = BB.Max[0] - BB.Min[0];
                double yDist = BB.Max[1] - BB.Min[1];
                double aspectRatio = xDist / yDist;

                MGViz viz = new MGViz(Mgop);
                DGField[] Fields = viz.ProlongateToDg(currentSol, "Error");

                for (int p = 0; p < Fields.Length; p++) {
                    var field = Fields[p];

                    int DOF = field.DOFLocal;
                    double N = Math.Sqrt(DOF);
                    int Nx = (int)Math.Round(Math.Sqrt(aspectRatio) * N);
                    int Ny = (int)Math.Round(1 / Math.Sqrt(aspectRatio) * N);

                    // make sure we have always uneven number of frequencys
                    if(Nx % 2 == 0) Nx++;
                    if (Ny % 2 == 0) Ny++;

                    SamplePoints = MultidimensionalArray.Create(Ny, Nx);

                    for (int i = 0; i < Nx; i++) {
                        MultidimensionalArray points = MultidimensionalArray.Create(Ny, 2);

                        for (int k = 0; k < Ny; k++) {
                            points[k, 0] = BB.Min[0] + (i + 1) * xDist / (Nx + 1);
                            points[k, 1] = BB.Min[1] + (k + 1) * yDist / (Ny + 1);
                        }

                        List<DGField> fields = new List<DGField>();
                        fields.Add(field);

                        FieldEvaluation FE = new FieldEvaluation(GD);

                        MultidimensionalArray Result = MultidimensionalArray.Create(Ny, 1);

                        FE.Evaluate(1.0, fields, points, 1.0, Result);

                        //points.SaveToTextFile("points_of_"+i+"th_col_"+ p + "th_field");

                        SamplePoints.ExtractSubArrayShallow(-1, i).Acc(1.0, Result.ExtractSubArrayShallow(-1, 0));
                    }

                    SamplePoints.SaveToTextFile("ResampleFFT_lvl" + Mgop.LevelIndex + "_" + iterIndex + "_" + component + "_" + field.Identification + ".txt");
                }

            }
            if (Mgop.GridData.SpatialDimension == 3 && Mgop.LevelIndex == 0) {
                MultidimensionalArray SamplePoints;

                GridData GD = (GridData)Mgop.Mapping.AggGrid.AncestorGrid;

                BoundingBox BB = GD.GlobalBoundingBox;

                double xDist = BB.Max[0] - BB.Min[0];
                double yDist = BB.Max[1] - BB.Min[1];
                double zDist = BB.Max[2] - BB.Min[2];

                MGViz viz = new MGViz(Mgop);
                DGField[] Fields = viz.ProlongateToDg(currentSol, "Error");
                double xy_ratio = yDist / xDist;
                double xz_ratio = zDist / xDist;

                for (int p = 0; p < Fields.Length; p++) {
                    var field = Fields[p];

                    int DOF = field.DOFLocal;
                    int Nx = (int)Math.Round(Math.Sqrt(1 / xy_ratio * 1 / xz_ratio * DOF));
                    int Ny = (int)Math.Round(Nx * xy_ratio * DOF);
                    int Nz = (int)Math.Round(Nx * xz_ratio * DOF);

                    SamplePoints = MultidimensionalArray.Create(Nx, Ny, Nz);

                    for (int x = 0; x < Nx; x++) {
                        for (int y = 0; y < Ny; y++) {
                            MultidimensionalArray points = MultidimensionalArray.Create(Nz, 3);
                            for (int z = 0; z < Nz; z++) {
                                points[z, 0] = BB.Min[0] + (x + 1) * xDist / (Nx + 1);
                                points[z, 1] = BB.Min[1] + (y + 1) * yDist / (Ny + 1);
                                points[z, 2] = BB.Min[2] + (z + 1) * yDist / (Nz + 1);
                            }

                            List<DGField> fields = new List<DGField>();
                            fields.Add(field);

                            FieldEvaluation FE = new FieldEvaluation(GD);

                            MultidimensionalArray Result = MultidimensionalArray.Create(Nz, 1);

                            FE.Evaluate(1.0, fields, points, 1.0, Result);

                            SamplePoints.ExtractSubArrayShallow(x, y, -1).Acc(1.0, Result.ExtractSubArrayShallow(-1, 0));
                        }
                    }
                    // Entweder neues Format nötig oder Frickel-Lsg alias wie hoch ist der Leidensdruck?
                    //SamplePoints.SaveToTextFile("ResampleFFT_lvl" + Mgop.LevelIndex + "_" + iterIndex + "_" + component + "_" + field.Identification + ".txt");
                }
            }
        }
    }
}
