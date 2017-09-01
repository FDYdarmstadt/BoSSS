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

namespace BoSSS.Solution.Multigrid {

    /// <summary>
    /// A helper class to analyze and visualize the convergence of an iterative solver algorithm.
    /// - writes tecplot files containing the residuals, error (if the solution is known), <see cref="TecplotOut"/>.
    /// - records convergence trends for different grid and \f$ p \f$-levels.
    /// </summary>
    public class ConvergenceObserver {

        public ConvergenceObserver(MultigridOperator muop, BlockMsrMatrix MassMatrix, double[] __ExactSolution) {
            if(__ExactSolution != null) {
                if(__ExactSolution.Length != muop.BaseGridProblemMapping.LocalLength)
                    throw new ArgumentException();
            }
            this.SolverOperator = muop;

            List<AggregationGridBasis[]> aggBasisSeq = new List<AggregationGridBasis[]>();
            for(var mo = muop; mo != null; mo = mo.CoarserLevel) {
                aggBasisSeq.Add(mo.Mapping.AggBasis);
            }
            
            this.ExactSolution = __ExactSolution;
            int[] Degrees = muop.BaseGridProblemMapping.BasisS.Select(b => b.Degree).ToArray();

            BlockMsrMatrix DummyOpMatrix = new BlockMsrMatrix(muop.BaseGridProblemMapping, muop.BaseGridProblemMapping);
            DummyOpMatrix.AccEyeSp(123);


            MultigridOperator.ChangeOfBasisConfig[][] config = new MultigridOperator.ChangeOfBasisConfig[1][];
            config[0] = new MultigridOperator.ChangeOfBasisConfig[muop.BaseGridProblemMapping.BasisS.Count];
            for(int iVar = 0; iVar < config[0].Length; iVar++) {
                config[0][iVar] = new MultigridOperator.ChangeOfBasisConfig() {
                    Degree = Degrees[iVar],
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { iVar }
                };
            }

            //this.DecompositionOperator = muop; this.DecompositionOperator_IsOrthonormal = false;
            this.DecompositionOperator = new MultigridOperator(aggBasisSeq, muop.BaseGridProblemMapping, DummyOpMatrix, MassMatrix, config); this.DecompositionOperator_IsOrthonormal = true;
            
            ResNormTrend = new Dictionary<Tuple<int, int, int>, List<double>>();
            ErrNormTrend = new Dictionary<Tuple<int, int, int>, List<double>>();
            for(var mgop = this.DecompositionOperator; mgop != null; mgop = mgop.CoarserLevel) {
                int[] _Degrees = mgop.Mapping.DgDegree;
                for(int iVar = 0; iVar < _Degrees.Length; iVar++) {
                    for(int p = 0; p <= _Degrees[iVar]; p++) {
                        ResNormTrend.Add(new Tuple<int, int, int>(mgop.LevelIndex, iVar, p), new List<double>());
                        ErrNormTrend.Add(new Tuple<int, int, int>(mgop.LevelIndex, iVar, p), new List<double>());
                    }
                }
            }

            
        }

        bool DecompositionOperator_IsOrthonormal {
            get;
            set;
        }


        public MultigridOperator DecompositionOperator;
        MultigridOperator SolverOperator;

        double[] ExactSolution;

        /// <summary>
        /// L2-norm of the residual per variable, multigrid level and DG polynomial degree;<br/>
        /// 1st index: multigrid level index <br/>
        /// 2nd index: variable index<br/>
        /// 3rd index: polynomial degree <br/>
        /// </summary>
        Dictionary<Tuple<int,int,int>,List<double>> ResNormTrend;

        /// <summary>
        /// L2-norm of the error per variable, multigrid level and DG polynomial degree;<br/>
        /// 1st index: multigrid level index <br/>
        /// 2nd index: variable index<br/>
        /// 3rd index: polynomial degree <br/>
        /// </summary>
        Dictionary<Tuple<int,int,int>,List<double>> ErrNormTrend;


        public void PlotTrend(bool ErrorOrResidual, bool SepPoly, bool SepLev) {
            using(var gp = new BoSSS.Solution.Gnuplot.Gnuplot()) {
                gp.Cmd("set logscale y");
                gp.Cmd("set title " + (ErrorOrResidual ? "\"Error trend\"" : "\"Residual trend\""));

                string[] Titels;
                MultidimensionalArray ConvTrendData;
                WriteTrendToTable(ErrorOrResidual, SepPoly, SepLev, out Titels, out ConvTrendData);
                double[] IterNo = ConvTrendData.GetLength(0).ForLoop(i => ((double)i));

                for(int iCol = 0; iCol < Titels.Length; iCol++) {
                    gp.PlotXY(IterNo, ConvTrendData.GetColumn(iCol), Titels[iCol],
                        new PlotFormat(lineColor: ((LineColors)(iCol + 1)), Style: Styles.Lines));
                }

                gp.Execute();

                Console.WriteLine("Press any key to continue...");
                Console.ReadKey(true);
                Console.WriteLine("killing gnuplot...");
            }
            Console.WriteLine("done.");
        }

        public void WriteTrendToCSV(bool ErrorOrResidual, bool SepPoly, bool SepLev, string name) {
            string[] Titels;
            MultidimensionalArray ConvTrendData;
            WriteTrendToTable(ErrorOrResidual, SepPoly, SepLev, out Titels, out ConvTrendData);
                

            using(StreamWriter stw = new StreamWriter(name)) {
                stw.WriteLine(Titels.CatStrings(" "));
                stw.Flush();
                ConvTrendData.SaveToStream(stw);
                stw.Flush();
            }

        }

        public void WriteTrendToTable(bool ErrorOrResidual, bool SepPoly, bool SepLev, out string[] Titels, out MultidimensionalArray ConvTrendData) {

            Dictionary<Tuple<int,int,int>,List<double>> data = ErrorOrResidual ? ErrNormTrend : ResNormTrend;

            
            List<string> titleS = new List<string>();
            int iCol = 0;
            Dictionary<Tuple<int,int,int>,int> ColumnIndex = new Dictionary<Tuple<int,int,int>,int>();
            foreach(var kv in data) {
                int iLevel = kv.Key.Item1;
                int pDG = kv.Key.Item3;
                int iVar = kv.Key.Item2;

                string title;
                if(SepPoly == false && SepLev == false) {
                    title = string.Format("(var#{0})", iVar);
                } else if(SepPoly == false && SepLev == true) {
                    title = string.Format("(var#{0},mg.lev.{1})", iVar, iLevel);
                } else if(SepPoly == true && SepLev == false) {
                    title = string.Format("(var#{0},p={2})", iVar, pDG);
                } else if(SepPoly == true && SepLev == true) {
                    title = string.Format("(var#{0},mg.lev.{1},p={2})", iVar, iLevel, pDG);
                } else {
                    throw new ApplicationException();
                }

                int iColKv = titleS.IndexOf(title);
                if(iColKv < 0) {
                    titleS.Add(title);
                    iColKv = iCol;
                    iCol++;
                }
                ColumnIndex.Add(kv.Key, iColKv);
            }
            Titels = titleS.ToArray();
            ConvTrendData = MultidimensionalArray.Create(data.First().Value.Count, titleS.Count);

            foreach(var kv in data) {
                double[] Column = kv.Value.ToArray();
                for(int l = 0; l < Column.Length; l++) {
                    Column[l] = Column[l].Pow2();
                }

                if(ConvTrendData.GetLength(0) != Column.Length)
                    throw new ApplicationException();
                
                int iColkv = ColumnIndex[kv.Key];
                ConvTrendData.ExtractSubArrayShallow(-1, iColkv).AccVector(1.0, Column);
                    
            }

            ConvTrendData.ApplyAll(x => Math.Sqrt(x));
        }


        /// <summary>
        /// Decomposition of some solution vector <paramref name="vec"/> into the different multigrid levels.
        /// </summary>
        public void PlotDecomposition<V>(V vec, string plotName) 
            where V : IList<double> //
        {
            int L0 = this.DecompositionOperator.Mapping.LocalLength;
            double[] Vec0 = new double[L0];
            this.DecompositionOperator.TransformSolInto(vec, Vec0);

            var Decomp = this.OrthonormalMultigridDecomposition(Vec0, false);

            List<DGField> DecompFields = new List<DGField>();

            MultigridOperator op = this.DecompositionOperator;
            for(int iLevel = 0; iLevel < Decomp.Count; iLevel++) {
                double[] vec_i = Decomp[iLevel];
                var DecompVec = this.InitProblemDGFields("Level" + iLevel);

                MultigridOperator opi = op;
                for(int k = iLevel; k > 0; k--) {
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


            Tecplot.Tecplot.PlotFields(DecompFields, plotName, 0.0, 4);
        }


        public IList<double[]> OrthonormalMultigridDecomposition(double[] Vec, bool decompose = true) {
            // vector length on level 0
            int L0 = DecompositionOperator.Mapping.LocalLength;
            if(Vec.Length != L0)
                throw new ArgumentException("Mismatch in vector length.", "Vec");


            List<double[]> OrthoVecs = new List<double[]>();
            OrthoVecs.Add(Vec.CloneAs());

            double l2pow2_Vec = OrthoVecs[0].L2NormPow2().MPISum();

            MultigridOperator coarsest = null;
            for(var mgop = this.DecompositionOperator.CoarserLevel; mgop != null; mgop = mgop.CoarserLevel) {
                int L = mgop.Mapping.LocalLength;
                int iLevel = mgop.LevelIndex;
                OrthoVecs.Add(new double[L]);

                mgop.Restrict(OrthoVecs[iLevel - 1], OrthoVecs[iLevel]);

                coarsest = mgop;
            }

            double Check_l2pow2_OrthoVecsTotal = 0.0;

            for(var mgop = this.DecompositionOperator.CoarserLevel; mgop != null; mgop = mgop.CoarserLevel) {
                int iLevel = mgop.LevelIndex;
                
                
                if(decompose)
                    mgop.Prolongate(-1.0, OrthoVecs[iLevel - 1], 1.0, OrthoVecs[iLevel]);


                coarsest = mgop;
            }
            foreach(double[] OrthoVec in OrthoVecs)
                Check_l2pow2_OrthoVecsTotal += OrthoVec.L2NormPow2().MPISum();

            // if the vectors are really orthonormal, their squared L2-norms must be additive!
            if(decompose && this.DecompositionOperator_IsOrthonormal) 
                Debug.Assert((Math.Abs(Check_l2pow2_OrthoVecsTotal - l2pow2_Vec) / l2pow2_Vec < 1.0e-8), "something wrong with orthonormal decomposition");
            
            return OrthoVecs;
        }

        /// <summary>
        /// Basis filename for the Tecplot output.
        /// </summary>
        public string TecplotOut = null;
        
        /// <summary>
        /// Callback routine, see <see cref="ISolverWithCallback.IterationCallback"/> or <see cref="NonlinearSolver.IterationCallback"/>.
        /// </summary>
        public void IterationCallback(int iter, double[] xI, double[] rI, MultigridOperator mgOp) {
            if(xI.Length != SolverOperator.Mapping.LocalLength)
                throw new ArgumentException();
            if(rI.Length != SolverOperator.Mapping.LocalLength)
                throw new ArgumentException();

            int Lorg = SolverOperator.BaseGridProblemMapping.LocalLength;

            // transform residual and solution back onto the orignal grid
            // ==========================================================

            double[] Res_Org = new double[Lorg];
            double[] Sol_Org = new double[Lorg];

            SolverOperator.TransformRhsFrom(Res_Org, rI);
            SolverOperator.TransformSolFrom(Sol_Org, xI);

            double[] Err_Org = Sol_Org.CloneAs();
            Err_Org.AccV(-1.0, this.ExactSolution);

            if(TecplotOut != null) {


                
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

            }


            // Console out
            // ===========
            double l2_RES = rI.L2NormPow2().MPISum().Sqrt();
            double l2_ERR =  Err_Org.L2NormPow2().MPISum().Sqrt();
            Console.WriteLine("Iter: {0}\tRes: {1:0.##E-00}\tErr: {2:0.##E-00}", iter, l2_RES, l2_ERR);


            // decompose error and residual into orthonormal vectors
            // =====================================================


            int L0 = DecompositionOperator.Mapping.LocalLength;
            double[] Err_0 = new double[L0], Res_0 = new double[L0];
            DecompositionOperator.TransformSolInto(Err_Org, Err_0);  
            DecompositionOperator.TransformRhsInto(Res_Org, Res_0);

            IList<double[]> Err_OrthoLevels = OrthonormalMultigridDecomposition(Err_0);
            IList<double[]> Res_OrthoLevels = OrthonormalMultigridDecomposition(Res_0);
            

            // compute L2 norms on each level
            // ==============================
            for(var mgop = this.DecompositionOperator; mgop != null; mgop = mgop.CoarserLevel) {
                int[] _Degrees = mgop.Mapping.DgDegree;
                
                double[] Resi = Res_OrthoLevels[mgop.LevelIndex];
                double[] Errr = Err_OrthoLevels[mgop.LevelIndex];
                int JAGG = mgop.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;


                for(int iVar = 0; iVar < _Degrees.Length; iVar++) {
                    for(int p = 0; p <= _Degrees[iVar]; p++) {
                        List<double> ResNorm = this.ResNormTrend[new Tuple<int, int, int>(mgop.LevelIndex, iVar, p)];
                        List<double> ErrNorm = this.ErrNormTrend[new Tuple<int, int, int>(mgop.LevelIndex, iVar, p)];

                        double ResNormAcc = 0.0;
                        double ErrNormAcc = 0.0;

                        for(int jagg = 0; jagg < JAGG; jagg++) {
                            int[] NN = mgop.Mapping.AggBasis[iVar].ModeIndexForDegree(jagg, p, _Degrees[iVar]);
                            
                            foreach(int n in NN) {
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
    }
}
