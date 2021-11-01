using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {


    public static class MGopExtensions {

        /// <summary>
        /// Converts a solution vector from some multigrid level into DG fields-
        /// </summary>
        public static DGField[] ProlongateSolToDg<T>(this MultigridOperator op, T V, string baseName)
            where T : IList<double> //
        {
            var viz = new MGViz(op);
            return viz.ProlongateSolToDg(V.ToArray(), baseName);
        }

        /// <summary>
        /// Converts a RHS/residual vector from some multigrid level into DG fields-
        /// </summary>
        public static DGField[] ProlongateRhsToDg<T>(this MultigridOperator op, T V, string baseName) 
            where T : IList<double> //
        {
            var viz = new MGViz(op);
            return viz.ProlongateRhsToDg(V.ToArray(), baseName);
        }


    }



    /// <summary>
    /// Utility class for visualization of intermediate results.
    /// </summary>
    internal class MGViz {

        public MGViz(MultigridOperator op) {
            m_op = op;
        }

        MultigridOperator m_op;

        public int FindLevel(int L) {
            int iLv = 0;
            for (var Op4Level = m_op.FinestLevel; Op4Level != null; Op4Level = Op4Level.CoarserLevel) {
                if (L == Op4Level.Mapping.LocalLength) {
                    Debug.Assert(Op4Level.LevelIndex == iLv);
                    return iLv;
                }
                iLv++;
            }
            return -1;
        }

        public DGField[] ProlongateSolToDg(double[] V, string name) {
            return ProlongateToDg(V, name, false);
        }

        public DGField[] ProlongateRhsToDg(double[] V, string name) {
            return ProlongateToDg(V, name, true);
        }

        public DGField[] ProlongateToDg(double[] V, string name, bool OnRes = false) {
            double[] Curr = ProlongateToTop(V);

            var gdat = m_op.BaseGridProblemMapping.GridDat;
            var basisS = m_op.BaseGridProblemMapping.BasisS;
            DGField[] dgCurrentSol = new DGField[basisS.Count];
            for (int i = 0; i < basisS.Count; i++) {
                var basis = basisS[i];
                if (basis is XDGBasis)
                    dgCurrentSol[i] = new XDGField((XDGBasis)basis, name + i);
                else
                    dgCurrentSol[i] = new SinglePhaseField(basis, name + i);
            }
            CoordinateVector cv = new CoordinateVector(dgCurrentSol);
            if(OnRes )
                m_op.FinestLevel.TransformSolFrom(cv, Curr);
            else 
                m_op.FinestLevel.TransformSolFrom(cv, Curr);
            return dgCurrentSol;
        }

        public double[] ProlongateToTop(double[] V) {
            int iLv = FindLevel(V.Length);

            MultigridOperator op_iLv = m_op.FinestLevel;
            for (int i = 0; i < iLv; i++)
                op_iLv = op_iLv.CoarserLevel;
            Debug.Assert(op_iLv.LevelIndex == iLv);
            Debug.Assert(V.Length == op_iLv.Mapping.LocalLength);

            double[] Curr = V;
            for (var Op4Level = op_iLv; Op4Level.FinerLevel != null; Op4Level = Op4Level.FinerLevel) {
                double[] Next = new double[Op4Level.FinerLevel.Mapping.LocalLength];
                Op4Level.Prolongate(1.0, Next, 0.0, Curr);
                Curr = Next;
            }

            return Curr;
        }

        static int counter = 0;
        public void PlotVectors(IEnumerable<double[]> VV, string[] names, string filename = null) {

            List<DGField> all = new List<DGField>();
            for (int i = 0; i < VV.Count(); i++) {
                all.AddRange(ProlongateToDg(VV.ElementAt(i), names[i]));
            }


            if(filename == null) {
                filename = "MGviz-" + counter;
                counter++;
            }

            Tecplot.Tecplot.PlotFields(all, filename, counter, 2);
            
        }
    }
}
