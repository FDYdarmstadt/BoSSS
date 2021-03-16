using AdvancedSolverTests.SubBlocking;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.AdvancedSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Solution.AdvancedSolvers.MultigridOperator;

namespace AdvancedSolverTests.SolverChooser
{
    class Utils
    {

        public static MgoSolverPair CreateTestMGOperator(XDGusage UseXdg = XDGusage.all, int DGOrder = 2, MatrixShape MShape = MatrixShape.laplace, int Resolution = 4) {
            //MultigridOperator retMGOp;
            //using (var solver = new SubBlockTestSolver2Var() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape, m_Res = Resolution }) {
            var solver = new SubBlockTestSolver2Var() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape, m_Res = Resolution };
            solver.Init(null);
            solver.RunSolverMode();
            //retMGOp = solver.MGOp;
            //MGSeq = solver.MgSeq;
            //}
            //return retMGOp;
            return new MgoSolverPair(solver);
        }

        public static ChangeOfBasisConfig[][] GetAllMGConfig(MultigridOperator mgo) {
            List<ChangeOfBasisConfig[]> cobc = new List<ChangeOfBasisConfig[]>();
            return MgDescend(mgo, cobc).ToArray();
        }

        public static AggregationGridBasis[][] GetAllAggGridBasis(MultigridOperator mgo) {
            List<AggregationGridBasis[]> agggb = new List<AggregationGridBasis[]>();
            return MgDescend(mgo, agggb).ToArray();
        }

        private static List<AggregationGridBasis[]> MgDescend(MultigridOperator mgo, List<AggregationGridBasis[]> agggb) {
            agggb.Add(mgo.Mapping.AggBasis);
            if (mgo.CoarserLevel != null)
                return MgDescend(mgo.CoarserLevel, agggb);
            else
                return agggb;
        }

        private static List<ChangeOfBasisConfig[]> MgDescend(MultigridOperator mgo, List<ChangeOfBasisConfig[]> cobc) {
            cobc.Add(mgo.Config);
            if (mgo.CoarserLevel != null)
                return MgDescend(mgo.CoarserLevel, cobc);
            else
                return cobc;
        }

        //private static int MgDescend(MultigridOperator mgo, int Depth) {
        //if (mgo.CoarserLevel != null)
        //    return MgDescend(mgo.CoarserLevel, Depth++);
        //else
        //    return Depth;
        //}

    }
}
