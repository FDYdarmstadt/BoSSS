using BoSSS.Solution.Control;
using ilPSP.LinSolvers;
using NUnit.Framework.Constraints;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {
    public class SchurComplement :ISolverSmootherTemplate {

        [Serializable]
        public class Config :ISolverFactory {
            public string Name => throw new NotImplementedException();

            public string Shortname => throw new NotImplementedException();


            [DataMember]
            public ISolverFactory ConvectionDiffusionSolver = new AdvancedSolvers.OrthoMGSchwarzConfig();

            public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                throw new NotImplementedException();
            }

            public bool Equals(ISolverFactory other) {
                throw new NotImplementedException();
            }
        }




        public int IterationsInNested {
            get;
            private set;
        }
            

        public int ThisLevelIterations {
            get;
            private set;
        }

        public bool Converged {
            get;
            private set;
        }

        public object Clone() {
            throw new NotImplementedException();
        }

        public void Dispose() {
            
        }

        public void Init(MultigridOperator op) {
            int D = op.GridData.SpatialDimension;
            int NoOfVars = op.Mapping.NoOfVariables;
            if(NoOfVars != D+1) {
                throw new NotSupportedException($"expecting {D+1} variables, got {NoOfVars}");
            }

            AggregationGridBasis[] GetCDBasesFromLevel(MultigridOperator op) {
                var cdBases = new AggregationGridBasis[op.Mapping.NoOfVariables - 1];
                for(int i = 0; i < cdBases.Length; i++) {
                    cdBases[i] = op.Mapping.AggBasis[i];
                }
                return cdBases;
            }

            AggregationGridBasis[][] GetCDBasesAllLevels(MultigridOperator op) {
                var cdBases = new AggregationGridBasis[op.NoOfLevels][];
                for(int i = 0; i < cdBases.Length; i++) {
                    cdBases[i] = GetCDBasesFromLevel(op.GetLevel(i));
                }
                return cdBases;
            }



            // extract the sub-matrices
            var VelIdx = op.Mapping.GetSubvectorIndices(D.ForLoop(d => d));
            var PIdx = op.Mapping.GetSubvectorIndices(new int[] { D });
            var CDmatrix = new BlockMsrMatrix()
            op.OperatorMatrix.WriteSubMatrixTo(

                VelIdx, default(long[]), VelIdx, default(long[]));

            var cdMultigridOp = new MultigridOperator(
               GetCDBasesAllLevels(op),
               new Foundation.UnsetteledCoordinateMapping(op.BaseGridProblemMapping.BasisS.Take(D)),
               ,
               ,
               ,
               op.AbstractOperator);
        }






        public void ResetStat() {
            IterationsInNested = 0;
            ThisLevelIterations = 0;
            Converged = false;
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            throw new NotImplementedException();
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}
