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
using System.Threading.Tasks;

using ilPSP;
using ilPSP.LinSolvers;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;


namespace BoSSS.Solution.XNSECommon.newXSpatialOperator {

    public abstract class XOperatorFactoryBase {

        protected XSpatialOperatorMk2 m_XOp;

        protected LevelSetTracker LsTrk;


        protected string[] DomName;
        protected string[] CodName;
        //protected string[] Params;

        protected int D;


        public XOperatorFactoryBase() {

        }


        /// <summary>
        /// base matrix assembly
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="OpMatrix"></param>
        /// <param name="OpAffine"></param>
        /// <param name="RowMapping"></param>
        /// <param name="ColMapping"></param>
        /// <param name="CurrentState"></param>
        /// <param name="ParamsMap"></param>
        /// <param name="AgglomeratedCellLengthScales"></param>
        /// <param name="time"></param>
        public void AssembleMatrix<T>(BlockMsrMatrix OpMatrix, double[] OpAffine,
            UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping,
            IEnumerable<T> CurrentState, IList<DGField> ParamsMap, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, 
            double time) where T : DGField {

            // checks
            if(ColMapping.BasisS.Count != m_XOp.DomainVar.Count)
                throw new ArgumentException();
            if(RowMapping.BasisS.Count != m_XOp.CodomainVar.Count)
                throw new ArgumentException();
            if(OpMatrix == null && CurrentState == null)
                throw new ArgumentException();

            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

            if(OpMatrix != null) {

                XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = m_XOp.GetMatrixBuilder(LsTrk, ColMapping, ParamsMap, RowMapping, SpcToCompute);

                foreach(var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                }

                mtxBuilder.time = time;

                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = m_XOp.GetEvaluatorEx(this.LsTrk,
                    CurrentState.ToArray(), ParamsMap, RowMapping,
                    SpcToCompute);

                foreach(var kv in AgglomeratedCellLengthScales) {
                    eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                }

                eval.time = time;

                eval.Evaluate(1.0, 1.0, OpAffine);

            }
        }


        public void AssembleMatrix<T>(BlockMsrMatrix OpMatrix, double[] OpAffine,
            UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping,
            IEnumerable<T> CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales,
            double time) where T : DGField {

            AssembleMatrix(OpMatrix, OpAffine, RowMapping, ColMapping, CurrentState, AgglomeratedCellLengthScales, time);
        }
    }
}
