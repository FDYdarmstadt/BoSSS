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
using System.Diagnostics;
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {
    public class GenericRestriction : ISolverSmootherTemplate {

        MultigridOperator m_OpThisLevel;

        public void Init(MultigridOperator op) {
            this.m_OpThisLevel = op;

            if(op.CoarserLevel == null) {
                throw new NotSupportedException("Multigrid algorithm cannot be used as a solver on the finest level.");
            }

            this.CoarserLevelSolver.Init(op.CoarserLevel);
        }

        public ISolverSmootherTemplate CoarserLevelSolver;

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            //
            using (var tr = new FuncTrace()) {


                int N = this.m_OpThisLevel.CoarserLevel.Mapping.LocalLength;
                double[] xc = new double[N], bc = new double[N];
                //this.RestrictionOperator.SpMVpara(1.0, X, 0.0, xc);
                //this.RestrictionOperator.SpMVpara(1.0, B, 0.0, bc);

                this.m_OpThisLevel.CoarserLevel.Restrict(B, bc);
                this.m_OpThisLevel.CoarserLevel.Restrict(X, xc);

                CoarserLevelSolver.Solve(xc, bc);

                //this.PrologateOperator.SpMV(1.0, xc, 0.0, X);
                this.m_OpThisLevel.CoarserLevel.Prolongate(1.0, X, 0.0, xc);
            }
        }
        

        public int IterationsInNested {
            get { 
                return CoarserLevelSolver.IterationsInNested; 
            }
        }

        public int ThisLevelIterations {
            get { return 0; }
        }

        public bool Converged {
            get { return CoarserLevelSolver.Converged; }
        }

        public void ResetStat() {
            CoarserLevelSolver.ResetStat();
        }
        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

    }
    
   
}
