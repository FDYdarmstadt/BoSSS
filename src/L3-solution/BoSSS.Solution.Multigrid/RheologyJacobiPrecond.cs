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

using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP;
using ilPSP.Connectors.Matlab;
using BoSSS.Solution.Multigrid;

namespace BoSSS.Solution.Multigrid {
    public class RheologyJacobiPrecond : ISolverSmootherTemplate, ISolverWithCallback {
        public int IterationsInNested {
            get {
                return 0;
                //throw new NotImplementedException();
            }
        }

        public int ThisLevelIterations {
            get {
                return 0;
                //throw new NotImplementedException();
            }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get {
                return null;
                //throw new NotImplementedException();
            }
            set {

                //throw new NotImplementedException();
            }
        }

        int D;
        string[] DomName;
        string[] CodName;
        string[] Params;

        public double m_We;

        int[] ConstEqIdx;

        SpatialOperator LocalOp;

        MsrMatrix LocalMatrix;
        BlockMsrMatrix P;

        public void Init(MultigridOperator op) {
            int D = op.GridData.SpatialDimension;

            CodName = new string[] { "momX", "momY", "div", "constitutiveXX", "constitutiveXY", "constitutiveYY" };
            Params = new string[] { "VelocityX_GradientX", "VelocityX_GradientY", "VelocityY_GradientX", "VelocityY_GradientY" }; 
            DomName = new string[] { "VelocityX","VelocityY","Pressure","StressXX", "StressXY", "StressYY" };

            LocalOp = new SpatialOperator(DomName, Params, CodName, (A, B, C) => 4);

            LocalOp.EquationComponents["constitutiveXX"].Add(new LocalJacobiFlux() { m_component = 0, We = m_We });
            LocalOp.EquationComponents["constitutiveXY"].Add(new LocalJacobiFlux() { m_component = 1, We = m_We });
            LocalOp.EquationComponents["constitutiveYY"].Add(new LocalJacobiFlux() { m_component = 2, We = m_We });

            LocalOp.Commit();

            var U0 = ((BoSSS.Foundation.CoordinateMapping)op.Mapping.ProblemMapping).Fields.GetSubVector(0, 2);
            Basis U0Basis = new Basis(op.Mapping.GridData, 0);

            VectorField <SinglePhaseField> VelocityXGradient = new VectorField<SinglePhaseField>(D, U0Basis, "VelocityX_Gradient", SinglePhaseField.Factory);
            VectorField<SinglePhaseField> VelocityYGradient = new VectorField<SinglePhaseField>(D, U0Basis, "VelocityY_Gradient", SinglePhaseField.Factory);

            VelocityXGradient.Clear();
            VelocityXGradient.GradientByFlux(1.0, U0[0]);
            VelocityYGradient.Clear();
            VelocityYGradient.GradientByFlux(1.0, U0[1]);

            var Parameters = ArrayTools.Cat<DGField>(VelocityXGradient, VelocityYGradient);

            LocalMatrix = LocalOp.ComputeMatrix(op.BaseGridProblemMapping, Parameters, op.BaseGridProblemMapping);

            ConstEqIdx = op.Mapping.ProblemMapping.GetSubvectorIndices(true, 3,4,5);

            P = (BlockMsrMatrix)op.MassMatrix.Clone();

            for (int i = ConstEqIdx[0]; i <= ConstEqIdx.Length; i++) {
                for (int j = ConstEqIdx[0]; j <= ConstEqIdx.Length; j++) {
                    if (LocalMatrix[i, j] != 0) {
                        P[i, j] = LocalMatrix[i, j];
                    }
                }
            }

            LocalMatrix.SaveToTextFileSparse("LocalMatrix");
            op.MassMatrix.SaveToTextFileSparse("MassMatrix");
            P.SaveToTextFileSparse("PrecondMatrix");
        }

        public void ResetStat() {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        bool m_Converged = false;
        int m_ThisLevelIterations;

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {

            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(P);
                solver.Solve(X, B);
            }

        }

    }


    public class LocalJacobiFlux : IVolumeForm, IEdgeForm {
        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.GradUxV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                switch (m_component)
                {
                    case 0:
                        return new string[] { VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressXX, VariableNames.StressXY };
                    case 1:
                        return new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.StressXX, VariableNames.StressXY };
                    case 2:
                        return new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.StressXY, VariableNames.StressYY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public int m_component;

        public IList<string> ParameterOrdering {
            get {
                switch (m_component)
                {
                    case 0:
                        return new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY };
                    case 1:
                        return new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                    case 2:
                        return new string[] { VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] Tin, double[] Tout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double Vin, double Vout, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double res = 0.0;
            res += (Tin[0] - Tout[0]) + (Tin[1] - Tout[1]) + (Tin[2] - Tout[2]) + (Tin[3] - Tout[3]);
            // return res * res * (Vin - Vout);
            return 0;
        }

        public double We;


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double acc = 0.0;

            acc += -We * cpv.Parameters[0] * V;


            return acc;
        }

    }
}
