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

namespace BoSSS.Solution.Multigrid
{
    public class LocalizedOperatorPrec : ISolverSmootherTemplate, ISolverWithCallback
    {
        public int IterationsInNested => throw new NotImplementedException();

        public int ThisLevelIterations => throw new NotImplementedException();

        public bool Converged => throw new NotImplementedException();

        public Action<int, double[], double[], MultigridOperator> IterationCallback { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        int D;
        string[] DomName;
        string[] CodName;
        string[] Params;

        string[] CodNameSelected = new string[0];
        string[] DomNameSelected = new string[0];

        public double m_dt = 1;
        public double m_muA;

        SpatialOperator LocalOp;

        public void Init(MultigridOperator op)
        {
            int D = op.GridData.SpatialDimension;

            CodName = (new string[] { "mom0", "mom1" });
            Params = ArrayTools.Cat(
                VariableNames.Velocity0Vector(D));
            DomName = ArrayTools.Cat(VariableNames.VelocityVector(D));

            LocalOp = new SpatialOperator(DomName, Params,CodName, (A, B, C) => 1);

            for (int d = 0; d < D; d++)
            {

                LocalOp.EquationComponents["mom" + d].Add(
                    new LocalDiffusiveFlux() { m_component = d, dt = m_dt, muA = m_muA });             
            }

            LocalOp.Commit();

            MsrMatrix LocalMatrix = op.MassMatrix.CloneAs().ToMsrMatrix();
            LocalMatrix.Clear();

            UnsetteledCoordinateMapping test = new UnsetteledCoordinateMapping(op.BaseGridProblemMapping.BasisS.GetSubVector(0,D));

            var U0_U0mean = new SinglePhaseField[D];

            LocalMatrix =  LocalOp.ComputeMatrix(test,U0_U0mean , test);

            LocalMatrix.SaveToTextFileSparse("LocalConvDiffMatrix");

        }

        public void ResetStat()
        {
            throw new NotImplementedException();
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>
        {
            throw new NotImplementedException();
        }

    }


    class LocalDiffusiveFlux : IVolumeForm, IEdgeForm
    {
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
                return new string[] { VariableNames.Velocity_d(m_component) };
            }
        }

        public int m_component;

        public IList<string> ParameterOrdering {
            get {
                return VariableNames.Velocity0Vector(2);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV | TermActivationFlags.UxV | TermActivationFlags.GradUxV;
            }
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            //throw new NotImplementedException();

            return 0.0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
        {
            int D = inp.D;

            double acc = 0;
            for (int d = 0; d < D; d++)
            {
                acc += (-_Grad_uIN[0, d] * inp.Normale[d] * _vIN - _Grad_uOUT[0, d] * inp.Normale[d] * _vOUT) * muA;
            }

            return acc;

        }

        public double muA;
        public double dt;


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            int D = cpv.D;

            // diffusive
            double acc = 0;
            for (int d = 0; d < D; d++)
            {
                acc += GradU[0, d] * GradV[d] * muA;
            }

            // temporal
            acc += (1 / dt) * U[0] * V;

            // convective
            for (int d = 0; d < D; d++)
            {
                acc += cpv.Parameters[d] * GradU[0, d] * V;
            }


            return acc;
        }

    }
}
