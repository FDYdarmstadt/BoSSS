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
using System.Diagnostics;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Solution.LevelSetTools.EllipticReInit {

    public static class DiffusionRates {
        public static double DoubleWell(double s, bool eval) {
            if (s >= 1) {
                return 1 / s;
            }
            else {
                return 3 * s - 2 * s * s;
            }
        }

        public static double SingleWell(double x, bool eval) {
            return 1 / x;
        }

        /// <summary>
        /// EllipticReinit Fluxes According to the potential function
        /// p_4= (1/2)*x^1*(x-1)^2
        /// </summary>
        public static double DoubleWellAlternative(double s, bool eval) {
            if (s >= 1) {
                return 1 / s;
            }
            else {
                return 3 - (3 / 2) * s - 1 / (2 * s);
                //(3 / 4) * s - 1 / (4 * s);
                //1-(5 / 4) * Math.Sqrt(s) + 3 / (2 *Math.Sqrt(s)) - 1 / (4 * s.Pow(3 / 2));

            }
        }

        public static double SingleWellNear(double x, bool eval) {
            if (eval) { return 1 / x; }
            else { return 0; }
        }

        public static double SingleWellOnCutDoubleWellElse(double x, bool eval) {
            if (eval) return SingleWell(x, eval);
            else return DoubleWell(x, eval);
        }
    }


    public abstract class RHSForm : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IObserver<LevelSetTracker.LevelSetRegions> {

        public RHSForm(double PenaltyBase, LevelSetTracker LSTrck) {
            this.LSTrck = LSTrck;
            NearFieldBitMask = LSTrck.Regions.GetNearFieldMask(1).GetBitMaskWithExternal();
            CutCells = LSTrck.Regions.GetNearFieldMask(0).GetBitMask();
        }

        LevelSetTracker LSTrck;
        internal System.Collections.BitArray NearFieldBitMask;
        internal System.Collections.BitArray CutCells;

        public TermActivationFlags VolTerms
        {
            get
            {
                return (TermActivationFlags.GradUxGradV);
            }
        }

        /// <summary>
        /// Diffusion Rate
        /// </summary>
        /// <param name="d">Abs(Grad(Phi))</param>
        /// <param name="eval">Evaluation Flag: True, if the Cell is in the Near Region</param>
        /// <returns></returns>
        public Func<double, bool, double> DiffusionRate;


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = cpv.D;
            Debug.Assert(GradU.GetLength(0) == 1);
            Debug.Assert(GradU.GetLength(1) == D);

            double AbsGradU = 0;
            double Acc = 0;
            for (int d = 0; d < D; d++) {
                AbsGradU += GradU[0, d] * GradU[0, d];
                Acc += GradU[0, d] * GradV[d];
            }
            AbsGradU = Math.Sqrt(AbsGradU);

            Acc *= DiffusionRate(AbsGradU, CutCells[cpv.jCell]);
#if DEBUG
            if (Acc.IsNaN()) throw new ArithmeticException();
#endif
            return - Acc;
        }


        public virtual IList<string> ParameterOrdering { get; }
        

        public TermActivationFlags BoundaryEdgeTerms
        {
            get
            {
                return TermActivationFlags.GradUxV;
            }
        }

        public abstract TermActivationFlags InnerEdgeTerms { get; }

        public IList<string> ArgumentOrdering
        {
            get
            {
                return new string[] { VariableNames.LevelSet };
            }
        }

  


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] uA, double[,] Grad_uA, double vA, double[] Grad_vA) {
            return 0;
        }

        public abstract double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT);

        public void OnNext(LevelSetTracker.LevelSetRegions value) {
            CutCells = LSTrck.Regions.GetCutCellMask().GetBitMask();
            NearFieldBitMask = LSTrck.Regions.GetNearFieldMask(1).GetBitMaskWithExternal();
        }

        public void OnError(Exception error) {
            // Do Nothing
        }

        public void OnCompleted() {
            // Do nothing
        }

    }
}
