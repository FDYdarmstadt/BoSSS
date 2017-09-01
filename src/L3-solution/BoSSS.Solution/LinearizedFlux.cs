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
using System.Runtime.CompilerServices;
using BoSSS.Foundation;
using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// A flux that represents the linearization of an arbitrary nonlinear
    /// flux. The linearization is computed via a Taylor expansion about the
    /// state given as a parameter (see <see cref="ParameterOrdering" />), while
    /// the involved gradient is approximated numerically by means of a finite
    /// difference.
    /// </summary>
    /// <remarks>
    /// This flux is purely experimental. Both, the accuracy and the efficiency
    /// have not really been tested. Use at own risk!
    /// </remarks>
    public class LinearizedFlux : LinearFlux {

        /// <summary>
        /// The flux to be linearized
        /// </summary>
        private INonlinearFlux baseFlux;

        /// <summary>
        /// Basis for the step sizes computed via <see cref="GetStepSize"/>.
        /// In short, it is given by the square root of the machine precision.
        /// For more information, see 
        /// "Numerical Recipes 3rd Edition: The Art of Scientific Computing"
        /// </summary>
        private static readonly double baseStepSize = Math.Sqrt(1.1102230246251566639e-16);

        /// <summary>
        /// Initializes a new instance of the <see cref="LinearizedFlux"/> class.
        /// </summary>
        /// <param name="baseFlux">The base flux.</param>
        public LinearizedFlux(INonlinearFlux baseFlux) {
            this.baseFlux = baseFlux;

            if (baseFlux.ParameterOrdering != null && baseFlux.ParameterOrdering.Count > 0) {
                throw new NotImplementedException(
                    "Linearized flux currently does not support fluxes with parameters");
            }
        }

        /// <summary>
        /// Computes a step size $h that guarantees that <paramref name="x"/>+h
        /// can exactly represented by a binary number.
        /// </summary>
        /// <param name="x">The point of evaluation</param>
        /// <returns>
        /// A step size $h that minimizes round-off errors.
        /// </returns>
        /// <remarks>
        /// For more information, see
        /// "Numerical Recipes 3rd Edition: The Art of Scientific Computing"
        /// </remarks>
        [MethodImpl(MethodImplOptions.NoOptimization)]
        private double GetStepSize(double x) {
            double temp = x + baseStepSize;
            return temp - x;
        }

        #region ILinearFlux Members

        /// <summary>
        /// Affine-linear flux function on a border edge; Boundary conditions
        /// are implemented here. An affine-linear flux must implemented as a
        /// matrix and an offset vector: The linear border flux, in a
        /// mathematical notation is defined as BorderFlux(U) =
        /// M_in*U + b, where U is a column vector containing function values
        /// of other fields. Which fields in which order is specified by
        /// <see cref="IEquationComponent.ArgumentOrdering" /> property. In
        /// this implementation, the resulting flux solely depends on the
        /// implementation of <see cref="INonlinearFlux.BorderEdgeFlux"/> that
        /// has been supplied to the constructor.
        /// </summary>
        /// <param name="inp">
        /// A set of input parameters such as time, coordinate and values of
        /// parameters. Given as reference for performance reasons; DO NOT
        /// WRITE to this structure.</param>
        /// <param name="Uin"></param>
        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            MultidimensionalArray refFlux = MultidimensionalArray.Create(1, 1);
            baseFlux.BorderEdgeFlux(
                0.0,
                inp.iEdge,
                MultidimensionalArray.CreateWrapper(inp.X, 1, 1, inp.X.Length),
                MultidimensionalArray.CreateWrapper(inp.Normale, 1, inp.Normale.Length),
                false,
                new byte[] { inp.EdgeTag },
                0,
                inp.Parameters_IN.Select(p =>
                    MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                0,
                1,
                refFlux);

            double[] FunctionMatrix = new double[inp.Parameters_IN.Length];
            for (int i = 0; i < inp.Parameters_IN.Length; i++) {
                double[] perturbedParameters = inp.Parameters_IN.CloneAs();
                double stepSize = GetStepSize(perturbedParameters[i]);
                perturbedParameters[i] += stepSize;

                MultidimensionalArray perturbedFlux = MultidimensionalArray.Create(1, 1);
                baseFlux.BorderEdgeFlux(
                    0.0,
                    inp.iEdge,
                    MultidimensionalArray.CreateWrapper(inp.X, 1, 1, inp.X.Length),
                    MultidimensionalArray.CreateWrapper(inp.Normale, 1, inp.Normale.Length),
                    false,
                    new byte[] { inp.EdgeTag },
                    0,
                    perturbedParameters.Select(p =>
                        MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                    0,
                    1,
                    perturbedFlux);

                FunctionMatrix[i] = (perturbedFlux[0, 0] - refFlux[0, 0]) / stepSize;
            }

            double result = refFlux[0, 0];
            for (int i = 0; i < inp.Parameters_IN.Length; i++) {
                result += FunctionMatrix[i] * (Uin[i] - inp.Parameters_IN[i]);
            }
            return result;
        }

        /// <summary>
        /// Affine-linear flux function on an inner edge. An affine-linear flux
        /// must implemented as two matrices and an offset vector: The linear
        /// flux, in a mathematical notation is defined as InnerFlux(Uin,Uout) =
        /// BorderFlux(U) =  M_in*Uin + M_out*Uout + b
        /// where Uin and Uout are column vectors containing function values of
        /// other fields on the respective side of the concerned edge. Which
        /// fields in which order is specified by
        /// <see cref="IEquationComponent.ArgumentOrdering" /> property. In
        /// this implementation, the resulting flux solely depends on the
        /// implementation of <see cref="INonlinearFlux.InnerEdgeFlux"/> that
        /// has been supplied to the constructor.
        /// </summary>
        /// <param name="inp">
        /// A set of input parameters such as time, coordinate and values of
        /// parameters. Given as reference for performance reasons; DO NOT
        /// WRITE to this structure.
        /// </param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            MultidimensionalArray refFlux = MultidimensionalArray.Create(1, 1);
            baseFlux.InnerEdgeFlux(
                0.0,
                inp.iEdge,
                MultidimensionalArray.CreateWrapper(inp.X, 1, 1, inp.X.Length),
                MultidimensionalArray.CreateWrapper(inp.Normale, 1, inp.Normale.Length),
                inp.Parameters_IN.Select(p =>
                    MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                inp.Parameters_OUT.Select(p =>
                    MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                0,
                1,
                refFlux);

            double[] FunctionMatrixIn = new double[inp.Parameters_IN.Length];
            for (int i = 0; i < inp.Parameters_IN.Length; i++) {
                double[] perturbedParameters = inp.Parameters_IN.CloneAs();
                double stepSize = GetStepSize(perturbedParameters[i]);
                perturbedParameters[i] += stepSize;

                MultidimensionalArray perturbedFlux = MultidimensionalArray.Create(1, 1);
                baseFlux.InnerEdgeFlux(
                    0.0,
                    inp.iEdge,
                    MultidimensionalArray.CreateWrapper(inp.X, 1, 1, inp.X.Length),
                    MultidimensionalArray.CreateWrapper(inp.Normale, 1, inp.Normale.Length),
                    perturbedParameters.Select(p =>
                        MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                    inp.Parameters_OUT.Select(p =>
                        MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                    0,
                    1,
                    perturbedFlux);

                FunctionMatrixIn[i] = (perturbedFlux[0, 0] - refFlux[0, 0]) / stepSize;
            }

            double[] FunctionMatrixOut = new double[inp.Parameters_OUT.Length];
            for (int i = 0; i < inp.Parameters_OUT.Length; i++) {
                double[] perturbedParameters = inp.Parameters_OUT.CloneAs();
                double stepSize = GetStepSize(perturbedParameters[i]);
                perturbedParameters[i] += stepSize;

                MultidimensionalArray perturbedFlux = MultidimensionalArray.Create(1, 1);
                baseFlux.InnerEdgeFlux(
                    0.0,
                    inp.iEdge,
                    MultidimensionalArray.CreateWrapper(inp.X, 1, 1, inp.X.Length),
                    MultidimensionalArray.CreateWrapper(inp.Normale, 1, inp.Normale.Length),
                    inp.Parameters_IN.Select(p =>
                        MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                    perturbedParameters.Select(p =>
                        MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                    0,
                    1,
                    perturbedFlux);

                FunctionMatrixOut[i] = (perturbedFlux[0, 0] - refFlux[0, 0]) / stepSize;
            }

            double result = refFlux[0, 0];
            for (int i = 0; i < inp.Parameters_IN.Length; i++) {
                result += FunctionMatrixIn[i] * (Uin[i] - inp.Parameters_IN[i])
                    + FunctionMatrixOut[i] * (Uout[i] - inp.Parameters_OUT[i]);
            }
            return result;
        }

        /// <summary>
        /// The linear flux of the underlying PDE, expressed as F(U) =
        /// M * U + n. The flux may depend on the space
        /// variable <see cref="CommonParamsVol.Xglobal"/> and on the field
        /// values of <see cref="CommonParamsVol.Parameters" /> that can
        /// be defined in <see cref="IEquationComponent.ParameterOrdering" />.
        /// In this implementation, the resulting flux solely depends on the
        /// implementation of <see cref="INonlinearFlux.Flux"/> that
        /// has been supplied to the constructor.
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="U"></param>
        /// <param name="output"></param>
        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            MultidimensionalArray refFlux = MultidimensionalArray.Create(1, 1, inp.Xglobal.Length);
            baseFlux.Flux(
                0.0,
                MultidimensionalArray.CreateWrapper(inp.Xglobal, 1, 1, inp.Xglobal.Length),
                inp.Parameters.Select(p =>
                    MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                0,
                1,
                refFlux);

            double[,] FunctionMatrix = new double[inp.Xglobal.Length, inp.Parameters.Length];
            for (int i = 0; i < inp.Parameters.Length; i++) {
                double[] perturbedParameters = inp.Parameters.CloneAs();
                double stepSize = GetStepSize(perturbedParameters[i]);
                perturbedParameters[i] += stepSize;

                MultidimensionalArray perturbedFlux = MultidimensionalArray.Create(1, 1, inp.Xglobal.Length);
                baseFlux.Flux(
                    0.0,
                    MultidimensionalArray.CreateWrapper(inp.Xglobal, 1, 1, inp.Xglobal.Length),
                    perturbedParameters.Select(p =>
                        MultidimensionalArray.CreateWrapper(new double[] { p }, 1, 1)).ToArray(),
                    0,
                    1,
                    perturbedFlux);

                for (int d = 0; d < inp.Xglobal.Length; d++) {
                    FunctionMatrix[d, i] = (perturbedFlux[0, 0, d] - refFlux[0, 0, d]) / stepSize;
                }
            }

            for (int d = 0; d < inp.Xglobal.Length; d++) {
                output[d] = refFlux[0, 0, d];
                for (int i = 0; i < inp.Parameters.Length; i++) {
                    output[d] += FunctionMatrix[d, i] * (U[i] - inp.Parameters[i]);
                }
            }
        }

        #endregion

        #region IEquationComponent Members

        /// <summary>
        /// <see cref="IEquationComponent.ArgumentOrdering"/>
        /// </summary>
        override public IList<string> ArgumentOrdering {
            get {
                return baseFlux.ArgumentOrdering;
            }
        }

        /// <summary>
        /// Linearization points for each argument defined in
        /// <see cref="ArgumentOrdering"/>.
        /// </summary>
        /// <remarks>
        /// Variable names are currently created by append "0" to each original
        /// argument name. Arguments of the original flux are ignored at the
        /// moment.
        /// </remarks>
        override public IList<string> ParameterOrdering {
            get {
                return baseFlux.ArgumentOrdering.Select(x => x + "0").ToList();
            }
        }

        #endregion
    }
}
