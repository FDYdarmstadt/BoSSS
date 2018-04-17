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
using System.Numerics;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.IntegralTransforms.Algorithms;
using MathNet.Numerics.Interpolation.Algorithms;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Statistic;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.FourierLevelSet {

    /// <summary>
    /// An explicit Interface description based on control points, which form a Fourier Series in the x-y space
    /// </summary>
    public class PlanarFourierLevSet : FourierLevSetBase {

        /// <summary>
        /// Ctr
        /// </summary>
        public PlanarFourierLevSet(FourierLevSetControl Control) : base(Control) {

            setProjectingFourierModes();

            // set the material sample points
            current_interfaceP = new MultidimensionalArray(2);
            current_interfaceP.Allocate(numFp, 2);
            for (int sp = 0; sp < numFp; sp++) {
                current_interfaceP[sp, 0] = FourierP[sp];
                current_interfaceP[sp, 1] = current_samplP[sp];
            }
            setInterfaceLength(current_interfaceP);
            InterfaceResolution = interfaceLength / numFp;

        }


        /// <summary>
        /// set the projecting Fourier modes for the Level set and curvature evaluation
        /// </summary>
        protected override void setProjectingFourierModes() {
            cf_end = (int)Math.Ceiling(DomainSize / (2.0 * h_min));
            if (mode > cf_end)
                Console.WriteLine("WARNING: projected Mode is smaller than grid resolution!");
        }


        /// <summary>
        /// computes the length of the samnple points polygon 
        /// </summary>
        /// <returns></returns>
        protected override void setInterfaceLength(MultidimensionalArray interP) {

            interfaceLength = 0.0;
            int numP = interP.Lengths[0];
            for (int sp = 0; sp < numP - 1; sp++) {
                interfaceLength += Math.Sqrt((interP[sp + 1, 0] - interP[sp, 0]).Pow2()
                    + (interP[sp + 1, 1] - interP[sp, 1]).Pow2());
            }
            interfaceLength += Math.Sqrt((interP[0, 0] + DomainSize - interP[numP - 1, 0]).Pow2()
                                            + (interP[0, 1] - interP[numP - 1, 1]).Pow2());
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="velocity"></param>
        /// <returns></returns>
        public override double[] ComputeChangerate(double dt, ConventionalDGField[] velocity, double[] current_FLSprop) {

            GridData grdat = (GridData)velocity[0].GridDat;
            FieldEvaluation fEval = new FieldEvaluation(grdat);

            MultidimensionalArray VelocityAtSamplePoints = MultidimensionalArray.Create(current_interfaceP.Lengths);

            int outP = fEval.Evaluate(1, velocity, current_interfaceP, 0, VelocityAtSamplePoints);

            if (outP != 0)
                throw new Exception("points outside the grid for fieldevaluation");

            // change rate for the material points is the velocity at the points
            if (FourierEvolve == Fourier_Evolution.MaterialPoints) {
                double[] velAtP = new double[2 * numFp];
                for (int sp = 0; sp < numFp; sp++) {
                    velAtP[sp * 2] = VelocityAtSamplePoints[sp, 0];
                    velAtP[sp * 2 + 1] = VelocityAtSamplePoints[sp, 1];
                }
                return velAtP;
            }

            // compute an infinitesimal change of sample points at the Fourier points/ change of Fourier modes
            MultidimensionalArray interfaceP_evo = current_interfaceP.CloneAs();
            double dt_infin = dt * 1e-3; 
            interfaceP_evo.Acc(dt_infin, VelocityAtSamplePoints);

            RearrangeOntoPeriodicDomain(interfaceP_evo);
            double[] samplP_change = current_samplP.CloneAs();
            InterpolateOntoFourierPoints(interfaceP_evo, samplP_change);

            if (FourierEvolve == Fourier_Evolution.FourierPoints) {
                samplP_change.AccV(-1.0, current_samplP);
                samplP_change.ScaleV(1.0 / dt_infin);
                return samplP_change;
            } else 
            if (FourierEvolve == Fourier_Evolution.FourierModes) {
                Complex[] samplP_complex = new Complex[numFp];
                for (int sp = 0; sp < numFp; sp++) {
                    samplP_complex[sp] = (Complex)samplP_change[sp];
                }
                Complex[] DFT_change = DFT.NaiveForward(samplP_complex, FourierOptions.Matlab);
                double[] DFT_change_double = new double[2 * numFp];
                for (int sp = 0; sp < numFp; sp++) {
                    DFT_change_double[sp * 2] = (DFT_change[sp].Real - DFT_coeff[sp].Real) / dt_infin;
                    DFT_change_double[sp * 2 + 1] = (DFT_change[sp].Imaginary - DFT_coeff[sp].Imaginary) / dt_infin;
                }
                return DFT_change_double;
            } else
                throw new ArgumentException();

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="current_FLSproperty"></param>
        public override void EvolveFourierLS(ref double[] current_FLSproperty) {

            switch (FourierEvolve) {
                case Fourier_Evolution.MaterialPoints:
                    current_interfaceP = MultidimensionalArray.CreateWrapper(current_FLSproperty, new int[] { numFp, 2 });
                    //Reseeding();
                    RearrangeOntoPeriodicDomain(current_interfaceP);
                    // interpolation onto the equidistant Fourier points
                    InterpolateOntoFourierPoints(current_interfaceP, current_samplP);
                    //setProjectingFourierModes();
                    setDFTcoeff();
                    SmoothSamplePoints();

                    current_FLSproperty = new double[numFp * 2];
                    for (int sp = 0; sp < numFp; sp++) {
                        current_interfaceP[sp, 0] = FourierP[sp];
                        current_interfaceP[sp, 1] = current_samplP[sp];
                        current_FLSproperty[sp * 2] = current_interfaceP[sp, 0];
                        current_FLSproperty[(sp * 2) + 1] = current_interfaceP[sp, 1];
                    }

                    break;
                case Fourier_Evolution.FourierPoints:
                    current_samplP = current_FLSproperty;
                    // set material points on Fourier points
                    for (int p = 0; p < numFp; p++) {
                        current_interfaceP[p, 1] = current_samplP[p];
                    }
                    setDFTcoeff();
                    break;
                case Fourier_Evolution.FourierModes:
                    for (int sp = 0; sp < numFp; sp++) {
                        DFT_coeff[sp] = new Complex(current_FLSproperty[sp * 2], current_FLSproperty[sp * 2 + 1]);
                    }
                    // project onto Fourier points
                    for (int fp = 0; fp < numFp; fp++) {
                        if (mode == -1) {
                            //constant values
                            current_samplP[fp] = DFT_coeff[0].Real;
                            if (cf_end == numFp / 2)
                                current_samplP[fp] += Complex.Multiply(DFT_coeff[numFp / 2],
                                                                       Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[fp] * (2 * Math.PI / DomainSize) * (numFp / 2)))
                                                                       ).Real;
                            for (int cf = 1; cf < cf_end; cf++) {
                                Complex cmplx_res = DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[fp] * (2 * Math.PI / DomainSize) * cf));
                                current_samplP[fp] += 2.0 * cmplx_res.Real;  //-f(x) for  k=1... numSp/2-1 , values for k=numSp/2+1...k= numSp-1 are complex conjugated, thus the same real part
                            }
                        } else {
                            if (mode == 0) {
                                current_samplP[fp] += DFT_coeff[0].Real;
                            } else {
                                current_samplP[fp] += Complex.Multiply(DFT_coeff[mode], Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[fp] * (2 * Math.PI / DomainSize) * mode))).Real;
                            }
                            if (mode > 0 && mode < numFp / 2)
                                current_samplP[fp] *= 2.0;
                        }

                        current_samplP[fp] /= (double)numFp;
                        current_interfaceP[fp, 1] = current_samplP[fp];
                    }
                    break;
                default:
                    throw new ArgumentException();
            }
        }


        /// <summary>
        /// LevelSet = y - f(x)
        /// </summary>
        /// <returns></returns>
        public override ScalarFunction PhiEvaluation(int mode) {

            return (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {

                //signal is real valued - i.e. for numSp sample points, numSp/2 coefficients are complex conjugated, thus they don't have to be evaluated individually
                Complex cmplx_res;
                for (int nd = 0; nd < nodes.Lengths[0]; nd++) {

                    if (mode == -1) {

                        //constant values
                        results[nd] -= DFT_coeff[0].Real;
                        if (cf_end == numFp / 2)
                            results[nd] -= Complex.Multiply(DFT_coeff[numFp / 2],
                                                                   Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * (numFp / 2)))
                                                                   ).Real;

                        for (int cf = 1; cf < cf_end; cf++) {
                            cmplx_res = DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * cf));
                            results[nd] -= 2.0 * cmplx_res.Real;  //-f(x) for  k=1... numSp/2-1 , values for k=numSp/2+1...k= numSp-1 are complex conjugated, thus the same real part
                        }

                    } else {

                        if (mode == 0) {
                            results[nd] -= DFT_coeff[0].Real;
                        } else {
                            results[nd] -= Complex.Multiply(DFT_coeff[mode], Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * mode))).Real;
                        }
                        if (mode > 0 && mode < numFp / 2)
                            results[nd] *= 2.0;
                    }

                    results[nd] /= (double)numFp;
                    results[nd] += nodes[nd, 1]; //+y

                }
            };
        }


        //double[] H_exact = new double[] { 0.01,0.00999999927546276,0.00999999710545042,0.00999999349347782,0.00999998844216346,0.00999998195370245,0.00999997403002475,0.00999996467287354,0.00999995388385115,
        //    0.00999994166444868,0.0099999280160665,0.00999991294002918,0.00999989643759664,0.00999987850997288,0.0099998591583129,0.00999983838372824,0.00999981618729166,
        //    0.00999979257004095,0.00999976753298216,0.00999974107709245,0.00999971320332241,0.00999968391259823,0.00999965320582345,0.00999962108388061,0.0099995875476327,
        //    0.00999955259792438,0.00999951623558319,0.00999947846142054,0.00999943927623265,0.00999939868080139,0.00999935667589509,0.0099993132622692,0.00999926844066695};


        //public void ProjectPrescribedLS(SinglePhaseField LevelSet, double time,  LevelSetTracker LsTrk = null) {

        //    CellMask VolMask;
        //    if (LsTrk == null) {
        //        VolMask = CellMask.GetFullMask(LevelSet.Basis.GridDat);
        //    } else {

        //        // set values in positive and negative FAR region to +1 and -1
        //        CellMask Near = LsTrk.GetNearMask4LevSet(0, 1);
        //        CellMask PosFar = LsTrk.GetLevelSetWing(0, +1).VolumeMask.Except(Near);
        //        CellMask NegFar = LsTrk.GetLevelSetWing(0, -1).VolumeMask.Except(Near);

        //        LevelSet.Clear(PosFar);
        //        LevelSet.AccConstant(1, PosFar);
        //        LevelSet.Clear(NegFar);
        //        LevelSet.AccConstant(-1, NegFar);

        //        // project Fourier levelSet to DGfield on near field
        //        VolMask = Near;
        //    }

        //    LevelSet.Clear(VolMask);
        //    LevelSet.ProjectField(1.0,
        //        (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {

        //            int timestepIndex = (int)(time / (1.25e-5 / 2.0));
        //            //Console.WriteLine("timestepIndex {0}: H_exact = {1}", timestepIndex, H_exact[timestepIndex]);
        //            for (int nd = 0; nd < nodes.Lengths[0]; nd++) {
        //                results[nd] = - H_exact[timestepIndex] * Math.Sin(nodes[nd, 0] * (2 * Math.PI / DomainSize));
        //                //results[nd] -= H_exact[timestepIndex] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize))).Imaginary;
        //                results[nd] += nodes[nd, 1]; //+y
        //            }
        //        },
        //        new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));

        //}


        //public void ProjectPrescribedCurv(SinglePhaseField Curvature, out VectorField<SinglePhaseField> LevSetGradient, double time, CellMask VolMask = null) {

        //    if (VolMask == null)
        //        VolMask = CellMask.GetFullMask(Curvature.Basis.GridDat);

        //    Curvature.Clear();
        //    Curvature.ProjectField(1.0,
        //        (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {

        //            for (int nd = 0; nd < nodes.Lengths[0]; nd++) {
        //                if (time <= 0.0) {
        //                    results[nd] = H_exact[0] * Math.Pow((2 * Math.PI / DomainSize), 2) * Math.Sin(nodes[nd, 0] * (2 * Math.PI / DomainSize));
        //                    //results[nd] = H_exact[0] * Math.Pow((2 * Math.PI / DomainSize), 2) * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize))).Imaginary;
        //                } else {
        //                    int timestepIndex = (int)(time / 1.25e-5);
        //                    results[nd] = H_exact[timestepIndex] * Math.Pow((2 * Math.PI / DomainSize), 2) * Math.Sin(nodes[nd, 0] * (2 * Math.PI / DomainSize));
        //                    //results[nd] = H_exact[timestepIndex] * Math.Pow((2 * Math.PI / DomainSize), 2) * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize))).Imaginary;
        //                }
        //            }
        //        },
        //        new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));

        //    LevSetGradient = null;

        //}


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override ScalarFunction CurvEvaluation(int mode) {

            return (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {

                // kappa = -F''(x) / (F'(x)^2 + 1)^(3/2)
                Complex cF_x;
                double F_x;

                Complex cF_xx;
                double F_xx;

                for (int nd = 0; nd < nodes.Lengths[0]; nd++) {

                    F_x = 0.0;
                    F_xx = 0.0;

                    if (mode == -1) {

                        // F_x
                        if (cf_end == numFp / 2) {
                            cF_x = Complex.ImaginaryOne * DFT_coeff[numFp / 2] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * (numFp / 2)));
                            F_x = (2.0 * Math.PI / DomainSize) * (double)(numFp / 2) * cF_x.Real / (double)numFp;
                        }

                        for (int cf = 1; cf < cf_end; cf++) {
                            cF_x = Complex.ImaginaryOne * DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * cf));
                            F_x += (2 * Math.PI / DomainSize) * cf * 2.0 * cF_x.Real / (double)numFp;
                        }

                        // F_xx
                        if (cf_end == numFp / 2) {
                            cF_xx = -DFT_coeff[numFp / 2] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * (numFp / 2)));
                            F_xx = Math.Pow((2 * Math.PI / DomainSize) * (double)(numFp / 2), 2) * cF_xx.Real / (double)numFp;
                        }

                        for (int cf = 1; cf < cf_end; cf++) {
                            cF_xx = -DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * cf));
                            F_xx += Math.Pow((2 * Math.PI / DomainSize) * cf, 2) * 2.0 * cF_xx.Real / (double)numFp;
                        }

                    } else {

                        if (mode != 0) {

                            cF_x = Complex.ImaginaryOne * DFT_coeff[mode] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * mode));
                            F_x += (2 * Math.PI / DomainSize) * mode * cF_x.Real / (double)numFp;

                            cF_xx = -DFT_coeff[mode] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, nodes[nd, 0] * (2 * Math.PI / DomainSize) * mode));
                            F_xx += Math.Pow((2 * Math.PI / DomainSize) * mode, 2) * cF_xx.Real / (double)numFp;

                            if (mode != numFp / 2) {
                                F_x *= 2.0;
                                F_xx *= 2.0;
                            }
                        }

                    }

                    results[nd] = -F_xx / Math.Pow(Math.Pow(F_x, 2) + 1, 1.5);
                }

            };
        }


    }
}
