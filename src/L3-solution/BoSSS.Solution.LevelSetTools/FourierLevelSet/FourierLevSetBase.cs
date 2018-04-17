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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.IO;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.IntegralTransforms.Algorithms;
using MathNet.Numerics.Interpolation.Algorithms;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Statistic;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.FourierLevelSet {


    /// <summary>
    /// specifies the coordoniate system for the Fourier-Series, e.g. y-x or cylinder
    /// </summary>
    public enum FourierType {

        /// <summary>
        /// 
        /// </summary>
        Planar,

        /// <summary>
        /// 
        /// </summary>
        Polar
    }

    /// <summary>
    /// Options for the interpolation on the Fourier points
    /// </summary>
    public enum Interpolationtype {

        /// <summary>
        /// 
        /// </summary>
        LinearSplineInterpolation,

        /// <summary>
        /// 
        /// </summary>
        CubicSplineInterpolation
    }

    /// <summary>
    /// options for the evolution of the Fourier Level-set
    /// </summary>
    public enum Fourier_Evolution {

        /// <summary>
        /// evolution by movement of the material sample points at the interface
        /// </summary>
        MaterialPoints,

        /// <summary>
        /// evolution by change rate of sample points at the static Fourier points
        /// </summary>
        FourierPoints,

        /// <summary>
        /// evolution by change rate of Fourier modes
        /// </summary>
        FourierModes,

    }


    public static class DiscreteFourierTransformation {

        /// <summary>
        /// class for the fast fourier transformation
        /// </summary>
        static DiscreteFourierTransform DFT = new DiscreteFourierTransform();

        /// <summary>
        /// Defines how to interpolate the evolved sample points onto equidistant spacing
        /// </summary>
        static Interpolationtype InterpolationType = Interpolationtype.LinearSplineInterpolation;

        /// <summary>
        /// Forward transformation for equidistant Fourier points
        /// </summary>
        /// <param name="samplFp"></param>
        /// <returns></returns>
        public static Complex[] TransformForward(double[] samplFp) {

            int numSFp = samplFp.Length;
            Complex[] invDFT_coeff = new Complex[numSFp];
            for (int sp = 0; sp < numSFp; sp++) {
                invDFT_coeff[sp] = (Complex)samplFp[sp];
            }
            Complex[] DFT_coeff = DFT.NaiveForward(invDFT_coeff, FourierOptions.Matlab);

            return DFT_coeff;
        }

        /// <summary>
        /// Forward transformation for nonequidistant sample points, will be interpolated on equidistant points on periodic domian
        /// </summary>
        /// <param name="samplP"></param>
        /// <param name="DomainSize">periodic domain</param>
        /// <param name="numSFp">number of interpolated sample points</param>
        /// <returns></returns>
        public static Complex[] TransformForward_nonequidistant(MultidimensionalArray samplP, double DomainSize, int numSFp = 0) {

            if (numSFp <= 0)
                numSFp = samplP.Lengths[0];

            double[] samplFp = new double[numSFp];
            InterpolateOntoFourierPoints(samplP, DomainSize, samplFp);

            Complex[] invDFT_coeff = new Complex[numSFp];
            for (int sp = 0; sp < numSFp; sp++) {
                invDFT_coeff[sp] = (Complex)samplFp[sp];
            }
            Complex[] DFT_coeff = DFT.NaiveForward(invDFT_coeff, FourierOptions.Matlab);

            return DFT_coeff;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="DFT_coeff"></param>
        /// <returns></returns>
        public static double[] TwoSidedPowerSpectrum(Complex[] DFT_coeff) {

            int numFp = DFT_coeff.Length;
            double[] P2 = DFT_coeff.Select(c => c.Magnitude / (double)numFp).ToArray();

            return P2;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="DFT_coeff"></param>
        /// <returns></returns>
        public static double[] SingleSidedPowerSpectrum(Complex[] DFT_coeff) {

            int numFp = DFT_coeff.Length;
            double[] P2 = DFT_coeff.Select(c => c.Magnitude / (double)numFp).ToArray();
            double[] P1 = new double[numFp / 2];
            P1[0] = P2[0];
            for (int i = 1; i < numFp / 2; i++) {
                P1[i] = 2.0 * P2[i];
            }

            return P1;

        }


        /// <summary>
        /// Interpolation of no equidistant sample point onto equidistant Fourier points
        /// </summary>
        public static void InterpolateOntoFourierPoints(MultidimensionalArray samplP, double DomainSize, double[] samplFp) {

            int numSp = samplP.Lengths[0];

            // set interpolation data (delete multiple independent values)
            ArrayList independentList = new ArrayList();
            ArrayList dependentList = new ArrayList();
            for (int sp = 0; sp < numSp; sp++) {
                if (independentList.Contains(samplP[sp, 0]) == false) {
                    independentList.Add(samplP[sp, 0]);
                    dependentList.Add(samplP[sp, 1]);
                }
            }
            // extend the interpolation data for sample points at the boundary of the domain
            independentList.Insert(0, samplP[numSp - 1, 0] - DomainSize);
            independentList.Insert(independentList.Count, samplP[0, 0] + DomainSize);
            dependentList.Insert(0, samplP[numSp - 1, 1]);
            dependentList.Insert(dependentList.Count, samplP[0, 1]);

            double[] independentVal = (double[])independentList.ToArray(typeof(double));
            double[] dependentVal = (double[])dependentList.ToArray(typeof(double));

            // set Fourier points
            int numSFp = samplFp.Length;
            double[] FourierP = new double[numSFp];
            for (int i = 0; i < numSFp; i++) {
                FourierP[i] = (DomainSize / numSFp) * i;
            }


            switch (InterpolationType) {
                case Interpolationtype.LinearSplineInterpolation:

                    LinearSplineInterpolation LinSpline = new LinearSplineInterpolation();
                    LinSpline.Initialize(independentVal, dependentVal);

                    for (int Fp = 0; Fp < numSFp; Fp++) {
                        samplFp[Fp] = LinSpline.Interpolate(FourierP[Fp]);
                    }

                    break;
                case Interpolationtype.CubicSplineInterpolation:

                    CubicSplineInterpolation CubSpline = new CubicSplineInterpolation();
                    CubSpline.Initialize(independentVal, dependentVal);

                    for (int Fp = 0; Fp < numSFp; Fp++) {
                        samplFp[Fp] = CubSpline.Interpolate(FourierP[Fp]);
                    }

                    break;
                default:
                    throw new NotImplementedException();
            }

        }

    }


    /// <summary>
    /// An explicit Interface description based on control points, which form a Fourier Series
    /// base class
    /// </summary>
    public abstract class FourierLevSetBase {

        readonly FourierLevSetControl control;

        /// <summary>
        /// class for the fast fourier transformation
        /// </summary>
        protected DiscreteFourierTransform DFT;

        /// <summary>
        /// Number of Fourier points for the fast Fourier transformation
        /// </summary>
        protected int numFp;

        /// <summary>
        /// cut-off for the projected Fourier modes
        /// </summary>
        protected int cf_end;

        /// <summary>
        /// If > -1: projection the given single mode 
        /// </summary>
        protected int mode;

        /// <summary>
        /// equidistant Fourier points, independent coordinate
        /// </summary>
        protected double[] FourierP;

        /// <summary>
        /// interpolated interface points at the Fourier points, dependent coordinate
        /// most recent iteration
        /// </summary>
        public double[] current_samplP;

        /// <summary>
        /// current 2-dimensional coordinates (x,y) of the interface points 
        /// </summary>
        public MultidimensionalArray current_interfaceP = new MultidimensionalArray(2);

        /// <summary>
        /// real-valued sample points
        /// </summary>
        protected Complex[] invDFT_coeff;

        /// <summary>
        /// complex Fourier-coefficients
        /// </summary>
        public Complex[] DFT_coeff;

        /// <summary>
        /// Size of the periodic domain of the independent variables
        /// </summary>
        protected double DomainSize;

        /// <summary>
        /// minimal gridSize: projected wavelengths are > 2*h_min
        /// </summary>
        protected double h_min;

        /// <summary>
        /// length of the material interface
        /// </summary>
        protected double interfaceLength;

        /// <summary>
        /// the resolution of the material interface
        /// </summary>
        protected double InterfaceResolution;

        /// <summary>
        /// Defines how to interpolate the evolved sample points onto equidistant spacing
        /// </summary>
        protected Interpolationtype InterpolationType;

        /// <summary>
        /// See <see cref="Fourier_Evolution"/>
        /// </summary>
        protected Fourier_Evolution FourierEvolve;


        /// <summary>
        /// Initiate the Series
        /// </summary>
        /// <param name="Control">Control Options for the FourierLevSet</param>
        protected FourierLevSetBase(FourierLevSetControl Control) {
            this.control = Control;

            this.FourierP = Control.FourierP;
            this.current_samplP = Control.samplP;
            this.numFp = FourierP.Length;
            this.DomainSize = Control.DomainSize;
            this.h_min = Control.h_min;
            this.mode = Control.mode;
            this.FourierEvolve = Control.FourierEvolve;

            // set up DFT
            invDFT_coeff = new Complex[numFp];
            for (int sp = 0; sp < numFp; sp++) {
                invDFT_coeff[sp] = (Complex)current_samplP[sp];
            }
            DFT = new DiscreteFourierTransform();
            DFT_coeff = DFT.NaiveForward(invDFT_coeff, FourierOptions.Matlab);

        }

        /// <summary>
        /// set the projecting Fourier modes for the Level set and curvature evaluation
        /// </summary>
        protected abstract void setProjectingFourierModes();


        /// <summary>
        /// computes the length of the samnple points polygon 
        /// </summary>
        protected abstract void setInterfaceLength(MultidimensionalArray interP);


        /// <summary>
        /// computes the changerate for the FourierLevelSet property according to <see cref="Fourier_Evolution"/>
        /// </summary>
        /// <param name="velocity"></param>
        /// <returns></returns>
        public abstract double[] ComputeChangerate(double dt, ConventionalDGField[] velocity, double[] current_FLSprop);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="current_FLSprop"></param>
        public abstract void EvolveFourierLS(ref double[] current_FLSprop);


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public virtual double[] GetFLSproperty() {

            switch (FourierEvolve) {
                case Fourier_Evolution.MaterialPoints:
                    double[] mpoints = new double[2 * numFp];
                    for (int sp = 0; sp < numFp; sp++) {
                        mpoints[sp * 2] = current_interfaceP[sp, 0];
                        mpoints[sp * 2 + 1] = current_interfaceP[sp, 1];
                    }
                    return mpoints;
                case Fourier_Evolution.FourierPoints:
                    return current_samplP;
                case Fourier_Evolution.FourierModes: {
                        double[] DFT_double = new double[2 * numFp];
                        for (int sp = 0; sp < numFp; sp++) {
                            DFT_double[sp * 2] = DFT_coeff[sp].Real;
                            DFT_double[sp * 2 + 1] = DFT_coeff[sp].Imaginary;
                        }
                        return DFT_double;
                    }
                default:
                    throw new ArgumentException();
            }
        }
 

        /// <summary>
        /// Since the Fourier representation is periodic, sample points, which moved outside the domain, have to be rearranged onto the domain
        /// </summary>
        protected void RearrangeOntoPeriodicDomain(MultidimensionalArray interP) {

            int numP = interP.Lengths[0];
            int[] rearranged = new int[numP];

            // periodic domain
            for (int sp = 0; sp < numP; sp++) {
                if (interP[sp, 0] >= DomainSize) {
                    interP[sp, 0] = interP[sp, 0] - DomainSize;
                    rearranged[sp] = 1;
                }
                if (interP[sp, 0] < 0) {
                    interP[sp, 0] = interP[sp, 0] + DomainSize;
                    rearranged[sp] = -1;
                }
            }

            int[] permutation = new int[numP];
            for (int ind = 0; ind < numP; ind++)
                permutation[ind] = ind;

            Array.Sort(permutation, (int i, int j) => Math.Sign(interP[i, 0] - interP[j, 0]));

            MultidimensionalArray interPcache = interP.CloneAs();

            for (int sp = 0; sp < numP; sp++) {
                interP[sp, 0] = interPcache[permutation[sp], 0];
                interP[sp, 1] = interPcache[permutation[sp], 1];
            }

        }


        /// <summary>
        /// Interpolation of the interface points onto the equidistant Fourier points
        /// </summary>
        protected void InterpolateOntoFourierPoints(MultidimensionalArray interP, double[] samplP) {

            int numP = interP.Lengths[0];

            // set interpolation data
            double[] independentVal = new double[numP + 2];
            double[] dependentVal = new double[numP + 2];
            for (int sp = 1; sp <= numP; sp++) {
                independentVal[sp] = interP[sp - 1, 0];
                dependentVal[sp] = interP[sp - 1, 1];
            }
            // extend the interpolation data for sample points at the boundary of the domain
            independentVal[0] = interP[numP - 1, 0] - DomainSize;
            dependentVal[0] = interP[numP - 1, 1];
            independentVal[numP + 1] = interP[0, 0] + DomainSize;
            dependentVal[numP + 1] = interP[0, 1];

            switch (this.InterpolationType) {
                case Interpolationtype.LinearSplineInterpolation:

                    LinearSplineInterpolation LinSpline = new LinearSplineInterpolation();
                    LinSpline.Initialize(independentVal, dependentVal);

                    for (int sp = 0; sp < numFp; sp++) {
                        samplP[sp] = LinSpline.Interpolate(FourierP[sp]);
                        //invDFT_coeff[sp] = (Complex)samplP[sp];
                    }

                    break;
                case Interpolationtype.CubicSplineInterpolation:

                    CubicSplineInterpolation CubSpline = new CubicSplineInterpolation();
                    CubSpline.Initialize(independentVal, dependentVal);

                    for (int sp = 0; sp < numFp; sp++) {
                        samplP[sp] = CubSpline.Interpolate(FourierP[sp]);
                        //invDFT_coeff[sp] = (Complex)samplP[sp];
                    }

                    break;
                default:
                    throw new NotImplementedException();
            }

        }


        /// <summary>
        /// smoothing of the Interface by limiting the evaluated DFT-coefficients
        /// </summary>
        protected void SmoothSamplePoints() {

            Complex cmplx_res;
            // smooth interpolated current sample points
            for (int sp = 0; sp < numFp; sp++) {

                current_samplP[sp] = DFT_coeff[0].Real;
                if (cf_end == numFp / 2)
                    current_samplP[sp] += Complex.Multiply(DFT_coeff[numFp / 2],
                                                           Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * (2 * Math.PI / DomainSize) * (numFp / 2)))
                                                           ).Real;

                for (int cf = 1; cf < cf_end; cf++) {
                    cmplx_res = DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * (2 * Math.PI / DomainSize) * cf));
                    current_samplP[sp] += 2.0 * cmplx_res.Real;
                }
                current_samplP[sp] /= (double)numFp;
            }

        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="interP"></param>
        /// <param name="p"></param>
        /// <param name="dir"></param>
        /// <returns></returns>
        protected double getDomainDistance(MultidimensionalArray interP, int p, int dir) {

            int numP = interP.Lengths[0];
            if (p < 0 || p >= numP)
                throw new ArgumentException();

            double domain_shift = DomainSize;
            if (this is PolarFourierLevSet)
                domain_shift = 0.0;

            double dist = 0.0;
            if (p + dir >= 0 && p + dir <= numP - 1) {
                dist = Math.Sqrt((interP[p + dir, 0] - interP[p, 0]).Pow2() + (interP[p + dir, 1] - interP[p, 1]).Pow2());
            } else {
                if (p + dir >= numP)
                    dist = Math.Sqrt((interP[p + dir - numP, 0] + domain_shift - interP[p, 0]).Pow2() + (interP[p + dir - numP, 1] - interP[p, 1]).Pow2());
                if (p + dir < 0)
                    dist = Math.Sqrt((interP[p + dir + numP, 0] + domain_shift - interP[p, 0]).Pow2() + (interP[p + dir + numP, 1] - interP[p, 1]).Pow2());
            }

            return dist;
        }


        /// <summary>
        /// Checks the minimal/maximal resolution of the material sample points and add/remove sample points if necessary
        /// </summary>
        protected bool Reseeding() {

            bool reseed = false;

            double domain_shift = DomainSize;
            if (this is PolarFourierLevSet)
                domain_shift = 0.0;

            int sp = 0;
            int numSp = current_interfaceP.Lengths[0];
            while (sp < numSp) {
                // determine distance between two successive points (independent var in domainsize)
                double dist = getDomainDistance(current_interfaceP, sp, +1);

                if (dist <= 2.0 * InterfaceResolution && dist >= InterfaceResolution / 2.0) {
                    sp += 1;
                } else {

                    // check for adding
                    if (dist > 2.0 * InterfaceResolution) {

                        Console.WriteLine("Reseeding: Adding after sp = {0}", sp);
                        reseed = true;

                        // add to current_interfaceP
                        MultidimensionalArray interP_Cache = current_interfaceP.CloneAs();
                        current_interfaceP = MultidimensionalArray.Create(numSp + 1, 2);
                        int dp = 0;
                        for (int p = 0; p < numSp + 1; p++) {
                            if (p == sp + 1) {
                                // interpolate the new sample point
                                if (p == numSp) {
                                    current_interfaceP[p, 0] = (interP_Cache[sp, 0] + domain_shift + interP_Cache[0, 0]) / 2;
                                    current_interfaceP[p, 1] = (interP_Cache[sp, 1] + interP_Cache[0, 1]) / 2;
                                } else {
                                    current_interfaceP[p, 0] = (interP_Cache[sp, 0] + interP_Cache[sp + 1, 0]) / 2;
                                    current_interfaceP[p, 1] = (interP_Cache[sp, 1] + interP_Cache[sp + 1, 1]) / 2;
                                }
                                dp = -1;
                            } else {
                                current_interfaceP[p, 0] = interP_Cache[p + dp, 0];
                                current_interfaceP[p, 1] = interP_Cache[p + dp, 1];
                            }
                        }
                        numSp = current_interfaceP.Lengths[0];

                        // add to interfaceP (History)
                        //History_MultiArray interP_historyCache = interfaceP.CloneAs();
                        //interfaceP.Clear();
                        //for (int i = 0; i <= interfaceP.m_capacity; i++) {
                        //    MultidimensionalArray item = MultidimensionalArray.Create(current_interfaceP.Lengths);
                        //    dp = 0;
                        //    for (int p = 0; p < numSp; p++) {
                        //        if (p == sp + 1) {
                        //            // interpolate the new sample point
                        //            if (p == numSp - 1) {
                        //                item[p, 0] = ((interP_historyCache[-i])[sp, 0] + domain_shift + (interP_historyCache[-i])[0, 0]) / 2;
                        //                item[p, 1] = ((interP_historyCache[-i])[sp, 1] + (interP_historyCache[-i])[0, 1]) / 2;
                        //            } else {
                        //                item[p, 0] = ((interP_historyCache[-i])[sp, 0] + (interP_historyCache[-i])[sp + 1, 0]) / 2;
                        //                item[p, 1] = ((interP_historyCache[-i])[sp, 1] + (interP_historyCache[-i])[sp + 1, 1]) / 2;
                        //            }
                        //            dp = -1;
                        //        } else {
                        //            item[p, 0] = (interP_historyCache[-i])[p + dp, 0];
                        //            item[p, 1] = (interP_historyCache[-i])[p + dp, 1];
                        //        }
                        //    }
                        //    interfaceP.setHistoryEntry(-i, item);
                        //}

                        // velocity hist
                        //History_MultiArray velAtInter_historyCache = velocityAtInterfaceP.CloneAs();
                        //velocityAtInterfaceP.Clear();
                        //for (int i = 0; i <= velocityAtInterfaceP.m_capacity; i++) {
                        //    MultidimensionalArray item = MultidimensionalArray.Create(current_interfaceP.Lengths);
                        //    dp = 0;
                        //    for (int p = 0; p < numSp; p++) {
                        //        if (p == sp + 1) {
                        //            // interpolate the new sample point
                        //            if (p == numSp - 1) {
                        //                item[p, 1] = ((velAtInter_historyCache[-i])[sp, 0] + (velAtInter_historyCache[-i])[0, 0]) / 2;
                        //                item[p, 1] = ((velAtInter_historyCache[-i])[sp, 1] + (velAtInter_historyCache[-i])[0, 1]) / 2;
                        //            } else {
                        //                item[p, 1] = ((velAtInter_historyCache[-i])[sp, 0] + (velAtInter_historyCache[-i])[sp + 1, 0]) / 2;
                        //                item[p, 1] = ((velAtInter_historyCache[-i])[sp, 1] + (velAtInter_historyCache[-i])[sp + 1, 1]) / 2;
                        //            }
                        //            dp = -1;
                        //        } else {
                        //            item[p, 0] = (velAtInter_historyCache[-i])[p + dp, 0];
                        //            item[p, 1] = (velAtInter_historyCache[-i])[p + dp, 1];
                        //        }
                        //    }
                        //    velocityAtInterfaceP.setHistoryEntry(-i, item);
                        //}

                        // set next sample point to check
                        sp += 2;
                    }

                    // check for removing
                    if (dist < (InterfaceResolution / 2.0)) {

                        // check which points to remove
                        double dist_left = getDomainDistance(current_interfaceP, sp, -1);
                        double dist_right = getDomainDistance(current_interfaceP, sp + 1, +1);
                        int pdel = 0;
                        if (dist_left <= dist_right) {
                            // remove current sample point sp
                            pdel = sp;
                        } else {
                            // remove the next one sp+1
                            pdel = sp + 1;
                        }
                        Console.WriteLine("Reseeding: Removing sp = {0}", pdel);
                        reseed = true;

                        MultidimensionalArray interP_Cache = current_interfaceP.CloneAs();
                        current_interfaceP = MultidimensionalArray.Create(numSp - 1, 2);
                        int dp = 0;
                        for (int p = 0; p < numSp - 1; p++) {
                            if (p == pdel) {
                                dp = +1;
                            }
                            current_interfaceP[p, 0] = interP_Cache[p + dp, 0];
                            current_interfaceP[p, 1] = interP_Cache[p + dp, 1];
                        }
                        numSp = current_interfaceP.Lengths[0];

                        // remove in interfaceP (History)
                        //History_MultiArray interP_historyCache = interfaceP.CloneAs();
                        //interfaceP.Clear();
                        //for (int i = 0; i <= interfaceP.m_capacity; i++) {
                        //    MultidimensionalArray item = MultidimensionalArray.Create(current_interfaceP.Lengths);
                        //    dp = 0;
                        //    for (int p = 0; p < numSp; p++) {
                        //        if (p == pdel) {
                        //            dp = +1;
                        //        }
                        //        item[p, 0] = (interP_historyCache[-i])[p + dp, 0];
                        //        item[p, 1] = (interP_historyCache[-i])[p + dp, 1];
                        //    }
                        //    interfaceP.setHistoryEntry(-i, item);
                        //}

                        // velocity hist
                        //History_MultiArray velAtInter_historyCache = velocityAtInterfaceP.CloneAs();
                        //velocityAtInterfaceP.Clear();
                        //for (int i = 0; i <= velocityAtInterfaceP.m_capacity; i++) {
                        //    MultidimensionalArray item = MultidimensionalArray.Create(current_interfaceP.Lengths);
                        //    dp = 0;
                        //    for (int p = 0; p < numSp; p++) {
                        //        if (p == pdel) {
                        //            dp = +1;
                        //        }
                        //        item[p, 0] = (velAtInter_historyCache[-i])[p + dp, 0];
                        //        item[p, 1] = (velAtInter_historyCache[-i])[p + dp, 1];
                        //    }
                        //    velocityAtInterfaceP.setHistoryEntry(-i, item);
                        //}

                    }
                }

            }
            return reseed;

        }


        /// <summary>
        /// Computes the DFT-coefficients for the current sample points 
        /// </summary>
        protected void setDFTcoeff() {
            for (int sp = 0; sp < numFp; sp++) {
                invDFT_coeff[sp] = (Complex)current_samplP[sp];
            }
            DFT_coeff = DFT.NaiveForward(invDFT_coeff, FourierOptions.Matlab);
        }


        /// <summary>
        /// Calculate A Level-Set field from the Explicit description
        /// </summary>
        /// <param name="LevelSet">Target Field</param>
        /// <param name="LsTrk"></param>
        public void ProjectToDGLevelSet(SinglePhaseField LevelSet, LevelSetTracker LsTrk = null) {

            CellMask VolMask;
            if (LsTrk == null) {
                VolMask = CellMask.GetFullMask(LevelSet.Basis.GridDat);
            } else {

                // set values in positive and negative FAR region to +1 and -1
                CellMask Near = LsTrk.Regions.GetNearMask4LevSet(0, 1);
                CellMask PosFar = LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask.Except(Near);
                CellMask NegFar = LsTrk.Regions.GetLevelSetWing(0, -1).VolumeMask.Except(Near);

                LevelSet.Clear(PosFar);
                LevelSet.AccConstant(1, PosFar);
                LevelSet.Clear(NegFar);
                LevelSet.AccConstant(-1, NegFar);

                // project Fourier levelSet to DGfield on near field
                VolMask = Near;
            }

            LevelSet.Clear(VolMask);
            // scalar function is already vectorized for parallel execution
            // nodes in global coordinates
            LevelSet.ProjectField(1.0,
                PhiEvaluation(mode),
                new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));

            // check the projection error
            projErr_phiDG = LevelSet.L2Error(
                PhiEvaluation(mode),
                new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));
            //if (projErr_phiDG >= 1e-5)
            //    Console.WriteLine("WARNING: LevelSet projection error onto PhiDG = {0}", projErr_phiDG);


            // project on higher degree field and take the difference
            //SinglePhaseField higherLevSet = new SinglePhaseField(new Basis(LevelSet.GridDat, LevelSet.Basis.Degree * 2), "higherLevSet");
            //higherLevSet.ProjectField(1.0, PhiEvaluator(), new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));
            //double higherProjErr = higherLevSet.L2Error(
            //    PhiEvaluator(),
            //    new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));
            //Console.WriteLine("LevelSet projection error onto higherPhiDG = {0}", higherProjErr.ToString());


            // check the projection from current sample points on the DGfield
            //MultidimensionalArray interP = MultidimensionalArray.Create(current_interfaceP.Lengths);
            //if (this is PolarFourierLevSet) {
            //    interP = ((PolarFourierLevSet)this).interfaceP_cartesian;
            //} else {
            //    interP = current_interfaceP;
            //}

            //projErr_interface = InterfaceProjectionError(LevelSet, current_interfaceP);
            //if (projErr_interface >= 1e-3)
            //    Console.WriteLine("WARNING: Interface projection Error onto PhiDG = {0}", projErr_interface);


        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        abstract public ScalarFunction PhiEvaluation(int mode);

        /// <summary>
        /// computes the L2-error of the interface points onto the LevelSet-DGfield
        /// </summary>
        /// <param name="LevelSet"></param>
        /// <param name="interP"></param>
        public double InterfaceProjectionError(SinglePhaseField LevelSet, MultidimensionalArray interP) {
            GridData grdat = (GridData)LevelSet.GridDat;
            FieldEvaluation fEval = new FieldEvaluation(grdat);

            MultidimensionalArray PhiAtSamplePoints = MultidimensionalArray.Create(interP.Lengths[0], 1);

            DGField[] DGLevSet = new DGField[] { LevelSet };
            int outP = fEval.Evaluate(1, DGLevSet, interP, 0, PhiAtSamplePoints);

            if (outP != 0)
                throw new Exception("points outside the grid for fieldevaluation");

            return PhiAtSamplePoints.L2Norm();
        }


        /// <summary>
        /// calculate the curvature corresponding to the Level-Set
        /// </summary>
        /// <param name="Curvature"></param>
        /// <param name="LevSetGradient"></param>
        /// <param name="VolMask"></param>
        public void ProjectToDGCurvature(SinglePhaseField Curvature, out VectorField<SinglePhaseField> LevSetGradient, CellMask VolMask = null) {

            if (VolMask == null)
                VolMask = CellMask.GetFullMask(Curvature.Basis.GridDat);

            Curvature.Clear();
            Curvature.ProjectField(1.0,
                CurvEvaluation(mode),
                new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));

            // check the projection error
            projErr_curv = Curvature.L2Error(
                CurvEvaluation(mode),
                new Foundation.Quadrature.CellQuadratureScheme(true, VolMask));
            //if (projErr_curv >= 1e-5)
            //    Console.WriteLine("WARNING: Curvature projection error onto PhiDG = {0}", projErr_curv);


            LevSetGradient = null;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        abstract public ScalarFunction CurvEvaluation(int mode);


        public double projErr_phiDG;

        public double projErr_interface;

        public double projErr_curv;


        /// <summary>
        /// L2-norm of the distance between current sampls points and the ones from last iteration (for InnerLoop)
        /// </summary>
        /// <returns></returns>
        //public double FourierResidual() {
        //    MultidimensionalArray InterfaceP_Residual = current_interfaceP.CloneAs();
        //    InterfaceP_Residual.Acc(-1.0, interfaceP[0]);

        //    return InterfaceP_Residual.L2Norm();
        //}


        /// <summary>
        /// returns the interpolated sample points as a string for the log file entry
        /// </summary>
        /// <returns></returns>
        protected virtual string samplP_LogFormat() {
            string samplP_DB = "";
            for (int sp = 0; sp < numFp; sp++) {
                samplP_DB = samplP_DB + "\t" + current_samplP[sp].ToString();
            }
            return samplP_DB;
        }

        /// <summary>
        /// returns the material interface points as a string for the log file entry
        /// </summary>
        /// <returns></returns>
        protected virtual string interfaceP_LogFormat() {
            string interP_DB = "";
            int numP = current_interfaceP.Lengths[0];
            for (int ip = 0; ip < numP; ip++) {
                interP_DB = interP_DB + "\t" + current_interfaceP[ip, 0].ToString();
                interP_DB = interP_DB + "\t" + current_interfaceP[ip, 1].ToString();
            }
            return interP_DB;
        }

        /// <summary>
        /// returns the projections erros as a string for the log file entry
        /// </summary>
        /// <returns></returns>
        private string projErr_LogFormat() {
            string projErr_DB = String.Format("\t{0}\t{1}\t{2}", projErr_phiDG, projErr_interface, projErr_curv);
            return projErr_DB;
        }

        /// <summary>
        /// returns necessary values for restarting the Fourier Level set
        /// </summary>
        /// <returns></returns>
        public virtual double[] getRestartInfo() {
            return current_samplP;
        }

        /// <summary>
        /// saves L2-error of the Fourier series onto the DG-Fields
        /// </summary>
        TextWriter Log_projErrFLS;

        /// <summary>
        /// saves interpolated sample points
        /// </summary>
        TextWriter Log_samplP;

        /// <summary>
        /// saves material interface points
        /// </summary>
        TextWriter Log_interP;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="DbDriver"></param>
        /// <param name="sessInfo"></param>
        /// <param name="timestep"></param>
        /// <param name="phystime"></param>
        public void setUpLogFiles(IDatabaseDriver DbDriver, SessionInfo sessInfo, int timestep, double phystime) {

            string timeHeader = String.Format("{0}\t{1}", "#timestep", "#time");
            string timeEntry = String.Format("{0}\t{1}", timestep, phystime);

            string firstline = null;
            string line = null;

            // l2 projection error
            //Log_projErrFLS = DbDriver.FsDriver.GetNewLog("Log_projErrFLS", sessInfo.ID);
            //string firstline = String.Format("{0}\t{1}\t{2}\t{3}", timeHeader, "#L2 error PhiDG", "#L2 error interface", "#L2 error curvature");
            //Log_projErrFLS.WriteLine(firstline);
            //string line = timeEntry + projErr_LogFormat();
            //Log_projErrFLS.WriteLine(line);
            //Log_projErrFLS.Flush();

            // interpolated sample points
            Log_samplP = DbDriver.FsDriver.GetNewLog("Log_samplP", sessInfo.ID);
            if (this is PolarFourierLevSet)
                firstline = String.Format("{0}\t{1}\t{2}", timeHeader, "center", "#samplP");
            else
                firstline = String.Format("{0}\t{1}", timeHeader, "#samplP");
            Log_samplP.WriteLine(firstline);
            line = timeEntry + samplP_LogFormat();
            Log_samplP.WriteLine(line);
            Log_samplP.Flush();

            // material interface points
            Log_interP = DbDriver.FsDriver.GetNewLog("Log_interP", sessInfo.ID);
            firstline = String.Format("{0}\t{1}", timeHeader, "#interP");
            Log_interP.WriteLine(firstline);
            line = timeEntry + interfaceP_LogFormat();
            Log_interP.WriteLine(line);
            Log_interP.Flush();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="TimestepNo"></param>
        /// <param name="phystime"></param>
        public void saveToLogFiles(int TimestepNo, double phystime) {

            string timeEntry = String.Format("{0}\t{1}", TimestepNo, phystime);
            //if (Log_projErrFLS != null) {
            //    string line = timeEntry + projErr_LogFormat();
            //    Log_projErrFLS.WriteLine(line);
            //    Log_projErrFLS.Flush();
            //}

            if (Log_samplP != null) {
                string line = timeEntry + samplP_LogFormat();
                Log_samplP.WriteLine(line);
                Log_samplP.Flush();
            }

            if (Log_interP != null) {
                string line = timeEntry + interfaceP_LogFormat();
                Log_interP.WriteLine(line);
                Log_interP.Flush();
            }
        }


    }




}
