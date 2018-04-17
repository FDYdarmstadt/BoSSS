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
using BoSSS.Solution.Timestepping;

namespace BoSSS.Solution.LevelSetTools.FourierLevelSet {

    /// <summary>
    /// Encapsulation of Options for FourierLevSet
    /// </summary>
    public class FourierLevSetControl {

        /// <summary>
        /// planar (x,y) or cylindric (r, phi) representation of the Fourier coeff
        /// </summary>
        public FourierType FType;

        /// <summary>
        /// number of sample Points
        /// </summary>
        public int numSp;

        /// <summary>
        /// Size of the periodic domain
        /// </summary>
        public double DomainSize;

        /// <summary>
        /// minmal size of the grid
        /// </summary>
        public double h_min;

        /// <summary>
        /// periodic function for sample point generation
        /// </summary>
        public Func<double, double> PeriodicFunc;

        /// <summary>
        /// periodic Fourier points for the transformation
        /// </summary>
        public double[] FourierP;

        /// <summary>
        /// corresponding sample points for the Fourier points
        /// </summary>
        public double[] samplP;

        /// <summary>
        /// coordinates of the center for the PolarFourierLevSet
        /// </summary>
        public double[] center = new double[] { 0.0, 0.0 };

        /// <summary>
        /// If > -1: projection the given single mode 
        /// </summary>
        public int mode = -1;

        /// <summary>
        /// See <see cref="Interpolationtype"/>
        /// </summary>
        public Interpolationtype InterpolationType = Interpolationtype.LinearSplineInterpolation;

        /// <summary>
        /// See <see cref="FourierLevelSet_Timestepper"/>
        /// </summary>
        public FourierLevelSet_Timestepper Timestepper = FourierLevelSet_Timestepper.ExplicitEuler;

        /// <summary>
        /// underrelaxation factor for the coupled/implicit level-set evolution
        /// </summary>
        public double UnderRelax = 1.0;

        /// <summary>
        /// See <see cref="Fourier_Evolution"/>
        /// </summary>
        public Fourier_Evolution FourierEvolve = Fourier_Evolution.MaterialPoints;

        /// <summary>
        /// See <see cref="CenterMovement"/>
        /// </summary>
        public CenterMovement centerMove = CenterMovement.None;

        /// <summary>
        /// See <see cref="PolarLevSetForm"/>
        /// </summary>
        public PolarLevSetForm PLevSetForm = PolarLevSetForm.signedDist;

        /// <summary>
        /// curvature on the DG field extended from interface or not 
        /// </summary>
        public bool curvComp_extended = true;


        /// <summary>
        /// ctr
        /// </summary>
        /// <param name="FType">Fourier Type</param>
        /// <param name="numSp">Number of Sample Points</param>
        /// <param name="DomainSize">Size of the Periodic Domain</param>
        /// <param name="PeriodicFunc">Function to Project</param>
        /// <param name="h_min">minimal gridSize</param>
        public FourierLevSetControl(FourierType FType, int numSp, double DomainSize, Func<double, double> PeriodicFunc, double h_min) {
            this.FType = FType;
            this.numSp = numSp;
            this.DomainSize = DomainSize;
            this.PeriodicFunc = PeriodicFunc;
            this.h_min = h_min;

            Initialize();

        }

        /// <summary>
        /// ctr
        /// </summary>
        /// <param name="FType">Fourier Type</param>
        /// <param name="DomainSize">Size of the Periodic Domain</param>
        /// <param name="PeriodicFunc">Function to Project</param>
        /// <param name="r_min"></param>
        /// <param name="h_min">minimal gridSize</param>
        public FourierLevSetControl(FourierType FType, double DomainSize, Func<double, double> PeriodicFunc, double r_min, double h_min) {
            this.FType = FType;
            this.DomainSize = DomainSize;
            this.PeriodicFunc = PeriodicFunc;
            this.h_min = h_min;

            Initialize(r_min);

        }

        /// <summary>
        /// ctr
        /// </summary>
        /// <param name="FType"></param>
        /// <param name="DomainSize"></param>
        /// <param name="FourierP"></param>
        /// <param name="samplP"></param>
        /// <param name="h_min">minimal gridSize</param>
        public FourierLevSetControl(FourierType FType, double DomainSize, double[] FourierP, double[] samplP, double h_min) {
            this.FType = FType;
            this.DomainSize = DomainSize;
            this.FourierP = FourierP;
            this.samplP = samplP;
            this.numSp = samplP.Length;
            this.h_min = h_min;

        }


        /// <summary>
        /// ctr - you must asign the public variables manually
        /// </summary>
        public FourierLevSetControl() { }


        /// <summary>
        /// 
        /// </summary>
        public void Initialize(double r_min = 0.0) {

            if (r_min != 0.0) {
                int cf_end = (int)Math.Ceiling(2.0 * Math.PI * r_min / (2.0 * h_min));
                this.numSp = cf_end * 2;
            }

            this.FourierP = new double[numSp];
            this.samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (DomainSize / (double)numSp);
                samplP[sp] = PeriodicFunc(FourierP[sp]);
            }

        }

    }


    /// <summary>
    /// A Factory Class Creating a derived class of <see cref="FourierLevSetBase"/>
    /// </summary>
    public static class FourierLevelSetFactory {

        /// <summary>
        ///  Main Method of the FactoryClass for <see cref="FourierLevSetBase"/>, creating the instance specified in the <paramref name="Control"/>
        /// </summary>
        /// <param name="Control"></param>
        /// <returns></returns>
        public static FourierLevSetBase Build(FourierLevSetControl Control) {

            switch (Control.FType) {
                case FourierType.Planar:
                    return new PlanarFourierLevSet(Control);
                case FourierType.Polar:
                    return new PolarFourierLevSet(Control);
                default: throw new ArgumentException();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Control"></param>
        /// <returns></returns>
        public static FourierLevSetTimestepper Build_Timestepper(FourierLevSetControl Control, double[] currentState, DelComputeChangerate _DelCompChange, DelEvolveFourier _DelEvolveFourier) {

            RungeKuttaScheme RKscheme;
            switch (Control.Timestepper) {
                case FourierLevelSet_Timestepper.ExplicitEuler:
                    RKscheme = RungeKuttaScheme.ExplicitEuler;
                    break;
                case FourierLevelSet_Timestepper.Middlepoint:
                    RKscheme = RungeKuttaScheme.Middlepoint;
                    break;
                case FourierLevelSet_Timestepper.RungeKutta1901:
                    RKscheme = RungeKuttaScheme.RungeKutta1901;
                    break;
                case FourierLevelSet_Timestepper.ThreeOverEight:
                    RKscheme = RungeKuttaScheme.ThreeOverEight;
                    break;
                case FourierLevelSet_Timestepper.Heun:
                    RKscheme = RungeKuttaScheme.Heun2;
                    break;
                case FourierLevelSet_Timestepper.Heun3:
                    RKscheme = RungeKuttaScheme.Heun;
                    break;
                case FourierLevelSet_Timestepper.TVD3:
                    RKscheme = RungeKuttaScheme.TVD3;
                    break;
                case FourierLevelSet_Timestepper.RungeKutta3:
                    RKscheme = RungeKuttaScheme.RungeKutta3;
                    break;
                default:
                    throw new ArgumentException();
            }
            return new FourierLevSetTimestepper(Control, RKscheme, currentState, _DelCompChange, _DelEvolveFourier);

        }
    }


}
