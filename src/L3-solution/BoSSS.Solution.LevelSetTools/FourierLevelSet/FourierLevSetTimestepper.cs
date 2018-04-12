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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Solution.LevelSetTools.FourierLevelSet {


    /// <summary>
    /// options for the time-discretization of the Fourier level-set evolution
    /// </summary>
    public enum FourierLevelSet_Timestepper {

        /// <summary>
        /// Explicit Euler - Rule, Order 1, 1 Stage;
        /// </summary>
        ExplicitEuler,

        /// <summary>
        /// Middlepoint - Rule, Order 2, 2 Stages;
        /// </summary>  
        Middlepoint,

        /// <summary>
        /// classical Method of Runge and Kutta, anno 1901, order 4, 4 stages;
        /// </summary>
        RungeKutta1901,

        /// <summary>
        /// 3/8 - Rule, order 4, 4 stages
        /// </summary>
        ThreeOverEight,

        /// <summary>
        /// Heun - Rule, order 2, 2 stages
        /// </summary>
        Heun,

        /// <summary>
        /// Heun - Rule, order 3, 3 stages
        /// </summary>
        Heun3,

        /// <summary>
        /// Runge Kutta - Rule, order 3, 3 stages
        /// </summary>
        RungeKutta3,

        /// <summary>
        /// 
        /// </summary>
        TVD3


    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="dt"></param>
    /// <param name="velocity"></param>
    /// <returns></returns>
    public delegate double[] DelComputeChangerate(double dt, ConventionalDGField[] velocity, double[] current_FLSprop);

    /// <summary>
    /// 
    /// </summary>
    /// <param name="current_FLSproperty"></param>
    public delegate void DelEvolveFourier(ref double[] current_FLSproperty);


    /// <summary>
    /// base class for all Fourier Timesteppers. implements standard Explicit Euler scheme
    /// </summary>
    public class FourierLevSetTimestepper {


        internal double[] current_FLSproperty;

        internal double[] previous_FLSproperty;

        /// <summary>
        /// See <see cref="FourierLevelSet_Timestepper"/>
        /// </summary>
        protected FourierLevelSet_Timestepper Timestepper;

        internal RungeKuttaScheme RKscheme;

        internal double[] changerateAtStage;

        /// <summary>
        /// underrelaxation factor for the Evolution
        /// </summary>
        protected double underrelaxation;

        protected DelComputeChangerate DelCompChange;

        protected DelEvolveFourier DelEvolveFourier;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Control"></param>
        /// <param name="currentState"></param>
        public FourierLevSetTimestepper(FourierLevSetControl Control, RungeKuttaScheme _RKscheme, double[] currentState, DelComputeChangerate _DelCompChange, DelEvolveFourier _DelEvolveFourier) {

            this.Timestepper = Control.Timestepper;
            this.RKscheme = _RKscheme;
            this.underrelaxation = Control.UnderRelax;
            this.DelCompChange = _DelCompChange;
            this.DelEvolveFourier = _DelEvolveFourier;

            current_FLSproperty = currentState;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="velocity"></param>
        public virtual void moveLevelSet(double dt, ConventionalDGField[] velocity) {

            int s = RKscheme.Stages;

            double[] FLSpropertyAtStage;
            ArrayList stage_changerates = new ArrayList();
            for (int j = 0; j < s; j++) {
                FLSpropertyAtStage = previous_FLSproperty.CloneAs();
                for (int l = 0; l < s; l++) {
                    if (RKscheme.a[j, l] != 0.0)
                        FLSpropertyAtStage.AccV(dt * RKscheme.a[j, l], (double[])stage_changerates[l]);
                }
                // compute the change rate at current stage
                changerateAtStage = DelCompChange(dt, velocity, FLSpropertyAtStage);
                stage_changerates.Insert(j, changerateAtStage);

            }

            current_FLSproperty = previous_FLSproperty.CloneAs();
            for (int j = 0; j < s; j++) {
                current_FLSproperty.AccV(dt * RKscheme.b[j], (double[])stage_changerates[j]);
            }

            DelEvolveFourier(ref current_FLSproperty);

            //// compute the current change rate at current state
            //current_changerate = DelCompChange(dt, velocity);

            //// Move FourierLevelSet property
            //double[] FLSproperty_evo = current_FLSproperty.CloneAs();
            //FLSproperty_evo.AccV(dt, current_changerate);

            //// Underrelaxation
            //current_FLSproperty.ScaleV(1.0 - underrelaxation);
            //current_FLSproperty.AccV(underrelaxation, FLSproperty_evo);

            //DelEvolveFourier(ref current_FLSproperty);

        }

        /// <summary>
        /// update the history of the interface points after convergence
        /// </summary>
        public virtual void updateFourierLevSet() {

            previous_FLSproperty = current_FLSproperty.CloneAs();

        }

    }

}
