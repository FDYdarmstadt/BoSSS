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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.LevelSetTools {

    /// <summary>
    /// Interface for Level-Set Motion Algorithms
    /// </summary>
    public interface ILevelSetAdvection {

        /// <summary>
        /// Move Level-Set with Timestep <paramref name="dt"/>
        /// </summary>
        /// <param name="dt">Timestep</param>
        void Advect(double dt);



        void FinishTimeStep();



    }

    /// <summary>
    /// Default implementation : Do nothing
    /// </summary>
    public abstract class LevelSetAdvection : ILevelSetAdvection {

        // A shared instance that can be used for comparisons
        public static readonly ILevelSetAdvection Null = new NullLevelSetAdvection();      


        private class NullLevelSetAdvection : ILevelSetAdvection {
            /// <summary>
            /// Do not move the Level-set
            /// </summary>
            /// <param name="dt">unnecessary, since no motion</param>
            public void Advect(double dt = 0) {
                /// Do nothing
            }

            public void FinishTimeStep() {
                /// Do nothing
            }
        }

        public abstract void Advect(double dt);
        public abstract void FinishTimeStep();

    }

    public class PrescribedMover : ILevelSetAdvection {

        public double time = 0.0;
        Func<double[], double, double> MotionFunc;

        DGField LevelSet;

        public PrescribedMover(Func<double[], double, double> MotionFunc, DGField LevelSet) {
            this.MotionFunc = MotionFunc;
            this.LevelSet = LevelSet;
        }

        public void Advect(double dt) {
            time += dt;
            LevelSet.ProjectField(MotionFunc.Vectorize(time));
        }

        public void FinishTimeStep() {
            /// Do nothing;
        }
    }

}
