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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using System.Globalization;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;


namespace BoSSS.Application.LsTest {

    /// <summary>
    /// Basic test class for the level-set unit tests
    /// </summary>
    abstract class LevelSetBaseTest : ILevelSetTest {

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetBaseTest(int spatDim, int LevelSetDegree) {
            this.SpatialDimension = spatDim;
            this.LevelsetPolynomialDegree = LevelSetDegree;
        }


        public double dt {
            get {
                return -1.0;    // will be set in LevelSetTest() according to level set cfl 
            }
        }


        /// <summary>
        /// computes the timestep size according to the level-set CFL condition
        /// </summary>
        /// <param name="Resolution"></param>
        /// <param name="LSdegree"></param>
        /// <returns></returns>
        public abstract double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel, int temporalResolution);

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public abstract double getEndTime();

       
        /// <summary>
        /// creates grid according to set resolution
        /// </summary>
        /// <param name="Resolution"></param>
        /// <returns></returns>
        public abstract GridCommons CreateGrid(int Resolution);

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public abstract IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig();

        /// <summary>
        /// Level-Set function
        /// </summary>
        public abstract Func<double[], double, double>[] GetPhi();

        /// <summary>
        /// velocity field
        /// </summary>
        public abstract Func<double[], double, double>[][] GetU();        

        public int LevelsetPolynomialDegree {
            get;
            private set;
        }

        public virtual double[,] AcceptableError {
            get;
        }


        /// <summary>
        /// returns spatial dimension
        /// </summary>
        public int SpatialDimension {
            get;
            private set;
        }

        public abstract int NoOfLevelsets { get; }
    }

}
