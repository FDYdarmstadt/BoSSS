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
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XdgNastyLevsetLocationTest {    
    
    interface ITest {
        GridCommons GetGrid();

        bool NextTestCase();

        void ResetTest();
        
        /// <summary>
        /// 
        /// </summary>
        bool VolumeTestSupported {
            get;
        }

        /// <summary>
        /// 
        /// </summary>
        bool EdgeTestSupported {
            get;
        }

        /// <summary>
        /// 
        /// </summary>
        bool LevelsetTestSupported {
            get;
        }

        /// <summary>
        /// 
        /// </summary>
        bool BoundaryPlusLevelsetTestSupported {
            get;
        }


        double GetLevelSet(double x, double y);

        /// <summary>
        /// returns the volume of cell <paramref name="jcell"/> for species <paramref name="species"/>
        /// </summary>
        /// <remarks>
        /// implementation is only required if <see cref="VolumeTestSupported"/> is true
        /// </remarks>
        double CellVolume(int jCell, string species, GridData g);
        
        /// <summary>
        /// returns the area edge <paramref name="jcell"/> for species <paramref name="species"/>
        /// </summary>
        /// <remarks>
        /// implementation is only required if <see cref="VolumeTestSupported"/> is true
        /// </remarks>
        double EdgeArea(int iEdge, string species, GridData g);

        /// <summary>
        /// returns the area of the level-set 
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        double LevelsetArea(int jCell, GridData g);
        
        /// <summary>
        /// returns the area of cell-boundary <paramref name="jcell"/> for species <paramref name="species"/>
        /// plus the area of the level-set surface.
        /// </summary>
        /// <remarks>
        /// implementation is only required if <see cref="BoundaryPlusLevelsetTestSupported"/> is true
        /// </remarks>
        double CellBoundaryPlusLevelsetArea(int jCell, string species, GridData g);
    }
}
