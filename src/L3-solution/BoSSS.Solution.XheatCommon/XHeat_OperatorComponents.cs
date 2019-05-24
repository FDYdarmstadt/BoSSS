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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Solution.XheatCommon {

    /// <summary>
    /// 
    /// </summary>
    public static partial class XOperatorComponentsFactory {

        //==============
        // Heat equation
        //==============

        public static void AddSpeciesHeatEq(XSpatialOperatorMk2 XOp) {

        }


        public static void AddInterfaceHeatEq(XSpatialOperatorMk2 XOp) {

        }

    }


    /// <summary>
    /// configuration options for the bulk NSE including continuity equation
    /// </summary>
    public interface IHeat_Configuration : ISolver_Configuration {

        /// <summary>
        /// thermal parameters
        /// </summary>
        ThermalParameters GetThermParams { get; }

        /// <summary>
        /// include transport operator
        /// </summary>
        bool isTransport { get; }

    }


    /// <summary>
    /// extended configuration options for interface discretizations
    /// </summary>
    public interface IXHeat_Configuration : IHeat_Configuration {

        /// <summary>
        /// switch for moving mesh flux discretizations
        /// </summary>
        bool isMovingMesh { get; }

        /// <summary>
        /// true if the interface is a material interface
        /// </summary>
        bool isMatInt { get; }

    }


}
