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

using BoSSS.Solution;
using NUnit.Framework;
using System;

namespace CNS.Tests.IsentropicVortex {

    /// <summary>
    /// Tests using the example of an isentropic vortex in a uniform
    /// background flow
    /// </summary>
    public class IsentropicVortexTest : TestProgram<VortexControl> {

        //public static void Main(string[] args) {
        //    IsentropicVortexIdealGasRusanovTest();
        //    //IsentropicVortexCovolumeGasRusanovTest();
        //}

        /// <summary>
        /// Tests the CNS solver with the <see cref="Convection.RusanovFlux"/>
        /// using the example moving isentropic vortex in an ideal gas.
        /// </summary>
        [Test]
        public static void IsentropicVortexIdealGasRusanovTest() {
            Program<VortexControl> p = null;
            Application<VortexControl>._Main(
                new string[] { @"-c cs:CNS.Tests.IsentropicVortex.ControlFiles.IsentropicVortexIdealGasRusanov()" },
                false,
                delegate() {
                    p = new IsentropicVortexTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 5e-3),
                Tuple.Create("L2ErrorPressure", 6e-3),
                Tuple.Create("L2ErrorEntropy", 4e-3));
        }


        /// <summary>
        /// Tests the CNS solver with the <see cref="Convection.HLLCFlux"/>
        /// using the example moving isentropic vortex in an ideal gas.
        /// </summary>
        [Test]
        public static void IsentropicVortexIdealGasHLLCTest() {
            Program<VortexControl> p = null;
            Application<VortexControl>._Main(
                new string[] { @"-c cs:CNS.Tests.IsentropicVortex.ControlFiles.IsentropicVortexIdealGasHLLC()" },
                false,
                delegate() {
                    p = new IsentropicVortexTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 4e-3),
                Tuple.Create("L2ErrorPressure", 5e-3),
                Tuple.Create("L2ErrorEntropy", 3.15e-3));
        }

        /// <summary>
        /// Tests the CNS solver with the <see cref="Convection.HLLCFlux"/>
        /// using the example moving isentropic vortex in an ideal gas.
        /// </summary>
        [Test]
        public static void IsentropicVortexIdealGasOptimizedHLLCTest() {
            Program<VortexControl> p = null;
            Application<VortexControl>._Main(
                new string[] { @"-c cs:CNS.Tests.IsentropicVortex.ControlFiles.IsentropicVortexIdealGasOptimizedHLLC()" },
                false,
                delegate() {
                    p = new IsentropicVortexTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 4e-3),
                Tuple.Create("L2ErrorPressure", 5e-3),
                Tuple.Create("L2ErrorEntropy", 3.15e-3));
        }

        /// <summary>
        /// Tests the CNS solver with the <see cref="Convection.RusanovFlux"/>
        /// using the example moving isentropic vortex in a stiffened gas.
        /// </summary>
        [Test]
        public static void IsentropicVortexStiffenedlGasRusanovTest() {
            Program<VortexControl> p = null;
            Application<VortexControl>._Main(
                new string[] { @"-c cs:CNS.Tests.IsentropicVortex.ControlFiles.IsentropicVortexStiffenedGasRusanov()" },
                false,
                delegate() {
                    p = new IsentropicVortexTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2e-2),
                Tuple.Create("L2ErrorPressure", 6e-2),
                Tuple.Create("L2ErrorEntropy", 2e-3));
        }

        /// <summary>
        /// Tests the CNS solver with the <see cref="Convection.RusanovFlux"/>
        /// using the example moving isentropic vortex in a covolume gas.
        /// </summary>
        [Test]
        public static void IsentropicVortexCovolumeGasRusanovTest() {
            Program<VortexControl> p = null;
            Application<VortexControl>._Main(
                new string[] { @"-c cs:CNS.Tests.IsentropicVortex.ControlFiles.IsentropicVortexCovolumeGasRusanov()" },
                false,
                delegate() {
                    p = new IsentropicVortexTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 3.6e-3),
                Tuple.Create("L2ErrorPressure", 5.1e-3),
                Tuple.Create("L2ErrorEntropy", 3.1e-3));
        }
    }
}
