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
using System.Linq;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;

namespace BoSSS.Foundation.XDG {

    public static class XDGUtils {


        /// <summary>
        /// Computes the Integral of a given function over the zero iso-contour of the Level-Set
        /// </summary>
        /// <param name="LsTrk">Level-Set tracker</param>
        /// <param name="func">function which is integrated</param>
        /// <param name="momentFittingVariant"></param>
        /// <param name="SchemeHelper">optional XQuadSchemeHelper</param>
        /// <param name="HMForder"></param>
        /// <returns>Integral of <param name="func">func</param> over all MPI processors</returns>
        static public double GetIntegralOverZeroLevelSet(LevelSetTracker LsTrk, ScalarFunctionEx func, XQuadFactoryHelper.MomentFittingVariants momentFittingVariant, int HMForder, XQuadSchemeHelper SchemeHelper = null) {
            using (new FuncTrace()) {
                if (LsTrk.LevelSets.Count != 1)
                    throw new NotImplementedException();


                if(SchemeHelper == null)
                    SchemeHelper = LsTrk.GetXDGSpaceMetrics(momentFittingVariant, HMForder, 1).XQuadSchemeHelper;
                        // new XQuadSchemeHelper(LsTrk, momentFittingVariant);

                // Classic HMF uses order+1 for Surface Integrals and additionally 1 order higher for the HMF system
                // e.g order-2 is the cached quad rule 
                if (momentFittingVariant == XQuadFactoryHelper.MomentFittingVariants.Classic)
                    HMForder -= 2;

                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                double force = 0;

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        func(i0, Length, QR.Nodes, EvalResult);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            force += ResultsOfIntegration[i, 0];
                    }
                ).Execute();

                double localForce = force;
                double globalForce;
                unsafe
                {
                    csMPI.Raw.Allreduce(
                        (IntPtr)(&localForce),
                        (IntPtr)(&globalForce),
                        1,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        csMPI.Raw._OP.SUM,
                        csMPI.Raw._COMM.WORLD);
                }

                return globalForce;
            }
        }
    }
}
