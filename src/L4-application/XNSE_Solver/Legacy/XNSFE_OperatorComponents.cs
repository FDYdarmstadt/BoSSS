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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSE_Solver.Legacy {
    public static partial class XOperatorComponentsFactory {

        //==============================================================
        // additional interface components for NSFE interface components
        //==============================================================


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE_withEvaporation_component(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config,
            int d, int D, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.MomentumEquationComponent(d);
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            double sigma = physParams.Sigma;


            // set components
            var comps = XOp.EquationComponents[CodName];


            if (config.isTransport) {
                if (!config.isMovingMesh) {
                    comps.Add(new MassFluxAtInterface(d, D, thermParams, sigma, config.isMovingMesh));
                    comps.Add(new ConvectionAtLevelSet_nonMaterialLLF(d, D, LsTrk, thermParams, sigma));
                    comps.Add(new ConvectionAtLevelSet_Consistency(d, D, LsTrk, -1, false, thermParams, sigma));
                } 
            } else {
                comps.Add(new MassFluxAtInterface(d, D, thermParams, sigma, config.isMovingMesh));
            }


            if (config.isViscous) {
                comps.Add(new ViscosityAtLevelSet_FullySymmetric_withEvap(LsTrk.GridDat.SpatialDimension, physParams.mu_A, physParams.mu_B, dntParams.PenaltySafety, d, thermParams, sigma));
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE_withEvaporation(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config, int D, LevelSetTracker LsTrk) {

            for (int d = 0; d < D; d++) {
                AddInterfaceNSE_withEvaporation_component(XOp, config, d, D, LsTrk);
            }
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceContinuityEq_withEvaporation(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config, int D, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.ContinuityEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");


            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set components
            var comps = XOp.EquationComponents[CodName];

            var divEvap = new DivergenceAtLevelSet_withEvaporation(D, LsTrk, -1, false, thermParams, config.getPhysParams.Sigma);
            comps.Add(divEvap);

        }

    }
}
