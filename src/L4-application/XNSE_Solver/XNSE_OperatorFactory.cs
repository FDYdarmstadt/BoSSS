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

using ilPSP.Utils;

using BoSSS.Foundation.XDG;

using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.newXSpatialOperator;


namespace BoSSS.Application.XNSE_Solver {

    public class XNSE_OperatorFactory : XOperatorFactoryBase {

        string[] CodNameSelected = new string[0];
        string[] DomNameSelected = new string[0];

        int HMFDegree;

        public XNSE_OperatorFactory(XNSE_OperatorConfiguration config, LevelSetTracker _LsTrk, int _HMFdegree, IncompressibleMultiphaseBoundaryCondMap BcMap) {

            LsTrk = _LsTrk;
            D = _LsTrk.GridDat.SpatialDimension;

            HMFDegree = _HMFdegree;

            // test input
            // ==========
            {
                if(config.DomBlocks.GetLength(0) != 2 || config.CodBlocks.GetLength(0) != 2)
                    throw new ArgumentException();

                if((config.physParams.mu_A <= 0) || (config.physParams.mu_A <= 0))
                        throw new ArgumentException();

                if((config.physParams.rho_A <= 0) || (config.physParams.rho_B <= 0))
                    throw new ArgumentException();

                if(_LsTrk.SpeciesNames.Count != 2)
                    throw new ArgumentException();
                if(!(_LsTrk.SpeciesNames.Contains("A") && _LsTrk.SpeciesNames.Contains("B")))
                    throw new ArgumentException();
            }

            // full operator:
            // ==============
            CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D)).Cat("div");
            Params = ArrayTools.Cat(
                VariableNames.Velocity0Vector(D),
                VariableNames.Velocity0MeanVector(D),
                "Curvature",
                (new string[] { "surfForceX", "surfForceY", "surfForceZ" }).GetSubVector(0, D),
                (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D),
                (new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D)),
                VariableNames.Temperature,
                "DisjoiningPressure");
            DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

            // selected part:
            if(config.CodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, D));
            if(config.CodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D, 1));

            if(config.DomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, D));
            if(config.DomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D, 1));


            // create Operator
            // ===============
            m_XOp = new XSpatialOperatorMk2(DomNameSelected, Params, CodNameSelected, (A, B, C) => _HMFdegree, null);

            // add components
            // ==============

            // species bulk components
            for(int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                // Navier Stokes equations
                XOperatorComponentsFactory.AddSpeciesNSE(m_XOp, CodName.GetSubVector(0, D), LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, config, LsTrk);

                // continuity equation
                XOperatorComponentsFactory.AddSpeciesContinuityEq(m_XOp, CodName[D], D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, config, LsTrk);
            }

            // interface components
            XOperatorComponentsFactory.AddInterfaceNSE(m_XOp, CodName.GetSubVector(0, D), BcMap, config, LsTrk);    // surface stress tensor
            XOperatorComponentsFactory.AddSurfaceTensionForce(m_XOp, CodName.GetSubVector(0, D), BcMap, config, LsTrk);     // surface tension force

            XOperatorComponentsFactory.AddInterfaceContinuityEq(m_XOp, CodName[D], D, BcMap, config, LsTrk);       // continuity equation



            m_XOp.Commit();
        }


        bool NormalsRequired = false;

        bool CurvatureRequired = false;

        bool U0meanrequired = false;


    }
}
