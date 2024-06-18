﻿/* =======================================================================
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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.XDG.Quadrature;
using IntersectingQuadrature;
using static BoSSS.Foundation.XDG.XQuadFactoryHelper;

namespace BoSSS.Foundation.XDG {




    /// <summary>
    /// Auxiliary class that helps with the creation of XDG-quadrature schemes;
    /// instances can be obtained via <see cref="LevelSetTracker.GetXQuadFactoryHelper"/>.
    /// </summary>
    public class XQuadFactoryHelperAlgoim : XQuadFactoryHelperBase {


        public XQuadFactoryHelperAlgoim(LevelSetTracker.LevelSetData[] lsDatas) : base(lsDatas) {
           
        }

        public override IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefVol) {
            throw new NotImplementedException();
        }

        public override IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            throw new NotImplementedException();
        }

        public override IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol) {
            throw new NotImplementedException();
        }

        public override IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            throw new NotImplementedException();
        }


        public override IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref) {
            var algoimFactory = new AlgoimFactories(m_LevelSetDatas[levSetIndex], Kref);
            
            return algoimFactory.GetSurfaceFactory(); 
        }

        public override IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref) {
            bool negativeVolume;
            if (jmp == JumpTypes.Heaviside)
                negativeVolume = false;
            else if (jmp == JumpTypes.OneMinusHeaviside)
                negativeVolume = true;
            else
                throw new NotImplementedException();

            var algoimFactory = new AlgoimFactories(m_LevelSetDatas[levSetIndex], Kref, negativeVolume);

            return algoimFactory.GetVolumeFactory();
        }
    }




}
