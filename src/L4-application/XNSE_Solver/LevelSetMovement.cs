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
using ilPSP.Tracing;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Foundation.SpecFEM;
using System.Diagnostics;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.Smoothing;
using BoSSS.Solution.Tecplot;
using ilPSP;
using BoSSS.Platform;
using ilPSP.Utils;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// options that control level-set evolution
    /// </summary>
    public enum LevelSetEvolution {

        /// <summary>
        /// turn level-set evolution off
        /// </summary>
        None,

        /// <summary>
        /// Prescribed level-set, the initial value is set every time.
        /// </summary>
        Prescribed,

        /// <summary>
        /// evolution is described by an explicit Fourier representation
        /// </summary>
        Fourier,

        /// <summary>
        /// Extension velocity based fast marching algorithm.
        /// </summary>
        FastMarching,

        /// <summary>
        /// Level set evolution where the level set field is treated by a scalar convection.
        /// </summary>
        ScalarConvection,

        /// <summary>
        /// The Level-Set is moved by an extension velocity idea using a PDE
        /// ReInitialization is done by fast marching, no ReInit on Cut Cells!
        /// -> TODO
        /// </summary>
        ExtensionVelocity
    }

    class LevelSetMover {

        DGField[] Velocity;
        LevelSetTracker LsTrk;
        ScalarFieldHistory<SinglePhaseField> DGLevSet;
        XVelocityProjection VelocityFilter;
        ExplicitNonconservativeAdvection Advection;
        VectorFieldHistory<SinglePhaseField> FilteredVelocity;

        public LevelSetMover(DGField[] _Velocity,
            VectorFieldHistory<SinglePhaseField> FilteredVelocity,
            LevelSetTracker _LsTrk, 
            XVelocityProjection.CutCellVelocityProjectiontype _CutCellVelocityProjectiontype,
            ScalarFieldHistory<SinglePhaseField> _DGLevSet,
            IncompressibleBoundaryCondMap BcMap) {
            
            this.Velocity = _Velocity;
            this.LsTrk = _LsTrk;
            this.DGLevSet = _DGLevSet;
            this.FilteredVelocity = FilteredVelocity;
                       

            VelocityFilter = new XVelocityProjection(_LsTrk, Velocity, FilteredVelocity.Current);
            VelocityFilter.Config.CutCellVelocityProjectiontype = _CutCellVelocityProjectiontype;
            //Advection = new ImplicitNonconservativeAdvection(DGLevSet.Current, FilteredVelocity.Current, BcMap);
            //Advection = new ImplicitNonconservativeAdvection(DGLevSet.Current, FilteredVelocity.Current, BcMap, AssumeDivergenceFreeVelocity: true);
            Advection = new ExplicitNonconservativeAdvection(DGLevSet.Current,  FilteredVelocity, BcMap, AssumeDivergenceFreeVelocity: true);
        }

        public void Advect(double dt) {
            using (new FuncTrace()) {
               
                // do the evo
                this.VelocityFilter.Perform(this.Velocity, this.FilteredVelocity.Current);
                Advection.Advect(dt);
                
            }
        }
    }
}
