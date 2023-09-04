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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XESF.Fluxes {
    public class GodunovFlux_Interface : ILevelSetForm,
        //INonlinLevelSetForm_V, 
        IEquationComponentChecking {

        #region ILevelSetForm members
        public int LevelSetIndex {
            get;
            private set;
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
        
        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            StateVector stateIn = new StateVector(uA, material);
            StateVector stateOut = new StateVector(uB, material);

            return this.bulkFlux.InnerEdgeFlux(null, int.MinValue, stateIn, stateOut, ref inp.Normal, int.MaxValue) * (vA - vB);
        }
        #endregion

        #region INonlinLevelSetForm_V members
        //public void NonlinInternalEdge_V(ref EdgeFormParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_Vin, MultidimensionalArray Koeff_Vot) {
        //    Debug.Assert(inp.Len == 1);
        //    Debug.Assert(Koeff_Vin.GetLength(0) == inp.Len);
        //    Debug.Assert(Koeff_Vot.GetLength(0) == inp.Len);

        //    // Split coefficient array into species A and B
        //    MultidimensionalArray fA = Koeff_Vin;//.ExtractSubArrayShallow(-1, -1, 0);
        //    MultidimensionalArray fB = Koeff_Vot;//.ExtractSubArrayShallow(-1, -1, 1);
        //    Debug.Assert(fA.Lengths.ListEquals(fB.Lengths, (a, b) => a == b));

        //    MultidimensionalArray tmp = MultidimensionalArray.Create(fA.Lengths);

        //    // Pass only those things that are really needed
        //    //bulkFlux.InnerEdgeFlux(
        //    //    time: double.NaN,
        //    //    stateIn: new StateVector(material, uA[], new Vector(uA[1], uA[2]), uA[3]);
        //    //    );

        //    fA.Acc(+1.0, tmp);
        //    fB.Acc(-1.0, tmp);
        //}
        #endregion

        #region IEquationComponentChecking
        public bool IgnoreVectorizedImplementation => false;
        #endregion

        //private readonly LevelSetTracker levelSetTracker;

        private readonly GodunovFlux bulkFlux;

        private readonly Material material;

        public GodunovFlux_Interface(LevelSetTracker levelSetTracker, CompressibleControl control, IBoundaryConditionMap boundaryMap, Material material, FluxVariables flux, int component = int.MinValue, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
            //this.levelSetTracker = levelSetTracker;
            this.material = material;

            switch (flux) {
                case FluxVariables.Density:
                    this.bulkFlux = new GodunovFlux(control, boundaryMap, new EulerDensityComponent(), material);
                    break;
                case FluxVariables.Momentum:
                    this.bulkFlux = new GodunovFlux(control, boundaryMap, new EulerMomentumComponent(component, material.EquationOfState.HeatCapacityRatio, control.MachNumber, levelSetTracker.GridDat.SpatialDimension), material);
                    break;
                case FluxVariables.Energy:
                    this.bulkFlux = new GodunovFlux(control, boundaryMap, new EulerEnergyComponent(), material);
                    break;
                default:
                    throw new NotImplementedException();
            }

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }
    }
}
