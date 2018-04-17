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

using System.Collections.Generic;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using System;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Advection {
    /// <summary>
    /// A class for obtaining the time integrated spatial discretization of the level set advection equation
    /// </summary>
    public class ScalarVelocityAdvection : ILevelSetAdvection{
        
        /// <summary>
        /// ctor 
        /// </summary>
        /// <param name="LevelSet"></param>
        /// <param name="ScalarExtVel"></param>
        /// <param name="e"><see cref="ExplicitEuler.ChangeRateCallback"/></param>
        public ScalarVelocityAdvection(LevelSetTracker LSTrk, SinglePhaseField LevelSet, SinglePhaseField ScalarExtVel,  ExplicitEuler.ChangeRateCallback e = null, bool nearfield = false) {
            GridDat = (GridData)( LevelSet.Basis.GridDat);
            D = GridDat.SpatialDimension;
            this.LSTrk = LSTrk;
            this.nearfield = nearfield;
            CreateAdvectionSpatialOperator(LevelSet, ScalarExtVel, e, nearfield ? LSTrk.Regions.GetNearFieldSubgrid(1) :null);
        }

        /// <summary>
        /// Internal context of the class which is connected to the outside via the constructor
        /// </summary>
        GridData GridDat;

        LevelSetTracker LSTrk;
        bool nearfield;

        ExplicitEuler TimeIntegrator;

        /// <summary>
        /// Creating the time integrated DG-FEM discretization of the level set advection equation
        /// </summary>
        /// <param name="LevelSet"></param>
        /// <param name="ExtensionVelocity"></param>
        /// <param name="e"></param>
        void CreateAdvectionSpatialOperator(SinglePhaseField LevelSet, SinglePhaseField ExtensionVelocity, ExplicitEuler.ChangeRateCallback e, SubGrid subGrid) {
            SpatialOperator SO;
            Func<int[], int[], int[], int> QuadOrderFunction = QuadOrderFunc.Linear();
            int D = LevelSet.GridDat.SpatialDimension;
            //FieldFactory<SinglePhaseField> fac = new FieldFactory<SinglePhaseField>(SinglePhaseField.Factory);
            //VectorField<SinglePhaseField> LevelSetGradient = new VectorField<SinglePhaseField>(D,
            //    LevelSet.Basis,fac);
            
            SO = new SpatialOperator(1, 1 , 1, QuadOrderFunction, new string[] { "LS", "S", "Result" });
            double PenaltyBase = ((double)((LevelSet.Basis.Degree + 1) * (LevelSet.Basis.Degree + D))) / ((double)D);
            SO.EquationComponents["Result"].Add(new ScalarVelocityAdvectionFlux(GridDat, PenaltyBase));
            SO.Commit();
            this.TimeIntegrator = new RungeKutta(RungeKuttaScheme.ExplicitEuler, SO, new CoordinateMapping(LevelSet), new CoordinateMapping(ExtensionVelocity), subGrid);

            // Performing the task e 
            if (e != null)
                this.TimeIntegrator.OnBeforeComputeChangeRate += e;
        }
            
        /// <summary>
        /// Spatial dimension
        /// </summary>
        protected int D;

        /// <summary>
        /// Perform Advection
        /// </summary>
        /// <param name="dt">TimestepSize</param>
        public void Advect(double dt) {
            if (nearfield) throw new NotImplementedException();
            TimeIntegrator.Perform(dt);
        }

        /// <summary>
        /// Do nothing
        /// </summary>
        public void FinishTimeStep() {
            //Do nothing;
        }
    }

    /// <summary>
    /// A class for calculating the numerical flux included in the spatial discretization of the level set advection equation.
    /// </summary>
    class ScalarVelocityAdvectionFlux : IEdgeForm, IVolumeForm {
        public ScalarVelocityAdvectionFlux(GridData GridDat,  double PenaltyBase) {
            D = GridDat.SpatialDimension;
            this.PenaltyLengthScales = GridDat.Cells.cj;
            this.PenaltyBase = PenaltyBase;
        }

        /// <summary>
        /// Spatial dimension
        /// </summary>
        protected int D;

        private MultidimensionalArray PenaltyLengthScales;

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "LS" };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { "S" };
            }
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV | TermActivationFlags.V;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.None;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.GradUxGradV;
            }
        }


        double PenaltyBase;


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            //penalty:
            // return GetPenalty(inp.jCellIn, inp.jCellOut, inp.GridDat.Cells.cj)*(_uA[0]-_uB[0])*(_vA-_vB);

            //do nothing:
            //return 0.0;

            //Bernauer, Herzog 2012, (4.1)
            double GradUTimesNormalIn = 0;
            double GradUTimesNormalOut = 0;
            double AbsGradUIn = 0;
            double AbsGradUOut = 0;
            for (int d = 0; d < D; d++) {
                GradUTimesNormalIn  += _Grad_uA[0, d] * inp.Normale[d];
                GradUTimesNormalOut += _Grad_uB[0, d] * inp.Normale[d];
                AbsGradUIn = _Grad_uA[0, d] * _Grad_uA[0, d];
                AbsGradUOut = _Grad_uB[0, d] * _Grad_uB[0, d];
            }
            AbsGradUIn = Math.Sqrt(AbsGradUIn);
            AbsGradUOut = Math.Sqrt(AbsGradUOut);

            //double DirectionIn = GradUTimesNormalIn  * inp.Parameters_IN[0];
            //double DirectionOut = GradUTimesNormalOut  * inp.Parameters_OUT[0];

            double DirectionIn = GradUTimesNormalIn  * _uA[0];
            double DirectionOut = GradUTimesNormalOut  * _uB[0];
            double Acc = 0;

            if (DirectionIn + DirectionOut > 0) {
                Acc =  (_uA[0] - _uB[0]) * (-_vB);
                return Acc;
            }
            else {
                Acc =   (_uA[0] - _uB[0]) * (_vA);
                return Acc;
            }

        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0.0;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            //return cpv / Math.Abs(GradU);
            double Acc = 0;
            for (int d = 0; d < D; d++) {
                Acc += GradU[0, d] * GradU[0, d];
            }
            return cpv.Parameters[0] * V * Math.Sqrt(Acc);
        }
    }
}